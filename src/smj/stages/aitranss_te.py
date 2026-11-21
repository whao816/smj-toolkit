"""aitranss_te：自洽场结果 → aitranss 透射谱（远端阶段）。

prepare 写 ``tcontrol``、``self.energy.in`` 与 ``sub_te.sh``：
- ``tcontrol`` 由 ``aims_scf`` 目录里的 ``mos.aims``/``geometry.in`` 与
  ``electrode_build/electrode_metadata.json`` 的表面索引生成；
- ``self.energy.in`` 按几何里实际 Au 数分层生成（左右以 metadata 的 tip 索引为准）。
aitranss 跑完产出 ``TE.dat``（collect 只做轻量核对，真正的量在 junction_analyze）。
"""

from __future__ import annotations

import json
from pathlib import Path

from ..junction.geometry import read_aims_geometry
from ..junction.inputs import make_aitranss_input
from ..junction.landauer import read_te_ordered
from ..junction.self_energy import make_self_energy
from ..junction.tcontrol import parse_efermi, parse_mos_nsaos
from ..errors import StageBlocked, StageError
from ..landauer_spin import read_te_channels
from ..provenance import (
    parse_aitranss_versions,
    record_observed_software,
    sha256_file,
    structural_fingerprint,
    verify_parent_geometry,
)
from ..spin import spin_state
from ._aims_artifacts import probe, problem
from ._helpers import est_record, submit_command
from .base import Collected, Prepared, Stage, register


@register
class AitranssTeStage(Stage):
    name = "aitranss_te"
    depends = ("aims_scf",)
    script = "sub_te.sh"
    input_hint = "消费 aims_scf/{aims.out,mos.aims,geometry.in}（--scf-dir 可显式覆盖）"
    outputs = ("TE.dat",)

    # -- prepare ----------------------------------------------------------- #
    def prepare(self, ctx, stage_dir: Path) -> Prepared:
        stage_dir = Path(stage_dir)
        state = spin_state(ctx.config)
        scf_dir = self._scf_dir(ctx)
        scf_geometry = scf_dir / "geometry.in"
        if not scf_geometry.is_file():
            raise StageBlocked(
                self.name, "缺 aims_scf/geometry.in（tcontrol 需要结几何）", [str(scf_geometry)])
        metadata = self._metadata(ctx)
        geometry = ctx.artifact("electrode_build", "geometry.in")
        parent_stage = "electrode_build"
        if not geometry.is_file():
            # 独立运行（--scf-dir）时上游目录可能不在：SCF 目录里的几何是同一份。
            geometry = scf_geometry
            parent_stage = "aims_scf"
        # P0-3：透射计算消费的几何必须与 metadata 记录的指纹一致；并且
        # tcontrol 的 $coord 与 self.energy.in 的电极划分必须来自同一份几何。
        self._verify_geometry_chain(geometry, scf_geometry, metadata)
        # 2026-11-20 复核 P1-4：老 run 的 collect 没有真实性校验，tcontrol 仍会
        # 直接读 basis-indices.out / omat.aims —— 在这里补一道前置检查。
        for name in ("basis-indices.out", "omat.aims"):
            path = scf_dir / name
            if not path.is_file():
                raise StageBlocked(self.name, "缺 aims_scf/%s" % name, [str(path)])
            reason = problem(probe(path))
            if reason is not None:
                raise StageBlocked(
                    self.name, "aims_scf/%s 不可用：%s（tcontrol 会读到空矩阵）"
                    % (name, reason), [str(path)])
        natoms = self.count_aims_atoms(scf_geometry)
        estimate = ctx.estimate(self.name, natoms)
        env = ctx.stage_env(self.name, estimate)
        transport = ctx.config.transport
        # eV 口径（默认 E_F ±2 eV）必须先拿到 SCF 的 efermi 才能平移到 Ha 绝对值；
        # Hartree 口径是旧的绝对能量三元组，不需要换算。
        efermi = None
        if transport.energy_window_unit == "eV":
            efermi = self._efermi(scf_dir)
        try:
            ener, estep, eend = transport.resolved_window(efermi)
        except ValueError as exc:
            raise StageBlocked(self.name, "能量窗无法换算成 Ha 绝对窗口: %s" % exc,
                               [str(scf_dir)]) from exc
        # 2026-11-20 复核 P1-2: 能量窗必须跨 E_F。整窗落在 E_F 一侧时费米窗权重全在同
        # 一侧，L0/L1/L2 不是 E_F 处的展开，S 的符号/量级都没有意义。这里在**生成输入前**
        # 拦住，不把无意义的谱跑出来再靠判据兜底。
        if efermi is not None and not (ener <= efermi <= eend):
            rel = transport.relative_window()
            raise StageBlocked(
                self.name,
                "能量窗未覆盖 E_F：E−E_F ∈ [%+.3f, %+.3f] eV，efermi=%r Ha。"
                "费米窗积分要求窗口跨越 E_F（默认 E_F±2 eV）；请改 "
                "transport.energy_window 或 energy_window_unit"
                % (rel[0], rel[2], efermi),
                [str(scf_dir)])
        nlayers = self._resolve_nlayers(ctx.config.electrode.nlayers, metadata)
        # 自旋极化 SCF 不写 mos.aims（改写成 alpha.aims/beta.aims），
        # make_tcontrol 没法自动探 nsaos，必须显式传。
        nsaos = self._spin_nsaos(scf_dir) if state.open_shell else None
        try:
            written = make_aitranss_input(
                scf_dir, stage_dir, metadata=metadata, nlayers=nlayers, nsaos=nsaos,
                ener=ener, estep=estep, eend=eend, efermi=efermi, env=env,
            )
        except (FileNotFoundError, KeyError, ValueError) as exc:
            raise StageBlocked(
                self.name, "aitranss_te 输入不完整（mos.aims / efermi / 表面索引）: %s" % exc,
                [str(scf_dir)]) from exc
        values = written.get("values") or {}
        spin_notes = self._patch_tcontrol(
            Path(written["tcontrol"]), state, stage_dir / "input_summary.json")
        # self.energy.in：左右电极划分以 metadata 的 tip 索引为准（不得按距离猜）。
        element = self._electrode_element(metadata, ctx)
        metal_count = sum(
            1 for atom in read_aims_geometry(geometry) if atom.element == element)
        if metal_count < 4 or metal_count % 2:
            raise StageBlocked(
                self.name, "结内 %s 数异常(需两侧对称且 >=4): %d" % (element, metal_count),
                [str(geometry)])
        try:
            self_energy = make_self_energy(
                geometry, stage_dir / "self.energy.in",
                au_per_electrode=metal_count // 2, metadata=metadata, element=element,
            )
        except ValueError as exc:
            raise StageError("aitranss_te 无法生成 self.energy.in: %s" % exc) from exc
        ctx.log("aitranss_te: 已写 tcontrol + self.energy.in（%d 原子，%d %s）"
                % (natoms, metal_count, element))
        if transport.energy_window_unit == "eV":
            r_start, r_step, r_end = transport.relative_window()
            window_note = (
                "能量窗 E−E_F ∈ [%+.3f, %+.3f] eV（步长 %g eV，%d 点）；efermi=%s Ha "
                "→ tcontrol ener=%s / estep=%s / eend=%s"
                % (r_start, r_end, r_step,
                   int(round((r_end - r_start) / r_step)) + 1, efermi,
                   values.get("ener", ener), values.get("estep", estep),
                   values.get("eend", eend)))
        else:
            window_note = ("能量窗（Ha 绝对）ener=%s / estep=%s / eend=%s"
                           % (values.get("ener", ener), values.get("estep", estep),
                              values.get("eend", eend)))
        return self.prepared(
            ctx, stage_dir,
            inputs=[Path(written["tcontrol"]), stage_dir / "self.energy.in",
                    Path(written["submit_script"]), stage_dir / "input_summary.json",
                    scf_geometry, Path(metadata)],
            resources=est_record(estimate),
            notes=[
                "%s，nlayers=%d" % (window_note, values.get("nlayers", nlayers)),
                "表面索引 l=%s r=%s（来自 %s）"
                % ([values.get("lsurc"), values.get("lsurx"), values.get("lsury")],
                   [values.get("rsurc"), values.get("rsurx"), values.get("rsury")],
                   Path(metadata).name),
                "self.energy.in：每侧 %d 个电极 %s（含 tip，tip 记 empty）"
                % (metal_count // 2, element),
                *spin_notes,
                "prepare 只写输入与作业脚本，未提交任何远端命令",
            ],
            commands=submit_command(self.script),
            parent_stage=parent_stage,
            parent_fingerprint=sha256_file(geometry),
        )

    # -- collect ----------------------------------------------------------- #
    def collect(self, ctx, stage_dir: Path) -> Collected:
        stage_dir = Path(stage_dir)
        te = ctx.require_stage_file(self.name, "TE.dat")
        try:
            rows, order = read_te_ordered(te)
        except ValueError as exc:
            raise StageError("aitranss_te 无法解析 TE.dat: %s" % exc) from exc

        notes: list[str] = []
        te_out = stage_dir / "te.out"
        require_te_out = bool(getattr(ctx.config.transport, "require_te_out", True))
        software: dict[str, object] = {}
        if not te_out.is_file():
            if require_te_out:
                raise StageError(
                    "aitranss_te 缺 te.out（aitranss 主输出）：TE.dat 无法复核。"
                    "若这是外部/历史 run，显式设 transport.require_te_out: false 豁免")
            notes.append("缺 te.out：已按 transport.require_te_out=false 显式豁免")
        else:
            try:
                te_text = te_out.read_text(encoding="utf-8", errors="replace")
            except OSError as exc:
                raise StageError("aitranss_te 无法读取 te.out: %s" % exc) from exc
            software = parse_aitranss_versions(te_text)
            if software:
                record_observed_software(ctx.layout.run_dir, self.name, software)
            else:
                notes.append("未从 te.out 解析到 aitranss 版本：远端软件版本未入档")
        if not (stage_dir / "self.energy.in").is_file():
            notes.append("缺 self.energy.in：结果不可复现")

        try:
            _, channels = read_te_channels(te)
        except ValueError:
            channels = 1
        # 2026-11-20 复审 P1-3: rows 已按 E−E_F 升序, order 是文件原始行序。
        if order != "ascending":
            notes.append(
                "TE.dat 能量列非升序(%s)：分析时已按 E-E_F 升序排序，S 符号不受行序影响"
                % order)
        # 2026-11-21 自由基支线 S3：开壳层必须真的有两列 α/β。只剩一列时
        # aitranss 已经把两个通道合并成均值，自旋极化不可能再恢复——默认阻断，
        # 只有显式 allow_single_channel 才放行（外部/历史数据）。
        state = spin_state(ctx.config)
        allow_single = bool(getattr(getattr(ctx.config, "transport", None),
                                    "allow_single_channel", False))
        single_channel_allowed = False
        if state.open_shell and channels < 2:
            if not allow_single:
                raise StageError(
                    "开壳层 TE.dat 只有 1 个自旋通道：α/β 已被 aitranss 合并成均值，"
                    "自旋极化无法恢复。请检查 tcontrol 是否把 $scfmo 换成 "
                    "$uhfmo_alpha/$uhfmo_beta、SCF 是否以 spin collinear 收敛、"
                    "nsaos 是否非零；若这份谱确实只能单列，显式设 "
                    "transport.allow_single_channel: true")
            single_channel_allowed = True
            notes.append(
                "开壳层 TE.dat 只有 1 个自旋通道：已按 transport.allow_single_channel=true "
                "放行，先导自旋极化不可用")
        results = {
            "aitranss_te": {
                "points": len(rows),
                "te_bytes": te.stat().st_size,
                "max_transmission": max(row[2] for row in rows),
                "transmission_channels": channels,
                "single_channel_allowed": single_channel_allowed,
                "energy_order": order,
                "has_te_out": te_out.is_file(),
                "require_te_out": require_te_out,
                "software": software or None,
            }
        }
        artifacts = ["TE.dat"] + (["te.out"] if te_out.is_file() else [])
        return self.collected(
            ctx, stage_dir, artifacts=artifacts, results=results, notes=notes)

    def _verify_geometry_chain(self, geometry: Path, scf_geometry: Path,
                               metadata: Path) -> None:
        """SCF 几何 == metadata/A2 几何 == self.energy.in 所依据的几何。"""
        verify_parent_geometry(self.name, scf_geometry, metadata)
        if Path(geometry).resolve() == Path(scf_geometry).resolve():
            return
        geo_hash = structural_fingerprint(geometry)
        scf_hash = structural_fingerprint(scf_geometry)
        if geo_hash is None or scf_hash is None:
            raise StageBlocked(
                self.name, "无法读取上游/SCF 几何，拒绝继续",
                [str(geometry), str(scf_geometry)])
        if geo_hash != scf_hash:
            raise StageBlocked(
                self.name,
                "SCF 几何与上游几何不一致（aims_scf/geometry.in=%s，%s=%s）；"
                "tcontrol 的 $coord 与 self.energy.in 的电极划分必须来自同一份几何"
                % (str(scf_hash)[:12], Path(geometry).name, str(geo_hash)[:12]),
                [str(scf_geometry), str(geometry)])

    # -- 内部 -------------------------------------------------------------- #
    def _spin_nsaos(self, scf_dir: Path) -> int:
        """开壳层 SCF 的 nsaos：FHI-aims 把 mos.aims 换成 alpha.aims/beta.aims。"""
        alpha = Path(scf_dir) / "alpha.aims"
        beta = Path(scf_dir) / "beta.aims"
        missing = [path.name for path in (alpha, beta) if not path.is_file()]
        if missing:
            raise StageBlocked(
                self.name,
                "开壳层结需要自旋极化 SCF 波函数 alpha.aims + beta.aims"
                "（FHI-aims 在 spin collinear 下不再写 mos.aims）；缺 %s"
                % "、".join(missing),
                [str(alpha), str(beta)])
        # 2026-11-20 复审 P1-4：is_file() 不等于可用；零字节/截断的通道会被
        # tcontrol 当成空矩阵读进去，这里按头部再验一遍（不依赖上游已验）。
        for path in (alpha, beta):
            info = probe(path)
            reason = problem(info)
            if reason is not None:
                raise StageBlocked(
                    self.name,
                    "开壳层 SCF 波函数不可用：%s；tcontrol 会读到空的自旋通道" % reason,
                    [str(path), str(alpha), str(beta)])
        try:
            return int(parse_mos_nsaos(alpha))
        except (OSError, ValueError) as exc:
            raise StageBlocked(
                self.name, "无法从 alpha.aims 解析 nsaos: %s" % exc, [str(alpha)]) from exc

    def _patch_tcontrol(self, control: Path, state, summary: Path) -> list[str]:
        """把 tcontrol 改成自旋极化口径；返回人读的 notes。

        * ``$scfmo file=mos.aims`` → ``$uhfmo_alpha file=alpha.aims`` +
          ``$uhfmo_beta file=beta.aims``（手册要求；自旋极化时 mos.aims 不存在）；
        * ``$valence_electrons`` 扣除结的净电荷（FHI-aims 的 element 价电子表
          只给中性原子求和，带电结不扣就不守恒）。
        """
        lines = control.read_text(encoding="utf-8").splitlines()
        out: list[str] = []
        notes: list[str] = []
        saw_scfmo = False
        valence_base: int | None = None
        adjusted: int | None = None
        for line in lines:
            stripped = line.strip()
            if state.open_shell and stripped.startswith("$scfmo"):
                saw_scfmo = True
                out.append("$uhfmo_alpha file=alpha.aims")
                out.append("$uhfmo_beta file=beta.aims")
                continue
            if state.charge and stripped.startswith("$valence_electrons"):
                fields = stripped.split()
                try:
                    valence_base = int(float(fields[-1]))
                except (IndexError, ValueError) as exc:
                    raise StageError(
                        "aitranss_te: 无法解析 tcontrol 的 $valence_electrons: %r"
                        % stripped) from exc
                adjusted = valence_base - state.charge
                if adjusted <= 0:
                    raise StageError(
                        "aitranss_te: $valence_electrons 扣除 charge=%+d 后为 %d ≤ 0"
                        "（中性求和 %d）；检查电极尺寸或电荷"
                        % (state.charge, adjusted, valence_base))
                out.append("$valence_electrons  %d" % adjusted)
                continue
            out.append(line)
        if state.open_shell and not saw_scfmo:
            raise StageError(
                "aitranss_te: tcontrol 缺 $scfmo 行，无法换成 $uhfmo_alpha/$uhfmo_beta")
        if out != lines:
            control.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")
        if state.open_shell:
            notes.append(
                "开壳层结：tcontrol 用 $uhfmo_alpha file=alpha.aims + "
                "$uhfmo_beta file=beta.aims 读两套自旋波函数（替代 $scfmo/mos.aims）")
        if adjusted is not None and valence_base is not None:
            notes.append("电荷 %+d：$valence_electrons 由 %d 调整为 %d"
                         % (state.charge, valence_base, adjusted))
        self._record_spin_summary(summary, state, adjusted)
        return notes

    @staticmethod
    def _record_spin_summary(path: Path, state, valence_electrons: int | None) -> None:
        try:
            data = json.loads(Path(path).read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return
        if not isinstance(data, dict):
            return
        data["spin_polarized"] = state.open_shell
        data["charge"] = float(state.charge)
        if state.open_shell:
            data["wavefunction_files"] = ["alpha.aims", "beta.aims"]
        if valence_electrons is not None:
            data["valence_electrons"] = int(valence_electrons)
        Path(path).write_text(
            json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")

    def _scf_dir(self, ctx) -> Path:
        override = ctx.override_path("scf_dir")
        if override is not None:
            if not override.is_dir():
                raise StageBlocked(
                    self.name, "显式 scf 目录不存在: %s" % override, [str(override)])
            return override
        return ctx.upstream(self.name, "aims.out").parent

    def _metadata(self, ctx) -> Path:
        override = ctx.override_path("metadata")
        if override is not None:
            if not override.is_file():
                raise StageBlocked(
                    self.name, "显式 metadata 不存在: %s" % override, [str(override)])
            return override
        # au_relax 已从主链路移除：metadata 以 electrode_build 为准。旧 run 可通过
        # --metadata 显式指定其 au_relax 快照。
        candidate = ctx.artifact("electrode_build", "electrode_metadata.json")
        if candidate.is_file():
            return candidate
        raise StageBlocked(
            self.name,
            "缺 electrode_build/electrode_metadata.json（tcontrol 需要六个表面索引，"
            "self.energy.in 需要左右 tip 索引）",
            [str(candidate)],
        )

    @classmethod
    def _efermi(cls, scf_dir: Path) -> float:
        """SCF 的费米能级（Ha）：eV 相对口径换算成 Ha 绝对窗口的锚点。

        缺了 fail-closed，不用占位值——否则整个透射窗会平移到错误位置。
        """
        try:
            return parse_efermi([scf_dir / "aims.out", scf_dir / "scf.out",
                                 scf_dir / "out"])
        except ValueError as exc:
            raise StageBlocked(
                cls.name,
                "缺 efermi：transport.energy_window 为 eV 相对口径（默认 E_F ±2 eV），"
                "必须用 SCF 的费米能级换算成 Ha 绝对窗口（aims_scf/aims.out 未收敛或"
                "未写费米能级？）",
                [str(scf_dir)]) from exc

    @staticmethod
    def _electrode_element(metadata: Path, ctx) -> str:
        """电极金属：以 electrode_metadata.json 的 electrode_element 为准，缺省 Au。

        旧 run 的 metadata 没有该键（0.2.0 之前只支持 Au），按 "Au" 处理；
        配置里的 electrode.element 只作为读不到 metadata 时的兜底。
        """
        try:
            data = json.loads(Path(metadata).read_text(encoding="utf-8"))
        except (OSError, ValueError):
            data = {}
        element = data.get("electrode_element") if isinstance(data, dict) else None
        return str(element or ctx.config.electrode.element)

    @staticmethod
    def _resolve_nlayers(configured: int | None, metadata: Path) -> int:
        """``electrode.nlayers`` 为 null 时用库文件建议值（手册 §5.3.2），再回退 3。"""
        if configured is not None:
            return int(configured)
        try:
            data = json.loads(Path(metadata).read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return 3
        suggested = (data.get("library_cluster") or {}).get("nlayers_suggested")
        try:
            return int(suggested) if suggested else 3
        except (TypeError, ValueError):
            return 3
