"""aims_scf：Au60 结的 FHI-aims 自洽场计算（远端阶段）。

prepare 把 ``electrode_build/geometry.in``（贴装后的结几何；``au_relax`` 已从主链路
移除）拷进本阶段目录，写 ``control.in`` + ``sub_aims.sh``；AIMS 跑完产出
``aims.out``、``basis-indices.out`` 与 ``mos.aims`` —— 后两者是 aitranss_te 的
必要输入。结几何在 SCF/输运阶段保持冻结，不再弛豫。
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

from ..junction.inputs import make_aims_scf_input
from ..errors import StageBlocked, StageError
from ..provenance import (
    parse_aims_versions,
    record_observed_software,
    sha256_file,
    verify_parent_geometry,
)
from ..spin import (DEFAULT_INITIAL_MOMENT, insert_control_lines, map_molecule_moments,
                    record_initial_charge, record_junction_spin, remove_control_lines,
                    spin_control_lines, spin_state, write_initial_charges,
                    write_initial_moments)
from ._aims_artifacts import probe, problem
from ._helpers import est_record, submit_command, vdw_mode
from .base import Collected, Prepared, Stage, register


@register
class AimsScfStage(Stage):
    name = "aims_scf"
    depends = ("electrode_build",)
    script = "sub_aims.sh"
    input_hint = "消费 electrode_build/geometry.in（冻结的结几何）"
    # aitranss 的 tcontrol 直接消费后三个文件；缺任何一个都跑不了透射。
    outputs = ("aims.out", "basis-indices.out", "omat.aims", "mos.aims")

    # -- prepare ----------------------------------------------------------- #
    def prepare(self, ctx, stage_dir: Path) -> Prepared:
        stage_dir = Path(stage_dir)
        # 结级自旋：SCF 必须与 tip 弛豫同一套 spin/charge/fixed_spin_moment，
        # 否则弛豫出的自旋态在 SCF 里被重新投影成闭壳层。
        state = spin_state(ctx.config)
        upstream = ctx.upstream(
            self.name, "geometry.in", override=ctx.override_path("geometry"))
        species = self._species_dir(ctx)
        geometry = stage_dir / "geometry.in"
        shutil.copy2(upstream, geometry)
        notes: list[str] = []
        # P0-3：SCF 必须用 electrode_build 产出并已由 A2 确认/锚定的那份几何。
        metadata_path = ctx.artifact("electrode_build", "electrode_metadata.json")
        parent_stage = "electrode_build"
        if metadata_path.is_file():
            verify_parent_geometry(self.name, upstream, metadata_path)
        else:
            notes.append("未找到 electrode_metadata.json：跳过父几何指纹校验（独立运行模式）")
        natoms = self.count_aims_atoms(geometry)
        estimate = ctx.estimate(self.name, natoms)
        env = ctx.stage_env(self.name, estimate)
        # 额外输出（默认空）：只有开壳层才允许要自旋密度——闭壳层没有自旋密度可写，
        # 静默接受会让人以为"自旋密度已产出"。
        transport = getattr(ctx.config, "transport", None)
        extra_outputs = tuple(getattr(transport, "extra_outputs", ()) or ())
        if "spin_density" in extra_outputs and not state.open_shell:
            raise StageError(
                "transport.extra_outputs 含 spin_density，但 molecule.multiplicity=%d "
                "是闭壳层：FHI-aims 没有共线自旋密度可写"
                % int(getattr(ctx.config.molecule, "multiplicity", 1) or 1))
        if extra_outputs:
            notes.append(
                "额外输出：%s（control.in 已写对应 output 行；collect 会核验产物真的出现，"
                "缺件拒绝 collected）" % " / ".join(extra_outputs))
        written = make_aims_scf_input(
            geometry, stage_dir, species_dir=species, env=env,
            spin=state.spin, charge=float(state.charge),
            fixed_spin_moment=state.fixed_spin_moment,
            vdw=vdw_mode(ctx.config.electrode.vdw),
            extra_outputs=extra_outputs,
            occupation_width=float(getattr(transport, "occupation_width", 0.01)),
            charge_mix_param=float(getattr(transport, "charge_mix_param", 0.2)),
            n_max_pulay=int(getattr(transport, "n_max_pulay", 10)),
        )
        # scf 模板已带 fixed_spin_moment 行；这里只补 FHI-aims 必需的初猜。
        control = Path(written["control"])
        molecule = getattr(ctx.config, "molecule", None)
        moments = dict(getattr(molecule, "initial_moments", None) or {})
        default_moment: float | None = DEFAULT_INITIAL_MOMENT
        written_moments: dict[int, float] = {}
        if moments:
            if not state.open_shell:
                raise StageError(
                    "molecule.initial_moments 只在开壳层（multiplicity>1）下有效")
            tip_summary = self._tip_summary(ctx)
            removed = (tip_summary.get("strip_thiol_h") or {}).get(
                "removed_h_indices") or []
            molecule_atoms = self._molecule_atoms(tip_summary)
            try:
                mapped = map_molecule_moments(
                    moments, removed, molecule_atoms=molecule_atoms)
                written_moments = write_initial_moments(geometry, mapped)
            except ValueError as exc:
                raise StageError("aims_scf 写初始磁矩失败: %s" % exc) from exc
            default_moment = None
        inserted = insert_control_lines(
            control, spin_control_lines(state, with_fixed_spin_moment=False,
                                        default_initial_moment=default_moment))
        if not state.fix_spin_moment and state.open_shell:
            # scf 模板没有"不写 fixed_spin_moment"的开关：写完后把这一行拿掉。
            remove_control_lines(control, ("fixed_spin_moment",))
        if inserted or written_moments or not state.fix_spin_moment:
            record_junction_spin(stage_dir / "input_summary.json", state,
                                 default_initial_moment=default_moment,
                                 initial_moments=written_moments)
        if state.charge:
            per_atom = write_initial_charges(geometry, float(state.charge))
            if per_atom is not None:
                record_initial_charge(
                    stage_dir / "input_summary.json", float(state.charge), per_atom)
        if state.open_shell:
            if state.fix_spin_moment:
                constraint = "fixed_spin_moment %g" % state.fixed_spin_moment
            else:
                constraint = "无 fixed_spin_moment（解除总自旋约束，对照解）"
            guess = ("default_initial_moment %g" % default_moment
                     if default_moment is not None else "per-atom initial_moment")
            notes.append(
                "开壳层结：spin collinear + %s + %s；FHI-aims 将把 mos.aims 写成 "
                "alpha.aims + beta.aims（自旋极化 SCF 不产出 mos.aims）"
                % (constraint, guess))
        if written_moments:
            notes.append(
                "初始磁矩（分子坐标系 %s）→ 结几何 1-based %s：%s"
                % (sorted(moments), sorted(written_moments),
                   "、".join("原子 %d=%.3g μB" % (index, value)
                            for index, value in written_moments.items())))
        scf_tuning = (float(getattr(transport, "occupation_width", 0.01)),
                      float(getattr(transport, "charge_mix_param", 0.2)),
                      int(getattr(transport, "n_max_pulay", 10)))
        if scf_tuning != (0.01, 0.2, 10):
            notes.append(
                "静态 SCF 收敛参数偏离默认：occupation_width=%g eV、"
                "charge_mix_param=%g、n_max_pulay=%d（默认 0.01 / 0.2 / 10）"
                % scf_tuning)
        if state.charge:
            notes.append(
                "电荷 %+d：geometry.in 已按原子数均分写 per-atom initial_charge"
                "（合计 = charge，满足 FHI-aims 一致性检查）" % state.charge)
        ctx.log("aims_scf: 已写 control.in（%d 原子，spin=%s）" % (natoms, state.spin))
        return self.prepared(
            ctx, stage_dir,
            inputs=[geometry, Path(written["control"]), Path(written["submit_script"]),
                    stage_dir / "input_summary.json"],
            resources=est_record(estimate),
            notes=[
                "%d 原子，ntasks=%s，partition=%s"
                % (natoms, env.get("ntasks"), env.get("partition")),
                *notes,
                "prepare 只写输入与作业脚本，未提交任何远端命令",
            ],
            commands=submit_command(self.script),
            parent_stage=parent_stage,
            parent_fingerprint=sha256_file(upstream),
        )

    # -- collect ----------------------------------------------------------- #
    def collect(self, ctx, stage_dir: Path) -> Collected:
        stage_dir = Path(stage_dir)
        output = ctx.require_stage_file(self.name, "aims.out")
        try:
            text = output.read_text(encoding="utf-8", errors="replace")
        except OSError as exc:
            raise StageError("aims_scf 无法读取 aims.out: %s" % exc) from exc
        if not text.strip():
            raise StageError("aims_scf/aims.out 为空：FHI-aims 没有产生可用输出")
        if "Have a nice day" not in text:
            raise StageError(
                "FHI-aims 未正常结束（aims.out 无 'Have a nice day'）；"
                "拒绝 collect 残缺的 SCF 结果")
        # 2026-11-20 复核 P1-4：日志写着正常结束，但 SCF 可能根本没收敛
        # （"SCF is not yet converged" / 迭代到上限后仍正常退出）。
        if "Self-consistency cycle converged" not in text:
            raise StageError(
                "aims_scf 的 SCF 未收敛：aims.out 里没有 'Self-consistency cycle converged'"
                "（可能命中 scf 迭代上限或提前退出）；tcontrol 会把未收敛的密度当输入")

        # tcontrol 的 $basis/$read_omat 两个口径共用；$scfmo 在自旋极化 SCF 里
        # 被 alpha.aims + beta.aims 取代——二选一缺了就 fail-closed。
        #
        # 2026-11-20 复审 P1-4（自旋通道）：必需波函数由**实际自旋模式**决定，
        # 并且逐个 probe。旧实现只看 is_file()、只 probe 一个波函数文件，
        # 零字节/截断的 beta.aims 因此被放行，tcontrol 读到的是一份空的自旋
        # 下降通道。配置声明 multiplicity>1 但磁盘上只有 mos.aims 同样拒绝——
        # 不得静默降级成闭壳层。
        present = {name: (stage_dir / name).is_file()
                   for name in ("mos.aims", "alpha.aims", "beta.aims")}
        spin_expected = bool(spin_state(ctx.config).open_shell
                             or present["alpha.aims"] or present["beta.aims"])
        if spin_expected:
            missing = [name for name in ("alpha.aims", "beta.aims")
                       if not present[name]]
            if missing:
                raise StageError(
                    "FHI-aims SCF 产物不完整：自旋极化 SCF 需要 alpha.aims + beta.aims "
                    "两份波函数（缺 %s）；aitranss_te 的 tcontrol 需要完整自旋通道"
                    % "、".join(missing))
            required_wf: tuple[str, ...] = ("alpha.aims", "beta.aims")
        elif present["mos.aims"]:
            required_wf = ("mos.aims",)
        else:
            raise StageError(
                "FHI-aims SCF 产物不完整：既没有 mos.aims（闭壳层），也没有 "
                "alpha.aims + beta.aims（自旋极化 SCF 会把 mos.aims 拆成两份）；"
                "aitranss_te 的 tcontrol 需要其中一组波函数文件")
        required = ("basis-indices.out", "omat.aims")
        missing = [name for name in required if not (stage_dir / name).is_file()]
        if missing:
            raise StageError(
                "FHI-aims SCF 产物不完整（缺 %s）；aitranss_te 的 tcontrol 需要 "
                "basis-indices.out + omat.aims" % "、".join(missing))
        # 存在 ≠ 可用：零字节/截断文件会让 tcontrol 读到空矩阵。
        observed: dict[str, dict] = {}
        for name in ("basis-indices.out", "omat.aims", *required_wf):
            info = probe(stage_dir / name)
            reason = problem(info)
            if reason is not None:
                raise StageError(
                    "aims_scf 产物不可用：%s；tcontrol 需要真实的 basis/omat/波函数文件"
                    % reason)
            observed[name] = info
        has_mos = "mos.aims" in observed
        has_alpha = "alpha.aims" in observed
        has_beta = "beta.aims" in observed

        # 额外输出（2026-11-20 复核 P1）：请求了就必须真的产出，不能"我开了开关"就算数。
        # spin_density 的产物文件名是 FHI-aims 文档化的 ``spin_density.cube``；
        # hirshfeld 打印在 aims.out（不是文件），这里只核验分析段确实出现。
        requested = tuple(getattr(getattr(ctx.config, "transport", None),
                                 "extra_outputs", ()) or ())
        extra = {"extra_outputs": list(requested)}
        if "spin_density" in requested:
            cube = stage_dir / "spin_density.cube"
            if not cube.is_file() or cube.stat().st_size <= 0:
                raise StageError(
                    "请求了 spin_density（control.in 有 'output cube spin_density'），"
                    "但 spin_density.cube 缺失或为空；拒绝把没有自旋密度证据的 run "
                    "标成 collected")
            extra["spin_density_cube_bytes"] = cube.stat().st_size
        if "hirshfeld" in requested:
            # 控制文件回显里也有 "vdw_correction_hirshfeld"，先扣掉再数。
            echo = text.lower().count("vdw_correction_hirshfeld")
            mentions = text.lower().count("hirshfeld") - echo
            if mentions <= 0:
                raise StageError(
                    "请求了 hirshfeld（control.in 有 'output hirshfeld'），但 aims.out "
                    "里没有 Hirshfeld 分析段；可能是该二进制不接受这个关键字")
            extra["hirshfeld_mentions"] = mentions

        notes: list[str] = []
        software = parse_aims_versions(text)
        if software:
            record_observed_software(ctx.layout.run_dir, self.name, software)
        else:
            notes.append("未从 aims.out 解析到 FHI-aims 版本/uuid：远端软件版本未入档")
        results = {
            "aims_scf": {
                "atom_count": self.count_aims_atoms(stage_dir / "geometry.in"),
                "normal_termination": "Have a nice day" in text,
                "scf_converged": "Self-consistency cycle converged" in text,
                "has_basis_indices": "basis-indices.out" in observed,
                "has_omat": "omat.aims" in observed,
                "has_mos": has_mos,
                "has_alpha": has_alpha,
                "has_beta": has_beta,
                "spin_unpolarized": has_mos,
                "basis_indices_bytes": observed["basis-indices.out"]["bytes"],
                "basis_functions": observed["basis-indices.out"]["size_value"],
                "omat_bytes": observed["omat.aims"]["bytes"],
                "omat_nsaos": observed["omat.aims"]["size_value"],
                "wavefunction_files": list(required_wf),
                "wavefunction_bytes": sum(observed[name]["bytes"]
                                          for name in required_wf),
                "log_bytes": len(text.encode("utf-8", "replace")),
                "software": software or None,
                **extra,
            }
        }
        artifacts = [name for name in ("aims.out", "control.in", "geometry.in",
                                       "basis-indices.out", "omat.aims",
                                       "mos.aims", "alpha.aims", "beta.aims",
                                       "spin_density.cube")
                     if (stage_dir / name).is_file()]
        return self.collected(
            ctx, stage_dir,
            artifacts=artifacts, results=results, notes=notes,
        )

    # -- 内部 -------------------------------------------------------------- #
    @staticmethod
    def _tip_summary(ctx) -> dict:
        """读 tip_relax 的 input_summary（去质子索引/分子原子数），缺失时返回空。"""
        try:
            path = ctx.artifact("aims_tip_relax", "input_summary.json")
        except Exception:  # noqa: BLE001 - artifact 只是布局查询，失败不影响主流程
            return {}
        try:
            data = json.loads(Path(path).read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return {}
        return data if isinstance(data, dict) else {}

    @staticmethod
    def _molecule_atoms(tip_summary: dict) -> int | None:
        """分子原子数 = tip_relax 的 atom_count − 两枚 tip；信息不全时返回 None。"""
        try:
            atom_count = int(tip_summary.get("atom_count"))
            tip_count = int(tip_summary.get("au_tip_count") or 2)
        except (TypeError, ValueError):
            return None
        if atom_count <= tip_count:
            return None
        return atom_count - tip_count

    def _species_dir(self, ctx) -> Path:
        species = ctx.species_dir()
        if species is None:
            raise StageBlocked(
                self.name,
                "缺少 FHI-aims species 目录：请在 smj.yaml 设 cluster.species_dir 或 paths.species_defaults")
        if not species.is_dir():
            raise StageBlocked(self.name, "FHI-aims species 目录不存在: %s" % species, [str(species)])
        return species
