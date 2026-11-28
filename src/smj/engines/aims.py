"""FHI-aims 引擎插件：PBE/NAO 全分子弛豫（默认引擎）。

写 ``geometry.in`` + ``control.in`` + ``sub_aims.sh``；collect 从
``geometry.in.next_step`` 生成 ``opt.xyz``，并核对 FHI-aims 的收敛标记与几何
合理性。SCF/弛豫精度、开壳层自旋固定等细节与 ``aims_tip_relax`` 同源。
"""

from __future__ import annotations

import json
from pathlib import Path

from .. import aims_out, cube, species as species_mod
from ..junction.geometry import read_aims_geometry, read_xyz, write_xyz
from ..junction.inputs import make_aims_molecule_opt_input
from ..checks import check_geometry, total_electrons
from ..errors import StageError
from ..provenance import parse_aims_versions, record_observed_software
from ..spin import describe_spin_constraint, record_initial_charge, write_initial_charges
from . import Engine, EngineCollected, EnginePrepared, register


@register
class AimsEngine(Engine):
    name = "aims"
    script = "sub_aims.sh"
    outputs = ("aims.out", "geometry.in", "control.in", "sub_aims.sh",
               "input_summary.json", "opt.xyz")
    resource_anchor = "aims_tip_relax"
    env_anchor = "aims_tip_relax"
    skip_probe_note = ("opt_engine=aims：分子优化走 FHI-aims，"
                       "跳过 Gaussian 路径/许可检查")
    remote_files = ("geometry.in", "control.in", "sub_aims.sh",
                    "input_summary.json", "aims.out",
                    "geometry.in.next_step", ".dsh_done")

    @staticmethod
    def _spin(molecule) -> str:
        """FHI-aims 的 ``spin`` 只接受 ``none`` / ``collinear``（不接受数字）。

        多重度 M 由 ``fixed_spin_moment = M - 1 = 2S`` 表达；早期实现把 M−1
        直接写进 ``spin``（如 ``spin 2``），FHI-aims 会因关键字取值非法而停止。
        ``molecule.broken_symmetry``（M=1 的 BS 单重态）同样要 ``collinear``：
        独立的自旋通道允许反号的局域磁矩。
        """
        multiplicity = int(getattr(molecule, "multiplicity", 1) or 1)
        broken = bool(getattr(molecule, "broken_symmetry", False))
        return "collinear" if (multiplicity > 1 or broken) else "none"

    @staticmethod
    def _initial_moments(molecule, natoms: int) -> dict[int, float]:
        """校验并返回 ``molecule.initial_moments``（1-based 原子序号 -> μB）。"""
        raw = getattr(molecule, "initial_moments", None) or {}
        moments = {int(key): float(value) for key, value in dict(raw).items()}
        if not moments:
            return {}
        multiplicity = int(getattr(molecule, "multiplicity", 1) or 1)
        broken = bool(getattr(molecule, "broken_symmetry", False))
        if multiplicity <= 1 and not broken:
            raise StageError(
                "molecule.initial_moments 只在开壳层时有意义：请设 multiplicity > 1，"
                "或 molecule.broken_symmetry=true（M=1 的 BS 单重态）；"
                "当前 multiplicity=1 且未开破缺对称")
        out_of_range = sorted(index for index in moments if not 1 <= index <= natoms)
        if out_of_range:
            raise StageError(
                "molecule.initial_moments 的原子序号越界：结构共 %d 原子（1-based），"
                "越界序号 %s" % (natoms, out_of_range))
        return moments

    @staticmethod
    def _write_initial_moments(geometry: Path, moments: dict[int, float]) -> None:
        """把 per-atom ``initial_moment`` 子 tag 插到对应 ``atom`` 行之后。

        FHI-aims 的 atom 子 tag 必须**独占一行**（与 ``constrain_relaxation``
        同一种层级写法）；``read_aims_geometry`` 只认 ``atom`` 行，所以写在末尾
        同一行会让 collect 解析不到原子，必须单独成行。
        """
        if not moments:
            return
        lines = Path(geometry).read_text(encoding="utf-8").splitlines()
        out: list[str] = []
        index = 0
        for line in lines:
            out.append(line)
            parts = line.split()
            if parts and parts[0] == "atom":
                index += 1
                value = moments.get(index)
                if value is not None:
                    out.append("  initial_moment  %g" % value)
        if index < max(moments):
            raise StageError(
                "molecule_opt: geometry.in 只有 %d 个 atom 行，但 initial_moments "
                "指定到第 %d 号" % (index, max(moments)))
        Path(geometry).write_text(
            "\n".join(out) + "\n", encoding="utf-8", newline="\n")

    @staticmethod
    def _record_initial_moments(summary_path: Path, moments: dict[int, float]) -> None:
        """把初始磁矩补进 ``input_summary.json``（引擎自己写的输入清单）。"""
        try:
            data = json.loads(Path(summary_path).read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return
        if not isinstance(data, dict):
            return
        data["initial_moments"] = {str(key): value for key, value in sorted(moments.items())}
        Path(summary_path).write_text(
            json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")

    @staticmethod
    def _write_initial_charges(geometry: Path, charge: float) -> float | None:
        """非零电荷时把总电荷均分成 per-atom ``initial_charge`` 写进 geometry.in。

        FHI-aims 要求 ``charge q`` 等于 geometry.in 里 ``initial_charge`` 之和（不写
        这类 per-atom 子 tag 时 q 必须为 0），否则直接停算；项目使用的 221103 不认
        ``override_initial_charge_check`` 开关（集群实测 parse_control 停算），所以
        补 per-atom ``initial_charge``。均分只是初始密度猜测，不影响 SCF 自洽结果。
        """
        return write_initial_charges(geometry, charge)

    @staticmethod
    def _record_initial_charges(summary_path: Path, charge: float, per_atom: float) -> None:
        """把初始电荷分布记进 ``input_summary.json``，便于 truth chain 核对。"""
        record_initial_charge(summary_path, charge, per_atom)

    @staticmethod
    def _vdw(ctx) -> str | None:
        # 与 stages._helpers.vdw_mode 同义；插件层不 import smj.stages（避免循环）。
        enabled = bool(getattr(getattr(ctx.config, "electrode", None), "vdw", True))
        return "hirshfeld" if enabled else None

    def prepare(self, ctx, stage_dir: Path, *, stage: str, molecule, xyz: Path,
                natoms: int, electron_note: str, estimate: dict,
                env: dict) -> EnginePrepared:
        stage_dir = Path(stage_dir)
        species = ctx.species_dir()
        multiplicity = int(getattr(molecule, "multiplicity", 1) or 1)
        broken = bool(getattr(molecule, "broken_symmetry", False))
        spin_moment = float(multiplicity - 1) if multiplicity > 1 else None
        if broken and getattr(molecule, "fixed_spin_moment", False):
            spin_moment = 0.0
        functional = getattr(molecule, "opt_functional", "pbe")
        basis = getattr(molecule, "opt_basis", "light")

        charge = float(getattr(molecule, "charge", 0) or 0)
        written = make_aims_molecule_opt_input(
            xyz, stage_dir, species,
            xc=functional,
            spin=self._spin(molecule),
            charge=charge,
            force_threshold=getattr(molecule, "opt_force_threshold", "1.e-2"),
            basis=basis,
            env=env,
            vdw=self._vdw(ctx),
            spin_moment=spin_moment,
        )
        geom = Path(written["geometry"])
        control = Path(written["control"])
        charge_note = None
        if charge:
            per_atom = self._write_initial_charges(geom, charge)
            self._record_initial_charges(stage_dir / "input_summary.json", charge, per_atom)
            charge_note = (
                "电荷 %+g：FHI-aims 要求 charge 等于 geometry.in 的 initial_charge "
                "之和，已按原子数均分写入 per-atom initial_charge（%.6g e/原子，"
                "合计 %+g）；这只是初始密度猜测，不影响 SCF 自洽结果"
                % (charge, per_atom, charge))
        moments = self._initial_moments(molecule, natoms)
        if moments:
            self._write_initial_moments(geom, moments)
            self._record_initial_moments(stage_dir / "input_summary.json", moments)
        species_json = stage_dir / "species.json"
        species_notes: list[str] = []
        species_artifacts: list[str] = []
        try:
            payload = species_mod.summary(
                species, read_xyz(xyz), label=basis)
            species_mod.write_species_json(payload, species_json)
            species_artifacts.append("species.json")
            species_notes.append(
                "species 实际档位 %s（%d 个元素，逐文件 sha256 已记录）；"
                "配置标签 molecule.opt_basis=%s 只是标签"
                % (payload["actual_tier"], len(payload["elements"]), basis))
            species_notes.extend(payload["warnings"])
        except (OSError, ValueError) as exc:
            species_notes.append(
                "species 审计失败（%s）：本阶段仍按 species 目录计算，但无法核对"
                "实际档位与文件指纹" % exc)

        ctx.log("molecule_opt: 已写 %s（%d 原子，FHI-aims %s/%s，spin=%s）"
                % (geom.name, natoms, functional, basis, self._spin(molecule)))
        if broken:
            shell_note = (
                "破缺对称单重态初猜（M=1）：spin collinear + 反号 per-atom initial_moment；"
                "%s。需检查收敛后的局域磁矩或自旋密度，确认反号分布是否保留。"
                % describe_spin_constraint(spin_moment))
        elif spin_moment:
            shell_note = (
                "开壳层：spin collinear + fixed_spin_moment %g（μB = 2S = M−1）把总自旋"
                "固定住，弛豫过程不允许在步间改变自旋态" % spin_moment)
        else:
            shell_note = "闭壳层：spin none，不写 fixed_spin_moment / initial_moment"
        moment_notes: list[str] = []
        if moments:
            moment_notes.append(
                "初始磁矩（geometry.in per-atom initial_moment，1-based）：%s；合计 %g μB"
                % ("，".join("原子 %d=%g" % (index, value)
                            for index, value in sorted(moments.items())),
                   sum(moments.values())))
            target = 0.0 if broken else (spin_moment if spin_moment is not None else 0.0)
            if abs(sum(moments.values()) - target) > 1e-6:
                if broken:
                    moment_notes.append(
                        "破缺对称初猜的 initial_moments 合计 %g μB ≠ 0，请检查初猜符号"
                        % sum(moments.values()))
                else:
                    moment_notes.append(
                        "initial_moments 合计 %g μB ≠ M−1 = %g μB：fixed_spin_moment 仍会"
                        "强制总自旋，初猜分布只影响收敛路径与落到哪个解"
                        % (sum(moments.values()), spin_moment))
        elif spin_moment:
            moment_notes.append(
                "未设 molecule.initial_moments：spin collinear 下 FHI-aims 用默认初始"
                "自旋密度，开壳层 SCF 可能收敛到错误自旋分布（或更慢）；建议把 μB 加在"
                "自由基中心原子上（molecule.initial_moments，1-based 序号）")
        return EnginePrepared(
            inputs=[geom, Path(written["control"]),
                    Path(written["submit_script"]),
                    stage_dir / "input_summary.json"]
                   + ([species_json] if species_artifacts else []),
            artifacts=["geometry.in", "control.in", "sub_aims.sh",
                       "input_summary.json"] + species_artifacts,
            notes=[
                "%d 原子，FHI-aims %s/%s，spin=%s，charge=%s，ntasks=%s，partition=%s"
                % (natoms, functional, basis, self._spin(molecule),
                   getattr(molecule, "charge", 0), env.get("ntasks"),
                   env.get("partition")),
                "aims 引擎：无虚频检查（FHI-aims 未提供一键 freq）；几何来自 "
                "relax_geometry bfgs，最大力阈值 %s eV/Å（FHI-aims 判据是最大力"
                "分量，单位不是 Hartree/Bohr）；collect 会核对 "
                "\"Present geometry is converged\" 并做几何合理性检查"
                % getattr(molecule, "opt_force_threshold", "1.e-2"),
                electron_note,
                "SCF 精度显式写进 control.in（与 control_scf.in 同源）："
                "rho 1E-5、eev 1E-3 eV、etot 1E-6 eV、forces 5E-3 eV/Å，"
                "最多 1000 次迭代；旧版依赖 FHI-aims 默认（forces 为 not_checked）",
                shell_note,
                *moment_notes,
                *([charge_note] if charge_note else []),
                "prepare 只写输入与作业脚本，未提交任何远端命令",
            ] + species_notes,
            script=self.script,
        )

    def collect(self, ctx, stage_dir: Path, *, stage: str) -> EngineCollected:
        stage_dir = Path(stage_dir)
        out = ctx.require_stage_file(stage, "aims.out")
        text = out.read_text(encoding="utf-8", errors="replace")
        if "Have a nice day" not in text:
            raise StageError(
                "molecule_opt(aims) 未正常结束：aims.out 缺 \"Have a nice day\""
                "（%d 字节），拒绝把残缺几何交给下游 aims_tip_relax"
                % out.stat().st_size)

        relax = aims_out.relax_summary(text)
        threshold = getattr(ctx.config.molecule, "opt_force_threshold", "1.e-2")
        if relax["scf_not_converged"]:
            raise StageError(
                "molecule_opt(aims) 有 SCF 未收敛的弛豫步：aims.out 出现 "
                "\"self-consistency cycle ... not ... converged\"。那一步的力没有自洽"
                "保证，BFGS 会沿着不可信的方向走；请收紧/放宽 sc_accuracy_* 或加大 "
                "sc_iter_limit 后重跑（%d 步，最后一步最大力 %s eV/Å）"
                % (relax["relax_steps"], _force_str(relax["max_force_ev_per_A"])))
        if not relax["converged"]:
            raise StageError(
                "molecule_opt(aims) 几何未收敛：aims.out 有正常结束标记，但没有 "
                "\"Present geometry is converged\"（最后一步最大力 %s eV/Å，阈值 %s eV/Å，"
                "共 %d 步）。FHI-aims 弛豫达到步数上限时仍正常退出并写出 "
                "geometry.in.next_step，所以必须显式拒绝——非驻点几何会以初猜身份进入 "
                "aims_tip_relax，误差被后续局部弛豫掩盖"
                % (_force_str(relax["max_force_ev_per_A"]), threshold,
                   relax["relax_steps"]))

        next_step = stage_dir / "geometry.in.next_step"
        if not next_step.exists():
            raise StageError(
                "molecule_opt(aims) 缺 geometry.in.next_step：FHI-aims 未写出优化后几何，"
                "拒绝继续")
        atoms = read_aims_geometry(next_step)

        geometry_in = stage_dir / "geometry.in"
        prev = read_aims_geometry(geometry_in) if geometry_in.exists() else None
        check = check_geometry(prev, atoms)
        if not check["ok"]:
            raise StageError(
                "molecule_opt(aims) 几何合理性检查未通过（%d 原子）：%s。"
                "优化后几何已不是输入的那个分子，拒绝交给 aims_tip_relax"
                % (len(atoms), "；".join(check["problems"])))

        opt_xyz = stage_dir / "opt.xyz"
        write_xyz(opt_xyz, atoms, comment="FHI-aims %s/%s relax" % (
            getattr(ctx.config.molecule, "opt_functional", "pbe"),
            getattr(ctx.config.molecule, "opt_basis", "light")))
        software = parse_aims_versions(text)
        run_dir = getattr(getattr(ctx, "layout", None), "run_dir", stage_dir.parent)
        record_observed_software(run_dir, stage, software)
        results = {
            "molecule_opt": {
                "engine": "aims",
                "normal_termination": True,
                "atom_count": len(atoms),
                "converged": True,
                "force_threshold_ev_per_A": _as_float(threshold),
                "max_force_ev_per_A": relax["max_force_ev_per_A"],
                "relax_steps": relax["relax_steps"],
                "geometry_check": check,
            },
            "software": software,
        }
        orbitals = aims_out.frontier_orbitals(text)
        if orbitals is not None:
            results["molecule_opt"]["frontier_orbitals"] = orbitals
        charge = int(getattr(ctx.config.molecule, "charge", 0) or 0)
        try:
            results["molecule_opt"]["electron_count"] = total_electrons(atoms, charge)
        except KeyError:
            pass

        try:
            weights = cube.frontier_weights(stage_dir)
        except (ValueError, OSError) as exc:
            weights = []
            results["molecule_opt"]["frontier_orbital_weights_error"] = str(exc)
        if weights:
            results["molecule_opt"]["frontier_orbital_weights"] = weights

        control_path = stage_dir / "control.in"
        scf = aims_out.scf_settings(
            control_path.read_text(encoding="utf-8", errors="replace")
            if control_path.exists() else "")
        scf["converged_steps"] = relax["scf_converged_steps"]
        scf["relax_steps"] = relax["relax_steps"]
        results["molecule_opt"]["scf"] = scf

        notes = [
            "几何弛豫收敛：%d 步，最后一步最大力 %s eV/Å < 阈值 %s eV/Å"
            % (relax["relax_steps"], _force_str(relax["max_force_ev_per_A"]), threshold),
        ]
        species_path = stage_dir / "species.json"
        if species_path.exists():
            try:
                payload = json.loads(species_path.read_text(encoding="utf-8"))
            except (OSError, ValueError) as exc:
                notes.append("species.json 读取失败：%s" % exc)
            else:
                results["molecule_opt"]["species"] = payload
                notes.append(
                    "species 实际档位 %s（%d 个元素；基组与文件指纹见 results.species）"
                    % (payload.get("actual_tier", "unknown"),
                       len(payload.get("elements") or [])))
                notes.extend(payload.get("warnings") or [])
        else:
            notes.append("未找到 species.json（旧 run 或 gaussian 引擎）：无法核对"
                         "实际 species 档位与文件指纹")
        if prev is None:
            notes.append("未找到 prepare 写的 geometry.in：只做重叠/孤立检查，"
                         "跳过键图对比与位移 RMSD")
        notes.extend(check["warnings"])
        if scf.get("sc_accuracy_forces"):
            notes.append(
                "SCF 精度（control.in 显式）：rho=%s，eev=%s eV，etot=%s eV，"
                "forces=%s eV/Å，最多 %s 次迭代；SCF 收敛步数 %d / 弛豫步数 %d"
                % (scf.get("sc_accuracy_rho", "?"), scf.get("sc_accuracy_eev", "?"),
                   scf.get("sc_accuracy_etot", "?"), scf["sc_accuracy_forces"],
                   scf.get("sc_iter_limit", "?"), scf["converged_steps"],
                   scf["relax_steps"]))
        else:
            notes.append(
                "control.in 未显式设置 sc_accuracy_forces：这一步的力精度由 FHI-aims "
                "默认值决定（旧版默认 not_checked），几何阈值附近的力可信度打折")
        if scf["relax_steps"] and scf["converged_steps"] < scf["relax_steps"]:
            scf["missing_steps"] = scf["relax_steps"] - scf["converged_steps"]
            notes.append(
                "SCF 收敛标记 %d 次 < 弛豫步数 %d：有 %d 步没有 "
                "\"Self-consistency cycle converged\" 行（可能被日志截断），"
                "请人工核对"
                % (scf["converged_steps"], scf["relax_steps"], scf["missing_steps"]))
        if orbitals is not None:
            if orbitals["gap_ev"] is not None:
                notes.append("气相 HOMO-LUMO gap %.3f eV（PBE；只能作后续输运的对照，"
                             "不能直接当器件能隙）" % orbitals["gap_ev"])
            if orbitals["fractional_occupation"]:
                notes.append("aims.out 有分数占据（金属性）：HOMO/LUMO/gap 只作参考，"
                             "不是整数占据的带隙")
        if weights:
            notes.append(
                "锚定轨道权重（半径 %.1f Å 球内 |ψ|² 占全箱的比例）：%s。"
                "权重高 = 该前线轨道在锚定原子上有振幅（有耦合通道），"
                "接近 0 = 这个轨道不经过锚定原子"
                % (weights[0]["radius_A"], "；".join(
                    "%s %.1f%%" % (_orbit_label(item),
                                   100.0 * (item["anchor_fraction"] or 0.0))
                    for item in weights)))
        elif results["molecule_opt"].get("frontier_orbital_weights_error"):
            notes.append("未计算锚定轨道权重：%s"
                         % results["molecule_opt"]["frontier_orbital_weights_error"])
        if not software:
            notes.append("aims.out 未解析到版本信息（不伪造）")
        return EngineCollected(
            artifacts=["opt.xyz", "aims.out", "geometry.in", "geometry.in.next_step",
                       "control.in", "sub_aims.sh", "input_summary.json"],
            results=results, notes=notes)


def _orbit_label(item: dict) -> str:
    """cube 的报告标签：优先 control.in 里的 HOMO/LUMO，否则退到本征态号。"""
    orbit = item.get("orbit")
    if orbit:
        return str(orbit)
    return "本征态 %d" % int(item.get("eigenstate") or 0)


def _force_str(value: float | None) -> str:
    return "未知" if value is None else "%.3e" % value


def _as_float(value) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None
