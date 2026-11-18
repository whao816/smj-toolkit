"""FHI-aims 引擎插件：PBE/NAO 全分子弛豫（默认引擎）。

写 ``geometry.in`` + ``control.in`` + ``sub_aims.sh``；collect 从
``geometry.in.next_step`` 生成 ``opt.xyz``，并核对 FHI-aims 的收敛标记与几何
合理性。SCF/弛豫精度、开壳层自旋固定等细节与 ``aims_tip_relax`` 同源。
"""

from __future__ import annotations

import json
from pathlib import Path

from .. import aims_out, cube, species as species_mod
from .._vendor.junction.geometry import read_aims_geometry, read_xyz, write_xyz
from .._vendor.junction.inputs import make_aims_molecule_opt_input
from ..checks import check_geometry, total_electrons
from ..errors import StageError
from ..provenance import parse_aims_versions, record_observed_software
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
        multiplicity = int(getattr(molecule, "multiplicity", 1) or 1)
        return "none" if multiplicity <= 1 else str(multiplicity - 1)

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
        spin_moment = float(multiplicity - 1) if multiplicity > 1 else None
        functional = getattr(molecule, "opt_functional", "pbe")
        basis = getattr(molecule, "opt_basis", "light")

        written = make_aims_molecule_opt_input(
            xyz, stage_dir, species,
            xc=functional,
            spin=self._spin(molecule),
            charge=float(getattr(molecule, "charge", 0) or 0),
            force_threshold=getattr(molecule, "opt_force_threshold", "1.e-2"),
            basis=basis,
            env=env,
            vdw=self._vdw(ctx),
            spin_moment=spin_moment,
        )
        geom = Path(written["geometry"])
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
                ("开壳层：已写 fixed_spin_moment %g（μB），弛豫过程总自旋固定，"
                 "不允许在步间改变自旋态" % spin_moment) if spin_moment
                else "闭壳层：不写 fixed_spin_moment",
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
