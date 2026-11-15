"""molecule_opt：分子 XYZ → 分子几何优化（远端阶段）。

引擎由 ``molecule.opt_engine`` 选择：

* ``aims``（默认）—— FHI-aims PBE/NAO 全分子优化，产物 ``geometry.in`` +
  ``control.in`` + ``sub_aims.sh``；collect 把 ``geometry.in.next_step`` 落成
  ``opt.xyz``。与后续 ``aims_tip_relax`` / ``aims_scf`` 同一套 FHI-aims 设置，
  不再依赖 Gaussian。
* ``gaussian``（兼容旧 run）—— 写 ``mol.gjf`` + ``sub_g16.sh``，collect 解析
  ``mol.log`` 的 Standard orientation。

两种引擎的契约相同：``collect`` 后 ``opt.xyz`` 是 ``aims_tip_relax`` 的唯一
分子输入；失败时宁可不产出，也不把残缺几何交给下游。
"""

from __future__ import annotations

import json
from pathlib import Path

from .. import aims_out, cube, species as species_mod
from .._vendor.junction.geometry import read_aims_geometry, read_xyz, write_xyz
from .._vendor.junction.inputs import make_aims_molecule_opt_input, make_gaussian_input
from .._vendor.junction.spectra import parse_gaussian_log
from ..checks import check_geometry, total_electrons
from ..errors import StageError
from ..provenance import parse_aims_versions, record_observed_software
from ._helpers import est_record, submit_command, vdw_mode
from .base import Collected, Prepared, Stage, register


@register
class MoleculeOptStage(Stage):
    name = "molecule_opt"
    script = "sub_g16.sh"
    input_hint = "--xyz <分子.xyz>（或 smj.yaml 的 molecule.xyz）"
    outputs = ("opt.xyz", "mol.log", "gaussian_summary.json", "aims.out", "input_summary.json")

    @staticmethod
    def _engine(ctx) -> str:
        engine = str(getattr(ctx.config.molecule, "opt_engine", "gaussian") or "gaussian")
        return engine.lower()

    @staticmethod
    def _spin(molecule) -> str:
        multiplicity = int(getattr(molecule, "multiplicity", 1) or 1)
        return "none" if multiplicity <= 1 else str(multiplicity - 1)

    @staticmethod
    def _electron_state(xyz: Path, molecule) -> dict[str, int]:
        """电子数 ↔ 多重度奇偶校验。

        偶数电子只能是奇多重度（单重态/三重态…），奇数电子必须是偶多重度。
        写错时 SCF 会为一个"不存在的电子态"收敛出漂亮的结果，链上没有任何
        别的环节会发现它，所以这里在写输入之前就拒绝。
        """
        charge = int(getattr(molecule, "charge", 0) or 0)
        multiplicity = int(getattr(molecule, "multiplicity", 1) or 1)
        if multiplicity < 1:
            raise StageError("molecule_opt: multiplicity=%d 非法（至少为 1）"
                             % multiplicity)
        try:
            electrons = total_electrons(read_xyz(xyz), charge)
        except KeyError as exc:
            raise StageError(
                "molecule_opt: 无法计算电子数——元素 %s 不在价电子表"
                "（smj.checks.ATOMIC_NUMBERS）里，无法做电荷/多重度校验；"
                "请补表或确认输入" % exc.args[0]) from exc
        if electrons <= 0:
            raise StageError(
                "molecule_opt: charge=%d 下电子数 %d ≤ 0，检查 charge 或元素组成"
                % (charge, electrons))
        if (electrons % 2) == (multiplicity % 2):
            raise StageError(
                "molecule_opt: 电子数 %d（%s）与 multiplicity=%d 矛盾：偶数电子要求"
                "奇多重度，奇数电子要求偶多重度。写错会被 SCF 当成另一个电子态的"
                "分子算出\"收敛\"结果，下游无法察觉"
                % (electrons, "偶数" if electrons % 2 == 0 else "奇数", multiplicity))
        return {"electron_count": electrons, "charge": charge,
                "multiplicity": multiplicity}

    # -- prepare ----------------------------------------------------------- #
    def prepare(self, ctx, stage_dir: Path) -> Prepared:
        stage_dir = Path(stage_dir)
        xyz = ctx.molecule_xyz()
        natoms = self.count_xyz_atoms(xyz)
        molecule = ctx.config.molecule
        engine = self._engine(ctx)
        electrons = self._electron_state(xyz, molecule)
        electron_note = ("电子数 %d（charge=%d，multiplicity=%d）：奇偶校验通过"
                         % (electrons["electron_count"], electrons["charge"],
                            electrons["multiplicity"]))

        if engine == "aims":
            estimate = ctx.estimate("aims_tip_relax", natoms)
            env = ctx.stage_env("aims_tip_relax", estimate)
            species = ctx.species_dir()
            multiplicity = int(getattr(molecule, "multiplicity", 1) or 1)
            spin_moment = float(multiplicity - 1) if multiplicity > 1 else None
            written = make_aims_molecule_opt_input(
                xyz, stage_dir, species,
                xc=getattr(molecule, "opt_functional", "pbe"),
                spin=self._spin(molecule),
                charge=float(getattr(molecule, "charge", 0) or 0),
                force_threshold=getattr(molecule, "opt_force_threshold", "1.e-2"),
                basis=getattr(molecule, "opt_basis", "light"),
                env=env,
                vdw=vdw_mode(getattr(ctx.config.electrode, "vdw", True)),
                spin_moment=spin_moment,
            )
            geom = Path(written["geometry"])
            species_json = stage_dir / "species.json"
            species_notes: list[str] = []
            species_artifacts: list[str] = []
            try:
                payload = species_mod.summary(
                    species, read_xyz(xyz), label=getattr(molecule, "opt_basis", ""))
                species_mod.write_species_json(payload, species_json)
                species_artifacts.append("species.json")
                species_notes.append(
                    "species 实际档位 %s（%d 个元素，逐文件 sha256 已记录）；"
                    "配置标签 molecule.opt_basis=%s 只是标签"
                    % (payload["actual_tier"], len(payload["elements"]),
                       getattr(molecule, "opt_basis", "")))
                species_notes.extend(payload["warnings"])
            except (OSError, ValueError) as exc:
                species_notes.append(
                    "species 审计失败（%s）：本阶段仍按 species 目录计算，但无法核对"
                    "实际档位与文件指纹" % exc)
            ctx.log("molecule_opt: 已写 %s（%d 原子，FHI-aims %s/%s，spin=%s）"
                    % (geom.name, natoms,
                       getattr(molecule, "opt_functional", "pbe"),
                       getattr(molecule, "opt_basis", "light"),
                       self._spin(molecule)))
            return self.prepared(
                ctx, stage_dir,
                inputs=[geom, Path(written["control"]),
                        Path(written["submit_script"]), stage_dir / "input_summary.json"]
                       + ([species_json] if species_artifacts else []),
                resources=est_record(estimate),
                script="sub_aims.sh",
                artifacts=["geometry.in", "control.in", "sub_aims.sh",
                           "input_summary.json"] + species_artifacts,
                notes=[
                    "%d 原子，FHI-aims %s/%s，spin=%s，charge=%s，ntasks=%s，partition=%s"
                    % (natoms, getattr(molecule, "opt_functional", "pbe"),
                       getattr(molecule, "opt_basis", "light"), self._spin(molecule),
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
                commands=submit_command("sub_aims.sh"),
            )

        estimate = ctx.estimate(self.name, natoms)
        env = ctx.stage_env(self.name, estimate)
        written = make_gaussian_input(
            xyz, stage_dir, name="mol",
            charge=molecule.charge, multiplicity=molecule.multiplicity,
            method=molecule.method, basis=molecule.basis, env=env,
        )
        gjf = Path(written["gjf"])
        ctx.log("molecule_opt: 已写 %s（%d 原子，%s/%s）"
                % (gjf.name, natoms, molecule.method, molecule.basis))
        return self.prepared(
            ctx, stage_dir,
            inputs=[gjf, Path(written["submit_script"]), stage_dir / "input_summary.json"],
            resources=est_record(estimate),
            notes=[
                "%d 原子，%s/%s，nproc=%s，partition=%s"
                % (natoms, molecule.method, molecule.basis,
                   env.get("nproc"), env.get("partition")),
                electron_note,
                "prepare 只写输入与作业脚本，未提交任何远端命令",
            ],
            commands=submit_command(self.script),
        )

    # -- collect ----------------------------------------------------------- #
    def collect(self, ctx, stage_dir: Path) -> Collected:
        stage_dir = Path(stage_dir)
        if self._engine(ctx) == "aims":
            return self._collect_aims(ctx, stage_dir)
        return self._collect_gaussian(ctx, stage_dir)

    def _collect_aims(self, ctx, stage_dir: Path) -> Collected:
        out = ctx.require_stage_file(self.name, "aims.out")
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
        record_observed_software(run_dir, self.name, software)
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
        return self.collected(
            ctx, stage_dir,
            artifacts=["opt.xyz", "aims.out", "geometry.in", "geometry.in.next_step",
                       "control.in", "sub_aims.sh", "input_summary.json"],
            results=results, notes=notes)

    def _collect_gaussian(self, ctx, stage_dir: Path) -> Collected:
        log = ctx.require_stage_file(self.name, "mol.log")
        xyz = stage_dir / "opt.xyz"
        summary_path = stage_dir / "gaussian_summary.json"
        try:
            summary = parse_gaussian_log(log, xyz_path=xyz, summary_path=summary_path)
        except ValueError as exc:
            raise StageError("molecule_opt 无法解析 mol.log: %s" % exc) from exc
        if not summary.get("normal_termination"):
            # 残缺日志里的"最后一步"几何可能是中途结构：宁可不产出，也不静默交给下游。
            xyz.unlink(missing_ok=True)
            raise StageError(
                "molecule_opt 未正常结束(mol.log 无 \"Normal termination\", %d 字节): "
                "拒绝把残缺日志解析出的中间几何交给下游 aims_tip_relax"
                % log.stat().st_size
            )

        prev = None
        source = None
        try:
            source = ctx.molecule_xyz()
        except Exception:  # pragma: no cover - 配置缺 xyz 时 prepare 已经失败
            source = None
        if source is not None and Path(source).exists():
            prev = read_xyz(source)
        check = check_geometry(prev, read_xyz(xyz))
        if not check["ok"]:
            xyz.unlink(missing_ok=True)
            summary_path.unlink(missing_ok=True)
            raise StageError(
                "molecule_opt 几何合理性检查未通过（%d 原子）：%s。"
                "拒绝把该几何交给下游 aims_tip_relax"
                % (check["atom_count"], "；".join(check["problems"])))

        imaginary = int(summary.get("imaginary_frequency_count") or 0)
        notes: list[str] = []
        if prev is None:
            notes.append("未找到输入 xyz：只做重叠/孤立检查，跳过键图对比与位移 RMSD")
        notes.extend(check["warnings"])
        if imaginary:
            notes.append("Gaussian 虚频 %d 个：优化可能停在鞍点，请人工确认" % imaginary)
        results = {
            "gaussian": {
                "normal_termination": True,
                "final_scf_energy_hartree": summary.get("final_scf_energy_hartree"),
                "imaginary_frequency_count": imaginary,
                "atom_count": summary.get("atom_count"),
                "geometry_check": check,
            }
        }
        if prev is not None:
            try:
                results["gaussian"]["electron_count"] = total_electrons(
                    prev, int(getattr(ctx.config.molecule, "charge", 0) or 0))
            except KeyError:
                pass
        return self.collected(
            ctx, stage_dir, artifacts=["opt.xyz", "mol.log", "gaussian_summary.json"], results=results, notes=notes)


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
