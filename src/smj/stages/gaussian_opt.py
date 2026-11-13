"""gaussian_opt：分子 XYZ → 分子几何优化（远端阶段）。

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

from pathlib import Path

from .._vendor.junction.geometry import read_aims_geometry, write_xyz
from .._vendor.junction.inputs import make_aims_molecule_opt_input, make_gaussian_input
from .._vendor.junction.spectra import parse_gaussian_log
from ..errors import StageError
from ..provenance import parse_aims_versions, record_observed_software
from ._helpers import est_record, submit_command, vdw_mode
from .base import Collected, Prepared, Stage, register


@register
class GaussianOptStage(Stage):
    name = "gaussian_opt"
    script = "sub_g16.sh"
    input_hint = "--xyz <分子.xyz>（或 smj.yaml 的 molecule.xyz）"
    outputs = ("opt.xyz", "mol.log", "gaussian_summary.json")

    @staticmethod
    def _engine(ctx) -> str:
        engine = str(getattr(ctx.config.molecule, "opt_engine", "gaussian") or "gaussian")
        return engine.lower()

    @staticmethod
    def _spin(molecule) -> str:
        multiplicity = int(getattr(molecule, "multiplicity", 1) or 1)
        return "none" if multiplicity <= 1 else str(multiplicity - 1)

    # -- prepare ----------------------------------------------------------- #
    def prepare(self, ctx, stage_dir: Path) -> Prepared:
        stage_dir = Path(stage_dir)
        xyz = ctx.molecule_xyz()
        natoms = self.count_xyz_atoms(xyz)
        molecule = ctx.config.molecule
        engine = self._engine(ctx)

        if engine == "aims":
            estimate = ctx.estimate("aims_tip_relax", natoms)
            env = ctx.stage_env("aims_tip_relax", estimate)
            species = ctx.species_dir()
            written = make_aims_molecule_opt_input(
                xyz, stage_dir, species,
                xc=getattr(molecule, "opt_functional", "pbe"),
                spin=self._spin(molecule),
                charge=float(getattr(molecule, "charge", 0) or 0),
                force_threshold=getattr(molecule, "opt_force_threshold", "1.e-2"),
                basis=getattr(molecule, "opt_basis", "light"),
                env=env,
                vdw=vdw_mode(getattr(ctx.config.electrode, "vdw", True)),
            )
            geom = Path(written["geometry"])
            ctx.log("gaussian_opt: 已写 %s（%d 原子，FHI-aims %s/%s，spin=%s）"
                    % (geom.name, natoms,
                       getattr(molecule, "opt_functional", "pbe"),
                       getattr(molecule, "opt_basis", "light"),
                       self._spin(molecule)))
            return self.prepared(
                ctx, stage_dir,
                inputs=[geom, Path(written["control"]),
                        Path(written["submit_script"]), stage_dir / "input_summary.json"],
                resources=est_record(estimate),
                script="sub_aims.sh",
                artifacts=["geometry.in", "control.in", "sub_aims.sh",
                           "input_summary.json"],
                notes=[
                    "%d 原子，FHI-aims %s/%s，spin=%s，charge=%s，ntasks=%s，partition=%s"
                    % (natoms, getattr(molecule, "opt_functional", "pbe"),
                       getattr(molecule, "opt_basis", "light"), self._spin(molecule),
                       getattr(molecule, "charge", 0), env.get("ntasks"),
                       env.get("partition")),
                    "aims 引擎：无虚频检查（FHI-aims 未提供一键 freq）；几何来自 "
                    "relax_geometry bfgs，收敛阈值 %s" % getattr(
                        molecule, "opt_force_threshold", "1.e-2"),
                    "prepare 只写输入与作业脚本，未提交任何远端命令",
                ],
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
        ctx.log("gaussian_opt: 已写 %s（%d 原子，%s/%s）"
                % (gjf.name, natoms, molecule.method, molecule.basis))
        return self.prepared(
            ctx, stage_dir,
            inputs=[gjf, Path(written["submit_script"]), stage_dir / "input_summary.json"],
            resources=est_record(estimate),
            notes=[
                "%d 原子，%s/%s，nproc=%s，partition=%s"
                % (natoms, molecule.method, molecule.basis,
                   env.get("nproc"), env.get("partition")),
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
                "gaussian_opt(aims) 未正常结束：aims.out 缺 \"Have a nice day\""
                "（%d 字节），拒绝把残缺几何交给下游 aims_tip_relax"
                % out.stat().st_size)
        next_step = stage_dir / "geometry.in.next_step"
        if not next_step.exists():
            raise StageError(
                "gaussian_opt(aims) 缺 geometry.in.next_step：FHI-aims 未写出优化后几何，"
                "拒绝继续")
        atoms = read_aims_geometry(next_step)
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
            },
            "software": software,
        }
        notes = []
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
            raise StageError("gaussian_opt 无法解析 mol.log: %s" % exc) from exc
        if not summary.get("normal_termination"):
            # 残缺日志里的"最后一步"几何可能是中途结构：宁可不产出，也不静默交给下游。
            xyz.unlink(missing_ok=True)
            raise StageError(
                "gaussian_opt 未正常结束(mol.log 无 \"Normal termination\", %d 字节): "
                "拒绝把残缺日志解析出的中间几何交给下游 aims_tip_relax"
                % log.stat().st_size
            )
        imaginary = int(summary.get("imaginary_frequency_count") or 0)
        notes: list[str] = []
        if imaginary:
            notes.append("Gaussian 虚频 %d 个：优化可能停在鞍点，请人工确认" % imaginary)
        results = {
            "gaussian": {
                "normal_termination": True,
                "final_scf_energy_hartree": summary.get("final_scf_energy_hartree"),
                "imaginary_frequency_count": imaginary,
                "atom_count": summary.get("atom_count"),
            }
        }
        return self.collected(
            ctx, stage_dir, artifacts=list(self.outputs), results=results, notes=notes)
