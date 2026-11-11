"""gaussian_opt：分子 XYZ → Gaussian 几何优化（远端阶段）。

prepare 只写 ``mol.gjf`` + ``sub_g16.sh``（不提交、不等待）；
collect 解析 ``mol.log``，把最后一步 Standard orientation 落成 ``opt.xyz`` —— 这是
``aims_tip_relax`` 的唯一分子输入。
"""

from __future__ import annotations

from pathlib import Path

from .._vendor.junction.inputs import make_gaussian_input
from .._vendor.junction.spectra import parse_gaussian_log
from ..errors import StageError
from ._helpers import est_record, submit_command
from .base import Collected, Prepared, Stage, register


@register
class GaussianOptStage(Stage):
    name = "gaussian_opt"
    script = "sub_g16.sh"
    input_hint = "--xyz <分子.xyz>（或 smj.yaml 的 molecule.xyz）"
    outputs = ("opt.xyz", "mol.log", "gaussian_summary.json")

    # -- prepare ----------------------------------------------------------- #
    def prepare(self, ctx, stage_dir: Path) -> Prepared:
        stage_dir = Path(stage_dir)
        xyz = ctx.molecule_xyz()
        natoms = self.count_xyz_atoms(xyz)
        estimate = ctx.estimate(self.name, natoms)
        env = ctx.stage_env(self.name, estimate)
        molecule = ctx.config.molecule
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
