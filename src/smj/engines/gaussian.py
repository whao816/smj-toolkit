"""Gaussian 引擎插件：兼容旧 run 的 ``opt freq`` 分子优化。

写 ``mol.gjf`` + ``sub_g16.sh``；collect 解析 ``mol.log`` 的 Standard orientation。
Gaussian 分支只保证"旧 run 还能跑通"，新流程默认走 FHI-aims。
"""

from __future__ import annotations

from pathlib import Path

from ..junction.geometry import read_xyz, write_xyz  # noqa: F401
from ..junction.inputs import make_gaussian_input
from ..junction.spectra import parse_gaussian_log
from ..checks import check_geometry, total_electrons
from ..errors import StageError
from . import Engine, EngineCollected, EnginePrepared, register


@register
class GaussianEngine(Engine):
    name = "gaussian"
    script = "sub_g16.sh"
    outputs = ("opt.xyz", "mol.log", "gaussian_summary.json",
               "mol.gjf", "sub_g16.sh", "input_summary.json")
    resource_anchor = "molecule_opt"
    env_anchor = "molecule_opt"
    probes = ("gaussian",)

    def prepare(self, ctx, stage_dir: Path, *, stage: str, molecule, xyz: Path,
                natoms: int, electron_note: str, estimate: dict,
                env: dict) -> EnginePrepared:
        stage_dir = Path(stage_dir)
        written = make_gaussian_input(
            xyz, stage_dir, name="mol",
            charge=molecule.charge, multiplicity=molecule.multiplicity,
            method=molecule.method, basis=molecule.basis, env=env,
        )
        gjf = Path(written["gjf"])
        ctx.log("molecule_opt: 已写 %s（%d 原子，%s/%s）"
                % (gjf.name, natoms, molecule.method, molecule.basis))
        return EnginePrepared(
            inputs=[gjf, Path(written["submit_script"]),
                    stage_dir / "input_summary.json"],
            artifacts=["mol.gjf", "sub_g16.sh", "input_summary.json"],
            notes=[
                "%d 原子，%s/%s，nproc=%s，partition=%s"
                % (natoms, molecule.method, molecule.basis,
                   env.get("nproc"), env.get("partition")),
                electron_note,
                "prepare 只写输入与作业脚本，未提交任何远端命令",
            ],
            script=self.script,
        )

    def collect(self, ctx, stage_dir: Path, *, stage: str) -> EngineCollected:
        stage_dir = Path(stage_dir)
        log = ctx.require_stage_file(stage, "mol.log")
        xyz = stage_dir / "opt.xyz"
        summary_path = stage_dir / "gaussian_summary.json"
        try:
            summary = parse_gaussian_log(log, xyz_path=xyz,
                                         summary_path=summary_path)
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
        return EngineCollected(
            artifacts=["opt.xyz", "mol.log", "gaussian_summary.json"],
            results=results, notes=notes)
