"""molecule_opt：分子 XYZ → 分子几何优化（远端阶段）。

引擎由 ``molecule.opt_engine`` 选择，具体实现在 :mod:`smj.engines`：

* ``aims``（默认）—— FHI-aims PBE/NAO 全分子优化，产物 ``geometry.in`` +
  ``control.in`` + ``sub_aims.sh``；collect 把 ``geometry.in.next_step`` 落成
  ``opt.xyz``。与后续 ``aims_tip_relax`` / ``aims_scf`` 同一套 FHI-aims 设置。
* ``gaussian``（兼容旧 run）—— 写 ``mol.gjf`` + ``sub_g16.sh``，collect 解析
  ``mol.log`` 的 Standard orientation。

本阶段只保留与引擎无关的公共流程：电子数/多重度奇偶校验、按注册表 dispatch、
把引擎的 prepare/collect 结果补成 ``Prepared``/``Collected``。两种引擎的契约
相同：``opt.xyz`` 是 ``aims_tip_relax`` 的唯一分子输入；失败时宁可不产出，也
不把残缺几何交给下游。
"""

from __future__ import annotations

from pathlib import Path

from ..junction.geometry import read_xyz
from ..checks import total_electrons
from ..engines import DEFAULT_ENGINE, all_engines, engine_for_config, get_engine
from ..errors import StageError
from ..structure_archive import archive_note, archive_structure
from ._helpers import est_record, submit_command
from .base import Collected, Prepared, Stage, register


@register
class MoleculeOptStage(Stage):
    name = "molecule_opt"
    # 元数据取缺省引擎；实际脚本由 prepare 按 molecule.opt_engine 写进 Prepared.script。
    script = get_engine(DEFAULT_ENGINE).script
    input_hint = "--xyz <分子.xyz>（或 smj.yaml 的 molecule.xyz）"
    outputs = tuple(dict.fromkeys(
        output for engine in all_engines() for output in engine.outputs))

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
        engine = engine_for_config(ctx.config)
        electrons = self._electron_state(xyz, molecule)
        electron_note = ("电子数 %d（charge=%d，multiplicity=%d）：奇偶校验通过"
                         % (electrons["electron_count"], electrons["charge"],
                            electrons["multiplicity"]))
        estimate = ctx.estimate(engine.resource_anchor, natoms)
        env = ctx.stage_env(engine.env_anchor, estimate)
        part = engine.prepare(
            ctx, stage_dir, stage=self.name, molecule=molecule, xyz=xyz,
            natoms=natoms, electron_note=electron_note,
            estimate=estimate, env=env)
        notes = list(part.notes)
        # 输入结构是"构建好"的那一版：先落盘归档，之后用户改 molecule.xyz
        # 也不影响本次 run 的审计链。
        notes.extend(self._archive(ctx, xyz, kind="input", engine=engine.name))
        return self.prepared(
            ctx, stage_dir,
            inputs=list(part.inputs),
            resources=est_record(estimate),
            script=part.script,
            artifacts=list(part.artifacts),
            notes=notes,
            commands=submit_command(part.script),
        )

    # -- collect ----------------------------------------------------------- #
    def collect(self, ctx, stage_dir: Path) -> Collected:
        stage_dir = Path(stage_dir)
        engine = engine_for_config(ctx.config)
        part = engine.collect(ctx, stage_dir, stage=self.name)
        notes = list(part.notes)
        opt_xyz = stage_dir / "opt.xyz"
        if opt_xyz.is_file():
            # 只有真正产出 opt.xyz 的 collect 才归档：失败分支宁可不留半截结构。
            notes.extend(self._archive(
                ctx, opt_xyz, kind="optimized",
                engine=engine.name, results=part.results))
        return self.collected(
            ctx, stage_dir,
            artifacts=list(part.artifacts),
            results=dict(part.results),
            notes=notes,
        )

    # -- 结构归档 ---------------------------------------------------------- #
    @staticmethod
    def _archive(ctx, source: Path, *, kind: str, engine: str,
                 results=None) -> list[str]:
        """归档一份 xyz；写盘失败只降级为 note，不阻断阶段（审计口径）。"""
        try:
            entry = archive_structure(
                ctx.layout.run_dir, source, stage="molecule_opt",
                kind=kind, engine=engine, results=results)
        except Exception as exc:  # noqa: BLE001 - 归档失败不能吞掉科学结果
            return ["结构归档失败（不影响本阶段）: %s" % exc]
        return [archive_note(entry)]
