"""aims_tip_relax：opt.xyz → 加两个电极金属 tip（默认 Au）→ FHI-aims BFGS 弛豫（远端阶段）。

prepare 只写 ``geometry.in``（初始几何 + 两枚 tip）、``control.in``、
``sub_aims.sh``；BFGS 结束后 FHI-aims 自己写出 ``geometry.in.next_step``。
collect 量**弛豫后**的接触键长，避免"目标 2.40 Å"停留在写输入那一刻。
tip 元素随 ``electrode.element``（Au/Ag/Cu/…）变化，默认 Au 时逐字节不变。
"""

from __future__ import annotations

import json
from pathlib import Path

from ..junction.geometry import read_xyz
from ..junction.inputs import make_tip_relax_input, measure_au_s_bonds
from ..errors import StageBlocked, StageError
from ..spin import (DEFAULT_INITIAL_MOMENT, insert_control_lines, map_molecule_moments,
                    record_initial_charge, record_junction_spin, spin_control_lines,
                    spin_state, write_initial_charges, write_initial_moments)
from ..structure_archive import safe_archive
from ._aims_log import describe, relax_status
from ._helpers import (est_record, resolve_anchor_indices,
                       resolve_geometry_anchor_indices, submit_command, vdw_mode)
from .base import Collected, Prepared, Stage, register


@register
class AimsTipRelaxStage(Stage):
    name = "aims_tip_relax"
    depends = ("molecule_opt",)
    script = "sub_aims.sh"
    input_hint = "消费 molecule_opt/opt.xyz"
    outputs = ("geometry.in.next_step",)

    # -- prepare ----------------------------------------------------------- #
    def prepare(self, ctx, stage_dir: Path) -> Prepared:
        stage_dir = Path(stage_dir)
        # 开壳层结：tip 弛豫同样要 spin collinear + fixed_spin_moment=M−1，
        # 否则自由基的未配对电子在弛豫阶段就被 SCF 抹成闭壳层。
        state = spin_state(ctx.config)
        # --geometry 允许跳过 molecule_opt，直接把已优化分子结构喂给本阶段。
        xyz = ctx.upstream(self.name, "opt.xyz", override=ctx.override_path("geometry"))
        species = self._species_dir(ctx)
        electrode = ctx.config.electrode
        left, right = resolve_anchor_indices(
            xyz, electrode.left_anchor, electrode.right_anchor, stage=self.name)
        left_dir, right_dir = self._tip_directions(xyz, left, right,
                                                   getattr(electrode, "tip_direction", "axis"))
        natoms = self.count_xyz_atoms(xyz) + 2
        estimate = ctx.estimate(self.name, natoms)
        env = ctx.stage_env(self.name, estimate)
        written = make_tip_relax_input(
            xyz, stage_dir,
            left_anchor=left, right_anchor=right,
            au_tip_distance=electrode.au_tip_distance,
            left_direction=left_dir, right_direction=right_dir,
            species_dir=species, env=env,
            spin=state.spin, charge=float(state.charge),
            strip_thiol_h=electrode.strip_thiol_h,
            anchor_au_s_target=electrode.anchor_au_s_target,
            pin_contact=electrode.pin_contact,
            vdw=vdw_mode(electrode.vdw),
            tip_element=electrode.element,
        )
        summary = written["summary"]
        control = Path(written["control"])
        geometry = Path(written["geometry"])
        # 初始磁矩：molecule.initial_moments 是分子坐标系的 1-based 序号；tips 追加在
        # 分子之后，去质子只会删掉锚定 S–H 的氢，所以映射 = 减掉前面被删的氢。
        molecule = getattr(ctx.config, "molecule", None)
        moments = dict(getattr(molecule, "initial_moments", None) or {})
        default_moment: float | None = DEFAULT_INITIAL_MOMENT
        written_moments: dict[int, float] = {}
        if moments:
            if not state.open_shell:
                raise StageError(
                    "molecule.initial_moments 只在开壳层下有效：请设 multiplicity>1，"
                    "或 molecule.broken_symmetry=true（M=1 的 BS 单重态）")
            strip_info = summary.get("strip_thiol_h") or {}
            removed = (strip_info.get("removed_h_indices") or []
                       if isinstance(strip_info, dict) else [])
            try:
                mapped = map_molecule_moments(
                    moments, removed, molecule_atoms=natoms - 2)
                written_moments = write_initial_moments(geometry, mapped)
            except ValueError as exc:
                raise StageError("aims_tip_relax 写初始磁矩失败: %s" % exc) from exc
            default_moment = None
        # relax 模板只有 spin_moment_line 占位（make_tip_relax_input 没传），
        # fixed_spin_moment / default_initial_moment 在这里补。
        inserted = insert_control_lines(
            control, spin_control_lines(state, with_fixed_spin_moment=True,
                                        default_initial_moment=default_moment))
        if inserted or written_moments:
            record_junction_spin(stage_dir / "input_summary.json", state,
                                 default_initial_moment=default_moment,
                                 initial_moments=written_moments)
        if state.charge:
            per_atom = write_initial_charges(
                Path(written["geometry"]), float(state.charge))
            if per_atom is not None:
                record_initial_charge(
                    stage_dir / "input_summary.json", float(state.charge), per_atom)
        element = electrode.element
        anchor_element = summary.get("anchor_element") or electrode.left_anchor
        notes = [
            "%d 原子（分子 %d + 2 %s tip），锚定原子 %d/%d（%s）"
            % (natoms, natoms - 2, element, summary.get("left_anchor"),
               summary.get("right_anchor"), electrode.left_anchor),
            "tip 全原子弛豫（pin_contact=%s，direction=%s），%s–%s 初始放置距离 %s Å，species 目录 %s"
            % (summary.get("contact_pinned"),
               getattr(electrode, "tip_direction", "axis"),
               element, anchor_element,
               summary.get("placement_distance_A"), species),
            "prepare 只写输入与作业脚本，未提交任何远端命令",
        ]
        strip = summary.get("strip_thiol_h") or {}
        removed = int(strip.get("removed") or 0) if isinstance(strip, dict) else 0
        if removed:
            notes.append("去质子：移除 %d 个 %s–H 氢（原几何存档 opt_neutral.xyz）"
                         % (removed, anchor_element))
        elif electrode.strip_thiol_h is None and anchor_element != "S":
            notes.append("锚定原子 %s 非 S：strip_thiol_h=null 默认保留 H"
                         "（确需去质子请设 strip_thiol_h: true）" % anchor_element)
        if state.open_shell:
            if state.broken_symmetry:
                constraint = ("破缺对称单重态（M=1，不写 fixed_spin_moment，总自旋不钉）")
            elif state.fix_spin_moment:
                constraint = ("fixed_spin_moment %g（=2S=M−1）；弛豫全程固定总自旋"
                              % state.fixed_spin_moment)
            else:
                constraint = ("未写 fixed_spin_moment：解除总自旋约束，让 SCF 自寻自旋态"
                              "（release-fixed-spin 对照）")
            guess = ("default_initial_moment %g" % default_moment
                     if default_moment is not None else "per-atom initial_moment")
            notes.append("开壳层结：control.in 写 spin collinear + %s + %s（初猜）"
                         % (constraint, guess))
        if written_moments:
            notes.append(
                "初始磁矩（分子坐标系 %s）→ 结几何 1-based %s：%s"
                % (sorted(moments), sorted(written_moments),
                   "、".join("原子 %d=%.3g μB" % (index, value)
                            for index, value in written_moments.items())))
        if state.charge:
            notes.append(
                "电荷 %+d：geometry.in 已按原子数均分写 per-atom initial_charge"
                "（合计 = charge，满足 FHI-aims 一致性检查；初始电荷分布交给 "
                "SCF 自洽）" % state.charge)
        ctx.log("aims_tip_relax: 已写 geometry.in（%d 原子）" % natoms)
        return self.prepared(
            ctx, stage_dir,
            inputs=[Path(written["geometry"]), Path(written["control"]),
                    Path(written["submit_script"]), stage_dir / "input_summary.json"],
            resources=est_record(estimate),
            notes=notes,
            commands=submit_command(self.script),
        )

    # -- tip 方向 ---------------------------------------------------------- #
    @staticmethod
    def _tip_directions(xyz, left: int, right: int, mode: str):
        """按 ``electrode.tip_direction`` 给出两枚 tip 的外推方向。

        ``axis``(默认) 返回 (None, None) —— 交给 ``add_tip_au`` 沿两锚原子连线外推，
        与旧行为逐字节一致。``lone_pair`` 返回锚原子孤对方向的单位向量（以 S 为例：
        C–S–C 角平分线的反向），用于 S···S 连线与 S–C 键近共线的间位取代分子。
        """
        if mode == "axis":
            return None, None
        if mode != "lone_pair":
            raise StageError("未知的 electrode.tip_direction: %r" % (mode,))
        atoms = read_xyz(xyz)
        out = []
        for label, idx in (("left", int(left)), ("right", int(right))):
            anchor = atoms[idx - 1]
            ax, ay, az = anchor.coord()
            nbrs = [a for j, a in enumerate(atoms)
                    if j != idx - 1
                    and ((a.coord()[0] - ax) ** 2 + (a.coord()[1] - ay) ** 2
                         + (a.coord()[2] - az) ** 2) ** 0.5 < 1.95]
            if len(nbrs) != 2:
                raise StageError(
                    "tip_direction=lone_pair 需要锚定原子恰有 2 个近邻（C–S–C），"
                    "%s 锚定原子 %d 实得 %d 个" % (label, idx, len(nbrs)))
            (x1, y1, z1), (x2, y2, z2) = nbrs[0].coord(), nbrs[1].coord()
            vx, vy, vz = 2.0 * ax - x1 - x2, 2.0 * ay - y1 - y2, 2.0 * az - z1 - z2
            norm = (vx * vx + vy * vy + vz * vz) ** 0.5
            if norm < 1e-6:
                raise StageError("tip_direction=lone_pair: %s 锚定原子的孤对方向退化" % label)
            out.append((vx / norm, vy / norm, vz / norm))
        return out[0], out[1]

    # -- collect ----------------------------------------------------------- #
    def collect(self, ctx, stage_dir: Path) -> Collected:
        stage_dir = Path(stage_dir)
        # 2026-11-20 复核 P1-3：先验收弛豫收敛，再谈"弛豫后的键长"。
        # 不读日志时，未收敛的 next_step 会被当成弛豫结果直接送去建电极。
        status = relax_status(stage_dir / "aims.out")
        if status is None:
            raise StageError(
                "aims_tip_relax 缺少 aims.out：无法验收弛豫是否收敛，拒绝采集（%s）。"
                "若日志未回取，请先回取 aims.out 再 collect。" % stage_dir)
        if not status["converged"]:
            raise StageError(
                "aims_tip_relax 未收敛，拒绝采集：%s（日志 %s）。接触几何不是弛豫结果，"
                "不能用于建电极；请调整弛豫步数/收敛阈值后重跑。"
                % (describe(status), status["path"]))
        next_step = stage_dir / "geometry.in.next_step"
        if not next_step.is_file():
            raise StageError(
                "aims_tip_relax 日志显示已收敛，但缺少 geometry.in.next_step：产物不完整"
                "（可能只回取了日志），请重新回取或重跑")
        geometry, source = next_step, "geometry.in.next_step"
        left, right = self._anchors(ctx, stage_dir, geometry)
        summary = self._input_summary(stage_dir)
        tip_element = summary.get("tip_element") or ctx.config.electrode.element
        anchor_element = summary.get("anchor_element") or "S"
        try:
            distances = measure_au_s_bonds(geometry, left, right)
        except ValueError as exc:
            raise StageError("aims_tip_relax 无法量弛豫后的 %s–%s 键长: %s"
                             % (tip_element, anchor_element, exc)) from exc
        target = summary.get("anchor_au_s_target")
        payload = {
            "source": source,
            "atom_count": self.count_aims_atoms(geometry),
            "relaxed_au_s_A": distances,
            "left_anchor": left,
            "right_anchor": right,
            "contact_pinned": summary.get("contact_pinned"),
            "tip_element": tip_element,
            "anchor_element": anchor_element,
            "relax_converged": True,
            "relax_steps": status["steps"],
            "relax_last_marker": status["last_marker"],
            "relax_log_bytes": status["bytes"],
        }
        notes: list[str] = []
        if target is not None:
            try:
                target_value = float(target)
            except (TypeError, ValueError):
                target_value = None
            if target_value is not None:
                payload["target_A"] = target_value
                payload["delta_A"] = {side: distances[side] - target_value for side in distances}
        # 未收敛的日志在上面已被拒绝，这里 source 恒为 next_step；断言式保留提示。
        if source != "geometry.in.next_step":
            notes.append("未找到 geometry.in.next_step：BFGS 未写出最终几何，按初始 geometry.in 记录")
        # 「分子 + 两枚 tip」是后续装配电极的中间结构，同样进 saved_structures，
        # 这样从优化分子到完整分子结的每一版几何都可追溯。
        notes.extend(safe_archive(
            ctx.layout.run_dir, geometry, stage=self.name, kind="molecule_tips",
            engine="aims", results={"aims_tip_relax": payload},
            fmt="aims", natoms=payload["atom_count"]))
        return self.collected(
            ctx, stage_dir, artifacts=[source],
            results={"aims_tip_relax": payload}, notes=notes,
        )

    # -- 内部 -------------------------------------------------------------- #
    def _species_dir(self, ctx) -> Path:
        species = ctx.species_dir()
        if species is None:
            raise StageBlocked(
                self.name,
                "缺少 FHI-aims species 目录：请在 smj.yaml 设 cluster.species_dir 或 paths.species_defaults")
        if not species.is_dir():
            raise StageBlocked(self.name, "FHI-aims species 目录不存在: %s" % species, [str(species)])
        return species

    @staticmethod
    def _input_summary(stage_dir: Path) -> dict:
        path = Path(stage_dir) / "input_summary.json"
        if not path.is_file():
            return {}
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return {}
        return data if isinstance(data, dict) else {}

    def _anchors(self, ctx, stage_dir: Path, geometry: Path) -> tuple[int, int]:
        summary = self._input_summary(stage_dir)
        left, right = summary.get("left_anchor"), summary.get("right_anchor")
        if isinstance(left, int) and isinstance(right, int) and left and right:
            return left, right
        electrode = ctx.config.electrode
        try:
            return resolve_geometry_anchor_indices(
                geometry, electrode.left_anchor, electrode.right_anchor, stage=self.name)
        except ValueError as exc:
            raise StageError("aims_tip_relax 无法定位锚定原子: %s" % exc) from exc
