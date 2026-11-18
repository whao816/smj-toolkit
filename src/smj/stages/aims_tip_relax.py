"""aims_tip_relax：opt.xyz → 加两个电极金属 tip（默认 Au）→ FHI-aims BFGS 弛豫（远端阶段）。

prepare 只写 ``geometry.in``（初始几何 + 两枚 tip）、``control.in``、
``sub_aims.sh``；BFGS 结束后 FHI-aims 自己写出 ``geometry.in.next_step``。
collect 量**弛豫后**的接触键长，避免"目标 2.40 Å"停留在写输入那一刻。
tip 元素随 ``electrode.element``（Au/Ag/Cu/…）变化，默认 Au 时逐字节不变。
"""

from __future__ import annotations

import json
from pathlib import Path

from .._vendor.junction.geometry import read_xyz
from .._vendor.junction.inputs import make_tip_relax_input, measure_au_s_bonds
from ..errors import StageBlocked, StageError
from ..structure_archive import safe_archive
from ._helpers import (est_record, resolve_anchor_indices, resolve_geometry_anchor_indices,
                       submit_command, vdw_mode)
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
            strip_thiol_h=electrode.strip_thiol_h,
            anchor_au_s_target=electrode.anchor_au_s_target,
            pin_contact=electrode.pin_contact,
            vdw=vdw_mode(electrode.vdw),
            tip_element=electrode.element,
        )
        summary = written["summary"]
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
        next_step = stage_dir / "geometry.in.next_step"
        if next_step.is_file():
            geometry, source = next_step, "geometry.in.next_step"
        else:
            # BFGS 未写出 next_step 时 FHI-aims 只留下初始几何：照样记录，但在 notes 里说清。
            geometry = ctx.require_stage_file(self.name, "geometry.in")
            source = "geometry.in"
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
