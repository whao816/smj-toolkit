"""aims_tip_relax：opt.xyz → 加两个 Au tip → FHI-aims BFGS 弛豫（远端阶段）。

prepare 只写 ``geometry.in``（初始几何 + 两枚 Au tip）、``control.in``、
``sub_aims.sh``；BFGS 结束后 FHI-aims 自己写出 ``geometry.in.next_step``。
collect 量**弛豫后**的 Au–S 键长，避免"目标 2.40 Å"停留在写输入那一刻。
"""

from __future__ import annotations

import json
from pathlib import Path

from .._vendor.junction.inputs import make_tip_relax_input, measure_au_s_bonds
from ..errors import StageBlocked, StageError
from ._helpers import (est_record, resolve_anchor_indices, resolve_geometry_anchor_indices,
                       submit_command, vdw_mode)
from .base import Collected, Prepared, Stage, register


@register
class AimsTipRelaxStage(Stage):
    name = "aims_tip_relax"
    depends = ("gaussian_opt",)
    script = "sub_aims.sh"
    input_hint = "消费 gaussian_opt/opt.xyz"
    outputs = ("geometry.in.next_step",)

    # -- prepare ----------------------------------------------------------- #
    def prepare(self, ctx, stage_dir: Path) -> Prepared:
        stage_dir = Path(stage_dir)
        # --geometry 允许跳过 gaussian_opt，直接把已优化分子结构喂给本阶段。
        xyz = ctx.upstream(self.name, "opt.xyz", override=ctx.override_path("geometry"))
        species = self._species_dir(ctx)
        electrode = ctx.config.electrode
        left, right = resolve_anchor_indices(
            xyz, electrode.left_anchor, electrode.right_anchor, stage=self.name)
        natoms = self.count_xyz_atoms(xyz) + 2
        estimate = ctx.estimate(self.name, natoms)
        env = ctx.stage_env(self.name, estimate)
        written = make_tip_relax_input(
            xyz, stage_dir,
            left_anchor=left, right_anchor=right,
            au_tip_distance=electrode.au_tip_distance,
            species_dir=species, env=env,
            strip_thiol_h=electrode.strip_thiol_h,
            anchor_au_s_target=electrode.anchor_au_s_target,
            vdw=vdw_mode(electrode.vdw),
        )
        summary = written["summary"]
        notes = [
            "%d 原子（分子 %d + 2 Au tip），锚定原子 %d/%d（%s）"
            % (natoms, natoms - 2, summary.get("left_anchor"), summary.get("right_anchor"),
               electrode.left_anchor),
            "Au–S 放置距离 %s Å（pin_contact=%s），species 目录 %s"
            % (summary.get("placement_distance_A"), summary.get("contact_pinned"), species),
            "prepare 只写输入与作业脚本，未提交任何远端命令",
        ]
        strip = summary.get("strip_thiol_h") or {}
        removed = int(strip.get("removed") or 0) if isinstance(strip, dict) else 0
        if removed:
            notes.append("去质子：移除 %d 个 S–H 氢（原几何存档 opt_neutral.xyz）" % removed)
        ctx.log("aims_tip_relax: 已写 geometry.in（%d 原子）" % natoms)
        return self.prepared(
            ctx, stage_dir,
            inputs=[Path(written["geometry"]), Path(written["control"]),
                    Path(written["submit_script"]), stage_dir / "input_summary.json"],
            resources=est_record(estimate),
            notes=notes,
            commands=submit_command(self.script),
        )

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
        try:
            distances = measure_au_s_bonds(geometry, left, right)
        except ValueError as exc:
            raise StageError("aims_tip_relax 无法量弛豫后的 Au–S 键长: %s" % exc) from exc
        summary = self._input_summary(stage_dir)
        target = summary.get("anchor_au_s_target")
        payload = {
            "source": source,
            "atom_count": self.count_aims_atoms(geometry),
            "relaxed_au_s_A": distances,
            "left_anchor": left,
            "right_anchor": right,
            "contact_pinned": summary.get("contact_pinned"),
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
