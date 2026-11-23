"""thermal_device（构建声子器件）：把分子贴到两片 Au(111) 周期板上（本地阶段）。

**与电子链的 ``electrode_build`` 是两条独立的几何**，两者互不覆盖：

* ``electrode_build`` 造的是三角金字塔团簇 —— aitranss 的 ``self.energy`` 与
  ``basis-indices`` 都按"逐原子"给，团簇足够；但声子 NEGF 需要**半无限 lead**的
  表面格林函数，团簇既没有横向元胞、也没有可迭代的层间耦合，做不出自能。
* 本阶段造的是带横向元胞的二维周期板：面内 ``in_plane`` 个原胞（3×3 → 9 原子/
  层）、每侧 ``layers`` 层、沿 c 轴 ``vacuum`` 埃真空。这样层间力常数 V₁/V₂ 与
  面内耦合都有定义，decimation 才有东西可迭。

**去质子在这一步做**：``molecule_opt`` 产出的是中性二硫醇（S 上带 H），
直接贴 Au 就是报告 §6 陷阱⑥（用 SH 建模，Au–S 键连不起来）。``build_device``
内部先调 ``strip_thiol_h``，本阶段再把"去掉几个 H、哪几个原子"写进 metadata。
"""

from __future__ import annotations

from pathlib import Path

from ..errors import StageError
from ..junction.geometry import read_aims_geometry, read_xyz
from ..phonon.device import DeviceSpec, build_device, write_device
from ._thermal_common import upstream_file
from .base import Collected, Prepared, Stage, register


def _read_structure(path: Path):
    """按扩展名读结构：``.xyz`` 走 read_xyz，其余当 FHI-aims geometry.in。"""
    path = Path(path)
    if path.suffix.lower() == ".xyz":
        return read_xyz(path)
    return read_aims_geometry(path)


def spec_from_config(thermal) -> DeviceSpec:
    """``thermal.electrode.*`` → :class:`DeviceSpec`（两处共用，避免口径漂移）。"""
    electrode = thermal.electrode
    return DeviceSpec(
        a_lattice_A=float(electrode.a_lattice),
        in_plane=(int(electrode.in_plane[0]), int(electrode.in_plane[1])),
        electrode_layers=int(electrode.layers),
        device_layers=int(electrode.device_layers),
        vacuum_A=float(electrode.vacuum),
        au_s_distance_A=float(electrode.au_s_distance),
        gap_A=float(electrode.gap),
        element=str(electrode.element),
        anchor_element=str(electrode.anchor_element),
        strip_thiol_h=bool(electrode.strip_thiol_h),
        freeze_outer_layers=int(electrode.freeze_outer_layers),
    )


@register
class ThermalDeviceStage(Stage):
    name = "thermal_device"
    depends = ("molecule_opt",)
    needs_submit = False
    input_hint = "消费 molecule_opt/geometry.in（中性分子）或 --geometry"
    outputs = ("geometry.in", "device_metadata.json")

    def prepare(self, ctx, stage_dir: Path) -> Prepared:
        stage_dir = Path(stage_dir)
        thermal = getattr(ctx.config, "thermal", None)
        if thermal is None:
            raise StageError("配置里没有 thermal 段：声子热导链未启用")
        source = upstream_file(
            ctx, self.name, "geometry.in", "molecule_opt",
            override=ctx.override_path("geometry"))
        try:
            molecule = _read_structure(Path(source))
        except (OSError, ValueError) as exc:
            raise StageError("thermal_device 无法读取分子几何 %s: %s" % (source, exc)) from exc
        spec = spec_from_config(thermal)
        try:
            atoms, meta = build_device(molecule, spec)
        except ValueError as exc:
            raise StageError("thermal_device 建器件失败: %s" % exc) from exc
        geometry, metadata_path = write_device(stage_dir, atoms, meta)

        removed = list(meta.removed_h)
        atoms_per_layer = int(meta.spec["in_plane"][0]) * int(meta.spec["in_plane"][1])
        layers_per_side = int(meta.spec["electrode_layers"])
        molecule_atoms = len(meta.molecule_indices)
        device_au = len(meta.device_indices) - molecule_atoms
        notes = [
            "分子 %d 原子 → 贴装前去掉 %d 个 S–H → %d 原子"
            % (len(molecule), len(removed), molecule_atoms),
            "器件总 %d 原子：Au %d（每侧 %d 层 × %d 原子/层）+ 分子 %d"
            % (meta.n_atoms, len(meta.lead_indices) + device_au,
               layers_per_side, atoms_per_layer, molecule_atoms),
            "器件区 %d 原子（将参与有限位移），冻结 lead %d 原子"
            % (len(meta.device_indices), len(meta.lead_indices)),
            "面内 %.4f × %.4f Å，c 轴 %.4f Å（含 %.1f Å 真空）"
            % (meta.lattice[0][0], meta.lattice[1][1], meta.lattice[2][2],
               spec.vacuum_A),
        ]
        if len(removed) != 2:
            notes.append(
                "警告：去掉的 S–H 是 %d 个（预期 2 个）——锚定元素是 %s 时每个分子端"
                "应恰好 1 个 H；请核对输入分子" % (len(removed), spec.anchor_element))
        ctx.log("thermal_device: %d → %d 原子（%d 层/侧，%d 原子/层）"
                % (len(molecule), meta.n_atoms, layers_per_side, atoms_per_layer))
        return self.prepared(
            ctx, stage_dir,
            inputs=[Path(source), geometry, metadata_path],
            notes=notes,
            commands=[],
            parent_stage="molecule_opt",
        )

    def collect(self, ctx, stage_dir: Path) -> Collected:
        stage_dir = Path(stage_dir)
        geometry = ctx.require_stage_file(self.name, "geometry.in")
        meta_path = ctx.require_stage_file(self.name, "device_metadata.json")
        import json

        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise StageError("thermal_device 无法解析 device_metadata.json: %s" % exc) from exc
        expected = meta.get("n_atoms")
        actual = self.count_aims_atoms(geometry)
        if expected is not None and int(expected) != actual:
            raise StageError(
                "thermal_device 原子数不一致：metadata 记 %s，geometry.in 实际 %d"
                % (expected, actual))
        device_indices = meta.get("device_indices") or []
        lead_indices = meta.get("lead_indices") or []
        if len(device_indices) + len(lead_indices) != actual:
            raise StageError(
                "thermal_device 分区不完整：器件区 %d + lead %d ≠ 总原子 %d"
                % (len(device_indices), len(lead_indices), actual))
        results = {
            "thermal_device": {
                "n_atoms": actual,
                "molecule_atoms": len(meta.get("molecule_indices") or []),
                "removed_h": meta.get("removed_h"),
                "atoms_per_layer": (int(meta["spec"]["in_plane"][0])
                                    * int(meta["spec"]["in_plane"][1])),
                "layers_per_side": meta["spec"]["electrode_layers"],
                "n_device_atoms": len(device_indices),
                "n_lead_atoms": len(lead_indices),
                "n_displacements": 6 * len(device_indices),
                "n_constrained": len(meta.get("constrained_indices") or []),
                "lattice": meta.get("lattice"),
            }
        }
        return self.collected(
            ctx, stage_dir, artifacts=list(self.outputs), results=results)
