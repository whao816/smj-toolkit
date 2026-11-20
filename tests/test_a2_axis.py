"""A2 表面审计的输运轴回归: tip 连线优先于 linker 连线。

2026-11-13 的 2CN-656 验收 run: 分子弛豫后 S-C 键相对两侧 tip 连线倾斜约 10°,
audit 用 tip->linker 定向, 把 15 原子表面层投影成 1 个原子, 拒绝签署。这里用
一个 tip 与 linker 明显不共线的合成结, 锁住"用 tip 轴"的行为。
"""
from __future__ import annotations

from pathlib import Path

from smj.junction.a2_viz import _axis, surface_audit
from smj.junction.geometry import read_aims_geometry

GEOMETRY = """\
atom -5.0  0.0  0.0 Au
atom  5.0  0.0  0.0 Au
atom -4.0  1.0  0.0 S
atom  4.0 -1.0  0.0 S
atom -6.0 -1.0 -1.0 Au
atom -6.0  1.0 -1.0 Au
atom -6.0  0.0  1.0 Au
atom  6.0 -1.0 -1.0 Au
atom  6.0  1.0 -1.0 Au
atom  6.0  0.0  1.0 Au
"""

METADATA = {
    "left_tip_index": 1,
    "right_tip_index": 2,
    "left_linker_index": 3,
    "right_linker_index": 4,
    "left_surface_indices": [5, 6, 7],
    "right_surface_indices": [8, 9, 10],
    "left_backplane_indices": [5, 6, 7],
    "right_backplane_indices": [8, 9, 10],
}


def test_axis_prefers_tips_over_linkers(workdir: Path):
    geom = workdir / "geometry.in"
    geom.write_text(GEOMETRY, encoding="utf-8")
    atoms = read_aims_geometry(geom)
    axis, _, _ = _axis(atoms, METADATA)
    # tip 连线是 +x; linker 连线是 (8,-2,0)/sqrt(68), 二者明显不同
    assert axis[0] > 0.99
    assert abs(axis[1]) < 0.01


def test_surface_audit_accepts_tilted_linker(workdir: Path):
    geom = workdir / "geometry.in"
    geom.write_text(GEOMETRY, encoding="utf-8")
    result = surface_audit(geom, METADATA)
    assert result["ok"], [
        (c["name"], c["detail"]) for c in result["checks"] if not c["ok"]
    ]
    # 表面层必须仍是 3 个三元组原子所在层, 而不是被 linker 轴投影成 1 个
    for side in ("left", "right"):
        same = [c for c in result["checks"] if c["name"] == "same_surface_layer"]
        assert all(c["ok"] for c in same)
