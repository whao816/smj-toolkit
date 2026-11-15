"""几何合理性检查：重叠 / 断键 / 解离 / 未知元素。

背景：``molecule_opt`` 的 collect 原来只核"进程正常结束 + 文件存在"，
键断裂、质子转移、原子塌缩都能一路进 ``aims_tip_relax``。``smj.checks``
是补上的那层廉价闸门，这里钉住它的判定口径。
"""
from __future__ import annotations

import math

import pytest

from smj._vendor.junction.geometry import Atom
from smj.checks import anchor_contacts, check_geometry, min_atom_distance


def atoms(*rows) -> list[Atom]:
    return [Atom(el, float(x), float(y), float(z)) for el, x, y, z in rows]


H2 = atoms(("H", 0.0, 0.0, 0.0), ("H", 0.74, 0.0, 0.0))


def test_h2_equilibrium_passes_and_counts_one_bond():
    check = check_geometry(H2, H2)
    assert check["ok"] is True
    assert check["problems"] == []
    assert check["atom_count"] == 2
    assert check["bond_count"] == [1, 1]
    assert check["min_distance_A"] == pytest.approx(0.74)
    assert check["rmsd_A"] == pytest.approx(0.0)
    assert check["bond_check_skipped"] is False


def test_overlap_is_rejected():
    check = check_geometry(None, atoms(("H", 0.0, 0.0, 0.0), ("H", 0.1, 0.0, 0.0)))
    assert check["ok"] is False
    assert any("重叠" in p for p in check["problems"])


def test_dissociation_leaves_isolated_atoms():
    check = check_geometry(H2, atoms(("H", 0.0, 0.0, 0.0), ("H", 5.0, 0.0, 0.0)))
    assert check["ok"] is False
    assert any("孤立原子" in p for p in check["problems"])


def test_heavy_bond_breaking_is_rejected():
    before = atoms(("C", -0.77, 0.0, 0.0), ("C", 0.77, 0.0, 0.0))
    after = atoms(("C", -2.0, 0.0, 0.0), ("C", 2.0, 0.0, 0.0))
    check = check_geometry(before, after)
    assert check["ok"] is False
    assert any("重原子键消失" in p for p in check["problems"])


def test_lost_h_bond_only_warns():
    """质子从 O1 迁移到 O2：O1–H 键消失，但只算警告（不阻断）。"""
    before = atoms(
        ("O", 0.0, 0.0, 0.0), ("H", 0.96, 0.0, 0.0), ("H", -0.24, 0.93, 0.0),
        ("O", 2.8, 0.0, 0.0), ("H", 3.76, 0.0, 0.0), ("H", 2.56, 0.93, 0.0),
    )
    after = atoms(
        ("O", 0.0, 0.0, 0.0), ("H", 2.0, -0.6, 0.0), ("H", -0.24, 0.93, 0.0),
        ("O", 2.8, 0.0, 0.0), ("H", 3.76, 0.0, 0.0), ("H", 2.56, 0.93, 0.0),
    )
    check = check_geometry(before, after)
    assert check["ok"] is True, check["problems"]
    assert any("含 H 的键消失" in w for w in check["warnings"])


def test_element_composition_change_is_rejected():
    check = check_geometry(H2, atoms(("H", 0.0, 0.0, 0.0), ("He", 0.74, 0.0, 0.0)))
    assert check["ok"] is False
    assert any("元素组成变化" in p for p in check["problems"])


def test_atom_count_change_is_rejected():
    check = check_geometry(H2, atoms(("H", 0.0, 0.0, 0.0)))
    assert check["ok"] is False
    assert any("原子数变化" in p for p in check["problems"])


def test_non_finite_coordinates_are_rejected():
    check = check_geometry(None, atoms(("H", math.nan, 0.0, 0.0), ("H", 0.74, 0.0, 0.0)))
    assert check["ok"] is False
    assert any("非有限值" in p for p in check["problems"])


def test_unknown_element_skips_bond_graph_instead_of_guessing():
    """未知共价半径的元素不参与键图——宁可不判，也不乱判。"""
    pair = atoms(("Xx", 0.0, 0.0, 0.0), ("Xx", 1.0, 0.0, 0.0))
    check = check_geometry(pair, pair)
    assert check["ok"] is True
    assert check["bond_check_skipped"] is True
    assert check["bond_count"] is None
    assert any("未知共价半径" in w for w in check["warnings"])


def test_min_atom_distance_of_single_atom():
    assert min_atom_distance(atoms(("H", 0.0, 0.0, 0.0))) == (None, None)

def test_anchor_contacts_reports_distance_and_angle():
    """S 在原点、C 沿 +x（1.82 Å）、Au 在 100° 方向（2.42 Å）。"""
    ang = math.radians(100.0)
    ax, ay = 2.42 * math.cos(ang), 2.42 * math.sin(ang)
    system = atoms(
        ("S", 0.0, 0.0, 0.0),
        ("C", 1.82, 0.0, 0.0),
        ("Au", ax, ay, 0.0),
    )
    contacts = anchor_contacts(system)
    assert len(contacts) == 1
    item = contacts[0]
    assert item["anchor"] == 1 and item["anchor_element"] == "S"
    assert item["metal"] == 3 and item["metal_element"] == "Au"
    assert item["metal_distance_A"] == pytest.approx(2.42)
    assert item["neighbor"] == 2 and item["neighbor_element"] == "C"
    assert item["angle_deg"] == pytest.approx(100.0, abs=1e-6)


def test_anchor_contacts_ignores_far_or_missing_metal():
    far = atoms(("S", 0.0, 0.0, 0.0), ("C", 1.82, 0.0, 0.0),
                ("Au", 5.0, 0.0, 0.0))
    assert anchor_contacts(far) == []
    assert anchor_contacts(atoms(("S", 0.0, 0.0, 0.0))) == []