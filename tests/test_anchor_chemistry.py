"""锚定化学：2.40–2.49 Å 是硫醇盐/表面配位带，不能套到硫醚 -SCH3 上。

2026-11-13 的验收 run 锚定基团是 -SCH3（S 与两个 C 成键的硫醚），
Au–S = 2.68 Å 被 `junction_geom_anchor_aus` 按硫醇盐带有源带判成偏离。
硫醚的 S–Au 是配位接触，不是硫醇盐的 Au–S 共价键，必须按锚定化学分流。
"""

from __future__ import annotations

from smj._vendor.junction.geometry import Atom, write_aims_geometry
from smj._vendor.junction.validate import validate_junction_geometry


def _geometry(workdir, *, chemistry: str):
    atoms = [
        Atom("S", 0.0, 0.0, 0.0),
        Atom("C", 1.80, 0.0, 0.0),
        Atom("Au", 0.0, 0.0, 2.68),   # 接触 Au
        Atom("Au", -5.0, 0.0, 0.0),   # 左表面
        Atom("Au", 5.0, 0.0, 0.0),    # 右表面
    ]
    if chemistry == "thioether":
        atoms.append(Atom("C", -1.80, 0.0, 0.0))       # S 与两个 C 成键
    elif chemistry == "thiol":
        atoms.append(Atom("H", -1.34, 0.0, 0.0))       # S-H 完整
    elif chemistry == "thiolate":
        pass                                          # 只与一个 C 成键
    else:
        raise ValueError(chemistry)

    geo = workdir / "geometry.in"
    write_aims_geometry(geo, atoms)
    metadata = {
        "left_surface_indices": [4],
        "right_surface_indices": [5],
        "left_linker_index": 1,
        "right_linker_index": 1,
        "left_tip_index": 1,
        "right_tip_index": 1,
    }
    return validate_junction_geometry(workdir, metadata=metadata, geometry=geo)


def _anchor_check(report):
    return next(c for c in report["checks"] if c["name"] == "junction_geom_anchor_aus")


def test_thioether_does_not_use_thiolate_band(workdir):
    report = _geometry(workdir, chemistry="thioether")
    check = _anchor_check(report)
    assert check["severity"] == "info"
    assert check["value"]["chemistry"] == "thioether"
    assert "2.40" not in check["detail"] or "不适用" in check["detail"]


def test_thiolate_still_uses_thiolate_band(workdir):
    report = _geometry(workdir, chemistry="thiolate")
    check = _anchor_check(report)
    assert check["severity"] == "warn"
    assert check["value"]["chemistry"] == "thiolate"


def test_intact_thiol_is_not_thiolate(workdir):
    report = _geometry(workdir, chemistry="thiol")
    check = _anchor_check(report)
    assert check["severity"] == "info"
    assert check["value"]["chemistry"] == "thiol"
