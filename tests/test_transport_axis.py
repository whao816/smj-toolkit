"""输运轴定向：左右电极质心优先，锚定原子连线只作兜底。

2026-11-20 ttf-release 实测：分子相对电极倾斜时，两锚定 S 的连线方向被分子内
几何主导（S6→S10 主要沿 y），而两电极实际沿 z 分离 19.7 Å。按锚定连线定向会把
合法的右表面索引（82 号，仅越过中线 0.25 Å）判到"错误一侧"。
"""
from __future__ import annotations

from smj.junction.geometry import Atom, write_aims_geometry
from smj.junction.validate import validate_junction_geometry


def _tilted_junction(workdir):
    """分子沿 y 倾斜，两电极沿 z 分开；表面索引本身完全合法。"""
    atoms = [
        Atom("S", 0.0, 0.0, 0.0),        # 1: 左锚定
        Atom("S", 0.0, 3.0, 0.0),        # 2: 右锚定（连线沿 y）
        Atom("Au", 0.0, 2.0, -8.0),      # 3
        Atom("Au", 0.0, 2.0, -8.5),      # 4: 左表面（沿 y 偏在中线右侧）
        Atom("Au", 0.0, 2.0, -9.0),      # 5
        Atom("Au", 0.0, 1.0, 8.0),       # 6
        Atom("Au", 0.0, 1.0, 8.5),       # 7: 右表面（沿 y 偏在中线左侧）
        Atom("Au", 0.0, 1.0, 9.0),       # 8
    ]
    geo = workdir / "geometry.in"
    write_aims_geometry(geo, atoms)
    metadata = {
        "left_surface_indices": [4],
        "right_surface_indices": [7],
        "left_linker_index": 1,
        "right_linker_index": 2,
    }
    return geo, metadata


def _surface_side(report):
    return next(c for c in report["checks"]
                if c["name"] == "junction_geom_surface_side")


def test_electrode_centroids_define_the_transport_axis(workdir):
    geo, metadata = _tilted_junction(workdir)
    metadata["left_electrode_indices"] = [3, 4, 5]
    metadata["right_electrode_indices"] = [6, 7, 8]

    check = _surface_side(validate_junction_geometry(
        workdir, metadata=metadata, geometry=geo))

    assert check["passed"] is True
    assert check["value"]["axis"] == "z"
    assert check["value"]["basis"] == "electrode_centroids"


def test_linkers_are_only_a_fallback(workdir):
    """缺电极索引时退回锚定连线——同一个几何在倾斜分子上会定向到 y。"""
    geo, metadata = _tilted_junction(workdir)

    check = _surface_side(validate_junction_geometry(
        workdir, metadata=metadata, geometry=geo))

    assert check["value"]["axis"] == "y"
    assert check["value"]["basis"] == "linkers"
    assert check["passed"] is False


def test_missing_both_sources_warns_instead_of_misjudging(workdir):
    geo, metadata = _tilted_junction(workdir)
    metadata.pop("left_linker_index")
    metadata.pop("right_linker_index")

    check = _surface_side(validate_junction_geometry(
        workdir, metadata=metadata, geometry=geo))

    assert check["severity"] == "warn"
    assert "left/right_electrode_indices" in check["detail"]
