"""2026-11-22 复审 P0：连接核对不能靠"表面索引没变"放行。

回归的真实缺陷：
1. ``surface_index_deviation`` 只比较记录/元数据里的索引集合；几何坐标被整体
   平移或某个表面原子被挪动时，集合完全不变 -> 偏差 0 -> 按 info 放行；
2. ``geometry.in`` 缺失时 ``check_confirmation`` 退回 metadata 里的快照指纹自比，
   旧记录被判有效；
3. ``connection_check_deviation`` 因此形同虚设，与 README 的"签名绑定几何
   sha256"承诺相反。
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from smj.config import ReviewConfig
from smj.junction.fingerprint import (
    CONFIRM_FILE,
    check_confirmation,
    record_confirmation,
)
from smj.review import (
    GEOMETRY_FILE,
    METADATA_FILE,
    ReviewResult,
    _connection_check_review,
    _connection_check_state,
    surface_coord_deviation,
)


def _geometry_text(z_shift: float = 0.0) -> str:
    atoms = [
        ("Au", 0.0, 0.0, 0.0),
        ("Au", 0.0, 0.0, 1.0),
        ("Au", 0.0, 0.0, 2.0),
        ("Au", 0.0, 0.0, 3.0),
        ("Au", 0.0, 0.0, 4.0),
        ("Au", 0.0, 0.0, 5.0 + z_shift),
    ]
    return "".join(
        "atom %.8f %.8f %.8f %s\n" % (x, y, z, element)
        for element, x, y, z in atoms)


def _metadata() -> dict:
    return {
        "left_surface_indices": [1, 2, 3],
        "right_surface_indices": [4, 5, 6],
        "left_anchor": 1,
        "right_anchor": 6,
    }


def _signed_build(root: Path) -> Path:
    build = root / "electrode_build"
    build.mkdir(parents=True)
    (build / GEOMETRY_FILE).write_text(_geometry_text(), encoding="utf-8")
    (build / METADATA_FILE).write_text(
        json.dumps(_metadata(), ensure_ascii=False), encoding="utf-8")
    record_confirmation(build, who="tester")
    return build


def _review(build: Path, **config) -> ReviewResult:
    result = ReviewResult(run_id="r")
    _connection_check_review(result, build, ReviewConfig(**config))
    return result.finalize()


def test_confirmation_records_surface_coordinates(workdir: Path):
    build = _signed_build(workdir)
    rec = json.loads((build / CONFIRM_FILE).read_text(encoding="utf-8"))
    assert rec["surface_coords"]["left"][0][:2] == [1, "Au"]
    assert check_confirmation(build)["valid"] is True


def test_moved_surface_atom_blocks_despite_identical_indices(workdir: Path):
    build = _signed_build(workdir)
    # 索引集合一个都没变，只是把最外层表面原子沿 z 挪了 0.05 Å。
    (build / GEOMETRY_FILE).write_text(_geometry_text(z_shift=0.05), encoding="utf-8")

    state = _connection_check_state(build)
    assert state["valid"] is False
    assert state["deviation"] == 0.0          # 旧的索引口径确实看不出问题
    assert state["deviation_geom_A"] == pytest.approx(0.05, abs=1e-9)

    result = _review(build, connection_check=True)
    assert result.status == "blocked"
    assert "connection_check_geometry_changed" in result.codes()


def test_missing_geometry_blocks_instead_of_self_comparing_metadata(workdir: Path):
    build = _signed_build(workdir)
    (build / GEOMETRY_FILE).unlink()

    checked = check_confirmation(build)
    assert checked["present"] is True
    assert checked["valid"] is False

    result = _review(build, connection_check=True)
    assert result.status == "blocked"
    assert "connection_check_stale_unquantified" in result.codes()


def test_legacy_record_without_coords_is_unquantified(workdir: Path):
    build = _signed_build(workdir)
    path = build / CONFIRM_FILE
    rec = json.loads(path.read_text(encoding="utf-8"))
    rec.pop("surface_coords")
    path.write_text(json.dumps(rec), encoding="utf-8")
    (build / GEOMETRY_FILE).write_text(_geometry_text(z_shift=0.05), encoding="utf-8")

    assert surface_coord_deviation(rec, build) is None
    result = _review(build, connection_check=True)
    assert result.status == "blocked"
    assert "connection_check_stale_unquantified" in result.codes()


def test_accept_deviation_is_the_only_waiver(workdir: Path):
    build = _signed_build(workdir)
    (build / GEOMETRY_FILE).write_text(_geometry_text(z_shift=0.05), encoding="utf-8")

    waived = _review(build, connection_check=True, accept_deviation=0.2)
    assert waived.accepted is True
    assert waived.status == "warning"
    assert "connection_check_deviation_accepted" in waived.codes()

    too_small = _review(build, connection_check=True, accept_deviation=0.01)
    assert too_small.status == "blocked"
    assert "connection_check_geometry_changed" in too_small.codes()


def test_untouched_geometry_still_passes(workdir: Path):
    build = _signed_build(workdir)
    result = _review(build, connection_check=True)
    assert result.status == "passed"
    assert "connection_check_confirmed" in result.codes()
