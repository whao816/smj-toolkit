"""声子热导链的两个 control.in 模板(thermal_relax / thermal_force)。

为什么值得单独钉死(方案附录 A-3):FHI-aims 里打印力的关键字是
``compute_forces``,写成 ``calculate_forces`` **不会报错**,只会静默不给力 ——
整个 324 次位移扫描算完才发现 aims.out 里没有力块。所以这里数的是关键字出现的
次数,不是「包含」。
"""
from __future__ import annotations

import pathlib

import pytest

from smj.junction.geometry import Atom, write_aims_geometry
from smj.junction.tcontrol import make_control

_ATOMS = [
    Atom("Au", 0.0, 0.0, 0.0),
    Atom("S", 0.0, 0.0, 2.4),
    Atom("C", 1.4, 0.0, 3.9),
]


@pytest.fixture()
def setup(workdir):
    """几何 + 假的 species 目录(每个元素一个 ``<El>.in``)。"""
    geometry = workdir / "geometry.in"
    write_aims_geometry(geometry, _ATOMS,
                        lattice=((4.07823, 0.0, 0.0), (0.0, 4.07823, 0.0),
                                 (0.0, 0.0, 25.0)))
    species = workdir / "species"
    species.mkdir()
    for element in ("Au", "S", "C"):
        (species / f"{element}.in").write_text(
            "# %s block\n  species %s\n" % (element, element), encoding="utf-8")
    return geometry, species


def _render(setup, mode: str, **kwargs) -> str:
    geometry, species = setup
    out = geometry.parent / ("control_%s.in" % mode)
    make_control(geometry, out, mode, species,
                 k_grid=(3, 3, 1) if mode.startswith("thermal") else None,
                 **kwargs)
    return out.read_text(encoding="utf-8")


def _active(text: str) -> str:
    """只留 FHI-aims 会读的行:去掉 ``#`` 注释。

    模板里的注释**故意**提到 ``calculate_forces`` / ``output aitranss`` 这些
    名字(解释为什么不写它们)。直接在全文里搜关键字会被自己的注释绊倒。
    """
    return "\n".join(
        line for line in text.splitlines()
        if line.strip() and not line.lstrip().startswith("#"))


def test_thermal_force_has_exactly_one_compute_forces(setup):
    text = _active(_render(setup, "thermal_force"))
    assert text.count("compute_forces .true.") == 1
    # 拼错的关键字会被 FHI-aims 静默忽略 —— 一个都不能出现
    assert "calculate_forces" not in text


def test_thermal_force_has_k_grid_and_no_aitranss_outputs(setup):
    text = _active(_render(setup, "thermal_force"))
    assert "k_grid 3 3 1" in text
    assert "output aitranss" not in text
    assert "output mulliken" not in text
    assert "restart aims.restart" not in text
    assert "final_forces_cleaned .true." in text
    assert "sc_accuracy_forces 5E-3" in text


def test_thermal_relax_adds_bfgs(setup):
    text = _active(_render(setup, "thermal_relax"))
    assert "relax_geometry bfgs 1.e-2" in text
    assert text.count("compute_forces .true.") == 1
    assert "k_grid 3 3 1" in text
    assert "output aitranss" not in text


def test_thermal_modes_require_k_grid(setup):
    geometry, species = setup
    with pytest.raises(ValueError):
        make_control(geometry, geometry.parent / "c.in", "thermal_force", species)


def test_kz_must_be_one(setup):
    geometry, species = setup
    with pytest.raises(ValueError):
        make_control(geometry, geometry.parent / "c.in", "thermal_force", species,
                     k_grid=(3, 3, 4))


def test_k_grid_must_be_three_positive_ints(setup):
    geometry, species = setup
    for bad in ((3, 3), (3, 3, 0), None):
        with pytest.raises(ValueError):
            make_control(geometry, geometry.parent / "c.in", "thermal_force",
                         species, k_grid=bad)


def test_existing_scf_mode_is_unchanged(setup):
    """既有 scf 口径不能因为新增 mode 而漂:aitranss/mulliken/restart 仍在。"""
    text = _render(setup, "scf")
    assert "output aitranss" in text
    assert "output mulliken" in text
    assert "restart aims.restart" in text
    assert "compute_forces" not in text
    assert "k_grid" not in text


def test_unknown_mode_still_rejected(setup):
    geometry, species = setup
    with pytest.raises(ValueError):
        make_control(geometry, geometry.parent / "c.in", "thermal", species)
