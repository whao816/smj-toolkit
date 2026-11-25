"""软模归属（``smj.phonon.modes``）的解析基准。

虚频只报一个波数是不负责任的：−21 cm⁻¹ 可能是分子骨架的浅弯曲（无所谓），
也可能是 lead 截断造出来的假模（说明力常数错了）。这两者的区别只能从
**本征矢**里读，所以这里针对四件事做固定：

* 负曲率确实落在"最软"的那根键/那个原子上；
* 质量加权 vs 物理位移的换算不可省（否则 H 永远占头条）；
* 刚体自由度（平移/转动）能与内部形变分开；
* 区域权重归一、键长截断生效。
"""

from __future__ import annotations

import numpy as np
import pytest

from smj.phonon.modes import (
    BOND_CUTOFF_A,
    MODE_CM1_TOL,
    bond_pairs,
    describe_soft_modes,
)


def _two_atom_matrix(k: float, m1: float, m2: float) -> np.ndarray:
    """沿 x 方向的一维双原子弹簧，``k < 0`` 时给一个虚模。"""
    block = np.array([[k, -k], [-k, k]]) / np.sqrt(np.outer([m1, m2], [m1, m2]))
    matrix = np.zeros((6, 6))
    for row, col in ((0, 0), (0, 3), (3, 0), (3, 3)):
        matrix[row, col] = block[row // 3, col // 3]
    return matrix


def test_bond_pairs_respect_the_cutoff():
    positions = [(0.0, 0.0, 0.0), (2.5, 0.0, 0.0), (2.5, 3.0, 0.0)]
    elements = ["C", "H", "H"]
    pairs = bond_pairs(positions, elements)
    labels = [p["label"] for p in pairs]
    assert labels.count("C-H") == 1  # 只有 2.5 Å 的那根算键
    assert BOND_CUTOFF_A == 2.6


def test_mass_weighting_decides_which_atom_moves():
    """本征矢在质量加权坐标里;不除 sqrt(m) 会让重原子看起来动得更多。

    这里的两个原子振幅相同（质量加权），去掉 sqrt(m) 后**轻的**原子必须占优。
    """
    matrix = _two_atom_matrix(-4.0, 12.011, 1.008)
    modes = describe_soft_modes(matrix, masses=[12.011, 1.008],
                                positions=[(0, 0, 0), (1.1, 0, 0)],
                                elements=["C", "H"])
    assert modes
    top = modes[0]["top_atoms"][0]
    assert top["element"] == "H"


def test_negative_mode_lands_on_its_own_bond():
    matrix = _two_atom_matrix(-4.0, 12.011, 1.008)
    modes = describe_soft_modes(matrix, masses=[12.011, 1.008],
                                positions=[(0, 0, 0), (1.1, 0, 0)],
                                elements=["C", "H"])
    assert modes[0]["cm1"] < 0.0
    assert modes[0]["bonds"]
    assert modes[0]["bonds"][0]["label"] == "C-H"
    assert abs(modes[0]["bonds"][0]["stretch"]) == pytest.approx(1.0, abs=0.2)


def test_region_weights_are_a_partition():
    matrix = _two_atom_matrix(-4.0, 12.011, 1.008)
    modes = describe_soft_modes(matrix, masses=[12.011, 1.008],
                                positions=[(0, 0, 0), (1.1, 0, 0)],
                                elements=["C", "H"],
                                regions={"molecule": [0, 1], "lead": []})
    regions = modes[0]["regions"]
    assert regions["molecule"] == pytest.approx(1.0)
    assert regions["lead"] == pytest.approx(0.0)
    assert sum(regions.values()) == pytest.approx(1.0)


def _rank_one_direction(direction: np.ndarray, masses) -> np.ndarray:
    """做「唯一负曲率沿 direction」的矩阵。

    注意方向必须先在**质量加权坐标**里给定（``eigh`` 的本征矢活在这一套坐标），
    否则 ``describe_soft_modes`` 去掉 ``√m`` 之后拿到的已经不是原方向——这正是
    "刚体占比"要把质量因子算清楚的原因。
    """
    weighted = np.asarray(direction, dtype=float) * np.repeat(np.sqrt(masses), 3)
    unit = weighted / np.linalg.norm(weighted)
    return -4.0 * np.outer(unit, unit)


MASSES = [12.011, 1.008]


def test_rigid_translation_is_flagged_as_rigid_body():
    """整体平移的刚体占比应为 1:这种模不是"分子在形变"。"""
    direction = np.zeros(6)
    direction[0] = 1.0
    direction[3] = 1.0
    matrix = _rank_one_direction(direction, MASSES)
    modes = describe_soft_modes(matrix, masses=MASSES,
                                positions=[(0, 0, 0), (1.1, 0, 0)],
                                elements=["C", "H"])
    assert modes[0]["rigid_body"]["fraction"] == pytest.approx(1.0, abs=1e-6)


def test_internal_deformation_is_not_rigid():
    direction = np.zeros(6)
    direction[0] = 1.0
    direction[3] = -1.0
    matrix = _rank_one_direction(direction, MASSES)
    modes = describe_soft_modes(matrix, masses=MASSES,
                                positions=[(0, 0, 0), (1.1, 0, 0)],
                                elements=["C", "H"])
    assert modes[0]["rigid_body"]["fraction"] < 0.05


def test_modes_are_sorted_and_capped():
    matrix = np.diag([-1.0, -4.0, -9.0, 0.0, 0.0, 0.0])
    modes = describe_soft_modes(matrix, masses=[1.0, 1.0],
                                positions=[(0, 0, 0), (1.0, 0, 0)],
                                elements=["C", "C"], max_modes=2)
    assert len(modes) == 2
    assert modes[0]["cm1"] < modes[1]["cm1"]


def test_positive_tolerance_is_rejected():
    with pytest.raises(ValueError):
        describe_soft_modes(np.eye(6), cm1_tol=1.0)
    assert MODE_CM1_TOL == -10.0


def test_matrix_shape_is_validated():
    with pytest.raises(ValueError):
        describe_soft_modes(np.zeros((5, 5)))   # 边长不是 3 的倍数
    with pytest.raises(ValueError):
        describe_soft_modes(np.zeros((6, 3)))


def test_without_structure_only_frequencies():
    """没给元素/坐标时只报频率与参与比,不得抛错（旧 thermal.json 兼容）。"""
    matrix = _two_atom_matrix(-4.0, 12.011, 1.008)
    modes = describe_soft_modes(matrix)
    assert modes[0].get("bonds", []) == []
    assert modes[0].get("regions", {}) == {}
    assert all(atom["element"] is None for atom in modes[0]["top_atoms"])
