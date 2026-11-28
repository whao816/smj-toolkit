"""动力学矩阵与虚频口径的解析基准。

模型:两个原子由一根沿 x 的弹簧相连(C–H 级别参数),解析结果完全已知 ——
6 个模里 5 个是 0(3 个平移 + 2 个无恢复力的 y/z 方向),1 个光学模
``ħω = EV_PER_SQRT·sqrt(k(1/m0+1/m1))``。

这组测试专门覆盖调研报告 §6 的第一号陷阱:把 ``sqrt(λ)`` 当 eV(漏掉
``EV_PER_SQRT``)。漏了会让 C–H 伸缩从 2960 cm⁻¹ 变成 19 cm⁻¹ 这种荒谬值。
"""
from __future__ import annotations

import numpy as np
import pytest

from smj.phonon import dynamical, units

K_SPRING = 30.0            # eV/Å²
M_C = 12.011               # amu
M_H = 1.008


def _spring_matrix(k: float = K_SPRING) -> np.ndarray:
    """两原子沿 x 的弹簧:K[0,0]=K[3,3]=k,K[0,3]=K[3,0]=-k,其余为零。"""
    mat = np.zeros((6, 6))
    mat[0, 0] = k
    mat[3, 3] = k
    mat[0, 3] = -k
    mat[3, 0] = -k
    return mat


def _hand_block(k: float = K_SPRING) -> np.ndarray:
    """x 方向 2×2 质量加权重块 eV/(Å²·amu) —— 与实现独立地手写一遍。"""
    return np.array([
        [k / M_C, -k / np.sqrt(M_C * M_H)],
        [-k / np.sqrt(M_C * M_H), k / M_H],
    ])


def test_two_atom_spring_eigenvalues():
    """5 个零模 + 1 个光学模,光学模与手写 2×2 的本征值一致。"""
    d = dynamical.dynamical_matrix(_spring_matrix(), [M_C, M_H])
    lam = np.linalg.eigvalsh(d)

    hand = np.linalg.eigvalsh(_hand_block()) * units.EV_PER_SQRT ** 2
    optical = np.sort(lam)[-1]
    assert optical == pytest.approx(np.sort(hand)[-1], rel=1.0e-12)
    # 5 个零模(平移 3 + 无恢复力 2)
    assert np.allclose(np.sort(lam)[:5], 0.0, atol=1.0e-12)


def test_physical_frequency_of_ch_stretch():
    """C–H 级别弹簧必须给出 2800–3100 cm⁻¹ 的伸缩频率(物理地标)。"""
    d = dynamical.dynamical_matrix(_spring_matrix(), [M_C, M_H])
    freqs = dynamical.mode_frequencies_cm1(d)
    assert 2800.0 < freqs[-1] < 3100.0
    # 期望值的解析式
    expect = units.ev_to_cm1(
        units.EV_PER_SQRT * np.sqrt(K_SPRING * (1.0 / M_C + 1.0 / M_H)))
    assert freqs[-1] == pytest.approx(expect, rel=1.0e-12)


def test_units_off_gives_eva2_amu():
    """scale_to_ev2=False 时本征值回到 eV/(Å²·amu),比值恰为 EV_PER_SQRT²。"""
    raw = dynamical.dynamical_matrix(_spring_matrix(), [M_C, M_H],
                                     scale_to_ev2=False)
    ev2 = dynamical.dynamical_matrix(_spring_matrix(), [M_C, M_H])
    assert np.allclose(ev2, raw * units.EV_PER_SQRT ** 2, rtol=1e-14)
    lam_raw = np.sort(np.linalg.eigvalsh(raw))[-1]
    assert lam_raw == pytest.approx(
        _hand_block().trace() - np.linalg.eigvalsh(_hand_block())[0], rel=1.0e-12)


def test_symmetrize_removes_machine_asymmetry():
    k = _spring_matrix()
    k[1, 2] = 1.0e-12          # 伪造一点非对称
    d = dynamical.dynamical_matrix(k, [M_C, M_H])
    assert np.allclose(d, d.T, atol=0.0, rtol=0.0)
    d_nosym = dynamical.dynamical_matrix(k, [M_C, M_H], symmetrize=False)
    # atol/rtol 必须为 0:1e-12 的非对称会被 allclose 的默认 atol=1e-8 吞掉
    assert not np.allclose(d_nosym, d_nosym.T, atol=0.0, rtol=0.0)


# --- 虚频口径 -------------------------------------------------------------

def _diag_ev2(cm1_values) -> np.ndarray:
    """按 cm⁻¹ 直接构造一个对角 D̃(负值 = 虚频)。"""
    lam = np.array([np.sign(v) * (abs(v) / units.EV_TO_CM1) ** 2 for v in cm1_values])
    return np.diag(lam)


def test_imaginary_modes_tolerance():
    """-20 cm⁻¹ 标记、-5 cm⁻¹ 放过 —— 容差默认 -10 cm⁻¹。"""
    d = _diag_ev2([-20.0, -5.0, 120.0, 400.0])
    info = dynamical.imaginary_modes(d)
    assert info["count"] == 1
    assert info["worst_cm1"] == pytest.approx(-20.0, rel=1.0e-9)
    assert info["modes"][0]["cm1"] == pytest.approx(-20.0, rel=1.0e-9)
    assert info["n_modes"] == 4
    assert info["tolerance_cm1"] == -10.0


def test_imaginary_modes_none_for_stable_system():
    d = dynamical.dynamical_matrix(_spring_matrix(), [M_C, M_H])
    info = dynamical.imaginary_modes(d)
    assert info["count"] == 0
    assert info["worst_cm1"] == 0.0


def test_imaginary_modes_sorted_and_negative():
    """多个虚频按从负到正排序,``worst_cm1`` 是最负的那个。"""
    d = _diag_ev2([-30.0, -50.0, -15.0, 10.0])
    info = dynamical.imaginary_modes(d)
    got = [m["cm1"] for m in info["modes"]]
    assert got == pytest.approx([-50.0, -30.0, -15.0], rel=1.0e-9)
    assert info["worst_cm1"] == pytest.approx(-50.0, rel=1.0e-9)


def test_imaginary_modes_rejects_positive_tolerance():
    """cm1_tol 必须是负阈值,传正数多半是把符号约定搞反了。"""
    with pytest.raises(ValueError):
        dynamical.imaginary_modes(np.eye(2), cm1_tol=10.0)


# --- ASR ------------------------------------------------------------------

def test_acoustic_sum_residual_zero_for_translation_invariant():
    assert dynamical.acoustic_sum_residual(_spring_matrix()) == pytest.approx(0.0)


def test_acoustic_sum_residual_reports_row_offset():
    """给某个自由度加一个未补偿的力 → 残差恰等于它。"""
    k = _spring_matrix()
    k[4, 1] += 0.05                      # y 方向多出一根「外部弹簧」
    assert dynamical.acoustic_sum_residual(k) == pytest.approx(0.05, rel=1.0e-12)
    rel = dynamical.acoustic_sum_relative_residual(k)
    assert rel == pytest.approx(0.05 / 30.0, rel=1.0e-12)


def test_acoustic_residual_is_computed_on_k_not_on_d():
    """残差在 K 上算:质量加权后每行归一化不同,不能在 D 上做平移检查。"""
    masses = [M_C, M_H]
    k_ok = _spring_matrix()
    d_ok = dynamical.dynamical_matrix(k_ok, masses)
    # D̃·1 并不为零(各原子分量的行和按 1/sqrt(m_i m_j) 缩放),这不是 ASR 违反
    assert np.max(np.abs(d_ok.sum(axis=1))) > 1.0e-3
    assert dynamical.acoustic_sum_residual(k_ok) == pytest.approx(0.0)


# --- 输入校验 -------------------------------------------------------------

def test_validation():
    with pytest.raises(ValueError):
        dynamical.dynamical_matrix(np.eye(5), [M_C, M_H])       # 维度不符
    with pytest.raises(ValueError):
        dynamical.dynamical_matrix(_spring_matrix(), [M_C, 0.0])   # 质量非正
    with pytest.raises(ValueError):
        k = _spring_matrix()
        k[0, 1] = np.nan
        dynamical.dynamical_matrix(k, [M_C, M_H])
    with pytest.raises(ValueError):
        dynamical.mode_frequencies_cm1(np.zeros((2, 3)))
