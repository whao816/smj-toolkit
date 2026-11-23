"""声子 Landauer 积分的解析基准:常数透射必须精确回到 G_Q 的整数倍。

这一组是 M1 的核心:``T_ph ≡ n`` ⇒ ``G = n·G_Q`` 把「前因子、ħ、k_B、单位制」
一次全部钉死。任何一处漏了 ``2π`` 或 ``k_BT``,这里的偏差都是数量级的,不会
被网格误差掩盖。调研报告 §1 的 ``ħ²/2πk_BT²`` 形式是对的,曾被怀疑的是它的
``(k_B/h)`` 简化写法(量纲少一个 k_BT)。
"""
from __future__ import annotations

import numpy as np
import pytest

from smj.phonon import landauer, units

TEMPERATURES = [4.0, 100.0, 300.0]


def _grid(temperature: float, x_max: float = 30.0, points: int = 2001):
    """覆盖到 x=ħω/k_BT=30 的网格 —— 更远的尾巴对积分贡献 <1e-10。"""
    kt = units.KB_EV_PER_K * temperature
    return np.linspace(0.0, x_max * kt, points)


@pytest.mark.parametrize("temperature", TEMPERATURES)
def test_n_channels_give_n_quantum(temperature):
    """``T_ph ≡ 3`` 必须是 ``3·G_Q(T)`` —— 精确到网格截断(1e-9 以内)。"""
    e = _grid(temperature)
    got = landauer.thermal_conductance(
        e, np.full_like(e, 3.0), temperature, band_top_mev=1.0e9)
    want = 3.0 * units.gq_pw_per_k(temperature)
    assert got == pytest.approx(want, rel=1.0e-9)


@pytest.mark.parametrize("temperature", TEMPERATURES)
def test_single_channel_is_the_quantum_bound(temperature):
    """单通道:``T_ph ≡ 1`` 恰好取到 G_Q —— 这就是本方法的天花板。"""
    e = _grid(temperature)
    got = landauer.thermal_conductance(
        e, np.full_like(e, 1.0), temperature, band_top_mev=1.0e9)
    assert got == pytest.approx(units.gq_pw_per_k(temperature), rel=1.0e-9)


def test_quantum_bound_scales_linearly_with_temperature():
    """常数透射下 G ∝ T(玻色因子已积掉)—— 用来验证没有把 T 多算或少算一次。"""
    g300 = landauer.thermal_conductance(
        _grid(300.0), np.full(2001, 2.0), 300.0, band_top_mev=1.0e9)
    g600 = landauer.thermal_conductance(
        _grid(600.0), np.full(2001, 2.0), 600.0, band_top_mev=1.0e9)
    assert g600 / g300 == pytest.approx(2.0, rel=1.0e-9)


def test_bose_kernel_identity():
    """∫₀^∞ x²eˣ/(eˣ-1)²dx = π²/3,且 x=0 处极限为 1(不能是 0 或 inf)。"""
    x = np.linspace(0.0, 40.0, 40001)
    assert landauer._trapz(landauer.bose_kernel(x), x) == pytest.approx(
        np.pi ** 2 / 3.0, rel=1.0e-8)
    assert landauer.bose_kernel(np.array([0.0]))[0] == 1.0
    # 大 x 处按 x²e^{-x} 衰减;sinh 写法不能溢出成 nan/inf
    assert np.all(np.isfinite(landauer.bose_kernel(np.array([700.0, 2000.0]))))
    assert landauer.bose_kernel(np.array([40.0]))[0] < 1.0e-13


def test_zero_transmission_gives_zero():
    e = _grid(300.0)
    assert landauer.thermal_conductance(e, np.zeros_like(e), 300.0) == 0.0


# --- 三条 fail-closed -----------------------------------------------------

def test_rejects_transmission_above_gold_band_top():
    """19.5 meV 之上不许有通道:超过 Au 最高支的弹性透射物理上不存在。"""
    e = np.linspace(0.0, 0.030, 301)          # 到 30 meV,越过带顶
    t = np.zeros_like(e)
    t[:] = np.where(e <= 0.0195, 1.0, 0.3)
    with pytest.raises(ValueError, match="带顶"):
        landauer.thermal_conductance(e, t, 300.0)
    # 把带顶调高就能通过 —— 证明拒绝的是「带顶之上非零」而不是别的
    landauer.thermal_conductance(e, t, 300.0, band_top_mev=1.0e9)


def test_rejects_more_acoustic_channels_than_exist():
    """T_ph(0) 超过声学支数:最常见原因是 k 求和没有除 N_k。"""
    e = _grid(300.0, x_max=6.0, points=301)
    t = np.full_like(e, 5.0)
    with pytest.raises(ValueError, match="声学支"):
        landauer.thermal_conductance(e, t, 300.0, band_top_mev=1.0e9)


def test_rejects_conductance_above_quantum_bound():
    """T_ph ≡ 10(但 t(0)=0 以绕过声学检查)必须被量子上限拦住。"""
    e = _grid(300.0, x_max=30.0)
    t = np.full_like(e, 10.0)
    t[0] = 0.0
    with pytest.raises(ValueError, match="G_Q"):
        landauer.thermal_conductance(e, t, 300.0, band_top_mev=1.0e9)


# --- 累积热导 -------------------------------------------------------------

def test_cumulative_equals_total_at_the_top():
    e = _grid(300.0, x_max=30.0, points=301)
    t = np.full_like(e, 2.0)
    cum = landauer.cumulative_thermal_conductance(e, t, 300.0)
    total = landauer.thermal_conductance(e, t, 300.0, band_top_mev=1.0e9)
    assert cum.shape == e.shape
    assert cum[0] == 0.0
    assert cum[-1] == pytest.approx(total, rel=1.0e-12)
    assert np.all(np.diff(cum) >= 0.0)


# --- 输入校验 -------------------------------------------------------------

def test_grid_validation():
    e = np.array([0.0, 0.01, 0.02])
    t = np.array([1.0, 1.0, 1.0])
    with pytest.raises(ValueError):
        landauer.thermal_conductance(e[::-1], t, 300.0)      # 未递增
    with pytest.raises(ValueError):
        landauer.thermal_conductance(np.array([-0.01, 0.01, 0.02]), t, 300.0)
    with pytest.raises(ValueError):
        landauer.thermal_conductance(e, np.array([1.0, 1.0]), 300.0)   # 长度不符
    with pytest.raises(ValueError):
        landauer.thermal_conductance(e, np.array([1.0, np.nan, 1.0]), 300.0)
    with pytest.raises(ValueError):
        landauer.thermal_conductance(e, t, 0.0)              # 温度必须为正
