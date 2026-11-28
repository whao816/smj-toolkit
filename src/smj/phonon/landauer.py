"""声子 Landauer 线性响应:``T_ph(ω)`` → ``G_ph(T)``(唯一公式源)。

**公式**(方案 §3.3,与调研报告 §1 等价,报告的前因子本身是对的)::

    G_ph(T) = (ħ² / 2π k_B T²) ∫₀^∞ ω² T_ph(ω) e^{ħω/k_BT}/(e^{ħω/k_BT}-1)² dω

把 ``ω`` 换成能量 ``E = ħω``(我们的网格单位),化成实现里用的形式::

    G_ph(T) = 1/(h k_B T²) ∫₀^∞ E² T_ph(E) K(E/k_BT) dE ,
    K(x) = e^x/(e^x-1)² = 1/(4 sinh²(x/2))

两者差一个 ``h = 2πħ``,可互换;本模块只实现后一种(避免每处都除 2π)。

**两个必须成立的解析结论**(测试钉住):

- ``T_ph ≡ n`` ⇒ ``G = n·G_Q(T)``,``G_Q = π²k_B²T/(3h)``(300 K 283.93 pW/K);
- 单通道结 ``T_ph ≡ 1`` 在最理想情形下取到 ``G_Q`` —— 这是硬上界,超过就是错的。

**为什么不用 ``(k_B/h)∫T_ph x²…dx``**:那个写法量纲是 J/K 而不是 W/K,少一个
``k_BT``。调研报告 §1 给的是 ``ħ²/2πk_BT²`` 形式(对的),不要「简化」成前者。
"""
from __future__ import annotations

import math

import numpy as np

from . import units

__all__ = ["bose_kernel", "thermal_conductance", "cumulative_thermal_conductance"]

DEFAULT_MAX_CHANNELS = 3
NEGATIVE_TRANSMISSION_TOL = 1.0e-12
MAX_ABOVE_BAND_T = 1.0e-6
ACOUSTIC_ENERGY_EV = 0.5e-3
ACOUSTIC_TOL = 1.0e-6


def above_band(energy_mev, band_top_mev):
    return np.asarray(energy_mev) > band_top_mev + 1.0e-9


def conductance_bound(temperature_k, max_channels):
    """Upper bound in pW/K when at most max_channels modes transmit at any energy."""
    if not math.isfinite(max_channels) or max_channels <= 0:
        raise ValueError("max_channels 必须为正有限值")
    return max_channels * units.gq_pw_per_k(temperature_k)


def _trapz(y, x):
    """复合梯形积分(自带,不用 ``np.trapz``)。

    numpy 2.0 把 ``np.trapz`` 改名成 ``np.trapezoid`` 并删掉旧名;本仓要同时支持
    1.24+ 与 2.x,**不能**依赖任一版本的存在性。电子链的
    ``smj.junction.landauer`` 出于同样的理由自带了一份,这里保持同一口径。
    """
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def bose_kernel(x):
    """``x²·e^x/(e^x-1)²``,在 ``x=0`` 处按极限取 1。

    直接用 ``x²eˣ/(eˣ-1)²`` 在 ``x>700`` 时溢出成 inf/inf=nan;这里改成
    ``((x/2)/sinh(x/2))²``,``x→0`` 自然趋于 1、``x→∞`` 按 ``x²e^{-x}`` 衰减,
    全程无奇点。
    """
    x = np.asarray(x, dtype=float)
    if np.any(x < 0):
        raise ValueError("x=ħω/k_BT 必须非负;负值通常意味着温度或能量维度搞反了")
    half = 0.5 * x
    out = np.ones_like(x)
    big = half > 1.0e-8
    hb = half[big]
    res = np.empty_like(hb)
    # sinh(half) 在 half≳710 溢出成 inf(→ 0/0 或 RuntimeWarning),故 300 以上
    # 直接用渐近式 ((2h)e^{-h}/2)² = 4h²e^{-2h},相对误差 O(e^{-600})。
    small = hb <= 300.0
    res[small] = (hb[small] / np.sinh(hb[small])) ** 2
    tail = hb[~small]
    res[~small] = 4.0 * tail ** 2 * np.exp(-2.0 * tail)
    out[big] = res
    return out


def validate_spectrum(omega_ev, t_ph):
    """Validate the energy axis and clamp negative transmission within roundoff."""
    e = np.asarray(omega_ev, dtype=float)
    t = np.asarray(t_ph, dtype=float)
    if e.ndim != 1 or t.ndim != 1:
        raise ValueError("omega_ev 与 t_ph 都必须是一维数组")
    if e.size != t.size:
        raise ValueError(
            "omega_ev(%d 个点)与 t_ph(%d 个点)长度不一致" % (e.size, t.size))
    if e.size < 2:
        raise ValueError("至少需要 2 个频点才能积分")
    if not np.all(np.isfinite(e)) or not np.all(np.isfinite(t)):
        raise ValueError("频率/透射含非有限值;损坏的谱不能当物理结果")
    if np.any(np.diff(e) <= 0):
        raise ValueError("omega_ev 必须严格递增")
    if e[0] < 0:
        raise ValueError("omega_ev 不能为负")
    if np.any(t < -NEGATIVE_TRANSMISSION_TOL):
        raise ValueError("T_ph 含显著负透射（最小值 %.3g）" % float(np.min(t)))
    return e, np.maximum(t, 0.0)


def thermal_conductance(
    omega_ev,
    t_ph,
    temperature_k: float,
    *,
    band_top_mev: float = units.AU_BAND_TOP_MEV,
    max_acoustic: float = 3.0,
    max_channels: float = DEFAULT_MAX_CHANNELS,
):
    """积分得 ``G_ph(T)``,单位 **pW/K**。

    Parameters
    ----------
    omega_ev : (N,) array
        声子能量 ħω,eV,**严格递增**,从 0(或很小的正值)开始。
    t_ph : (N,) array
        声子透射 ``T_ph(E)``,无量纲。``E > band_top_mev`` 处必须为 0。
    temperature_k : float
        温度 K。与电子链共用 ``transport.temperature``(方案 D8)。
    band_top_mev : float
        电极声子带顶(默认 fcc-Au 19.5 meV)。**积分上限不能取分子的最高振动
        频率**:分子内 100–400 meV 的模在 Au 里没有可透射通道,把它们算进去
        等于给弹性透射编造出不存在的通道(调研报告 §6 陷阱⑨)。
    max_acoustic : float
        声学支数上限，用于检查低频极限。
    max_channels : float
        全能谱上的传输通道数上限。默认三通道协议取 3。含额外通道的
        多模电极应按其色散及横向归一化提供适当上限。

    Returns
    -------
    float
        ``G_ph``,pW/K。

    Raises
    ------
    ValueError
        网格/数值不合法,或三条物理上界之一被突破(见 Notes)。

    Notes
    -----
    三条 fail-closed 检查都对应「最贵的错」:

    1. 带顶之上非零 → 电极里没有那条支,透射是算错的(常见于自能构造错误);
    2. ``T_ph(0) > max_acoustic`` → 声学支计数错(常见于 k 求和未除 ``N_k``);
    3. ``G > max_channels·G_Q`` → 超过给定传输通道数可承载的热导。
    """
    if not math.isfinite(temperature_k) or temperature_k <= 0:
        raise ValueError("温度必须为正:%r K" % (temperature_k,))
    e, t = validate_spectrum(omega_ev, t_ph)

    # 1) 带顶之上必须无通道
    above = above_band(e * 1.0e3, band_top_mev)
    if np.any(above) and float(np.max(np.abs(t[above]))) > MAX_ABOVE_BAND_T:
        worst = float(e[above][int(np.argmax(np.abs(t[above])))] * 1e3)
        raise ValueError(
            "T_ph 在电极声子带顶(%.2f meV)之上仍非零(最大在 %.2f meV);"
            "弹性声子透射没有高于 lead 最高支的通道,积分上限与自能都要查"
            % (band_top_mev, worst))

    # 2) 声学极限
    if e[0] <= ACOUSTIC_ENERGY_EV and t[0] > max_acoustic + ACOUSTIC_TOL:
        raise ValueError(
            "T_ph(ω→0)=%.4f 超过声学支数 %.1f;常见原因是横向 k 求和未除 N_k"
            % (t[0], max_acoustic))

    kt = units.KB_EV_PER_K * temperature_k
    integrand = kt ** 2 * bose_kernel(e / kt) * t
    integral = _trapz(integrand, e)

    # G = 1/(h k_B T²) ∫ E² T_ph K dE;用 (kT)²·w 代替 E²K 后等价于对
    # (kT)²wT dE 积分,故预因子里的 T² 用 kt·T 抵消。
    g_w_per_k = integral / (units.H_EV_S * kt * temperature_k) * units.E_CHARGE
    g_pw_per_k = g_w_per_k * 1.0e12

    # 3) 量子上限
    bound = conductance_bound(temperature_k, max_channels)
    if g_pw_per_k > bound:
        raise ValueError(
            "G_ph=%.3f pW/K 超过给定通道数的上限 %g×G_Q(=%.1f pW/K)"
            % (g_pw_per_k, max_channels, bound))

    return g_pw_per_k


def cumulative_thermal_conductance(omega_ev, t_ph, temperature_k: float):
    """累积热导 ``G(ω) = G_ph({E ≤ ω})``,pW/K,与 :func:`thermal_conductance` 同口径。

    用于出图(方案 §6.6「累积热导」):把 ``G_ph`` 拆成「哪些频率贡献了多少」,
    这是判断结果是否被界面低频模主导的主要证据。
    """
    if not math.isfinite(temperature_k) or temperature_k <= 0:
        raise ValueError("温度必须为正:%r K" % (temperature_k,))
    e, t = validate_spectrum(omega_ev, t_ph)
    kt = units.KB_EV_PER_K * temperature_k
    integrand = kt ** 2 * bose_kernel(e / kt) * t
    # 从 0 开始累积:trapz 在每段的贡献 = 0.5*(f_i+f_{i+1})*ΔE
    seg = 0.5 * (integrand[1:] + integrand[:-1]) * np.diff(e)
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    scale = units.E_CHARGE / (units.H_EV_S * kt * temperature_k) * 1.0e12
    return cum * scale
