"""NEGF 内核的解析基准:每个函数都要能对上一个闭式解或大体系直接求逆。

为什么必须这样做(方案 R7):本仓做的是 FHI-aims 基的分子结声子输运,**没有任何
第三方实现可以逐数对账**。唯一可靠的锚是「同一模型的两种独立算法」——
1D 链的半无限表面格林函数有闭式二次根;有限长链可以直接对 200×200 求逆。
两边都过的代码,才有资格去接真实的 Au 电极矩阵。

模型:每层一个自由度、d0=2k、v=-k 的紧束缚链,色散 ω²=2k(1-cos q),
通带 ω∈[0, 2√k]。取 k=1 → 通带 [0,2] eV。
"""
from __future__ import annotations

import numpy as np
import pytest

from smj.phonon import landauer, negf, units

K_CHAIN = 1.0
D_CHAIN = 2.0 * K_CHAIN
V_CHAIN = -K_CHAIN
ETA = 1.0e-9


def _closed_form_surface_gf(omega: float, d0: float, v: float, eta: float) -> complex:
    """半无限链表面 GF 的闭式解:g = (z-d0 ∓ √((z-d0)²-4v²)) / 2v²,取 Im g<0。"""
    z = omega ** 2 + 1j * eta
    a = z - d0
    disc = np.sqrt(a * a - 4.0 * v * v)
    g_minus = (a - disc) / (2.0 * v * v)
    g_plus = (a + disc) / (2.0 * v * v)
    return g_minus if g_minus.imag < 0 else g_plus


def _finite_stack_surface_gf(omega: float, d0: float, v: float, eta: float,
                             layers: int = 200) -> complex:
    """同参数 M 层有限长链的表面元 —— M 足够大时与半无限极限一致。"""
    z = omega ** 2 + 1j * eta
    ham = np.zeros((layers, layers), dtype=complex)
    np.fill_diagonal(ham, d0)
    idx = np.arange(layers - 1)
    ham[idx, idx + 1] = v
    ham[idx + 1, idx] = v
    return np.linalg.inv(z * np.eye(layers) - ham)[0, 0]


# --- T4:decimation vs 直接求逆 --------------------------------------------

@pytest.mark.parametrize("omega", [2.2, 2.5, 3.0, 4.0])
def test_outside_band_matches_finite_stack(omega):
    """带外(渐逝)区:半无限 decimation 必须与 200 层直接求逆逐位一致。

    只在这里比有限链:带**内**的有限链本征值是离散的,η→0 时逆矩阵在
    本征值处出极,有限链与半无限铅不可比(这不是 bug,是物理)。
    """
    got = negf.surface_gf(omega, D_CHAIN, V_CHAIN, ETA)[0, 0]
    ref = _finite_stack_surface_gf(omega, D_CHAIN, V_CHAIN, ETA)
    assert abs(got - ref) < 1.0e-8


@pytest.mark.parametrize("omega", [0.2, 0.5, 1.0, 1.5, 1.9])
def test_inside_band_matches_closed_form(omega):
    """通带内:与闭式二次根一致(这是带内唯一正确的独立基准)。"""
    got = negf.surface_gf(omega, D_CHAIN, V_CHAIN, ETA)[0, 0]
    ref = _closed_form_surface_gf(omega, D_CHAIN, V_CHAIN, ETA)
    assert abs(got - ref) < 1.0e-10


# --- 物理不变式 -----------------------------------------------------------

def test_causality_imag_negative():
    """推迟传播子:Im g ≤ 0,全频段(含带外与 ω=0)。"""
    omega = np.linspace(0.0, 3.0, 61)
    g = negf.surface_gf(omega, D_CHAIN, V_CHAIN, ETA)
    assert g.shape == (61, 1, 1)
    assert np.all(g[:, 0, 0].imag <= 0.0)
    assert not np.isnan(g).any()


def test_batch_matches_per_point():
    """批量 ω 与逐点调用必须逐位一致 —— 生产一次 251 个频点走的是批量路径。"""
    omega = np.linspace(0.05, 2.5, 25)
    batch = negf.surface_gf(omega, D_CHAIN, V_CHAIN, ETA)[:, 0, 0]
    single = np.array([negf.surface_gf(float(w), D_CHAIN, V_CHAIN, ETA)[0, 0]
                       for w in omega])
    assert np.max(np.abs(batch - single)) < 1.0e-14


def test_two_channel_block_diagonalizes():
    """两自由度/层:当 v 与 d 对易时问题按本征通道解耦,结果必须等于逐通道闭式解。

    这条专门验证矩阵版 decimation 没有把「按标量写的公式」硬套到矩阵上。
    """
    d = np.array([[2.0, -0.5], [-0.5, 2.0]])
    v = -1.0 * np.eye(2)
    omega = 1.0
    got = negf.surface_gf(omega, d, v, ETA)

    evals, evecs = np.linalg.eigh(d)
    ref = np.zeros((2, 2), dtype=complex)
    for lam, u in zip(evals, evecs):
        g_ch = _closed_form_surface_gf(omega, float(lam), -1.0, ETA)
        ref += g_ch * np.outer(u, u.conj())
    assert np.max(np.abs(got - ref)) < 1.0e-10
    assert got.shape == (2, 2)


# --- 形状与 fail-closed ---------------------------------------------------

def test_shape_conventions():
    assert negf.surface_gf(1.0, D_CHAIN, V_CHAIN, ETA).shape == (1, 1)
    assert negf.surface_gf([1.0, 1.5], D_CHAIN, V_CHAIN, ETA).shape == (2, 1, 1)


def test_rejects_mismatched_shapes():
    with pytest.raises(ValueError):
        negf.surface_gf(1.0, np.eye(2), np.eye(3), ETA)
    with pytest.raises(ValueError):
        negf.surface_gf(1.0, np.eye(2), np.eye(2), -1.0e-9)


def test_unconverged_raises():
    """不收敛必须报错,不能静默返回一个像解的数。"""
    with pytest.raises(RuntimeError):
        negf.surface_gf(1.0, D_CHAIN, V_CHAIN, ETA, max_iter=1, tol=1.0e-30)


# --- T2:透射的解析基准 ----------------------------------------------------

CHANNEL_OMEGA = [0.2, 0.5, 1.0, 1.5, 1.9]


def test_perfect_chain_transmission_is_one():
    """完美周期链(器件层与铅层同参数)通带内 T=1,带外 T=0。

    偏差 ~1e-9 来自 η 对带边的正则化,不是误差累积。
    """
    for omega in CHANNEL_OMEGA:
        t = negf.transmission(omega, D_CHAIN, D_CHAIN, D_CHAIN, V_CHAIN, V_CHAIN, ETA)
        assert t == pytest.approx(1.0, abs=1.0e-6)
    for omega in (2.2, 2.5, 3.0):
        t = negf.transmission(omega, D_CHAIN, D_CHAIN, D_CHAIN, V_CHAIN, V_CHAIN, ETA)
        assert abs(t) < 1.0e-9


@pytest.mark.parametrize("eps_c", [2.0, 2.3, 3.0])
def test_single_site_impurity_matches_closed_form(eps_c):
    """单个「杂质层」夹在两个完美铅之间,有闭式解可对:

    G = 1/(z - ε_c - 2Σ),  Σ = V² g_lead,  Γ = -2 Im Σ,  T = |G|² Γ_L Γ_R。

    参考值只用到 ``_closed_form_surface_gf``,与 decimation 无关。
    """
    omega = 1.0
    z = omega ** 2 + 1j * ETA
    g_lead = _closed_form_surface_gf(omega, D_CHAIN, V_CHAIN, ETA)
    sigma = V_CHAIN ** 2 * g_lead
    gamma = -2.0 * sigma.imag
    g_c = 1.0 / (z - eps_c - 2.0 * sigma)
    ref = abs(g_c) ** 2 * gamma * gamma

    got = negf.transmission(omega, eps_c, D_CHAIN, D_CHAIN, V_CHAIN, V_CHAIN, ETA)
    assert got == pytest.approx(ref, rel=1.0e-9)


def test_no_rectification_when_swapping_leads():
    """线性响应下 L↔R 交换透射不变(方案 A-9)—— 本方法给不出热整流。

    两侧故意取不同的 d 与耦合,确保「对称」不是因为两边本来就一样。
    """
    for omega in (0.4, 1.0, 1.7):
        t_lr = negf.transmission(omega, 2.0, 2.0, 2.5, -1.0, -0.8, ETA)
        t_rl = negf.transmission(omega, 2.0, 2.5, 2.0, -0.8, -1.0, ETA)
        assert t_lr == pytest.approx(t_rl, abs=1.0e-12)


def test_transmission_bounded_by_channel_count():
    """单通道链:0 ≤ T ≤ 1,任何频点都不许越界(越界说明迹或 Γ 写错)。"""
    omega = np.linspace(0.02, 2.6, 40)
    t = negf.transmission(omega, D_CHAIN, D_CHAIN, D_CHAIN, V_CHAIN, V_CHAIN, ETA)
    assert t.shape == (40,)
    assert np.all(t >= -1.0e-12)
    assert np.all(t <= 1.0 + 1.0e-6)


# --- k 权重归一化(方案 D10 的契约) ---------------------------------------

def test_k_weights_must_be_normalised():
    """Σw≠1 必须报错:忘了除 N_k 会把 G_ph 放大 N_k 倍,是最隐蔽的量级错误。"""
    with pytest.raises(ValueError):
        negf.transmission(1.0, [[2.0], [3.0]], [[2.0], [2.0]], [[2.0], [2.0]],
                          -1.0, -1.0, ETA, k_weights=[1.0, 1.0])


def test_k_weighted_average_matches_manual():
    """k 平均 = 逐 k 计算后按权重求和。"""
    d_c = np.array([[[2.0]], [[3.0]]])
    d_lead = np.array([[[2.0]], [[2.0]]])
    t_k = negf.transmission(1.0, d_c, d_lead, d_lead, -1.0, -1.0, ETA,
                            k_weights=[0.5, 0.5])
    manual = 0.5 * negf.transmission(1.0, 2.0, 2.0, 2.0, -1.0, -1.0, ETA) \
        + 0.5 * negf.transmission(1.0, 3.0, 2.0, 2.0, -1.0, -1.0, ETA)
    assert t_k == pytest.approx(manual, rel=1.0e-12)


def test_single_k_is_default():
    """N_k=1 时缺省与显式 w=[1.0] 必须逐位一致(生产路径不走特殊分支)。"""
    a = negf.transmission(1.0, D_CHAIN, D_CHAIN, D_CHAIN, V_CHAIN, V_CHAIN, ETA)
    b = negf.transmission(1.0, [[D_CHAIN]], [[D_CHAIN]], [[D_CHAIN]],
                          V_CHAIN, V_CHAIN, ETA, k_weights=[1.0])
    assert a == pytest.approx(b, rel=1.0e-14)


def test_vector_matches_loop():
    omega = np.linspace(0.1, 2.3, 12)
    batch = negf.transmission(omega, D_CHAIN, D_CHAIN, D_CHAIN, V_CHAIN, V_CHAIN, ETA)
    single = np.array([
        negf.transmission(float(w), D_CHAIN, D_CHAIN, D_CHAIN, V_CHAIN, V_CHAIN, ETA)
        for w in omega])
    assert np.max(np.abs(batch - single)) < 1.0e-14


def test_transmission_rejects_bad_shapes():
    with pytest.raises(ValueError):
        # 器件 2 自由度、耦合只有 3 列
        negf.transmission(1.0, np.eye(2), [[2.0]], [[2.0]], np.ones((2, 3)),
                          -1.0, ETA)


# --- T6:端到端(transmission → thermal_conductance) -------------------------
#
# 单独测每个函数还不够:单位要靠**串联**才暴露问题(D_c 是 eV²、T_ph 无量纲、
# G 是 pW/K,任何一环漏了 EV_PER_SQRT 或 E_CHARGE 都会在这里塌掉)。所以这里
# 从动力学矩阵级别的 1D 链一路算到热导,对准两个物理锚:
#
#   1. 带宽远大于 k_BT 的单通道 = 一个量子 G_Q(283.93 pW/K @300 K);
#   2. 带顶压到 Au 声子带顶量级时,与「T_ph 取解析台阶」的同一个积分对齐。

#: 贴近 fcc-Au 最高支(19.5 meV),故意取非整数值让它落在网格点之间。
BAND_TOP_MEV = 19.37


def _chain_for_band_top(band_top_ev: float):
    """把完美链的色散 ``ω² = d0(1+cos q)`` 的通带上沿压到 ``band_top_ev``。

    通带是 ``[0, sqrt(2·d0)]``,所以 ``d0 = top²/2``、``v = -d0/2``。
    """
    d0 = band_top_ev ** 2 / 2.0
    return d0, -0.5 * d0


def test_wide_band_chain_approaches_quantum_bound():
    """带宽(2 eV)≫k_BT 时,单个完美通道必须给出一个量子 G_Q。

    这是 M1 的端到端物理锚:从 D_c(eV²) 到 pW/K 的整条链只要有一处单位错,
    结果就会差几个数量级,或者差一个 E_CHARGE 倍(1.6e19)。
    """
    omega = np.linspace(0.0, 2.5, 20001)
    t_ph = negf.transmission(omega, D_CHAIN, D_CHAIN, D_CHAIN, V_CHAIN, V_CHAIN, 1.0e-9)
    g = landauer.thermal_conductance(omega, t_ph, 300.0, band_top_mev=2100.0)
    assert g == pytest.approx(units.gq_pw_per_k(300.0), rel=2.0e-3)


def test_chain_end_to_end_matches_band_limited_reference():
    """带顶压到 Au 量级:半无限铅 + 中心区链的热导与解析台阶透射一致。

    参考值用同一公式、40 倍密的网格对 ``T_ph = Θ(top − E)`` 积分 —— 它检验的是
    「NEGF 在生产网格(251 点)上给出的透射是否足以支撑这个积分」,而不是再抄
    一遍公式。251 点网格的梯形误差约 0.2%(带边落在网格点之间)。
    """
    band_top = BAND_TOP_MEV * 1.0e-3
    d0, v = _chain_for_band_top(band_top)

    omega = np.linspace(0.0, 0.025, 251)
    t_ph = negf.transmission(omega, d0, d0, d0, v, v, 1.0e-12)
    got = landauer.thermal_conductance(omega, t_ph, 300.0)

    ref_grid = np.linspace(0.0, 0.025, 10001)
    ref_t = (ref_grid < band_top).astype(float)
    ref = landauer.thermal_conductance(ref_grid, ref_t, 300.0)

    assert got == pytest.approx(ref, rel=5.0e-3)
    # 物理量级:19.37 meV 的通道在 300 K 下只剩约 22% 的量子
    assert got / units.gq_pw_per_k(300.0) == pytest.approx(0.2237, abs=2.0e-3)
