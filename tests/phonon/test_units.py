"""换算常量的独立验证:用物理地标而不是「和文档对一遍」钉死。

为什么这组测试值得单独存在(方案 A-1):调研报告的换算系数错了 3.9 倍,而错法
非常隐蔽 —— 它同时给了一组「看似自洽」的示例(λ=0.010 → 25.4 meV 「≈Au 带顶」),
只有拿真实的 Au/C-C 地标去对,才能发现它把带顶抬到了 57.6 meV、把 C-C 伸缩抬到
X-H 区。所以这里一律用**独立可查的物理量**作锚:

- Au 带顶 19.5 meV(4.71 THz,中子散射实测)
- Au θ_D=165 K → 14.2 meV(德拜平均,不是带顶)
- C-C 伸缩 k=30–50 eV/Å²、m=12.011 → 824–1064 cm⁻¹
- 单通道量子热导 300 K = 283.93 pW/K
"""
from __future__ import annotations

import math

import pytest

from smj.phonon import units


def test_codata_selfconsistent():
    """两套口径(ev / SI)必须是同一组 CODATA,不能各写一个数。"""
    assert units.KB_EV_PER_K * units.E_CHARGE == pytest.approx(units.KB_J_PER_K, rel=1e-15)
    assert units.H_EV_S * units.E_CHARGE == pytest.approx(units.H_J_S, rel=1e-15)
    assert units.HBAR_J_S == pytest.approx(units.H_J_S / (2.0 * math.pi), rel=1e-9)


def test_matches_phonopy_constant():
    """与 phonopy 的 VaspToTHz 交叉核对 —— 两条独立实现必须一致。

    phonopy 的 15.633302 是它自己从 CODATA 推的,我们推的是 15.633304。相对差
    1.4e-7 来自两边取的 CODATA 版本尾数不同,不影响任何物理结论。
    """
    rel = units.THZ_PER_SQRT / units.PHONOPY_VASP_TO_THZ - 1.0
    assert abs(rel) < 1e-6
    assert abs(rel) > 0.0  # 不是抄的同一个数


def test_lambda_to_mev_landmark():
    """λ=0.1 是 Au 带顶量级:20.445 meV(报告给的 40 meV 是 0.2543 的后果)。"""
    assert units.lambda_to_ev(0.1) * 1e3 == pytest.approx(20.4454, abs=2e-3)


def test_gold_band_top_lambda():
    """Au 带顶 19.5 meV 对应 λ≈0.0910;θ_D 的 14.2 meV 对应 λ≈0.0482。"""
    for e_mev, lam in ((19.5, 0.090965), (14.2, 0.048237)):
        got = (e_mev * 1e-3 / units.EV_PER_SQRT) ** 2
        assert got == pytest.approx(lam, rel=2e-4)


def test_cc_stretch_landmark():
    """萘的 C-C 伸缩必须落在 800–1100 cm⁻¹,不能落进 X-H 区(>2800)。"""
    for k_ev_a2, expected_cm1 in ((30.0, 824.14), (50.0, 1063.96)):
        lam = k_ev_a2 / 12.011
        cm1 = units.ev_to_cm1(units.lambda_to_ev(lam))
        assert cm1 == pytest.approx(expected_cm1, abs=0.05)
        assert 700.0 < cm1 < 1200.0


def test_report_coefficient_is_rejected():
    """回归:0.2543 必须被判为错 —— 它会把 Au 带顶抬到 57.6 meV。"""
    wrong = 0.2543
    assert wrong / units.EV_PER_SQRT == pytest.approx(3.9332, abs=1e-3)
    au_top_with_wrong = wrong * math.sqrt(0.090965) * 1e3
    assert au_top_with_wrong > 70.0  # 报告口径给出 76.7 meV,物理上不可能


def test_quantum_conductance():
    """G_Q(T)=π²k_B²T/(3h):300 K 283.93 pW/K,且严格正比于 T。"""
    assert units.gq_pw_per_k(300.0) == pytest.approx(283.929345, rel=1e-6)
    assert units.gq_pw_per_k(100.0) == pytest.approx(94.643115, rel=1e-6)
    assert units.gq_pw_per_k(4.0) == pytest.approx(3.785725, rel=1e-6)
    for T in (1.0, 77.0, 300.0, 1000.0):
        assert units.gq_pw_per_k(2.0 * T) == pytest.approx(2.0 * units.gq_pw_per_k(T), rel=1e-12)


def test_ev_to_cm1():
    assert units.ev_to_cm1(0.1) == pytest.approx(806.5544, abs=1e-3)
    # 1 eV = 8065.544 cm⁻¹;不要误用 8065.73(那是把 c 取成 2.9979e10 的旧尾数)
    assert abs(units.EV_TO_CM1 - 8065.73) > 0.1


def test_invalid_inputs_raise():
    with pytest.raises(ValueError):
        units.lambda_to_ev(-0.01)          # 虚频不走这个函数
    with pytest.raises(ValueError):
        units.gq_pw_per_k(0.0)
