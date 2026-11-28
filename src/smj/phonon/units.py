"""声子输运的唯一换算源:力常数 → 频率、能量 → 波数、单通道量子热导。

**为什么单独建一个模块(方案 A-1)**:调研报告 §3.4 给的是
``ħω [eV] = 0.2543 * sqrt(λ [eV/(Å²·amu)])``,这个系数大了 **3.9332 倍**
—— ``0.2543² = 0.0646685``,等于把正确值又开了一次方。用它算出的 Au 声子带顶
是 57.6 meV(≈670 K),物理上不可能;萘的 C-C 伸缩会落到 3242–4189 cm⁻¹(那是
X-H 区)。正确系数 ``0.0646542`` 由下面的 CODATA 常数直接推出,并用两个独立的
物理地标钉死(``tests/phonon/test_units.py``):

- **Au**:θ_D=165 K → 14.2 meV;λ=0.0514 → 14.65 meV。实测带顶 19.5 meV
  (4.71 THz, Lynn/Smith/Nicklow PRB 8, 3493 (1973))对应 λ=0.0910。
- **C-C**:k=30–50 eV/Å²、m=12.011 → λ=2.50–4.16 → 824–1064 cm⁻¹。

**口径说明**:本模块一律用 CODATA 2018。``smj.junction.landauer`` 里的
``KB_EV_PER_K``/``H_EV_S`` 是更早的一组值,被电子链的冻结锚点
(``log10_G_over_G0_like``/``k_e_pW_per_K``)绑住,**不要改动,也不要在这里
import 它们**。两组值相对差 <1e-8,对热导的影响远低于任何数值误差;但混用会让
"同一个 k_B 有两个值"变成事故来源。声子侧的单通道量子热导必须用 **SI 的 k_B/h**
再换 pW —— 用 eV 版会差一个 ``E_CHARGE``。
"""
from __future__ import annotations

import math

# --- CODATA 2018 基本常数 -------------------------------------------------
E_CHARGE = 1.602176634e-19          # C
KB_J_PER_K = 1.380649e-23           # J/K   (2019 SI 定义值,精确)
H_J_S = 6.62607015e-34              # J·s   (2019 SI 定义值,精确)
HBAR_J_S = 1.054571817e-34          # J·s
AMU_KG = 1.66053906660e-27          # kg
ANGSTROM_M = 1.0e-10                # m
SPEED_OF_LIGHT_M_S = 2.99792458e8   # m/s

# eV 口径(仅用于能量→温度/波数的换算,不用于量子热导)
KB_EV_PER_K = KB_J_PER_K / E_CHARGE         # 8.617333262e-5 eV/K
H_EV_S = H_J_S / E_CHARGE                   # 4.135667696e-15 eV·s
HBAR_EV_S = HBAR_J_S / E_CHARGE             # 6.582119569e-16 eV·s

# --- 力常数 → 频率 --------------------------------------------------------
# λ = K_ij / sqrt(m_i m_j),单位 eV/(Å²·amu)。角频率
#   ω [rad/s] = sqrt(1 eV/(Å²·amu)) * sqrt(λ)
#             = sqrt(E_CHARGE / (AMU_KG * ANGSTROM_M²)) * sqrt(λ)
# 再乘 ħ 折成 eV。两个因子都是常数,先合并成一个系数,避免每处各算一遍。
_OMEGA_PER_SQRT_RAD_S = math.sqrt(E_CHARGE / (AMU_KG * ANGSTROM_M ** 2))
EV_PER_SQRT = HBAR_J_S * _OMEGA_PER_SQRT_RAD_S / E_CHARGE      # 0.06465421... eV

# ω → ν 只差 2π:ν [THz] = ħω [eV] / h [eV·s] / 1e12
THZ_PER_SQRT = EV_PER_SQRT / H_EV_S / 1.0e12                   # 15.633306... THz

# 1 eV 对应多少 cm⁻¹:1/hc,eV 与 cm⁻¹ 都是能量口径(与 amu 无关)
EV_TO_CM1 = E_CHARGE / (H_J_S * SPEED_OF_LIGHT_M_S * 100.0)    # 8065.544... cm⁻¹

#: phonopy/ASE 口径的等价常量,用于交叉核对(phonopy.VaspToTHz)
#: 单位 THz·sqrt(eV/Å²/amu);与本模块的 ``THZ_PER_SQRT`` 必须一致到 1e-6。
PHONOPY_VASP_TO_THZ = 15.633302

#: fcc-Au 声子带顶(实测最高支 L 点 4.71 THz = 19.5 meV)。
#: 透射在此之上必须严格为零 —— 弹性声子透射没有大于电极最高支的通道。
AU_BAND_TOP_MEV = 19.5


def lambda_to_ev(lam_ev_per_a2_amu: float) -> float:
    """质量加权力常数 λ → ħω [eV]。

    λ = K_ij / sqrt(m_i m_j)(eV/Å²/amu)。λ<0 表示动力学矩阵负本征值(虚频),
    本函数会抛 ValueError —— 虚频走 :func:`smj.phonon.dynamical.imaginary_modes`。
    """
    if lam_ev_per_a2_amu < 0:
        raise ValueError(
            "λ<0(虚频)不能直接开方:%r;虚频请用 dynamical.imaginary_modes"
            % (lam_ev_per_a2_amu,))
    return EV_PER_SQRT * math.sqrt(lam_ev_per_a2_amu)


def ev_to_cm1(energy_ev: float) -> float:
    """能量 [eV] → 波数 [cm⁻¹](带符号,虚频可传负值)。"""
    return energy_ev * EV_TO_CM1


def gq_pw_per_k(temperature_k: float) -> float:
    """单量子通道热导上限 ``G_Q(T) = π²k_B²T/(3h)``,单位 pW/K。

    这是 Landauer 积分的硬上界:任一体系在任一温度下
    ``G_ph <= (声学支数)*G_Q``,单通道结 ``<= G_Q``。300 K 约 283.9 pW/K。

    必须用 SI 的 ``KB_J_PER_K``/``H_J_S`` 算完再乘 1e12 —— 若用 eV 口径的
    k_B/h 会差一个 ``E_CHARGE``(约 1.6e19 倍)。
    """
    if temperature_k <= 0:
        raise ValueError("温度必须为正:%r K" % (temperature_k,))
    g_w_per_k = math.pi ** 2 * KB_J_PER_K ** 2 * temperature_k / (3.0 * H_J_S)
    return g_w_per_k * 1.0e12


#: 元素的标准原子量(amu)。只收本仓真正会用到的元素:锚定基团 S/Se、
#: 电极金属 Au/Ag/Cu/Pt/Pd/Ni/Al、分子骨架 C/H/N/O 与常见取代卤素。
#: 来源:CIAAW 2021 标准原子量(常规值),与 FHI-aims species 默认值一致到 1e-4。
ELEMENT_MASS_AMU = {
    "H": 1.008, "He": 4.002602,
    "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998403163,
    "Si": 28.085, "P": 30.973761998, "S": 32.06, "Cl": 35.45,
    "Se": 78.971, "Br": 79.904,
    "Al": 26.9815385, "Ni": 58.6934, "Cu": 63.546,
    "Pd": 106.42, "Ag": 107.8682,
    "Pt": 195.084, "Au": 196.966569,
}
