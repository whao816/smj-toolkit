"""动力学矩阵、虚频检查、声学求和(ASR)残差。

**单位链**(方案 §3.4,最容易错的一环):力常数 ``K`` 是 eV/Å²;质量加权后
``D_ij = K_ij / sqrt(m_i m_j)`` 的单位是 eV/(Å²·amu),**不能直接开方当 eV**。
必须先乘 ``EV_PER_SQRT²`` 才把本征值变成 ``(ħω)²``、单位 eV²::

    D̃ = EV_PER_SQRT² · K_ij / sqrt(m_i m_j)          单位 eV²
    λ  = eig(D̃)                                     (ħω)²
    ħω = sqrt(λ)                                    eV
    ν̃  = sign(λ)·sqrt(|λ|)·8065.544                 cm⁻¹

漏掉 ``EV_PER_SQRT`` 会让频率偏 15.5 倍(把 ħω 当成了 ω),这是调研报告 §6
列的第一号陷阱。

**为什么要单独查 ASR**:有限位移得到的 ``K`` 只在截断半径内求和,``Σ_j K_ij``
一般不为零 —— 不为零就意味着整个体系「凭空受力」,低频声学支会漂。注意残差必须
在 **K** 上算:质量加权后的 ``D̃·1`` 每个原子分量不同,拿它当平移模式残差会把
质量差异误判成 ASR 违反。
"""
from __future__ import annotations

import numpy as np

from . import units

__all__ = [
    "dynamical_matrix",
    "mode_frequencies_cm1",
    "imaginary_modes",
    "acoustic_sum_residual",
    "acoustic_sum_relative_residual",
    "acoustic_sum_residuals",
]


def dynamical_matrix(force_constants, masses, *, scale_to_ev2: bool = True,
                     symmetrize: bool = True):
    """``K``(eV/Å²)→ 质量加权动力学矩阵,缺省直接给出 eV² 口径。

    Parameters
    ----------
    force_constants : (3N,3N) array
        力常数矩阵 ``K[iα, jβ] = ∂²E/∂u_iα∂u_jβ``,eV/Å²。
    masses : (N,) array
        每原子质量,amu。三个笛卡尔自由度共用同一质量。
    scale_to_ev2 : bool
        True(缺省)乘 ``EV_PER_SQRT²``,本征值即 ``(ħω)²``、单位 eV² —— 正是
        ``negf.transmission`` 要的输入。False 时返回 eV/(Å²·amu),只在做单位
        推导或对账时用。
    symmetrize : bool
        中心差分算出 ``K`` 可能在机器精度上不对称;取 ``(K+K^T)/2`` 保证
        ``eigvalsh`` 的结果是实数。不要靠它「修」真正的非对称物理错误。

    Notes
    -----
    不在这里做 ASR 投影。方案 §5.1 明确要求**先测残差再投影**:若一开始就构造性
    满足 ``Σ_j K_ij = 0``,``acoustic_sum_residual`` 会退化成永远为 0 的空检查。
    """
    k = np.asarray(force_constants, dtype=float)
    m = np.asarray(masses, dtype=float).reshape(-1)
    natoms = m.size
    if k.shape != (3 * natoms, 3 * natoms):
        raise ValueError(
            "force_constants 形状 %r 与 %d 个原子(期望 (%d,%d))不匹配"
            % (k.shape, natoms, 3 * natoms, 3 * natoms))
    if np.any(~np.isfinite(k)):
        raise ValueError("force_constants 含非有限值")
    if np.any(m <= 0.0):
        raise ValueError("质量必须为正:%r" % (m[m <= 0.0].tolist(),))

    if symmetrize:
        k = 0.5 * (k + k.T)
    m3 = np.repeat(m, 3)
    d = k / np.sqrt(np.outer(m3, m3))
    if scale_to_ev2:
        d = d * units.EV_PER_SQRT ** 2
    return d


def mode_frequencies_cm1(dynamical_ev2):
    """``D̃``(eV²)的本征频率,单位 cm⁻¹;**虚频返回负值**。

    排序从小到大(虚频在前),便于直接看最低模。
    """
    d = np.asarray(dynamical_ev2, dtype=float)
    if d.ndim != 2 or d.shape[0] != d.shape[1]:
        raise ValueError("dynamical 必须是方阵,得到 %r" % (d.shape,))
    lam = np.linalg.eigvalsh(0.5 * (d + d.T))
    # λ≥0 → 实频;λ<0 → 虚频(记负号),两者都按 |λ|^{1/2} 换算
    return np.sign(lam) * units.ev_to_cm1(np.sqrt(np.abs(lam)))


def imaginary_modes(dynamical_ev2, *, cm1_tol: float = -10.0):
    """列出**超出容差**的虚频模。

    有限位移 + 几何未完全收敛时,总会有几个 ``-1~-10 cm⁻¹`` 的噪声模;方案 §6
    允许到 ``-10 cm⁻¹``。这里有三个用途不同的口径:

    - ``count``:被标记的模数(应判为不通过时 >0);
    - ``worst_cm1``:最负的那个(越小越糟);
    - ``modes``:逐条 ``{index, cm1}``,按 cm1 升序。

    ``cm1_tol`` 是**负阈值**:``cm1 < cm1_tol`` 才标记。默认 -10 表示
    ``-10 cm⁻¹`` 以内的噪声不算。
    """
    if cm1_tol > 0:
        raise ValueError(
            "cm1_tol 应为负数(如 -10.0),表示允许到 -10 cm⁻¹;得到 %r" % (cm1_tol,))
    freqs = mode_frequencies_cm1(dynamical_ev2)
    flagged = [
        {"index": int(i), "cm1": float(v)}
        for i, v in enumerate(freqs) if v < cm1_tol
    ]
    flagged.sort(key=lambda item: item["cm1"])
    return {
        "count": len(flagged),
        "worst_cm1": flagged[0]["cm1"] if flagged else 0.0,
        "modes": flagged,
        "tolerance_cm1": float(cm1_tol),
        "n_modes": int(freqs.size),
    }


def acoustic_sum_residuals(force_constants, rows=None) -> tuple[float, float]:
    """逐平移方向求和,返回绝对残差 max_iαβ |Σ_j K_iα,jβ| (eV/Å²) 和相对残差。

    对完全收敛的有限位移集,这个数应在 ``1e-4`` 量级以下;若与 ``max|K|`` 同量级,
    说明截断半径太小或位移-力配对错位。

    ``rows`` 只在"行不完整"的矩阵上才需要给: 拼接后的完整矩阵里,lead 行来自
    被截断成两个主层的 lead 力常数,它们与第三层及更远的耦合根本不在矩阵里,
    行和按构造就不可能闭合。只有器件区的行其邻居集合是完整的,判据必须限定在
    这些行上(见 :func:`smj.phonon.fc.asr_diagnostics`)。
    """
    k = np.asarray(force_constants, dtype=float)
    if k.ndim != 2 or k.shape[0] != k.shape[1] or not k.size or k.shape[1] % 3:
        raise ValueError("force_constants 必须是非空 3N×3N 方阵,得到 %r" % (k.shape,))
    if rows is not None:
        index = np.asarray(rows, dtype=int).ravel()
        if index.size == 0:
            raise ValueError("rows 不能为空")
        if index.min() < 0 or index.max() >= k.shape[0]:
            raise ValueError("rows 越界: 最大 %d,矩阵只有 %d 行"
                             % (index.max(), k.shape[0]))
        k = k[index, :]
    sums = k.reshape(k.shape[0], -1, 3).sum(axis=1)
    absolute = float(np.max(np.abs(sums)))
    scale = float(np.max(np.abs(k)))
    return absolute, absolute / scale if scale else 0.0


def acoustic_sum_residual(force_constants, rows=None) -> float:
    """ASR 绝对残差 max_iαβ |Σ_j K_iα,jβ| (eV/Å²)。"""
    return acoustic_sum_residuals(force_constants, rows)[0]


def acoustic_sum_relative_residual(force_constants, rows=None) -> float:
    """ASR 残差除以受检行的 max|K| (无量纲)。"""
    return acoustic_sum_residuals(force_constants, rows)[1]


def masses_from_elements(elements) -> "np.ndarray":
    """元素符号列表 → 原子量数组(amu)。

    质量表在 :mod:`smj.phonon.units` 里(与其它常量同处),未知元素直接报错 ——
    猜一个质量会让动力学矩阵的每一行都错,而且错得很像"物理"。
    """
    import numpy as np

    from .units import ELEMENT_MASS_AMU

    masses = []
    for element in elements:
        try:
            masses.append(ELEMENT_MASS_AMU[element])
        except KeyError:
            raise ValueError("未知元素 %r:请先把它的原子量加进 units.ELEMENT_MASS_AMU"
                             % (element,)) from None
    return np.asarray(masses, dtype=float)
