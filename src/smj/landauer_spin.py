"""自旋极化 ``TE.dat`` 的双通道解析（不改 vendored ``landauer``）。

aitranss 在 ``nspin == 2`` 时把 ``TE.dat`` 写成两列透射（α 与 β 通道），而
vendored ``landauer.read_te`` 为了保持旧口径对两列取均值——总 G/S 不错，但
自旋极化信息在那一步就丢了。本模块按列读出 α/β，按同一套 L0 积分给出
每通道电导与 E_F 处的自旋极化率；总电导口径仍是
``G/G0 = (1/2) ∫ (Tα + Tβ) dE``（与 vendored 的均值一致）。
"""

from __future__ import annotations

from pathlib import Path

from .junction.landauer import (
    energy_order, minus_fermi_derivative, sort_by_energy, transport_moments, trapz)

#: ``P_S = (Sα−Sβ)/(Sα+Sβ)`` 的分母相对下限（2026-11-22 复审 P2）。
#: 分母是两通道 Seebeck 之差：Sα≈−Sβ 时比值对数值噪声极端敏感——实测
#: 分母 ~1e-6、|P_S|~2e6，打印出来是 ``+199999900%``，会被误读成"巨大自旋极化"。
#: 判据不是"分母严格等于 0"，而是"分母相对通道尺度近抵消"，故用 max(|Sα|,|Sβ|)
#: 作尺度；低于该阈值时 P_S 不可定义（None），并把原因写进 result。
SPIN_POLARIZATION_RTOL = 1e-6


def read_te_channels(path: str | Path) -> tuple[list[tuple[float, float, float, float]], int]:
    """``TE.dat`` -> ``[(E_Ha, E−E_F_eV, T_alpha, T_beta)]`` 与通道数（1 或 2）。

    单列文件里 α=β，通道数记 1；两列文件通道数记 2。非数据行（表头/空行）跳过。
    """
    rows: list[tuple[float, float, float, float]] = []
    channels = 1
    for line in Path(path).read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 3:
            continue
        try:
            energy = float(fields[0].replace("D", "E"))
            relative = float(fields[1].replace("D", "E"))
            alpha = float(fields[2].replace("D", "E"))
        except ValueError:
            continue
        beta = alpha
        if len(fields) >= 4:
            try:
                beta = float(fields[3].replace("D", "E"))
            except ValueError:
                beta = alpha
            else:
                channels = 2
        rows.append((energy, relative, alpha, beta))
    if len(rows) < 2:
        raise ValueError("Not enough TE data in %s" % path)
    return rows, channels


def _max_transmission(rows, index: int) -> dict:
    row = max(rows, key=lambda item: item[index])
    return {
        "E_Hartree": row[0],
        "E_minus_Ef_eV": row[1],
        "T": max(row[index], 0.0),
    }


def analyze_te_spin(path: str | Path, temperature: float = 300.0) -> dict:
    """双通道 Landauer 量：每通道 L0 / G/G0、E_F 处与积分后的自旋极化率。"""
    rows, channels = read_te_channels(path)
    # 2026-11-20 复审 P1-2: 行序不是物理量。旧实现按文件行序做 trapz, 乱序 TE.dat
    # 会得到与升序完全不同的 L0α/L0β(实测 0.541327 vs 0.249782)而总 L0 恰好不变,
    # 于是 results.json 里 "total 自洽、通道不自洽"。与 analyze_te 共用同一套排序。
    input_order = energy_order(rows)
    if input_order != "ascending":
        rows = sort_by_energy(rows)
    e_ev = [row[1] for row in rows]
    alpha = [max(row[2], 0.0) for row in rows]
    beta = [max(row[3], 0.0) for row in rows]
    weights = [minus_fermi_derivative(energy, temperature) for energy in e_ev]
    l0_alpha = abs(trapz([t * w for t, w in zip(alpha, weights)], e_ev))
    l0_beta = abs(trapz([t * w for t, w in zip(beta, weights)], e_ev))
    # 2026-11-21 批 B (S2): 每通道也积 L1/L2, 用与 analyze_te **同一** helper
    # 导出 Sα/Sβ、k_e、ZT —— 不复制公式, 否则"加权平均 == 总 S"无从校验。
    l1_alpha = trapz([e * t * w for e, t, w in zip(e_ev, alpha, weights)], e_ev)
    l1_beta = trapz([e * t * w for e, t, w in zip(e_ev, beta, weights)], e_ev)
    l2_alpha = abs(trapz([e * e * t * w for e, t, w in zip(e_ev, alpha, weights)], e_ev))
    l2_beta = abs(trapz([e * e * t * w for e, t, w in zip(e_ev, beta, weights)], e_ev))
    moments_alpha = transport_moments(l0_alpha, l1_alpha, l2_alpha, temperature)
    moments_beta = transport_moments(l0_beta, l1_beta, l2_beta, temperature)
    moments_total = transport_moments(
        0.5 * (l0_alpha + l0_beta),
        0.5 * (l1_alpha + l1_beta),
        0.5 * (l2_alpha + l2_beta),
        temperature)
    s_alpha = moments_alpha["S_uV_per_K"]
    s_beta = moments_beta["S_uV_per_K"]
    spin_polarization_s = None
    spin_polarization_s_note = None
    if s_alpha is not None and s_beta is not None:
        denominator = s_alpha + s_beta
        scale = max(abs(s_alpha), abs(s_beta), 1e-300)
        if abs(denominator) <= SPIN_POLARIZATION_RTOL * scale:
            spin_polarization_s_note = (
                "P_S 不可定义：Sα+Sβ=%.6g 相对通道尺度近抵消"
                "（|Sα+Sβ| ≤ %.0e·max(|Sα|,|Sβ|)），比值对数值噪声不稳定"
                % (denominator, SPIN_POLARIZATION_RTOL))
        else:
            spin_polarization_s = (s_alpha - s_beta) / denominator
    nearest = min(range(len(rows)), key=lambda i: abs(e_ev[i]))
    t_alpha_ef = alpha[nearest]
    t_beta_ef = beta[nearest]
    sum_ef = t_alpha_ef + t_beta_ef
    sum_l0 = l0_alpha + l0_beta
    result = {
        "source": str(path),
        "channels": channels,
        "temperature_K": temperature,
        "points": len(rows),
        "input_energy_order": input_order,
        "L0_alpha": l0_alpha,
        "L0_beta": l0_beta,
        "L1_alpha": l1_alpha,
        "L1_beta": l1_beta,
        "L2_alpha": l2_alpha,
        "L2_beta": l2_beta,
        "g_over_g0_alpha": 0.5 * l0_alpha,
        "g_over_g0_beta": 0.5 * l0_beta,
        "g_over_g0_total": 0.5 * sum_l0,
        "S_alpha_uV_per_K": s_alpha,
        "S_beta_uV_per_K": s_beta,
        "S_total_uV_per_K": moments_total["S_uV_per_K"],
        "k_e_alpha_pW_per_K": moments_alpha["k_e_pW_per_K"],
        "k_e_beta_pW_per_K": moments_beta["k_e_pW_per_K"],
        "ZT_alpha_electronic": moments_alpha["ZT_electronic"],
        "ZT_beta_electronic": moments_beta["ZT_electronic"],
        "spin_polarization_S": spin_polarization_s,
        "spin_polarization_S_note": spin_polarization_s_note,
        "T_alpha_at_ef": t_alpha_ef,
        "T_beta_at_ef": t_beta_ef,
        "spin_polarization_at_ef": (
            (t_alpha_ef - t_beta_ef) / sum_ef if sum_ef > 0 else None),
        "spin_polarization_L0": (
            (l0_alpha - l0_beta) / sum_l0 if sum_l0 > 0 else None),
        "max_transmission_alpha": _max_transmission(rows, 2),
        "max_transmission_beta": _max_transmission(rows, 3),
    }
    if channels == 1:
        result["note"] = (
            "TE.dat 为单列（α=β）：自旋极化恒为 0；若结确为开壳层，"
            "检查 FHI-aims 是否真的以 spin collinear 跑完")
    return result
