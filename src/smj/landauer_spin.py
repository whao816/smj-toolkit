"""自旋极化 ``TE.dat`` 的双通道解析（不改 vendored ``landauer``）。

aitranss 在 ``nspin == 2`` 时把 ``TE.dat`` 写成两列透射（α 与 β 通道），而
vendored ``landauer.read_te`` 为了保持旧口径对两列取均值——总 G/S 不错，但
自旋极化信息在那一步就丢了。本模块按列读出 α/β，按同一套 L0 积分给出
每通道电导与 E_F 处的自旋极化率；总电导口径仍是
``G/G0 = (1/2) ∫ (Tα + Tβ) dE``（与 vendored 的均值一致）。
"""

from __future__ import annotations

from pathlib import Path

from ._vendor.junction.landauer import minus_fermi_derivative, trapz


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
    e_ev = [row[1] for row in rows]
    alpha = [max(row[2], 0.0) for row in rows]
    beta = [max(row[3], 0.0) for row in rows]
    weights = [minus_fermi_derivative(energy, temperature) for energy in e_ev]
    l0_alpha = abs(trapz([t * w for t, w in zip(alpha, weights)], e_ev))
    l0_beta = abs(trapz([t * w for t, w in zip(beta, weights)], e_ev))
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
        "L0_alpha": l0_alpha,
        "L0_beta": l0_beta,
        "g_over_g0_alpha": 0.5 * l0_alpha,
        "g_over_g0_beta": 0.5 * l0_beta,
        "g_over_g0_total": 0.5 * sum_l0,
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
