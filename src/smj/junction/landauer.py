"""单分子结 Landauer 后处理:TE.dat 解析、L0/L1/L2 积分、G/S/k_e/ZT。

来源: 第三方工具 qmtransport/transport.py(74 行, tmp/自动化/, 只读来源不删除不入库)
按本项目协议适配迁入(AGD-2026-09-22-SMJ, P2):改名 landauer 避免与既有 transport
性质(体材料载流子迁移率)冲突;公式逐字保留 —— 五值管线回归锚点
(log10_G_over_G0_like 等)已冻结于 agenticdft/tests/data/junction_anchors.json(P0)。
**免责(原样保留)**:锚点来自 7 个能量点的示例数据, 只能验证解析→积分→命名管线的
稳定性, 不能作为生产级科学结论; `_like` 后缀表示定义近似, 不得当真实量报告。
纯本地工具: 只做后处理, 不做真实计算。
"""
from __future__ import annotations

from pathlib import Path
import json
import math
from typing import Iterable


KB_EV_PER_K = 8.6173324e-5
H_EV_S = 4.135667662e-15

# 求积充分性只在 E_F ± FERMI_WINDOW_KBT·kBT 内检查（2026-11-21 复审 R1）：
# 费米窗权重的 99.9% 落在 ±4 kBT 内，检查范围再放大会被远离 E_F 的粗网格稀释。
FERMI_WINDOW_KBT = 4.0


def read_te(path: str | Path) -> list[tuple[float, float, float]]:
    """TE.dat -> [(E_Ha, E-EF_eV, T)]。

    P2-1b(自旋口径): aitranss 非极化时写**单列** T(per spin), 极化(nspin==2)时写
    **两列** T(E).alpha / T(E).beta(见 trans.f 的表头分支)。原实现固定取 fields[2],
    对极化文件会**静默丢掉 β 通道**(G 少算约 2x) —— 故此处改为按通道取均值。

    口径依据(已核实的 aitranss 自述): `G = T(EFermi)*[e^2/h] per spin channel` 与
    G0 = 2e^2/h  =>  G/G0 = (1/2) * int (T_alpha + T_beta) w dE。
    退化情形(T_alpha = T_beta = T)自动回到单列结果, 与改动前**逐位一致**。
    """
    rows = []
    for lineno, line in enumerate(
            Path(path).read_text(encoding="utf-8", errors="ignore").splitlines(), 1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        if len(fields) >= 3:
            try:
                e_ha = float(fields[0].replace("D", "E"))
                e_ev = float(fields[1].replace("D", "E"))
                chans = [float(x.replace("D", "E")) for x in fields[2:4]]
            except ValueError:
                continue
            # 2026-11-22 复审 P0-4：NaN/inf 会一路穿过 max()、trapz 和
            # transport_moments 的 `if l0` 守卫，把 G/S/k_e/ZT 变成 NaN，
            # 最后让机械判据输出"否证 H0"。损坏的证据必须在这里拒绝，
            # 由调用方归到 te_present(missing=True) -> inconclusive。
            if not all(math.isfinite(v) for v in (e_ha, e_ev, *chans)):
                raise ValueError(
                    "TE.dat 含非有限数值(第 %d 行: %s)；损坏的谱不能当物理结果"
                    % (lineno, " ".join(fields[:4])))
            rows.append((e_ha, e_ev, sum(chans) / len(chans)))
    if len(rows) < 2:
        raise ValueError(f"Not enough TE data in {path}")
    return rows


def trapz(y: list[float], x: list[float]) -> float:
    return sum(0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]) for i in range(1, len(x)))


def energy_order(rows: list[tuple[float, float, float]]) -> str:
    """E-E_F(eV) 列的行序: ``ascending`` / ``descending`` / ``unsorted``。

    TE.dat 的行序**不是物理量**: 积分只依赖 (E, T) 集合。但梯形求积的符号随 x 的
    单调方向翻转, 而旧实现直接按行序积分、且只对 L0/L2 取 abs —— 于是**倒序文件会
    静默把 Seebeck 反号**(L1 是唯一带符号的矩, 不受 abs 保护)。这里只做观测,
    真正的修正在 :func:`analyze_te` 里先排序。
    """
    e = [r[1] for r in rows]
    if all(b > a for a, b in zip(e, e[1:])):
        return "ascending"
    if all(b < a for a, b in zip(e, e[1:])):
        return "descending"
    return "unsorted"


def sort_by_energy(rows: Iterable[tuple]) -> list[tuple]:
    """按 ``E−E_F``（第 2 列）升序排序，3 列/4 列行元组通用。

    单一出口：积分、峰检测、绘图都从同一套有序行出发 —— 2026-11-20 复审 P1-3
    的根因正是"积分自己排了序，峰检测又按文件行序读了一遍"。
    """
    return sorted(rows, key=lambda row: row[1])


def read_te_ordered(path: str | Path) -> tuple[list[tuple[float, float, float]], str]:
    """TE.dat -> (按 ``E−E_F`` 升序的行, 原始行序)。

    行序不是物理量（见 :func:`energy_order`），但峰检测、绘图、窗口判据都依赖
    "相邻行 = 相邻能量"这一前提，故统一从这里取数。返回的 ``order`` 是**文件
    原始**行序，供下游记录"是否被纠正过"。
    """
    rows = read_te(path)
    order = energy_order(rows)
    return (sort_by_energy(rows) if order != "ascending" else rows), order


def minus_fermi_derivative(e_ev: float, temperature: float) -> float:
    x = e_ev / (2.0 * KB_EV_PER_K * temperature)
    if abs(x) > 50:
        return 0.0
    c = math.cosh(x)
    return 1.0 / (4.0 * KB_EV_PER_K * temperature * c * c)


def fermi(e_ev: float, temperature: float) -> float:
    """费米-狄拉克占据 ``f(E)``，``E`` 为 ``E − E_F``（eV）。"""
    x = e_ev / (KB_EV_PER_K * temperature)
    if x > 50.0:
        return 0.0
    if x < -50.0:
        return 1.0
    return 1.0 / (1.0 + math.exp(x))


def thermal_window(rows: Iterable[tuple], temperature: float) -> dict:
    """费米窗覆盖：端点占据差 ``f(E_min) − f(E_max)`` 与端点相对 kBT 的位置。

    这是"截窗积分有没有代表性"的唯一直接量（2026-11-20 复审 P1-1）：
    覆盖≈1 时 ``∫T·(−∂f/∂E)dE`` 才是 E_F 附近的展开；覆盖 << 1 时 G 被系统性
    低估、S 的符号可被端点截断带偏。判据不是它本身，而是 ``coverage`` 是否
    超过 ``criteria.MIN_THERMAL_COVERAGE``。
    """
    energies = [r[1] for r in rows]
    e_min, e_max = min(energies), max(energies)
    kbt = KB_EV_PER_K * temperature
    return {
        "e_min_ev": e_min,
        "e_max_ev": e_max,
        "kbt_ev": kbt,
        "lower_kbt": e_min / kbt,
        "upper_kbt": e_max / kbt,
        "coverage": fermi(e_min, temperature) - fermi(e_max, temperature),
    }


def fermi_quadrature(rows: Iterable[tuple], temperature: float) -> dict:
    """费米窗求积充分性：权重积分/精确覆盖 与 权重一阶矩（2026-11-21 复审 R1）。

    :func:`thermal_window` 只看"窗口端点是否把权重算全"，与网格密度无关；真正的
    积分误差取决于网格在 E_F 附近的分辨率。12 点 / 0.2 eV 步长、300 K 时 ±kBT
    内只有 1 个采样点：梯形求积把 ``∫w dE`` 低估 42%、一阶矩偏到 −1.96 kBT，于是
    常数谱也报出 S=+169 uV/K。两个比值都是无量纲量，与体系无关。
    """
    ordered = sort_by_energy(rows)
    e_ev = [r[1] for r in ordered]
    kbt = KB_EV_PER_K * temperature
    weights = [minus_fermi_derivative(e, temperature) for e in e_ev]
    integral = trapz(weights, e_ev)
    exact = fermi(e_ev[0], temperature) - fermi(e_ev[-1], temperature)
    moment = trapz([e * w for e, w in zip(e_ev, weights)], e_ev)
    in_window = [e for e in e_ev if abs(e) <= FERMI_WINDOW_KBT * kbt]
    steps = [b - a for a, b in zip(in_window, in_window[1:])]
    return {
        "kbt_ev": kbt,
        "window_kbt": FERMI_WINDOW_KBT,
        "weight_integral": integral,
        "exact_coverage": exact,
        "ratio": integral / exact if exact > 0 else None,
        "first_moment_kbt": moment / (kbt * integral) if integral else None,
        "points_in_window": len(in_window),
        "min_step_in_window_ev": min(steps) if steps else None,
    }


def transport_moments(l0: float, l1: float, l2: float, temperature: float) -> dict:
    """由 ``L0/L1/L2`` 导出的输运系数（``analyze_te`` 与分通道共用同一套公式）。

    2026-11-21 批 B (S2): 分通道 S/k_e/ZT 必须与总口径完全同形，否则"加权平均
    等于总 S"的自洽检查就失去意义。这里只做公式搬运，不做任何近似。
    """
    if not all(math.isfinite(v) for v in (l0, l1, l2)):
        # 理论上 read_te 已拒绝非有限输入；这里的守卫是给直接调用者/未来
        # 数据源（spin 分通道）兜底，避免 NaN 被 `if l0` 当成真值继续传播。
        return {
            "log10_G_over_G0_like": None,
            "S_uV_per_K": None,
            "k_e_pW_per_K": None,
            "ZT_electronic": None,
        }
    return {
        "log10_G_over_G0_like": math.log10(l0) if l0 > 0 else None,
        "S_uV_per_K": 1e6 * (-(l1) / (l0 * temperature)) if l0 else None,
        "k_e_pW_per_K": (1.60217733e-7 * (l2 - (l1 * l1) / l0) * 2.0
                         / (H_EV_S * temperature)) if l0 else None,
        "ZT_electronic": 1.0 / ((l0 * l2) / (l1 * l1) - 1.0) if l1 else None,
    }


def analyze_te(path: str | Path, temperature: float = 300.0, output: str | Path | None = None) -> dict:
    rows, order = read_te_ordered(path)
    warnings: list[str] = []
    if order != "ascending":
        # 行序不是物理量: 排序后积分, 让升序/降序/乱序给出**同一**物理结果。
        # 不排序时梯形积分随 x 方向整体变号, 只有带符号的 L1 没被 abs 吸收,
        # 于是倒序 TE.dat 会静默把 S 反号(2026-11-20 复核 P1-1)。
        warnings.append(
            "TE.dat 的能量列不是升序(%s); 已按 E-E_F 升序排序后积分" % order)
    e_ev = [r[1] for r in rows]
    raw_trans = [r[2] for r in rows]
    negative = [t for t in raw_trans if t < 0.0]
    if negative:
        warnings.append(
            "TE.dat 有 %d 个负透射值(最小 %.6g)：按物理口径 clamp 为 0 后积分，"
            "max_transmission 仍记录原始峰值" % (len(negative), min(negative)))
    trans = [max(t, 0.0) for t in raw_trans]
    weights = [minus_fermi_derivative(e, temperature) for e in e_ev]
    l0_terms = [t * w for t, w in zip(trans, weights)]
    l1_terms = [e * t * w for e, t, w in zip(e_ev, trans, weights)]
    l2_terms = [e * e * t * w for e, t, w in zip(e_ev, trans, weights)]
    l0 = abs(trapz(l0_terms, e_ev))
    l1 = trapz(l1_terms, e_ev)
    l2 = abs(trapz(l2_terms, e_ev))
    moments = transport_moments(l0, l1, l2, temperature)
    result = {
        "source": str(path),
        "temperature_K": temperature,
        "points": len(rows),
        "L0": l0,
        "L1": l1,
        "L2": l2,
        **moments,
        "fermi_nearest": nearest_energy(rows, 0.0),
        "thermal_window": thermal_window(rows, temperature),
        "fermi_quadrature": fermi_quadrature(rows, temperature),
        "max_transmission": max(({"E_Hartree": r[0], "E_minus_Ef_eV": r[1], "T": r[2]} for r in rows), key=lambda x: x["T"]),
        # 行序观测: 让下游能区分"本来就是升序"与"被脚本纠正过"。
        "input_energy_order": order,
        "reordered": order != "ascending",
        "negative_transmission_points": len(negative),
        "warnings": warnings,
    }
    if output:
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        Path(output).write_text(json.dumps(result, indent=2, allow_nan=False),
                                encoding="utf-8")
    return result


def nearest_energy(rows: Iterable[tuple[float, float, float]], target_e_ev: float) -> dict:
    row = min(rows, key=lambda r: abs(r[1] - target_e_ev))
    return {"E_Hartree": row[0], "E_minus_Ef_eV": row[1], "T": row[2]}