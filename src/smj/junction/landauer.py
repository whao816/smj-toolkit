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
    for line in Path(path).read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        if len(fields) >= 3:
            try:
                chans = [float(x.replace("D", "E")) for x in fields[2:4]]
                rows.append((float(fields[0]), float(fields[1]), sum(chans) / len(chans)))
            except ValueError:
                pass
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


def minus_fermi_derivative(e_ev: float, temperature: float) -> float:
    x = e_ev / (2.0 * KB_EV_PER_K * temperature)
    if abs(x) > 50:
        return 0.0
    c = math.cosh(x)
    return 1.0 / (4.0 * KB_EV_PER_K * temperature * c * c)


def analyze_te(path: str | Path, temperature: float = 300.0, output: str | Path | None = None) -> dict:
    rows = read_te(path)
    order = energy_order(rows)
    warnings: list[str] = []
    if order != "ascending":
        # 行序不是物理量: 排序后积分, 让升序/降序/乱序给出**同一**物理结果。
        # 不排序时梯形积分随 x 方向整体变号, 只有带符号的 L1 没被 abs 吸收,
        # 于是倒序 TE.dat 会静默把 S 反号(2026-11-20 复核 P1-1)。
        rows = sorted(rows, key=lambda r: r[1])
        warnings.append(
            "TE.dat 的能量列不是升序(%s); 已按 E-E_F 升序排序后积分" % order)
    e_ev = [r[1] for r in rows]
    trans = [max(r[2], 0.0) for r in rows]
    weights = [minus_fermi_derivative(e, temperature) for e in e_ev]
    l0_terms = [t * w for t, w in zip(trans, weights)]
    l1_terms = [e * t * w for e, t, w in zip(e_ev, trans, weights)]
    l2_terms = [e * e * t * w for e, t, w in zip(e_ev, trans, weights)]
    l0 = abs(trapz(l0_terms, e_ev))
    l1 = trapz(l1_terms, e_ev)
    l2 = abs(trapz(l2_terms, e_ev))
    result = {
        "source": str(path),
        "temperature_K": temperature,
        "points": len(rows),
        "L0": l0,
        "L1": l1,
        "L2": l2,
        "log10_G_over_G0_like": math.log10(l0) if l0 > 0 else None,
        "S_uV_per_K": 1e6 * (-(l1) / (l0 * temperature)) if l0 else None,
        "k_e_pW_per_K": 1.60217733e-7 * (l2 - (l1 * l1) / l0) * 2.0 / (H_EV_S * temperature) if l0 else None,
        "ZT_electronic": 1.0 / ((l0 * l2) / (l1 * l1) - 1.0) if l1 else None,
        "fermi_nearest": nearest_energy(rows, 0.0),
        "max_transmission": max(({"E_Hartree": r[0], "E_minus_Ef_eV": r[1], "T": r[2]} for r in rows), key=lambda x: x["T"]),
        # 行序观测: 让下游能区分"本来就是升序"与"被脚本纠正过"。
        "input_energy_order": order,
        "reordered": order != "ascending",
        "warnings": warnings,
    }
    if output:
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        Path(output).write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def nearest_energy(rows: Iterable[tuple[float, float, float]], target_e_ev: float) -> dict:
    row = min(rows, key=lambda r: abs(r[1] - target_e_ev))
    return {"E_Hartree": row[0], "E_minus_Ef_eV": row[1], "T": row[2]}