# -*- coding: utf-8 -*-
"""smj 独立复核: 脱离管线锚点重算 Landauer L0/L1/L2(2026-10-14)。

动机(承 docs/taskcard-smj-bdt-benchmark.md §8 与 report §3 "未机械化部分"):
旗舰算例 `bdt-trans-smj-aae585` 的 L0/L1/L2 一直只有**管线自身**的实现值
(`tools/junction/landauer.py`)。科学结论要可证伪, 必须有一条**独立**算路:
独立数值解析 + 不同的求积(Simpson vs 梯形) + CODATA 常数, 若两者一致, 才排除
"解析→积分→命名"链上的单位/口径事故(本项目历史上出过 27x Hartree/eV 同类事故)。

独立性**边界**(如实声明, 不夸大)
----------------------------------
- 物理定义**相同**(同一 Landauer 公式), 本模块不是第二套物理;
- 独立之处: ① 独立解析数值并归并通道(与主管线共用文件列数校验); ② Simpson 求积(vs 梯形);
  ③ CODATA 2018 `k_B`; ④ 全程显式单位账(不做 `_like` 近似命名)。
- 因此它能抓的是**实现级**事故(单位、常数、符号、权重、通道), 抓不到"物理模型选错"。

数值口径
--------
T(E) 取 `eV` 能量轴(E - E_F), 费米窗 `w(E) = -df/dE = 1/(4kT cosh^2(E/2kT))`:
    L0 = ∫ T(E) w(E) dE          (量子电导单位: G/G0 = L0, G0 = 2e^2/h)
    L1 = ∫ E T(E) w(E) dE
    L2 = ∫ E^2 T(E) w(E) dE
恒等校验(闭式): T=const=c 时 L0 = c, L1 = 0, L2 = c * pi^2 (k_B T)^2 / 3。
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

if __package__:
    from .landauer import te_fields
else:  # 支持直接运行 independent_check.py
    from landauer import te_fields

#: CODATA 2018 玻尔兹曼常数(eV/K) —— 与管线用的 8.6173324e-5 独立取源
KB_EV_PER_K = 8.617333262e-5
#: 对比容差(相对): 求积方案差异 + 常数 1e-7 量级 —— 1e-3 足以抓单位/口径事故
DEFAULT_RTOL = 1e-3


def parse_te(path) -> list:
    """独立解析 TE.dat -> [(E_Ha, E_minus_Ef_eV, T)]; 表头以 # 起首, 数据列空格分隔。"""
    rows = []
    for lineno, fields in te_fields(path):
        try:
            vals = [float(x.replace('D', 'E').replace('d', 'e')) for x in fields]
        except ValueError as exc:
            raise ValueError(f'TE.dat 第 {lineno} 行含无效数值') from exc
        if not all(math.isfinite(x) for x in vals):
            raise ValueError(f'TE.dat 第 {lineno} 行含非有限数值')
        # 自旋口径: aitranss 非极化写 1 列, 极化写 2 列 -> 取通道均值(与物理定义一致)
        chans = vals[2:4]
        rows.append((vals[0], vals[1], sum(chans) / len(chans)))
    if len(rows) < 2:
        raise ValueError('TE.dat 数据行不足(<2): %s' % path)
    return rows


def fermi_window(e_ev: float, temperature: float) -> float:
    """-df/dE @ E(相对 E_F, eV) —— 解析式(与 k_B 独立取源)。"""
    x = e_ev / (2.0 * KB_EV_PER_K * temperature)
    if abs(x) > 50.0:
        return 0.0
    c = math.cosh(x)
    return 1.0 / (4.0 * KB_EV_PER_K * temperature * c * c)


def _simpson(y: list, x: list) -> float:
    """复合 Simpson(奇数点用 1/3 律, 偶数点末段用梯形收尾)。"""
    n = len(x) - 1
    if n < 1:
        return 0.0
    h = x[1] - x[0]
    uniform = all(abs((x[i] - x[i - 1]) - h) < 1e-9 * max(1.0, abs(h)) for i in range(1, len(x)))
    if not uniform or n < 2:
        return sum(0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]) for i in range(1, len(x)))
    m = n if n % 2 == 0 else n - 1
    total = y[0] + y[m] + 4.0 * sum(y[i] for i in range(1, m, 2)) + 2.0 * sum(y[i] for i in range(2, m, 2))
    total *= h / 3.0
    if m != n:                                    # 末段奇数: 梯形收尾
        total += 0.5 * (y[n] + y[m]) * (x[n] - x[m])
    return total


def l0_l1_l2(rows: list, temperature: float = 300.0) -> dict:
    """独立积分: L0/L1/L2(不调用管线 landauer 的任何函数)。

    先按 E-E_F 升序排序, 与 :func:`landauer.analyze_te` 同口径 —— 行序不是物理量,
    但不排序时求积结果随 x 方向整体变号, 只有 L1 带符号, 倒序 TE.dat 会静默反号。
    """
    rows = sorted(rows, key=lambda r: r[1])
    e = [r[1] for r in rows]
    t = [max(r[2], 0.0) for r in rows]
    w = [fermi_window(x, temperature) for x in e]
    l0 = _simpson([ti * wi for ti, wi in zip(t, w)], e)
    l1 = _simpson([ei * ti * wi for ei, ti, wi in zip(e, t, w)], e)
    l2 = _simpson([ei * ei * ti * wi for ei, ti, wi in zip(e, t, w)], e)
    return {'points': len(rows), 'L0': abs(l0), 'L1': l1, 'L2': abs(l2)}


def derived(l: dict, temperature: float) -> dict:
    """由 L0/L1/L2 出可复核量(单位显式; 不做 _like 近似命名)。"""
    l0, l1, l2 = l['L0'], l['L1'], l['L2']
    out = {'G_over_G0': l0,
           'S_uV_per_K': (1e6 * (-l1 / (l0 * temperature))) if l0 else None,
           'ZT_electronic': (1.0 / ((l0 * l2) / (l1 * l1) - 1.0)) if l1 else None}
    return out


def compare(pipeline: dict, independent: dict, rtol: float = DEFAULT_RTOL,
            temperature: float = 300.0) -> dict:
    """逐量比管线值与独立值; 返回 {'ok', 'rtol', 'diffs'}。

    L1 是**抵消主导**量(对称谱下解析恒为 0): 两个实现各自 ~1e-9 的舍入差会让相对差
    失去判别力。故 L1 用积分自然尺度兜底 —— |ΔL1| <= rtol * L0 * k_B T 判为一致
    (k_B T 是费米窗的能量尺度); L0/L2 用相对差。兜底阈值在 diffs 里显式回报。"""
    diffs = []
    l0 = pipeline.get('L0') if pipeline.get('L0') is not None else independent.get('L0')
    l1_floor = rtol * abs(l0 or 0.0) * KB_EV_PER_K * temperature
    for k in ('L0', 'L1', 'L2'):
        a, b = pipeline.get(k), independent.get(k)
        if a is None or b is None:
            diffs.append({'key': k, 'pipeline': a, 'independent': b, 'rel': None, 'ok': False})
            continue
        denom = max(abs(a), abs(b), 1e-300)
        rel = abs(a - b) / denom
        if k == 'L1':
            ok = abs(a - b) <= max(rtol * denom, l1_floor)
            diffs.append({'key': k, 'pipeline': a, 'independent': b, 'rel': rel,
                          'atol_floor': l1_floor, 'ok': ok})
        else:
            diffs.append({'key': k, 'pipeline': a, 'independent': b, 'rel': rel, 'ok': rel <= rtol})
    return {'ok': all(d['ok'] for d in diffs), 'rtol': rtol, 'diffs': diffs}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog='junction_independent_check',
                                 description='smj 独立复核: 重算 Landauer L0/L1/L2(与管线实现对比)')
    ap.add_argument('--te', required=True, help='TE.dat 路径')
    ap.add_argument('--temperature', type=float, default=300.0, help='温度 K(缺省 300)')
    ap.add_argument('--pipeline-json', default=None,
                    help='管线结果 JSON(含 L0/L1/L2; 如 run 的 junction 字段), 给了才做对比')
    ap.add_argument('--rtol', type=float, default=DEFAULT_RTOL, help='对比相对容差(缺省 1e-3)')
    ap.add_argument('--json', action='store_true')
    args = ap.parse_args(argv)
    rows = parse_te(args.te)
    ind = l0_l1_l2(rows, args.temperature)
    out = {'te': args.te, 'temperature_K': args.temperature, 'independent': ind,
           'derived': derived(ind, args.temperature)}
    if args.pipeline_json:
        pipe = json.loads(Path(args.pipeline_json).read_text(encoding='utf-8'))
        pipe = pipe.get('junction', pipe)
        out['pipeline'] = {k: pipe.get(k) for k in ('L0', 'L1', 'L2')}
        out['comparison'] = compare(pipe, ind, args.rtol, args.temperature)
    if args.json:
        print(json.dumps(out, ensure_ascii=False, indent=1, sort_keys=True))
    else:
        print('独立复核 %s(%d 点, %.0f K)' % (args.te, ind['points'], args.temperature))
        print('  L0=%.9f  L1=%.9e  L2=%.9e' % (ind['L0'], ind['L1'], ind['L2']))
        if 'comparison' in out:
            c = out['comparison']
            print('  与管线对比: %s(rtol=%g)' % ('一致' if c['ok'] else '**不一致**', c['rtol']))
            for d in c['diffs']:
                print('    %-3s 管线=%.9g 独立=%.9g 相对差=%s' % (
                    d['key'], d['pipeline'] if d['pipeline'] is not None else float('nan'),
                    d['independent'], ('%.2e' % d['rel']) if d['rel'] is not None else 'n/a'))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
