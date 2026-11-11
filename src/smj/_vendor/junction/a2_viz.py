"""A2 表面三元组审计 + 可视化 + 人工签署登记(AGD-2026-09-22-SMJ P0-4/P0-6 配套)。

为什么需要这个工具(而不是"跑个临时脚本出张图"):
- 首猜算法 electrode.surface_indices 只在"沿 tip->linker 最远投影层"里任选一个
  非共线三元组; 132 原子 FCC(111) 表面层有 21 个等价三角 -> **具体索引不唯一**,
  人眼要判断的是"这三个原子是否真的构成电极与分子之间的表面平面", 光看索引数字看不出来。
- 引擎 run 的产物按阶段子目录存放(au60_build/geometry.in), A2 确认必须绑定几何指纹
  (fingerprint.py): 审计、出图、登记三件事绑成一条命令, 才不会出现"图是一份几何、
  签名是另一份几何"的静默错位。

审计判据(纯本地可复算, 不消耗机时):
① 索引: 每侧恰好 3 个, 在 1..n, 元素为 Au, 左右不相交;
② 表面层: 三个原子同在最远投影层(沿 tip->linker, 容差 1e-3 A);
③ 三角形: 非共线(面积 > 0.5 A^2), 法向与输运轴夹角 < 1 度(否则 aitranss rotate 拿不到面);
④ 平面边界: 平面外侧无其它电极原子, 且 tip 在平面内侧(平面确实切在电极端头);
⑤ 与实际提交一致: 有 aitranss_te/tcontrol 时, 比对 $lsurc/$lsurx/$lsury/$rsurc/$rsurx/$rsury;
⑥ 索引可复现(W0-3, 2026-11-05): metadata 带 `surface_rule='det-v1'` 的 run 上,
   "按几何重算 == metadata" 是硬检查(不一致 = 索引不是几何的确定性函数);
   早于 det-v1 的 run 没有该键, 只把重算结果记进 detail, 不判死(它们是 legacy)。
   ①-⑤ 共 15 项(每侧 7 + tcontrol_match), 加上 ⑥ 后每侧 8 项、合计 17 项。

纯文件操作 + 只读几何, 不做真实计算。
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

from .electrode import surface_indices
from .fingerprint import (CONFIRM_FILE, attach_fingerprint, record_confirmation,
                         surface_fingerprint)
from .geometry import cross, dot, norm, read_aims_geometry, sub, unit

TOL_LAYER = 1e-3        # 同层判据: 与最远投影的差 < 1e-3 A(与 surface_indices 同口径)
MIN_AREA = 0.5          # 三角形面积下限 A^2: 低于此 aitranss rotate 视为退化
MAX_AXIS_ANGLE = 1.0    # 三角形法向与输运轴夹角上限(度)


def _triangle_metrics(atoms, chosen, away):
    """三元组几何: 三边长/面积/法向与 away 的夹角(度)。"""
    a, b, c = (atoms[i].coord() for i in chosen)
    edges = [norm(sub(b, a)), norm(sub(c, a)), norm(sub(c, b))]
    nrm = cross(sub(b, a), sub(c, a))
    area = 0.5 * norm(nrm)
    if norm(nrm) < 1e-9:
        return {'edges': edges, 'area': 0.0, 'normal_angle_deg': 90.0, 'normal': (0.0, 0.0, 0.0)}
    ang = math.degrees(math.acos(max(-1.0, min(1.0, abs(dot(unit(nrm), away))))))
    return {'edges': edges, 'area': area, 'normal_angle_deg': ang, 'normal': unit(nrm)}


def _axis(atoms, metadata):
    """输运轴 = 左锚 -> 右锚 单位方向(不预设 x/y/z, 与 validate_junction_geometry 同口径)。"""
    lk = int(metadata['left_linker_index']) - 1
    rk = int(metadata['right_linker_index']) - 1
    return unit(sub(atoms[rk].coord(), atoms[lk].coord())), atoms[lk].coord(), atoms[rk].coord()


def _cluster(atoms, metadata, side, axis, mid):
    """该侧电极 Au(不含对侧 tip): 用输运轴中点分侧, 不靠索引段猜。"""
    other_tip = int(metadata['right_tip_index' if side == 'left' else 'left_tip_index']) - 1
    out = []
    for i, atom in enumerate(atoms):
        if atom.element != 'Au' or i == other_tip:
            continue
        s = dot(sub(atom.coord(), mid), axis)      # >0 偏向右侧
        if (side == 'left' and s < 0) or (side == 'right' and s > 0):
            out.append(i)
    return out


def _side_audit(atoms, metadata, side, axis, mid):
    """单侧审计: 返回 checks/metrics(不抛错, 缺字段按未通过记录)。"""
    tip = int(metadata['%s_tip_index' % side]) - 1
    linker = int(metadata['%s_linker_index' % side]) - 1
    raw = list(metadata.get('%s_surface_indices' % side) or [])
    chosen = [int(i) - 1 for i in raw if str(i).strip()]
    away = unit(sub(atoms[tip].coord(), atoms[linker].coord()))
    cluster = _cluster(atoms, metadata, side, axis, mid)
    proj = {i: dot(atoms[i].coord(), away) for i in range(len(atoms))}
    top = max((proj[i] for i in cluster), default=0.0)
    layer = [i for i in cluster if top - proj[i] <= TOL_LAYER]
    checks = []

    def add(name, ok, detail, value=None):
        checks.append({'name': name, 'ok': bool(ok), 'detail': detail, 'value': value})

    in_range = [i for i in chosen if 0 <= i < len(atoms)]
    elems = [atoms[i].element for i in in_range]
    add('indices_shape', len(chosen) == 3 and len(in_range) == 3,
        '%s 表面索引 %s: 须恰好 3 个且在 1..%d' % (side, raw, len(atoms)), raw)
    add('indices_are_au', len(elems) == 3 and all(e == 'Au' for e in elems),
        '%s 表面索引元素: %s' % (side, elems), elems)
    others = set(int(i) for i in (metadata.get('right_surface_indices' if side == 'left'
                                              else 'left_surface_indices') or []))
    add('sides_disjoint', not (set(int(i) for i in raw) & others),
        '%s 与对侧表面索引不相交' % side, sorted(set(int(i) for i in raw) & others))

    # 表面层: 三个原子都必须在最远投影层(容差内)
    gaps = [top - proj[i] for i in chosen if i in proj]
    same_layer = bool(chosen) and all(g <= TOL_LAYER for g in gaps)
    add('same_surface_layer', same_layer,
        '%s 三个原子与最远投影层的最大偏差 %.6f A(容差 %g); 该层共 %d 个 Au'
        % (side, max(gaps, default=float('nan')), TOL_LAYER, len(layer)),
        {'layer_size': len(layer), 'layer_indices': [i + 1 for i in layer],
         'max_gap_A': max(gaps, default=None)})

    tm = _triangle_metrics(atoms, chosen, away) if len(in_range) == 3 else {}
    add('non_collinear', tm.get('area', 0.0) > MIN_AREA,
        '%s 三角形面积 %.4f A^2(下限 %g)' % (side, tm.get('area', 0.0), MIN_AREA),
        {'area_A2': tm.get('area'), 'edges_A': tm.get('edges')})
    add('plane_normal_along_axis', tm.get('normal_angle_deg', 90.0) < MAX_AXIS_ANGLE,
        '%s 三角形法向与输运轴夹角 %.3f 度(上限 %g)'
        % (side, tm.get('normal_angle_deg', 90.0), MAX_AXIS_ANGLE),
        {'angle_deg': tm.get('normal_angle_deg'), 'normal': tm.get('normal')})

    # 平面边界: 平面外侧(更靠 away)无电极原子, tip 在平面内侧
    plane = atom_plane = None
    if len(chosen) == 3:
        atom_plane = atoms[chosen[0]].coord()
        # 平面外侧 = 比最远层更靠外(正常应为空); 内侧 = tip 与分子都必须落在平面内
        beyond = [i + 1 for i in cluster if proj[i] - top > TOL_LAYER]
        tip_inside = dot(sub(atoms[tip].coord(), atom_plane), away) < 0
        mol_out = [i + 1 for i, a in enumerate(atoms)
                   if a.element != 'Au' and dot(sub(a.coord(), atom_plane), away) > 0]
        add('plane_is_boundary', not beyond and not mol_out and tip_inside,
            '%s 平面外侧电极原子 %s / 平面外分子原子 %s; tip 在内侧=%s'
            % (side, beyond[:6], mol_out[:6], tip_inside),
            {'beyond': beyond[:8], 'molecule_outside': mol_out[:8],
             'tip_inside': tip_inside})
        plane = {'point': atom_plane, 'normal': away}

    metrics = {
        'indices': raw, 'elements': elems,
        'tip_index': tip + 1, 'linker_index': linker + 1,
        'away': away, 'layer_size': len(layer), 'layer_indices': [i + 1 for i in layer],
        'cluster_atoms': len(cluster),
        'edges_A': tm.get('edges'), 'area_A2': tm.get('area'),
        'normal_angle_deg': tm.get('normal_angle_deg'),
        'tip_to_plane_A': (dot(sub(atoms[tip].coord(), atom_plane), away) if atom_plane else None),
        'equivalent_triple': None, 'equiv_same_set': None,
    }
    if cluster and len(chosen) == 3:
        # W0-3(2026-11-05): 索引规则改成确定性的 det-v1(surface_index.pick_surface_triple)。
        # 对 det-v1 的 run, "重算 == metadata" 是**硬检查**; 早于 det-v1 的 run 没有
        # surface_rule, 只记录不判死(它们已被 §11-S8 标成 legacy)。
        equiv = surface_indices(atoms, cluster, tip, linker)
        same_set = set(equiv) == set(i + 1 for i in chosen)
        rule = str(metadata.get('surface_rule') or 'legacy')
        metrics['equivalent_triple'] = equiv
        metrics['equiv_same_set'] = same_set
        metrics['equiv_same_layer'] = all(i in layer for i in (x - 1 for x in equiv))
        metrics['index_rule'] = rule
        if rule == 'det-v1':
            add('deterministic_indices', same_set,
                '%s det-v1 重算 %s vs metadata %s(不一致=索引不是几何的确定性函数)'
                % (side, equiv, sorted(i + 1 for i in chosen)), equiv)
        else:
            add('deterministic_indices', True,
                '%s metadata 无 surface_rule(legacy, 早于 det-v1): 重算 %s 与 metadata %s'
                ' 仅记录, 不判死 —— 新 run 一律带 det-v1 标记'
                % (side, equiv, sorted(i + 1 for i in chosen)), equiv)
    return {'side': side, 'checks': checks, 'metrics': metrics,
            'plane': plane, 'ok': all(c['ok'] for c in checks)}


def _tcontrol_values(path):
    """解析 tcontrol 里 aitranss 的六个表面原子索引; 缺键返回 None。"""
    vals = {}
    for raw in Path(path).read_text(encoding='utf-8', errors='ignore').splitlines():
        parts = raw.split()
        if len(parts) >= 2 and parts[0].startswith('$') and parts[0][1:] in (
                'lsurc', 'lsurx', 'lsury', 'rsurc', 'rsurx', 'rsury'):
            try:
                vals[parts[0][1:]] = int(parts[1])
            except ValueError:
                pass
    need = ('lsurc', 'lsurx', 'lsury', 'rsurc', 'rsurx', 'rsury')
    return vals if all(k in vals for k in need) else None


def surface_audit(geometry_path, metadata, tcontrol=None) -> dict:
    """审计一侧几何的 6 个表面索引是否构成合法表面三元组(不通过不抛错, 结果里给原因)。"""
    geometry_path = Path(geometry_path)
    atoms = read_aims_geometry(geometry_path)
    if not isinstance(metadata, dict):
        raise ValueError('metadata 必须是 dict(electrode_metadata.json 内容), 实为 %r'
                         % type(metadata).__name__)
    n = len(atoms)
    for k in ('left_tip_index', 'right_tip_index', 'left_linker_index', 'right_linker_index'):
        if not metadata.get(k):
            raise ValueError('metadata 缺 %s, 无法定向输运轴' % k)
    axis, lpos, rpos = _axis(atoms, metadata)
    mid = tuple((lpos[i] + rpos[i]) / 2.0 for i in range(3))
    sides = {s: _side_audit(atoms, metadata, s, axis, mid) for s in ('left', 'right')}

    checks = []
    both = [c for s in sides.values() for c in s['checks']]
    checks.extend(both)
    # tcontrol 一致性: 若本次提交真跑过, 签名的索引必须就是 aitranss 用的那 6 个
    tc_path = Path(tcontrol) if tcontrol else None
    tc = None
    if tc_path is not None and tc_path.is_file():
        vals = _tcontrol_values(tc_path)
        if vals is None:
            checks.append({'name': 'tcontrol_match', 'ok': False,
                           'detail': 'tcontrol 缺表面索引键: %s' % tc_path, 'value': None})
        else:
            ok = True
            detail = []
            for side, pre in (('left', 'l'), ('right', 'r')):
                got = [vals[pre + k] for k in ('surc', 'surx', 'sury')]
                exp = list(metadata.get('%s_surface_indices' % side) or [])
                same_set = set(got) == set(int(i) for i in exp)
                ok = ok and same_set
                detail.append('%s tcontrol=%s metadata=%s' % (side, got, exp))
            checks.append({'name': 'tcontrol_match', 'ok': ok,
                           'detail': 'tcontrol(%s) 与 metadata 表面索引一致: %s'
                                     % (tc_path, '; '.join(detail)),
                           'value': vals})
            tc = {'path': str(tc_path), 'values': vals, 'consistent': ok}
    return {
        'geometry': str(geometry_path),
        'fingerprint': surface_fingerprint(geometry_path),
        'atom_count': n,
        'axis': axis,
        'sides': sides,
        'tcontrol': tc,
        'checks': checks,
        'ok': all(c['ok'] for c in checks),
    }


def _plane_basis(normal):
    """构造与 normal 正交的一组基(e1, e2), 用于沿轴看表面层的平面视图。"""
    n = unit(normal)
    seed = (1.0, 0.0, 0.0) if abs(n[0]) < 0.9 else (0.0, 1.0, 0.0)
    e1 = unit(sub(seed, tuple(dot(seed, n) * n[i] for i in range(3))))
    e2 = cross(n, e1)
    return e1, e2


def render_figure(geometry_path, metadata, out_png, audit=None, title=None,
                  dpi=150) -> str:
    """出 A2 可视核对图: 侧视全景 + 左右表面沿轴视图(带索引标注)。

    matplotlib 是可选依赖(不在 pyproject 必装里): 只在真正出图时导入,
    缺失时抛 RuntimeError 并提示 --no-figure —— 审计与签署不依赖绘图。
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except Exception as exc:   # pragma: no cover - 环境缺 matplotlib 时走这里
        raise RuntimeError('出图需要 matplotlib(可选依赖): %s; 可加 --no-figure 只审计/签署'
                           % exc)
    geometry_path = Path(geometry_path)
    atoms = read_aims_geometry(geometry_path)
    audit = audit or surface_audit(geometry_path, metadata)
    axis, lpos, rpos = _axis(atoms, metadata)
    mid = tuple((lpos[i] + rpos[i]) / 2.0 for i in range(3))
    # 侧视的第二个基: 取电极质心相对中点的垂直分量(分子的 V 形张开方向), 退化则退到 z 轴
    au = [a.coord() for a in atoms if a.element == 'Au']
    cent = tuple(sum(c[i] for c in au) / len(au) for i in range(3))
    d = sub(cent, mid)
    v = sub(d, tuple(dot(d, axis) * axis[i] for i in range(3)))
    if norm(v) < 1e-6:
        seed = (0.0, 0.0, 1.0)
        v = sub(seed, tuple(dot(seed, axis) * axis[i] for i in range(3)))
    v = unit(v)

    colors = {'Au': '#9a9a9a', 'C': '#333333', 'H': '#c8c8c8', 'S': '#e0b000',
              'N': '#2e6fd0', 'O': '#d04040', 'F': '#40a040'}
    sizes = {'Au': 46, 'C': 30, 'H': 14, 'S': 40, 'N': 30, 'O': 30}
    surf_color = {'left': '#d02020', 'right': '#1040d0'}

    fig = plt.figure(figsize=(15.0, 9.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.55, 1.0])

    def _scatter(ax, xs, ys, atom_idx, scale=1.0, alpha=1.0):
        ax.scatter([xs[i] for i in atom_idx], [ys[i] for i in atom_idx],
                   s=[sizes.get(atoms[i].element, 24) * scale for i in atom_idx],
                   c=[colors.get(atoms[i].element, '#888888') for i in atom_idx],
                   edgecolors='#404040', linewidths=0.3, alpha=alpha, zorder=2)

    # --- (a) 侧视: 横轴 = 输运轴, 纵轴 = 张开方向 ---
    ax = fig.add_subplot(gs[0, :])
    xs = [dot(a.coord(), axis) for a in atoms]
    ys = [dot(a.coord(), v) for a in atoms]
    _scatter(ax, xs, ys, list(range(len(atoms))))
    for side in ('left', 'right'):
        idx = [int(i) - 1 for i in metadata.get('%s_surface_indices' % side) or []]
        col = surf_color[side]
        if len(idx) == 3:
            tri = [idx[0], idx[1], idx[2], idx[0]]
            ax.plot([xs[i] for i in tri], [ys[i] for i in tri], '-', color=col,
                    lw=1.6, zorder=4)
        ax.scatter([xs[i] for i in idx], [ys[i] for i in idx], s=300, facecolors='none',
                   edgecolors=col, linewidths=2.0, zorder=5, label='%s-surface' % side)
        for i in idx:
            ax.annotate(str(i + 1), (xs[i], ys[i]), textcoords='offset points',
                        xytext=(7, 5), fontsize=11, color=col, fontweight='bold', zorder=6)
    for key, label in (('%s_tip_index', 'tip'), ('%s_linker_index', 'linker')):
        for side, pre in (('left', 'L-'), ('right', 'R-')):
            k = key % side
            if metadata.get(k):
                i = int(metadata[k]) - 1
                ax.annotate(pre + label, (xs[i], ys[i]), textcoords='offset points',
                            xytext=(-38, -12), fontsize=9, color='#006060', zorder=6)
    for side, alpha in (('left', 0.55), ('right', 0.55)):
        idx = [int(i) - 1 for i in metadata.get('%s_surface_indices' % side) or []]
        if len(idx) == 3:
            pt = atoms[idx[0]].coord()
            c = dot(pt, axis)
            ax.axvline(c, color=surf_color[side], ls='--', lw=1.0, alpha=alpha, zorder=1)
    ax.set_xlabel('transport axis (A, left->right)')
    ax.set_ylabel('spread direction (A)')
    ax.set_title('(a) full junction side view: red/blue rings = confirmed surface triples',
                 fontsize=11)
    ax.set_aspect('equal', adjustable='datalim')
    ax.grid(alpha=0.25)
    ax.legend(loc='upper left', fontsize=9, frameon=True)

    # --- (b)(c) 沿轴看表面层: 平面内两个正交方向 ---
    for col_i, side in enumerate(('left', 'right')):
        ax = fig.add_subplot(gs[1, col_i])
        tip = int(metadata['%s_tip_index' % side]) - 1
        linker = int(metadata['%s_linker_index' % side]) - 1
        away = unit(sub(atoms[tip].coord(), atoms[linker].coord()))
        e1, e2 = _plane_basis(away)
        cluster = _cluster(atoms, metadata, side, axis, mid)
        layer = set(int(i) - 1 for i in meta_layer(audit, side))   # 审计里是 1-based
        chosen = [int(i) - 1 for i in metadata.get('%s_surface_indices' % side) or []]
        px = [dot(atoms[i].coord(), e1) for i in cluster]
        py = [dot(atoms[i].coord(), e2) for i in cluster]
        ax.scatter(px, py, s=26, c='#bbbbbb', edgecolors='#777777', linewidths=0.3,
                   zorder=1, label='electrode Au')
        if layer:
            lx = [dot(atoms[i].coord(), e1) for i in layer]
            ly = [dot(atoms[i].coord(), e2) for i in layer]
            ax.scatter(lx, ly, s=70, c='#7a7a7a', edgecolors='#333333',
                       linewidths=0.6, zorder=2, label='outermost layer')
        if len(chosen) == 3:
            tri = [chosen[0], chosen[1], chosen[2], chosen[0]]
            ax.plot([dot(atoms[i].coord(), e1) for i in tri],
                    [dot(atoms[i].coord(), e2) for i in tri], '-',
                    color=surf_color[side], lw=1.8, zorder=3)
        ax.scatter([dot(atoms[i].coord(), e1) for i in chosen],
                   [dot(atoms[i].coord(), e2) for i in chosen],
                   s=260, facecolors='none', edgecolors=surf_color[side], linewidths=2.0,
                   zorder=4)
        for i in chosen:
            ax.annotate(str(i + 1), (dot(atoms[i].coord(), e1), dot(atoms[i].coord(), e2)),
                        textcoords='offset points', xytext=(8, 6), fontsize=11,
                        color=surf_color[side], fontweight='bold', zorder=5)
        m = audit['sides'][side]['metrics']
        ax.set_title('(%s) %s surface seen along axis: layer=%d atoms, area=%.2f A^2, '
                     'normal-vs-axis=%.2f deg'
                     % ('b' if side == 'left' else 'c', side, m.get('layer_size') or 0,
                        m.get('area_A2') or 0.0, m.get('normal_angle_deg') or 0.0),
                     fontsize=10)
        ax.set_aspect('equal', adjustable='datalim')
        ax.set_xlabel('in-plane axis 1 (A)')
        ax.set_ylabel('in-plane axis 2 (A)')
        ax.grid(alpha=0.25)
        ax.legend(loc='upper left', fontsize=8)

    fig.suptitle(title or ('A2 surface-triple check: %s (n=%d, fingerprint %s)'
                           % (geometry_path.parent.name, audit['atom_count'],
                              audit['fingerprint'][:12])), fontsize=12)
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    return str(out_png)


def meta_layer(audit, side):
    """审计结果里该侧最远层的索引(1-based); 供绘图复用, 不重算。"""
    return list(audit['sides'][side]['metrics'].get('layer_indices') or [])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog='python -m agenticdft.tools.junction.a2_viz',
        description='A2 表面三元组: 几何审计 + 出图 + 人工签署登记(绑定几何指纹)')
    ap.add_argument('--stage', required=True,
                    help='au60_build/au10_build 阶段目录(含 geometry.in + electrode_metadata.json)')
    ap.add_argument('--png', default=None, help='输出图(缺省 <stage>/a2_surface_viz.png)')
    ap.add_argument('--json', default=None, help='输出审计 json(缺省 <stage>/a2_surface_audit.json)')
    ap.add_argument('--tcontrol', default=None,
                    help='aitranss_te/tcontrol(缺省按 <stage>/../aitranss_te/tcontrol 找)')
    ap.add_argument('--no-figure', action='store_true', help='只审计, 不出图(无 matplotlib 时用)')
    ap.add_argument('--bind-fingerprint', action='store_true',
                    help='把 surface_fingerprint 写进 electrode_metadata.json(不写确认记录)')
    ap.add_argument('--sign-who', default=None,
                    help='登记 A2 确认(who); 审计未通过时拒绝签署并返回 2')
    ap.add_argument('--sign-when', default=None, help='确认时间(缺省 UTC now)')
    args = ap.parse_args(argv)

    stage = Path(args.stage)
    geo = stage / 'geometry.in'
    mp = stage / 'electrode_metadata.json'
    if not geo.is_file():
        print('[A2] 缺几何: %s' % geo)
        return 1
    if not mp.is_file():
        print('[A2] 缺 metadata: %s' % mp)
        return 1
    metadata = json.loads(mp.read_text(encoding='utf-8'))
    tc = Path(args.tcontrol) if args.tcontrol else (stage.parent / 'aitranss_te' / 'tcontrol')
    audit = surface_audit(geo, metadata, tcontrol=tc if tc.is_file() else None)

    out_json = Path(args.json) if args.json else stage / 'a2_surface_audit.json'
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf-8')

    print('[A2] %s | n=%d | fingerprint=%s' % (geo, audit['atom_count'],
                                              audit['fingerprint']))
    for side in ('left', 'right'):
        m = audit['sides'][side]['metrics']
        print('  [%s] indices=%s layer=%d/%d atoms | edges=%s | area=%.4f A^2 | '
              'angle=%.3f deg | tip-plane=%.3f A | equivalent=%s'
              % (side, m['indices'], m['layer_size'], m['cluster_atoms'],
                 ['%.3f' % e for e in (m.get('edges_A') or [])],
                 m.get('area_A2') or 0.0, m.get('normal_angle_deg') or 0.0,
                 m.get('tip_to_plane_A') or 0.0, m.get('equivalent_triple')))
    for c in audit['checks']:
        if not c['ok']:
            print('  [FAIL] %s: %s' % (c['name'], c['detail']))
    if audit.get('tcontrol'):
        print('[A2] tcontrol: %s -> %s' % (audit['tcontrol']['path'],
                                           audit['tcontrol']['values']))
    print('[A2] audit %s -> %s' % ('PASS' if audit['ok'] else 'FAIL', out_json))

    if args.bind_fingerprint:
        # 只绑指纹(证据层), 不写 a2_confirmation.json —— 确认记录只能由人工签署产生
        meta2 = attach_fingerprint(stage)
        print('[A2] 已绑定几何指纹 -> %s (surface_fingerprint=%s, surface_indices=%s)'
              % (mp, meta2['surface_fingerprint'][:12], meta2['surface_indices']))

    png = Path(args.png) if args.png else stage / 'a2_surface_viz.png'
    if not args.no_figure:
        try:
            render_figure(geo, metadata, png, audit=audit,
                          title='A2 surface-triple check: %s' % stage.parent.name)
            print('[A2] figure -> %s' % png)
        except RuntimeError as exc:
            if args.sign_who:
                print('[A2] 拒绝在无图情况下按"人工目视"签署: %s' % exc)
                return 2
            print('[A2] 未出图: %s' % exc)

    if args.sign_who:
        if not audit['ok']:
            print('[A2] 审计未通过, 拒绝签署(须先修正几何/索引或重新建结)')
            return 2
        rec = record_confirmation(stage, who=args.sign_who, when=args.sign_when)
        print('[A2] 已登记人工确认 -> %s' % (stage / CONFIRM_FILE))
        print('     who=%s when=%s fingerprint=%s' % (rec['who'], rec['when'],
                                                     rec['fingerprint']))
        print('     indices=%s' % rec['indices'])
    return 0 if audit['ok'] else 1


if __name__ == '__main__':   # pragma: no cover
    sys.exit(main())
