"""确定性表面三角(W0-3, 2026-11-05)。

aitranss 的 ``$lsurc/$lsurx/$lsury``(及右侧三个)要求每侧给 **3 个电极表面原子**,
它们定义一个平面与坐标系。旧实现 ``electrode.surface_indices`` 有两个不确定性来源:

1. 层判据是**绝对**容差 ``1e-3 Å`` —— 与几何尺度无关, 电极一加厚/一倾斜就可能
   把两层并成一层或把一层拆成两层;
2. 层内用"三边平方和最小"挑三角, 平局由**遍历顺序**决定。FCC(111) 外层是一整片
   三角格, 与接触点等价的单位三角往往不止一个(实测 6 个, 质心差 < 1e-8 Å),
   于是"同一份几何、不同次运行"可能给出不同索引 —— 三种 ``opt`` 设置曾各给一套,
   且三套都能通过审计。

本模块把选择改成**几何的确定性函数**: 同样的几何坐标 -> 同样的三个原子(与原子在
文件里的先后顺序无关)。规则:

a) 侧向成员: 与 tip 同元素的 Au、排除对侧 tip、按输运轴中点分侧(与 ``a2_viz`` 同口径);
b) 外层: 沿 ``away``(tip - linker)投影最远的层, 容差 = ``max(1e-5, 0.02 x 层间距)``;
c) 参考点 q: tip 沿 ``away`` 投影到外层平面上(即"接触点正对的那一小块表面");
d) 三角排序键: ``(桶(最长边), 桶(|质心 - q|), 桶(面积), 顶点坐标字典序)`` ——
   前三项按 ``BUCKET_A``(=0.05 Å, 约 2% 晶格常数)分桶, 于是"对称等价的多个三角"
   落在同一个桶里, 胜负由坐标字典序定, 而不是由浮点噪声或遍历顺序定; 第四项用
   坐标而不是索引, 所以换一个原子顺序不会换一组原子;
e) 桶化的必要性: 实测对称几何里多个单位三角的 gap 差 < 1e-8 Å, 而 SCF 收敛后
   坐标自带 ~1e-3 Å 噪声 —— 用 1e-6 Å 精度比较, 1e-4 Å 的抖动就能换掉三个原子
   (W0-3 的抖动测试就是这样失败的)。分桶把这个自由度吃掉;
f) 若与最优同分(同桶)的有多个, 记 ``tie_count``: 审计据此知道"这里是对称退化,
   选择等价但必须一致"。

退化(R7)与回退: 外层不足 3 个原子(极薄电极)时按"离 q 最近的 3 个且非共线"取,
并在结果里标 ``fallback=True``, 不静默。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .geometry import cross, dot, norm, sub, unit

LAYER_TOL_FRAC = 0.02      # 层容差 = 层间距 x 该比例(下限 1e-5 A)
MIN_LAYER_TOL = 1e-5
SCORE_NDIGITS = 6          # (保留给调用方做人类可读比较)
BUCKET_A = 0.05            # 判分量化的桶宽(A): 0.05 A ≈ 2% 晶格常数
COORD_NDIGITS = 3          # 破平局用的坐标量化位数(0.001 A)


def _bucket(x: float, width: float = BUCKET_A) -> int:
    """把连续判分量放进宽度 width 的桶。

    为什么要桶: "质心到接触点投影的距离"这类量在**对称几何**下会有多个三角落在
    1e-8 A 内(平局), 而在**近似对称**几何(SCF 收敛到 1e-3 A 量级)下顺序会翻 ——
    用 1e-6 A 精度比较等于让浮点噪声决定选哪个原子。按 0.05 A 分桶后, 噪声只影响
    "桶内谁更近"这种本来就等价的选择, 最终由坐标字典序定案, 结果是几何的稳定函数。
    """
    return int(round(float(x) / width))


@dataclass
class SurfaceTriple:
    """一侧表面三角的选择结果(索引一律 1-based, 与 metadata 同口径)。"""

    indices: list[int]
    side: str = ''
    layer_size: int = 0
    layer_indices: list[int] = field(default_factory=list)
    spacing_A: float = 0.0
    tolerance_A: float = 0.0
    centroid_gap_A: float = 0.0
    area_A2: float = 0.0
    longest_edge_A: float = 0.0
    tie_count: int = 1
    fallback: bool = False
    candidates: int = 0
    rule: str = 'det-v1'

    @property
    def symmetric_tie(self) -> bool:
        """是否属于"多个等价三角同分"的对称退化。"""
        return self.tie_count > 1

    def as_dict(self) -> dict:
        return {
            'indices': list(self.indices), 'side': self.side,
            'layer_size': self.layer_size, 'spacing_A': round(self.spacing_A, 6),
            'tolerance_A': round(self.tolerance_A, 6),
            'centroid_gap_A': round(self.centroid_gap_A, 6),
            'area_A2': round(self.area_A2, 6),
            'longest_edge_A': round(self.longest_edge_A, 6), 'tie_count': self.tie_count,
            'fallback': self.fallback, 'candidates': self.candidates, 'rule': self.rule,
        }


def _coords(atoms, i):
    return atoms[i].coord()


def _triangle(atoms, i: int, j: int, k: int) -> tuple[float, tuple[float, float, float], float]:
    """返回 (面积, 质心, 最长边)。

    最长边用于"紧凑度": 外层是一整片三角格时, 与接触点同心的三角形可以有很多
    (实测 7 个, 质心差 < 1e-8 A), 其中既有单位三角也有横跨 10 A 的近共线三角 ——
    它们的面积读数都可能是同一量化值, 所以**必须先按最长边收窄到近邻尺度**,
    否则会选中近乎共线的三角形(aitranss 的表面旋转会退化)。
    """
    a, b, c = _coords(atoms, i), _coords(atoms, j), _coords(atoms, k)
    u = sub(b, a)
    v = sub(c, a)
    area = 0.5 * norm(cross(u, v))
    centroid = ((a[0] + b[0] + c[0]) / 3.0,
                (a[1] + b[1] + c[1]) / 3.0,
                (a[2] + b[2] + c[2]) / 3.0)
    longest = max(norm(u), norm(v), norm(sub(c, b)))
    return area, centroid, longest


def junction_frame(atoms, left_linker: int, right_linker: int) -> dict:
    """输运轴 = 左锚 -> 右锚单位向量; 返回轴、左锚坐标、轴中点(与 a2_viz 同口径)。"""
    lpos, rpos = _coords(atoms, left_linker), _coords(atoms, right_linker)
    axis = unit(sub(rpos, lpos))
    mid = tuple((lpos[t] + rpos[t]) / 2.0 for t in range(3))
    return {'axis': axis, 'left': lpos, 'right': rpos, 'mid': mid}


def electrode_members(atoms, metadata: dict, side: str) -> list[int]:
    """该侧电极原子(tip 元素, 不含对侧 tip), 用轴中点分侧 —— 不靠索引段猜。"""
    if side not in ('left', 'right'):
        raise ValueError('side 只能是 left/right, 实为 %r' % (side,))
    tip = int(metadata['%s_tip_index' % side]) - 1
    other = int(metadata['right_tip_index' if side == 'left' else 'left_tip_index']) - 1
    elem = atoms[tip].element
    fr = junction_frame(atoms, int(metadata['left_linker_index']) - 1,
                        int(metadata['right_linker_index']) - 1)
    out = []
    for i, atom in enumerate(atoms):
        if atom.element != elem or i == other:
            continue
        s = dot(sub(atom.coord(), fr['mid']), fr['axis'])   # >0 偏右
        if (side == 'left' and s < 0) or (side == 'right' and s > 0):
            out.append(i)
    return out


def surface_layer(atoms, members: list[int], away: tuple[float, float, float]) -> dict:
    """沿 away 投影最远的一层; 容差按层间距相对给出。"""
    if not members:
        raise ValueError('members 为空, 无法定外层')
    proj = {i: dot(_coords(atoms, i), away) for i in members}
    top = max(proj.values())
    diffs = sorted({round(top - p, 6) for p in proj.values() if top - p > 1e-6})
    spacing = 0.0
    if diffs:
        mid = len(diffs) // 2
        spacing = diffs[mid] if len(diffs) % 2 else 0.5 * (diffs[mid - 1] + diffs[mid])
    tol = max(MIN_LAYER_TOL, LAYER_TOL_FRAC * spacing)
    layer = sorted(i for i in members if top - proj[i] <= tol)
    return {'layer': layer, 'top': top, 'spacing': spacing, 'tolerance': tol}


def _axis_point(tip_coord, top: float, away) -> tuple[float, float, float]:
    """tip 沿 away 投影到"外层平面"上的点。"""
    t = top - dot(tip_coord, away)
    return (tip_coord[0] + t * away[0],
            tip_coord[1] + t * away[1],
            tip_coord[2] + t * away[2])


def _tie_key(atoms, i: int, j: int, k: int) -> tuple:
    """破对称平局用: 三个顶点的坐标(量化后升序), 与原子顺序无关。"""
    pts = sorted(tuple(round(v, COORD_NDIGITS) for v in _coords(atoms, x))
                 for x in (i, j, k))
    return tuple(pts)


def _score_key(longest: float, gap: float, area: float, tie: tuple) -> tuple:
    """三角排序键(升序取第一个)。

    前三分量**分桶**(宽 BUCKET_A ≈ 2% 晶格常数), 末分量是顶点坐标字典序:

    * 第 1 桶 = 最长边: 优先"紧邻"的三个原子 —— 边长 2.88 A 落在同一个桶,
      跨越倍数的组合被排到后面(避免把 `electrode_members` 里隔几层的不相邻原子配成"三角");
    * 第 2 桶 = 质心到接触点投影的距离: 优先正对接触点的位置(occ 锚的最先近邻选择);
    * 第 3 桶 = 面积;
    * 第 4 项 = 顶点坐标排序后的字典序 —— **与遍历顺序无关**, 也吸收 ~1e-2 A 的抖动
      (坐标量化到 1e-3 A)。

    分桶的意义: 对称/近对称几何里多个三角的 gap 差可以小到 1e-8(纯平局) 甚至
    1e-3(SCF 噪声), 若用 1e-6 精度比较, 选哪个由浮点噪声决定; 分桶后它们同属
    一个桶, 胜负交给坐标字典序, 结果是几何的稳定函数(见模块 docstring "为什么").
    """
    return (_bucket(longest), _bucket(gap), _bucket(area), tie)


def pick_surface_triple(atoms, members: list[int], tip_index: int,
                        side: str = '', linker_index: int | None = None,
                        away: tuple[float, float, float] | None = None) -> SurfaceTriple:
    """在 members 的最外层里选 3 个非共线原子(确定性, 见模块 docstring)。

    away: 可显式给出“离开结区”的轴(通常两侧 tip 连线); 缺省仍由 tip->linker 定向。
    分子弛豫后 S-C 键可能倾斜, 此时 tip->linker 不再等于输运轴, 审计/复算应传 away。
    """
    if away is None:
        if linker_index is None:
            raise ValueError('pick_surface_triple 需要 linker_index 定向输运轴')
        away = unit(sub(_coords(atoms, tip_index), _coords(atoms, linker_index)))
    lay = surface_layer(atoms, members, away)
    layer = lay['layer']
    q = _axis_point(_coords(atoms, tip_index), lay['top'], away)

    best = []          # [(score_key, i, j, k, gap, area, longest)]
    n = len(layer)
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(b + 1, n):
                i, j, k = layer[a], layer[b], layer[c]
                area, centroid, longest = _triangle(atoms, i, j, k)
                if area <= 1e-6:
                    continue        # 共线: aitranss rotate 会退化
                gap = norm(sub(centroid, q))
                key = _score_key(longest, gap, area, _tie_key(atoms, i, j, k))
                best.append((key, i, j, k, gap, area, longest))

    fallback = False
    if not best:
        # 外层退化(不足 3 个非共线): 按离 q 的距离取最近的 3 个非共线组合
        fallback = True
        ranked = sorted(layer, key=lambda i: (round(norm(sub(_coords(atoms, i), q)),
                                                SCORE_NDIGITS), _tie_key(atoms, i, i, i)))
        pool = ranked[:max(3, min(len(ranked), 8))]
        for a in range(len(pool)):
            for b in range(a + 1, len(pool)):
                for c in range(b + 1, len(pool)):
                    i, j, k = pool[a], pool[b], pool[c]
                    area, centroid, longest = _triangle(atoms, i, j, k)
                    if area <= 1e-6:
                        continue
                    gap = norm(sub(centroid, q))
                    best.append((_score_key(longest, gap, area,
                                            _tie_key(atoms, i, j, k)),
                                 i, j, k, gap, area, longest))
        if not best:
            raise ValueError('表面三角无法构造: members=%d layer=%d(共线?)'
                             % (len(members), len(layer)))

    best.sort(key=lambda t: t[0])
    key0, i0, j0, k0, gap0, area0, long0 = best[0]
    tie = sum(1 for t in best if t[0][:3] == key0[:3])
    return SurfaceTriple(
        indices=sorted([i0 + 1, j0 + 1, k0 + 1]), side=side,
        layer_size=len(layer), layer_indices=[i + 1 for i in layer],
        spacing_A=lay['spacing'], tolerance_A=lay['tolerance'],
        centroid_gap_A=gap0, area_A2=area0, longest_edge_A=long0, tie_count=tie,
        fallback=fallback, candidates=len(best))


def surface_indices(atoms, electrode_indices: list[int], tip_index: int,
                    linker_index: int, side: str = '') -> list[int]:
    """旧签名兼容: 只返回 1-based 索引(顺序为升序)。

    老调用方(`append_au60`)传 [tip] + 本侧新增原子; 新实现需要 linker 定向轴,
    由参数给出。需要 tie/回退等诊断时用 :func:`pick_surface_triple`。
    """
    members = sorted(set(int(i) for i in electrode_indices))
    return list(pick_surface_triple(atoms, members, tip_index, side=side,
                                    linker_index=linker_index).indices)


def verify_surface_indices(atoms, metadata: dict, side: str) -> dict:
    """重算并核对 metadata 里的三个表面索引(不抛错, 结论在 ok/reasons)。"""
    reasons = []
    given = [int(i) for i in (metadata.get('%s_surface_indices' % side) or [])]
    meta = dict(metadata)
    try:
        members = electrode_members(atoms, meta, side)
        tip = int(meta['%s_tip_index' % side]) - 1
        linker = int(meta['%s_linker_index' % side]) - 1
        got = pick_surface_triple(atoms, members, tip, side=side, linker_index=linker)
    except Exception as exc:                                   # noqa: BLE001
        return {'ok': False, 'side': side, 'given': given, 'expected': None,
                'reasons': ['重算失败: %s' % exc], 'tie_count': None}

    if len(given) != 3:
        reasons.append('metadata 表面索引不是 3 个: %r' % (given,))
    if len(set(given)) != len(given):
        reasons.append('metadata 表面索引有重复: %r' % (given,))
    for i in given:
        if not (1 <= i <= len(atoms)):
            reasons.append('索引 %d 越界(1..%d)' % (i, len(atoms)))
        elif atoms[i - 1].element != atoms[int(meta['%s_tip_index' % side]) - 1].element:
            reasons.append('索引 %d 不是电极元素(期望 %s)'
                           % (i, atoms[int(meta['%s_tip_index' % side]) - 1].element))
    if set(given) != set(got.indices):
        reasons.append('与确定性规则重算结果不一致: metadata=%s 重算=%s'
                       % (sorted(given), got.indices))
    rule = str(metadata.get('surface_rule') or 'legacy')
    return {'ok': not reasons, 'side': side, 'given': sorted(given),
            'expected': got.indices, 'rule': rule,
            'tie_count': got.tie_count, 'symmetric_tie': got.symmetric_tie,
            'layer_size': got.layer_size, 'spacing_A': round(got.spacing_A, 6),
            'fallback': got.fallback, 'reasons': reasons}
