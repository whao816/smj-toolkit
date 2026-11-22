"""接触构型枚举(W3-2)。

一条偶氮苯结的"接触"不是单一量: 锚定原子落在表面的**哪种位置**(顶位/桥位/空心位)、
分子轴相对表面法向**倾多少**、锚–表面**间距**多少, 三者共同决定 T(E_F) 量级的差异
(基准里 6 层电极的 BDT 结 T(E_F)=0.849 G₀ 与文献 0.03–0.40 G₀ 的差距, 第一顺位归因
就是接触几何)。生产上过去只用一个 `au_tip_distance` 标量 + BFGS 自由弛豫, 等于把
"接触"交给优化器碰运气 —— 本模块把它的**候选集**显式枚举出来。

枚举维度(框架 §4.2 的 V7 + §4.2.3):

* 吸附位 ``atop / bridge / hollow`` —— 由外层原子的三角格直接几何构造;
* 间距 —— ``2.20 … 2.80 Å``、步长 ``0.05 Å``(锚–最近表面原子距离);
* 倾斜 —— ``0 / 5 / 10 / 15°``(分子轴与表面法向的夹角);
* 方位 —— ``0 / 90°``(绕法向; 只取两个等价类代表, 不做 360° 密采样)。

打分只看**几何**、不看电子结构(那要 SCF, 太贵): 先剔原子重叠, 再按锚–表面距离偏离
建结目标带(``2.40–2.49 Å``, W0-1 的裁决目标; **不是** au_s_bond 判据带)排序, 最后用"非锚原子到表面的最近距离"作次级
判据(越大越不挤)。同一份几何的枚举结果**逐位可复现**: 站点、位姿、排序都有确定性键。

刻意不做的事:

* 不替代 ``aims_tip_relax`` —— 枚举出的位姿是**初猜**, 之后仍由 BFGS 弛豫;
* 不做"两侧镜像/非对称"的组合 —— 那是 junction 级别的选择(左右各自取一个位姿),
  由调用方组合, 本模块只回答"单侧有哪些候选";
* 不引入随机化 ±0.1 Å 破简并变体 —— 该变体属于"同一候选的扰动", 由调用方在写出
  geometry 时按种子生成(见 ``perturb_pose``), 以免枚举表本身变成随机量。
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field
from pathlib import Path

from .geometry import Atom, add, dot, norm, read_aims_geometry, read_xyz, scale, sub, unit
from .surface_index import electrode_members, surface_layer

#: 吸附位种类(顺序即排序权重: 顶位最常用, 空心位最"躺")
SITE_KINDS = ('atop', 'bridge', 'hollow')

#: 锚–表面**建结目标带**(Au–S, 源自 W0-1 裁决)。这是搜索/初猜排序的目标，
#: **不是** au_s_bond 判据带；判据带见 criteria.AU_S_BOND_BAND (2.40–2.44 Å)。
TARGET_BAND_A = (2.40, 2.49)
#: 判定"进带"的容差; 超出此值直接淘汰
BAND_TOLERANCE_A = 0.10
#: 非锚原子与电极原子的最近距离下限(A)。2.2 约等于 Au/S 的强排斥起点
MIN_SEPARATION_A = 2.20
#: 间距扫描
DISTANCE_SCAN_A = tuple(round(2.20 + 0.05 * i, 2) for i in range(13))   # 2.20 … 2.80
#: 倾斜角(度)
TILT_SCAN_DEG = (0.0, 5.0, 10.0, 15.0)
#: 方位角(度)。只取两个等价类代表, 不做密采样
AZIMUTH_SCAN_DEG = (0.0, 90.0)

_TOL_POINT = 1e-3          # 两个"点"被认为同一位置的容差(A)


@dataclass
class Site:
    """一个候选吸附位(由外层原子直接构造, 不含分子)。"""

    kind: str
    members: tuple[int, ...]          # 0-based 原子索引(外层)
    point: tuple[float, float, float]
    extent_A: float = 0.0             # 定义该位点的几何尺度(bridge: 键长; hollow: 最长边)
    label: str = ''

    def as_dict(self) -> dict:
        return {'kind': self.kind, 'members': [m + 1 for m in self.members],
                'point': [round(v, 4) for v in self.point],
                'extent_A': round(self.extent_A, 4), 'label': self.label}


@dataclass
class Pose:
    """分子在给定吸附位上的一个具体位姿。"""

    site: Site
    distance_A: float
    tilt_deg: float
    azimuth_A_deg: float
    atoms: list[Atom] = field(default_factory=list)     # 已摆好的分子原子
    anchor_position: tuple[float, float, float] = (0.0, 0.0, 0.0)

    def as_dict(self) -> dict:
        return {'site': self.site.as_dict(), 'distance_A': round(self.distance_A, 3),
                'tilt_deg': round(self.tilt_deg, 2),
                'azimuth_deg': round(self.azimuth_A_deg, 2),
                'anchor_position': [round(v, 4) for v in self.anchor_position]}


# ---------------------------------------------------------------------------
# 几何工具
# ---------------------------------------------------------------------------

def _dist(a, b) -> float:
    return norm(sub(a, b))


def face_normal(points) -> tuple[float, float, float]:
    """由三点定单位法向; 共线时抛 ValueError(不猜)。"""
    a, b, c = points
    n = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]), \
        (b[2] - a[2]) * (c[0] - a[0]) - (b[0] - a[0]) * (c[2] - a[2]), \
        (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    if norm(n) < 1e-8:
        raise ValueError('三点共线, 无法定法向')
    return unit(n)


def rotate_about(axis, angle_rad, vec):
    """Rodrigues 旋转: 把 vec 绕单位轴 axis 转 angle_rad。"""
    ax = unit(axis)
    c, s = math.cos(angle_rad), math.sin(angle_rad)
    return (vec[0] * c + (ax[1] * vec[2] - ax[2] * vec[1]) * s + ax[0] * dot(ax, vec) * (1 - c),
            vec[1] * c + (ax[2] * vec[0] - ax[0] * vec[2]) * s + ax[1] * dot(ax, vec) * (1 - c),
            vec[2] * c + (ax[0] * vec[1] - ax[1] * vec[0]) * s + ax[2] * dot(ax, vec) * (1 - c))


def rotation_matrix_to(from_dir, to_dir):
    """把 from_dir 转到 to_dir 的最小旋转(3x3), 用 Rodrigues 公式。"""
    a, b = unit(from_dir), unit(to_dir)
    v = (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
    c = max(-1.0, min(1.0, dot(a, b)))
    s = norm(v)
    if s < 1e-9:
        if c > 0:
            return ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
        # 反向: 任意垂直轴转 pi
        helper = (1.0, 0.0, 0.0) if abs(a[0]) < 0.9 else (0.0, 1.0, 0.0)
        axis = unit((a[1] * helper[2] - a[2] * helper[1],
                     a[2] * helper[0] - a[0] * helper[2],
                     a[0] * helper[1] - a[1] * helper[0]))
        return tuple(tuple(2.0 * axis[i] * axis[j] - (1.0 if i == j else 0.0)
                           for j in range(3)) for i in range(3))
    k = scale(v, 1.0 / s)
    kx = ((0.0, -k[2], k[1]), (k[2], 0.0, -k[0]), (-k[1], k[0], 0.0))
    kk = tuple(tuple(k[i] * k[j] for j in range(3)) for i in range(3))
    # Rodrigues: R = c*I + s*K + (1-c)*(k k^T)  —— 注意对角项是 c 而不是 1
    return tuple(tuple((c if i == j else 0.0) + s * kx[i][j] + (1.0 - c) * kk[i][j]
                       for j in range(3)) for i in range(3))


def apply_matrix(m, v):
    return (sum(m[0][j] * v[j] for j in range(3)),
            sum(m[1][j] * v[j] for j in range(3)),
            sum(m[2][j] * v[j] for j in range(3)))


# ---------------------------------------------------------------------------
# 站点构造
# ---------------------------------------------------------------------------

def nearest_neighbour_A(atoms, layer: list[int]) -> float:
    """外层点阵的最近邻距中位数(A)。单原子层没有邻居时抛错。"""
    best = []
    for i in layer:
        d = sorted(_dist(atoms[i].coord(), atoms[j].coord()) for j in layer if j != i)
        if d:
            best.append(d[0])
    if not best:
        raise ValueError('外层不足 2 个原子, 无法定最近邻距')
    best.sort()
    return best[len(best) // 2]


def _nn_pairs(atoms, layer: list[int], a_nn: float, tol_frac: float = 0.15) -> list[tuple[int, int]]:
    lo, hi = a_nn * (1.0 - tol_frac), a_nn * (1.0 + tol_frac)
    out = []
    for i, j in itertools.combinations(layer, 2):
        d = _dist(atoms[i].coord(), atoms[j].coord())
        if lo <= d <= hi:
            out.append((i, j))
    return out


def enumerate_sites(atoms: list[Atom], layer: list[int],
                    kinds: tuple[str, ...] = SITE_KINDS) -> list[Site]:
    """从外层原子直接构造 atop/bridge/hollow 站点(确定性排序, 可复现)。

    规则:

    * ``atop``  —— 每个外层原子一个;
    * ``bridge`` —— 最近邻原子对的键中点;
    * ``hollow`` —— 三个**互为最近邻**的原子构成的三角, 取质心(FCC(111) 的三重空心位)。

    同一 kind 下坐标相距 < 1e-3 Å 的站点视为同一个(取成员索引字典序最小的那组)。
    """
    if not layer:
        raise ValueError('layer 为空')
    for k in kinds:
        if k not in SITE_KINDS:
            raise ValueError('未知吸附位种类 %r, 只支持 %s' % (k, '/'.join(SITE_KINDS)))
    a_nn = nearest_neighbour_A(atoms, layer) if len(layer) > 1 else 0.0
    sites: list[Site] = []

    if 'atop' in kinds:
        for i in layer:
            p = atoms[i].coord()
            sites.append(Site(kind='atop', members=(i,), point=p,
                              extent_A=0.0,
                              label='atop:%d' % (i + 1)))

    pairs = _nn_pairs(atoms, layer, a_nn) if a_nn and 'bridge' in kinds else []
    if 'bridge' in kinds:
        for i, j in pairs:
            p = scale(add(atoms[i].coord(), atoms[j].coord()), 0.5)
            sites.append(Site(kind='bridge', members=(i, j), point=p,
                              extent_A=a_nn, label='bridge:%d-%d' % (i + 1, j + 1)))

    if 'hollow' in kinds and pairs:
        pair_set = {tuple(sorted(p)) for p in pairs}
        for tri in itertools.combinations(layer, 3):
            if all(tuple(sorted((a, b))) in pair_set for a, b in itertools.combinations(tri, 2)):
                pts = [atoms[i].coord() for i in tri]
                p = scale((sum(q[0] for q in pts), sum(q[1] for q in pts),
                           sum(q[2] for q in pts)), 1.0 / 3.0)
                longest = max(_dist(pts[0], pts[1]), _dist(pts[1], pts[2]), _dist(pts[0], pts[2]))
                sites.append(Site(kind='hollow', members=tuple(sorted(tri)), point=p,
                                  extent_A=longest,
                                  label='hollow:%d-%d-%d' % tuple(sorted(x + 1 for x in tri))))

    # 去重 + 确定性排序: kind 权重 -> 坐标(量化) -> 成员索引
    uniq: list[Site] = []
    seen: set[tuple] = set()
    for s in sorted(sites, key=lambda x: (SITE_KINDS.index(x.kind),
                                          round(x.point[0], 3), round(x.point[1], 3),
                                          round(x.point[2], 3), x.members)):
        key = (s.kind, round(s.point[0], 3), round(s.point[1], 3), round(s.point[2], 3))
        if key in seen:
            continue
        seen.add(key)
        uniq.append(s)
    # 二次合并: 同 kind 且相距 < _TOL_POINT 的(量化边界可能拆成两条)
    merged: list[Site] = []
    for s in uniq:
        if any(m.kind == s.kind and _dist(m.point, s.point) < _TOL_POINT for m in merged):
            continue
        merged.append(s)
    return merged


# ---------------------------------------------------------------------------
# 位姿
# ---------------------------------------------------------------------------

def place_anchor(molecule: list[Atom], anchor_index: int, site: Site,
                 into_molecule: tuple[float, float, float], distance_A: float,
                 tilt_deg: float = 0.0, azimuth_deg: float = 0.0) -> Pose:
    """把分子摆到站点上: 锚原子落在 ``site.point + into*distance``。

    分子自身的参考方向取"锚原子 -> 分子质心"; 先把它转到 ``into_molecule``(tilt=0),
    再绕表面法向转 ``azimuth``, 最后按 ``tilt`` 由法向偏出。**不翻转手性**: 用的是
    最小旋转, 不会引入镜像。
    """
    if not 0 <= anchor_index < len(molecule):
        raise ValueError('anchor_index 越界: %d' % anchor_index)
    if distance_A <= 0:
        raise ValueError('distance_A 必须为正')
    anchor = molecule[anchor_index].coord()
    rest = [a for k, a in enumerate(molecule) if k != anchor_index]
    if not rest:
        raise ValueError('分子只有锚原子, 无法定向')
    centroid = (sum(a.coord()[0] for a in rest) / len(rest) + anchor[0] / len(rest),
                sum(a.coord()[1] for a in rest) / len(rest) + anchor[1] / len(rest),
                sum(a.coord()[2] for a in rest) / len(rest) + anchor[2] / len(rest))
    axis = unit(sub(centroid, anchor))
    n = unit(into_molecule)
    # 目标方向: 从法向出发, 绕"垂直于法向且在方位角方向"的轴偏 tilt
    helper = (0.0, 0.0, 1.0) if abs(n[2]) < 0.9 else (1.0, 0.0, 0.0)
    e1 = unit((n[1] * helper[2] - n[2] * helper[1],
               n[2] * helper[0] - n[0] * helper[2],
               n[0] * helper[1] - n[1] * helper[0]))
    e2 = unit((n[1] * e1[2] - n[2] * e1[1],
               n[2] * e1[0] - n[0] * e1[2],
               n[0] * e1[1] - n[1] * e1[0]))
    az = math.radians(azimuth_deg)
    inplane = (e1[0] * math.cos(az) + e2[0] * math.sin(az),
               e1[1] * math.cos(az) + e2[1] * math.sin(az),
               e1[2] * math.cos(az) + e2[2] * math.sin(az))
    tilt = math.radians(tilt_deg)
    target = (n[0] * math.cos(tilt) + inplane[0] * math.sin(tilt),
              n[1] * math.cos(tilt) + inplane[1] * math.sin(tilt),
              n[2] * math.cos(tilt) + inplane[2] * math.sin(tilt))
    rot = rotation_matrix_to(axis, unit(target))
    anchor_pos = add(site.point, scale(n, distance_A))
    placed = []
    for a in molecule:
        rel = sub(a.coord(), anchor)
        r = apply_matrix(rot, rel)
        placed.append(Atom(element=a.element, x=anchor_pos[0] + r[0],
                           y=anchor_pos[1] + r[1], z=anchor_pos[2] + r[2]))
    return Pose(site=site, distance_A=distance_A, tilt_deg=tilt_deg,
                azimuth_A_deg=azimuth_deg, atoms=placed, anchor_position=anchor_pos)


def perturb_pose(pose: Pose, seed: int, amplitude_A: float = 0.1,
                 skip: tuple[int, ...] = ()) -> Pose:
    """±amplitude 的确定性抖动(框架 §4.2.3 的破简并变体)。

    用 ``random.Random(seed)`` 逐原子各三个分量; 只抖**非锚**原子, 锚与冻结原子保持
    原位, 否则"接触距离"这个自变量就被抖掉了。
    """
    import random
    rng = random.Random(seed)
    anchor_pos = pose.anchor_position
    out = []
    for k, a in enumerate(pose.atoms):
        # 锚原子永远不抖(它是接触点本身), skip 供调用方额外指定冻结原子
        if k in skip or a.coord() == anchor_pos:
            out.append(Atom(element=a.element, x=a.x, y=a.y, z=a.z))
            continue
        out.append(Atom(element=a.element,
                        x=a.x + rng.uniform(-amplitude_A, amplitude_A),
                        y=a.y + rng.uniform(-amplitude_A, amplitude_A),
                        z=a.z + rng.uniform(-amplitude_A, amplitude_A)))
    return Pose(site=pose.site, distance_A=pose.distance_A, tilt_deg=pose.tilt_deg,
                azimuth_A_deg=pose.azimuth_A_deg, atoms=out,
                anchor_position=pose.anchor_position)


# ---------------------------------------------------------------------------
# 打分与枚举
# ---------------------------------------------------------------------------

def score_contact(electrode: list[Atom], pose: Pose,
                  target_band: tuple[float, float] = TARGET_BAND_A,
                  min_separation: float = MIN_SEPARATION_A) -> dict:
    """纯几何打分: 锚–表面距离是否进带 + 非锚原子到表面的最近距离。"""
    anchor_pos = pose.anchor_position
    d_anchor = min(_dist(anchor_pos, e.coord()) for e in electrode)
    min_sep = float('inf')
    clashes = 0
    for k, a in enumerate(pose.atoms):
        if a.coord() == anchor_pos:
            continue
        d = min(_dist(a.coord(), e.coord()) for e in electrode)
        if k and d < min_sep:
            min_sep = d
        if d < min_separation:
            clashes += 1
    lo, hi = target_band
    dev = 0.0 if lo <= d_anchor <= hi else min(abs(d_anchor - lo), abs(d_anchor - hi))
    return {'au_s_A': d_anchor, 'dev_A': dev, 'min_separation_A': min_sep,
            'clashes': clashes,
            'in_band': lo <= d_anchor <= hi,
            'bond_ok': dev <= BAND_TOLERANCE_A}


def enumerate_contacts(atoms: list[Atom], metadata: dict, molecule: list[Atom],
                       anchor_index: int, side: str = 'left',
                       kinds: tuple[str, ...] = SITE_KINDS,
                       distances: tuple[float, ...] = DISTANCE_SCAN_A,
                       tilts: tuple[float, ...] = TILT_SCAN_DEG,
                       azimuths: tuple[float, ...] = AZIMUTH_SCAN_DEG,
                       top_k: int = 3,
                       min_separation: float = MIN_SEPARATION_A) -> list[dict]:
    """枚举一侧的接触候选并按几何特征排序(确定性)。

    淘汰规则(硬): 有原子重叠(``clashes > 0``); 锚–表面距离偏离目标带超过
    ``BAND_TOLERANCE_A``。
    排序键: ``(dev_A, -min_separation_A, kind 权重, 坐标, 倾斜, 间距)`` —— 距带内越近
    越优先, 同等接近时"更不挤"的优先, 其余按坐标字典序破平局。
    """
    tip_index = int(metadata['%s_tip_index' % side]) - 1
    tip = atoms[tip_index]
    linker = atoms[int(metadata['%s_linker_index' % side]) - 1]
    into_molecule = unit(sub(linker.coord(), tip.coord()))
    # 只按**电极表面**枚举: tip 由 tip_relax 阶段沿锚–表面连线插入, 不参与本步打分,
    # 否则"锚–最近 Au 距离"会被已有的 tip 抢占, 得到的 d_elec 就不再是同一个量。
    members = [i for i in electrode_members(atoms, metadata, side) if i != tip_index]
    layer = surface_layer(atoms, members, into_molecule)['layer']
    sites = enumerate_sites(atoms, layer, kinds)
    electrode = [atoms[i] for i in members]
    out: list[dict] = []
    for site in sites:
        for d in distances:
            for tilt in tilts:
                for az in azimuths:
                    pose = place_anchor(molecule, anchor_index, site, into_molecule,
                                        d, tilt, az)
                    sc = score_contact(electrode, pose, min_separation=min_separation)
                    if sc['clashes'] or not sc['bond_ok']:
                        continue
                    rec = {'side': side, 'anchor_index': anchor_index + 1,
                           'pose': pose.as_dict(), 'score': {k: round(v, 4)
                                                             if isinstance(v, float) else v
                                                             for k, v in sc.items()},
                           'sort_key': [round(sc['dev_A'], 4), -round(sc['min_separation_A'], 4),
                                        SITE_KINDS.index(site.kind),
                                        round(site.point[0], 3), round(site.point[1], 3),
                                        round(site.point[2], 3), tilt, d, az]}
                    out.append(rec)
    out.sort(key=lambda r: r['sort_key'])
    return out[:top_k]


def load_molecule(path: str | Path) -> list[Atom]:
    """读分子几何(``.xyz`` 或 FHI-aims ``geometry.in``, 按扩展名判, 不猜)。"""
    p = Path(path)
    if p.suffix.lower() == '.xyz':
        return read_xyz(p)
    if p.name.endswith('geometry.in') or p.suffix == '.in':
        return read_aims_geometry(p)
    raise ValueError('无法判定的分子几何格式: %s' % p)
