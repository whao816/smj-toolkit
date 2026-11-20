"""阶段产物几何的合理性检查（纯 python，无第三方依赖）。

``molecule_opt`` 原来只核"进程正常结束 + 文件存在"，键断裂、质子转移、原子重叠
都能一路进链，而下游 ``aims_tip_relax`` 会把它当"已优化的分子"用。这里补一组
廉价检查：

* 原子重叠：最近原子对距离 < 0.7 Å；
* 键图对比（输入几何 vs 输出几何）：重原子键消失 = 断键/解离（拒绝），
  含 H 的键消失 = 警告（质子可能迁移，人工确认但不阻断）；
* 孤立原子：整枚原子没有任何一根键（拒绝）；
* 位移：平移对齐后的 RMSD（只报告，不作为判据——优化本来就会挪动原子）。

键判定用共价半径和 ``d < bond_factor * (r_i + r_j)``。未知元素不参与键图
（``bond_check_skipped=True``），只做重叠与有限值检查——宁可不判，也不乱判。

约定：``problems`` 非空即 ``ok=False``，调用方应 raise；``warnings`` 写进 notes。
"""

from __future__ import annotations

import math
from itertools import combinations
from typing import Any, Iterable, Sequence

from .junction.geometry import Atom
from .junction.molecules import rmsd as _rmsd

__all__ = ["ATOMIC_NUMBERS", "COVALENT_RADII", "anchor_contacts", "check_geometry",
           "min_atom_distance", "total_electrons"]

#: 原子序数（只覆盖本仓库会用到的元素；未知元素让调用方显式报错，不猜）。
ATOMIC_NUMBERS: dict[str, int] = {
    "H": 1, "Li": 3, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Na": 11,
    "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "K": 19,
    "Ca": 20, "Ti": 22, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29,
    "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Ag": 47,
    "Cd": 48, "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I": 53, "Pt": 78,
    "Au": 79,
}

#: 共价半径（Å），Cordero 2008 量级；只用于"两个原子是否成键"的粗判。
COVALENT_RADII: dict[str, float] = {
    "H": 0.31, "Li": 1.28, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66,
    "F": 0.57, "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07,
    "S": 1.05, "Cl": 1.02, "K": 2.03, "Ca": 1.76, "Ti": 1.60, "Mn": 1.39,
    "Fe": 1.32, "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Ga": 1.22,
    "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20, "Ag": 1.45, "Cd": 1.44,
    "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "Pt": 1.36,
    "Au": 1.36,
}

#: 最近原子对距离低于该值判为重叠/塌缩（H2 平衡键长 0.74 Å 仍可通过）。
DEFAULT_MIN_DISTANCE = 0.7
#: 键判定阈值系数：d < factor * (r_i + r_j)。
DEFAULT_BOND_FACTOR = 1.3
#: 锚定原子–金属接触判定上限（Au–S 常见 2.35–2.55 Å，放到 3.2 Å 兜住弛豫后拉长）。
DEFAULT_CONTACT_DISTANCE = 3.2
#: 参与接触报告的锚定元素（与 electrode.left_anchor/right_anchor 的默认口径一致）。
DEFAULT_ANCHOR_ELEMENTS: tuple[str, ...] = ("S", "N", "Se")


def _distance(a: Atom, b: Atom) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def min_atom_distance(atoms: Sequence[Atom]) -> tuple[float | None, tuple[int, int] | None]:
    """最近原子对的距离与（0 基）下标；原子数 < 2 时返回 ``(None, None)``。"""
    best: float | None = None
    pair: tuple[int, int] | None = None
    for i, j in combinations(range(len(atoms)), 2):
        d = _distance(atoms[i], atoms[j])
        if best is None or d < best:
            best, pair = d, (i, j)
    return best, pair


def _bond_graph(
    atoms: Sequence[Atom], factor: float
) -> tuple[dict[tuple[int, int], tuple[str, str, float]], set[str]]:
    """共价半径判键。返回 ``({(i, j): (elem_i, elem_j, d)}, 未知元素集合)``。"""
    graph: dict[tuple[int, int], tuple[str, str, float]] = {}
    unknown: set[str] = set()
    for i, j in combinations(range(len(atoms)), 2):
        ei, ej = atoms[i].element, atoms[j].element
        ri, rj = COVALENT_RADII.get(ei), COVALENT_RADII.get(ej)
        if ri is None or rj is None:
            unknown.update(el for el, r in ((ei, ri), (ej, rj)) if r is None)
            continue
        d = _distance(atoms[i], atoms[j])
        if d < factor * (ri + rj):
            graph[(i, j)] = (ei, ej, d)
    return graph, unknown


def _label(atoms: Sequence[Atom], key: tuple[int, int]) -> str:
    i, j = key
    return "%s%d-%s%d" % (atoms[i].element, i + 1, atoms[j].element, j + 1)


def _is_heavy(atoms: Sequence[Atom], key: tuple[int, int]) -> bool:
    return all(atoms[k].element != "H" for k in key)


def _angle_deg(a: Atom, b: Atom, c: Atom) -> float:
    """顶点 ``b`` 处的 a–b–c 夹角（度）；退化时返回 NaN。"""
    v1 = (a.x - b.x, a.y - b.y, a.z - b.z)
    v2 = (c.x - b.x, c.y - b.y, c.z - b.z)
    n1 = math.sqrt(sum(v * v for v in v1))
    n2 = math.sqrt(sum(v * v for v in v2))
    if n1 == 0.0 or n2 == 0.0:
        return float("nan")
    cos = sum(x * y for x, y in zip(v1, v2)) / (n1 * n2)
    return math.degrees(math.acos(max(-1.0, min(1.0, cos))))


def anchor_contacts(
    atoms: Iterable[Atom],
    *,
    metal_element: str = "Au",
    anchor_elements: Sequence[str] = DEFAULT_ANCHOR_ELEMENTS,
    max_distance: float = DEFAULT_CONTACT_DISTANCE,
) -> list[dict[str, Any]]:
    """贴着金属的锚定原子的接触几何（M–X 距离 + C–X–M 夹角）。

    只报告、不判定：合格区间取决于锚定化学（硫醇盐/硫醚）与表面位点，本函数
    给出原始量（``metal_distance_A``、``angle_deg``）供 review/report 使用。
    ``neighbor`` 取锚定原子最近的"非金属、非 H"原子，预期是与它成键的 C。
    """
    atoms = list(atoms)
    metal = [i for i, atom in enumerate(atoms) if atom.element == metal_element]
    contacts: list[dict[str, Any]] = []
    for idx, atom in enumerate(atoms):
        if atom.element not in anchor_elements:
            continue
        nearest: tuple[int, float] | None = None
        for j in metal:
            d = _distance(atom, atoms[j])
            if nearest is None or d < nearest[1]:
                nearest = (j, d)
        if nearest is None or nearest[1] > max_distance:
            continue
        neighbor: tuple[int, float] | None = None
        for k, other in enumerate(atoms):
            if k == idx or other.element in (metal_element, "H"):
                continue
            d = _distance(atom, other)
            if neighbor is None or d < neighbor[1]:
                neighbor = (k, d)
        angle = None
        if neighbor is not None:
            angle = _angle_deg(atoms[neighbor[0]], atom, atoms[nearest[0]])
        contacts.append({
            "anchor": idx + 1,
            "anchor_element": atom.element,
            "metal": nearest[0] + 1,
            "metal_element": atoms[nearest[0]].element,
            "metal_distance_A": round(nearest[1], 4),
            "neighbor": neighbor[0] + 1 if neighbor is not None else None,
            "neighbor_element": atoms[neighbor[0]].element if neighbor is not None else None,
            "neighbor_distance_A": round(neighbor[1], 4) if neighbor is not None else None,
            "angle_deg": round(angle, 2) if angle is not None else None,
        })
    return contacts


def check_geometry(
    before: Iterable[Atom] | None,
    after: Iterable[Atom],
    *,
    min_distance: float = DEFAULT_MIN_DISTANCE,
    bond_factor: float = DEFAULT_BOND_FACTOR,
) -> dict[str, Any]:
    """比较输入/输出几何，给出 ``ok`` 与逐条问题。

    Args:
        before: 优化前的几何（可 ``None``：只做绝对检查，如 ``geometry.in`` 缺失）；
        after: 优化后的几何；
        min_distance: 原子重叠阈值（Å）；
        bond_factor: 键判定系数。

    Returns:
        dict，键：``ok`` / ``atom_count`` / ``problems`` / ``warnings`` /
        ``min_distance_A`` / ``bond_count`` / ``lost_bonds`` / ``new_bonds`` /
        ``isolated_atoms`` / ``rmsd_A`` / ``bond_check_skipped``。
    """
    atoms = list(after)
    prev = list(before) if before is not None else None
    problems: list[str] = []
    warnings: list[str] = []

    finite = True
    for idx, atom in enumerate(atoms):
        if not all(math.isfinite(v) for v in (atom.x, atom.y, atom.z)):
            finite = False
            problems.append("第 %d 个原子(%s)坐标含非有限值" % (idx + 1, atom.element))

    dist: float | None = None
    if len(atoms) > 1:
        dist, pair = min_atom_distance(atoms)
        if dist is not None and pair is not None and dist < min_distance:
            problems.append(
                "原子重叠：%s 与 %s 相距 %.3f Å（阈值 %.2f Å）"
                % (atoms[pair[0]].element, atoms[pair[1]].element, dist, min_distance))

    after_graph, unknown = _bond_graph(atoms, bond_factor) if finite else ({}, set())
    skipped = bool(unknown)
    if skipped:
        warnings.append("未知共价半径的元素 %s：跳过键图对比（只做重叠/有限值检查）"
                        % "/".join(sorted(unknown)))

    degree: dict[int, int] = {}
    for i, j in after_graph:
        degree[i] = degree.get(i, 0) + 1
        degree[j] = degree.get(j, 0) + 1
    isolated: list[str] = []
    if len(atoms) > 1 and not skipped and finite:
        for idx, atom in enumerate(atoms):
            if degree.get(idx, 0) == 0:
                isolated.append("%s%d" % (atom.element, idx + 1))
        if isolated:
            problems.append("孤立原子（无任何键）：%s；几何可能已解离"
                            % ", ".join(isolated))

    lost: list[str] = []
    new: list[str] = []
    bond_count: list[int] | None = None
    rmsd_value: float | None = None
    if prev is not None:
        if len(prev) != len(atoms):
            problems.append("原子数变化：%d → %d" % (len(prev), len(atoms)))
        elif [a.element for a in prev] != [a.element for a in atoms]:
            if sorted(a.element for a in prev) == sorted(a.element for a in atoms):
                warnings.append("元素顺序与输入不一致（如 Gaussian Standard orientation "
                                "重排）：跳过键图对比")
            else:
                problems.append("元素组成变化：与输入不是同一个分子")
        elif not skipped and finite:
            before_graph, _ = _bond_graph(prev, bond_factor)
            bond_count = [len(before_graph), len(after_graph)]
            for key in sorted(set(before_graph) - set(after_graph)):
                lost.append(_label(prev, key))
                if _is_heavy(prev, key):
                    problems.append("重原子键消失：%s（%.3f Å）——断键/解离"
                                    % (_label(prev, key), before_graph[key][2]))
                else:
                    warnings.append("含 H 的键消失：%s（%.3f Å）——质子可能迁移，请人工确认"
                                    % (_label(prev, key), before_graph[key][2]))
            for key in sorted(set(after_graph) - set(before_graph)):
                new.append(_label(atoms, key))
                if _is_heavy(atoms, key):
                    warnings.append("新增重原子键：%s（%.3f Å）——请确认不是成环/聚合"
                                    % (_label(atoms, key), after_graph[key][2]))
        if len(prev) == len(atoms) and [a.element for a in prev] == [a.element for a in atoms]:
            try:
                rmsd_value = float(_rmsd(prev, atoms))
            except ValueError:  # pragma: no cover - 长度已校验，理论上不会发生
                rmsd_value = None

    return {
        "ok": not problems,
        "atom_count": len(atoms),
        "problems": problems,
        "warnings": warnings,
        "min_distance_A": None if dist is None else round(dist, 4),
        "bond_count": bond_count,
        "lost_bonds": lost,
        "new_bonds": new,
        "isolated_atoms": isolated,
        "rmsd_A": None if rmsd_value is None else round(rmsd_value, 4),
        "bond_check_skipped": skipped,
    }


def total_electrons(atoms, charge: int = 0) -> int:
    """原子序数之和减电荷。未知元素抛 KeyError —— 调用方应显式报错而非猜。"""
    total = 0
    for atom in atoms:
        total += ATOMIC_NUMBERS[atom.element.capitalize()]
    return total - int(charge)