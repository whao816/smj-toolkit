"""单分子结几何工具:Atom、XYZ/FHI-aims geometry.in 读写、向量代数、Au tip 添加。

来源: 第三方工具 qmtransport/geometry.py(134 行, tmp/自动化/, 只读来源不删除不入库)
按本项目协议适配迁入(AGD-2026-09-22-SMJ, P2):保留全部公式与语义;承载层改为
本仓库头注/命名/门控约定。纯本地工具, 不做真实计算。
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math
import re
from typing import Iterable


@dataclass
class Atom:
    element: str
    x: float
    y: float
    z: float

    def coord(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)


_NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
_XYZ_RE = re.compile(rf"^\s*([A-Za-z][A-Za-z]?)\d*\s+({_NUM})\s+({_NUM})\s+({_NUM})\s*$")
_AIMS_RE = re.compile(rf"^\s*atom\s+({_NUM})\s+({_NUM})\s+({_NUM})\s+([A-Za-z][A-Za-z]?)\s*$")


def read_xyz(path: str | Path) -> list[Atom]:
    atoms: list[Atom] = []
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    start = 0
    if lines and lines[0].strip().isdigit():
        start = 2
    for line in lines[start:]:
        match = _XYZ_RE.match(line)
        if match:
            element, x, y, z = match.groups()
            atoms.append(Atom(element, float(x), float(y), float(z)))
    if not atoms:
        raise ValueError(f"No XYZ atoms found in {path}")
    return atoms


def read_aims_geometry(path: str | Path) -> list[Atom]:
    atoms: list[Atom] = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        match = _AIMS_RE.match(line)
        if match:
            x, y, z, element = match.groups()
            atoms.append(Atom(element, float(x), float(y), float(z)))
    if not atoms:
        raise ValueError(f"No FHI-aims atoms found in {path}")
    return atoms


def write_xyz(path: str | Path, atoms: Iterable[Atom], comment: str = "") -> None:
    atoms = list(atoms)
    text = [str(len(atoms)), comment]
    text.extend(f"{a.element:2s} {a.x:16.8f} {a.y:16.8f} {a.z:16.8f}" for a in atoms)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("\n".join(text) + "\n", encoding="utf-8", newline="\n")


def write_aims_geometry(path: str | Path, atoms: Iterable[Atom],
                        constrained: Iterable[int] = ()) -> None:
    """写 FHI-aims geometry.in。

    `constrained` = 需要冻结的原子序号(**1-based**, 与 left_anchor/right_anchor 同一计数),
    对应 atom 行之后追加一行 `constrain_relaxation .true.`(FHI-aims 语法)。
    缺省空集合 = 输出与旧版逐字节一致(W0-1 兼容约束: 既有 25 条 smj 记录的锚点回归靠它)。
    """
    atoms = list(atoms)
    fixed = {int(i) for i in constrained}
    if fixed and (min(fixed) < 1 or max(fixed) > len(atoms)):
        raise ValueError('constrained 原子序号越界(1..%d): %r' % (len(atoms), sorted(fixed)))
    text: list[str] = []
    for idx, a in enumerate(atoms, start=1):
        text.append(f"atom  {a.x:16.8f}  {a.y:16.8f}  {a.z:16.8f}  {a.element}")
        if idx in fixed:
            text.append("constrain_relaxation .true.")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("\n".join(text) + "\n", encoding="utf-8", newline="\n")


def elements_in_geometry(atoms: Iterable[Atom]) -> list[str]:
    return sorted({a.element for a in atoms}, key=lambda x: (x != "H", x))


def sub(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def add(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def scale(a: tuple[float, float, float], s: float) -> tuple[float, float, float]:
    return (a[0] * s, a[1] * s, a[2] * s)


def dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def norm(a: tuple[float, float, float]) -> float:
    return math.sqrt(dot(a, a))


def unit(a: tuple[float, float, float]) -> tuple[float, float, float]:
    n = norm(a)
    if n == 0:
        raise ValueError("Zero-length vector cannot be normalized")
    return scale(a, 1.0 / n)


def parse_vector(text: str) -> tuple[float, float, float]:
    parts = [float(x.strip()) for x in text.split(",")]
    if len(parts) != 3:
        raise ValueError("Vector must have form x,y,z")
    return (parts[0], parts[1], parts[2])


def add_tip_au(
    atoms: list[Atom],
    left_anchor: int,
    right_anchor: int,
    distance: float,
    left_direction: tuple[float, float, float] | None = None,
    right_direction: tuple[float, float, float] | None = None,
    element: str = "Au",
) -> list[Atom]:
    """在两个锚定原子外侧各放一个电极金属 tip。

    ``element`` 默认 Au（旧行为不变）；构建 Ag/Cu/… 电极时由调用方给出，
    tip 元素必须与电极金属一致，否则后续 append_au60 找不到成对的 tip 原子。
    """
    left = atoms[left_anchor - 1]
    right = atoms[right_anchor - 1]
    if left_direction is None:
        left_direction = sub(left.coord(), right.coord())
    if right_direction is None:
        right_direction = sub(right.coord(), left.coord())
    left_pos = add(left.coord(), scale(unit(left_direction), distance))
    right_pos = add(right.coord(), scale(unit(right_direction), distance))
    return atoms + [Atom(element, *left_pos), Atom(element, *right_pos)]


def strip_thiol_h(atoms: list[Atom], anchors: tuple[int, int] | None = None,
                  sh_cutoff: float = 1.70,
                  element: str | None = None) -> tuple[list[Atom], dict]:
    """移除锚定原子上的 H(锚-H < sh_cutoff A) —— "中性优化、去质子建结"口径。

    历史路径(二硫醇盐, docs/taskcard-smj-bdt-benchmark.md §1.4.1): Gaussian 优化中性
    HS-R-SH, 建结时移除两个 S-H 的 H, 让两端 S 直接与电极 tip 成键。本函数已推广到
    任意锚定元素: 分子里没有 S 时(N/Se/…锚定), 由 ``anchors`` 显式给两个 1-based
    索引(或由 ``element`` 指定元素), 只移除这两个锚原子上的 H。

    只做"移除"这一件事, 归因清楚: 是否该去由调用方(``strip_thiol_h`` 配置)显式决定;
    化学歧义一律 fail-closed(不静默挑原子) —— 每个锚原子必须恰好 0 或 1 个 H,
    总共 0 或 2 个才移除; 例如 -NH2 带 2 个 H 会直接报错, 不会猜着只去一个。

    ``anchors=(left, right)`` 为原始 1-based 索引时, 返回的 info 里给出重映射后的索引
(去 H 后索引前移, 不静默沿用旧索引)。返回 (kept_atoms, info)。
    """
    if anchors and anchors[0] and anchors[1]:
        left_in, right_in = int(anchors[0]), int(anchors[1])
        for label, idx in (("left", left_in), ("right", right_in)):
            if not 1 <= idx <= len(atoms):
                raise ValueError("strip_thiol_h: %s 锚定索引 %d 越界(分子 %d 原子)"
                                 % (label, idx, len(atoms)))
        anchor_idx = [left_in - 1, right_in - 1]
    elif element is not None:
        anchor_idx = [i for i, a in enumerate(atoms) if a.element == element]
        if len(anchor_idx) != 2:
            raise ValueError("strip_thiol_h 需恰好 2 个 %s 锚原子, 实得 %d"
                             % (element, len(anchor_idx)))
    else:
        anchor_idx = [i for i, a in enumerate(atoms) if a.element == "S"]
        if len(anchor_idx) != 2:
            raise ValueError("strip_thiol_h 需恰好 2 个 S 锚原子, 实得 %d(无 S 分子请显式给出 anchors)"
                             % len(anchor_idx))
    anchor_element = atoms[anchor_idx[0]].element
    hydrogen_anchors = [i + 1 for i in anchor_idx if atoms[i].element == "H"]
    if hydrogen_anchors:
        raise ValueError("strip_thiol_h: 锚定原子不能是 H(%r)" % hydrogen_anchors)
    removed: list[int] = []
    for j in anchor_idx:
        hs = [i for i, a in enumerate(atoms)
              if a.element == "H" and math.dist(a.coord(), atoms[j].coord()) < sh_cutoff]
        if len(hs) > 1:
            raise ValueError("strip_thiol_h: %s 锚定原子(第 %d 个)带 %d 个 H, 去质子会改变"
                             "化学态; 请显式建模或设 strip_thiol_h: false"
                             % (anchor_element, j + 1, len(hs)))
        removed.extend(hs)
    if len(removed) not in (0, 2):
        raise ValueError("strip_thiol_h 期望移除 0 或 2 个 %s-H 氢(%s-H < %.2f A), "
                         "实得 %d 个: %r" % (anchor_element, anchor_element, sh_cutoff,
                                            len(removed), [i + 1 for i in removed]))
    drop = set(removed)
    kept = [a for i, a in enumerate(atoms) if i not in drop]
    remap: dict[int, int] = {}
    shift = 0
    for i in range(len(atoms)):
        if i in drop:
            shift += 1
        else:
            remap[i + 1] = i + 1 - shift
    info = {"removed": len(removed),
            "removed_h_indices": [i + 1 for i in removed],
            "sh_cutoff_A": sh_cutoff,
            "anchor_element": anchor_element,
            "anchor_indices": [remap[i + 1] for i in anchor_idx],
            "atom_count": len(kept),
            "route": ("thiolate" if anchor_element == "S" else "anchor_h_removed")
                     if removed else "no_anchor_h"}
    if anchor_element == "S":
        info["sulfur_indices"] = [remap[i + 1] for i in anchor_idx]
    if anchors and anchors[0] and anchors[1]:
        info["left_anchor_in"] = left_in
        info["right_anchor_in"] = right_in
        info["left_anchor"] = remap[left_in]
        info["right_anchor"] = remap[right_in]
    return kept, info


def count_atoms(path: str | Path) -> int:
    return len(read_aims_geometry(path))