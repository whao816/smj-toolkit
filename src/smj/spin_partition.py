"""Hirshfeld 自旋矩解析与分区（2026-11-21 自由基支线 S1）。

FHI-aims 在 ``output hirshfeld`` 下会把每个原子的 Hirshfeld 布居打印进
``aims.out``（``Performing Hirshfeld analysis ...`` 段），其中
``Hirshfeld spin moment`` 就是该原子分到的自旋磁矩（μB）。总磁矩守恒不能
说明"自由基还在分子上"——必须看分子/电极各自的份额。

这里只做**机械解析与分区**，不做任何物理判定：

* 取**最后一段**分析（前面的是迭代中间态，不能当终态）；
* 缺 ``Hirshfeld spin moment`` 行的原子记 ``None``，不冒充 0；
* 分区由调用方给（分子 = 前 N 个原子；接触区 = 距锚原子 ≤cutoff 的电极原子）；
* 比值/占比是派生量，供报告与复核陈述使用；判据不使用未校准的阈值。
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

__all__ = [
    "AtomMoment",
    "HirshfeldTable",
    "contact_atoms",
    "default_groups",
    "parse_hirshfeld",
    "parse_hirshfeld_text",
    "read_geometry",
    "summarize",
]

_ATOM_RE = re.compile(r"^\s*\|\s*Atom\s+(\d+):\s*(\S+)\s*$")
_SPIN_RE = re.compile(
    r"^\s*\|?\s*Hirshfeld spin moment\s*:\s*"
    r"(-?\d+(?:\.\d*)?(?:[eE][+-]?\d+)?)\s*$")
_SECTION_RE = re.compile(r"Performing Hirshfeld analysis", re.IGNORECASE)


@dataclass(frozen=True)
class AtomMoment:
    """单个原子的 Hirshfeld 自旋矩（μB）；缺字段时 ``spin_moment is None``。"""

    index: int
    element: str
    spin_moment: float | None


@dataclass
class HirshfeldTable:
    """一次 ``aims.out`` 里最后一段 Hirshfeld 分析的解析结果。"""

    atoms: list[AtomMoment] = field(default_factory=list)
    analyses: int = 0

    @property
    def total(self) -> float | None:
        """参与分析的原子的自旋矩之和；一个可用值都没有时返回 ``None``。"""
        values = [a.spin_moment for a in self.atoms if a.spin_moment is not None]
        return sum(values) if values else None

    @property
    def missing(self) -> int:
        return sum(1 for a in self.atoms if a.spin_moment is None)

    def as_dict(self) -> dict:
        return {
            "analyses": self.analyses,
            "atoms": len(self.atoms),
            "missing_moments": self.missing,
            "total": self.total,
            "atom_moments": [
                {"index": a.index, "element": a.element, "spin_moment": a.spin_moment}
                for a in self.atoms
            ],
        }


def parse_hirshfeld_text(text: str) -> HirshfeldTable:
    """解析 ``aims.out`` 文本，只保留最后一段 Hirshfeld 分析。"""
    table = HirshfeldTable()
    current: list[AtomMoment] = []
    for line in text.splitlines():
        if _SECTION_RE.search(line):
            table.analyses += 1
            current = []
            continue
        atom = _ATOM_RE.match(line)
        if atom is not None:
            current.append(AtomMoment(int(atom.group(1)), atom.group(2), None))
            continue
        moment = _SPIN_RE.match(line)
        if moment is not None and current:
            value = float(moment.group(1))
            current[-1] = AtomMoment(current[-1].index, current[-1].element, value)
    # 循环结束时 ``current`` 就是最后一段（空段也保留 analyses 计数）。
    table.atoms = current
    return table


def parse_hirshfeld(path: str | Path) -> HirshfeldTable:
    """读 ``aims.out`` 并解析最后一段 Hirshfeld 分析。"""
    return parse_hirshfeld_text(
        Path(path).read_text(encoding="utf-8", errors="replace"))


def read_geometry(path: str | Path) -> list[tuple[str, float, float, float]]:
    """``geometry.in`` -> ``[(element, x, y, z)]``（1-based 顺序与文件一致）。"""
    atoms: list[tuple[str, float, float, float]] = []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split()
        if len(parts) < 5 or parts[0] != "atom":
            continue
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            continue
        atoms.append((parts[4], x, y, z))
    return atoms


def contact_atoms(
    geometry: Sequence[tuple[str, float, float, float]],
    *,
    molecule_atoms: int,
    anchors: Sequence[int],
    cutoff: float = 3.0,
) -> list[int]:
    """距任一锚原子 ≤``cutoff`` Å 的电极原子（1-based，升序）。

    ``anchors`` 是分子内的 1-based 原子序号；越界或非正数被忽略（不抛错，
    锚定信息来自配置，坏值不该让整份报告崩掉）。
    """
    if molecule_atoms <= 0 or molecule_atoms >= len(geometry):
        return []
    points: list[tuple[float, float, float]] = []
    for index in anchors:
        if isinstance(index, int) and 1 <= index <= molecule_atoms:
            _, x, y, z = geometry[index - 1]
            points.append((x, y, z))
    if not points:
        return []
    limit = float(cutoff) ** 2
    found: list[int] = []
    for position in range(molecule_atoms, len(geometry)):
        _, x, y, z = geometry[position]
        if any((x - ax) ** 2 + (y - ay) ** 2 + (z - az) ** 2 <= limit
               for ax, ay, az in points):
            found.append(position + 1)
    return found


def default_groups(
    n_molecule: int,
    n_atoms: int,
    *,
    contact: Sequence[int] = (),
) -> dict[str, list[int]]:
    """分子 / 接触区 / 电极的 1-based 分区（互不重叠、覆盖全部原子）。

    接触区从电极里拆出来单独显示；不传 ``contact`` 时只给分子+电极两组。
    """
    molecule = list(range(1, max(0, n_molecule) + 1))
    electrode = list(range(max(0, n_molecule) + 1, max(0, n_atoms) + 1))
    contact_set = {index for index in contact if index in set(electrode)}
    if not contact_set:
        return {"molecule": molecule, "electrode": electrode}
    return {
        "molecule": molecule,
        "contact": [index for index in electrode if index in contact_set],
        "electrode": [index for index in electrode if index not in contact_set],
    }


def _moment_stats(values: Sequence[float]) -> dict[str, float | int]:
    """一组自旋矩的计数/净矩/最大绝对值（正负阈值与 ``localized`` 一致）。"""
    moments = [float(value) for value in values]
    return {
        "n_positive": sum(1 for value in moments if value > 1e-8),
        "n_negative": sum(1 for value in moments if value < -1e-8),
        "net": sum(moments),
        "max_abs": max((abs(value) for value in moments), default=0.0),
    }


def summarize(
    table: HirshfeldTable,
    groups: Mapping[str, Sequence[int]],
) -> dict | None:
    """按分区汇总自旋矩；表为空时返回 ``None``（调用方据此跳过展示）。

    分区必须互不重叠且只引用表里存在的原子，否则抛 ``ValueError``——重复
    计数会让"占比"变成假数字，宁可显式失败。

    破缺对称（BS）单重态的证据是"分子内正负自旋矩共存、净矩近零"：为此除
    全表的 ``n_positive``/``n_negative``/``net``/``max_abs`` 外，还给出分子组
    的对应量与 ``localized_antiparallel``（正负共存、单个原子至少 0.05 μB
    且 ``|net| < 0.05*max_abs``）。0.05 μB 的下限用来挡住数值噪声（``ifs-bs``
    的局部矩只有 ~5e-5 μB，净矩/幅度比虽小但毫无意义）。没有 ``molecule``
    分区时这些键为 ``None``——缺证据不冒充"未检出"。
    """
    if not table.atoms:
        return None
    by_index = {atom.index: atom for atom in table.atoms}
    owner: dict[int, str] = {}
    by_group: dict[str, dict] = {}
    for name, members in groups.items():
        indices = list(members)
        values: list[float] = []
        for index in indices:
            atom = by_index.get(index)
            if atom is None:
                raise ValueError("分区 %r 引用了不存在的原子 %d" % (name, index))
            if index in owner:
                raise ValueError("分区重复：原子 %d 同时出现在 %r 与 %r"
                                 % (index, owner[index], name))
            owner[index] = name
            if atom.spin_moment is not None:
                values.append(atom.spin_moment)
        stats = _moment_stats(values)
        by_group[name] = {
            "atoms": len(indices),
            "spin_moment": stats["net"],
            "fraction": None,
            **stats,
        }
    total = table.total or 0.0
    if total:
        for data in by_group.values():
            data["fraction"] = data["spin_moment"] / total
    moments = [atom.spin_moment or 0.0 for atom in table.atoms]
    overall = _moment_stats(moments)
    molecule = by_group.get("molecule")
    if molecule is not None:
        molecule_net = molecule["net"]
        molecule_max_abs = molecule["max_abs"]
        localized_antiparallel = bool(
            molecule["n_positive"] > 0 and molecule["n_negative"] > 0
            and molecule_max_abs >= 0.05
            and abs(molecule_net) < 0.05 * molecule_max_abs)
        molecule_positive = molecule["n_positive"]
        molecule_negative = molecule["n_negative"]
    else:
        molecule_net = molecule_max_abs = None
        molecule_positive = molecule_negative = None
        localized_antiparallel = None
    return {
        "total": table.total,
        "atoms": len(table.atoms),
        "missing": table.missing,
        "groups": by_group,
        "n_positive": overall["n_positive"],
        "n_negative": overall["n_negative"],
        "net": overall["net"],
        "max_abs": overall["max_abs"],
        "localized": overall["n_positive"] > 0 and overall["n_negative"] > 0,
        "molecule_net": molecule_net,
        "molecule_max_abs": molecule_max_abs,
        "molecule_n_positive": molecule_positive,
        "molecule_n_negative": molecule_negative,
        "localized_antiparallel": localized_antiparallel,
    }
