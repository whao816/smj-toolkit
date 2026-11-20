"""Materials Studio / Tripos ``.mol2`` -> 标准 xyz(供 smj 管线使用)。

只取 ``@<TRIPOS>ATOM`` 段的元素与坐标: **坐标原样写出, 不旋转/不平移/不居中** ——
键长/环几何的弛豫是管线第一步 ``molecule_opt``(FHI-aims)的职责, 这里只做格式转换,
避免"转换时顺手改了结构"导致结果不可追溯。

元素取自 Tripos 原子类型(``C.2`` -> ``C``, ``N.am`` -> ``N``); 类型缺失/不认识时
退回原子名(``C1`` -> ``C``)的字母前缀。xyz 注释行写分子名/化学式/锚点序号, 见
``--anchors``。

用法::

    python scripts/mol2_to_xyz.py INPUT.mol2 OUTPUT.xyz --anchors N
    python scripts/mol2_to_xyz.py INPUT.mol2 OUTPUT.xyz --anchors S --anchors N

``--anchors`` 给元素符号时, 额外打印该元素全部原子的 1-based 序号(与 smj.yaml 的
``left_anchor``/``right_anchor`` 同一计数)与两两距离, 便于判左右与初始 Au–锚距离。
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from smj.junction.geometry import Atom, write_xyz  # noqa: E402

#: 双字母元素符号白名单(单字母直接取首字母)。只列分子结里可能出现的重元素。
_TWO_LETTER = frozenset({
    "Cl", "Br", "Si", "Se", "Li", "Na", "Mg", "Al", "Ca", "Ti", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Kr", "Rb", "Sr", "Zr",
    "Mo", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "Xe", "Cs",
    "Ba", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
})


def element_of(atom_name: str, atom_type: str) -> str:
    """mol2 原子类型优先(``C.2``/``N.am``/``Cl``), 退回原子名的字母前缀。"""
    for token in (atom_type, atom_name):
        if not token:
            continue
        letters = "".join(ch for ch in token.split(".")[0] if ch.isalpha())
        if not letters:
            continue
        head = letters[:2].capitalize()
        if len(letters) >= 2 and head in _TWO_LETTER:
            return head
        return letters[0].upper()
    raise ValueError("无法从原子名/类型判定元素: name=%r type=%r" % (atom_name, atom_type))


def parse_mol2(path: str | Path) -> tuple[str, list[tuple[str, float, float, float]],
                                             list[tuple[int, int]]]:
    """读 .mol2, 返回 ``(name, [(element, x, y, z)], [(i, j)])``; 键为 1-based。"""
    name: str | None = None
    atoms: list[tuple[str, float, float, float]] = []
    bonds: list[tuple[int, int]] = []
    section = ""
    for raw in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("@<TRIPOS>"):
            section = line[len("@<TRIPOS>"):].strip().upper()
            continue
        if section == "MOLECULE":
            if name is None:
                name = line.strip().strip('"')
            continue
        if section == "ATOM":
            fields = line.split()
            if len(fields) < 6:
                raise ValueError("mol2 ATOM 行字段不足(需 >=6): %r" % raw)
            atoms.append((element_of(fields[1], fields[5]),
                          float(fields[2]), float(fields[3]), float(fields[4])))
        elif section == "BOND":
            fields = line.split()
            if len(fields) < 3:
                raise ValueError("mol2 BOND 行字段不足(需 >=3): %r" % raw)
            bonds.append((int(fields[1]), int(fields[2])))
    if not atoms:
        raise ValueError("mol2 里没有 @<TRIPOS>ATOM 段或没有原子: %s" % path)
    return (name or Path(path).stem), atoms, bonds


def hill_formula(elements: list[str]) -> str:
    """Hill 式化学式(C 在前, H 随其后, 其余按字母序)。"""
    counts: dict[str, int] = {}
    for el in elements:
        counts[el] = counts.get(el, 0) + 1
    if "C" in counts:
        order = ["C"] + (["H"] if "H" in counts else [])
        order += sorted(el for el in counts if el not in ("C", "H"))
    else:
        order = sorted(counts)
    return "".join(el + (str(counts[el]) if counts[el] > 1 else "") for el in order)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="mol2 -> xyz(坐标原样迁移)")
    ap.add_argument("mol2", help="输入 .mol2 文件")
    ap.add_argument("xyz", help="输出 .xyz 文件")
    ap.add_argument("--anchors", action="append", default=[],
                    metavar="ELEMENT", help="要报告的锚元素, 可重复(如 --anchors N)")
    ap.add_argument("--name", default=None, help="覆盖注释行里的分子名")
    args = ap.parse_args(argv)

    name, atoms, bonds = parse_mol2(args.mol2)
    formula = hill_formula([a[0] for a in atoms])

    parts = [args.name or name, formula, "%d atoms" % len(atoms),
             "from " + Path(args.mol2).name, "coords unchanged"]
    report: list[str] = []
    for el in args.anchors:
        idx = [i + 1 for i, a in enumerate(atoms) if a[0].capitalize() == el.capitalize()]
        if not idx:
            raise ValueError("分子里没有 %s 原子, 无法作锚(现有元素: %s)"
                             % (el, ", ".join(sorted({a[0] for a in atoms}))))
        parts.append("anchor %s: %s" % (el.capitalize(), ",".join(str(i) for i in idx)))
        report.append("  anchor %-2s -> 1-based indices %s" % (el.capitalize(), idx))
        for i in range(len(idx)):
            for j in range(i + 1, len(idx)):
                a, b = atoms[idx[i] - 1], atoms[idx[j] - 1]
                d = math.dist(a[1:], b[1:])
                parts.append("%s%d-%s%d %.3f A" % (el, idx[i], el, idx[j], d))
                report.append("    d(%s%d-%s%d) = %.4f A" % (el, idx[i], el, idx[j], d))

    write_xyz(args.xyz, [Atom(e, x, y, z) for e, x, y, z in atoms],
              comment=" | ".join(parts))
    print("%s -> %s" % (args.mol2, args.xyz))
    print("  name      : %s" % (args.name or name))
    print("  formula   : %s (%d atoms)" % (formula, len(atoms)))
    print("  bonds     : %d (仅用于校验, 不写进 xyz)" % len(bonds))
    for line in report:
        print(line)
    print("  comment   : %s" % " | ".join(parts))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
