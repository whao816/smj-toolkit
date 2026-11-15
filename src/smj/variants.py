"""分子化学态对照臂：硫醇盐（去质子）与金帽（S–Au）模型。

``molecule_opt`` 优化的是**中性**硫醇/硫醚，而器件里锚定原子是与 Au 成键的。
两个 demo 同为 2CN-656 硫醚，建成结后的 S–Au 接触却给出 2.6877 Å / 141.5°
（硫醚）与 2.4267 Å / 111.4°（硫醇盐）——差 0.26 Å / 30°。这条差异不该只靠
``aims_tip_relax`` 的局部弛豫去掩盖，所以提供两条**对照臂**，让"质子化态
有没有系统性偏差"可以被量化而不是被假设：

* ``thiolate``：把 S–H 的 H 去掉（每个去质子 +1 净负电），中性分子 → 阴离子。
  电子数不变（去掉 1 个 H⁺，核电荷 −1、电荷 −1 相抵），因此奇偶校验与多重度
  不变。
* ``aucap``：把每个 S–H 的 H 换成 Au，放在原 S–H 方向上 ``au_distance``（默认
  2.40 Å，与 demo 硫醇盐接触 2.43 Å 一致）。电荷不变，中性。

两条臂都**只生成输入**（xyz + manifest），真正的判定仍由下游
``molecule_opt → aims_tip_relax`` 跑出来；本模块不替用户选构象、不声称哪条臂
"更对"。``variants.json`` 里给出每条臂的 ``config_hint``（charge / multiplicity /
xyz），便于直接跑一次对照。

口径边界（写进 manifest 的 ``caveat``）：

* 都是气相、0 K、无溶剂/无电极极化；
* 硫醇盐臂是**裸阴离子**，PBE/light 下 HOMO 可能高于真空能级（无束缚态），
  只适合看 C–S 键长/接触几何的**趋势**，不能当成器件里的化学态；
* 金帽臂的 Au 位置是几何启发（沿 S–H 方向 2.40 Å），不是能量极小，需要下游
  弛豫；单原子 Au 也不是真实电极表面。
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from ._vendor.junction.geometry import Atom, read_xyz, write_xyz
from .checks import total_electrons

__all__ = [
    "DEFAULT_ANCHOR_ELEMENTS",
    "DEFAULT_AU_DISTANCE",
    "DEFAULT_SH_CUTOFF",
    "Variant",
    "build_variants",
    "thiolate_variant",
    "aucap_variant",
    "write_variants",
]

DEFAULT_ANCHOR_ELEMENTS = ("S", "Se")
#: 锚–H 距离小于该值（Å）视为可去质子/可替换的氢（S–H ≈ 1.34，Se–H ≈ 1.47）。
DEFAULT_SH_CUTOFF = 1.70
#: 金帽沿原 S–H 方向放 Au 的距离（Å）；demo 硫醇盐 S–Au 接触 2.4267 Å。
DEFAULT_AU_DISTANCE = 2.40

CAVEAT = (
    "两条臂都是气相/0 K/无溶剂无电极极化；硫醇盐臂是裸阴离子，PBE/light 下 "
    "HOMO 可能高于真空能级（无束缚态），只适合比较 C–S 键长与接触几何的趋势；"
    "金帽臂的 Au 位置是沿原 S–H 方向的几何启发（默认 2.40 Å），不是能量极小，"
    "单原子 Au 也不是真实电极表面。真正的判定要用下游 molecule_opt → aims_tip_relax "
    "跑出来，不能只看本文件。"
)


@dataclass(frozen=True)
class Variant:
    """一条对照臂：几何 + 电荷/多重度 + 为什么这样做。"""

    name: str
    atoms: tuple[Atom, ...]
    charge: int
    multiplicity: int
    removed: tuple[int, ...]      # 被移除的 H 的 1-based 序号（原始几何）
    note: str

    @property
    def atom_count(self) -> int:
        return len(self.atoms)

    def config_hint(self) -> dict:
        """可直接抄进 ``smj.yaml`` 的 ``molecule`` 段。"""
        return {"charge": self.charge, "multiplicity": self.multiplicity}


def _distance(a: Atom, b: Atom) -> float:
    return math.dist(a.coord(), b.coord())


def _anchor_hydrogens(
    atoms: Sequence[Atom],
    anchor_elements: Iterable[str],
    cutoff: float,
) -> list[tuple[int, int]]:
    """返回 ``[(h_index, anchor_index), ...]``（0 基）：锚–H 距离 < ``cutoff``。"""
    anchors = set(anchor_elements)
    pairs: list[tuple[int, int]] = []
    for h, atom in enumerate(atoms):
        if atom.element != "H":
            continue
        best: tuple[float, int] | None = None
        for j, other in enumerate(atoms):
            if other.element not in anchors:
                continue
            d = _distance(atom, other)
            if d < cutoff and (best is None or d < best[0]):
                best = (d, j)
        if best is not None:
            pairs.append((h, best[1]))
    return pairs


def _drop(atoms: Sequence[Atom], drop: set[int]) -> tuple[Atom, ...]:
    return tuple(a for i, a in enumerate(atoms) if i not in drop)


def thiolate_variant(
    atoms: Sequence[Atom],
    *,
    charge: int = 0,
    multiplicity: int = 1,
    anchor_elements: Iterable[str] = DEFAULT_ANCHOR_ELEMENTS,
    sh_cutoff: float = DEFAULT_SH_CUTOFF,
    name: str = "thiolate",
) -> Variant:
    """去掉锚定原子上的 H（每个 S–H → S⁻，总电荷 −1）。"""
    pairs = _anchor_hydrogens(atoms, anchor_elements, sh_cutoff)
    if not pairs:
        raise ValueError(
            "thiolate_variant: %s 上没有找到锚–H（距离 < %.2f Å），无法生成硫醇盐臂；"
            "硫醚/无 S–H 的分子本来就没有可对照的质子化态"
            % ("/".join(anchor_elements), sh_cutoff))
    removed = sorted({h for h, _ in pairs})
    new_charge = int(charge) - len(removed)
    kept = _drop(atoms, set(removed))
    return Variant(
        name=name, atoms=kept, charge=new_charge, multiplicity=int(multiplicity),
        removed=tuple(i + 1 for i in removed),
        note=("去掉 %d 个锚–H（S–H/Se–H < %.2f Å）：电荷 %d → %d。电子数不变"
              "（每去一个 H⁺，核电荷与电荷各减 1），多重度校验不受影响"
              % (len(removed), sh_cutoff, charge, new_charge)))


def aucap_variant(
    atoms: Sequence[Atom],
    *,
    charge: int = 0,
    multiplicity: int = 1,
    anchor_elements: Iterable[str] = DEFAULT_ANCHOR_ELEMENTS,
    sh_cutoff: float = DEFAULT_SH_CUTOFF,
    au_distance: float = DEFAULT_AU_DISTANCE,
    element: str = "Au",
    name: str = "aucap",
) -> Variant:
    """把锚–H 换成 ``element``（默认 Au），沿原 H 方向放在 ``au_distance`` 处。"""
    pairs = _anchor_hydrogens(atoms, anchor_elements, sh_cutoff)
    if not pairs:
        raise ValueError(
            "aucap_variant: %s 上没有找到锚–H（距离 < %.2f Å），无法放金帽"
            % ("/".join(anchor_elements), sh_cutoff))
    removed = sorted({h for h, _ in pairs})
    caps: list[Atom] = []
    for h, anchor in sorted(pairs, key=lambda p: p[0]):
        a = atoms[anchor]
        b = atoms[h]
        d = _distance(a, b)
        direction = ((b.x - a.x) / d, (b.y - a.y) / d, (b.z - a.z) / d)
        caps.append(Atom(element,
                         a.x + direction[0] * au_distance,
                         a.y + direction[1] * au_distance,
                         a.z + direction[2] * au_distance))
    kept = list(_drop(atoms, set(removed))) + caps
    return Variant(
        name=name, atoms=tuple(kept), charge=int(charge),
        multiplicity=int(multiplicity), removed=tuple(i + 1 for i in removed),
        note=("把 %d 个锚–H 换成 %s（沿原 H 方向 %.2f Å，电荷不变）：这是局部模型，"
              "不是电极表面；Au 位置是几何启发，必须由下游弛豫找极小"
              % (len(removed), element, au_distance)))


def build_variants(
    atoms: Sequence[Atom],
    *,
    charge: int = 0,
    multiplicity: int = 1,
    kinds: Iterable[str] = ("thiolate", "aucap"),
    anchor_elements: Iterable[str] = DEFAULT_ANCHOR_ELEMENTS,
    sh_cutoff: float = DEFAULT_SH_CUTOFF,
    au_distance: float = DEFAULT_AU_DISTANCE,
) -> list[Variant]:
    """按 ``kinds`` 生成对照臂（顺序保持给定顺序）。"""
    out: list[Variant] = []
    for kind in kinds:
        if kind == "thiolate":
            out.append(thiolate_variant(
                atoms, charge=charge, multiplicity=multiplicity,
                anchor_elements=anchor_elements, sh_cutoff=sh_cutoff))
        elif kind == "aucap":
            out.append(aucap_variant(
                atoms, charge=charge, multiplicity=multiplicity,
                anchor_elements=anchor_elements, sh_cutoff=sh_cutoff,
                au_distance=au_distance))
        else:
            raise ValueError("未知对照臂 %r：只支持 thiolate / aucap" % (kind,))
    return out


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_variants(
    variants: Iterable[Variant],
    outdir: str | Path,
    *,
    source: str | Path | None = None,
) -> dict:
    """把每条臂写成 ``<name>.xyz`` 并汇总 ``variants.json``。"""
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    entries: list[dict] = []
    for variant in variants:
        xyz = out / ("%s.xyz" % variant.name)
        write_xyz(xyz, variant.atoms,
                  comment="%s variant of %s (%d atoms, charge %+d)"
                          % (variant.name, Path(source).name if source else "input",
                             variant.atom_count, variant.charge))
        entry = {
            "name": variant.name,
            "xyz": str(xyz),
            "charge": variant.charge,
            "multiplicity": variant.multiplicity,
            "atom_count": variant.atom_count,
            "electron_count": total_electrons(variant.atoms, variant.charge),
            "removed_h_indices": list(variant.removed),
            "note": variant.note,
            "config_hint": variant.config_hint(),
        }
        entries.append(entry)
    manifest = {
        "source": str(source) if source else None,
        "source_sha256": _sha256(Path(source)) if source else None,
        "caveat": CAVEAT,
        "variants": entries,
    }
    path = out / "variants.json"
    path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8")
    return {"manifest": str(path), "variants": entries, "data": manifest}


def load_source(xyz_path: str | Path) -> list[Atom]:
    """读输入 xyz（保留顺序，便于把 removed 索引指回原文件）。"""
    return list(read_xyz(xyz_path))
