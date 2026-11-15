"""构象筛选前端：SMILES → 3D 构象系综 → 排序 → 落盘。

为什么需要：``molecule_opt`` 只吃**一份** xyz，PBE 弛豫会落到离初猜最近的局部极小。
柔性分子（甲硫基可旋转、长链二硫醇）的构象直接改变锚–锚距离与 S–Au 接触几何，
进而改变 T(E)；"用哪份初猜"不该是手工默认值。

流程（全在 RDKit 内完成，不需要集群）：

1. ETKDGv3 生成 ``n_conformers`` 个候选（``randomSeed=seed``，可复现）；
2. MMFF94s 优化（RDKit），拿到每个构象的力场能量；
3. **先按能量升序**，再按 RMSD 贪心去重——每个构象簇保留最低能代表；
   去重用 ``rdMolAlign.GetBestRMS``（旋转 + 取代基对称等价），不是平移 RMSD：
   ETKDG 给每个构象的取向不同，平移 RMSD 会把"旋转过的同一个构象"当成两个；
   比较的是**重原子骨架**（S–H、甲基氢的朝向差异不构成新构象——这些氢在
   ``aims_tip_relax`` 里会被剥掉）；
4. 取 top-``keep`` 写 ``conformer_<rank>.xyz``，第 1 个同时写成 ``--out``
   （默认 ``molecule.xyz``，可直接交给 ``smj init --xyz``）。

实测注意（RDKit 2026.09.1）：``randomSeed=0`` 会让每个构象用同一初猜，
6 个构象两两 RMS 全为 0.0；因此默认 ``seed=42``。换 seed 会换整批候选，
换的是"粗筛池"，不影响排序规则。

口径（同时写进 ``conformers.json`` 的 ``caveat``）：

* MMFF94s 是力场，不是电子结构判据。它只负责把 20 个候选压到 3 个，
  **不是**"分子在结里就是这个构象"的结论；top-k 应各跑一次下游链
  （或至少 ``molecule_opt`` + ``aims_tip_relax``）再比较。
* MMFF 不含金属参数；本模块只处理裸分子（无 Au/表面），所以适用。
* 能量单位 kcal/mol；不同 SMILES/电荷态之间不可直接比较。
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from ._vendor.junction.geometry import Atom, write_xyz
from ._vendor.junction.molecules import _rdkit

__all__ = ["DEFAULT_RMSD_TOL", "DEFAULT_SEED", "Conformer",
           "select_conformers", "write_conformers"]

#: 构象去重阈值（Å）。0.25 是 RDKit 常识值，与 ``molecules.DEFAULT_RMSD_TOL`` 对齐。
DEFAULT_RMSD_TOL = 0.25
#: 默认随机种子；**不要用 0**（见模块 docstring：RDKit 里 0 会退化成同一初猜）。
DEFAULT_SEED = 42
#: MMFF94s 优化步数上限（与 ``molecule_from_smiles`` 一致）。
_MMFF_MAX_ITERS = 2000


@dataclass(frozen=True)
class Conformer:
    """一个保留的构象；``rank`` 从 1 开始，按 MMFF 能量升序。"""

    rank: int
    energy_kcal: float
    rmsd_to_first_A: float
    mmff_converged: bool
    atoms: tuple[Atom, ...]

    @property
    def atom_count(self) -> int:
        return len(self.atoms)


def _atoms_of(mol, cid: int) -> tuple[Atom, ...]:
    conformer = mol.GetConformer(int(cid))
    return tuple(
        Atom(element=atom.GetSymbol(),
             x=conformer.GetAtomPosition(idx).x,
             y=conformer.GetAtomPosition(idx).y,
             z=conformer.GetAtomPosition(idx).z)
        for idx, atom in enumerate(mol.GetAtoms())
    )


def select_conformers(
    smiles: str,
    *,
    n_conformers: int = 20,
    seed: int = DEFAULT_SEED,
    rmsd_tol: float = DEFAULT_RMSD_TOL,
    keep: int = 3,
) -> list[Conformer]:
    """SMILES → top-``keep`` 构象（MMFF94s 能量升序 + RMSD 去重）。

    Raises:
        ValueError: SMILES 非法、MMFF 缺参数、ETKDG 生成不出构象、参数越界。
        MissingDependency: 未安装 RDKit。
    """
    if n_conformers < 1:
        raise ValueError("n_conformers 至少为 1")
    if keep < 1:
        raise ValueError("keep 至少为 1")
    if rmsd_tol <= 0:
        raise ValueError("rmsd_tol 必须为正数")
    if int(seed) == 0:
        # 不静默改参数：用户显式给 0 时说明后果，让他自己选。
        raise ValueError("seed 不能为 0（RDKit 里 randomSeed=0 会生成同一初猜的重复构象），"
                         "请给 1..2**31-1 之间的值，默认 %d" % DEFAULT_SEED)

    Chem, AllChem = _rdkit()
    from rdkit.Chem import rdMolAlign  # 函数内 import：rdkit 是可选依赖

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("SMILES 解析失败: %r" % (smiles,))
    mol = Chem.AddHs(mol)
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        raise ValueError("MMFF94s 缺该分子的参数（含金属/非常见元素时不可用）: %r"
                         % (smiles,))

    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.useSmallRingTorsions = True
    conf_ids = list(AllChem.EmbedMultipleConfs(
        mol, numConfs=int(n_conformers), params=params))
    if not conf_ids:
        raise ValueError("ETKDGv3 没能为该 SMILES 生成构象: %r" % (smiles,))
    optimized = AllChem.MMFFOptimizeMoleculeConfs(
        mol, maxIters=_MMFF_MAX_ITERS, mmffVariant="MMFF94s")
    # 去重只看重原子骨架：RemoveHs 保留 conformer id，可直接复用 cid 比较。
    heavy = Chem.RemoveHs(mol)

    candidates = sorted(
        ((float(energy), int(cid), int(status) == 0)
         for cid, (status, energy) in zip(conf_ids, optimized)),
        key=lambda item: (item[0], item[1]),
    )

    kept: list[tuple[float, int, bool]] = []
    for energy, cid, converged in candidates:
        if any(rdMolAlign.GetBestRMS(heavy, heavy, cid, other) < rmsd_tol
               for _, other, _ in kept):
            continue
        kept.append((energy, cid, converged))
        if len(kept) >= keep:
            break

    out: list[Conformer] = []
    for rank, (energy, cid, converged) in enumerate(kept, start=1):
        to_first = 0.0 if rank == 1 else float(
            rdMolAlign.GetBestRMS(heavy, heavy, kept[0][1], cid))
        out.append(Conformer(rank=rank, energy_kcal=energy,
                             rmsd_to_first_A=to_first,
                             mmff_converged=converged, atoms=_atoms_of(mol, cid)))
    return out


def write_conformers(
    conformers: Iterable[Conformer],
    outdir: str | Path,
    *,
    out: str = "molecule.xyz",
    manifest: str = "conformers.json",
    smiles: str = "",
    n_conformers: int = 0,
    seed: int = DEFAULT_SEED,
    rmsd_tol: float = DEFAULT_RMSD_TOL,
) -> dict[str, Any]:
    """把选中的构象写成 xyz + 清单；返回 ``{"out", "manifest", "conformers", "data"}``。"""
    items = list(conformers)
    if not items:
        raise ValueError("没有可写的构象")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    entries: list[dict[str, Any]] = []
    paths: list[str] = []
    for item in items:
        name = "conformer_%02d.xyz" % item.rank
        path = outdir / name
        write_xyz(path, item.atoms,
                  comment="MMFF94s %.4f kcal/mol (rank %d)" % (item.energy_kcal, item.rank))
        entries.append({
            "rank": item.rank,
            "file": name,
            "energy_kcal": round(item.energy_kcal, 6),
            "rmsd_to_first_A": round(item.rmsd_to_first_A, 4),
            "mmff_converged": item.mmff_converged,
            "atom_count": item.atom_count,
        })
        paths.append(str(path))

    out_path = outdir / out
    write_xyz(out_path, items[0].atoms,
              comment="MMFF94s rank 1, %.4f kcal/mol" % items[0].energy_kcal)

    try:  # 只用于记录，不参与判定
        from rdkit import rdBase

        rdkit_version: str | None = rdBase.rdkitVersion
    except Exception:  # pragma: no cover - 只有 rdkit 缺失才会到这里
        rdkit_version = None

    data: dict[str, Any] = {
        "generator": "smj conformers",
        "smiles": smiles,
        "seed": int(seed),
        "requested_conformers": int(n_conformers),
        "rmsd_tolerance_A": float(rmsd_tol),
        "kept": len(entries),
        "energy_unit": "kcal/mol (MMFF94s)",
        "rdkit_version": rdkit_version,
        "selection": "MMFF94s energy ascending, then greedy GetBestRMS dedupe "
                     "(lowest-energy representative per cluster)",
        "caveat": "MMFF94s/ETKDGv3 只作构象粗筛，不是电子结构判据；"
                  "top-k 应各跑一次下游链再比较",
        "conformers": entries,
    }
    manifest_path = outdir / manifest
    manifest_path.write_text(
        json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return {"out": str(out_path), "manifest": str(manifest_path),
            "conformers": paths, "data": data}
