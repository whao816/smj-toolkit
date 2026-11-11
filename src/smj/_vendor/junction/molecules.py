"""分子前端(W3-1): SMILES → 3D、E/Z 判定、构象采样、去重、锚定原子定位。

**RDKit 是可选依赖**(`pyproject.toml` 的 `dft` extra)。凡是需要 RDKit 的入口都会先
调 :func:`_rdkit`, 缺依赖时抛 :class:`MissingDependency`, 消息里写清装法与降级路径
(给现成 xyz / geometry.in 就能跳过整个分子前端)。CI 侧由 ``ci/check_no_rdkit.py``
守住"工具模块顶层不许裸 import rdkit"。

不依赖 RDKit 的能力(手上有几何时仍可用):

* :func:`load_molecule` —— 读 xyz / geometry.in(RDKit 无关)
* :func:`anchor_candidates` —— 按元素找候选锚原子
* :func:`dedupe_conformers` —— 同一分子的多构象按 RMSD 去重
* :func:`write_xyz` —— 落盘
* :func:`molecule_signature` —— 组成指纹(元素计数 + 原子数), 进 provenance
"""

from __future__ import annotations

import itertools
import math
from pathlib import Path
from typing import Iterable

from .geometry import Atom, read_aims_geometry, read_xyz

#: 做去重时的 RMSD 阈值(A)。经验值: 小于 0.25 A 视为同一构象。
DEFAULT_RMSD_TOL = 0.25


class MissingDependency(RuntimeError):
    """可选依赖缺失(RDKit)。"""


# ---------------------------------------------------------------------------
# 依赖探测
# ---------------------------------------------------------------------------

def rdkit_available() -> bool:
    """本机是否装了 RDKit(不抛异常, 供 CLI 决定是否走分子前端)。"""
    try:
        import rdkit  # noqa: F401
    except Exception:
        return False
    return True


def _rdkit():
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception as exc:  # pragma: no cover - 依赖包缺
        raise MissingDependency(
            '需要 RDKit 才能从 SMILES 生成几何。安装: pip install "agenticdft[dft]" '
            '(或 pip install rdkit)。若已有 xyz/geometry.in, 直接给几何即可跳过分子前端。'
        ) from exc
    return Chem, AllChem


# ---------------------------------------------------------------------------
# 几何 I/O(无 RDKit 依赖)
# ---------------------------------------------------------------------------

def load_molecule(path: str | Path) -> list[Atom]:
    """读 xyz 或 geometry.in(按扩展名分派)。

    注意: ``geometry.in`` 是 FHI-aims 格式, 其中 ``#`` 行是注释, 与 xyz 不同。
    """
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError('找不到分子几何文件: %s' % p)
    suffix = p.suffix.lower()
    if suffix == '.xyz':
        return read_xyz(p)
    if suffix == '.in' or p.name.endswith('geometry.in'):
        return read_aims_geometry(p)
    raise ValueError('不认识的分子几何格式: %s(支持 .xyz / .in)' % p.suffix)


def write_xyz(atoms: Iterable[Atom], path: str | Path, comment: str = '') -> Path:
    """写标准 xyz(元素 4 位左对齐 + 三列坐标); 返回写出的路径。"""
    atoms = list(atoms)
    if not atoms:
        raise ValueError('原子列表为空, 不写 xyz')
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    lines = ['%d' % len(atoms), comment]
    for a in atoms:
        lines.append('%-4s %14.8f %14.8f %14.8f' % (a.element, a.x, a.y, a.z))
    p.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return p


def molecule_signature(atoms: Iterable[Atom]) -> dict:
    """组成指纹: 原子数 + 各元素计数(字典按元素排序, 便于直接写进 JSON)。"""
    counts: dict[str, int] = {}
    total = 0
    for a in atoms:
        counts[a.element] = counts.get(a.element, 0) + 1
        total += 1
    if not total:
        raise ValueError('原子列表为空')
    return {'atom_count': total, 'elements': dict(sorted(counts.items()))}


def anchor_candidates(atoms: Iterable[Atom], element: str = 'S') -> list[int]:
    """按元素找候选锚原子(1-based), 例如硫醇结的两端 S。

    找不到就返回空列表 —— 调用方(构型枚举)据此报"这个分子没有可锚原子"。
    """
    out = []
    for k, a in enumerate(atoms):
        if a.element.strip().capitalize() == element.strip().capitalize():
            out.append(k + 1)
    return out


# ---------------------------------------------------------------------------
# 构象去重(无 RDKit 依赖)
# ---------------------------------------------------------------------------

def _centered(atoms: list[Atom]) -> list[tuple[float, float, float]]:
    n = len(atoms)
    cx = sum(a.x for a in atoms) / n
    cy = sum(a.y for a in atoms) / n
    cz = sum(a.z for a in atoms) / n
    return [(a.x - cx, a.y - cy, a.z - cz) for a in atoms]


def rmsd(a: list[Atom], b: list[Atom]) -> float:
    """同一原子顺序下的直接 RMSD(只做平移对齐; 不做旋转/对称匹配)。

    构象来自同一 SMILES、原子顺序一致, 所以到这里足够了; 真要对称匹配得接 RDKit。
    """
    if len(a) != len(b):
        raise ValueError('原子数不同, 无法比较: %d vs %d' % (len(a), len(b)))
    if not a:
        raise ValueError('空分子')
    ca, cb = _centered(a), _centered(b)
    s = sum((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2
            for p, q in zip(ca, cb))
    return math.sqrt(s / len(ca))


def dedupe_conformers(conformers: Iterable[Iterable[Atom]],
                      tol: float = DEFAULT_RMSD_TOL) -> list[list[Atom]]:
    """按顺序贪心去重: 与已保留的每个构象 RMSD 都大于 ``tol`` 才保留。"""
    if tol <= 0:
        raise ValueError('tol 必须为正')
    kept: list[list[Atom]] = []
    for conf in conformers:
        atoms = list(conf)
        if not atoms:
            raise ValueError('构象里没有原子')
        if any(rmsd(atoms, prev) <= tol for prev in kept):
            continue
        kept.append(atoms)
    return kept


# ---------------------------------------------------------------------------
# SMILES → 3D(RDKit)
# ---------------------------------------------------------------------------

def smile_to_mol(smiles: str):
    """SMILES → RDKit Mol(带显式氢); 非法 SMILES 抛 ValueError。"""
    Chem, _ = _rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError('SMILES 解析失败: %r' % (smiles,))
    return Chem.AddHs(mol)


def isomer_of(smiles: str) -> str:
    """从 SMILES 判 E/Z。

    先看 RDKit 标出的双键立体化学(``STEREOE``/``STEREOZ``); 没标就用 2D 坐标算
    双键两侧取代基的对位关系。判不出来时返回 ``''``(调用方必须显式给出异构体,
    不允许"猜一个")。
    """
    Chem, _ = _rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError('SMILES 解析失败: %r' % (smiles,))
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        stereo = bond.GetStereo()
        if stereo == Chem.BondStereo.STEREOE:
            return 'E'
        if stereo == Chem.BondStereo.STEREOZ:
            return 'Z'
    return ''


def molecule_from_smiles(smiles: str, *, seed: int = 0,
                         n_conformers: int = 1,
                         optimize: bool = True,
                         dedupe_tol: float | None = None) -> list[list[Atom]]:
    """SMILES → 一组大致合理的 3D 构象(ETKDGv3 + MMFF94s)。

    返回 ``list[list[Atom]]``: 每个元素是一个构象。``optimize=False`` 时只做 ETKDG
    嵌入(快, 用于粗筛)。去重阈值缺省 ``DEFAULT_RMSD_TOL``, 传 ``None`` 则不去重。
    """
    if n_conformers < 1:
        raise ValueError('n_conformers 至少为 1')
    Chem, AllChem = _rdkit()
    mol = smile_to_mol(smiles)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.useSmallRingTorsions = True
    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=int(n_conformers), params=params))
    if not conf_ids:
        raise ValueError('ETKDG 没能为该 SMILES 生成构象: %r' % (smiles,))
    if optimize:
        AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000, mmffVariant='MMFF94s')

    out: list[list[Atom]] = []
    conf = mol.GetConformer()
    for cid in conf_ids:
        atoms = []
        for idx, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(idx) if conf.GetNumAtoms() == mol.GetNumAtoms() \
                else mol.GetConformer(int(cid)).GetAtomPosition(idx)
            atoms.append(Atom(element=atom.GetSymbol(), x=pos.x, y=pos.y, z=pos.z))
        out.append(atoms)
    if dedupe_tol is None:
        return out
    return dedupe_conformers(out, tol=dedupe_tol)


def smiles_pair(smiles_trans: str, smiles_cis: str) -> dict[str, str]:
    """把用户给的两个 SMILES 归档成 ``{'E': ..., 'Z': ...}``(判不出来就抛错)。

    框架的输入表里 E/Z 是两个字符串; 这里只做**校验**(标注的异构必须与 SMILES 的
    立体化学一致), 不替用户"反推"另一个异构体 —— 手性写错在单分子结里代价太大。
    """
    out = {}
    for smiles in (smiles_trans, smiles_cis):
        tag = isomer_of(smiles)
        if tag not in ('E', 'Z'):
            raise ValueError('SMILES 没有标出 E/Z 立体化学: %r(请用 C/N=N/C 这类写法)'
                             % (smiles,))
        if tag in out:
            raise ValueError('两个 SMILES 都判成 %s, 需要各给一个 E 和一个 Z' % tag)
        out[tag] = smiles
    return out
