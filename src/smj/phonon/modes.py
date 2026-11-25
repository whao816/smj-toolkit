"""虚频（软模）的归属：哪些原子在动、动的是哪类键.

一个虚频只报 ``-21.3 cm⁻¹`` 没有意义 —— 必须回答"谁在动、动什么"才能判断
它是物理软模（分子在接触上摇摆、环扭转）还是数值假象（lead 截断、基组、
力噪声）。本模块把本征矢拆成三样东西：

* **原子权重**：``w_a = Σ_k |v_{3a+k}|²``（本征矢已按质量加权归一），配合
  参与比 ``1/Σ w_a²`` 判断是"整体刚体式摆动"还是"个别键在振"；
* **区域权重**：分子 / 器件区 Au / lead Au 的占比 —— 落在 lead 上说明这个
  软模来自被截断的 lead 力常数，不是分子本身的自由度；
* **键伸缩幅度**：对每根近邻键算沿键方向的反相位移 ``(v_j-v_i)·u``，按幅值
  排序。若最强的是 Au–S，就是接触软模；若是环内 C–C/C–H，就是分子内扭转。

这些数字只依赖力常数矩阵，不碰 DFT，因此可以随时复算（``smj thermal
--analyze-only`` 也会重算并写进 ``thermal.json`` 的 ``soft_modes``）。
"""
from __future__ import annotations

from typing import Iterable, Mapping, Sequence

import numpy as np

from .units import EV_TO_CM1

__all__ = ["describe_soft_modes", "bond_pairs"]

#: 默认的虚频容限（cm⁻¹）。FHI-aims 有限差分在 -10 cm⁻¹ 量级的负本征值属于
#: 噪声（方案 §6 陷阱④、附录 A-5 的同一口径）。
MODE_CM1_TOL = -10.0

#: 近邻键的判定距离（Å）。要覆盖 Au–S(2.35)、S–C(1.77)、C–C(1.39)、C–H(1.09)，
#: 同时排除 Au–Au 最近邻(2.884)。
BOND_CUTOFF_A = 2.6


def _as_positions(positions: np.ndarray) -> np.ndarray:
    array = np.asarray(positions, dtype=float)
    if array.ndim != 2 or array.shape[1] != 3:
        raise ValueError("positions 形状应为 (N,3)，得到 %r" % (array.shape,))
    return array


def bond_pairs(positions: np.ndarray, elements: Sequence[str], *,
               cutoff_A: float = BOND_CUTOFF_A) -> list[dict]:
    """列出近邻键：``{"i","j","label"}``，``label`` 是排好序的元素对。

    只按距离判定（不查成键表）—— 对 Au–分子–Au 这种体系足够：2.6 Å 以内
    只剩真实化学键，Au–Au 的最近邻 2.884 Å 被排除在外。
    """
    pos = _as_positions(positions)
    if len(elements) != pos.shape[0]:
        raise ValueError("elements 长度 %d 与坐标行数 %d 不一致"
                         % (len(elements), pos.shape[0]))
    pairs: list[dict] = []
    cutoff2 = float(cutoff_A) ** 2
    for i in range(pos.shape[0] - 1):
        delta = pos[i + 1:] - pos[i]
        dist2 = np.einsum("ij,ij->i", delta, delta)
        for offset in np.nonzero(dist2 <= cutoff2)[0]:
            j = int(i + 1 + offset)
            label = "-".join(sorted((str(elements[i]).capitalize(),
                                     str(elements[j]).capitalize())))
            pairs.append({"i": i, "j": j, "label": label})
    return pairs


def _region_weights(weights: np.ndarray, regions: Mapping[str, Iterable[int]]) -> dict:
    out: dict[str, float] = {}
    total = float(weights.sum())
    if total <= 0.0:
        return out
    for name, atoms in regions.items():
        index = np.fromiter((int(a) for a in atoms), dtype=int)
        if index.size == 0:
            out[str(name)] = 0.0
            continue
        out[str(name)] = float(weights[index].sum() / total)
    return out


def _bond_stretch(vector: np.ndarray, positions: np.ndarray,
                  pairs: Sequence[Mapping[str, object]]) -> list[dict]:
    """沿键方向的反相位移 ``(v_j - v_i)·u``，按 |stretch|² 排序。"""
    out: list[dict] = []
    for pair in pairs:
        i, j = int(pair["i"]), int(pair["j"])
        delta = positions[j] - positions[i]
        norm = float(np.linalg.norm(delta))
        if norm <= 0.0:
            continue
        unit = delta / norm
        stretch = float((vector[3 * j:3 * j + 3] - vector[3 * i:3 * i + 3]) @ unit)
        out.append({"i": i, "j": j, "label": str(pair["label"]),
                    "stretch": stretch, "weight": stretch * stretch})
    out.sort(key=lambda item: -item["weight"])
    return out


def _rigid_body_fraction(vector: np.ndarray, positions: np.ndarray) -> dict:
    """把一个位移场分解成"整体平移+整体转动"与"内部形变"两部分。

    低频软模最常见的两种身份是**刚体式摆动**（分子绕锚点铰链整体摇）与**分子
    内部扭转**。两者的物理含义完全不同：前者多半来自接触势能的非谐性，后者
    可能预示真实的结构失稳。用 6 个刚体自由度做最小二乘投影即可分辨：

    ``d ≈ B c``，``B`` 的列是 3 个平移与 3 个绕质心的转动，返回
    ``{"fraction": 刚体占比, "axis": 主导转轴, "rotation": 各轴幅度}``。
    """
    pos = positions - positions.mean(axis=0)
    b = np.zeros((vector.size, 6))
    for axis in range(3):
        b[axis::3, axis] = 1.0
    for axis in range(3):
        unit = np.zeros(3)
        unit[axis] = 1.0
        rotation = np.cross(unit, pos)            # (N,3)：绕 axis 转动的位移场
        for component in range(3):
            b[component::3, 3 + axis] = rotation[:, component]
    norms = np.linalg.norm(b, axis=0)
    norms[norms == 0.0] = 1.0
    b = b / norms
    coeffs, *_ = np.linalg.lstsq(b, vector, rcond=None)
    fitted = b @ coeffs
    total = float(vector @ vector)
    fraction = float(fitted @ fitted) / total if total > 0.0 else 0.0
    rotation = np.abs(coeffs[3:])
    axis = "xyz"[int(np.argmax(rotation))] if rotation.size else None
    return {"fraction": fraction, "axis": axis,
            "rotation": {name: float(value) for name, value in
                         zip("xyz", np.abs(coeffs[3:]))},
            "translation": float(np.linalg.norm(coeffs[:3]))}


def describe_soft_modes(
    dynamical_ev2: np.ndarray,
    *,
    masses: Sequence[float] | None = None,
    elements: Sequence[str] | None = None,
    positions: np.ndarray | None = None,
    regions: Mapping[str, Iterable[int]] | None = None,
    cm1_tol: float = MODE_CM1_TOL,
    max_modes: int = 6,
    top_atoms: int = 6,
    top_bonds: int = 5,
    bond_cutoff_A: float = BOND_CUTOFF_A,
) -> list[dict]:
    """把软模（本征值 < 0 且低于 ``cm1_tol``）拆成原子/区域/键三张表。

    返回按频率升序（最软在前）的 list；每一项::

        {"index", "cm1", "participation_ratio",
         "regions": {name: 权重}, "top_atoms": [{"atom","element","weight"}],
         "bonds": [{"pair","label","stretch"}]}

    ``masses`` 给了就把本征矢换成**物理位移** ``v_a/√m_a`` 再统计。这一步不能
    省：``eigh`` 出来的本征矢是质量加权坐标（``D = K/√(m_i m_j)``）下的分量，
    H 的相对振幅被放大 ``√14 ≈ 3.7`` 倍。直接用它会得出"每个低频模都是 C–H
    在伸缩"的错误结论（实测 −21 cm⁻¹ 模的头几根键全是 C–H），而那些模真正
    动的是重金属与骨架。

    没有元素/坐标时只给原子与区域权重。``cm1_tol`` 取正会抛 ``ValueError``
    （正频率不是虚频，调用方多半把符号弄反了）。
    """
    if cm1_tol > 0:
        raise ValueError("cm1_tol 应为非正值（虚频记作负波数），得到 %r" % (cm1_tol,))
    matrix = np.asarray(dynamical_ev2, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("dynamical 应为方阵，得到 %r" % (matrix.shape,))
    if matrix.shape[0] % 3 != 0:
        raise ValueError("dynamical 边长 %d 不是 3 的倍数(应为 3N×3N)"
                         % matrix.shape[0])
    if not np.all(np.isfinite(matrix)):
        raise ValueError("dynamical 含非有限数值")

    matrix = 0.5 * (matrix + matrix.T)
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    n_atoms = matrix.shape[0] // 3

    sqrt_mass = None
    if masses is not None:
        values = np.asarray(masses, dtype=float)
        if values.shape != (n_atoms,):
            raise ValueError("masses 长度 %d 与原子数 %d 不一致"
                             % (values.size, n_atoms))
        if np.any(values <= 0.0) or not np.all(np.isfinite(values)):
            raise ValueError("masses 必须全为正的有限数")
        sqrt_mass = np.repeat(np.sqrt(values), 3)

    if elements is not None and len(elements) != n_atoms:
        raise ValueError("elements 长度 %d 与原子数 %d 不一致" % (len(elements), n_atoms))
    pos = None if positions is None else _as_positions(positions)
    if pos is not None and pos.shape[0] != n_atoms:
        raise ValueError("positions 行数 %d 与原子数 %d 不一致" % (pos.shape[0], n_atoms))

    pairs: list[dict] = []
    if pos is not None:
        # 元素未知时用占位标签，至少能报出"第 i 与第 j 个原子之间的键"
        labels = list(elements) if elements is not None else ["X"] * n_atoms
        pairs = bond_pairs(pos, labels, cutoff_A=bond_cutoff_A)

    out: list[dict] = []
    for index in range(eigenvalues.size):
        value = float(eigenvalues[index])
        if value >= 0.0:
            break                       # 升序，遇到实频即可停
        cm1 = float(-np.sqrt(-value) * EV_TO_CM1)
        if cm1 >= cm1_tol:
            continue
        vector = eigenvectors[:, index]
        if sqrt_mass is not None:
            # 质量加权坐标 -> 物理位移；再归一化，权重才可跨模比较
            vector = vector / sqrt_mass
            norm = float(np.linalg.norm(vector))
            if norm > 0.0:
                vector = vector / norm
        weights = np.zeros(n_atoms)
        for atom in range(n_atoms):
            block = vector[3 * atom:3 * atom + 3]
            weights[atom] = float(block @ block)
        total = float(weights.sum())
        if total > 0.0:
            weights = weights / total
        ratio = float(1.0 / np.sum(weights ** 2)) if total > 0.0 else 0.0
        item: dict = {
            "index": int(index),
            "cm1": cm1,
            "participation_ratio": ratio,
        }
        if regions:
            item["regions"] = _region_weights(weights, regions)
        order = np.argsort(-weights)[:max(int(top_atoms), 1)]
        item["top_atoms"] = [
            {"atom": int(atom),
             "element": str(elements[atom]) if elements is not None else None,
             "weight": float(weights[atom])}
            for atom in order if weights[atom] > 0.0]
        if pairs:
            bonds = _bond_stretch(vector, pos, pairs)[:max(int(top_bonds), 1)]
            item["bonds"] = [{"pair": [b["i"], b["j"]], "label": b["label"],
                              "stretch": b["stretch"]} for b in bonds]
        if pos is not None:
            item["rigid_body"] = _rigid_body_fraction(vector, pos)
        out.append(item)
        if len(out) >= max(int(max_modes), 1):
            break
    return out
