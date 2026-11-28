"""器件和体相电极共用的面内周期距离。"""
from __future__ import annotations

import numpy as np


def atom_distances(positions, lattice=None) -> np.ndarray:
    """原子对的最短距离（Å）。只允许沿前两条晶格矢量平移。"""
    xyz = np.asarray(positions, dtype=float)
    delta = xyz[:, None, :] - xyz[None, :, :]
    if lattice is None:
        return np.linalg.norm(delta, axis=-1)

    vectors = np.asarray(lattice, dtype=float)
    if vectors.ndim != 2 or vectors.shape[0] < 2 or vectors.shape[1] != 3:
        raise ValueError("lattice 必须含至少两条三维晶格矢量")
    a, b = vectors[:2]
    if not np.isfinite([a, b]).all() or np.linalg.norm(np.cross(a, b)) < 1.0e-12:
        raise ValueError("面内晶格矢量必须有限且线性无关")

    # Gauss 约化保留同一格点集，使 |a| ≤ |b| 且 |a·b| ≤ |a|²/2。
    while True:
        if a @ a > b @ b:
            a, b = b, a
        multiple = np.rint((a @ b) / (a @ a))
        if multiple == 0:
            break
        b = b - multiple * a

    b_perp = b - (a @ b) / (a @ a) * a
    nearest_b = np.rint((delta @ b_perp) / (b_perp @ b_perp))
    # 每条平行于 a 的格线上，最优 a 系数由一维投影给出。约化后
    # |b_perp| ≥ √3 |a|/2，更远的格线无法优于最近三条中的候选。
    best_squared = np.full(delta.shape[:2], np.inf)
    for offset in (-1, 0, 1):
        remainder = delta - (nearest_b + offset)[..., None] * b
        nearest_a = np.rint((remainder @ a) / (a @ a))
        image = remainder - nearest_a[..., None] * a
        best_squared = np.minimum(best_squared, np.einsum("ijk,ijk->ij", image, image))
    return np.sqrt(best_squared)
