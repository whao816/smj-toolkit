"""FHI-aims cube 文件解析与"锚定轨道权重"（纯 numpy）。

``control_relax.in`` 请求了 4 个前线轨道的 cube（``output cube eigenstate
HOMO-1/HOMO/LUMO/LUMO+1``），FHI-aims 写成
``cube_<序号>_eigenstate_<本征态>_spin_<自旋>.cube``。原先这些文件既不收集也
不解析，白花 IO；这里把它们变成可解释的报告量：

* **锚定轨道权重** = 以锚定原子（S/N/Se）为球心、半径 ``radius`` 的球内
  |ψ|² 之和 / 全箱 |ψ|² 之和。权重高说明该前线轨道在锚定原子上有振幅，
  与金属耦合有通道；权重接近 0 说明这个轨道根本不经过锚定原子。
* 只算**占比**，所以不受 FHI-aims 那条 "内容单位是 1/Å³、体素单位是原子单位"
  的警告影响（分子分母同单位抵消）。

两个文件格式上的坑，代码里显式处理：

1. cube 坐标与步长按标准 cube 约定是 **bohr**，这里乘 ``BOHR`` 换成 Å；
2. 文件不说明存的是 ψ 还是 |ψ|²：出现负值即按带符号的 ψ 处理（平方后再积分），
   全非负则视为已是密度。两种情况的占比都可比较，但**同一批文件应同源**。
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ._vendor.junction.geometry import Atom
from .checks import DEFAULT_ANCHOR_ELEMENTS

__all__ = [
    "BOHR",
    "DEFAULT_ANCHOR_RADIUS",
    "Cube",
    "anchor_weight",
    "frontier_weights",
    "read_cube",
]

#: bohr → Å（CODATA 2018）。
BOHR = 0.529177210903

#: 锚定球半径（Å）：覆盖 S/N/Se 及其第一近邻的一圈，足以代表"锚定原子上的振幅"。
DEFAULT_ANCHOR_RADIUS = 1.6

_CUBE_NAME_RE = re.compile(r"^cube_(\d+)_eigenstate_(\d+)_spin_(\d+)\.cube$")
_ORBIT_RE = re.compile(r"^\s*output\s+cube\s+eigenstate\s+(\S+)\s*$", re.IGNORECASE)


@dataclass(frozen=True)
class Cube:
    """一个标准 cube 文件。``values`` 形状 ``(nx, ny, nz)``，**x 下标变化最快**。"""

    atoms: tuple[Atom, ...]
    origin: tuple[float, float, float]  # Å
    steps: tuple[float, float, float]   # 每个下标轴的步长，Å
    shape: tuple[int, int, int]
    values: np.ndarray
    comments: tuple[str, ...]

    @property
    def signed(self) -> bool:
        """文件里存的是带符号的 ψ（有负值）还是 |ψ|²。"""
        return bool(np.any(self.values < 0.0))


def _axis_vectors(lines: list[str]) -> tuple[list[int], list[tuple[float, float, float]]]:
    counts: list[int] = []
    vectors: list[tuple[float, float, float]] = []
    for line in lines:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError("cube 头部步长行残缺：%r" % line)
        counts.append(int(parts[0]))
        vectors.append(tuple(float(v) for v in parts[1:4]))
    return counts, vectors


def read_cube(path: str | Path) -> Cube:
    """读取标准 Gaussian cube（FHI-aims 输出）。体素网格必须是正交/轴对齐的。"""
    path = Path(path)
    if not path.exists():
        raise ValueError("cube 文件不存在：%s" % path)
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 7:
        raise ValueError("cube 文件太短（%d 行）：%s" % (len(lines), path.name))
    comments = (lines[0], lines[1])

    head = lines[2].split()
    if len(head) < 4:
        raise ValueError("cube 第 3 行应为 \"natoms ox oy oz\"：%r" % lines[2])
    natoms = int(head[0])
    origin_au = tuple(float(v) for v in head[1:4])
    counts, vectors = _axis_vectors(lines[3:6])
    if not (0 <= natoms <= len(lines) - 6):
        raise ValueError("cube 声明 %d 个原子，但文件只剩 %d 行" % (natoms, len(lines) - 6))

    atom_lines = lines[6:6 + natoms]
    atoms: list[Atom] = []
    for line in atom_lines:
        parts = line.split()
        if len(parts) < 5:
            raise ValueError("cube 原子行残缺：%r" % line)
        atoms.append(Atom(
            element=_element_symbol(int(parts[0]), parts),
            x=float(parts[2]) * BOHR,
            y=float(parts[3]) * BOHR,
            z=float(parts[4]) * BOHR,
        ))

    # 只支持轴对齐网格：FHI-aims 默认如此；非对角时宁可报错也不乱积分。
    matrix = np.array(vectors, dtype=float)  # 每行 = 一个下标轴的步长向量（bohr）
    off = matrix - np.diag(np.diag(matrix))
    if np.max(np.abs(off)) > 1e-6:
        raise ValueError("cube 网格不是轴对齐的（步长矩阵有非对角项），无法按球积分：%s"
                         % path.name)
    steps = tuple(float(v) * BOHR for v in np.diag(matrix))

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for _ in range(6 + natoms):
            handle.readline()
        values = _read_values(handle, path)

    shape = (counts[0], counts[1], counts[2])
    expected = shape[0] * shape[1] * shape[2]
    if values.size != expected:
        raise ValueError("cube 体素数不符：声明 %d×%d×%d=%d，实际 %d"
                         % (shape[0], shape[1], shape[2], expected, values.size))
    # x 下标变化最快 → 按 Fortran 顺序 reshape，使 values[i, j, k] 的 i 就是 x 轴。
    grid = values.reshape(shape, order="F")
    return Cube(
        atoms=tuple(atoms),
        origin=tuple(float(v) * BOHR for v in origin_au),
        steps=steps,
        shape=shape,
        values=grid,
        comments=comments,
    )


#: 每次转 float 的 token 数：cube 动辄 7e6 个体素，整表 split 后再转会多占几百 MB。
_CHUNK = 200_000


def _read_values(handle, path: Path) -> np.ndarray:
    """流式读取体素值（cube 的最后一行常不满 6 列，loadtxt 会因此报错）。"""
    parts: list[np.ndarray] = []
    buffer: list[str] = []
    for line in handle:
        tokens = line.split()
        if not tokens:
            continue
        buffer.extend(tokens)
        if len(buffer) >= _CHUNK:
            parts.append(_floats(buffer, path))
            buffer = []
    if buffer:
        parts.append(_floats(buffer, path))
    return np.concatenate(parts) if parts else np.empty(0, dtype=float)


def _floats(tokens: list[str], path: Path) -> np.ndarray:
    try:
        return np.array(tokens, dtype=float)
    except ValueError as exc:
        raise ValueError("cube 体素值无法解析（%s）：%s" % (path.name, exc)) from exc


def _element_symbol(atomic_number: int, parts: list[str]) -> str:
    """cube 原子行的第 1 列是原子序数；第 6 列（如果有）通常是元素名。"""
    table = _SYMBOLS.get(atomic_number)
    if table is not None:
        return table
    if len(parts) >= 6 and parts[5].isalpha():
        return parts[5].capitalize()
    raise ValueError("无法识别 cube 原子序数 %d" % atomic_number)


def _axes(cube: Cube) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ox, oy, oz = cube.origin
    dx, dy, dz = cube.steps
    return (ox + np.arange(cube.shape[0]) * dx,
            oy + np.arange(cube.shape[1]) * dy,
            oz + np.arange(cube.shape[2]) * dz)


def _scores(cube: Cube) -> np.ndarray:
    """按体素给出可加和的"密度"：带符号 ψ → ψ²，非负 → 原值。"""
    return np.square(cube.values) if cube.signed else cube.values


def anchor_weight(cube: Cube, *, anchor_elements: tuple[str, ...] = DEFAULT_ANCHOR_ELEMENTS,
                  radius: float = DEFAULT_ANCHOR_RADIUS) -> dict[str, object]:
    """锚定球内的轨道密度占比。

    返回 ``{"fraction", "total", "signed_cube", "radius_A", "anchors": [...]}``；
    ``fraction`` 为 ``None`` 表示全箱密度为 0（该轨道在本网格上无振幅），
    调用方应把它当"异常但可报告"处理，不要除零。
    """
    if radius <= 0:
        raise ValueError("radius 必须为正（Å）")
    wanted = {str(e).capitalize() for e in anchor_elements}
    score = _scores(cube)
    total = float(score.sum())

    xs, ys, zs = _axes(cube)
    r2 = float(radius) ** 2
    union = np.zeros(cube.shape, dtype=bool)
    anchors: list[dict[str, object]] = []
    for index, atom in enumerate(cube.atoms, start=1):
        if atom.element.capitalize() not in wanted:
            continue
        members = _sphere_slices(xs, ys, zs, (atom.x, atom.y, atom.z), radius)
        if members is None:
            # 球完全在网格外：权重为 0，但仍要报告（说明网格没罩住这个原子）。
            anchors.append({"index": index, "element": atom.element,
                            "weight": 0.0, "fraction": 0.0 if total else None})
            continue
        (ix0, ix1), (iy0, iy1), (iz0, iz1), mask = members
        weight = float(score[ix0:ix1, iy0:iy1, iz0:iz1][mask].sum())
        union[ix0:ix1, iy0:iy1, iz0:iz1] |= mask
        anchors.append({
            "index": index,
            "element": atom.element,
            "weight": weight,
            "fraction": (weight / total) if total else None,
        })
    union_weight = float(score[union].sum()) if anchors else 0.0
    return {
        "fraction": (union_weight / total) if total else None,
        "total": total,
        "signed_cube": cube.signed,
        "radius_A": float(radius),
        "anchors": anchors,
    }


def _sphere_slices(xs: np.ndarray, ys: np.ndarray, zs: np.ndarray,
                   center: tuple[float, float, float],
                   radius: float):
    """球在网格内的包围盒与局部 mask；完全在网格外返回 None。"""
    cx, cy, cz = center
    ix0, ix1 = _window(xs, cx, radius)
    iy0, iy1 = _window(ys, cy, radius)
    iz0, iz1 = _window(zs, cz, radius)
    if ix0 >= ix1 or iy0 >= iy1 or iz0 >= iz1:
        return None
    gx, gy, gz = xs[ix0:ix1], ys[iy0:iy1], zs[iz0:iz1]
    d2 = ((gx[:, None, None] - cx) ** 2
          + (gy[None, :, None] - cy) ** 2
          + (gz[None, None, :] - cz) ** 2)
    mask = d2 <= radius * radius
    if not mask.any():
        return None
    return (ix0, ix1), (iy0, iy1), (iz0, iz1), mask


def _window(axis: np.ndarray, center: float, radius: float) -> tuple[int, int]:
    """网格点落在 [center-radius, center+radius] 的半开下标区间。"""
    lo = int(np.searchsorted(axis, center - radius, side="left"))
    hi = int(np.searchsorted(axis, center + radius, side="right"))
    return max(lo, 0), min(hi, axis.size)


def requested_orbitals(control_text: str) -> list[str]:
    """按出现顺序返回 ``output cube eigenstate <Q>`` 里的 ``Q``（HOMO-1/HOMO/…）。"""
    return [m.group(1) for m in (_ORBIT_RE.match(line) for line in control_text.splitlines())
            if m is not None]


def frontier_weights(stage_dir: str | Path, *,
                     anchor_elements: tuple[str, ...] = DEFAULT_ANCHOR_ELEMENTS,
                     radius: float = DEFAULT_ANCHOR_RADIUS) -> list[dict[str, object]]:
    """扫描阶段目录里的 cube，逐个给出锚定权重（按 cube 序号排序）。

    cube 与 ``control.in`` 里请求的顺序一一对应，所以标签取
    ``requested_orbitals()``；文件数与请求数不符时标签留 ``None``（只报本征态号）。
    没有任何 cube 时返回空列表——旧 run 没有 cube 不该报错。
    """
    stage_dir = Path(stage_dir)
    cubes = [p for p in stage_dir.glob("cube_*.cube") if _CUBE_NAME_RE.match(p.name)]
    if not cubes:
        return []

    def sort_key(p: Path) -> tuple[int, int]:
        match = _CUBE_NAME_RE.match(p.name)
        assert match is not None
        return int(match.group(1)), int(match.group(3))

    control = stage_dir / "control.in"
    labels: list[str | None] = []
    if control.exists():
        labels = list(requested_orbitals(
            control.read_text(encoding="utf-8", errors="replace")))
    aligned = len(labels) == len(cubes)

    report: list[dict[str, object]] = []
    for order, path in enumerate(sorted(cubes, key=sort_key)):
        match = _CUBE_NAME_RE.match(path.name)
        assert match is not None
        cube = read_cube(path)
        weights = anchor_weight(cube, anchor_elements=anchor_elements, radius=radius)
        report.append({
            "file": path.name,
            "orbit": labels[order] if aligned else None,
            "eigenstate": int(match.group(2)),
            "spin": int(match.group(3)),
            "atom_count": len(cube.atoms),
            "signed_cube": weights["signed_cube"],
            "radius_A": weights["radius_A"],
            "anchor_fraction": weights["fraction"],
            "anchors": weights["anchors"],
        })
    return report


#: 原子序数 → 元素符号（cube 原子行只给 Z；覆盖到 Xe + 常见重元素）。
_SYMBOLS: dict[int, str] = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F",
    10: "Ne", 11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S",
    17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 21: "Sc", 22: "Ti", 23: "V",
    24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu", 30: "Zn",
    31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br", 36: "Kr", 37: "Rb",
    38: "Sr", 39: "Y", 40: "Zr", 41: "Nb", 42: "Mo", 43: "Tc", 44: "Ru",
    45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn", 51: "Sb",
    52: "Te", 53: "I", 54: "Xe", 78: "Pt", 79: "Au", 80: "Hg",
}
