"""有限位移力常数:位移集、中心差分组装、声学求和规则、lead 力常数拼接。

几条刻意的取舍(方案 §5.1、§5.4、D11):

* **中心差分**而不是单边差分:单边差分的截断误差是 O(δ),中心是 O(δ²),
  且中心差分天然消掉偶次非谐项。代价是 6N 次 SCF 而不是 3N。
* **不在这里做 ASR 对角补齐**。Mosso 等人把对角块用 ``K_ii = -Σ_{j≠i} K_ij``
  补齐,这会让 ``acoustic_sum_residual`` 变成恒为零的空检查。本模块的做法是
  **先**用实测 K 算出残差(作为独立的力常数质量指标),**再**投影。
* lead FC 是**冻结入库**的(``lead_fc.npy`` + 带 sha256 的 ``lead_fc.json``),
  不与器件一起重算 —— 器件区的化学环境与体相 lead 不同,混在一起算会两头不靠。
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .dynamical import acoustic_sum_residual, acoustic_sum_residuals
from .units import EV_PER_SQRT

__all__ = [
    "LeadFC",
    "FullFCMeta",
    "acoustic_residual_of",
    "apply_asr",
    "asr_diagnostics",
    "assemble_fc",
    "assemble_full_fc",
    "lead_device_conflicts",
    "load_lead_fc",
    "make_displacements",
    "save_lead_fc",
]

_AXIS_NAMES = ("x", "y", "z")


def make_displacements(n_atoms: int, device_indices, amplitude_A: float) -> list[dict]:
    """生成 ``6 * N_dev`` 条位移记录,顺序固定(先原子、再轴、再符号)。

    固定顺序是为了让 ``disp_index`` 与 Slurm 数组下标一一对应 —— 数组作业
    重启、部分重跑时不会把 ``disp_0007`` 的力塞进 ``disp_0008`` 的槽位。
    符号顺序是 ``+`` 在前、``-`` 在后。
    """
    n_atoms = int(n_atoms)
    if n_atoms < 1:
        raise ValueError("n_atoms 必须 ≥1,得到 %d" % n_atoms)
    amplitude = float(amplitude_A)
    if not np.isfinite(amplitude) or amplitude <= 0.0:
        raise ValueError("位移幅度必须为正的有限值,得到 %r" % (amplitude_A,))
    indices = sorted({int(i) for i in device_indices})
    if not indices:
        raise ValueError("device_indices 不能为空")
    if indices[0] < 0 or indices[-1] >= n_atoms:
        raise ValueError(
            "device_indices 越界: 取值范围 0..%d,得到 %d..%d"
            % (n_atoms - 1, indices[0], indices[-1]))
    records: list[dict] = []
    slot = 0
    for atom in indices:
        for axis in range(3):
            for sign in (+1, -1):
                records.append({
                    "slot": slot,
                    "atom": atom,
                    "axis": axis,
                    "axis_name": _AXIS_NAMES[axis],
                    "sign": sign,
                    "delta_A": sign * amplitude,
                    "label": "%04d_%s_%s" % (atom, _AXIS_NAMES[axis],
                                             "p" if sign > 0 else "m"),
                })
                slot += 1
    return records


def _device_atoms(disps: list[dict]) -> list[int]:
    atoms = sorted({int(d["atom"]) for d in disps})
    if not atoms:
        raise ValueError("位移集为空")
    return atoms


def assemble_fc(disps: list[dict], forces_by_disp, n_atoms: int) -> np.ndarray:
    """中心差分组装器件行,返回 ``(3*N_dev, 3*N)``,单位 eV/Å²。

    ``K[3i+a, 3j+b] = -(F_i^a(+δ) - F_i^a(-δ)) / (2δ)``,其中 ``i`` 枚举
    ``disps`` 里的器件原子(升序),``j`` 枚举全部原子。只有器件原子的行有值;
    lead 行由 :func:`assemble_full_fc` 从冻结的 lead 力常数填。
    """
    n_atoms = int(n_atoms)
    atoms = _device_atoms(disps)
    lookup = {atom: pos for pos, atom in enumerate(atoms)}
    expect = 6 * len(atoms)
    if len(disps) != expect:
        raise ValueError(
            "位移记录数应为 6*%d=%d,得到 %d" % (len(atoms), expect, len(disps)))
    if len(forces_by_disp) != len(disps):
        raise ValueError(
            "力的个数(%d)与位移记录数(%d)不一致"
            % (len(forces_by_disp), len(disps)))

    plus: dict[tuple[int, int], np.ndarray] = {}
    minus: dict[tuple[int, int], np.ndarray] = {}
    delta: dict[tuple[int, int], float] = {}
    for record, forces in zip(disps, forces_by_disp):
        key = (int(record["atom"]), int(record["axis"]))
        array = np.asarray(forces, dtype=float)
        if array.shape != (n_atoms, 3):
            raise ValueError(
                "力数组形状应为 (%d, 3),得到 %r" % (n_atoms, array.shape))
        step = float(record["delta_A"])
        if step == 0.0:
            raise ValueError("位移幅度为 0(原子 %d,%s 轴)" % key)
        target = plus if step > 0 else minus
        if key in target:
            raise ValueError("重复的位移记录: 原子 %d,%s 轴,%s" % (
                key[0], _AXIS_NAMES[key[1]], "正" if step > 0 else "负"))
        target[key] = array
        delta[key] = abs(step)

    fc = np.zeros((3 * len(atoms), 3 * n_atoms), dtype=float)
    for key in sorted(plus):
        if key not in minus:
            raise ValueError(
                "位移原子 %d,%s 轴缺少负向记录 —— 中心差分必须成对"
                % (key[0], _AXIS_NAMES[key[1]]))
        row_atom, axis = key
        row = 3 * lookup[row_atom] + axis
        fc[row, :] = -(plus[key] - minus[key]).reshape(-1) / (2.0 * delta[key])
    return fc


def _rebuild_force_constants(dynamical_ev2, masses) -> np.ndarray:
    matrix = np.asarray(dynamical_ev2, dtype=float)
    mass = np.asarray(masses, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("dynamical_matrix 必须是方阵,得到 %r" % (matrix.shape,))
    n = matrix.shape[0]
    reps = n // 3
    if n % 3 or mass.shape != (reps,):
        raise ValueError(
            "质量数组长度应为 %d,得到 %d" % (reps, mass.shape[0] if mass.ndim else 1))
    if not np.all(np.isfinite(matrix)):
        raise ValueError("dynamical_matrix 含非有限数值")
    if not np.all(np.isfinite(mass)) or np.any(mass <= 0):
        raise ValueError("质量必须是有限正数")
    scale = np.repeat(mass, 3)
    return matrix * np.sqrt(np.outer(scale, scale)) / (EV_PER_SQRT ** 2)


def _translation_basis(mass: np.ndarray) -> np.ndarray:
    return (np.sqrt(mass)[:, None, None] * np.eye(3)).reshape(-1, 3)


def apply_asr(dynamical_ev2, masses) -> tuple[np.ndarray, float]:
    """把三条质量加权平移模投影掉,返回 ``(D_proj, residual)``。

    ``residual`` 是投影**之前**在力常数矩阵上算的 ``max_iαβ |Σ_j K_iα,jβ|``
    (eV/Å²)。先测残差再投影,否则这个指标会退化成恒为零 —— 详见方案 §5.1。
    """
    matrix = np.asarray(dynamical_ev2, dtype=float)
    mass = np.asarray(masses, dtype=float)
    force_constants = _rebuild_force_constants(matrix, mass)
    residual = acoustic_sum_residual(force_constants)
    basis = _translation_basis(mass)
    complement = np.eye(matrix.shape[0]) - basis @ basis.T / mass.sum()
    projected = complement @ matrix @ complement
    return projected, residual


def asr_diagnostics(dynamical_ev2, masses, *, threshold: float = 1.0e-3,
                    rows=None) -> dict:
    """先测 ASR 残差,超过门限再投影,并把两边的数都记下来(D13)。

    返回的 dict 同时带"投影前"和"投影后"两套指标:

    ``residual_ev_per_A2`` / ``relative_residual``
        投影前在**力常数**矩阵上算的 ``max_iαβ |Σ_j K_iα,jβ|`` 及其相对值
        (``/ max|K|``)。判据 ``ph_asr_residual`` 判的就是这两个数 ——
        投影之后行和恒为机器零,拿投影后的值判等于没判。
    ``projected`` / ``projection_norm``
        是否做了投影,以及相对 Frobenius 投影量
        ``||D_proj - D|| / ||D||``(投影把三条平移模抹掉,量级应当很小)。
    ``residual_after_ev_per_A2``
        投影后重算的残差,用来证明投影确实生效。

    ``threshold`` 是**相对**门限:有限差分与基组噪声本身就在 1e-3~1e-1
    eV/Å² 的绝对量级,拿绝对门限判会让每个 run 都 inconclusive。

    ``rows`` 是自由度下标,限定在**邻居集合完整**的区域上(拼接矩阵里只有
    器件区满足;lead 行来自被截断成两个主层的 lead 力常数,行和按构造不闭合)。
    给了 ``rows`` 时,只修改器件子块,计入固定区外耦合的平移响应。
    如果无法用对称器件块抵消该响应,返回原矩阵并记录 projection_reason。
    """
    matrix = np.asarray(dynamical_ev2, dtype=float)
    mass = np.asarray(masses, dtype=float)
    force_constants = _rebuild_force_constants(matrix, mass)

    index = None if rows is None else np.asarray(rows, dtype=int).ravel()
    if index is not None:
        index = np.unique(index)
    residual, relative = acoustic_sum_residuals(force_constants, rows=index)

    report = {
        "residual_ev_per_A2": residual,
        "relative_residual": relative,
        "projected": False,
        "projection_reason": None,
        "projection_norm": None,
        "residual_after_ev_per_A2": None,
        "dynamical": matrix,
    }
    if relative <= threshold:
        return report

    if index is None:
        projected, _ = apply_asr(matrix, mass)
    else:
        basis = _translation_basis(mass)
        outside = np.setdiff1d(np.arange(matrix.shape[0]), index)
        device_basis = basis[index]
        inverse = np.linalg.pinv(device_basis)
        complement = np.eye(index.size) - device_basis @ inverse
        # For D_cc T_c = -D_cl T_l, the action on translations is prescribed.
        # Keep the optical-optical block and complete the cross terms by transpose.
        required = -(matrix[np.ix_(index, outside)] @ basis[outside]) @ inverse
        block = matrix[np.ix_(index, index)]
        sub = complement @ block @ complement + required + required.T @ complement
        projected = matrix.copy()
        projected[np.ix_(index, index)] = (sub + sub.T) * 0.5
        # An asymmetric translation-translation target admits no reciprocal solution.
        error = acoustic_residual_of(projected, mass, rows=index)
        scale = float(np.max(np.abs(force_constants[index, :])))
        if error > 32 * np.finfo(float).eps * matrix.shape[0] * scale:
            report["projection_reason"] = "incompatible_fixed_coupling"
            return report
    norm = float(np.linalg.norm(matrix))
    report.update(
        projected=True,
        projection_norm=float(np.linalg.norm(projected - matrix) / norm) if norm else None,
        residual_after_ev_per_A2=acoustic_residual_of(projected, mass, rows=index),
        dynamical=projected,
    )
    return report


def acoustic_residual_of(dynamical_ev2, masses, *, rows=None) -> float:
    """投影后的残差:把 ``D`` 还原成 ``K`` 再取行和的最大绝对值(eV/Å²)。"""
    force_constants = _rebuild_force_constants(dynamical_ev2, masses)
    return acoustic_sum_residual(force_constants, rows)



@dataclass(frozen=True)
class LeadFC:
    """冻结的体相 lead 力常数(2D 周期板,``layers`` 个原子层 × ``in_plane`` 胞)。"""

    element: str
    a_lattice_A: float
    n_atoms_per_layer: int
    layers: int
    in_plane: tuple[int, int]
    fc_cutoff_A: float
    force_constants: np.ndarray
    sha256: str
    source: str = ""

    @property
    def n_atoms(self) -> int:
        return self.n_atoms_per_layer * self.layers


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def save_lead_fc(directory, force_constants, *, element: str, a_lattice_A: float,
                 n_atoms_per_layer: int, layers: int, in_plane,
                 fc_cutoff_A: float, source: str = "") -> LeadFC:
    """写 ``lead_fc.npy`` + ``lead_fc.json``(含 sha256)。仅供生成方使用。"""
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=True)
    array = np.ascontiguousarray(np.asarray(force_constants, dtype=float))
    n = int(n_atoms_per_layer) * int(layers)
    if array.shape != (3 * n, 3 * n):
        raise ValueError("lead 力常数形状应为 (%d, %d),得到 %r"
                         % (3 * n, 3 * n, array.shape))
    npy = root / "lead_fc.npy"
    np.save(npy, array)
    digest = _sha256_file(npy)
    payload = {
        "schema": "smj.lead_fc.v1",
        "element": element,
        "a_lattice_A": float(a_lattice_A),
        "n_atoms_per_layer": int(n_atoms_per_layer),
        "layers": int(layers),
        "in_plane": [int(v) for v in in_plane],
        "fc_cutoff_A": float(fc_cutoff_A),
        "units": "eV/Ang^2",
        "force_constants_file": npy.name,
        "sha256": digest,
        "source": source,
    }
    (root / "lead_fc.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return LeadFC(
        element=str(element), a_lattice_A=float(a_lattice_A),
        n_atoms_per_layer=int(n_atoms_per_layer), layers=int(layers),
        in_plane=(int(in_plane[0]), int(in_plane[1])),
        fc_cutoff_A=float(fc_cutoff_A), force_constants=array,
        sha256=digest, source=source)


def load_lead_fc(path) -> LeadFC:
    """读冻结 lead 力常数并校验 sha256;任何不符都 raise,不做"尽力而为"。"""
    root = Path(path)
    if root.is_dir():
        root = root / "lead_fc.json"
    if not root.exists():
        raise FileNotFoundError("找不到 lead 力常数元数据: %s" % root)
    payload = json.loads(root.read_text(encoding="utf-8"))
    for key in ("element", "a_lattice_A", "n_atoms_per_layer", "layers",
                "in_plane", "fc_cutoff_A", "force_constants_file", "sha256"):
        if key not in payload:
            raise ValueError("lead 力常数元数据缺少字段 %r: %s" % (key, root))
    npy = root.parent / payload["force_constants_file"]
    if not npy.exists():
        raise FileNotFoundError("找不到 lead 力常数数组: %s" % npy)
    actual = _sha256_file(npy)
    if actual != payload["sha256"]:
        raise ValueError(
            "lead 力常数 sha256 不符: 元数据 %s,实测 %s (%s)"
            % (payload["sha256"][:16], actual[:16], npy))
    array = np.load(npy)
    n = int(payload["n_atoms_per_layer"]) * int(payload["layers"])
    if array.shape != (3 * n, 3 * n):
        raise ValueError("lead 力常数形状应为 (%d, %d),得到 %r"
                         % (3 * n, 3 * n, array.shape))
    if not np.all(np.isfinite(array)):
        raise ValueError("lead 力常数含非有限数值: %s" % npy)
    return LeadFC(
        element=str(payload["element"]),
        a_lattice_A=float(payload["a_lattice_A"]),
        n_atoms_per_layer=int(payload["n_atoms_per_layer"]),
        layers=int(payload["layers"]),
        in_plane=(int(payload["in_plane"][0]), int(payload["in_plane"][1])),
        fc_cutoff_A=float(payload["fc_cutoff_A"]),
        force_constants=array,
        sha256=str(payload["sha256"]),
        source=str(payload.get("source", "")))


@dataclass(frozen=True)
class FullFCMeta:
    """把器件行/lead 行映射到完整 ``(3N, 3N)`` 矩阵所需的索引。"""

    n_atoms: int
    device_indices: tuple[int, ...]
    left_indices: tuple[int, ...]
    right_indices: tuple[int, ...]

    def __post_init__(self):
        n = int(self.n_atoms)
        seen: set[int] = set()
        for name in ("device_indices", "left_indices", "right_indices"):
            values = tuple(int(v) for v in getattr(self, name))
            if any(v < 0 or v >= n for v in values):
                raise ValueError("%s 越界(0..%d)" % (name, n - 1))
            if len(set(values)) != len(values):
                raise ValueError("%s 有重复下标" % name)
            overlap = seen & set(values)
            if overlap:
                raise ValueError(
                    "器件区与 lead 区下标重叠: %s" % sorted(overlap)[:5])
            seen |= set(values)


def _place_lead_block(fc: np.ndarray, lead: LeadFC, indices, *, reflect_z=False) -> None:
    indices = tuple(indices)
    expected = lead.n_atoms
    if len(indices) != expected:
        raise ValueError(
            "lead 区原子数 %d 与 lead 力常数原子数 %d 不符"
            % (len(indices), expected))
    block = lead.force_constants
    if reflect_z:
        # device.build_device mirrors z for every right-hand atom. Both Cartesian
        # indices of K transform with Q = diag(1, 1, -1) under that reflection.
        signs = np.tile([1., 1., -1.], expected)
        block = block * signs[:, None] * signs[None, :]
    for a_pos, atom_a in enumerate(indices):
        for b_pos, atom_b in enumerate(indices):
            fc[3 * atom_a:3 * atom_a + 3, 3 * atom_b:3 * atom_b + 3] = \
                block[3 * a_pos:3 * a_pos + 3, 3 * b_pos:3 * b_pos + 3]


def assemble_full_fc(k_dev: np.ndarray, lead: LeadFC, meta: FullFCMeta) -> np.ndarray:
    """拼出完整 ``(3N, 3N)`` 力常数矩阵,单位 eV/Å²。

    lead 区四块来自冻结 lead 力常数;器件行随后被器件侧实测值**整体覆盖**
    (交界处优先信实测)。覆盖造成的差异用 :func:`lead_device_conflicts` 单独统计。

    **器件↔lead 的交叉块必须补上转置**:器件行给出 ``K[device, lead]``,
    而 lead 力常数只覆盖 lead 原子本身,``K[lead, device]`` 不会被填。
    互易性要求两者互为转置;漏掉时矩阵不对称度约 8%,后果有两个:
    ① ASR 残差从 1e-13 涨到 2.7 eV/Å²(看似"力常数不合声学求和规则",
    实际是矩阵缺块);② ``dynamical_matrix(symmetrize=True)`` 把交叉块
    平均成一半,``v_l``(器件↔主层耦合)直接少一半,而透射按 ``|v|⁴`` 变化 ——
    弱耦合处可差一个数量级,强耦合处因共振饱和不易察觉,属"看不出错但数值错"
    的一类(2026-11-25 实测:补转置后 ASR 残差从 2.7 降到 2.6e-13)。
    """
    n = int(meta.n_atoms)
    matrix = np.zeros((3 * n, 3 * n), dtype=float)
    _place_lead_block(matrix, lead, meta.left_indices)
    _place_lead_block(matrix, lead, meta.right_indices, reflect_z=True)

    array = np.asarray(k_dev, dtype=float)
    expected = (3 * len(meta.device_indices), 3 * n)
    if array.shape != expected:
        raise ValueError("器件行形状应为 %r,得到 %r" % (expected, array.shape))
    for pos, atom in enumerate(meta.device_indices):
        matrix[3 * atom:3 * atom + 3, :] = array[3 * pos:3 * pos + 3, :]

    lead_atoms = list(meta.left_indices) + list(meta.right_indices)
    if lead_atoms:
        device_dofs = _dofs(meta.device_indices)
        lead_dofs = _dofs(lead_atoms)
        matrix[np.ix_(lead_dofs, device_dofs)] = \
            matrix[np.ix_(device_dofs, lead_dofs)].T
    return matrix


def _dofs(atoms) -> np.ndarray:
    """把原子下标展开成自由度下标 ``[3a, 3a+1, 3a+2]``。"""
    return np.array([3 * int(a) + k for a in atoms for k in range(3)],
                    dtype=int)


def lead_device_conflicts(k_dev: np.ndarray, lead: LeadFC, meta: FullFCMeta,
                          *, tol: float = 1.0e-3) -> int:
    """统计器件行覆盖 lead 力常数时``|Δ| > tol`` 的元素个数(只查 lead 列)。

    不为零说明器件交界处的力常数与体相 lead 确有差异 —— 这是**预期**的,
    不是错误;报告里应把它当作"界面与体相不一致"的量化,而不是吞掉。
    """
    reference = assemble_full_fc(
        np.zeros((3 * len(meta.device_indices), 3 * meta.n_atoms)),
        lead, meta)
    device_rows = assemble_full_fc(np.asarray(k_dev, dtype=float), lead, meta)
    lead_columns = list(meta.left_indices) + list(meta.right_indices)
    if not lead_columns:
        return 0
    columns = np.concatenate([np.arange(3 * a, 3 * a + 3) for a in lead_columns])
    count = 0
    for atom in meta.device_indices:
        rows = np.arange(3 * atom, 3 * atom + 3)
        delta = np.abs(device_rows[np.ix_(rows, columns)]
                       - reference[np.ix_(rows, columns)])
        count += int(np.count_nonzero(delta > tol))
    return count
