"""周期板器件构造:Au(111) 3×3 × layers + 分子,按角色划分器件区/lead 区。

为什么不复用电子链的金字塔团簇(方案 §4.1):金字塔团簇没有面内周期性,
无法定义横向 k 采样,也无法给半无限 lead 一个确定的表面元胞;声子 NEGF 需要的是
"二维周期板 + 沿 z 的层堆叠",而不是一个尖朝分子的原子堆。

坐标约定:

* ``z`` 是传输方向,分子两个 S 沿 ``z`` 对齐,中点落在 ``z=0``;
* 左侧电极在 ``z<0``,右侧在 ``z>0``;
* ``layer_index`` 每侧独立编号,``0 = 最外层(远离分子)``,向分子递增。

**第一步必须是** :func:`smj.junction.geometry.strip_thiol_h`:``molecule_opt``
产出的是中性二硫醇(含 S–H),直接贴 Au 会得到"S 上还挂着 H、同时又离 Au 2.35 Å"
的非物理接触 —— 报告 §6 陷阱⑥,也是 T5b 存在的理由。
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path

from ..junction.geometry import Atom, strip_thiol_h, write_aims_geometry

__all__ = [
    "DeviceMeta",
    "DeviceSpec",
    "build_device",
    "fcc111_in_plane_vectors",
    "write_device",
]

DELTA_Z_PER_LAYER = 2.35454990  # a/sqrt(3),与 junction/electrode.py 同一约定


def fcc111_in_plane_vectors(a_lattice_A: float, in_plane=(3, 3)
                            ) -> tuple[tuple[float, float, float],
                                       tuple[float, float, float]]:
    """(111) 面内二维胞矢量,单位 Å。最近邻 = ``a/sqrt(2)``。"""
    nearest = float(a_lattice_A) / math.sqrt(2.0)
    a1 = (in_plane[0] * nearest, 0.0, 0.0)
    a2 = (in_plane[1] * nearest * 0.5,
          in_plane[1] * nearest * math.sqrt(3.0) / 2.0, 0.0)
    return a1, a2


@dataclass(frozen=True)
class DeviceSpec:
    """器件几何参数。M3 的 ``ThermalConfig`` 会映射到这个对象。"""

    a_lattice_A: float = 4.078230
    in_plane: tuple[int, int] = (3, 3)
    electrode_layers: int = 4
    device_layers: int = 2
    vacuum_A: float = 20.0
    au_s_distance_A: float = 2.35
    gap_A: float = 8.4
    element: str = "Au"
    anchor_element: str = "S"
    strip_thiol_h: bool = True
    freeze_outer_layers: int = 1

    def __post_init__(self):
        if not 1 <= self.device_layers <= self.electrode_layers - 2:
            raise ValueError(
                "每侧需要至少一层器件和两层冻结 lead：electrode_layers=%d, device_layers=%d"
                % (self.electrode_layers, self.device_layers))
        if self.vacuum_A <= 0:
            raise ValueError("vacuum_A 必须为正")
        if self.in_plane[0] < 1 or self.in_plane[1] < 1:
            raise ValueError("in_plane 两个方向都必须 ≥1")


@dataclass
class DeviceMeta:
    schema: str
    n_atoms: int
    lattice: tuple[tuple[float, float, float], ...]
    roles: list[str]
    layer_index: list[int]
    side: list[str]
    molecule_indices: list[int]
    device_indices: list[int]
    lead_indices: list[int]
    left_device_indices: list[int]
    right_device_indices: list[int]
    left_lead_indices: list[int]
    right_lead_indices: list[int]
    left_surface_indices: list[int]
    right_surface_indices: list[int]
    constrained_indices: list[int]
    anchors: dict
    sbond: list[float]
    spec: dict = field(default_factory=dict)
    removed_h: list[int] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "schema": self.schema,
            "n_atoms": self.n_atoms,
            "lattice": [list(v) for v in self.lattice],
            "roles": list(self.roles),
            "layer_index": list(self.layer_index),
            "side": list(self.side),
            "molecule_indices": list(self.molecule_indices),
            "device_indices": list(self.device_indices),
            "lead_indices": list(self.lead_indices),
            "left_device_indices": list(self.left_device_indices),
            "right_device_indices": list(self.right_device_indices),
            "left_lead_indices": list(self.left_lead_indices),
            "right_lead_indices": list(self.right_lead_indices),
            "left_surface_indices": list(self.left_surface_indices),
            "right_surface_indices": list(self.right_surface_indices),
            "constrained_indices": list(self.constrained_indices),
            "anchors": dict(self.anchors),
            "sbond": list(self.sbond),
            "removed_h": list(self.removed_h),
            "spec": dict(self.spec),
        }


def _rotation_to_z(atoms: list[Atom], left: int, right: int) -> list[Atom]:
    """把 ``right - left`` 旋到 ``+z``:Rodrigues 公式,不引入 scipy。"""
    a = atoms[left]
    b = atoms[right]
    axis = (b.x - a.x, b.y - a.y, b.z - a.z)
    norm = math.sqrt(sum(v * v for v in axis))
    if norm < 1e-8:
        raise ValueError("两个锚定原子重合,无法定义分子轴")
    target = (0.0, 0.0, 1.0)
    u = tuple(v / norm for v in axis)
    # 旋转轴 = u × target;若已平行则无需旋转
    cross = (u[1] * target[2] - u[2] * target[1],
             u[2] * target[0] - u[0] * target[2],
             u[0] * target[1] - u[1] * target[0])
    sine = math.sqrt(sum(v * v for v in cross))
    cosine = sum(u[i] * target[i] for i in range(3))
    if sine < 1e-12:
        if cosine > 0:
            return atoms
        # 反向:绕 x 轴转 180°
        return [Atom(at.element, at.x, -at.y, -at.z) for at in atoms]
    k = (cross[0] / sine, cross[1] / sine, cross[2] / sine)

    def rotate(atom: Atom) -> Atom:
        x, y, z = atom.coord()
        # v_rot = v cosθ + (k×v) sinθ + k (k·v)(1−cosθ)
        kx, ky, kz = k
        dot = kx * x + ky * y + kz * z
        cross_v = (ky * z - kz * y, kz * x - kx * z, kx * y - ky * x)
        return Atom(
            atom.element,
            x * cosine + cross_v[0] * sine + kx * dot * (1.0 - cosine),
            y * cosine + cross_v[1] * sine + ky * dot * (1.0 - cosine),
            z * cosine + cross_v[2] * sine + kz * dot * (1.0 - cosine),
        )

    return [rotate(atom) for atom in atoms]


def _centre(atoms: list[Atom], left: int, right: int) -> list[Atom]:
    mid = ((atoms[left].x + atoms[right].x) / 2.0,
           (atoms[left].y + atoms[right].y) / 2.0,
           (atoms[left].z + atoms[right].z) / 2.0)
    return [Atom(a.element, a.x - mid[0], a.y - mid[1], a.z - mid[2])
            for a in atoms]


def _in_plane_sites(a1, a2, in_plane, offset) -> list[tuple[float, float]]:
    """一层原子的面内位置(等间距网格 + ABC 堆垛偏移)。"""
    nx, ny = in_plane
    sites = []
    for i in range(nx):
        for j in range(ny):
            x = (i / nx) * a1[0] + (j / ny) * a2[0] + offset[0]
            y = (i / nx) * a1[1] + (j / ny) * a2[1] + offset[1]
            sites.append((x, y))
    return sites


def build_device(molecule_atoms: list[Atom], spec: DeviceSpec | None = None
                 ) -> tuple[list[Atom], DeviceMeta]:
    """构造器件:分子(去 S–H)+ 两侧 Au(111) 板,返回 ``(atoms, meta)``。"""
    spec = spec or DeviceSpec()
    molecule = list(molecule_atoms)
    removed: list[int] = []
    if spec.strip_thiol_h:
        molecule, info = strip_thiol_h(molecule, element=spec.anchor_element)
        removed = list(info.get("removed_h_indices", []))

    anchors = [i for i, a in enumerate(molecule) if a.element == spec.anchor_element]
    if len(anchors) != 2:
        raise ValueError(
            "去质子后应有恰好 2 个 %s 锚原子,实得 %d"
            % (spec.anchor_element, len(anchors)))
    left_anchor, right_anchor = anchors
    if molecule[left_anchor].z > molecule[right_anchor].z:
        left_anchor, right_anchor = right_anchor, left_anchor

    molecule = _rotation_to_z(molecule, left_anchor, right_anchor)
    molecule = _centre(molecule, left_anchor, right_anchor)
    s_left = molecule[left_anchor].coord()
    s_right = molecule[right_anchor].coord()

    a1, a2 = fcc111_in_plane_vectors(spec.a_lattice_A, spec.in_plane)
    delta_z = DELTA_Z_PER_LAYER * (spec.a_lattice_A / 4.078230)
    # ABC 堆垛按两条原胞基矢之和的 1/3 平移，与超胞的面内尺寸无关。
    nx, ny = spec.in_plane
    stack = tuple((a1[k] / nx + a2[k] / ny) / 3.0 for k in (0, 1))

    layers = spec.electrode_layers
    surface_z = spec.au_s_distance_A  # 表面 Au 平面到 S 的距离
    half_gap = spec.gap_A / 2.0

    atoms: list[Atom] = list(molecule)
    roles: list[str] = ["molecule"] * len(molecule)
    layer_index: list[int] = [-1] * len(molecule)
    side: list[str] = ["none"] * len(molecule)

    left_device: list[int] = []
    right_device: list[int] = []
    left_lead: list[int] = []
    right_lead: list[int] = []
    left_surface: list[int] = []
    right_surface: list[int] = []
    constrained: list[int] = []
    sbond: list[float] = []

    for which in ("left", "right"):
        # 平移量必须**每个电极只算一次**,用表面层的格点定:层与层之间的 ABC
        # 堆垛就体现在 offset 的差上,如果每层各自 snap 到 S 正下方,所有层会
        # 对齐到同一个原子(AA 堆垛),相邻层 Au-Au 只剩 delta_z=2.355 Å,
        # 而正常的最近邻应是 a/sqrt(2)=2.885 Å。2026-11-23 就是这样在 relax
        # 里出现 5 eV/Å 的初始力。
        anchor_xy = (s_left[0], s_left[1]) if which == "left" else \
            (s_right[0], s_right[1])
        surface_offset = (stack[0] * (layers - 1), stack[1] * (layers - 1))
        surface_sites = _in_plane_sites(a1, a2, spec.in_plane, surface_offset)
        nearest = min(surface_sites,
                      key=lambda s: (s[0] - anchor_xy[0]) ** 2
                      + (s[1] - anchor_xy[1]) ** 2)
        shift = (anchor_xy[0] - nearest[0], anchor_xy[1] - nearest[1])
        for layer in range(layers):
            # layer 0 = 最外层(远离分子);layers-1 = 表面层
            distance_from_surface = (layers - 1) - layer
            z = -(half_gap + surface_z) - distance_from_surface * delta_z
            if which == "right":
                z = -z
            offset = (stack[0] * layer, stack[1] * layer)
            sites = _in_plane_sites(a1, a2, spec.in_plane, offset)
            layer_start = len(atoms)
            for (x, y) in sites:
                atoms.append(Atom(spec.element, x + shift[0], y + shift[1], z))
                layer_index.append(layer)
                side.append(which)
                roles.append("device_au" if layer >= layers - spec.device_layers
                             else "lead_au")
            layer_indices = list(range(layer_start, len(atoms)))
            if layer == layers - 1:
                (left_surface if which == "left" else right_surface).extend(
                    layer_indices)
                for index in layer_indices:
                    atom = atoms[index]
                    anchor = s_left if which == "left" else s_right
                    sbond.append(math.dist(atom.coord(), anchor))
            if layer >= layers - spec.device_layers:
                (left_device if which == "left" else right_device).extend(
                    layer_indices)
            else:
                (left_lead if which == "left" else right_lead).extend(
                    layer_indices)
            if layer < spec.freeze_outer_layers:
                constrained.extend(layer_indices)

    # 器件区 = 分子 + 每侧靠近分子的 device_layers 层;位移集就取这个集合
    device_indices = sorted(list(range(len(molecule))) + left_device + right_device)
    lead_indices = sorted(left_lead + right_lead)
    z_values = [a.z for a in atoms]
    c_length = (max(z_values) - min(z_values)) + spec.vacuum_A
    lattice = (a1, a2, (0.0, 0.0, c_length))

    meta = DeviceMeta(
        schema="smj.thermal_device.v1",
        n_atoms=len(atoms),
        lattice=lattice,
        roles=roles,
        layer_index=layer_index,
        side=side,
        molecule_indices=list(range(len(molecule))),
        device_indices=device_indices,
        lead_indices=lead_indices,
        left_device_indices=left_device,
        right_device_indices=right_device,
        left_lead_indices=left_lead,
        right_lead_indices=right_lead,
        left_surface_indices=left_surface,
        right_surface_indices=right_surface,
        constrained_indices=sorted(constrained),
        anchors={"left": left_anchor, "right": right_anchor},
        sbond=sbond,
        removed_h=removed,
        spec={
            "a_lattice_A": spec.a_lattice_A,
            "in_plane": list(spec.in_plane),
            "electrode_layers": spec.electrode_layers,
            "device_layers": spec.device_layers,
            "vacuum_A": spec.vacuum_A,
            "au_s_distance_A": spec.au_s_distance_A,
            "gap_A": spec.gap_A,
            "element": spec.element,
            "anchor_element": spec.anchor_element,
            "strip_thiol_h": spec.strip_thiol_h,
            "freeze_outer_layers": spec.freeze_outer_layers,
            "delta_z_A": delta_z,
        },
    )
    return atoms, meta


def write_device(directory, atoms: list[Atom], meta: DeviceMeta
                 ) -> tuple[Path, Path]:
    """写 ``geometry.in``(带 lattice)与 ``device_metadata.json``。"""
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=True)
    geometry = root / "geometry.in"
    write_aims_geometry(geometry, atoms,
                        constrained=tuple(meta.constrained_indices),
                        lattice=meta.lattice)
    metadata = root / "device_metadata.json"
    metadata.write_text(
        json.dumps(meta.as_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8")
    return geometry, metadata
