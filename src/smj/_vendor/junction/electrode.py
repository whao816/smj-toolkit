"""单分子结电极构建:Electrode111 晶格常数、四元数贴装、表面索引、append_au60。

来源: 第三方工具 qmtransport/electrode.py(145 行, tmp/自动化/, 只读来源不删除不入库)
按本项目协议适配迁入(AGD-2026-09-22-SMJ, P2):保留全部常数与公式 ——
delta_z=2.35454990 / delta_y=1.66491820 / delta_x=1.44186146(6 层体 + 4 额外原子);
123 原子复算已冻结为回归锚点(agenticdft/tests/data/junction_anchors.json, P0)。
纯本地工具: 只做结构生成, 不做真实计算。
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
import math
from .geometry import Atom, add, cross, dot, norm, read_aims_geometry, scale, sub, unit, write_aims_geometry
from .surface_index import pick_surface_triple as _pick_surface_triple


Vec = tuple[float, float, float]


def _qmul(q1: tuple[float, float, float, float], q2: tuple[float, float, float, float]) -> tuple[float, float, float, float]:
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return (
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
    )


def _qconj(q: tuple[float, float, float, float]) -> tuple[float, float, float, float]:
    return (q[0], -q[1], -q[2], -q[3])


def rotation_quaternion(source: Vec, target: Vec) -> tuple[float, float, float, float]:
    s = unit(source)
    t = unit(target)
    d = max(-1.0, min(1.0, dot(s, t)))
    if d > 0.999999:
        return (1.0, 0.0, 0.0, 0.0)
    if d < -0.999999:
        axis = unit(cross(s, (1.0, 0.0, 0.0)) if abs(s[0]) < 0.9 else cross(s, (0.0, 1.0, 0.0)))
        return (0.0, axis[0], axis[1], axis[2])
    axis = unit(cross(s, t))
    theta = math.acos(d)
    return (math.cos(theta / 2.0), *(scale(axis, math.sin(theta / 2.0))))


def rotate(v: Vec, q: tuple[float, float, float, float]) -> Vec:
    out = _qmul(_qmul(q, (0.0, *v)), _qconj(q))
    return (out[1], out[2], out[3])


@dataclass
class Electrode111:
    layers: int = 6
    extras: int = 4
    delta_z: float = 2.35454990
    delta_y: float = 1.66491820
    delta_x: float = 1.44186146
    au_per_electrode: int | None = None  # 每侧 Au 总数(含 tip 重合原子); None=全量(60)

    def body(self) -> list[Vec]:
        coords: list[Vec] = [(0.0, 0.0, 0.0)]
        for z in range(1, self.layers):
            cz0 = self.delta_z * z
            cx0 = -self.delta_x * z
            cy0 = -0.5 * self.delta_y * z
            for y in range(z + 1):
                for x in range(y, z + 1):
                    coords.append((cx0 + self.delta_x * (2 * x - y), cy0 + self.delta_y * 1.5 * y, cz0))
        return coords

    def extra_pool(self) -> list[Vec]:
        pool: list[Vec] = []
        for z in range(max(2, self.layers - 4), self.layers - 1):
            cz0 = self.delta_z * z
            cx0 = -self.delta_x * z
            cy0 = -0.5 * self.delta_y * z
            for i in range(z):
                for x, y in ((i, -1), (i, i + 1), (z + 1, i + 1)):
                    pool.append((cx0 + self.delta_x * (2 * x - y), cy0 + self.delta_y * 1.5 * y, cz0))
        return pool

    def new_electrode(self) -> list[Vec]:
        full = self.body() + self.extra_pool()[: self.extras]
        if self.au_per_electrode is None or self.au_per_electrode >= len(full):
            return full
        # 按层截取: body 第 z 层原子数 = 三角数 (z+1)(z+2)/2 (z=0 层 1 个, z=1 层 3 个, ...)
        # 逐层累加至 >= au_per_electrode, 保持每层完整(表面层不残缺, 供非共线三角选择)。
        total = 0
        cut = len(full)
        for z in range(self.layers):
            total += (z + 1) * (z + 2) // 2
            if total >= self.au_per_electrode:
                cut = total
                break
        return full[:cut]


def nearest_linker(atoms: list[Atom], au_index: int, tip_element: str = "Au") -> int:
    au = atoms[au_index].coord()
    best_index = -1
    best_distance = float("inf")
    for i, atom in enumerate(atoms):
        if atom.element in {tip_element, "H"}:
            continue
        d = norm(sub(atom.coord(), au))
        if d < best_distance:
            best_distance = d
            best_index = i
    if best_index < 0:
        raise ValueError("Cannot find non-H linker atom for Au tip")
    return best_index


# ---------------------------------------------------------------------------
# P2-1 (AGD-2026-09-25-SMJ): 电极参数化 —— 元素与晶格常数
#
# 勘察结论: Electrode111 已参数化 layers/extras/au_per_electrode, **真正硬编码的只有
# 晶格常量 delta_z/y/x**(= Au(111)); 而"用户给定电极"要求能表达其它金属/晶面。
# 晶格导出关系(数值已验证, 见 tests/test_smj_automation.py):
#   delta_z = a/sqrt(3);  delta_y = delta_z/sqrt(2);  delta_x = delta_z*sqrt(3)/(2*sqrt(2))
#   a(Au) = 4.07823 A 时三者与 Electrode111 的字面常量在 8 位小数内一致。
#
# **刻意保守**: 不传 a_lattice 时**完全继承** Electrode111 的字面常量与几何 ->
# 默认路径逐字节不变(由 test_au60_build_golden_regression 与既有 123/23 原子冻结测试守住)。
# 只有显式给出晶格常数时才按公式重算 —— 避免浮点末位差异打破黄金回归。
# ---------------------------------------------------------------------------
@dataclass
class ElectrodeFCC(Electrode111):
    """通用 fcc(111) 电极构造器: 参数化**元素**与**晶格常数**(P2-1)。

    element: 写入 geometry.in 的元素符号(默认 'Au', 与改动前一致)。
    a_lattice: 晶格常数(A); None 表示沿用 Electrode111 的字面 delta(逐字节不变)。
    """
    element: str = 'Au'
    a_lattice: float | None = None

    def __post_init__(self):
        if self.a_lattice:
            dz = float(self.a_lattice) / math.sqrt(3.0)
            self.delta_z = dz
            self.delta_y = dz / math.sqrt(2.0)
            self.delta_x = dz * math.sqrt(3.0) / (2.0 * math.sqrt(2.0))


# 常用金属的 fcc 晶格常数(A) —— 供上层给出"金属"而非"晶格常数"时换算。
FCC_LATTICE_A = {'Au': 4.07823, 'Ag': 4.08530, 'Cu': 3.61496, 'Ni': 3.52360,
                 'Pt': 3.92360, 'Pd': 3.89070, 'Al': 4.04950}


def electrode_for_metal(element: str, layers: int = 6, extras: int = 4,
                        au_per_electrode: int | None = None) -> ElectrodeFCC:
    """按金属名构造 fcc(111) 电极; 未知金属抛错(不猜晶格常数)。"""
    a = FCC_LATTICE_A.get(element)
    if a is None:
        raise ValueError('未知金属: %r(已登记: %s)' % (element, sorted(FCC_LATTICE_A)))
    return ElectrodeFCC(element=element, layers=layers, extras=extras,
                        au_per_electrode=au_per_electrode, a_lattice=a)


def mount_electrode(tip: Atom, linker: Atom, constructor: Electrode111) -> list[Atom]:
    bond_to_molecule = sub(linker.coord(), tip.coord())
    q = rotation_quaternion((0.0, 0.0, -1.0), bond_to_molecule)
    mounted = []
    for coord in constructor.new_electrode()[1:]:
        pos = add(rotate(coord, q), tip.coord())
        mounted.append(Atom(getattr(constructor, 'element', 'Au'), *pos))
    return mounted


def surface_indices(atoms: list[Atom], electrode_indices: list[int], tip_index: int, linker_index: int) -> list[int]:
    """选电极表面(沿 tip→linker 方向最远投影层)的 3 个 Au 作 aitranss 表面。

    W0-3(2026-11-05): 选择规则迁到 :mod:`surface_index`(确定性 det-v1),
    这里的签名只为兼容既有调用方。任何"同一几何两次运行给不同索引"的怀疑,
    都应看 :func:`surface_index.pick_surface_triple` 的 ``tie_count``。
    """
    return _pick_surface_triple(atoms, sorted(set(int(i) for i in electrode_indices)),
                                tip_index, linker_index=linker_index).indices


def append_au60(input_path: str | Path, output_path: str | Path,
                metadata_path: str | Path | None = None,
                au_per_electrode: int | None = None,
                element: str = "Au", a_lattice: float | None = None,
                allowed_root: str | Path | None = None,
                layers: int | None = None,
                extras: int | None = None) -> dict:
    """贴装电极并落盘 geometry.in / electrode_metadata.json; 返回 metadata dict。

    ROB-08(2026-10-19 P2): allowed_root 给定时对**两个写点**(output_path 与
    metadata_path)做第二道路径门禁 —— 必须落在 allowed_root 内, 越界抛 ValueError。
    MCP 写类工具传 settings.run.work_root(第一道门禁在 mcp_server._write_path_error);
    直接调用方(stages 等)不传则行为完全不变。

    W0-2(2026-11-04): layers/extras 显式给出时才改电极厚度; 缺省(None)逐字节回旧行为。
    """
    if allowed_root is not None:
        _root = Path(allowed_root).resolve()
        for _label, _p in (("output_path", output_path), ("metadata_path", metadata_path)):
            if _p is None or str(_p).strip() == "":
                continue
            try:
                _ok = Path(_p).resolve().is_relative_to(_root)
            except Exception:
                _ok = False
            if not _ok:
                raise ValueError("ROB-08 路径门禁: %s 必须落在 %s 内: %s"
                                 % (_label, _root, _p))
    atoms = read_aims_geometry(input_path)
    au_indices = [i for i, atom in enumerate(atoms) if atom.element == element]
    if len(au_indices) != 2:
        raise ValueError(f"Expected exactly two {element} tip atoms, found {len(au_indices)}")
    # W0-2: 只有显式传了 layers/extras 才进入 kwargs(dataclass 默认仍是 6/4)
    thick = {}
    if layers is not None:
        thick["layers"] = layers
    if extras is not None:
        thick["extras"] = extras
    # P2-1: element='Au' 且未显式给晶格常数 -> 走 Electrode111 原路径(默认逐字节不变)
    if element == "Au" and a_lattice is None:
        constructor = Electrode111(au_per_electrode=au_per_electrode, **thick)
    else:
        _a = a_lattice if a_lattice is not None else FCC_LATTICE_A.get(element)
        if _a is None:
            raise ValueError('未知金属: %r(请显式给出 a_lattice; 已登记: %s)'
                             % (element, sorted(FCC_LATTICE_A)))
        constructor = ElectrodeFCC(element=element, a_lattice=_a,
                                   au_per_electrode=au_per_electrode, **thick)
    original_count = len(atoms)
    left_tip, right_tip = au_indices
    left_linker = nearest_linker(atoms, left_tip, tip_element=element)
    right_linker = nearest_linker(atoms, right_tip, tip_element=element)
    left_added = mount_electrode(atoms[left_tip], atoms[left_linker], constructor)
    left_added_indices = list(range(len(atoms), len(atoms) + len(left_added)))
    atoms.extend(left_added)
    right_added = mount_electrode(atoms[right_tip], atoms[right_linker], constructor)
    right_added_indices = list(range(len(atoms), len(atoms) + len(right_added)))
    atoms.extend(right_added)
    write_aims_geometry(output_path, atoms)
    # W0-3: 两侧各选一次, 索引与"同分等价三角个数"一起进 metadata —— 审计据此判断
    # 该选择是几何确定的, 还是对称退化后被规则打破的。
    _left_tri = _pick_surface_triple(atoms, sorted(set([left_tip] + left_added_indices)),
                                     left_tip, side='left', linker_index=left_linker)
    _right_tri = _pick_surface_triple(atoms, sorted(set([right_tip] + right_added_indices)),
                                      right_tip, side='right', linker_index=right_linker)
    metadata = {
        "input": str(input_path),
        "output": str(output_path),
        "au_per_electrode": au_per_electrode if au_per_electrode is not None else len(constructor.new_electrode()),
        "original_atom_count": original_count,
        "left_tip_index": left_tip + 1,
        "right_tip_index": right_tip + 1,
        "left_linker_index": left_linker + 1,
        "right_linker_index": right_linker + 1,
        "left_surface_indices": list(_left_tri.indices),
        "right_surface_indices": list(_right_tri.indices),
        # W0-3: 规则名与对称退化计数(>1 = 有多个几何等价的三角同分, 规则按坐标字典序取一个)
        "surface_rule": "det-v1",
        "surface_tie_counts": {"left": _left_tri.tie_count, "right": _right_tri.tie_count},
        # W0-2: 精确的左右归属(1-based, 含 tip)。厚电极的贴装轴随分子倾斜,
        # 远层原子可能离对侧 tip 更近 —— 按距离划分会错分, 这里直接给出真值。
        "left_electrode_indices": [left_tip + 1] + [i + 1 for i in left_added_indices],
        "right_electrode_indices": [right_tip + 1] + [i + 1 for i in right_added_indices],
        "final_atom_count": len(atoms),
    }
    if thick:
        # 仅在显式加厚时写入(既有调用方 metadata 逐字节不变)
        metadata["electrode_layers"] = getattr(constructor, "layers", None)
        metadata["electrode_extras"] = getattr(constructor, "extras", None)
    if metadata_path:
        Path(metadata_path).parent.mkdir(parents=True, exist_ok=True)
        Path(metadata_path).write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata