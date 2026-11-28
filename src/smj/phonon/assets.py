"""随包发布的电极力常数资产（lead FC）。

为什么要有这个模块
------------------
``thermal_analyze`` 需要一份 Au(111) 半无限电极的力常数（``lead_fc.npy`` +
``lead_fc.json``）。它只能由一次体相 DFT 有限位移得到（3×3×6 层、Γ、54 个
SCF，见 ``scripts/gen_au_lead_fc.py``），**不是每次算结都要重做的量**：电极
是金，与分子无关。所以本仓把一份验证过的结果随包发布，让第二个使用者不必
先搭一套体相计算才能跑分子。

资产位置与用法
--------------
``thermal.lead_fc`` 的默认值是 ``"auto"``，解析到本模块的 ``lead_fc_dir()``；
写具体路径（目录或 ``lead_fc.json``）就覆盖它。``smj lead-fc --show`` 打印
这份资产的身份（sha256、晶格常数、面内胞、层数），``--export DIR`` 把它复制
出去做成体相计算或替换成自己算的版本。

已知边界
--------
* 资产是**体相** Au(111) 的头两个主层（ABC 堆垛），不是"板的表面层"。这与
  Pauly 组做法一致：decimation 描述半无限 lead，器件里那两层只是它的开头。
* 表面弛豫没有显式包含（体相几何）。对声子带顶的影响在 1% 量级（实测带顶
  19.85 meV vs 实验 19.5 meV），且远小于 η 带来的系统不确定性。
"""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import numpy as np

from .device import fcc111_in_plane_vectors

__all__ = [
    "BUNDLED_LEAD_FC_DIR",
    "lead_fc_dir",
    "lead_fc_info",
    "export_lead_fc",
    "resolve_lead_fc_path",
    "lead_boundary",
]

#: 随包发布的资产目录（``lead_fc.npy`` + ``lead_fc.json``）。
BUNDLED_LEAD_FC_DIR = Path(__file__).resolve().parent / "data" / "au111_3x3_2layer"


def lead_fc_dir() -> Path:
    """返回随包资产的目录。

    目录不存在时抛 ``FileNotFoundError`` 并给出补救办法——这通常意味着安装
    时漏了 ``package-data``（``pip`` 的 ``--no-binary`` 或手工拷贝时会踩到）。
    """
    if not (BUNDLED_LEAD_FC_DIR / "lead_fc.json").is_file():
        raise FileNotFoundError(
            "缺少随包的电极力常数资产 %s；重新安装 smj-toolkit，或用 "
            "smj lead-fc --export DIR 从一份已有结果复制，再在 smj.yaml 里"
            "用 thermal.lead_fc 指向它" % BUNDLED_LEAD_FC_DIR)
    return BUNDLED_LEAD_FC_DIR


def resolve_lead_fc_path(value: str | Path | None) -> Path:
    """Resolve the bundled asset, an exported directory or a metadata filename."""
    path = lead_fc_dir() if value is None or str(value) == "auto" else Path(value)
    return path / "lead_fc.json" if path.is_dir() else path


def read_device_structure(geometry: Path):
    elements, positions, lattice = [], [], []
    for line in geometry.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if parts[:1] == ["atom"]:
            elements.append(parts[4])
            positions.append([float(v) for v in parts[1:4]])
        elif parts[:1] == ["lattice_vector"]:
            lattice.append([float(v) for v in parts[1:4]])
    return elements, np.asarray(positions), lattice


def lead_boundary(lead, meta: dict, elements: list[str], lattice, cutoff):
    """Validate the bulk asset and return (outer, principal) for each electrode."""
    per = lead.n_atoms_per_layer
    if (lead.layers != 2 or min(lead.in_plane) < 1 or per != np.prod(lead.in_plane)
            or not np.isfinite([lead.a_lattice_A, lead.fc_cutoff_A]).all()
            or min(lead.a_lattice_A, lead.fc_cutoff_A) <= 0):
        raise ValueError("lead FC 资产须含两个完整原子层、有限正晶格常数及截断半径")
    if cutoff is not None and not np.isclose(cutoff, lead.fc_cutoff_A, rtol=0, atol=1.e-6):
        raise ValueError(f"thermal.fc.cutoff={cutoff} Å 与资产截断 {lead.fc_cutoff_A} Å 不符")

    spec = meta.get("spec") or {}
    if ("element" in spec and spec["element"] != lead.element
            or "in_plane" in spec and tuple(spec["in_plane"]) != lead.in_plane
            or "a_lattice_A" in spec and not np.isclose(
                spec["a_lattice_A"], lead.a_lattice_A, rtol=0, atol=1.e-6)):
        raise ValueError("lead FC 的元素、面内尺寸或晶格常数与器件规格不符")
    expected_cell = fcc111_in_plane_vectors(lead.a_lattice_A, lead.in_plane)
    for cell in (meta.get("lattice"), lattice):
        if cell is not None and len(cell):
            vectors = np.asarray(cell, dtype=float)
            if vectors.shape != (3, 3) or not np.allclose(vectors[:2], expected_cell,
                                                          atol=2.e-6, rtol=0):
                raise ValueError("lead FC 的面内超胞晶格矢量与器件不符")

    layers, sides = meta["layer_index"], meta["side"]
    n = meta["n_atoms"]
    indices = meta["lead_indices"]
    if (len(layers) != n or len(sides) != n or len(elements) != n
            or len(set(indices)) != len(indices)
            or any(i < 0 or i >= n or sides[i] not in ("left", "right")
                   or elements[i] != lead.element for i in indices)):
        raise ValueError("器件电极原子编号、元素或左右分区与 lead FC 不符")
    boundaries = []
    for side in ("left", "right"):
        frozen = sorted(i for i in indices if sides[i] == side)
        levels = sorted({layers[i] for i in frozen})
        if len(levels) < 2 or levels != list(range(len(levels))):
            raise ValueError(f"{side} 侧须有至少两个连续冻结层")
        if ("electrode_layers" in spec and "device_layers" in spec
                and len(levels) != spec["electrode_layers"] - spec["device_layers"]):
            raise ValueError(f"{side} 侧冻结层数与器件规格不符")
        # Atom IDs follow the generated in-plane i,j site order. build_device
        # advances the same unwrapped ABC shift from every layer to the next.
        pair = [[i for i in frozen if layers[i] == level] for level in levels[-2:]]
        if any(len(group) != per for group in pair):
            raise ValueError(f"{side} 侧主层或外邻原子数与 lead FC 不符")
        boundaries.append(pair)
    return boundaries


def lead_fc_info(directory: str | Path | None = None) -> dict:
    """读 ``lead_fc.json`` 并**重新校验 sha256**，返回元数据。

    校验放在这里而不是只信任 json 里的字段：力常数被截断/被 git 的 CRLF
    转换动过时，只有实算哈希才能发现。
    """
    path = Path(directory) if directory is not None else lead_fc_dir()
    meta_path = path / "lead_fc.json" if path.is_dir() else path
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    npy_path = meta_path.parent / str(meta.get("force_constants_file", "lead_fc.npy"))
    digest = hashlib.sha256(npy_path.read_bytes()).hexdigest()
    expected = str(meta.get("sha256", ""))
    if expected and digest != expected:
        raise ValueError(
            "电极力常数 sha256 不符：%s（记录 %s，实际 %s）"
            % (npy_path, expected[:12], digest[:12]))
    info = dict(meta)
    info["sha256_actual"] = digest
    info["path"] = str(npy_path)
    return info


def export_lead_fc(target: str | Path) -> list[Path]:
    """把随包资产复制到 ``target`` 目录，返回落盘的两个文件路径。"""
    source = lead_fc_dir()
    target = Path(target)
    target.mkdir(parents=True, exist_ok=True)
    written = []
    for name in ("lead_fc.npy", "lead_fc.json"):
        destination = target / name
        shutil.copyfile(source / name, destination)
        written.append(destination)
    return written


def describe_lines() -> list[str]:
    """``smj lead-fc --show`` 的人类可读输出。"""
    info = lead_fc_info()
    in_plane = info.get("in_plane") or [0, 0]
    layers = int(info.get("layers", 0))
    per_layer = int(info.get("n_atoms_per_layer", 0))
    return [
        "电极力常数资产：%s" % info["path"],
        "  元素/晶格常数：%s / %.6f Å"
        % (info.get("element", "?"), float(info.get("a_lattice_A", 0.0))),
        "  面内胞：%d×%d（每层 %d 个原子）" % (in_plane[0], in_plane[1], per_layer),
        "  层数：%d（%d 个原子，decimation 用前两个主层）" % (layers, layers * per_layer),
        "  力常数截断：%.1f Å" % float(info.get("fc_cutoff_A", 0.0)),
        "  sha256：%s" % str(info.get("sha256", ""))[:16],
        "  来源：%s" % str(info.get("source", "")),
    ]
