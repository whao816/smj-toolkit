"""lead 力常数块与器件层号的对应关系。

``scripts/gen_au_lead_fc.py`` 生成的 ``lead_fc.npy`` 按「参考层 0 → 参考层 1」
排列:前 3N 个自由度是层 0,后 3N 个是层 1;``K(0,1)`` 存在左上-右下的块里。
器件里 ``layer 0`` 是**最外层**、``layer 1`` 是紧挨器件区的**主层**,所以装配时
必须按 ``outer + principal`` 的顺序传下标。反过来传(``principal + outer``)会把
``K(0,1)`` 当成 ``K(1,0)``,即丢掉一次转置——矩阵仍然对称(块本身是 ``[[K00,
K01], [K01^T, K00]]``),但方向错了,结果不会报错只会偏。
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from smj.config import ThermalConfig
from smj.phonon.fc import save_lead_fc
from smj.stages.thermal_analyze import _atom_distances, _build


def _asymmetric_lead_block() -> np.ndarray:
    """2 层 × 1 原子的 lead 块,层间耦合刻意非对称以暴露转置错误。"""
    coupling = np.array([[0.20, 0.07, 0.00],
                         [0.03, 0.15, 0.05],
                         [0.01, 0.02, 0.10]])
    block = np.zeros((6, 6))
    block[:3, :3] = 0.50 * np.eye(3)
    block[3:, 3:] = 0.50 * np.eye(3)
    block[:3, 3:] = coupling
    block[3:, :3] = coupling.T
    return block


def test_lead_block_is_ordered_outer_then_principal(workdir) -> None:
    root = Path(workdir).resolve() / "run"
    fc_dir = root / "thermal_fc"
    dev_dir = root / "thermal_device"
    lead_dir = root / "lead"
    for directory in (fc_dir, dev_dir, lead_dir):
        directory.mkdir(parents=True)

    # 6 原子:器件区 0/1;左侧 outer=2、principal=3;右侧 principal=4、outer=5。
    n_atoms = 6
    np.save(fc_dir / "force_constants.npy",
            np.zeros((3 * 2, 3 * n_atoms)))
    (fc_dir / "geometry.in").write_text(
        "".join("atom 0.0 0.0 %.3f  Au\n" % (1.5 * i) for i in range(n_atoms)),
        encoding="utf-8")
    block = _asymmetric_lead_block()
    save_lead_fc(lead_dir, block, element="Au", a_lattice_A=4.07823,
                 n_atoms_per_layer=1, layers=2, in_plane=(1, 1),
                 fc_cutoff_A=10.0)
    (dev_dir / "device_metadata.json").write_text(json.dumps({
        "n_atoms": n_atoms,
        "device_indices": [0, 1],
        "lead_indices": [2, 3, 4, 5],
        "layer_index": [2, 3, 0, 1, 1, 0],
        "side": ["left", "left", "left", "left", "right", "right"],
        "molecule_indices": [0, 1],
    }), encoding="utf-8")

    thermal = ThermalConfig()
    thermal.lead_fc = str(lead_dir)
    built = _build(fc_dir / "force_constants.npy",
                   dev_dir / "device_metadata.json",
                   fc_dir / "geometry.in", thermal)

    # 左侧:layer 0 = atom 2(outer)、layer 1 = atom 3(principal)。
    # 正确的 K(principal, outer) = K(0,1)^T。
    principal_rows = slice(3 * 3, 3 * 3 + 3)
    outer_columns = slice(3 * 2, 3 * 2 + 3)
    expected = block[:3, 3:].T
    actual = built["full"][principal_rows, outer_columns]
    assert np.allclose(actual, expected), (
        "lead 块顺序错了:期望 K(principal, outer)=K(0,1)^T,实际拿到 K(0,1)")


def test_atom_distances_wrap_in_plane_periodic_images() -> None:
    """lead 块来自周期体相:按器件坐标的裸距离会把跨边界邻居当成远邻。

    nap-sh 的 lead 对 (18,35) 裸距 11.89 Å、最小像 3.11 Å。截断敏感性如果
    用裸距,会在名义 10 Å 处把真实存在的近邻切掉,造出 4% 的伪差异。
    """
    positions = np.array([[0.0, 0.0, 0.0], [8.60, 0.0, 0.0]])
    lattice = ((8.65, 0.0, 0.0), (4.325, 7.49, 0.0), (0.0, 0.0, 40.0))
    raw = _atom_distances(positions)[0, 1]
    wrapped = _atom_distances(positions, lattice)[0, 1]
    assert raw > 8.0
    assert wrapped < 0.1
