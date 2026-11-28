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
from types import SimpleNamespace

import numpy as np
import pytest

from smj.config import ThermalConfig
from smj.errors import StageError
from smj.junction.geometry import Atom
from smj.phonon import assets
from smj.phonon.device import DeviceSpec, build_device, write_device
from smj.phonon.fc import save_lead_fc
from smj.review import _thermal_lead_match
from smj.stages.thermal_analyze import _atom_distances, _build


def _layer_case(root, layers=6):
    molecule = [Atom("S", 0, 0, -3), Atom("S", 0, 0, 3)]
    atoms, meta = build_device(molecule, DeviceSpec(electrode_layers=layers,
                                                     strip_thiol_h=False))
    geometry, metadata = write_device(root / "thermal_device", atoms, meta)
    fc_path = root / "force_constants.npy"
    np.save(fc_path, np.zeros((3 * len(meta.device_indices), 3 * len(atoms))))
    return fc_path, metadata, geometry


def test_six_layers_registers_closest_frozen_pair_with_bulk_ab(workdir):
    fc, metadata, geometry = _layer_case(workdir)
    built = _build(fc, metadata, geometry, ThermalConfig())
    meta = json.loads(metadata.read_text())
    for side in ("left", "right"):
        principal = built[side + "_pl"]
        outer = built[side + "_outer"]
        assert {meta["layer_index"][i] for i in principal} == {3}
        assert {meta["layer_index"][i] for i in outer} == {2}
        # Even across the C→A boundary a common translation takes this
        # generated unwrapped slab pair onto bulk A→B in the same site order.
        positions = built["positions"]
        a1, a2 = np.asarray(meta["lattice"])[:2]
        reference = np.array([i*a1/3 + j*a2/3 for i in range(3) for j in range(3)])
        np.testing.assert_allclose(positions[outer, :2] - positions[outer[0], :2],
                                   reference[:, :2], atol=1.e-7)
        expected_shift = (a1 + a2)[:2] / 9
        np.testing.assert_allclose(positions[principal, :2] - positions[outer, :2],
                                   np.tile(expected_shift, (9, 1)), atol=1.e-7)


@pytest.mark.parametrize(("change", "value"), [
    ("in_plane", [1, 9]), ("a_lattice_A", 5.0), ("element", "Ag"),
    ("layers", 1), ("cutoff", 7.0), ("lattice_vector", 1.2),
])
def test_review_and_build_reject_incompatible_asset(workdir, change, value):
    fc, metadata, geometry = _layer_case(workdir)
    cfg = ThermalConfig()
    assert _thermal_lead_match(workdir, SimpleNamespace(thermal=cfg)) is True
    target = workdir / "lead"
    assets.export_lead_fc(target)
    cfg.lead_fc = str(target)
    asset_file = target / "lead_fc.json"
    asset = json.loads(asset_file.read_text())
    if change == "cutoff":
        cfg.fc.cutoff = value
    elif change == "lattice_vector":
        lines = geometry.read_text().splitlines()
        fields = lines[0].split()
        fields[1] = str(float(fields[1]) * value)
        lines[0] = " ".join(fields)
        geometry.write_text("\n".join(lines))
    else:
        asset[change] = value
        asset_file.write_text(json.dumps(asset))
    with pytest.raises((ValueError, StageError)):
        _build(fc, metadata, geometry, cfg)
    assert _thermal_lead_match(workdir, SimpleNamespace(thermal=cfg)) is False


def test_unspecified_cutoff_follows_custom_asset_but_explicit_value_is_checked(workdir):
    fc, metadata, geometry = _layer_case(workdir, layers=4)
    lead_dir = workdir / "lead"
    assets.export_lead_fc(lead_dir)
    path = lead_dir / "lead_fc.json"
    record = json.loads(path.read_text())
    record["fc_cutoff_A"] = 8.0
    path.write_text(json.dumps(record))
    thermal = ThermalConfig(lead_fc=str(path))
    assert _build(fc, metadata, geometry, thermal)["lead"].fc_cutoff_A == 8.0
    thermal.fc.cutoff = 10.0
    with pytest.raises(StageError, match="cutoff"):
        _build(fc, metadata, geometry, thermal)


def test_review_checks_matrix_identity_against_archived_result(workdir):
    _layer_case(workdir, layers=4)
    cfg = SimpleNamespace(thermal=ThermalConfig())
    output = workdir / "thermal_analyze" / "thermal.json"
    output.parent.mkdir()
    output.write_text(json.dumps({"lead_fc_sha256": "0" * 64}))
    assert _thermal_lead_match(workdir, cfg) is False
    assert _thermal_lead_match(workdir, None) is None
    output.write_text(json.dumps({"lead_fc_sha256": assets.lead_fc_info()["sha256"]}))
    assert _thermal_lead_match(workdir, cfg) is True


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

    reflection = np.diag([1.0, 1.0, -1.0])
    right_principal_rows = slice(3 * 4, 3 * 5)
    right_outer_columns = slice(3 * 5, 3 * 6)
    expected_right = reflection @ block[:3, 3:].T @ reflection
    actual_right = built["full"][right_principal_rows, right_outer_columns]
    assert np.allclose(actual_right, expected_right)


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


@pytest.mark.parametrize(("delta", "lattice"), [
    ([2.88374409, 2.49739564, 0.0],
     [[8.65123226, 0.0, 0.0], [4.32561613, 7.49218692, 0.0]]),
    ([1.47, 0.49 * np.sqrt(3.0), 0.0],
     [[2.0, 0.0, 0.0], [1.0, np.sqrt(3.0), 0.0]]),
])
def test_atom_distances_find_nearest_image_in_skew_cell(delta, lattice) -> None:
    a, b = np.asarray(lattice)
    positions = np.array([np.zeros(3), delta])
    distance = _atom_distances(positions, lattice)[0, 1]
    expected = min(np.linalg.norm(delta - i*a - j*b)
                   for i in (-1, 0, 1) for j in (-1, 0, 1))
    assert distance == pytest.approx(expected)
