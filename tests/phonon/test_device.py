"""``device.build_device``:周期板器件几何与角色划分(T5b)。

夹具 ``tests/data/nap-dithiol.xyz`` 是真实 run 里 ``molecule_opt/geometry.in``
的逐原子坐标(20 原子,中性 HS–R–SH)。``thermal_device`` 阶段接的正是这个产物,
所以"去 S–H"这一步必须在这里被测到 —— 直接用 ``aims_scf`` 那份(已去质子)
会把最容易出错的一步跳过去。
"""
from __future__ import annotations

import json
import math
import pathlib

import numpy as np
import pytest

from smj.junction.geometry import read_xyz
from smj.phonon.device import (
    DeviceSpec,
    build_device,
    fcc111_in_plane_vectors,
    write_device,
)

DATA = pathlib.Path(__file__).resolve().parent.parent / "data"
MOLECULE = read_xyz(DATA / "nap-dithiol.xyz")


def _dist(a, b) -> float:
    return math.dist(a.coord(), b.coord())


@pytest.fixture(scope="module")
def device4():
    return build_device(MOLECULE, DeviceSpec(electrode_layers=4))


def test_molecule_is_deprotonated(device4):
    atoms, meta = device4
    molecule = [atoms[i] for i in meta.molecule_indices]
    assert len(molecule) == 18, "中性二硫醇 20 原子,去两个 S–H 后应为 18"
    sulfur = [a for a in molecule if a.element == "S"]
    assert len(sulfur) == 2
    # 锚定 S 上不能还挂着 H(报告 §6 陷阱⑥)
    for s in sulfur:
        neighbours = [a for a in molecule
                      if a.element == "H" and _dist(a, s) < 1.70]
        assert neighbours == []


def test_total_atom_count(device4):
    atoms, meta = device4
    assert meta.n_atoms == 18 + 2 * 4 * 9
    assert meta.n_atoms == 90
    assert len(atoms) == 90


def test_each_sulfur_bonds_exactly_one_surface_gold(device4):
    atoms, meta = device4
    for anchor in (meta.anchors["left"], meta.anchors["right"]):
        sulfur = atoms[anchor]
        gold = [a for a in atoms
                if a.element == "Au" and _dist(a, sulfur) < 2.75]
        assert len(gold) == 1, "每个 S 只应有一个最近的表面 Au"
        assert 2.2 < _dist(gold[0], sulfur) < 2.6


def test_layer_index_is_monotonic_and_complete(device4):
    atoms, meta = device4
    left = [i for i in range(meta.n_atoms) if meta.side[i] == "left"]
    right = [i for i in range(meta.n_atoms) if meta.side[i] == "right"]
    assert len(left) == len(right) == 4 * 9
    assert all(atoms[i].z < 0 for i in left)
    assert all(atoms[i].z > 0 for i in right)
    # 层号 0 = 最外层 → 层号越大离分子(以及 z=0)越近
    for indices in (left, right):
        by_layer: dict[int, list[float]] = {}
        for index in indices:
            by_layer.setdefault(meta.layer_index[index], []).append(
                abs(atoms[index].z))
        assert sorted(by_layer) == [0, 1, 2, 3]
        means = [sum(by_layer[k]) / len(by_layer[k]) for k in range(4)]
        assert means == sorted(means, reverse=True)
        assert all(len(by_layer[k]) == 9 for k in range(4))
    for index in meta.molecule_indices:
        assert meta.layer_index[index] == -1
        assert meta.side[index] == "none"


def test_device_and_lead_partition(device4):
    atoms, meta = device4
    assert len(meta.molecule_indices) == 18
    assert len(meta.device_indices) == 18 + 2 * 2 * 9 == 54
    assert len(meta.lead_indices) == 2 * 2 * 9 == 36
    assert set(meta.device_indices) | set(meta.lead_indices) == set(
        range(meta.n_atoms))
    assert not (set(meta.device_indices) & set(meta.lead_indices))
    assert meta.roles[0] == "molecule"
    assert sorted(meta.device_indices) == meta.device_indices


def test_constrained_atoms_are_the_outer_layers(device4):
    atoms, meta = device4
    assert len(meta.constrained_indices) == 2 * 9
    for index in meta.constrained_indices:
        assert meta.layer_index[index] == 0
        assert meta.roles[index] == "lead_au"
    assert not (set(meta.constrained_indices) & set(meta.device_indices))


def test_surface_indices_are_closest_to_molecule(device4):
    atoms, meta = device4
    assert len(meta.left_surface_indices) == 9
    assert len(meta.right_surface_indices) == 9
    for index in meta.left_surface_indices:
        assert meta.layer_index[index] == 3
        assert meta.side[index] == "left"
        assert atoms[index].z < 0
    for index in meta.right_surface_indices:
        assert atoms[index].z > 0
    assert len(meta.sbond) == 18
    assert min(meta.sbond) == pytest.approx(2.35, abs=0.15)


def test_lattice_and_vacuum(device4):
    atoms, meta = device4
    a1, a2 = fcc111_in_plane_vectors(4.078230, (3, 3))
    np.testing.assert_allclose(meta.lattice[0], a1)
    np.testing.assert_allclose(meta.lattice[1], a2)
    assert meta.lattice[2][0] == 0.0 and meta.lattice[2][1] == 0.0
    z = [a.z for a in atoms]
    c_length = meta.lattice[2][2]
    assert c_length - (max(z) - min(z)) == pytest.approx(20.0, abs=1e-9)
    # 面内最近邻 = a/sqrt(2)
    assert math.dist((a1[0], a1[1], 0), (0, 0, 0)) == pytest.approx(3 * 2.8845,
                                                                  abs=0.01)


def test_adjacent_layers_keep_abc_stacking(device4):
    """相邻层 Au-Au 必须是 a/sqrt(2),不能是层间距 delta_z。

    2026-11-23 的 bug:平移量按**每一层**各自 snap 到 S 正下方,于是所有层
    对齐到同一个原子(AA 堆垛),相邻层 Au-Au 只有 2.355 A(=a/sqrt(3))。
    relax 一开始就报 5 eV/A 的力。这条测试直接钉住堆垛。
    """
    atoms, meta = device4
    coord = np.array([a.coord() for a in atoms])
    gold = [i for i, a in enumerate(atoms) if a.element == "Au"]
    nn = 4.078230 / math.sqrt(2)
    worst = math.inf
    for k, i in enumerate(gold):
        for j in gold[k + 1:]:
            if meta.side[i] != meta.side[j]:
                continue
            worst = min(worst, float(np.linalg.norm(coord[i] - coord[j])))
    assert worst == pytest.approx(nn, abs=1e-3), "相邻层 Au-Au %.4f A" % worst
    # ABC 堆垛:任一原子都不会与上下层原子在同一条 z 轴上
    for k, i in enumerate(gold):
        for j in gold[k + 1:]:
            if meta.side[i] != meta.side[j]:
                continue
            if abs(meta.layer_index[i] - meta.layer_index[j]) != 1:
                continue
            planar = float(np.linalg.norm(coord[i][:2] - coord[j][:2]))
            assert planar > 1.0, "层间原子在同一 xy 柱上(AA 堆垛)"


def test_abc_offset_uses_primitive_vectors_in_rectangular_supercell():
    atoms, meta = build_device(MOLECULE, DeviceSpec(in_plane=(2, 4)))
    left = meta.left_lead_indices
    surface_shift = np.array(atoms[left[8]].coord()) - atoms[left[0]].coord()
    a, b = np.asarray(meta.lattice[:2])
    np.testing.assert_allclose(surface_shift[:2], ((a / 2 + b / 4) / 3)[:2])


def test_six_layers_gives_126_atoms():
    atoms, meta = build_device(MOLECULE, DeviceSpec(electrode_layers=6))
    assert meta.n_atoms == 18 + 2 * 6 * 9 == 126
    assert len(meta.device_indices) == 18 + 2 * 2 * 9


def test_spec_validation():
    with pytest.raises(ValueError):
        DeviceSpec(electrode_layers=1)
    with pytest.raises(ValueError):
        DeviceSpec(electrode_layers=4, device_layers=4)
    with pytest.raises(ValueError):
        DeviceSpec(electrode_layers=4, device_layers=3)
    with pytest.raises(ValueError):
        DeviceSpec(vacuum_A=0.0)
    with pytest.raises(ValueError):
        DeviceSpec(in_plane=(0, 3))


def test_write_device(workdir, device4):
    atoms, meta = device4
    geometry, metadata = write_device(workdir / "thermal", atoms, meta)
    text = geometry.read_text(encoding="utf-8")
    assert text.count("lattice_vector") == 3
    assert text.count("constrain_relaxation") == len(meta.constrained_indices)
    payload = json.loads(metadata.read_text(encoding="utf-8"))
    assert payload["schema"] == "smj.thermal_device.v1"
    assert payload["n_atoms"] == 90
    assert payload["spec"]["electrode_layers"] == 4
    assert len(payload["device_indices"]) == 54


def test_rejects_molecule_without_two_anchors():
    """只有一个 S 的分子必须报错,不能猜一个"差不多的"锚定原子。"""
    with pytest.raises(ValueError):
        build_device([a for a in MOLECULE if a.element != "S"] + [
            type(MOLECULE[0])("S", 0.0, 0.0, -4.0)], DeviceSpec())
