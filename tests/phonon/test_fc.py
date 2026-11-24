"""``fc.py``:位移集、中心差分、ASR 投影、lead 力常数装载与拼接。

中心差分的两条性质在这里被钉死:①对**线性**体系逐位复原 K;②对含三次非谐项
的体系误差随 δ 呈 O(δ²) 下降(δ 减半 → 误差降到 1/4)。第二条是"位移幅度选
0.01 Å"的真正依据 —— 不是拍脑袋,是误差阶数允许的最大步长。
"""
from __future__ import annotations

import json

import numpy as np
import pytest

from smj.phonon.dynamical import dynamical_matrix
from smj.phonon.fc import (
    FullFCMeta,
    acoustic_residual_of,
    apply_asr,
    asr_diagnostics,
    assemble_fc,
    assemble_full_fc,
    lead_device_conflicts,
    load_lead_fc,
    make_displacements,
    save_lead_fc,
)
from smj.phonon.units import EV_PER_SQRT

N = 3
MASSES = np.array([12.011, 12.011, 12.011])
SPRING = 10.0  # eV/Ang^2


def _chain_k() -> np.ndarray:
    """三原子一维链的解析力常数矩阵(只有 x-x 耦合)。"""
    k = SPRING
    block = np.array([[k, -k, 0.0], [-k, 2 * k, -k], [0.0, -k, k]])
    matrix = np.zeros((3 * N, 3 * N))
    for i in range(N):
        for j in range(N):
            matrix[3 * i, 3 * j] = block[i, j]
    return matrix


def _forces(matrix: np.ndarray, atom: int, axis: int, delta: float,
            anharmonic: float = 0.0) -> np.ndarray:
    """原子 ``atom`` 沿 ``axis`` 位移 ``delta`` 时全部原子受的力。"""
    displacement = np.zeros(3 * N)
    displacement[3 * atom + axis] = delta
    forces = -matrix @ displacement
    if anharmonic:
        forces[3 * atom + axis] -= anharmonic * delta ** 3
    return forces.reshape(N, 3)


def _displacements(amplitude: float = 0.01) -> list[dict]:
    return make_displacements(N, range(N), amplitude)


def test_make_displacements_order_and_count():
    records = _displacements(0.01)
    assert len(records) == 6 * N
    assert [r["slot"] for r in records] == list(range(6 * N))
    assert [r["label"] for r in records][:2] == ["0000_x_p", "0000_x_m"]
    assert [r["label"] for r in records][-2:] == ["0002_z_p", "0002_z_m"]
    assert records[0]["delta_A"] == pytest.approx(0.01)
    assert records[1]["delta_A"] == pytest.approx(-0.01)


def test_make_displacements_only_device_atoms():
    records = make_displacements(5, [1, 3], 0.01)
    assert len(records) == 12
    assert {r["atom"] for r in records} == {1, 3}


@pytest.mark.parametrize("bad", [0.0, -0.01, float("nan"), float("inf")])
def test_make_displacements_rejects_bad_amplitude(bad):
    with pytest.raises(ValueError):
        make_displacements(N, range(N), bad)


def test_make_displacements_rejects_out_of_range():
    with pytest.raises(ValueError):
        make_displacements(N, [0, N], 0.01)
    with pytest.raises(ValueError):
        make_displacements(N, [], 0.01)


def test_assemble_fc_recovers_analytic_matrix():
    matrix = _chain_k()
    disps = _displacements()
    forces = [_forces(matrix, d["atom"], d["axis"], d["delta_A"]) for d in disps]
    fc = assemble_fc(disps, forces, N)
    assert fc.shape == (3 * N, 3 * N)
    np.testing.assert_allclose(fc, matrix, atol=1e-12)


def test_central_difference_error_is_second_order():
    """含三次非谐项时误差 ∝ δ²:δ 减半,误差降到约 1/4。"""
    matrix = _chain_k()
    anharmonic = 500.0

    def error(amplitude: float) -> float:
        disps = _displacements(amplitude)
        forces = [_forces(matrix, d["atom"], d["axis"], d["delta_A"], anharmonic)
                  for d in disps]
        fc = assemble_fc(disps, forces, N)
        rows = np.arange(3 * N)
        return float(np.max(np.abs(fc[np.ix_(rows, rows)] - matrix)))

    coarse = error(0.02)
    fine = error(0.01)
    assert coarse > fine > 0.0
    assert 3.5 < coarse / fine < 4.5


def test_assemble_fc_rejects_unpaired_displacement():
    matrix = _chain_k()
    disps = _displacements()
    forces = [_forces(matrix, d["atom"], d["axis"], d["delta_A"]) for d in disps]
    with pytest.raises(ValueError):
        assemble_fc(disps[:-1], forces[:-1], N)


def test_assemble_fc_rejects_bad_force_shape():
    disps = _displacements()
    forces = [np.zeros((N + 1, 3)) for _ in disps]
    with pytest.raises(ValueError):
        assemble_fc(disps, forces, N)


def test_apply_asr_removes_translation_modes():
    matrix = _chain_k() + 0.05  # 给每一行加常数 → 明显违反 ASR
    dynamical = dynamical_matrix(matrix, MASSES)
    projected, residual = apply_asr(dynamical, MASSES)

    assert residual == pytest.approx(0.05 * 9, rel=1e-6)  # 每行 9 列 × 0.05
    eigenvalues = np.linalg.eigvalsh(projected)
    assert np.max(np.abs(eigenvalues[:3])) < 1e-10
    # 投影不能把光学模一起抹掉;D 的单位是 eV²(本征值即 (ħω)²),
    # C–C 弹簧 10 eV/Å² 的量级就在 1e-2 eV²,不是 1。
    assert np.max(eigenvalues) > 1.0e-3


def test_apply_asr_is_noop_for_a_clean_matrix():
    matrix = _chain_k()
    projected, residual = apply_asr(dynamical_matrix(matrix, MASSES), MASSES)
    assert residual < 1e-12
    np.testing.assert_allclose(
        np.linalg.eigvalsh(projected),
        np.linalg.eigvalsh(dynamical_matrix(matrix, MASSES)), atol=1e-10)


def test_apply_asr_validates_inputs():
    with pytest.raises(ValueError):
        apply_asr(np.zeros((4, 4)), MASSES)
    with pytest.raises(ValueError):
        apply_asr(dynamical_matrix(_chain_k(), MASSES), np.array([0.0, 1.0, 1.0]))


def _lead(workdir, layers: int = 2, per_layer: int = 2):
    n = per_layer * layers
    rng = np.random.default_rng(7)
    block = rng.normal(size=(3 * n, 3 * n))
    block = (block + block.T) / 2
    return save_lead_fc(workdir / "lead", block, element="Au",
                        a_lattice_A=4.07823, n_atoms_per_layer=per_layer,
                        layers=layers, in_plane=(3, 3), fc_cutoff_A=10.0,
                        source="test")


def test_lead_fc_roundtrip(workdir):
    lead = _lead(workdir)
    loaded = load_lead_fc(workdir / "lead")
    assert loaded.n_atoms == 4
    assert loaded.in_plane == (3, 3)
    np.testing.assert_allclose(loaded.force_constants, lead.force_constants)
    assert loaded.sha256 == lead.sha256
    payload = json.loads((workdir / "lead" / "lead_fc.json").read_text())
    assert payload["units"] == "eV/Ang^2"


def test_lead_fc_detects_tampering(workdir):
    lead = _lead(workdir)
    array = np.load(workdir / "lead" / "lead_fc.npy")
    np.save(workdir / "lead" / "lead_fc.npy", array + 1.0)
    with pytest.raises(ValueError) as excinfo:
        load_lead_fc(workdir / "lead")
    assert "sha256" in str(excinfo.value)


def test_lead_fc_missing_metadata(workdir):
    with pytest.raises(FileNotFoundError):
        load_lead_fc(workdir / "nope")


def test_assemble_full_fc_places_blocks(workdir):
    lead = _lead(workdir)
    # lead FC 是 2 层 × 2 原子 = 4 原子/侧,所以每侧要给 4 个下标
    meta = FullFCMeta(n_atoms=12, device_indices=(4, 5),
                      left_indices=(0, 1, 2, 3), right_indices=(8, 9, 10, 11))
    k_dev = np.arange(3 * 2 * 3 * 12, dtype=float).reshape(6, 36)
    full = assemble_full_fc(k_dev, lead, meta)
    assert full.shape == (36, 36)
    np.testing.assert_allclose(full[0:12, 0:12], lead.force_constants)
    np.testing.assert_allclose(full[24:36, 24:36], lead.force_constants)
    np.testing.assert_allclose(full[12:18, :], k_dev)
    # 未分配原子(6,7)的**行**是 0;器件行仍会测到它们的力系数,所以列不为零
    assert np.all(full[18:24, :] == 0.0)
    assert np.all(full[0:12, 18:24] == 0.0)
    assert np.all(full[24:36, 18:24] == 0.0)


def test_lead_device_conflicts_counts_differences(workdir):
    lead = _lead(workdir)
    meta = FullFCMeta(n_atoms=12, device_indices=(4, 5),
                      left_indices=(0, 1, 2, 3), right_indices=(8, 9, 10, 11))
    zeros = np.zeros((6, 36))
    assert lead_device_conflicts(zeros, lead, meta) == 0
    assert lead_device_conflicts(np.full((6, 36), 1.0), lead, meta) > 0


def test_full_fc_meta_rejects_overlap():
    with pytest.raises(ValueError):
        FullFCMeta(n_atoms=6, device_indices=(2,), left_indices=(2,),
                   right_indices=(4,))
    with pytest.raises(ValueError):
        FullFCMeta(n_atoms=6, device_indices=(2,), left_indices=(0, 0),
                   right_indices=(4,))
    with pytest.raises(ValueError):
        FullFCMeta(n_atoms=6, device_indices=(6,), left_indices=(0,),
                   right_indices=(4,))


def test_lead_cell_mismatch_is_rejected(workdir):
    lead = _lead(workdir)
    meta = FullFCMeta(n_atoms=6, device_indices=(2,),
                      left_indices=(0,), right_indices=(4, 5))
    with pytest.raises(ValueError):
        assemble_full_fc(np.zeros((3, 18)), lead, meta)


def test_rebuild_roundtrip_through_dynamical_matrix():
    matrix = _chain_k()
    dynamical = dynamical_matrix(matrix, MASSES)
    from smj.phonon.fc import _rebuild_force_constants
    np.testing.assert_allclose(_rebuild_force_constants(dynamical, MASSES),
                               matrix, atol=1e-12)
    assert EV_PER_SQRT == pytest.approx(0.064654151301)


def _with_row_shift(shift: float) -> np.ndarray:
    """在第 0 行对角元上抬 ``shift``,制造一个已知大小的 ASR 违反。

    只动对角元是为了保持矩阵对称 —— ``dynamical_matrix`` 默认会对
    ``D`` 做对称化,若只抬一整行会被折半,预期值就对不上了。
    """
    matrix = _chain_k().copy()
    matrix[0, 0] += shift
    return matrix


def test_asr_diagnostics_measures_before_projecting():
    """低于门限时只测不投影 —— 投影是"超阈值才做"的动作(D13)。"""
    shift = 0.01
    matrix = _with_row_shift(shift)
    dynamical = dynamical_matrix(matrix, MASSES)
    report = asr_diagnostics(dynamical, MASSES, threshold=1.0e-3)

    assert report["projected"] is False
    assert report["projection_norm"] is None
    assert report["residual_after_ev_per_A2"] is None
    assert report["residual_ev_per_A2"] == pytest.approx(shift, abs=1e-12)
    assert report["relative_residual"] == pytest.approx(
        shift / np.max(np.abs(matrix)), rel=1e-12)
    # 没投影时返回的矩阵必须原样
    np.testing.assert_array_equal(report["dynamical"], dynamical)


def test_asr_diagnostics_projects_when_above_threshold():
    matrix = _with_row_shift(5.0)      # 相对残差 ~5/10 = 0.5 >> 1e-3
    dynamical = dynamical_matrix(matrix, MASSES)
    report = asr_diagnostics(dynamical, MASSES, threshold=1.0e-3)

    assert report["projected"] is True
    assert report["relative_residual"] > 1.0e-3
    # 判据用的是投影前的值,所以它必须仍是那个大数
    assert report["residual_ev_per_A2"] == pytest.approx(5.0, abs=1e-12)
    # 投影后行和回到机器零,且投影本身是个小量
    assert report["residual_after_ev_per_A2"] < 1e-10
    assert 0.0 < report["projection_norm"] < 0.5
    assert not np.array_equal(report["dynamical"], dynamical)
    assert acoustic_residual_of(report["dynamical"], MASSES) < 1e-10


def test_asr_diagnostics_threshold_zero_always_projects():
    """门限取 0 时任何非零残差都要投影(便于确认投影确实能消掉残差)。"""
    matrix = _with_row_shift(1.0e-9)
    report = asr_diagnostics(dynamical_matrix(matrix, MASSES), MASSES,
                             threshold=0.0)
    assert report["projected"] is True
    assert report["residual_after_ev_per_A2"] < 1e-10

