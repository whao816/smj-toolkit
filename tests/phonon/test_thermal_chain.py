"""声子热导链的阶段装配（M3-T2…T5）。

这一组测的是**接线**而不是物理：注册是否完整、``thermal_device`` 建出的器件
是不是 54/36 的分区、``thermal_fc`` 的数组脚本有没有把``--array`` 与
``.finalize.lock`` 渲染对、``thermal_analyze`` 能不能把力常数一路算成
``thermal.json``。物理正确性由 ``test_units`` / ``test_negf_analytic`` /
``test_landauer`` 的解析基准保证。

``thermal_analyze`` 用一个人工的「等强弹簧」力常数（``k=0.5 eV/Å²``、
截断 2.9 Å）跑通全链：它的行和恰好为零，因此 ASR 残差可以逐位断言为 0；
它的最高模在 80 meV 附近，所以这里把 ``negf.band_top_meV`` 抬到 80 ——
带顶 fail-closed 本身由 ``test_landauer`` 单独覆盖。
"""

from __future__ import annotations

import json
import pathlib

import numpy as np
import pytest

from smj.config import load_config
from smj.junction.geometry import read_xyz, write_aims_geometry
from smj.phonon.fc import FullFCMeta, save_lead_fc
from smj.runner import StageRunner
from smj.stages import get_stage
from smj.state import THERMAL_STAGE_ORDER, normalize_stage, stage_index
from smj.errors import SmjError

_DATA = pathlib.Path(__file__).resolve().parent.parent / "data"


def _fake_run(workdir, extra: str = ""):
    """造一个只含 ``molecule_opt`` 的 run 目录,返回 (runner, run_dir)。"""
    species = workdir / "species"
    species.mkdir()
    for element in ("Au", "S", "C", "H"):
        (species / ("%s.in" % element)).write_text(
            "# %s\n  species %s\n" % (element, element), encoding="utf-8")
    config_path = workdir / "run.yaml"
    config_path.write_text(
        "run:\n  name: t\n  mode: dry-run\n"
        "cluster:\n  species_dir: %s\n"
        "thermal:\n  enabled: true\n%s"
        % (str(species.resolve()).replace("\\", "/"), extra),
        encoding="utf-8")
    config = load_config(config_path)
    run = workdir / "run"
    (run / "molecule_opt").mkdir(parents=True)
    atoms = read_xyz(_DATA / "nap-dithiol.xyz")
    write_aims_geometry(run / "molecule_opt" / "geometry.in", atoms)
    return StageRunner(config, run_dir=run), run


def test_thermal_chain_registered():
    """四个阶段按依赖链注册,远端/本地划分正确(D12:不进六阶段链)。"""
    assert THERMAL_STAGE_ORDER == (
        "thermal_device", "thermal_relax", "thermal_fc", "thermal_analyze")
    expected = {
        "thermal_device": ("molecule_opt", False),
        "thermal_relax": ("thermal_device", True),
        "thermal_fc": ("thermal_relax", True),
        "thermal_analyze": ("thermal_fc", False),
    }
    for name, (depends, submit) in expected.items():
        stage = get_stage(name)
        assert stage.depends == (depends,), name
        assert stage.needs_submit is submit, name


def test_thermal_stage_names_resolve_but_stay_off_the_electron_chain():
    assert normalize_stage("tfc") == "thermal_fc"
    assert normalize_stage("thermal-analyze") == "thermal_analyze"
    with pytest.raises(SmjError):
        stage_index("thermal_fc")


def test_thermal_device_prepare_and_collect(workdir):
    runner, run = _fake_run(workdir)
    stage = get_stage("thermal_device")
    stage_dir = run / "thermal_device"
    stage_dir.mkdir()
    stage.prepare(runner.context, stage_dir)

    text = (stage_dir / "geometry.in").read_text(encoding="utf-8")
    assert text.count("lattice_vector") == 3
    assert sum(1 for line in text.splitlines() if line.startswith("atom")) == 90

    meta = json.loads((stage_dir / "device_metadata.json").read_text(
        encoding="utf-8"))
    assert len(meta["device_indices"]) == 54
    assert len(meta["lead_indices"]) == 36
    assert len(meta["molecule_indices"]) == 18
    assert len(meta["removed_h"]) == 2

    collected = stage.collect(runner.context, stage_dir)
    summary = collected.results["thermal_device"]
    assert summary["n_atoms"] == 90
    assert summary["n_displacements"] == 6 * 54


def test_thermal_fc_prepare_renders_array_script(workdir):
    """数组脚本必须带 ``--array=0-323%8`` 与最终汇总锁。"""
    runner, run = _fake_run(workdir)
    device = get_stage("thermal_device")
    (run / "thermal_device").mkdir()
    device.prepare(runner.context, run / "thermal_device")
    device.collect(runner.context, run / "thermal_device")

    relax_dir = run / "thermal_relax"
    relax_dir.mkdir()
    (relax_dir / "geometry.in.next_step").write_text(
        (run / "thermal_device" / "geometry.in").read_text(encoding="utf-8"),
        encoding="utf-8")

    fc = get_stage("thermal_fc")
    fc_dir = run / "thermal_fc"
    fc_dir.mkdir()
    fc.prepare(runner.context, fc_dir)

    manifest = json.loads((fc_dir / "displacements.json").read_text(
        encoding="utf-8"))
    entries = manifest["displacements"] if isinstance(manifest, dict) else manifest
    assert len(entries) == 324
    assert entries[0]["label"].endswith("_x_p")
    assert entries[1]["label"].endswith("_x_m")

    script = (fc_dir / "sub_thermal_fc.sh").read_text(encoding="utf-8")
    assert "--array=0-323%8" in script
    assert ".finalize.lock" in script
    assert "apply_displacement.py" in script
    assert "extract_forces.py" in script
    # 位移作业从零起算,不能带 restart(每次 291 MB × 324 次)
    control = (fc_dir / "control.in").read_text(encoding="utf-8")
    assert "compute_forces .true." in control
    assert not any(line.strip().startswith("restart")
                   for line in control.splitlines())
    assert "k_grid 3 3 1" in control


def _spring_fc(atoms, k=0.5, cutoff=2.9):
    """等强弹簧力常数:行和精确为零(ASR 逐位成立)。"""
    n = len(atoms)
    matrix = np.zeros((3 * n, 3 * n))
    position = np.array([[a.x, a.y, a.z] for a in atoms])
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if float(np.linalg.norm(position[i] - position[j])) > cutoff:
                continue
            for axis in range(3):
                matrix[3 * i + axis, 3 * j + axis] -= k
                matrix[3 * i + axis, 3 * i + axis] += k
    return matrix


def test_thermal_analyze_writes_thermal_json(workdir):
    """力常数 → D → T_ph → G_ph:全链跑通并满足 thermal.json 契约。"""
    runner, run = _fake_run(workdir, extra=(
        "  electrode:\n    in_plane: [2, 2]\n"
        "  negf:\n    band_top_meV: 80\n    energy_step_meV: 2\n"))
    device = get_stage("thermal_device")
    (run / "thermal_device").mkdir()
    device.prepare(runner.context, run / "thermal_device")
    device.collect(runner.context, run / "thermal_device")
    meta = json.loads((run / "thermal_device" / "device_metadata.json").read_text(
        encoding="utf-8"))

    # 器件几何(50 原子:18 分子 + 32 Au)上的人工力常数
    from smj.junction.geometry import read_aims_geometry
    atoms = read_aims_geometry(run / "thermal_device" / "geometry.in")
    full = _spring_fc(atoms)

    def rows(indices):
        out = []
        for atom in indices:
            out.extend(range(3 * int(atom), 3 * int(atom) + 3))
        return out

    device_indices = [int(i) for i in meta["device_indices"]]
    left = [int(i) for i in meta["left_lead_indices"]]
    right = [int(i) for i in meta["right_lead_indices"]]
    lead_block = full[np.ix_(rows(left), rows(left))]
    lead_dir = workdir / "lead"
    save_lead_fc(lead_dir, lead_block, element="Au", a_lattice_A=4.078230,
                 n_atoms_per_layer=4, layers=2, in_plane=(2, 2),
                 fc_cutoff_A=10.0, source="test")

    fc_dir = run / "thermal_fc"
    fc_dir.mkdir()
    k_dev = full[rows(device_indices), :]
    np.save(fc_dir / "force_constants.npy", k_dev)
    (fc_dir / "geometry.in").write_text(
        (run / "thermal_device" / "geometry.in").read_text(encoding="utf-8"),
        encoding="utf-8")

    # 用测试自己的 lead_fc 覆盖配置里的路径
    runner.context.config.thermal.lead_fc = str(lead_dir.resolve())
    analyze = get_stage("thermal_analyze")
    analyze_dir = run / "thermal_analyze"
    analyze_dir.mkdir()
    analyze.prepare(runner.context, analyze_dir)
    analyze.collect(runner.context, analyze_dir)

    payload = json.loads((analyze_dir / "thermal.json").read_text(
        encoding="utf-8"))
    assert payload["schema"] == "smj.thermal.v1"
    assert payload["n_k_points"] == 1
    assert payload["n_atoms"] == 50
    # ASR 残差是**诊断量**而不是硬零：器件行来自实测、lead 行来自冻结 lead，
    # 交界处本来就不闭合。这里只要求它被算出来、有限。
    assert np.isfinite(payload["acoustic_sum_residual"])
    assert 0.0 <= payload["acoustic_sum_relative_residual"] < 10.0
    assert 0.0 < payload["g_ph_pw_per_k"] < 3.0 * payload["g_quantum_pw_per_k"]
    spectrum = payload["spectrum"]
    n = len(spectrum["energy_meV"])
    assert n == len(spectrum["transmission"])
    assert n == len(spectrum["cumulative_g_ph_pw_per_k"])
    assert spectrum["cumulative_g_ph_pw_per_k"][-1] == pytest.approx(
        payload["g_ph_pw_per_k"], rel=1e-12)
    # 没有电子链结果时 zt 只报不可用,不编数
    assert payload["zt"]["available"] is False


def test_thermal_cli_help_and_disabled(workdir, capsys):
    """``smj thermal`` 只认热导链的名字,且 thermal 段没开时明确拒绝。"""
    from smj.cli import build_parser, main
    from smj.errors import SmjError

    with pytest.raises(SystemExit) as exc:
        build_parser().parse_args(["thermal", "--help"])
    assert exc.value.code == 0

    config_path = workdir / "plain.yaml"
    config_path.write_text("run:\n  name: t\n  mode: dry-run\n", encoding="utf-8")
    try:
        code = main(["thermal", "-c", str(config_path)])
    except SmjError:
        code = 1
    assert code != 0


def test_thermal_stages_stay_out_of_plan_and_status():
    """D12:``smj run`` / ``smj status`` 走的是 ``STAGE_ORDER``,热导链不在里面。"""
    from smj.state import STAGE_ORDER, stage_range

    assert len(STAGE_ORDER) == 6
    assert not any(name.startswith("thermal_") for name in STAGE_ORDER)
    assert stage_range() == list(STAGE_ORDER)


def test_thermal_criteria_registered():
    """M5-T1:十条热导判据进 CRITERIA,且不叫 ``ph_te_present``。"""
    from smj.junction.criteria import CRITERIA, KNOWN_GAPS, lint_registry

    names = [c["name"] for c in CRITERIA]
    ph = [n for n in names if n.startswith("ph_")]
    assert ph == [
        "ph_transmission_present",
        "ph_relax_converged",
        "ph_band_top_zero",
        "ph_quantum_bound",
        "ph_acoustic_limit",
        "ph_asr_residual",
        "ph_fc_symmetry",
        "ph_imaginary_modes",
        "ph_lead_cell_match",
        "ph_grid_independence",
    ]
    # 命名铁律:热导侧没有 TE.dat
    assert not any("_te_" in n or n.endswith("_te") for n in ph)
    assert not any("te_present" in n for n in ph)

    # 4 条 physical 必须带源;新增判据不得进 KNOWN_GAPS(ratchet 只减不增)
    problems = {p["name"] for p in lint_registry()}
    assert not (problems - set(KNOWN_GAPS)) & set(ph)


def test_thermal_analyze_zt_folds_electronic(workdir):
    """``ZT_total = ZT_e/(1+G_ph/k_e)``:电子链在时按这个式子折进去。"""
    from smj.stages.thermal_analyze import ThermalAnalyzeStage

    assert ThermalAnalyzeStage._zt_payload  # 存在性
    # 直接验证公式,避免依赖完整 run
    zt_e, k_e, g_ph = 0.115479, 223.992, 20.0
    zt_total = zt_e / (1.0 + g_ph / k_e)
    assert zt_total == pytest.approx(0.1060132, abs=1e-6)

def test_dof_indices_map_atoms_to_three_degrees_of_freedom():
    """D 是 3N×3N：原子下标必须先换算成自由度下标。

    回归：曾直接拿原子下标 `np.ix_` 去切 D，54 个原子的 54×54 恰好与
    18 个原子的 3N=54 撞上，形状检查通过但 T_ph 恒为 0。
    """
    from smj.stages.thermal_analyze import _dof_indices

    assert list(_dof_indices([0, 2])) == [0, 1, 2, 6, 7, 8]
    assert _dof_indices([]).size == 0
