"""2026-11-20 复审验收套件入库（对应附件 ``test_smj_review_20261120.py`` 8 条）。

附件由复审者在 0.2.2rc1 / 04739baa 上运行，8 条全部失败。这里按原样移植
（夹具改仓库 ``workdir``；``tmp_path`` 在 DSH 沙箱下会因 0700 目录不可访问而
连带失败），作为这一轮 7 项问题的固定回归：

| # | 断言 | 修法批次 |
| --- | --- | --- |
| 1 | 窄费米窗的常数谱不得判 passed | 批 1 / R1 |
| 2 | 乱序双通道 ``g_over_g0_total`` 必须与主积分一致 | 批 1 / R3 |
| 3 | 峰检测与行序无关 | 批 1 / R2 |
| 4 | 零字节 ``beta.aims`` 必须让 ``aims_scf.collect`` 拒绝 | 批 2 / R4 |
| 5 | ``finished`` 不等于已验收，manifest 不能算完成 | 批 3 / R6 |
| 6 | 物理失败不改变工程判定 | 批 3 / R7 |
| 7 | 警告数在报告生成后不丢 | 批 3 / R7 |
| 8 | 清单损坏后归档不得覆盖既有结构 | 批 2 / R5 |
"""
from __future__ import annotations

import json
import math
import random
from pathlib import Path
from types import SimpleNamespace

import pytest

from smj.errors import StageError
from smj.junction.criteria import evaluate_run_criteria, physics_verdict
from smj.junction.landauer import analyze_te
from smj.junction.validate import validate_junction_run
from smj.landauer_spin import analyze_te_spin
from smj.report import build_manifest
from smj.stages.aims_scf import AimsScfStage
from smj.state import STAGE_ORDER, RunState
from smj.structure_archive import archive_structure

# 与附件一致：E_F = −0.188452434085 Ha，1 Ha = 27.2114 eV。
EF_HARTREE = -0.188452434085
H_EV = 27.2114


def _write_te(folder: Path, rows, *, header: str = "# efermi = -0.18") -> Path:
    """附件口径：第 1 列 Ha、第 2 列 E−E_F (eV)，其后每列一个自旋通道。"""
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / "TE.dat"
    lines = [header]
    for energy, *channels in rows:
        fields = [EF_HARTREE + energy / H_EV, energy, *channels]
        lines.append(" ".join(str(value) for value in fields))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _collected_run(workdir: Path) -> Path:
    run_dir = workdir / "run"
    run_dir.mkdir(exist_ok=True)
    stages = {name: {"status": "collected"} for name in STAGE_ORDER}
    (run_dir / "state.json").write_text(
        json.dumps({"schema": "smj.state.v1", "stages": stages}),
        encoding="utf-8")
    return run_dir


def test_constant_spectrum_with_truncated_thermal_window(workdir: Path):
    """#1：跨零但只覆盖 ~23% 费米导数权重，截窗积分 G=0.0586 不得判 passed。"""
    rows = [(-0.025 + 0.026 * i / 199, 0.25) for i in range(200)]
    _write_te(workdir, rows)
    checks = evaluate_run_criteria(workdir, temperature=300, system="BDT")
    assert physics_verdict(checks)["physics_verdict"] != "passed"


def test_spin_sum_matches_total_for_unsorted_spectrum(workdir: Path):
    """#2：α/β 分通道积分必须按能量排序，不能按文件行序。"""
    rows = []
    for i in range(201):
        energy = -0.2 + i * 0.002
        delta = 0.05 * math.tanh(energy / 0.04)
        rows.append((energy, 0.25 + delta, 0.25 - delta))
    random.Random(7).shuffle(rows)
    te = _write_te(workdir, rows)
    spin = analyze_te_spin(te)
    assert spin["g_over_g0_total"] == pytest.approx(analyze_te(te)["L0"])


def test_peak_check_is_invariant_to_row_order(workdir: Path):
    """#3：同一单调谱，升序与乱序必须得到同一个物理判定。"""
    rows = []
    for i in range(201):
        energy = -0.5 + i * 0.005
        rows.append((energy, 0.25 - 0.15 * math.tanh(energy / 0.12)))
    _write_te(workdir, rows)
    first = evaluate_run_criteria(workdir, system="BDT")
    random.Random(7).shuffle(rows)
    _write_te(workdir, rows)
    second = evaluate_run_criteria(workdir, system="BDT")
    assert (physics_verdict(first)["physics_verdict"]
            == physics_verdict(second)["physics_verdict"])


def test_scf_rejects_empty_beta(workdir: Path):
    """#4：零字节 beta.aims 在开壳层 SCF 里必须报错，不能当"存在"。"""
    stage = workdir / "aims_scf"
    stage.mkdir()
    files = {
        "aims.out": "Self-consistency cycle converged.\nHave a nice day.\n",
        "basis-indices.out": "$basis.indices.aims\n# amount of basis functions used: 4\n",
        "omat.aims": "$omat.aims format(4d20.14) nsaos=4\n",
        "alpha.aims": "$scfmo.aims nsaos=4\n",
        "beta.aims": "",
        "geometry.in": "atom 0 0 0 S\natom 3 0 0 S\n",
    }
    for name, contents in files.items():
        (stage / name).write_text(contents, encoding="utf-8")
    ctx = SimpleNamespace(
        layout=SimpleNamespace(run_dir=workdir),
        config=SimpleNamespace(),
        require_stage_file=lambda stage_name, name: workdir / stage_name / name,
    )
    with pytest.raises(StageError, match=r"beta\.aims"):
        AimsScfStage().collect(ctx, stage)


def test_finished_without_collect_is_incomplete(workdir: Path):
    """#5：aims_scf 只到 finished（产物取回、未验收），run 不算完成。"""
    run_dir = _collected_run(workdir)
    state = RunState.load(run_dir / "state.json")
    state.set("aims_scf", status="finished")
    _write_te(run_dir / "aitranss_te",
              [(-0.4 + i * 0.002, 0.25 - i * 0.0001) for i in range(401)])
    manifest = build_manifest(run_dir, run_id="t", state=state,
                              system="BDT", include_checks=True)
    assert not manifest["verdict"]["completeness"]["is_done"]


def test_physical_failure_preserves_engineering_verdict(workdir: Path):
    """#6：物理判据失败时工程判定仍为通过，manifest 与 validate 一致。"""
    run_dir = _collected_run(workdir)
    _write_te(run_dir / "aitranss_te",
              [(-0.4 + i * 0.002, 0.25 + i * 0.0001) for i in range(401)])
    raw = validate_junction_run(run_dir, system="BDT")
    assert raw["physics_verdict"] == "failed" and raw["engineering_passed"]
    manifest = build_manifest(run_dir, run_id="t", include_checks=True)
    assert manifest["verdict"]["engineering_passed"] == raw["engineering_passed"]


def test_warning_count_survives_report_generation(workdir: Path):
    """#7：validate 的 warn 级警告必须原样进 manifest，不能被丢掉。"""
    run_dir = _collected_run(workdir)
    _write_te(run_dir / "aitranss_te",
              [(-0.4 + i * 0.002, 0.25 - i * 0.0001) for i in range(401)],
              header="# missing Fermi energy")
    raw = validate_junction_run(run_dir, system="BDT")
    assert raw["warnings"] > 0
    manifest = build_manifest(run_dir, run_id="t", include_checks=True)
    assert manifest["verdict"]["warnings"] == raw["warnings"]


def test_archive_recovery_preserves_old_structure(workdir: Path):
    """#8：清单写坏后重新归档，既有快照字节不得被覆盖。"""
    source = workdir / "input.xyz"
    source.write_bytes(b"1\nfirst\nH 0 0 0\n")
    entry = archive_structure(workdir, source, stage="molecule_opt", kind="input")
    archive = workdir / "saved_structures"
    original = (archive / entry["file"]).read_bytes()
    (archive / "manifest.json").write_text("{", encoding="utf-8")
    source.write_bytes(b"1\nsecond\nH 0 1 0\n")
    archive_structure(workdir, source, stage="molecule_opt", kind="input")
    assert (archive / entry["file"]).read_bytes() == original
