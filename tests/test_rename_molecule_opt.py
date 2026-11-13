"""阶段改名：``gaussian_opt`` -> ``molecule_opt``。

规范名是 ``molecule_opt``；``gaussian_opt`` 只作为旧 run / 旧命令的别名，
并且旧 run 目录（``gaussian_opt/``）必须继续可读。
"""

from __future__ import annotations

import json

from smj.report import discover_stages
from smj.state import (
    RunLayout,
    RunState,
    STAGE_ORDER,
    normalize_stage,
)
from smj.stages import get_stage


def test_canonical_stage_name_is_molecule_opt():
    assert STAGE_ORDER[0] == "molecule_opt"
    assert "gaussian_opt" not in STAGE_ORDER
    assert normalize_stage("gaussian_opt") == "molecule_opt"
    assert normalize_stage("gaussian") == "molecule_opt"
    assert normalize_stage("g16") == "molecule_opt"
    assert normalize_stage("opt") == "molecule_opt"
    assert normalize_stage("molecule-opt") == "molecule_opt"


def test_legacy_state_key_migrates_to_molecule_opt(workdir):
    run_dir = workdir / "run"
    run_dir.mkdir()
    (run_dir / "state.json").write_text(
        json.dumps(
            {
                "schema": "smj.state.v1",
                "stages": {
                    "gaussian_opt": {"status": "collected", "artifacts": ["opt.xyz"]},
                    "aims_tip_relax": {"status": "pending"},
                },
            }
        ),
        encoding="utf-8",
    )
    state = RunState.load(run_dir / "state.json")
    assert state.stage("gaussian_opt").status == "collected"
    assert state.stage("molecule_opt").artifacts == ["opt.xyz"]
    payload = state.to_dict()
    assert "molecule_opt" in payload["stages"]
    assert "gaussian_opt" not in payload["stages"]


def test_layout_falls_back_to_legacy_gaussian_opt_dir(workdir):
    run_dir = workdir / "run"
    (run_dir / "gaussian_opt").mkdir(parents=True)
    layout = RunLayout(run_dir)
    assert layout.stage_dir("molecule_opt") == (run_dir / "gaussian_opt").resolve()
    assert layout.stage_dir("gaussian_opt") == (run_dir / "gaussian_opt").resolve()


def test_new_run_uses_molecule_opt_dir(workdir):
    run_dir = workdir / "run"
    layout = RunLayout(run_dir)
    assert layout.stage_dir("molecule_opt") == (run_dir / "molecule_opt").resolve()
    legacy = layout.ensure_stage("molecule_opt")
    assert legacy == (run_dir / "molecule_opt").resolve()
    assert legacy.is_dir()


def test_report_discovers_legacy_dir_as_molecule_opt(workdir):
    run_dir = workdir / "run"
    (run_dir / "gaussian_opt").mkdir(parents=True)
    assert discover_stages(run_dir) == ["molecule_opt"]


def test_stage_registry_exposes_canonical_name():
    stage = get_stage("gaussian_opt")
    assert stage.name == "molecule_opt"
    tip = get_stage("aims_tip_relax")
    assert tip.depends == ("molecule_opt",)
    assert "molecule_opt" in tip.input_hint
