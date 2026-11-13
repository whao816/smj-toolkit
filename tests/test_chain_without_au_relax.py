"""2026-11-13：au_relax 从主链路移除后的状态机/报告契约。"""

from __future__ import annotations

import pytest

from smj.errors import SmjError
from smj.report import discover_stages
from smj.state import (
    LEGACY_STAGES,
    STAGE_ORDER,
    RunState,
    next_stage,
    normalize_stage,
    previous_stage,
    stage_index,
)


def test_stage_order_is_six_stages_without_au_relax():
    assert STAGE_ORDER == (
        "gaussian_opt",
        "aims_tip_relax",
        "au60_build",
        "aims_scf",
        "aitranss_te",
        "junction_analyze",
    )
    assert "au_relax" not in STAGE_ORDER
    assert "au_relax" in LEGACY_STAGES


def test_legacy_au_relax_keeps_aliases_and_neighbours():
    assert normalize_stage("au-relax") == "au_relax"
    assert previous_stage("au_relax") == "au60_build"
    assert next_stage("au_relax") == "aims_scf"
    with pytest.raises(SmjError, match="已从主链路移除"):
        stage_index("au_relax")


def test_run_state_creates_legacy_stage_only_on_demand():
    fresh = RunState(run="fresh")
    assert "au_relax" not in fresh.to_dict()["stages"]

    legacy = RunState(run="legacy")
    stage = legacy.stage("au_relax")
    stage.status = "collected"
    assert legacy.to_dict()["stages"]["au_relax"]["status"] == "collected"


def test_discover_stages_ignores_legacy_au_relax_dir(workdir):
    run_dir = workdir / "run"
    for name in STAGE_ORDER + ("au_relax",):
        (run_dir / name).mkdir(parents=True)
    assert discover_stages(run_dir) == list(STAGE_ORDER)
