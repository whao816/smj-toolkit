"""2026-11-13：第三段规范名 au60_build → electrode_build 的兼容契约。"""

from __future__ import annotations

from smj.state import (
    STAGE_ORDER,
    RunLayout,
    RunState,
    normalize_stage,
)


def test_canonical_name_and_aliases():
    assert STAGE_ORDER[2] == "electrode_build"
    assert normalize_stage("electrode_build") == "electrode_build"
    assert normalize_stage("au60_build") == "electrode_build"
    assert normalize_stage("au20_build") == "electrode_build"
    assert normalize_stage("au60-build") == "electrode_build"
    assert normalize_stage("build") == "electrode_build"
    assert normalize_stage("electrode") == "electrode_build"


def test_old_au60_directory_is_read_through_layout(workdir):
    run_dir = workdir / "run"
    legacy = (run_dir / "au60_build").resolve()
    legacy.mkdir(parents=True)
    (legacy / "geometry.in").write_text("legacy\n", encoding="utf-8")

    layout = RunLayout(run_dir)
    assert layout.stage_dir("electrode_build") == legacy
    assert layout.electrode_geometry().read_text(encoding="utf-8") == "legacy\n"


def test_old_auN_directory_is_read_through_layout(workdir):
    run_dir = workdir / "run"
    legacy = (run_dir / "au20_build").resolve()
    legacy.mkdir(parents=True)
    (legacy / "geometry.in").write_text("small electrode\n", encoding="utf-8")

    layout = RunLayout(run_dir)
    assert layout.stage_dir("au20_build") == legacy
    assert layout.stage_dir("electrode_build") == legacy


def test_old_state_key_migrates_to_canonical(workdir):
    run_dir = workdir / "run"
    run_dir.mkdir(parents=True)
    state = RunState.from_dict({
        "run": "legacy",
        "stages": {"au60_build": {"status": "collected"}},
    })
    assert state.status("electrode_build") == "collected"
    assert "electrode_build" in state.to_dict()["stages"]
    assert "au60_build" not in state.to_dict()["stages"]


def test_a2_stage_argument_accepts_legacy_names(workdir):
    from smj.a2 import _resolve_build_dir

    run_dir = workdir / "run"
    legacy = (run_dir / "au60_build").resolve()
    legacy.mkdir(parents=True)

    assert _resolve_build_dir(run_dir, "au60_build").resolve() == legacy
    assert _resolve_build_dir(run_dir, "au20_build").resolve() == legacy
