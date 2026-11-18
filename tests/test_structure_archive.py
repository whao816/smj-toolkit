"""结构归档：输入/优化后的 xyz 与结果快照、manifest 的联动。"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

from smj.config import ElectrodeConfig, MoleculeConfig
from smj.engines import EngineCollected, EnginePrepared
from smj.stages import molecule_opt as molecule_opt_module
from smj.stages.molecule_opt import MoleculeOptStage
from smj.structure_archive import (
    ARCHIVE_DIRNAME,
    archive_structure,
    summarize_results,
)

XYZ = """2
H2 input
H 0.0 0.0 0.0
H 0.74 0.0 0.0
"""

OPT = """2
H2 optimized
H 0.0 0.0 0.0
H 0.75 0.0 0.0
"""


class _FakeEngine:
    name = "aims"
    script = "sub_aims.sh"
    outputs = ("opt.xyz",)
    resource_anchor = "aims_tip_relax"
    env_anchor = "aims_tip_relax"

    def prepare(self, ctx, stage_dir: Path, *, stage, molecule, xyz, natoms,
                electron_note, estimate, env) -> EnginePrepared:
        Path(stage_dir, "geometry.in").write_text("fake\n", encoding="utf-8")
        return EnginePrepared(
            inputs=[], artifacts=["geometry.in"], notes=["fake prepare"],
            script=self.script)

    def collect(self, ctx, stage_dir: Path, *, stage) -> EngineCollected:
        Path(stage_dir, "opt.xyz").write_text(OPT, encoding="utf-8")
        return EngineCollected(
            artifacts=["opt.xyz"],
            results={"molecule_opt": {
                "engine": "aims", "converged": True,
                "max_force_ev_per_A": 0.004, "relax_steps": 3,
            }},
            notes=["fake collect"])


def _ctx(root: Path) -> SimpleNamespace:
    run_dir = root / "run"
    run_dir.mkdir(parents=True, exist_ok=True)
    xyz = root / "molecule.xyz"
    xyz.write_text(XYZ, encoding="utf-8")
    molecule = MoleculeConfig(xyz=str(xyz), charge=0, multiplicity=1)
    config = SimpleNamespace(
        molecule=molecule, electrode=ElectrodeConfig(vdw=False))
    return SimpleNamespace(
        config=config,
        layout=SimpleNamespace(run_dir=run_dir),
        molecule_xyz=lambda: xyz,
        estimate=lambda *a, **k: {
            "walltime": "00:10:00", "time_limit": "00:20:00",
            "mem": "2G", "ntasks": 1, "cores": 1, "capped": [],
        },
        stage_env=lambda *a, **k: {"partition": "main", "ntasks": 1},
        log=lambda *_: None,
    )


def _entries(root: Path) -> list[dict]:
    path = root / "run" / ARCHIVE_DIRNAME / "manifest.json"
    return json.loads(path.read_text(encoding="utf-8"))["entries"]


def test_prepare_archives_input_structure(workdir, monkeypatch):
    ctx = _ctx(workdir)
    stage_dir = workdir / "run" / "molecule_opt"
    stage_dir.mkdir(parents=True, exist_ok=True)
    monkeypatch.setattr(
        molecule_opt_module, "engine_for_config", lambda cfg: _FakeEngine())
    prepared = MoleculeOptStage().prepare(ctx, stage_dir)

    entry = _entries(workdir)[0]
    assert entry["stage"] == "molecule_opt"
    assert entry["kind"] == "input"
    assert entry["engine"] == "aims"
    assert entry["natoms"] == 2
    assert entry["results_file"] is None
    archived = workdir / "run" / ARCHIVE_DIRNAME / entry["file"]
    source = workdir / "molecule.xyz"
    assert archived.read_bytes() == source.read_bytes()
    assert entry["sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert any("结构已归档" in note for note in prepared.notes)


def test_collect_archives_optimized_with_results_snapshot(workdir, monkeypatch):
    ctx = _ctx(workdir)
    stage_dir = workdir / "run" / "molecule_opt"
    stage_dir.mkdir(parents=True, exist_ok=True)
    monkeypatch.setattr(
        molecule_opt_module, "engine_for_config", lambda cfg: _FakeEngine())

    collected = MoleculeOptStage().collect(ctx, stage_dir)

    entry = _entries(workdir)[-1]
    assert entry["kind"] == "optimized"
    assert entry["results_file"]
    assert entry["summary"] == {
        "engine": "aims", "converged": True,
        "max_force_ev_per_A": 0.004, "relax_steps": 3,
        "results_key": "molecule_opt",
    }
    snapshot_path = workdir / "run" / ARCHIVE_DIRNAME / entry["results_file"]
    snapshot = json.loads(snapshot_path.read_text(encoding="utf-8"))
    assert snapshot["molecule_opt"]["converged"] is True
    assert entry["results_sha256"] == hashlib.sha256(
        snapshot_path.read_bytes()).hexdigest()
    assert any("结果快照" in note for note in collected.notes)


def test_archive_failure_only_adds_note(workdir, monkeypatch):
    ctx = _ctx(workdir)
    stage_dir = workdir / "run" / "molecule_opt"
    stage_dir.mkdir(parents=True, exist_ok=True)
    monkeypatch.setattr(
        molecule_opt_module, "engine_for_config", lambda cfg: _FakeEngine())

    def boom(*args, **kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(molecule_opt_module, "archive_structure", boom)
    collected = MoleculeOptStage().collect(ctx, stage_dir)

    # 归档是审计动作：失败不能吞掉已经算出来的科学结果。
    assert collected.results["molecule_opt"]["converged"] is True
    assert any("结构归档失败" in note for note in collected.notes)


def test_archive_appends_entries_and_keeps_bytes(workdir):
    run = workdir / "run"
    first = workdir / "a.xyz"
    second = workdir / "b.xyz"
    first.write_text(XYZ, encoding="utf-8")
    second.write_text(OPT, encoding="utf-8")

    e1 = archive_structure(run, first, stage="molecule_opt", kind="input")
    e2 = archive_structure(
        run, second, stage="molecule_opt", kind="optimized",
        results={"molecule_opt": {"converged": True}})

    assert (e1["seq"], e2["seq"]) == (1, 2)
    root = run / ARCHIVE_DIRNAME
    assert (root / e1["file"]).read_bytes() == first.read_bytes()
    assert e1["sha256"] != e2["sha256"]
    assert e2["results_file"] and (root / e2["results_file"]).is_file()
    assert len(_entries(workdir)) == 2


def test_summarize_results_keeps_scalars_only():
    summary = summarize_results({"gaussian": {
        "normal_termination": True,
        "final_scf_energy_hartree": -1.5,
        "imaginary_frequency_count": 0,
        "atom_count": 2,
        "geometry_check": {"ok": True},
    }})

    assert summary["results_key"] == "gaussian"
    assert summary["normal_termination"] is True
    assert summary["final_scf_energy_hartree"] == -1.5
    # 非标量块留给结果快照，manifest 只放可检索的摘要。
    assert "geometry_check" not in summary
    assert "atom_count" not in summary
