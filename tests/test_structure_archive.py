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
    assert summary["atom_count"] == 2


GEOMETRY = """# junction geometry
atom 0.0 0.0 0.0 Au
atom 0.0 0.0 2.5 S
atom 0.0 0.0 4.0 Au
"""


def test_archive_aims_geometry_uses_native_format(workdir):
    run = workdir / "run"
    source = workdir / "geometry.in"
    source.write_text(GEOMETRY, encoding="utf-8")

    entry = archive_structure(
        run, source, stage="electrode_build", kind="junction", fmt="aims",
        results={"electrode_build": {"final_atom_count": 3,
                                     "surface_fingerprint": "abc123"}})

    assert entry["format"] == "aims"
    assert entry["file"].endswith(".in")
    assert entry["natoms"] == 3
    root = run / ARCHIVE_DIRNAME
    assert (root / entry["file"]).read_bytes() == source.read_bytes()
    assert entry["summary"] == {
        "final_atom_count": 3, "surface_fingerprint": "abc123",
        "results_key": "electrode_build",
    }


def test_summarize_results_prefers_stage_block():
    summary = summarize_results(
        {"electrode_build": {"final_atom_count": 7},
         "molecule_opt": {"atom_count": 2}},
        stage="electrode_build")
    assert summary["results_key"] == "electrode_build"
    assert summary["final_atom_count"] == 7


def test_safe_archive_failure_degrades_to_note(workdir, monkeypatch):
    from smj import structure_archive

    source = workdir / "a.xyz"
    source.write_text(XYZ, encoding="utf-8")

    def boom(*args, **kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(structure_archive, "archive_structure", boom)
    notes = structure_archive.safe_archive(
        workdir / "run", source, stage="molecule_opt", kind="input")

    assert notes and "结构归档失败" in notes[0]


def test_electrode_build_collect_archives_junction(workdir):
    from smj.stages.electrode_build import ElectrodeBuildStage

    run = workdir / "run"
    stage_dir = run / "electrode_build"
    stage_dir.mkdir(parents=True)
    (stage_dir / "geometry.in").write_text(GEOMETRY, encoding="utf-8")
    (stage_dir / "electrode_metadata.json").write_text(json.dumps({
        "final_atom_count": 3,
        "original_atom_count": 3,
        "surface_fingerprint": "deadbeef",
        "surface_indices": {"left_surface_indices": [1, 2, 3],
                            "right_surface_indices": [1, 2, 3]},
        "left_surface_indices": [1, 2, 3],
        "right_surface_indices": [1, 2, 3],
        "electrode_element": "Au",
        "electrode_a_lattice": 4.078,
    }), encoding="utf-8")
    ctx = SimpleNamespace(
        layout=SimpleNamespace(run_dir=run),
        config=SimpleNamespace(),
        require_stage_file=lambda name, filename: stage_dir / filename,
    )

    collected = ElectrodeBuildStage().collect(ctx, stage_dir)

    entry = _entries(workdir)[-1]
    assert entry["stage"] == "electrode_build"
    assert entry["kind"] == "junction"
    assert entry["format"] == "aims"
    assert entry["natoms"] == 3
    assert entry["summary"]["final_atom_count"] == 3
    assert entry["summary"]["surface_fingerprint"] == "deadbeef"
    assert (run / ARCHIVE_DIRNAME / entry["file"]).read_bytes() == (
        stage_dir / "geometry.in").read_bytes()
    assert any("结构已归档" in note for note in collected.notes)


# --------------------------------------------------------------------------- #
# 2026-11-20 复审 P1-5：manifest 损坏后不得从头编号、更不得覆盖既有结构
# --------------------------------------------------------------------------- #
def test_corrupt_manifest_does_not_overwrite_first_structure(workdir):
    run = workdir / "run"
    first_src = workdir / "first.xyz"
    first_src.write_text(XYZ, encoding="utf-8")
    second_src = workdir / "second.xyz"
    second_src.write_text("1\nsecond\nH 0 1 0\n", encoding="utf-8")

    first = archive_structure(run, first_src, stage="molecule_opt", kind="input")
    root = run / ARCHIVE_DIRNAME
    first_path = root / first["file"]
    before = first_path.read_bytes()

    (root / "manifest.json").write_text("not json\n", encoding="utf-8")
    second = archive_structure(run, second_src, stage="molecule_opt", kind="input")

    assert first_path.read_bytes() == before, "首份结构被覆盖了"
    assert second["file"] != first["file"]
    assert second["seq"] > first["seq"], "清单损坏后必须从现存快照续号"

    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    seqs = [entry["seq"] for entry in manifest["entries"]]
    assert first["seq"] in seqs, "恢复时应把现存快照重新登记"
    assert second["seq"] in seqs


def test_archive_never_overwrites_existing_structure_file(workdir):
    """清单缺条目但文件已占号（如崩溃残留）：换号，不覆盖。"""
    run = workdir / "run"
    src = workdir / "a.xyz"
    src.write_text(XYZ, encoding="utf-8")
    first = archive_structure(run, src, stage="molecule_opt", kind="input")
    root = run / ARCHIVE_DIRNAME

    occupied = root / "0002_molecule_opt_input.xyz"
    occupied.write_text("sentinel\n", encoding="utf-8")

    second = archive_structure(run, src, stage="molecule_opt", kind="input")
    assert occupied.read_text(encoding="utf-8") == "sentinel\n"
    assert second["seq"] == 3 and second["file"].startswith("0003_")
    assert first["file"].startswith("0001_")
