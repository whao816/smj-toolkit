"""Regression tests for the standalone Au lead generator and assembly CLI."""
import itertools
import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from scripts import gen_au_lead_fc as lead
from smj.errors import StageError


def _json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def _case(path, layer, in_plane=(1, 1)):
    atoms, lattice = lead.build_bulk_cell(6, in_plane)
    lead.write_aims_geometry(path / "geometry.in", atoms, lattice=lattice)
    (path / "control.in").write_text("xc pbe\n", encoding="utf-8")
    ids = lead.layer_indices(in_plane, 6, layer)
    disps = lead.make_displacements(len(atoms), ids, 0.01)
    meta = dict(layers=6, n_atoms=len(atoms), in_plane=list(in_plane),
                n_displacements=len(disps), a_lattice_A=lead.A_LATTICE_AU,
                surface_layer_indices=ids)
    (path / "lead_case.json").write_text(json.dumps(meta), encoding="utf-8")
    (path / "displacements.json").write_text(json.dumps(disps), encoding="utf-8")
    rows = np.ones((3 * len(ids), 3 * len(atoms)))
    if len(ids) == 1:
        rows[:] = 0.0
        rows[:, :3] = np.eye(3) * (1 if layer == 0 else 7)
        rows[:, 3:6] = np.eye(3) * (3 if layer == 0 else 11)
        # Individually large x/y residuals cancel in a scalar row sum.
        rows[0, 6:8] = [20, -20]
    for rec in disps:
        row = 3 * ids.index(rec["atom"]) + rec["axis"]
        force = -rows[row].reshape(-1, 3) * rec["delta_A"]
        slot = path / "fc" / ("%04d" % rec["slot"])
        slot.mkdir(parents=True)
        (slot / "forces.json").write_text(json.dumps({"forces": force.tolist()}),
                                          encoding="utf-8")
    return rows


def _assemble(out, *cases, cutoff=10):
    return lead.main(["--out", str(out), "--cutoff", str(cutoff),
                      "--assemble", *map(str, cases)])


def test_joint_assembly_uses_both_measured_blocks(workdir):
    a, b, out = [workdir / name for name in ("a", "b", "out")]
    _case(a, 0)
    _case(b, 1)
    assert _assemble(out, b, a) == 0
    block = np.load(out / "lead_fc.npy")
    assert np.allclose(block, np.kron([[1, 5], [5, 11]], np.eye(3)))
    report = _json(out / "lead_asr.json")
    assert report["raw_reciprocity_max_ev_per_A2"] == pytest.approx(4.0)
    assert report["raw_intralayer_AA_reciprocity_max_ev_per_A2"] == 0.0
    assert report["raw_intralayer_BB_reciprocity_max_ev_per_A2"] == 0.0
    assert report["raw_interlayer_reciprocity_max_ev_per_A2"] == pytest.approx(4.0)
    assert report["asr_residual_full_cell_ev_per_A2"] == pytest.approx(38.0)
    assert report["measured_layers"] == [0, 1]


def test_only_a_retains_approximation_and_directional_asr(workdir):
    _case(workdir / "a", 0)
    out = workdir / "result"
    _assemble(out, workdir / "a")
    assert np.allclose(np.load(out / "lead_fc.npy"),
                       np.kron([[1, 3], [3, 1]], np.eye(3)))
    assert "approximation" in _json(out / "lead_fc.json")["source"]
    report = _json(out / "lead_asr.json")
    assert report["raw_reciprocity_max_ev_per_A2"] is None
    assert report["raw_interlayer_reciprocity_max_ev_per_A2"] is None
    assert report["raw_intralayer_reciprocity_max_ev_per_A2"] is None
    assert report["asr_residual_full_cell_ev_per_A2"] == pytest.approx(24.0)


def test_b_alone_fails_without_asset(workdir):
    _case(workdir / "b", 1)
    with pytest.raises((ValueError, SystemExit)):
        _assemble(workdir / "out", workdir / "b")
    assert not (workdir / "out").exists()


def test_asr_keeps_cartesian_translation_directions_separate():
    matrix = np.array([[1.0, -1.0, 0.0],
                       [-1.0, 1.0, 0.0],
                       [0.0, 0.0, 0.0]])
    assert np.max(np.abs(matrix.sum(axis=1))) == 0.0
    assert lead._directional_asr(matrix) == 1.0


@pytest.mark.parametrize("fault", ["control", "geometry", "duplicate_layer", "slot",
                                   "metadata", "missing_axis", "sign", "amplitude"])
def test_rejects_incompatible_or_malformed_cases(workdir, fault):
    a, b, out = [workdir / name for name in ("a", "b", "out")]
    _case(a, 0)
    _case(b, 0 if fault == "duplicate_layer" else 1)
    if fault == "control":
        (b / "control.in").write_text("xc pw-lda\n", encoding="utf-8")
    elif fault == "geometry":
        geom = b / "geometry.in"
        geom.write_text(geom.read_text().replace("-5.88641797", "-6.88641797"))
    elif fault == "metadata":
        meta = _json(b / "lead_case.json")
        meta["surface_layer_indices"] = [0]
        (b / "lead_case.json").write_text(json.dumps(meta))
    elif fault != "duplicate_layer":
        manifest = _json(b / "displacements.json")
        if fault == "slot":
            manifest[0]["slot"] = 1
        elif fault == "missing_axis":
            manifest[0]["axis"] = 2
        elif fault == "sign":
            manifest[0]["sign"] = -1
        elif fault == "amplitude":
            for rec in manifest:
                rec["delta_A"] *= 2
        (b / "displacements.json").write_text(json.dumps(manifest))
    with pytest.raises((ValueError, SystemExit)):
        _assemble(out, a, b)
    assert not out.exists()


def test_cutoff_uses_actual_two_layer_positions_and_hexagonal_images(workdir):
    path, out = workdir / "a", workdir / "out"
    _case(path, 0, (3, 3))
    geometry = path / "geometry.in"
    geometry.write_text(geometry.read_text().replace("-3.53185078", "-3.1"))
    _assemble(out, path, cutoff=4.5)
    lattice, all_positions = lead._geometry(geometry)
    positions = all_positions[:18]
    shifts = np.array(list(itertools.product(range(-2, 3), repeat=2))) @ lattice[:2]
    delta = positions[:, None, :] - positions[None, :, :]
    expected = np.linalg.norm(delta[:, :, None, :] + shifts, axis=-1).min(axis=-1) <= 4.5
    result = np.load(out / "lead_fc.npy").reshape(18, 3, 18, 3)
    np.testing.assert_array_equal(result, np.broadcast_to(expected[:, None, :, None], result.shape))


def test_prepare_b_ships_all_three_helpers(workdir, capsys):
    species = Path(__file__).resolve().parents[2] / "demo/azulene-26-dithiol/species"
    assert lead.main(["--out", str(workdir), "--species-dir", str(species),
                      "--reference-layer", "1", "--in-plane", "1", "1"]) == 0
    assert 'cd "$HOME/smj-runs/au-lead-b/lead_fc" && bash sub_lead_fc.sh --submit' in capsys.readouterr().out
    assert {d["atom"] for d in _json(workdir / "displacements.json")} == {1}
    assert (workdir / "concat_forces.py").is_file()
    for name in ["apply_displacement.py", "extract_forces.py", "concat_forces.py"]:
        assert (workdir / name).read_bytes() == (lead._ROOT / "src/smj/phonon" / name).read_bytes()


@pytest.mark.parametrize("location,allowed", [
    ("smj-runs/Au-A/lead_fc", True), ("archive/Au-A/lead_fc", False),
    ("smj-runs/Au-A/extra/lead_fc", False), ("smj-runs/Au-A/thermal_fc", False),
    ("smj-runs/Au--A/lead_fc", False), ("smj-runs/" + "A" * 41 + "/lead_fc", False),
    ("smj-runs/Au-B/lead_fc", False),  # fc is a link to an external directory
    ("smj-runs/link/lead_fc", False),  # whole run is linked outside smj-runs
])
def test_lead_cluster_workdir(workdir, location, allowed):
    bash = shutil.which("bash")
    if os.name == "nt" and shutil.which("git"):
        bash = Path(shutil.which("git")).parents[1] / "bin/bash.exe"
    if not bash or not Path(bash).is_file():
        pytest.skip("requires bash to execute the Slurm script")
    home = workdir.resolve()
    stage = home / location
    if "Au-B" in location or "/link/" in location:
        target = home / "external"
        target.mkdir()
        link = stage / "fc" if "Au-B" in location else stage.parent
        link.parent.mkdir(parents=True)
        if os.name == "nt":
            import _winapi
            _winapi.CreateJunction(str(target), str(link))
        else:
            link.symlink_to(target, target_is_directory=True)
    stage.mkdir(parents=True, exist_ok=True)
    species = lead._ROOT / "demo/azulene-26-dithiol/species"
    lead.main(["--out", str(stage), "--species-dir", str(species), "--in-plane", "1", "1",
               "--module-load", "echo worker-started; exit 0"])
    # Use a stand-in sbatch and stop the worker at module load, before running any calculation.
    stub = home / "sbatch"
    stub.write_bytes(b"#!/bin/bash\necho scheduler-called\n")
    stub.chmod(0o755)
    env = dict(os.environ, HOME=home.as_posix(), PATH=str(home) + os.pathsep + os.environ["PATH"])
    for arg, marker in [("--submit", "scheduler-called"), ("", "worker-started")]:
        result = subprocess.run(
            [bash, "-c", f'export SLURM_SUBMIT_DIR="$PWD"; bash sub_lead_fc.sh {arg}'],
            cwd=stage, env=env, capture_output=True, text=True)
        assert (result.returncode == 0) == allowed, result.stderr
        assert (marker in result.stdout) == allowed
        if not allowed:
            assert "lead_fc" in result.stderr
    assert not (stage / "fc/0000").exists()
    if allowed:
        result = subprocess.run([bash, "sub_lead_fc.sh", "--submit", "--chdir=/tmp"],
                                cwd=stage, env=env, capture_output=True, text=True)
        assert result.returncode != 0 and "scheduler-called" not in result.stdout


def test_rectangular_cell_scales_abc_shift_and_layer_spacing():
    atoms, lattice = lead.build_bulk_cell(6, (2, 3), a_lattice=5.0)
    a1, a2 = lead.in_plane_vectors(5.0, (2, 3))
    expected_shift = a1 / 6 + a2 / 9
    assert np.allclose(atoms[6].coord(), np.r_[expected_shift[:2], atoms[0].z + 5 / np.sqrt(3)])
    assert lattice[2][2] == pytest.approx(6 * 5 / np.sqrt(3))


def test_assemble_accepts_finalized_bundle_and_keeps_scf_provenance(workdir):
    case = workdir / "a"
    rows = _case(case, 0)
    forces = []
    for rec in _json(case / "displacements.json"):
        path = case / "fc" / ("%04d" % rec["slot"]) / "forces.json"
        forces.append(_json(path)["forces"])
        path.unlink()
    slots = [{"slot": slot, "converged": True, "n_cycles": 1}
             for slot in range(len(forces))]
    summary = {"count": len(forces), "slots": slots,
               "missing": [], "not_converged": [], "unknown": []}
    (case / "forces_all.json").write_text(json.dumps(
        {"n_atoms": 6, "count": len(forces), "forces": forces, "summary": summary}),
        encoding="utf-8")
    _assemble(workdir / "out", case)
    report = _json(workdir / "out" / "lead_asr.json")
    assert report["force_inputs_by_layer"]["0"]["scf_evidence"] == "converged (summary)"
    assert report["raw_interlayer_reciprocity_max_ev_per_A2"] is None
    assert report["raw_intralayer_reciprocity_max_ev_per_A2"] is None


def test_assemble_rejects_failed_scf_summary(workdir):
    case = workdir / "a"
    _case(case, 0)
    forces = [_json(case / "fc" / ("%04d" % slot) / "forces.json")["forces"]
              for slot in range(6)]
    summary = {"count": 6,
               "slots": [{"slot": slot, "converged": slot != 2}
                         for slot in range(6)],
               "missing": [], "not_converged": [2], "unknown": []}
    (case / "forces_all.json").write_text(json.dumps(
        {"n_atoms": 6, "count": 6, "forces": forces, "summary": summary}),
        encoding="utf-8")
    with pytest.raises(StageError):
        _assemble(workdir / "out", case)
    assert not (workdir / "out").exists()


def test_assemble_reads_mixed_legacy_force_json_and_aims_logs(workdir):
    case = workdir / "a"
    rows = _case(case, 0)
    for rec in _json(case / "displacements.json"):
        path = case / "fc" / ("%04d" % rec["slot"])
        forces = _json(path / "forces.json")["forces"]
        data = ["Self-consistency cycle converged",
                "Total atomic forces (unitary forces cleaned) [eV/Ang]:"]
        data.extend("| %d % .14e % .14e % .14e" %
                    (i + 1, *force) for i, force in enumerate(forces))
        (path / "aims.out").write_text("\n".join(data), encoding="utf-8")
        if rec["slot"] % 2:
            (path / "forces.json").unlink()
    _assemble(workdir / "out", case)
    report = _json(workdir / "out" / "lead_asr.json")
    assert report["force_inputs_by_layer"]["0"]["source"] == "mixed aims.out+forces.json"
    assert report["force_inputs_by_layer"]["0"]["scf_evidence"] == "converged (aims.out)"
