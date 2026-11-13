"""分子预优化引擎：FHI-aims 默认、Gaussian 兼容。

用户口径（2026-11-13）：``gaussian_opt`` 默认改为 FHI-aims PBE/light，让
分子优化、tip 弛豫、SCF 在同一软件栈内完成；旧的 Gaussian 分支保留，但不再是
默认。两种引擎都必须把优化后几何落成 ``opt.xyz`` 交给 ``aims_tip_relax``。
"""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from smj.config import ElectrodeConfig, MoleculeConfig
from smj.errors import StageError
from smj.stages.gaussian_opt import GaussianOptStage

XYZ = """2
H2
H 0.0 0.0 0.0
H 0.74 0.0 0.0
"""

AIMS_OUT = """\
FHI-aims version      : 221103
aims_uuid             : TEST-UUID
Commit number         : TEST-COMMIT

Have a nice day.
"""

AIMS_GEOM = """\
  atom  0.00000000  0.00000000  0.00000000 H
  atom  0.74000000  0.00000000  0.00000000 H
"""


def _species(root: Path) -> Path:
    species = root / "species"
    species.mkdir(parents=True, exist_ok=True)
    for element in ("H", "C", "N", "S", "Au"):
        (species / ("%s_default" % element)).write_text(
            "  species %s\n    # dummy\n" % element, encoding="utf-8")
    return species


def _ctx(root: Path, *, engine: str = "aims", multiplicity: int = 1) -> SimpleNamespace:
    run_dir = root / "run"
    run_dir.mkdir(parents=True, exist_ok=True)
    xyz = root / "molecule.xyz"
    xyz.write_text(XYZ, encoding="utf-8")
    species = _species(root)
    molecule = MoleculeConfig(
        xyz=str(xyz), opt_engine=engine, opt_functional="pbe",
        opt_basis="light", opt_force_threshold="1.e-2",
        charge=0, multiplicity=multiplicity,
    )
    config = SimpleNamespace(molecule=molecule, electrode=ElectrodeConfig(vdw=True))
    return SimpleNamespace(
        config=config,
        layout=SimpleNamespace(run_dir=run_dir),
        molecule_xyz=lambda: xyz,
        species_dir=lambda: species,
        estimate=lambda *a, **k: {
            "walltime": "00:30:00", "time_limit": "00:40:00",
            "mem": "4G", "ntasks": 2, "cores": 2, "capped": [],
        },
        stage_env=lambda *a, **k: {
            "partition": "main", "ntasks": 2, "aims_bin": "/bin/true",
            "aims_mem": "4G", "aims_walltime": "00:40:00",
        },
        require_stage_file=lambda stage, name: (root / "run" / stage / name),
        log=lambda *_: None,
    )


def test_default_engine_is_aims():
    assert MoleculeConfig().opt_engine == "aims"
    assert MoleculeConfig().opt_functional == "pbe"
    assert MoleculeConfig().opt_basis == "light"


def test_aims_prepare_writes_three_files_without_gaussian(workdir):
    ctx = _ctx(workdir)
    stage_dir = workdir / "run" / "gaussian_opt"
    prepared = GaussianOptStage().prepare(ctx, stage_dir)

    assert prepared.script == "sub_aims.sh"
    assert (stage_dir / "geometry.in").exists()
    assert (stage_dir / "control.in").exists()
    assert (stage_dir / "sub_aims.sh").exists()
    assert not (stage_dir / "mol.gjf").exists()
    control = (stage_dir / "control.in").read_text(encoding="utf-8")
    assert "relax_geometry" in control and "bfgs 1.e-2" in control
    assert "pbe" in control
    assert "spin" in control and "none" in control
    summary = json.loads((stage_dir / "input_summary.json").read_text(encoding="utf-8"))
    assert summary["engine"] == "aims"
    assert summary["xc"] == "pbe"
    assert summary["spin"] == "none"


def test_aims_prepare_maps_multiplicity_to_spin(workdir):
    ctx = _ctx(workdir, multiplicity=3)
    stage_dir = workdir / "run" / "gaussian_opt"
    GaussianOptStage().prepare(ctx, stage_dir)
    control = (stage_dir / "control.in").read_text(encoding="utf-8")
    assert "spin" in control and "2" in control


def test_aims_collect_writes_opt_xyz_and_records_software(workdir):
    ctx = _ctx(workdir)
    stage_dir = workdir / "run" / "gaussian_opt"
    stage_dir.mkdir(parents=True, exist_ok=True)
    (stage_dir / "aims.out").write_text(AIMS_OUT, encoding="utf-8")
    (stage_dir / "geometry.in.next_step").write_text(AIMS_GEOM, encoding="utf-8")

    collected = GaussianOptStage().collect(ctx, stage_dir)

    opt = stage_dir / "opt.xyz"
    assert opt.exists()
    lines = opt.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "2"
    observed = json.loads((ctx.layout.run_dir / "software_observed.json").read_text(
        encoding="utf-8"))
    assert observed["gaussian_opt"]["version"] == "221103"
    assert collected.results["software"]["version"] == "221103"


def test_aims_collect_fails_closed_on_missing_next_step(workdir):
    ctx = _ctx(workdir)
    stage_dir = workdir / "run" / "gaussian_opt"
    stage_dir.mkdir(parents=True, exist_ok=True)
    (stage_dir / "aims.out").write_text(AIMS_OUT, encoding="utf-8")

    with pytest.raises(StageError, match="next_step"):
        GaussianOptStage().collect(ctx, stage_dir)


def test_aims_collect_fails_closed_on_abnormal_termination(workdir):
    ctx = _ctx(workdir)
    stage_dir = workdir / "run" / "gaussian_opt"
    stage_dir.mkdir(parents=True, exist_ok=True)
    (stage_dir / "aims.out").write_text(
        "FHI-aims version : 221103\n", encoding="utf-8")
    (stage_dir / "geometry.in.next_step").write_text(AIMS_GEOM, encoding="utf-8")

    with pytest.raises(StageError, match="Have a nice day"):
        GaussianOptStage().collect(ctx, stage_dir)
