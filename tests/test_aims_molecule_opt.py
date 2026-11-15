"""分子预优化引擎：FHI-aims 默认、Gaussian 兼容。

用户口径（2026-11-13）：``gaussian_opt`` 默认改为 FHI-aims PBE/light，让
分子优化、tip 弛豫、SCF 在同一软件栈内完成；旧的 Gaussian 分支保留，但不再是
默认。两种引擎都必须把优化后几何落成 ``opt.xyz`` 交给 ``aims_tip_relax``。
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from types import SimpleNamespace

import pytest

from smj.config import ElectrodeConfig, MoleculeConfig, load_config
from smj.errors import ConfigError, StageError
from smj.stages.molecule_opt import MoleculeOptStage

XYZ = """2
H2
H 0.0 0.0 0.0
H 0.74 0.0 0.0
"""

AIMS_OUT = """\
FHI-aims version      : 221103
aims_uuid             : TEST-UUID
Commit number         : TEST-COMMIT

  Convergence accuracy for geometry relaxation: Maximum force < 0.100000E-01 eV/A.
  Self-consistency cycle converged.
  Relaxation step number     1: Predicting new coordinates.
  Self-consistency cycle converged.
  Relaxation step number     2: Predicting new coordinates.
  Self-consistency cycle converged.
  Net remaining forces (excluding translations, rotations) in present geometry:
  Maximum force component is  0.500000E-02 eV/A.
  Present geometry is converged.

  Writing Kohn-Sham eigenvalues.
  State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
      1       2.00000          -0.500000         -13.60569
      2       0.00000          -0.100000          -2.72114

Have a nice day.
"""

#: 未跑完弛豫就正常退出（FHI-aims 达到步数上限时的真实形态：仍写 next_step）。
AIMS_OUT_NOT_CONVERGED = AIMS_OUT.replace(
    "  Present geometry is converged.",
    "  Present geometry is not yet converged.")

#: 有一步 SCF 没有自洽收敛：那一步的力不可信（旧版 relax 模板不会写 sc_accuracy_*）。
AIMS_OUT_SCF_FAIL = AIMS_OUT.replace(
    "  Self-consistency cycle converged.",
    "  Self-consistency cycle NOT converged.", 1)

#: 新版 relax 模板显式写出的 SCF 精度块。
CONTROL_WITH_SCF = """\
  xc                 pbe
  spin               none
  charge             0.0

  sc_accuracy_rho    1E-5
  sc_accuracy_eev    1E-3
  sc_accuracy_etot   1E-6
  sc_accuracy_forces 5E-3
  sc_iter_limit      1000
"""

AIMS_GEOM = """\
  atom  0.00000000  0.00000000  0.00000000 H
  atom  0.74000000  0.00000000  0.00000000 H
"""

#: 优化后一枚 H 飞走（H-H 键断了，且该原子没有任何键）。
AIMS_GEOM_DISSOCIATED = """\
  atom  0.00000000  0.00000000  0.00000000 H
  atom  5.00000000  0.00000000  0.00000000 H
"""

#: 两枚 H 塌到一起（距离 0.1 Å，低于 0.7 Å 重叠阈值）。
AIMS_GEOM_OVERLAP = """\
  atom  0.00000000  0.00000000  0.00000000 H
  atom  0.10000000  0.00000000  0.00000000 H
"""


def _species(root: Path) -> Path:
    species = root / "species"
    species.mkdir(parents=True, exist_ok=True)
    for element in ("H", "C", "N", "S", "Au"):
        (species / ("%s_default" % element)).write_text(
            "  species %s\n    # dummy\n" % element, encoding="utf-8")
    return species


def _ctx(root: Path, *, engine: str = "aims", multiplicity: int = 1,
         charge: int = 0) -> SimpleNamespace:
    run_dir = root / "run"
    run_dir.mkdir(parents=True, exist_ok=True)
    xyz = root / "molecule.xyz"
    xyz.write_text(XYZ, encoding="utf-8")
    species = _species(root)
    molecule = MoleculeConfig(
        xyz=str(xyz), opt_engine=engine, opt_functional="pbe",
        opt_basis="light", opt_force_threshold="1.e-2",
        charge=charge, multiplicity=multiplicity,
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
    stage_dir = workdir / "run" / "molecule_opt"
    prepared = MoleculeOptStage().prepare(ctx, stage_dir)

    assert prepared.script == "sub_aims.sh"
    assert (stage_dir / "geometry.in").exists()
    assert (stage_dir / "control.in").exists()
    assert (stage_dir / "sub_aims.sh").exists()
    assert not (stage_dir / "mol.gjf").exists()
    control = (stage_dir / "control.in").read_text(encoding="utf-8")
    assert "relax_geometry" in control and "bfgs 1.e-2" in control
    assert "pbe" in control
    assert "spin" in control and "none" in control
    # 与 control_scf.in 同源的 SCF 精度必须显式写出来（旧版没有）。
    assert "sc_accuracy_rho    1E-5" in control
    assert "sc_accuracy_etot   1E-6" in control
    assert "sc_accuracy_forces 5E-3" in control
    assert "sc_iter_limit      1000" in control
    # 闭壳层不写 fixed_spin_moment（与旧版行为一致）。
    assert "fixed_spin_moment" not in control
    summary = json.loads((stage_dir / "input_summary.json").read_text(encoding="utf-8"))
    assert summary["engine"] == "aims"
    assert summary["xc"] == "pbe"
    assert summary["spin"] == "none"
    assert summary["spin_moment"] is None
    assert summary["scf"]["sc_accuracy_forces"] == "5E-3"
    assert summary["scf"]["sc_iter_limit"] == 1000


def test_aims_prepare_maps_multiplicity_to_spin(workdir):
    ctx = _ctx(workdir, multiplicity=3)
    stage_dir = workdir / "run" / "molecule_opt"
    MoleculeOptStage().prepare(ctx, stage_dir)
    control = (stage_dir / "control.in").read_text(encoding="utf-8")
    assert "spin" in control and "2" in control
    # 开壳层把总自旋固定住（μB = M-1），避免弛豫步间自旋态漂移。
    assert "fixed_spin_moment  2" in control
    summary = json.loads((stage_dir / "input_summary.json").read_text(encoding="utf-8"))
    assert summary["spin_moment"] == 2.0


def test_aims_collect_writes_opt_xyz_and_records_software(workdir):
    ctx = _ctx(workdir)
    stage_dir = workdir / "run" / "molecule_opt"
    stage_dir.mkdir(parents=True, exist_ok=True)
    (stage_dir / "geometry.in").write_text(AIMS_GEOM, encoding="utf-8")
    (stage_dir / "control.in").write_text(CONTROL_WITH_SCF, encoding="utf-8")
    (stage_dir / "aims.out").write_text(AIMS_OUT, encoding="utf-8")
    (stage_dir / "geometry.in.next_step").write_text(AIMS_GEOM, encoding="utf-8")

    collected = MoleculeOptStage().collect(ctx, stage_dir)

    opt = stage_dir / "opt.xyz"
    assert opt.exists()
    lines = opt.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "2"
    observed = json.loads((ctx.layout.run_dir / "software_observed.json").read_text(
        encoding="utf-8"))
    assert observed["molecule_opt"]["version"] == "221103"
    assert collected.results["software"]["version"] == "221103"

    result = collected.results["molecule_opt"]
    assert result["converged"] is True
    assert result["relax_steps"] == 2
    assert result["force_threshold_ev_per_A"] == pytest.approx(0.01)
    assert result["max_force_ev_per_A"] == pytest.approx(0.005)
    check = result["geometry_check"]
    assert check["ok"] is True
    assert check["atom_count"] == 2
    assert check["min_distance_A"] == pytest.approx(0.74)
    assert check["bond_count"] == [1, 1]
    assert check["lost_bonds"] == [] and check["new_bonds"] == []
    orbitals = result["frontier_orbitals"]
    assert orbitals["n_states"] == 2
    assert orbitals["homo_ev"] == pytest.approx(-13.60569)
    assert orbitals["lumo_ev"] == pytest.approx(-2.72114)
    assert orbitals["gap_ev"] == pytest.approx(10.88455)
    assert orbitals["fractional_occupation"] is False
    assert any("10.885 eV" in note for note in collected.notes)
    scf = result["scf"]
    assert scf["sc_accuracy_forces"] == "5E-3"
    assert scf["sc_accuracy_rho"] == "1E-5"
    assert scf["sc_iter_limit"] == "1000"
    assert scf["relax_steps"] == 2
    assert scf["converged_steps"] == 3
    assert "missing_steps" not in scf
    assert any("SCF 精度" in note and "5E-3 eV/Å" in note for note in collected.notes)


def test_aims_collect_rejects_unconverged_scf_step(workdir):
    ctx = _ctx(workdir)
    stage_dir = workdir / "run" / "molecule_opt"
    stage_dir.mkdir(parents=True, exist_ok=True)
    (stage_dir / "geometry.in").write_text(AIMS_GEOM, encoding="utf-8")
    (stage_dir / "control.in").write_text(CONTROL_WITH_SCF, encoding="utf-8")
    (stage_dir / "aims.out").write_text(AIMS_OUT_SCF_FAIL, encoding="utf-8")
    (stage_dir / "geometry.in.next_step").write_text(AIMS_GEOM, encoding="utf-8")

    with pytest.raises(StageError, match="SCF 未收敛"):
        MoleculeOptStage().collect(ctx, stage_dir)
    assert not (stage_dir / "opt.xyz").exists()


def test_aims_collect_flags_old_control_without_scf_accuracy(workdir):
    """旧 run 的 control.in 没有 sc_accuracy_*：不假装有，给出显式提醒。"""
    ctx = _ctx(workdir)
    stage_dir = workdir / "run" / "molecule_opt"
    stage_dir.mkdir(parents=True, exist_ok=True)
    (stage_dir / "geometry.in").write_text(AIMS_GEOM, encoding="utf-8")
    (stage_dir / "aims.out").write_text(AIMS_OUT, encoding="utf-8")
    (stage_dir / "geometry.in.next_step").write_text(AIMS_GEOM, encoding="utf-8")

    collected = MoleculeOptStage().collect(ctx, stage_dir)
    scf = collected.results["molecule_opt"]["scf"]
    assert "sc_accuracy_forces" not in scf
    assert any("未显式设置 sc_accuracy_forces" in note for note in collected.notes)


def test_aims_collect_rejects_unconverged_relaxation(workdir):
    ctx = _ctx(workdir)
    stage_dir = workdir / "run" / "molecule_opt"
    stage_dir.mkdir(parents=True, exist_ok=True)
    (stage_dir / "aims.out").write_text(AIMS_OUT_NOT_CONVERGED, encoding="utf-8")
    (stage_dir / "geometry.in.next_step").write_text(AIMS_GEOM, encoding="utf-8")

    with pytest.raises(StageError, match="未收敛"):
        MoleculeOptStage().collect(ctx, stage_dir)
    assert not (stage_dir / "opt.xyz").exists()


def test_aims_collect_rejects_dissociated_geometry(workdir):
    ctx = _ctx(workdir)
    stage_dir = workdir / "run" / "molecule_opt"
    stage_dir.mkdir(parents=True, exist_ok=True)
    (stage_dir / "geometry.in").write_text(AIMS_GEOM, encoding="utf-8")
    (stage_dir / "aims.out").write_text(AIMS_OUT, encoding="utf-8")
    (stage_dir / "geometry.in.next_step").write_text(
        AIMS_GEOM_DISSOCIATED, encoding="utf-8")

    with pytest.raises(StageError, match="孤立原子"):
        MoleculeOptStage().collect(ctx, stage_dir)
    assert not (stage_dir / "opt.xyz").exists()


def test_aims_collect_rejects_overlapping_geometry(workdir):
    ctx = _ctx(workdir)
    stage_dir = workdir / "run" / "molecule_opt"
    stage_dir.mkdir(parents=True, exist_ok=True)
    (stage_dir / "geometry.in").write_text(AIMS_GEOM, encoding="utf-8")
    (stage_dir / "aims.out").write_text(AIMS_OUT, encoding="utf-8")
    (stage_dir / "geometry.in.next_step").write_text(
        AIMS_GEOM_OVERLAP, encoding="utf-8")

    with pytest.raises(StageError, match="重叠"):
        MoleculeOptStage().collect(ctx, stage_dir)
    assert not (stage_dir / "opt.xyz").exists()


def test_aims_collect_without_geometry_in_skips_bond_graph(workdir):
    ctx = _ctx(workdir)
    stage_dir = workdir / "run" / "molecule_opt"
    stage_dir.mkdir(parents=True, exist_ok=True)
    (stage_dir / "aims.out").write_text(AIMS_OUT, encoding="utf-8")
    (stage_dir / "geometry.in.next_step").write_text(AIMS_GEOM, encoding="utf-8")

    collected = MoleculeOptStage().collect(ctx, stage_dir)

    check = collected.results["molecule_opt"]["geometry_check"]
    assert check["ok"] is True
    assert check["bond_count"] is None and check["rmsd_A"] is None
    assert any("未找到 prepare 写的 geometry.in" in note for note in collected.notes)


def test_aims_collect_fails_closed_on_missing_next_step(workdir):
    ctx = _ctx(workdir)
    stage_dir = workdir / "run" / "molecule_opt"
    stage_dir.mkdir(parents=True, exist_ok=True)
    (stage_dir / "aims.out").write_text(AIMS_OUT, encoding="utf-8")

    with pytest.raises(StageError, match="next_step"):
        MoleculeOptStage().collect(ctx, stage_dir)


def test_aims_collect_fails_closed_on_abnormal_termination(workdir):
    ctx = _ctx(workdir)
    stage_dir = workdir / "run" / "molecule_opt"
    stage_dir.mkdir(parents=True, exist_ok=True)
    (stage_dir / "aims.out").write_text(
        "FHI-aims version : 221103\n", encoding="utf-8")
    (stage_dir / "geometry.in.next_step").write_text(AIMS_GEOM, encoding="utf-8")

    with pytest.raises(StageError, match="Have a nice day"):
        MoleculeOptStage().collect(ctx, stage_dir)

def _write_config(root: Path, threshold) -> Path:
    cfg = root / "run.yaml"
    cfg.write_text(
        "run:\n  name: t\n  dir: .\n  mode: dry-run\n"
        "molecule:\n  xyz: molecule.xyz\n  opt_force_threshold: %r\n" % threshold,
        encoding="utf-8")
    return cfg


def test_config_accepts_default_force_threshold(workdir):
    cfg = load_config(_write_config(workdir, "1.e-2"))
    assert cfg.molecule.opt_force_threshold == "1.e-2"


def test_config_rejects_force_threshold_above_cap(workdir):
    """0.2 eV/Å 比默认松 20 倍，几乎等于不收敛——只能在 (0, 0.1] 内。"""
    with pytest.raises(ConfigError, match="opt_force_threshold"):
        load_config(_write_config(workdir, "0.2"))


def test_config_rejects_non_numeric_force_threshold(workdir):
    with pytest.raises(ConfigError, match="opt_force_threshold"):
        load_config(_write_config(workdir, "loose"))
def _fake_gaussian_log(monkeypatch, geometry: str):
    """把 parse_gaussian_log 换掉：只验证 collect 的几何闸门接线，不碰真日志。"""
    from smj.stages import molecule_opt as module

    def fake_parse(log, xyz_path, summary_path):
        Path(xyz_path).write_text(geometry, encoding="utf-8")
        Path(summary_path).write_text("{}", encoding="utf-8")
        return {"normal_termination": True, "imaginary_frequency_count": 0,
                "atom_count": len(geometry.strip().splitlines()) - 2,
                "final_scf_energy_hartree": -1.0}

    monkeypatch.setattr(module, "parse_gaussian_log", fake_parse)


def test_gaussian_collect_runs_geometry_check(workdir, monkeypatch):
    ctx = _ctx(workdir, engine="gaussian")
    stage_dir = workdir / "run" / "molecule_opt"
    stage_dir.mkdir(parents=True, exist_ok=True)
    (stage_dir / "mol.log").write_text("fake\n", encoding="utf-8")
    _fake_gaussian_log(monkeypatch, XYZ)

    collected = MoleculeOptStage().collect(ctx, stage_dir)

    check = collected.results["gaussian"]["geometry_check"]
    assert check["ok"] is True
    assert check["bond_count"] == [1, 1]


def test_gaussian_collect_rejects_dissociated_geometry(workdir, monkeypatch):
    ctx = _ctx(workdir, engine="gaussian")
    stage_dir = workdir / "run" / "molecule_opt"
    stage_dir.mkdir(parents=True, exist_ok=True)
    (stage_dir / "mol.log").write_text("fake\n", encoding="utf-8")
    _fake_gaussian_log(monkeypatch, "2\nH2\nH 0.0 0.0 0.0\nH 5.0 0.0 0.0\n")

    with pytest.raises(StageError, match="孤立原子"):
        MoleculeOptStage().collect(ctx, stage_dir)
    assert not (stage_dir / "opt.xyz").exists()

def test_prepare_rejects_charge_multiplicity_parity_mismatch(workdir):
    """H2 带 +1 电荷 = 1 个电子，却写 multiplicity=1（要求偶多重度）。"""
    ctx = _ctx(workdir, charge=1, multiplicity=1)
    with pytest.raises(StageError, match="多重度"):
        MoleculeOptStage().prepare(ctx, workdir / "run" / "molecule_opt")


def test_prepare_accepts_open_shell_cation(workdir):
    """同一个 H2+ 写 multiplicity=2 就成立，并把电子数记进 notes。"""
    ctx = _ctx(workdir, charge=1, multiplicity=2)
    prepared = MoleculeOptStage().prepare(ctx, workdir / "run" / "molecule_opt")
    assert any("电子数 1（charge=1，multiplicity=2）" in note
               for note in prepared.notes)
    control = (workdir / "run" / "molecule_opt" / "control.in").read_text(
        encoding="utf-8")
    spin_line = [line for line in control.splitlines()
                 if line.strip().startswith("spin")][0]
    assert spin_line.split()[-1] == "1"


def _write_cube(path: Path, *, counts=(11, 11, 11), step_a=0.2, sigma=0.45):
    """写一个 S 原子在原点、轨道密度为高斯球的 cube（坐标按 bohr）。"""
    nx, ny, nz = counts
    origin = (-1.0, -1.0, -1.0)
    lines = ["cube test", "generated by tests"]
    lines.append("1 %.6f %.6f %.6f" % tuple(v / 0.529177210903 for v in origin))
    for axis, n in enumerate(counts):
        vector = [0.0, 0.0, 0.0]
        vector[axis] = step_a / 0.529177210903
        lines.append("%d %.6f %.6f %.6f" % (n, *vector))
    lines.append("16 0.000000 %10.6f %10.6f %10.6f" % tuple(
        v / 0.529177210903 for v in (0.0, 0.0, 0.0)))
    chunk = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                x = origin[0] + i * step_a
                y = origin[1] + j * step_a
                z = origin[2] + k * step_a
                rho = math.exp(-(x * x + y * y + z * z) / (2.0 * sigma ** 2))
                chunk.append("%14.8E" % rho)
                if len(chunk) == 6:
                    lines.append(" ".join(chunk))
                    chunk = []
    if chunk:
        lines.append(" ".join(chunk))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_aims_collect_records_anchor_orbital_weights(workdir):
    """有 cube 时 collect 把锚定原子(S)上的轨道占比写进 results 与 notes。"""
    ctx = _ctx(workdir)
    (workdir / "molecule.xyz").write_text("1\nS\nS 0.0 0.0 0.0\n", encoding="utf-8")
    stage_dir = workdir / "run" / "molecule_opt"
    stage_dir.mkdir(parents=True, exist_ok=True)
    geom = "  atom  0.00000000  0.00000000  0.00000000 S\n"
    (stage_dir / "geometry.in").write_text(geom, encoding="utf-8")
    (stage_dir / "geometry.in.next_step").write_text(geom, encoding="utf-8")
    (stage_dir / "aims.out").write_text(AIMS_OUT, encoding="utf-8")
    (stage_dir / "control.in").write_text(
        "xc pbe\noutput cube eigenstate HOMO\n", encoding="utf-8")
    _write_cube(stage_dir / "cube_001_eigenstate_00008_spin_1.cube")

    collected = MoleculeOptStage().collect(ctx, stage_dir)

    result = collected.results["molecule_opt"]
    assert result["electron_count"] == 16
    weights = result["frontier_orbital_weights"]
    assert weights[0]["orbit"] == "HOMO"
    assert weights[0]["anchor_fraction"] > 0.9
    assert weights[0]["anchors"][0]["element"] == "S"
    assert any("锚定轨道权重" in note and "HOMO" in note for note in collected.notes)