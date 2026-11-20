"""开壳层（自由基）结链：配置解析 + 自旋输入接线。

2026-11-19：FHI-aims 的 ``spin`` 只接受 ``none``/``collinear``，多重度由
``fixed_spin_moment = M−1``（=2S）表达；结链路（tip_relax/aims_scf/aitranss_te）
必须把 ``molecule.charge/multiplicity`` 一路带下去，否则自由基静默退化成闭壳层结。
"""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from smj.config import ElectrodeConfig, MoleculeConfig, load_config
from smj.errors import ConfigError
from smj.stages.aims_scf import AimsScfStage
from smj.stages.aims_tip_relax import AimsTipRelaxStage


def _config(root: Path, molecule_yaml: str):
    cfg = root / "run.yaml"
    cfg.write_text(
        "run:\n  name: t\n  mode: dry-run\n"
        "molecule:\n  xyz: molecule.xyz\n" + molecule_yaml,
        encoding="utf-8")
    return load_config(cfg)


def test_initial_moments_parse_one_based_indices(workdir):
    config = _config(workdir, "  multiplicity: 2\n"
                              "  initial_moments:\n"
                              "    7: 1.0\n"
                              '    "9": 0.5\n')
    assert config.molecule.initial_moments == {7: 1.0, 9: 0.5}


def test_initial_moments_absent_is_none(workdir):
    assert _config(workdir, "").molecule.initial_moments is None


def test_initial_moments_reject_zero_index(workdir):
    with pytest.raises(ConfigError, match="1-based"):
        _config(workdir, "  initial_moments:\n    0: 1.0\n")


def test_initial_moments_reject_non_mapping(workdir):
    with pytest.raises(ConfigError, match="molecule.initial_moments"):
        _config(workdir, "  initial_moments: 1.0\n")


# --------------------------------------------------------------------------- #
# tip_relax：开壳层/带电 -> control.in 关键字
# --------------------------------------------------------------------------- #
def _ctx(workdir: Path, molecule: MoleculeConfig, electrode=None):
    species = workdir / "species"
    species.mkdir(exist_ok=True)
    xyz = workdir / "opt.xyz"
    xyz.write_text("2\nTTF\nS 0.0 0.0 0.0\nS 0.0 0.0 3.0\n", encoding="utf-8")
    return SimpleNamespace(
        config=SimpleNamespace(
            molecule=molecule,
            electrode=electrode or ElectrodeConfig()),
        override_path=lambda key: None,
        upstream=lambda stage, name, override=None: xyz,
        artifact=lambda stage, name: workdir / "run" / stage / name,
        require_stage_file=lambda stage, name: workdir / "run" / stage / name,
        species_dir=lambda: species,
        estimate=lambda name, natoms: {"atoms": natoms},
        stage_env=lambda name, estimate: {},
        log=lambda message: None,
    )


def _fake_tip_input(captured: dict):
    def fake(xyz_path, outdir, **kwargs):
        outdir = Path(outdir)
        captured.update(kwargs)
        geometry = outdir / "geometry.in"
        geometry.write_text(
            "atom 0.0 0.0 0.0 S\natom 0.0 0.0 3.0 S\n"
            "atom 0.0 0.0 -2.35 Au\natom 0.0 0.0 5.35 Au\n",
            encoding="utf-8")
        control = outdir / "control.in"
        control.write_text("xc pbe\nspin collinear\n  species S\n", encoding="utf-8")
        summary = outdir / "input_summary.json"
        summary.write_text(json.dumps({
            "left_anchor": 1, "right_anchor": 2, "anchor_element": "S",
            "contact_pinned": False, "placement_distance_A": 2.35,
            "strip_thiol_h": {"removed": 0}, "tip_element": "Au",
        }), encoding="utf-8")
        return {"geometry": str(geometry), "control": str(control),
                "submit_script": str(outdir / "sub_aims.sh"), "summary": json.loads(
                    summary.read_text(encoding="utf-8"))}
    return fake


def test_tip_relax_open_shell_writes_spin_keywords(workdir, monkeypatch):
    captured: dict = {}
    monkeypatch.setattr("smj.stages.aims_tip_relax.make_tip_relax_input",
                        _fake_tip_input(captured))
    ctx = _ctx(workdir, MoleculeConfig(multiplicity=2, charge=0))
    stage_dir = workdir / "run" / "aims_tip_relax"
    stage_dir.mkdir(parents=True)
    prepared = AimsTipRelaxStage().prepare(ctx, stage_dir)
    assert captured["spin"] == "collinear"
    assert captured["charge"] == 0.0
    control = (stage_dir / "control.in").read_text(encoding="utf-8")
    assert "fixed_spin_moment  1" in control
    assert "default_initial_moment" in control
    assert "override_initial_charge_check" not in control
    summary = json.loads((stage_dir / "input_summary.json").read_text(encoding="utf-8"))
    assert summary["spin"] == "collinear"
    assert summary["fixed_spin_moment"] == 1.0
    assert any("开壳层结" in note for note in prepared.notes)


def test_tip_relax_charged_writes_initial_charge(workdir, monkeypatch):
    captured: dict = {}
    monkeypatch.setattr("smj.stages.aims_tip_relax.make_tip_relax_input",
                        _fake_tip_input(captured))
    ctx = _ctx(workdir, MoleculeConfig(multiplicity=2, charge=-1))
    stage_dir = workdir / "run" / "aims_tip_relax"
    stage_dir.mkdir(parents=True)
    AimsTipRelaxStage().prepare(ctx, stage_dir)
    assert captured["charge"] == -1.0
    control = (stage_dir / "control.in").read_text(encoding="utf-8")
    assert "override_initial_charge_check" not in control
    geometry = (stage_dir / "geometry.in").read_text(encoding="utf-8")
    charge_lines = [line for line in geometry.splitlines()
                    if line.split()[:1] == ["initial_charge"]]
    assert charge_lines
    assert sum(float(line.split()[1]) for line in charge_lines) == pytest.approx(-1.0)
    summary = json.loads((stage_dir / "input_summary.json").read_text(encoding="utf-8"))
    assert summary["initial_charge_total"] == pytest.approx(-1.0)


# --------------------------------------------------------------------------- #
# aims_scf：开壳层/带电 -> control.in 关键字
# --------------------------------------------------------------------------- #
def _fake_scf_input(captured: dict):
    def fake(geometry_in, outdir, **kwargs):
        outdir = Path(outdir)
        captured.update(kwargs)
        control = outdir / "control.in"
        control.write_text(
            "xc pbe\nspin %s\nfixed_spin_moment %g\ncharge %g\n  species S\n"
            % (kwargs["spin"], kwargs["fixed_spin_moment"], kwargs["charge"]),
            encoding="utf-8")
        geometry = outdir / "geometry.in"
        geometry.write_text("atom 0.0 0.0 0.0 S\n", encoding="utf-8")
        summary = outdir / "input_summary.json"
        summary.write_text(json.dumps({"spin": kwargs["spin"],
                                       "charge": kwargs["charge"],
                                       "fixed_spin_moment": kwargs["fixed_spin_moment"]}),
                           encoding="utf-8")
        return {"geometry": str(geometry), "control": str(control),
                "submit_script": str(outdir / "sub_aims.sh"),
                "summary": json.loads(summary.read_text(encoding="utf-8"))}
    return fake


def test_scf_open_shell_passes_spin_and_adds_initial_moment(workdir, monkeypatch):
    captured: dict = {}
    monkeypatch.setattr("smj.stages.aims_scf.make_aims_scf_input",
                        _fake_scf_input(captured))
    ctx = _ctx(workdir, MoleculeConfig(multiplicity=3, charge=0))
    stage_dir = workdir / "run" / "aims_scf"
    stage_dir.mkdir(parents=True)
    # aims_scf 需要上游 geometry.in：把 tip 弛豫的几何放进 workdir 供 fake 用
    prepared = AimsScfStage().prepare(ctx, stage_dir)
    assert captured["spin"] == "collinear"
    assert captured["fixed_spin_moment"] == 2.0
    control = (stage_dir / "control.in").read_text(encoding="utf-8")
    assert "default_initial_moment" in control


def test_scf_collect_accepts_alpha_beta_instead_of_mos(workdir):
    ctx = _ctx(workdir, MoleculeConfig(multiplicity=2))
    stage_dir = workdir / "run" / "aims_scf"
    stage_dir.mkdir(parents=True)
    (stage_dir / "aims.out").write_text("Have a nice day.\n", encoding="utf-8")
    (stage_dir / "basis-indices.out").write_text("1\n", encoding="utf-8")
    (stage_dir / "omat.aims").write_text("omat\n", encoding="utf-8")
    (stage_dir / "geometry.in").write_text(
        "atom 0.0 0.0 0.0 S\natom 3.0 0.0 0.0 S\n", encoding="utf-8")
    (stage_dir / "alpha.aims").write_text("nsaos = 10\n", encoding="utf-8")
    (stage_dir / "beta.aims").write_text("nsaos = 10\n", encoding="utf-8")
    collected = AimsScfStage().collect(ctx, stage_dir)
    payload = collected.results["aims_scf"]
    assert payload["has_alpha"] is True and payload["has_beta"] is True
    assert payload["has_mos"] is False
    assert payload["spin_unpolarized"] is False


# --------------------------------------------------------------------------- #
# 总自旋约束开关 + 结级初始磁矩
# --------------------------------------------------------------------------- #
def test_fixed_spin_moment_config(workdir):
    assert _config(workdir, "").molecule.fixed_spin_moment is True
    assert _config(workdir, "  fixed_spin_moment: false\n").molecule.fixed_spin_moment is False


def test_map_molecule_moments_shifts_around_removed_hydrogens():
    from smj.spin import map_molecule_moments
    # 分子原子 2 的 H 被移除：原 3、5 号原子在结几何里变成 2、4 号。
    assert map_molecule_moments({3: 1.0, 5: 2.0}, [2]) == {2: 1.0, 4: 2.0}
    with pytest.raises(ValueError, match="被移除"):
        map_molecule_moments({2: 1.0}, [2])
    with pytest.raises(ValueError, match="超出分子原子数"):
        map_molecule_moments({7: 1.0}, [], molecule_atoms=6)


def test_tip_relax_writes_per_atom_initial_moments(workdir, monkeypatch):
    captured: dict = {}
    monkeypatch.setattr("smj.stages.aims_tip_relax.make_tip_relax_input",
                        _fake_tip_input(captured))
    ctx = _ctx(workdir, MoleculeConfig(multiplicity=2, initial_moments={2: 0.5}))
    stage_dir = workdir / "run" / "aims_tip_relax"
    stage_dir.mkdir(parents=True)
    prepared = AimsTipRelaxStage().prepare(ctx, stage_dir)
    geometry = (stage_dir / "geometry.in").read_text(encoding="utf-8").splitlines()
    position = next(i for i, line in enumerate(geometry)
                    if line.split()[:1] == ["atom"] and line.split()[-1] == "S"
                    and i > 0)
    assert geometry[position + 1].split() == ["initial_moment", "0.5"]
    control = (stage_dir / "control.in").read_text(encoding="utf-8")
    assert "default_initial_moment" not in control
    summary = json.loads((stage_dir / "input_summary.json").read_text(encoding="utf-8"))
    assert summary["initial_moments"] == {"2": 0.5}
    assert any("初始磁矩" in note for note in prepared.notes)


def test_tip_relax_release_omits_fixed_spin_moment(workdir, monkeypatch):
    captured: dict = {}
    monkeypatch.setattr("smj.stages.aims_tip_relax.make_tip_relax_input",
                        _fake_tip_input(captured))
    ctx = _ctx(workdir, MoleculeConfig(multiplicity=2, fixed_spin_moment=False))
    stage_dir = workdir / "run" / "aims_tip_relax"
    stage_dir.mkdir(parents=True)
    AimsTipRelaxStage().prepare(ctx, stage_dir)
    control = (stage_dir / "control.in").read_text(encoding="utf-8")
    assert "fixed_spin_moment" not in control
    assert "default_initial_moment" in control
    summary = json.loads((stage_dir / "input_summary.json").read_text(encoding="utf-8"))
    assert summary["spin_constraint"] == "released"
    assert summary["fixed_spin_moment"] is None


def test_scf_release_removes_template_fixed_spin_moment(workdir, monkeypatch):
    captured: dict = {}
    monkeypatch.setattr("smj.stages.aims_scf.make_aims_scf_input",
                        _fake_scf_input(captured))
    ctx = _ctx(workdir, MoleculeConfig(multiplicity=2, fixed_spin_moment=False))
    stage_dir = workdir / "run" / "aims_scf"
    stage_dir.mkdir(parents=True)
    AimsScfStage().prepare(ctx, stage_dir)
    control = (stage_dir / "control.in").read_text(encoding="utf-8")
    assert "fixed_spin_moment" not in control
    assert "default_initial_moment" in control
