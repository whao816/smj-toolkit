"""molecule.family：判据标定族的配置解析与"不适用"判定。"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from smj.cli import main
from smj.config import load_config


TE = "".join([
    "#non-spin-polarized calculation\n",
    "#bias =  0.000000 V    efermi =  -0.188452434085 H\n",
    "#   E [Hartree]       E-EF [eV]       T(E) per spin\n",
    "  -0.1894000000     -0.02578138      0.7000000000E+00\n",
    "  -0.1886000000     -0.00401290      0.7400000000E+00\n",
    "  -0.1884000000      0.00142992      0.7500000000E+00\n",
    "  -0.1874000000      0.02864548      0.7600000000E+00\n",
])


def _project(workdir: Path, family: str | None) -> Path:
    project = (workdir / "proj").resolve()
    project.mkdir(parents=True)
    (project / "molecule.xyz").write_text("2\n\nC 0 0 0\nH 1 0 0\n", encoding="utf-8")
    lines = [
        "run:", "  name: t", "  mode: dry-run", "  dir: %s" % project.as_posix(),
        "molecule:", "  xyz: %s" % (project / "molecule.xyz").as_posix(),
    ]
    if family:
        lines.append("  family: %s" % family)
    lines += ["electrode:", "  source: builtin", "  au_per_electrode: 4"]
    cfg = project / "smj.yaml"
    cfg.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return cfg


def test_config_reads_molecule_family(workdir):
    cfg = _project(workdir, "AZ")
    config = load_config(cfg)
    assert config.molecule.family == "AZ"


def test_config_family_absent_is_none(workdir):
    cfg = _project(workdir, None)
    assert load_config(cfg).molecule.family is None


def test_family_marks_bdt_bands_not_applicable(workdir, capsys):
    project = (workdir / "run").resolve()
    stage = project / "aitranss_te"
    stage.mkdir(parents=True)
    (stage / "TE.dat").write_text(TE, encoding="utf-8")

    from smj.junction.validate import validate_junction_run
    verdict = validate_junction_run(project, system="AZ")
    names = {c["name"] for c in verdict["checks"] if c.get("applicable") is False}
    assert {"g_in_theory_band", "seebeck_sign"} <= names
    assert "g_in_theory_band" not in (verdict.get("failed_criteria") or [])

    unknown = validate_junction_run(project, system="NOPE")
    assert all(c.get("applicable") is not False for c in unknown["checks"]
               if c["name"] in ("g_in_theory_band",))


def test_bpy_family_registered():
    from smj.junction.criteria import system_family

    assert system_family("bpy") == "Au-(4,4'-bipyridine)-Au"
    assert system_family("BPY") == "Au-(4,4'-bipyridine)-Au"


def test_bpy_family_marks_bdt_bands_not_applicable(workdir):
    """BPY 是吡啶 N 锚定、非硫醇盐体系，BDT 标定的两条判据必须记"不适用"。"""
    from smj.junction.validate import validate_junction_run

    project = workdir / "bpy-family"
    stage = project / "aitranss_te"
    stage.mkdir(parents=True, exist_ok=True)
    (stage / "TE.dat").write_text(TE, encoding="utf-8")

    verdict = validate_junction_run(project, system="BPY")
    names = {c["name"] for c in verdict["checks"] if c.get("applicable") is False}
    assert {"g_in_theory_band", "seebeck_sign"} <= names
    assert not {"g_in_theory_band", "seebeck_sign"} & set(verdict.get("failed_criteria") or [])


def test_sme_family_registered():
    from smj.junction.criteria import system_family

    assert system_family("sme") == "Au-(1,4-bis(methylthio)benzene)-Au"
    assert system_family("SME") == "Au-(1,4-bis(methylthio)benzene)-Au"


def test_sme_family_marks_bdt_bands_not_applicable(workdir):
    """SMe 是硫醚配位接触（S 上无 H，不脱质子），BDT 标定的两条判据必须记"不适用"。"""
    from smj.junction.validate import validate_junction_run

    project = workdir / "sme-family"
    stage = project / "aitranss_te"
    stage.mkdir(parents=True, exist_ok=True)
    (stage / "TE.dat").write_text(TE, encoding="utf-8")

    verdict = validate_junction_run(project, system="SME")
    names = {c["name"] for c in verdict["checks"] if c.get("applicable") is False}
    assert {"g_in_theory_band", "seebeck_sign"} <= names
    assert not {"g_in_theory_band", "seebeck_sign"} & set(verdict.get("failed_criteria") or [])


def test_alkane_family_registered():
    """HDT/ODT/DDT（C6/C8/C10 烷烃二硫醇）登记为同一判据标定族。"""
    from smj.junction.criteria import system_family

    for key in ("hdt", "HDT", "odt", "ODT", "ddt", "DDT", "c6", "C8", "C10"):
        assert system_family(key) == "Au-alkane-dithiol-Au"


def test_alkane_family_marks_bdt_bands_not_applicable(workdir):
    """饱和烷链走 σ 隧穿通道，BDT 标定的两条判据必须记"不适用"而非"否证"。"""
    from smj.junction.validate import validate_junction_run

    project = workdir / "alkane-family"
    stage = project / "aitranss_te"
    stage.mkdir(parents=True, exist_ok=True)
    (stage / "TE.dat").write_text(TE, encoding="utf-8")

    for system in ("HDT", "ODT", "DDT"):
        verdict = validate_junction_run(project, system=system)
        names = {c["name"] for c in verdict["checks"] if c.get("applicable") is False}
        assert {"g_in_theory_band", "seebeck_sign"} <= names
        assert not {"t_ef_magnitude", "g_in_theory_band",
                    "seebeck_sign"} & set(verdict.get("failed_criteria") or [])
