"""species 目录审计（smj.species）：档位、基组规模、指纹与标签核对。"""

from __future__ import annotations

import json
from pathlib import Path

from smj import species

C_SPECIES = """\
################################################################################
#
#  Suggested "light" defaults for C atom (to be pasted into control.in file)
#
################################################################################
     ionic 4 d auto
     hydro 2 p 1.7
     hydro 3 s auto
#     hydro 4 f 9.8
#
"""

S_SPECIES = """\
################################################################################
#
#  Suggested "tight" defaults for S atom (to be pasted into control.in file)
#
################################################################################
     ionic 3 d auto
     hydro 2 p 1.8
     hydro 4 f 7
     hydro 3 s auto
     hydro 3 d 3.9
#
"""


def _species_dir(root: Path) -> Path:
    directory = root / "species"
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "06_C_default").write_text(C_SPECIES, encoding="utf-8")
    (directory / "16_S_default").write_text(S_SPECIES, encoding="utf-8")
    return directory


def test_parse_species_reports_tier_and_basis_size():
    parsed = species.parse_species(C_SPECIES)
    assert parsed["tier"] == "light"
    assert parsed["basis_functions"] == 3  # 注释掉的 f 行不算
    assert parsed["max_angular"] == "d"  # ionic 4 d 是最高角动量
    assert species.parse_species(S_SPECIES)["max_angular"] == "f"
    assert species.parse_species(S_SPECIES)["basis_functions"] == 5


def test_species_info_records_sha256(workdir):
    directory = _species_dir(workdir)
    info = species.species_info(directory, "S")
    assert info.element == "S"
    assert len(info.sha256) == 64
    assert info.tier == "tight"
    assert info.path.endswith("16_S_default")


def test_audit_species_dir_scans_defaults(workdir):
    directory = _species_dir(workdir)
    infos = species.audit_species_dir(directory)
    assert [info.element for info in infos] == ["C", "S"]


def test_summary_flags_label_mismatch(workdir):
    directory = _species_dir(workdir)
    payload = species.summary(directory, elements=["C"], label="tight")
    assert payload["actual_tier"] == "light"
    assert payload["declared_label"] == "tight"
    assert any("只是配置标签" in warning for warning in payload["warnings"])
    assert payload["elements"][0]["basis_functions"] == 3


def test_summary_accepts_matching_label(workdir):
    directory = _species_dir(workdir)
    assert species.summary(directory, elements=["C"], label="light")["warnings"] == []


def test_summary_detects_mixed_tiers(workdir):
    directory = _species_dir(workdir)
    payload = species.summary(directory, elements=["C", "S"])
    assert payload["actual_tier"] == "mixed"
    assert any("混用了多个档位" in warning for warning in payload["warnings"])


def test_elements_of_keeps_geometry_order(workdir):
    atoms = [
        type("A", (), {"element": "S"})(),
        type("A", (), {"element": "C"})(),
        type("A", (), {"element": "S"})(),
    ]
    assert species.elements_of(atoms) == ["S", "C"]


def test_cli_species_prints_table_and_json(workdir, capsys):
    from smj.cli import main

    directory = _species_dir(workdir)
    assert main(["species", "--dir", str(directory), "--elements", "C",
                 "--label", "tight"]) == 0
    out = capsys.readouterr().out
    assert "实际档位: light" in out
    assert "只是配置标签" in out

    out_file = workdir / "species.json"
    assert main(["species", "--dir", str(directory), "--elements", "S",
                 "--out", str(out_file), "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["elements"][0]["element"] == "S"
    assert json.loads(out_file.read_text(encoding="utf-8"))["actual_tier"] == "tight"
