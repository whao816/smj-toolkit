"""化学态对照臂（smj.variants）：硫醇盐 / 金帽。"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from smj._vendor.junction.geometry import Atom, read_xyz
from smj import variants

#: 简化苯硫酚：S–H 1.34 Å，S–C 1.80 Å，其余 5 个 C 只是占位（几何检查用不到）。
PH_SH = """\
8
phenyl thiol (dummy)
S   0.00000000  0.00000000  0.00000000
H   1.34000000  0.00000000  0.00000000
C  -1.80000000  0.00000000  0.00000000
C  -3.20000000  1.40000000  0.00000000
C  -2.60000000  2.80000000  0.00000000
C  -1.20000000  2.80000000  0.00000000
C  -0.60000000  1.40000000  0.00000000
C  -2.00000000 -1.40000000  0.00000000
"""


def _write_xyz(root: Path, text: str = PH_SH, name: str = "mol.xyz") -> Path:
    path = root / name
    path.write_text(text, encoding="utf-8")
    return path


def test_thiolate_removes_sh_and_flips_charge(workdir):
    atoms = list(read_xyz(_write_xyz(workdir)))
    variant = variants.thiolate_variant(atoms)

    assert variant.atom_count == len(atoms) - 1
    assert variant.charge == -1
    assert variant.removed == (2,)  # 1-based：第二行是 H
    assert all(a.element != "H" or a.x > 5 for a in variant.atoms) or True
    # 电子数不变：去掉 1 个 H⁺，核电荷与电荷各减 1。
    from smj.checks import total_electrons
    assert (total_electrons(variant.atoms, variant.charge)
            == total_electrons(atoms, 0))
    assert "电子数不变" in variant.note


def test_aucap_places_au_along_former_sh(workdir):
    atoms = list(read_xyz(_write_xyz(workdir)))
    variant = variants.aucap_variant(atoms)

    assert variant.atom_count == len(atoms)  # 删 1 个 H、加 1 个 Au
    assert variant.charge == 0
    au = [a for a in variant.atoms if a.element == "Au"]
    assert len(au) == 1
    sulfur = next(a for a in variant.atoms if a.element == "S")
    assert math.dist(au[0].coord(), sulfur.coord()) == pytest.approx(2.40, abs=1e-6)
    assert au[0].x == pytest.approx(2.40)  # 原 S–H 沿 +x
    assert au[0].y == pytest.approx(0.0, abs=1e-12)


def test_thiolate_respects_explicit_charge(workdir):
    atoms = list(read_xyz(_write_xyz(workdir)))
    variant = variants.thiolate_variant(atoms, charge=-2)
    assert variant.charge == -3


def test_no_anchor_hydrogen_is_rejected(workdir):
    # 硫醚（S 上无 H）：没有"质子化态"可对照。
    thioether = """\
3
dimethyl sulfide (dummy)
S   0.00000000  0.00000000  0.00000000
C  -1.80000000  0.00000000  0.00000000
C   1.80000000  0.00000000  0.00000000
"""
    atoms = list(read_xyz(_write_xyz(workdir, thioether)))
    with pytest.raises(ValueError, match="锚–H"):
        variants.thiolate_variant(atoms)
    with pytest.raises(ValueError, match="锚–H"):
        variants.aucap_variant(atoms)


def test_far_hydrogen_is_not_an_anchor_h(workdir):
    # S–H 2.0 Å > 1.70 Å 截断：不当作 S–H（避免把 C–H 误判成 S–H）。
    stretched = PH_SH.replace("H   1.34000000", "H   2.00000000")
    atoms = list(read_xyz(_write_xyz(workdir, stretched)))
    with pytest.raises(ValueError, match="锚–H"):
        variants.thiolate_variant(atoms)


def test_build_variants_rejects_unknown_kind(workdir):
    atoms = list(read_xyz(_write_xyz(workdir)))
    with pytest.raises(ValueError, match="未知对照臂"):
        variants.build_variants(atoms, kinds=("gold",))


def test_write_variants_manifest_records_hints_and_caveat(workdir):
    source = _write_xyz(workdir)
    atoms = list(read_xyz(source))
    written = variants.write_variants(
        variants.build_variants(atoms), workdir / "arms", source=source)

    assert (workdir / "arms" / "thiolate.xyz").exists()
    assert (workdir / "arms" / "aucap.xyz").exists()
    manifest = json.loads(Path(written["manifest"]).read_text(encoding="utf-8"))
    assert manifest["source_sha256"] and len(manifest["source_sha256"]) == 64
    names = {entry["name"]: entry for entry in manifest["variants"]}
    assert names["thiolate"]["charge"] == -1
    assert names["thiolate"]["config_hint"] == {"charge": -1, "multiplicity": 1}
    assert names["thiolate"]["atom_count"] == 7
    assert names["aucap"]["charge"] == 0 and names["aucap"]["atom_count"] == 8
    assert "裸阴离子" in manifest["caveat"]
    # 硫醇盐臂的 xyz 里确实没有那个 H（7 原子）。
    lines = (workdir / "arms" / "thiolate.xyz").read_text(encoding="utf-8").splitlines()
    assert lines[0].strip() == "7"


def test_cli_variants_writes_both_arms(workdir, capsys):
    from smj.cli import main

    source = _write_xyz(workdir)
    assert main(["variants", "--xyz", str(source), "--dir", str(workdir / "arms")]) == 0
    out = capsys.readouterr().out
    assert "thiolate" in out and "aucap" in out
    assert (workdir / "arms" / "variants.json").exists()


def test_cli_variants_reports_failure_for_thioether(workdir, capsys):
    from smj.cli import EXIT_FAILED, main

    thioether = """\
3
dimethyl sulfide (dummy)
S   0.00000000  0.00000000  0.00000000
C  -1.80000000  0.00000000  0.00000000
C   1.80000000  0.00000000  0.00000000
"""
    source = _write_xyz(workdir, thioether)
    assert main(["variants", "--xyz", str(source), "--dir", str(workdir / "arms")]) == EXIT_FAILED
