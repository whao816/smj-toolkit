"""0.2.x 回归：分子里没有 S 时的锚定/去质子路径。

背景：锚定原子不再限定硫醇盐的 S —— 任意元素符号（N/Se/…）或 1-based 序号都能
指定 ``electrode.left_anchor/right_anchor``。去质子口径改为三态：
``strip_thiol_h: null``（缺省）= 只有锚定原子是 S 才自动移除 S–H；``true`` = 强制去，
``false`` = 从来不去。非 S 锚默认保留 H，且校验只报告接触距离、不套 2.40–2.49 Å
硫醇盐有源带（该带只对 S/Au 标定）。

本组测试锁住：
* 配置三态解析与模板缺省 null；
* 无 S 分子按元素符号/序号解析锚，``make_tip_relax_input`` 不报"需 2 个 S"；
* 无 S 结的 ``junction_geom_anchor_aus`` 是 info（不 warn），value 里写明锚定元素；
* S 分子在 null 口径下仍然自动去质子（旧行为不回归）；
* 强制去质子遇到一个锚原子带 2 个 H（-NH2）时 fail-closed。
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from smj.junction.electrode import append_au60
from smj.junction.geometry import read_aims_geometry
from smj.junction.inputs import make_tip_relax_input
from smj.junction.validate import validate_junction_geometry
from smj.config import default_config_text, load_config
from smj.stages._helpers import resolve_anchor_indices

#: 吡嗪（1,4-二氮杂苯）：两个 N 在 1/4 位，沿 x 轴相对；分子里没有 S。
PYR_XYZ = """6
pyrazine (N at 1/4, no sulfur)
N -1.400  0.000  0.000
C -0.700  1.200  0.000
C  0.700  1.200  0.000
N  1.400  0.000  0.000
C  0.700 -1.200  0.000
C -0.700 -1.200  0.000
"""

#: 两个 -NH2（每个 N 带 2 个 H）：强制去质子应当 fail-closed，而不是猜着去一个。
NH2_XYZ = """6
diamino
N -3.000  0.000  0.000
H -2.220  0.900  0.000
H -2.220 -0.900  0.000
N  3.000  0.000  0.000
H  2.220  0.900  0.000
H  2.220 -0.900  0.000
"""

BDT_XYZ = """14
benzene-1,4-dithiol (BDT)
C  1.390  0.000  0.000
C  0.695  1.204  0.000
C -0.695  1.204  0.000
C -1.390  0.000  0.000
C -0.695 -1.204  0.000
C  0.695 -1.204  0.000
S  3.160  0.000  0.000
S -3.160  0.000  0.000
H  1.235  2.140  0.000
H -1.235  2.140  0.000
H -1.235 -2.140  0.000
H  1.235 -2.140  0.000
H  4.490  0.000  0.000
H -4.490  0.000  0.000
"""


def _species(root: Path) -> Path:
    species = root / "species"
    species.mkdir()
    for el in ("Au", "C", "H", "N", "S"):
        (species / ("%s_default" % el)).write_text(
            "  species %s\n    # dummy\n" % el, encoding="utf-8")
    return species


def _env() -> dict:
    return {"aims_bin": "/bin/true", "partition": "main", "ntasks": 1}


def _config(workdir: Path, strip_line: str | None) -> Path:
    project = (workdir / "proj").resolve()
    project.mkdir(parents=True)
    xyz = project / "molecule.xyz"
    xyz.write_text(PYR_XYZ, encoding="utf-8")
    species = _species(project)
    lines = ["run:", "  name: t", "  mode: dry-run",
             "  dir: %s" % project.as_posix(),
             "molecule:", "  xyz: %s" % xyz.as_posix(),
             "electrode:", "  left_anchor: N", "  right_anchor: N"]
    if strip_line is not None:
        lines.append("  strip_thiol_h: %s" % strip_line)
    lines += ["cluster:", "  species_dir: %s" % species.as_posix()]
    cfg = project / "smj.yaml"
    cfg.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return cfg


def _tips(workdir: Path, xyz_text: str, *, left: int, right: int,
          strip) -> dict:
    xyz = workdir / "molecule.xyz"
    xyz.write_text(xyz_text, encoding="utf-8")
    info = make_tip_relax_input(
        xyz, workdir / "tips", left_anchor=left, right_anchor=right,
        au_tip_distance=2.35, species_dir=_species(workdir),
        strip_thiol_h=strip, anchor_au_s_target=None, pin_contact=False, env=_env())
    return info


# --------------------------------------------------------------------------- #
# 1) 配置：strip_thiol_h 三态
# --------------------------------------------------------------------------- #
def test_strip_thiol_h_defaults_to_auto(workdir):
    assert load_config(_config(workdir, None)).electrode.strip_thiol_h is None


def test_strip_thiol_h_parses_tristate(workdir):
    assert load_config(_config(workdir / "on", "true")).electrode.strip_thiol_h is True
    assert load_config(_config(workdir / "off", "false")).electrode.strip_thiol_h is False
    assert load_config(_config(workdir / "auto", "null")).electrode.strip_thiol_h is None


def test_default_template_strip_is_null():
    text = default_config_text(name="t", xyz="molecule.xyz")
    assert "strip_thiol_h: null" in text


# --------------------------------------------------------------------------- #
# 2) 无 S 分子：锚解析与加 tip
# --------------------------------------------------------------------------- #
def test_anchor_resolves_by_element_without_s(workdir):
    xyz = workdir / "pyr.xyz"
    xyz.write_text(PYR_XYZ, encoding="utf-8")
    assert resolve_anchor_indices(xyz, "N", "N") == (1, 4)


def test_anchor_resolves_by_index_without_s(workdir):
    xyz = workdir / "pyr.xyz"
    xyz.write_text(PYR_XYZ, encoding="utf-8")
    assert resolve_anchor_indices(xyz, "1", "4") == (1, 4)


def test_no_s_tip_prepare_keeps_h_by_default(workdir):
    info = _tips(workdir, PYR_XYZ, left=1, right=4, strip=None)
    summary = info["summary"]
    assert summary["anchor_element"] == "N"
    assert summary["strip_thiol_h"] is None       # null = 非 S 锚保留 H，不做任何移除
    atoms = read_aims_geometry(workdir / "tips" / "geometry.in")
    assert sum(1 for a in atoms if a.element == "Au") == 2
    assert not any(a.element == "S" for a in atoms)


# --------------------------------------------------------------------------- #
# 3) 无 S 结：校验给 info，不套硫醇盐带
# --------------------------------------------------------------------------- #
def test_no_s_junction_anchor_check_is_info(workdir):
    info = _tips(workdir, PYR_XYZ, left=1, right=4, strip=None)
    src = Path(info["geometry"])
    expanded = workdir / "expanded.in"
    meta_path = workdir / "electrode_metadata.json"
    meta = append_au60(src, expanded, meta_path, au_per_electrode=4)
    report = validate_junction_geometry(workdir, metadata=meta, geometry=expanded)

    check = next(c for c in report["checks"] if c["name"] == "junction_geom_anchor_aus")
    assert check["severity"] == "info"
    assert check["value"]["anchor_element"] == "N"
    assert check["value"]["chemistry"] == "non_sulfur"
    assert check["value"]["distance_A"] == pytest.approx(2.35, abs=1e-3)
    assert "不适用" in check["detail"]
    assert "缺 S" not in check["detail"]


# --------------------------------------------------------------------------- #
# 4) S 分子：null 口径仍然自动去质子（旧行为不回归）
# --------------------------------------------------------------------------- #
def test_sulfur_anchor_still_auto_strips(workdir):
    info = _tips(workdir, BDT_XYZ, left=7, right=8, strip=None)
    summary = info["summary"]
    assert summary["anchor_element"] == "S"
    assert summary["strip_thiol_h"]["removed"] == 2
    assert summary["strip_thiol_h"]["route"] == "thiolate"
    assert summary["strip_thiol_h"]["sulfur_indices"] == [7, 8]


# --------------------------------------------------------------------------- #
# 5) 强制去质子：-NH2 带 2 个 H 时 fail-closed
# --------------------------------------------------------------------------- #
def test_force_strip_on_nh2_fails_closed(workdir):
    with pytest.raises(ValueError) as excinfo:
        _tips(workdir, NH2_XYZ, left=1, right=4, strip=True)
    message = str(excinfo.value)
    assert "带 2 个 H" in message
    assert "strip_thiol_h: false" in message


def test_cli_tip_prepare_without_s(workdir, capsys):
    from smj.cli import main

    cfg = _config(workdir, "null")
    project = cfg.parent
    xyz = project / "molecule.xyz"
    rc = main(["stage", "aims_tip_relax", "--prepare", "--geometry", str(xyz),
               "-c", str(cfg)])
    assert rc == 0, capsys.readouterr()
    summary = json.loads((project / "t" / "aims_tip_relax" / "input_summary.json")
                         .read_text(encoding="utf-8"))
    assert summary["anchor_element"] == "N"
    assert summary["tip_element"] == "Au"
    assert summary["strip_thiol_h"] is None
    assert summary["au_tip_count"] == 2
