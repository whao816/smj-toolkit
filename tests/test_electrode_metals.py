"""0.2.x 回归：electrode_build 的「构建方式 + 电极类型」参数化。

背景：阶段 4 从「本地构建 Au(111) 电极」推广为「构建电极」——
是本地三角岛（``source: builtin``）还是 FHI-aims ``electrodes.library``
预定义团簇（``source: library``），以及电极金属（``electrode.element`` /
``electrode.a_lattice``，fcc(111) 单晶面）都成为本阶段参数。

本组测试锁住：
* ``electrode.element`` 的默认值、大小写归一化、未登记元素必须显式给晶格常数；
* tip 元素跟随 ``electrode.element``（``make_tip_relax_input(tip_element=...)``）；
* ``source: builtin`` 用非 Au 元素装配（ElectrodeFCC 晶格常数推导）；
* ``source: library`` 允许任意单一元素团簇，但必须与 ``electrode.element`` 一致
  （不一致 fail-closed）；
* 默认 Au / builtin 的 metadata 逐字节不变（不写 ``electrode_element``）。
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from smj.cli import main
from smj.config import load_config
from smj.errors import ConfigError
from smj.junction.electrode import append_au60, FCC_LATTICE_A
from smj.junction.geometry import read_aims_geometry, write_aims_geometry
from smj.junction.inputs import make_tip_relax_input

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

# 与 tests/test_library_au_relax.py 同一份 18 原子 Au(111) 团簇。
AU_LIB_TEXT = """#########################################################
#
#  number of au atoms = 18
#
#  execute "tcontrol.aims.x" script with option:
#  -nlayers 2
#
#########################################################
  atom     0.00000000     0.00000000     0.00000000    Au
  atom    -1.44186146    -0.83245910     2.35454990    Au
  atom     1.44186146    -0.83245910     2.35454990    Au
  atom     0.00000000     1.66491820     2.35454990    Au
  atom    -2.88372292    -1.66491820     4.70909981    Au
  atom     0.00000000    -1.66491820     4.70909981    Au
  atom     2.88372292    -1.66491820     4.70909981    Au
  atom    -1.44186146     0.83245910     4.70909981    Au
  atom     1.44186146     0.83245910     4.70909981    Au
  atom     0.00000000     3.32983641     4.70909981    Au
  atom    -1.44186146    -2.49737730     7.06364971    Au
  atom     1.44186146    -2.49737730     7.06364971    Au
  atom    -2.88372292     0.00000000     7.06364971    Au
  atom     0.00000000     0.00000000     7.06364971    Au
  atom     2.88372292     0.00000000     7.06364971    Au
  atom    -1.44186146     2.49737730     7.06364971    Au
  atom     1.44186146     2.49737730     7.06364971    Au
  atom     2.88372292     3.32983641     4.70909981    Au
#end
"""


def _tips_geometry(element: str = "Au") -> str:
    """分子 + 两枚 tip（M–S = 2.35 Å，沿 x 轴），模拟 aims_tip_relax 的产物。"""
    lines = ["# BDT + two %s tips" % element]
    for line in BDT_XYZ.splitlines()[2:]:
        parts = line.split()
        if len(parts) != 4:
            continue
        lines.append("  atom  %s %s %s %s" % (parts[1], parts[2], parts[3], parts[0]))
    lines.append("  atom  5.510000 0.000000 0.000000 %s" % element)
    lines.append("  atom  -5.510000 0.000000 0.000000 %s" % element)
    return "\n".join(lines) + "\n"


def _scaled_library(element: str, lattice: float, text: str = AU_LIB_TEXT) -> str:
    """把 Au 团簇按晶格常数缩放成另一元素（fcc(111) 层内平移等比）。"""
    factor = lattice / FCC_LATTICE_A["Au"]
    out = []
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 5 and parts[0] == "atom":
            x, y, z = (float(parts[i]) * factor for i in (1, 2, 3))
            out.append("  atom %15.8f %15.8f %15.8f    %s" % (x, y, z, element))
        else:
            out.append(line)
    return "\n".join(out) + "\n"


def _project(workdir: Path, *, source: str = "builtin", element: str | None = None,
             a_lattice: float | None = None, au_per_electrode: int = 4,
             library_cluster: str | None = None,
             library_element: str = "Au", tips_element: str = "Au") -> tuple[Path, Path]:
    """最小工程：smj.yaml + 分子 + species + 上游 tips 几何（+ 可选库文件）。"""
    project = (workdir / "proj").resolve()
    project.mkdir(parents=True)
    (project / "molecule.xyz").write_text(BDT_XYZ, encoding="utf-8")
    species = project / "species"
    species.mkdir()
    for el in ("Au", "Ag", "C", "H", "S"):
        (species / ("%s_default" % el)).write_text(
            "  species %s\n    # dummy\n" % el, encoding="utf-8")

    run_dir = project / "t"
    tips = run_dir / "aims_tip_relax"
    tips.mkdir(parents=True)
    (tips / "geometry.in.next_step").write_text(
        _tips_geometry(tips_element), encoding="utf-8")

    lines = [
        "run:", "  name: t", "  mode: dry-run", "  dir: %s" % project.as_posix(),
        "molecule:", "  xyz: %s" % (project / "molecule.xyz").as_posix(),
        "electrode:", "  source: %s" % source,
        "  au_per_electrode: %d" % au_per_electrode,
    ]
    if element is not None:
        lines.append("  element: %s" % element)
    if a_lattice is not None:
        lines.append("  a_lattice: %s" % a_lattice)
    if library_cluster is not None:
        lib = project / "lib"
        target = lib / library_cluster
        target.parent.mkdir(parents=True)
        if library_element == "Au":
            target.write_text(AU_LIB_TEXT, encoding="utf-8")
        else:
            target.write_text(_scaled_library(library_element,
                                              FCC_LATTICE_A[library_element]),
                              encoding="utf-8")
        lines.append("  library_dir: %s" % lib.as_posix())
        lines.append("  library_cluster: %s" % library_cluster)
    lines += ["cluster:", "  species_dir: %s" % species.as_posix()]
    cfg = project / "smj.yaml"
    cfg.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return cfg, run_dir


def _run(*args: str) -> int:
    return main(list(args))


def _meta(run_dir: Path) -> dict:
    return json.loads((run_dir / "electrode_build" / "electrode_metadata.json").read_text(
        encoding="utf-8"))


# --------------------------------------------------------------------------- #
# 1) 配置：元素 + 晶格常数
# --------------------------------------------------------------------------- #
def test_electrode_metal_defaults_to_au(workdir):
    cfg, _ = _project(workdir)
    config = load_config(cfg)
    assert config.electrode.element == "Au"
    assert config.electrode.a_lattice is None


def test_electrode_metal_normalizes_case(workdir):
    cfg, _ = _project(workdir, element="ag")
    assert load_config(cfg).electrode.element == "Ag"


def test_electrode_metal_rejects_bad_symbol(workdir):
    cfg, _ = _project(workdir, element="Auu")
    with pytest.raises(ConfigError) as excinfo:
        load_config(cfg)
    assert "electrode.element" in str(excinfo.value)


def test_unregistered_metal_requires_explicit_lattice(workdir):
    cfg, _ = _project(workdir, element="Xx")
    with pytest.raises(ConfigError) as excinfo:
        load_config(cfg)
    assert "a_lattice" in str(excinfo.value)


def test_unregistered_metal_with_lattice_ok(workdir):
    cfg, _ = _project(workdir, element="Xx", a_lattice=4.2)
    config = load_config(cfg)
    assert config.electrode.element == "Xx"
    assert config.electrode.a_lattice == pytest.approx(4.2)


def test_lattice_out_of_range(workdir):
    cfg, _ = _project(workdir, element="Ag", a_lattice=9.9)
    with pytest.raises(ConfigError) as excinfo:
        load_config(cfg)
    assert "2.0–6.0" in str(excinfo.value)


# --------------------------------------------------------------------------- #
# 2) tip 元素跟随 electrode.element
# --------------------------------------------------------------------------- #
def _species_dir(root: Path) -> Path:
    species = root / "species"
    species.mkdir()
    for el in ("Au", "Ag", "C", "H", "S"):
        (species / ("%s_default" % el)).write_text(
            "  species %s\n    # dummy\n" % el, encoding="utf-8")
    return species


def test_tip_element_follows_electrode_metal(workdir):
    xyz = workdir / "molecule.xyz"
    xyz.write_text(BDT_XYZ, encoding="utf-8")
    species = _species_dir(workdir)
    outdir = workdir / "tips"
    info = make_tip_relax_input(
        xyz, outdir, left_anchor=7, right_anchor=8, au_tip_distance=2.35,
        species_dir=species, anchor_au_s_target=None, pin_contact=False,
        env={"aims_bin": "/bin/true", "partition": "main", "ntasks": 1},
        tip_element="Ag")
    atoms = read_aims_geometry(outdir / "geometry.in")
    tips = [atom for atom in atoms if atom.element == "Ag"]
    assert len(tips) == 2
    assert not any(atom.element == "Au" for atom in atoms)
    assert info["summary"]["tip_element"] == "Ag"
    assert info["summary"]["au_tip_count"] == 2


# --------------------------------------------------------------------------- #
# 3) builtin：非 Au 三角岛
# --------------------------------------------------------------------------- #
def test_builtin_ag_electrode_build(workdir, capsys):
    cfg, run_dir = _project(workdir, source="builtin", element="Ag",
                            au_per_electrode=4, tips_element="Ag")
    assert _run("stage", "electrode_build", "--prepare", "-c", str(cfg)) == 0, \
        capsys.readouterr()

    meta = _meta(run_dir)
    assert meta["electrode_element"] == "Ag"
    assert meta["electrode_a_lattice"] == pytest.approx(4.0853)
    # 14 分子 + 2 tip + 2 侧各 3 个补全原子 = 22
    assert meta["final_atom_count"] == 22
    atoms = read_aims_geometry(run_dir / "electrode_build" / "geometry.in")
    assert sum(1 for atom in atoms if atom.element == "Ag") == 8
    assert not any(atom.element == "Au" for atom in atoms)


def test_builtin_custom_lattice_used(workdir, capsys):
    cfg, run_dir = _project(workdir, source="builtin", element="Xx",
                            a_lattice=4.0, au_per_electrode=4,
                            tips_element="Xx")
    assert _run("stage", "electrode_build", "--prepare", "-c", str(cfg)) == 0, \
        capsys.readouterr()
    meta = _meta(run_dir)
    assert meta["electrode_element"] == "Xx"
    assert meta["electrode_a_lattice"] == pytest.approx(4.0)


def test_au_builtin_metadata_stays_legacy(workdir):
    """默认 Au / builtin：metadata 不新增 electrode_element 键（与 0.2.0 一致）。"""
    from smj.junction.geometry import Atom

    atoms = [Atom("C", 0, 0, 0), Atom("S", 0, 1.8, 0), Atom("S", 0, -1.8, 0),
             Atom("Au", 0, 3.8, 0), Atom("Au", 0, -3.8, 0)]
    src = workdir / "tips.in"
    write_aims_geometry(src, atoms)
    meta = append_au60(src, workdir / "out.in", workdir / "meta.json",
                       au_per_electrode=4)
    assert "electrode_element" not in meta
    assert "electrode_a_lattice" not in meta


# --------------------------------------------------------------------------- #
# 4) library：任意单一元素，但必须与配置一致
# --------------------------------------------------------------------------- #
def test_library_ag_cluster_build(workdir, capsys):
    cfg, run_dir = _project(
        workdir, source="library", element="Ag",
        library_cluster="fcc111.ag.clusters.18_atoms/natoms18.x1.aims",
        library_element="Ag", tips_element="Ag")
    assert _run("stage", "electrode_build", "--prepare", "-c", str(cfg)) == 0, \
        capsys.readouterr()

    meta = _meta(run_dir)
    assert meta["electrode_source"] == "library"
    assert meta["electrode_element"] == "Ag"
    # 18 原子团簇去掉 on-top 接触原子 -> 每侧 17 个；14 + 2 + 34 = 50
    assert meta["final_atom_count"] == 50
    atoms = read_aims_geometry(run_dir / "electrode_build" / "geometry.in")
    assert sum(1 for atom in atoms if atom.element == "Ag") == 36
    assert not any(atom.element == "Au" for atom in atoms)


def test_library_element_mismatch_fails_closed(workdir, capsys):
    cfg, _ = _project(
        workdir, source="library", element="Ag",
        library_cluster="fcc111.au.clusters.18_atoms/natoms18.x1.aims",
        library_element="Au", tips_element="Ag")
    rc = _run("stage", "electrode_build", "--prepare", "-c", str(cfg))
    assert rc == 1
    assert "电极元素不一致" in capsys.readouterr().err
