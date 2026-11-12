"""0.1.2 回归：库电极（electrodes.library）+ 背层固定弛豫阶段 ``au_relax``。

背景（2026-11-11，FHI-aims 手册 §5.3.1）：官方流程是用 ``electrodes.library`` 里
预定义的 Au(111) 团簇装配电极，装配后**固定两侧背层做一次几何弛豫**，再进 SCF。
本组测试锁住：

* 库文件解析的几条硬规则（接触层恰 1 原子、背层 >=3、层内共面、单一元素、
  头部 ``-nlayers``）；
* ``au60_build`` 用库团簇装配后 metadata 带 ``electrode_source`` /
  ``library_cluster`` / 背层索引；
* ``au_relax`` prepare 把**整层背层**写成 ``constrain_relaxation .true.``；
* ``au_relax`` collect 的三条 fail-closed：没有 "Have a nice day"、缺
  ``geometry.in.next_step``、被冻结原子发生位移、表面三元组与 au60 不一致；
* ``electrode.relax=false`` 是本地直通（local_only），collect 不要求 aims.out；
* ``electrode.nlayers: null`` 时从库文件建议值取 nlayers。
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from smj.cli import main
from smj.electrode_library import parse_cluster_text
from smj.stages.aitranss_te import AitranssTeStage

# --------------------------------------------------------------------------- #
# 固件：FHI-aims 库里的 18 原子 Au(111) 团簇（头部建议 nlayers=2）
# 层结构 1 / 3 / 7 / 7，层距 2.35455 Å；与 vendored electrode.py 同源晶格。
# --------------------------------------------------------------------------- #
DZ = 2.35454990
LIB_TEXT = """#########################################################
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


def _tips_geometry() -> str:
    """分子 + 两枚 Au tip（S–Au = 2.35 Å，沿 x 轴），模拟 aims_tip_relax 的产物。"""
    lines = ["# BDT + two Au tips"]
    for line in BDT_XYZ.splitlines()[2:]:
        parts = line.split()
        if len(parts) != 4:
            continue
        lines.append("  atom  %s %s %s %s" % (parts[1], parts[2], parts[3], parts[0]))
    lines.append("  atom  5.510000 0.000000 0.000000 Au")
    lines.append("  atom  -5.510000 0.000000 0.000000 Au")
    return "\n".join(lines) + "\n"


def _write_project(workdir: Path, *, relax: bool = True, relax_freeze: str = "backplane",
                   cluster_relpath: str = "fcc111.au.clusters.18_atoms/natoms18.x1.aims",
                   source: str = "library") -> tuple[Path, Path]:
    """写一个最小工程：smj.yaml + 库文件 + species + 上游 tips 几何。"""
    project = (workdir / "proj").resolve()
    (project / "lib" / Path(cluster_relpath).parent).mkdir(parents=True)
    (project / "lib" / cluster_relpath).write_text(LIB_TEXT, encoding="utf-8")
    species = project / "species"
    species.mkdir()
    for element in ("Au", "C", "H", "S"):
        (species / ("%s_default" % element)).write_text(
            "  species %s\n    # dummy\n" % element, encoding="utf-8")
    (project / "molecule.xyz").write_text(BDT_XYZ, encoding="utf-8")

    run_dir = project / "t"
    tips = run_dir / "aims_tip_relax"
    tips.mkdir(parents=True)
    (tips / "geometry.in.next_step").write_text(_tips_geometry(), encoding="utf-8")

    cfg = project / "smj.yaml"
    cfg.write_text(
        "run:\n"
        "  name: t\n"
        "  mode: dry-run\n"
        "  dir: %s\n"
        "molecule:\n"
        "  xyz: %s\n"
        "electrode:\n"
        "  source: %s\n"
        "  library_dir: %s\n"
        "  library_cluster: %s\n"
        "  relax: %s\n"
        "  relax_freeze: %s\n"
        "cluster:\n"
        "  species_dir: %s\n"
        % (project.as_posix(), (project / "molecule.xyz").as_posix(), source,
           (project / "lib").as_posix(), cluster_relpath,
           "true" if relax else "false", relax_freeze, species.as_posix()),
        encoding="utf-8",
    )
    return cfg, run_dir


def _state(run_dir: Path) -> dict:
    return json.loads((run_dir / "state.json").read_text(encoding="utf-8"))


def _run(*args: str) -> int:
    return main(list(args))


# --------------------------------------------------------------------------- #
# 1) 库文件解析
# --------------------------------------------------------------------------- #
def test_parse_library_cluster_basic():
    cluster = parse_cluster_text(LIB_TEXT, source="mem", relpath="x/natoms18.x1.aims")
    assert cluster.natoms == 18
    assert cluster.nlayers_suggested == 2
    assert cluster.element == "Au"
    assert cluster.layers() == [(0.0, 1), (2.3545, 3), (4.7091, 7), (7.0636, 7)]
    assert cluster.translated is False
    assert cluster.sha256 == hashlib.sha256(LIB_TEXT.encode("utf-8")).hexdigest()
    assert cluster.provenance()["nlayers_suggested"] == 2


def test_parse_translates_non_origin_contact():
    import re

    def shift(match: re.Match) -> str:
        x, y, z = (float(match.group(i)) + d for i, d in ((1, 1.5), (2, 2.0), (3, 3.0)))
        return "atom %.8f %.8f %.8f Au" % (x, y, z)

    shifted = re.sub(r"atom\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+Au", shift, LIB_TEXT)
    cluster = parse_cluster_text(shifted, source="mem", relpath="x.aims")
    assert cluster.translated is True
    first = cluster.atoms[0]
    assert abs(first.x) < 1e-9 and abs(first.y) < 1e-9 and abs(first.z) < 1e-9


def _without_contact(text: str) -> str:
    """删掉 z 全零的接触原子行（保持其他行原样）。"""
    out = []
    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 5 and parts[0] == "atom" and parts[1:] == ["0.00000000", "0.00000000", "0.00000000", "Au"]:
            continue
        out.append(line)
    return "\n".join(out) + "\n"


def _with_z_jitter(text: str, old: str, new: str) -> str:
    assert old in text
    return text.replace(old, new, 1)


@pytest.mark.parametrize("text,message", [
    ("atom 0.0 0.0 0.0 Au\natom 2.88 0.0 0.0 Au\n#end\n", "只有一层"),
    (_without_contact(LIB_TEXT), "接触层"),
    (LIB_TEXT + "atom 9.0 9.0 0.0 Ag\n#end\n", "多种元素"),
    (_with_z_jitter(LIB_TEXT, "     4.70909981    Au", "     4.70959981    Au"), "层不共面"),
    ("# empty\n#end\n", "没有 atom 行"),
])
def test_parse_rejects_bad_library(text, message):
    with pytest.raises(Exception) as excinfo:  # SmjError，这里只锁文案
        parse_cluster_text(text, source="mem", relpath="x.aims")
    assert message in str(excinfo.value)


# --------------------------------------------------------------------------- #
# 2) au60_build（库团簇）+ au_relax prepare
# --------------------------------------------------------------------------- #
def test_au60_build_with_library_cluster(workdir, capsys):
    cfg, run_dir = _write_project(workdir)
    rc = _run("stage", "au60_build", "--prepare", "-c", str(cfg))
    assert rc == 0, capsys.readouterr()

    meta = json.loads((run_dir / "au60_build" / "electrode_metadata.json").read_text(
        encoding="utf-8"))
    assert meta["electrode_source"] == "library"
    assert meta["library_cluster"]["relpath"].endswith("natoms18.x1.aims")
    assert meta["library_cluster"]["natoms"] == 18
    assert meta["library_cluster"]["sha256"] == hashlib.sha256(
        LIB_TEXT.encode("utf-8")).hexdigest()
    # 18 原子团簇去掉 on-top 接触原子 -> 每侧 17 个；分子 14 + 2 tips + 34 = 50
    assert meta["final_atom_count"] == 50
    assert meta["au_per_electrode"] == 18
    assert len(meta["left_backplane_indices"]) == 7
    assert len(meta["right_backplane_indices"]) == 7
    assert meta["surface_fingerprint"]


def test_au_relax_prepare_freezes_whole_backplane(workdir, capsys):
    cfg, run_dir = _write_project(workdir)
    assert _run("stage", "au60_build", "--prepare", "-c", str(cfg)) == 0
    rc = _run("stage", "au_relax", "--prepare", "-c", str(cfg))
    assert rc == 0, capsys.readouterr()

    meta = json.loads((run_dir / "au60_build" / "electrode_metadata.json").read_text(
        encoding="utf-8"))
    frozen = sorted(meta["left_backplane_indices"] + meta["right_backplane_indices"])
    geometry = (run_dir / "au_relax" / "geometry.in").read_text(encoding="utf-8")
    constrained_lines = [line.strip() for line in geometry.splitlines()
                         if "constrain_relaxation" in line]
    assert len(constrained_lines) == len(frozen) == 14
    assert all(".true." in line for line in constrained_lines)

    summary = json.loads((run_dir / "au_relax" / "input_summary.json").read_text(
        encoding="utf-8"))
    assert summary["constrained"] == frozen
    assert summary["constrained_count"] == 14

    control = (run_dir / "au_relax" / "control.in").read_text(encoding="utf-8")
    assert "relax_geometry" in control and "bfgs" in control


# --------------------------------------------------------------------------- #
# 3) au_relax collect 的 fail-closed
# --------------------------------------------------------------------------- #
def test_au_relax_prepare_can_freeze_all_gold(workdir, capsys):
    """relax_freeze=electrodes：冻结全部 Au（含 tip），只让分子骨架动。"""
    cfg, run_dir = _write_project(workdir, relax_freeze="electrodes")
    assert _run("stage", "au60_build", "--prepare", "-c", str(cfg)) == 0
    rc = _run("stage", "au_relax", "--prepare", "-c", str(cfg))
    assert rc == 0, capsys.readouterr()

    geometry = (run_dir / "au_relax" / "geometry.in").read_text(encoding="utf-8")
    lines = [line for line in geometry.splitlines() if line.strip()]
    atom_lines = [line for line in lines if line.split()[0] == "atom"]
    gold = [i + 1 for i, line in enumerate(atom_lines) if line.split()[4].lower() == "au"]
    constrained = [line.strip() for line in lines if "constrain_relaxation" in line]
    assert len(constrained) == len(gold) > 0
    assert "relax_freeze=electrodes" in capsys.readouterr().out

    summary = json.loads((run_dir / "au_relax" / "input_summary.json").read_text(
        encoding="utf-8"))
    assert summary["constrained"] == gold


def _prepare_relax(workdir) -> tuple[Path, Path]:
    cfg, run_dir = _write_project(workdir)
    assert _run("stage", "au60_build", "--prepare", "-c", str(cfg)) == 0
    assert _run("stage", "au_relax", "--prepare", "-c", str(cfg)) == 0
    return cfg, run_dir


def test_au_relax_collect_requires_normal_termination(workdir, capsys):
    cfg, run_dir = _prepare_relax(workdir)
    (run_dir / "au_relax" / "aims.out").write_text("converged but no goodbye\n",
                                                   encoding="utf-8")
    rc = _run("stage", "au_relax", "--collect", "-c", str(cfg))
    assert rc == 1
    assert "没有正常结束" in capsys.readouterr().err


def test_au_relax_collect_requires_next_step(workdir, capsys):
    cfg, run_dir = _prepare_relax(workdir)
    (run_dir / "au_relax" / "aims.out").write_text("Have a nice day.\n", encoding="utf-8")
    rc = _run("stage", "au_relax", "--collect", "-c", str(cfg))
    assert rc == 1
    assert "缺 geometry.in.next_step" in capsys.readouterr().err


def test_au_relax_collect_rejects_moved_frozen_atoms(workdir, capsys):
    cfg, run_dir = _prepare_relax(workdir)
    (run_dir / "au_relax" / "aims.out").write_text("Have a nice day.\n", encoding="utf-8")
    meta = json.loads((run_dir / "au60_build" / "electrode_metadata.json").read_text(
        encoding="utf-8"))
    lines = (run_dir / "au_relax" / "geometry.in").read_text(encoding="utf-8").splitlines()
    idx = meta["left_backplane_indices"][0] - 1  # 0-based
    atom_lines = [i for i, line in enumerate(lines) if line.startswith("atom")]
    parts = lines[atom_lines[idx]].split()
    parts[1] = "%.8f" % (float(parts[1]) + 0.2)
    lines[atom_lines[idx]] = "  " + " ".join(parts)
    (run_dir / "au_relax" / "geometry.in.next_step").write_text(
        "\n".join(lines) + "\n", encoding="utf-8")
    rc = _run("stage", "au_relax", "--collect", "-c", str(cfg))
    assert rc == 1
    assert "被冻结的原子发生了位移" in capsys.readouterr().err


def test_au_relax_collect_rejects_surface_mismatch(workdir, capsys):
    cfg, run_dir = _prepare_relax(workdir)
    (run_dir / "au_relax" / "aims.out").write_text("Have a nice day.\n", encoding="utf-8")
    (run_dir / "au_relax" / "geometry.in.next_step").write_text(
        (run_dir / "au_relax" / "geometry.in").read_text(encoding="utf-8"), encoding="utf-8")
    meta_path = run_dir / "au60_build" / "electrode_metadata.json"
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    meta["left_surface_indices"] = [1, 2, 3]  # 与 det-v1 重算结果必然不符
    meta_path.write_text(json.dumps(meta), encoding="utf-8")
    rc = _run("stage", "au_relax", "--collect", "-c", str(cfg))
    assert rc == 1
    assert "表面三元组" in capsys.readouterr().err


def test_au_relax_collect_ok_writes_relaxed_metadata(workdir, capsys):
    cfg, run_dir = _prepare_relax(workdir)
    (run_dir / "au_relax" / "aims.out").write_text("Have a nice day.\n", encoding="utf-8")
    geometry = (run_dir / "au_relax" / "geometry.in").read_text(encoding="utf-8")
    (run_dir / "au_relax" / "geometry.in.next_step").write_text(geometry, encoding="utf-8")
    rc = _run("stage", "au_relax", "--collect", "-c", str(cfg))
    assert rc == 0, capsys.readouterr()

    state = _state(run_dir)
    assert state["stages"]["au_relax"]["status"] == "collected"
    assert state["stages"]["au_relax"]["results"]["au_relax"]["enabled"] is True
    assert state["stages"]["au_relax"]["results"]["au_relax"]["surface_indices_verified"] is True

    meta = json.loads((run_dir / "au_relax" / "electrode_metadata.json").read_text(
        encoding="utf-8"))
    assert meta["relax"]["enabled"] is True
    assert meta["relax"]["frozen_atoms"] == 14
    assert len(meta["relax"]["input_fingerprint"]) == 64
    assert len(meta["relax"]["relaxed_fingerprint"]) == 64


# --------------------------------------------------------------------------- #
# 4) relax=false：本地直通
# --------------------------------------------------------------------------- #
def test_relax_disabled_is_local_only(workdir, capsys):
    cfg, run_dir = _write_project(workdir, relax=False)
    assert _run("stage", "au60_build", "--prepare", "-c", str(cfg)) == 0
    rc = _run("stage", "au_relax", "--prepare", "-c", str(cfg))
    assert rc == 0, capsys.readouterr()

    stage_state = _state(run_dir)["stages"]["au_relax"]
    assert stage_state["status"] == "prepared"
    assert stage_state["local_only"] is True
    assert stage_state["script"] in (None, "")
    assert (run_dir / "au_relax" / "relax_skipped.json").is_file()
    assert not (run_dir / "au_relax" / "control.in").exists()

    # 直通阶段没有 aims.out 也能 collect
    rc = _run("stage", "au_relax", "--collect", "-c", str(cfg))
    assert rc == 0, capsys.readouterr()
    results = _state(run_dir)["stages"]["au_relax"]["results"]["au_relax"]
    assert results["enabled"] is False


# --------------------------------------------------------------------------- #
# 5) nlayers: null -> 库文件建议值
# --------------------------------------------------------------------------- #
def test_nlayers_auto_from_library_metadata(workdir):
    meta_path = workdir / "electrode_metadata.json"
    meta_path.write_text(json.dumps({"library_cluster": {"nlayers_suggested": 3}}),
                         encoding="utf-8")
    resolve = AitranssTeStage._resolve_nlayers
    assert resolve(None, meta_path) == 3
    assert resolve(2, meta_path) == 2
    empty = workdir / "empty.json"
    empty.write_text("{}", encoding="utf-8")
    assert resolve(None, empty) == 3
    assert resolve(None, workdir / "missing.json") == 3
