"""验收协议：tip 全原子弛豫 -> 扩展成完整团簇 -> 几何冻结。

用户口径（2026-11-13）：``aims_tip_relax`` 不做任何 Au–S 约束，全原子优化得到
可靠的 Au–S 键长/连接角；``electrode_build`` 只把每个 tip 扩展成完整金字塔团簇，
沿用已弛豫的接触几何；此后 SCF/输运不再弛豫。

本组测试锁住：
* 配置默认值 ``anchor_au_s_target=None`` / ``pin_contact=False`` / ``relax=False``；
* ``make_tip_relax_input`` 在默认口径下写出的 geometry 没有 ``constrain_relaxation``；
* ``aims_tip_relax.prepare`` 把配置里的 pin/target 原样传给输入生成器；
* 默认配置模板不再把 2.40 Å 或弛豫写成默认行为。
"""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

from smj.config import ElectrodeConfig, default_config_text, load_config
from smj.stages import aims_tip_relax as tip_stage
from smj.stages.aims_tip_relax import AimsTipRelaxStage
from smj._vendor.junction.electrode import append_au60
from smj._vendor.junction.inputs import make_tip_relax_input
from smj._vendor.junction.geometry import read_aims_geometry

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


def _species_dir(root: Path) -> Path:
    species = root / "species"
    species.mkdir(parents=True, exist_ok=True)
    for element in ("Au", "C", "H", "S"):
        (species / ("%s_default" % element)).write_text(
            "  species %s\n    # dummy\n" % element, encoding="utf-8")
    return species


def test_default_protocol_is_unpinned_and_unrelaxed_after_build():
    cfg = ElectrodeConfig()
    assert cfg.anchor_au_s_target is None
    assert cfg.pin_contact is False
    assert cfg.relax is False


def test_default_template_does_not_pin_or_relax():
    text = default_config_text(name="t", xyz="molecule.xyz")
    assert "anchor_au_s_target: null" in text
    assert "pin_contact: false" in text
    assert "relax: false" in text


def test_pinned_contact_is_flagged_as_control(workdir):
    """pinned 是明确要求才跑的对照口径，不是生产默认。"""
    cfg_path = workdir / "smj.yaml"
    text = default_config_text(name="t", xyz="molecule.xyz").replace(
        "pin_contact: false", "pin_contact: true")
    cfg_path.write_text(text, encoding="utf-8")
    cfg = load_config(cfg_path)
    assert cfg.electrode.pin_contact is True
    assert any("pinned 对照" in warning for warning in cfg.warnings)


def test_tip_input_without_target_has_no_constraints(workdir):
    xyz = workdir / "molecule.xyz"
    xyz.write_text(BDT_XYZ, encoding="utf-8")
    species = _species_dir(workdir)
    written = make_tip_relax_input(
        xyz, workdir / "tip", left_anchor=7, right_anchor=8,
        au_tip_distance=2.35, species_dir=species,
        anchor_au_s_target=None, pin_contact=False,
        env={"aims_bin": "/bin/true", "partition": "main", "ntasks": 1},
    )
    summary = written["summary"]
    assert summary["contact_pinned"] is False
    assert summary["anchor_au_s_target"] is None
    assert summary["placement_distance_A"] == 2.35
    geometry = Path(written["geometry"]).read_text(encoding="utf-8")
    assert "constrain_relaxation" not in geometry
    assert "Au" in geometry


def test_cluster_expansion_keeps_relaxed_contact_geometry(workdir):
    """扩展只加团簇体，不能移动已弛豫的 tip/S，也不能改 S–Au 距离。"""
    import math

    molecule = [line.split() for line in BDT_XYZ.splitlines()[2:] if line.strip()]
    tip_geometry = workdir / "tips.in"
    lines = []
    for element, x, y, z in molecule:
        lines.append("  atom %s %s %s %s" % (x, y, z, element))
    # 左侧 S 在 +x，右侧 S 在 -x；Au–S 故意放在 2.70 A（不是 2.40 目标）。
    lines.append("  atom  5.860000 0.000000 0.000000 Au")
    lines.append("  atom -5.860000 0.000000 0.000000 Au")
    tip_geometry.write_text("\n".join(lines) + "\n", encoding="utf-8")

    before = read_aims_geometry(tip_geometry)
    expanded = workdir / "expanded.in"
    append_au60(tip_geometry, expanded, au_per_electrode=18)
    after = read_aims_geometry(expanded)

    assert len(after) > len(before)
    for i, atom in enumerate(before):
        assert after[i].element == atom.element
        for got, want in zip(after[i].coord(), atom.coord()):
            assert abs(got - want) < 1e-6
    left_s, right_s = 6, 7
    left_tip, right_tip = 14, 15
    assert abs(math.dist(after[left_s].coord(), after[left_tip].coord()) - 2.70) < 1e-6
    assert abs(math.dist(after[right_s].coord(), after[right_tip].coord()) - 2.70) < 1e-6


def test_tip_stage_passes_unpinned_defaults(workdir, monkeypatch):
    xyz = workdir / "opt.xyz"
    xyz.write_text(BDT_XYZ, encoding="utf-8")
    species = _species_dir(workdir)
    captured: dict = {}

    def fake_make(*args, **kwargs):
        captured.update(kwargs)
        stage_dir = Path(args[1])
        stage_dir.mkdir(parents=True, exist_ok=True)
        (stage_dir / "geometry.in").write_text("# fake\n", encoding="utf-8")
        (stage_dir / "control.in").write_text("# fake\n", encoding="utf-8")
        (stage_dir / "sub_aims.sh").write_text("#!/bin/bash\n", encoding="utf-8")
        other = stage_dir.parent / "aims_tip_relax"
        other.mkdir(parents=True, exist_ok=True)
        (other / "geometry.in.next_step").write_text(
            "  atom 0 0 0 Au\n", encoding="utf-8")
        return {
            "geometry": str(stage_dir / "geometry.in"),
            "control": str(stage_dir / "control.in"),
            "submit_script": str(stage_dir / "sub_aims.sh"),
            "summary": {"left_anchor": 7, "right_anchor": 8,
                        "placement_distance_A": 2.35, "contact_pinned": False,
                        "anchor_au_s_target": None},
        }

    monkeypatch.setattr(tip_stage, "make_tip_relax_input", fake_make)
    stage = AimsTipRelaxStage()
    ctx = SimpleNamespace(
        config=SimpleNamespace(electrode=ElectrodeConfig()),
        upstream=lambda *a, **k: xyz,
        override_path=lambda *_: None,
        species_dir=lambda: species,
        estimate=lambda *a, **k: {},
        stage_env=lambda *a, **k: {"aims_bin": "/bin/true"},
        log=lambda *_: None,
    )

    prepared = stage.prepare(ctx, workdir / "aims_tip_relax")

    assert captured["pin_contact"] is False
    assert captured["anchor_au_s_target"] is None
    assert prepared.local_only is False
    assert prepared.notes
