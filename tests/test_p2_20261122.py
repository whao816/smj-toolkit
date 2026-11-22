"""2026-11-22 复审 P2 批的回归：口径收窄与三处"静默失效"。

覆盖：
* 自旋 Seebeck 极化率在分母近抵消时判不可定义（而不是 +199999900%）；
* Au–S 判据带只有一个来源，建结目标带不再冒充判据带；
* cube 负步长轴归一化后锚定窗口不再静默为空；
* plot 的 --labels/--spin 不再静默截断/丢弃；
* verify 真的检查退役阶段模块；
* 独立积分复核进入 validate 的 checks；
* junction.review 的结构检查点被 engine 写出、被 review 消费；
* junction.manifest 委托 report.build_manifest（不再有第二套 verdict）。
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from smj.cube import BOHR, anchor_weight, read_cube
from smj.junction import criteria as crit
from smj.junction.independent_check import l0_l1_l2, parse_te
from smj.junction.manifest import build_manifest as junction_manifest
from smj.junction.review import record_confirmation, write_review_xyz
from smj.landauer_spin import analyze_te_spin
from smj.review import ReviewResult, _structure_checkpoint_issues

# --------------------------------------------------------------------------- #
# 自旋 P_S
# --------------------------------------------------------------------------- #
def _write_spin_te(path: Path, rows) -> None:
    lines = ["# E_Ha  E-Ef_eV  T_alpha  T_beta"]
    lines.extend("%.10f %.10f %.10f %.10f" % row for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_spin_polarization_near_cancel_is_undefined_not_huge(workdir):
    """Sα≈−Sβ：比值对噪声极端敏感，必须 None + 解释，不能报 1e6。"""
    rows = []
    n = 101
    for i in range(n):
        x = -0.5 + 1.0 * i / (n - 1)
        rows.append((0.0, x, 1.0 + 0.20 * x, 1.0 - 0.20 * x))
    te = workdir / "TE.dat"
    _write_spin_te(te, rows)
    out = analyze_te_spin(te)
    assert out["spin_polarization_S"] is None
    assert "近抵消" in (out["spin_polarization_S_note"] or "")


def test_spin_polarization_defined_for_ordinary_asymmetry(workdir):
    """两通道同号但不等：分母有实质大小，P_S 正常给出。"""
    rows = []
    n = 101
    for i in range(n):
        x = -0.5 + 1.0 * i / (n - 1)
        rows.append((0.0, x, 1.0 + 0.30 * x, 1.0 + 0.10 * x))
    te = workdir / "TE.dat"
    _write_spin_te(te, rows)
    out = analyze_te_spin(te)
    assert out["spin_polarization_S"] is not None
    assert -1.0 < out["spin_polarization_S"] < 1.0
    assert out["spin_polarization_S_note"] is None


# --------------------------------------------------------------------------- #
# Au–S 口径
# --------------------------------------------------------------------------- #
def test_au_s_criterion_band_has_one_source():
    from smj.junction import contacts, validate

    assert validate.AU_S_BAND == crit.AU_S_BOND_BAND
    assert crit.AU_S_BOND_BAND == (2.40, 2.44)
    entry = next(c for c in crit.CRITERIA if c["name"] == "au_s_bond")
    assert entry["band"] == list(crit.AU_S_BOND_BAND)
    # 建结目标带仍是 2.40–2.49，但它是另一个名字/另一件事
    assert contacts.TARGET_BAND_A == (2.40, 2.49)
    assert "不是" in contacts.__doc__
    readme = Path(__file__).resolve().parents[1] / "README.md"
    text = readme.read_text(encoding="utf-8")
    assert "2.40–2.44" in text and "不是判据带" in text


# --------------------------------------------------------------------------- #
# cube 负步长
# --------------------------------------------------------------------------- #
def _write_cube(path: Path, *, origin, steps, counts, gaussian_center=None):
    nx, ny, nz = counts
    atoms = [(16, 0.1, 0.1, 0.1)]   # 物理位置固定，不随网格 origin 移动
    lines = ["cube", "test"]
    lines.append("%d %.8f %.8f %.8f" % (len(atoms), *(c / BOHR for c in origin)))
    for axis, n in enumerate(counts):
        vector = [0.0, 0.0, 0.0]
        vector[axis] = steps[axis]
        lines.append("%d %.8f %.8f %.8f" % (n, *(v / BOHR for v in vector)))
    for z, x, y, zc in atoms:
        lines.append("%d 0.0 %.8f %.8f %.8f" % (z, x / BOHR, y / BOHR, zc / BOHR))
    chunk = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                x = origin[0] + i * steps[0]
                y = origin[1] + j * steps[1]
                z = origin[2] + k * steps[2]
                if gaussian_center is None:
                    value = 1.0 / (1.0 + (i + j + k))
                else:
                    r2 = (x - gaussian_center[0]) ** 2 + (y - gaussian_center[1]) ** 2 \
                        + (z - gaussian_center[2]) ** 2
                    value = pow(2.718281828, -r2 / 0.5)
                chunk.append("%14.8E" % value)
                if len(chunk) == 6:
                    lines.append(" ".join(chunk))
                    chunk = []
    if chunk:
        lines.append(" ".join(chunk))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_negative_step_axis_is_flipped_and_anchor_window_matches(workdir):
    counts = (6, 7, 5)
    origin = (0.0, 0.0, 0.0)
    steps = (0.2, 0.2, 0.2)
    pos = workdir / "pos.cube"
    _write_cube(pos, origin=origin, steps=steps, counts=counts,
                gaussian_center=(0.4, 0.2, 0.1))
    neg_origin = (0.0, origin[1] + (counts[1] - 1) * steps[1], 0.0)
    neg = workdir / "neg.cube"
    _write_cube(neg, origin=neg_origin, steps=(steps[0], -steps[1], steps[2]),
                counts=counts, gaussian_center=(0.4, 0.2, 0.1))
    a = read_cube(pos)
    b = read_cube(neg)
    assert all(s > 0 for s in b.steps)
    assert b.origin == pytest.approx(a.origin, abs=1e-7)
    assert b.steps == pytest.approx(a.steps, abs=1e-9)
    # 网格点集合与体素值都应与正步长版一致（不只是"不崩"）
    import numpy as np
    assert np.allclose(b.values, a.values)
    assert anchor_weight(b)["fraction"] == pytest.approx(anchor_weight(a)["fraction"],
                                                         rel=1e-9)


# --------------------------------------------------------------------------- #
# plot 保护
# --------------------------------------------------------------------------- #
def _te(path: Path) -> Path:
    path.write_text("# E T\n0.0 1.0\n0.01 0.9\n0.02 1.1\n", encoding="utf-8")
    return path


def test_plot_labels_must_match_inputs(workdir, capsys):
    from smj.cli import main

    a, b = _te(workdir / "a.dat"), _te(workdir / "b.dat")
    rc = main(["plot", str(a), str(b), "--labels", "only-one",
               "--out", str(workdir / "x.png")])
    assert rc == 1
    assert "--labels" in capsys.readouterr().err


def test_plot_overlay_refuses_spin_instead_of_silently_dropping(workdir, capsys):
    from smj.cli import main

    a, b = _te(workdir / "a.dat"), _te(workdir / "b.dat")
    rc = main(["plot", str(a), str(b), "--spin", "mean", "--out", str(workdir / "x.png")])
    assert rc == 1
    assert "不支持 --spin" in capsys.readouterr().err


# --------------------------------------------------------------------------- #
# verify legacy stages
# --------------------------------------------------------------------------- #
def test_verify_checks_legacy_stage_modules(capsys):
    from smj.cli import main

    main(["verify", "--json"])
    payload = json.loads(capsys.readouterr().out)
    rows = {row["check"]: row for row in payload["checks"]}
    assert "legacy-stages" in rows
    assert rows["legacy-stages"]["severity"] == "warning"
    assert rows["legacy-stages"]["ok"] is True


# --------------------------------------------------------------------------- #
# 独立复核进入 checks
# --------------------------------------------------------------------------- #
def _smooth_te(path: Path) -> Path:
    lines = []
    n = 101
    for i in range(n):
        x = -0.5 + 1.0 * i / (n - 1)
        lines.append("%.10f %.10f %.10f" % (0.0, x, 1.0 + 0.5 * x * x))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_independent_recheck_passes_and_flags_mismatch(workdir):
    from smj.junction.validate import _independent_recheck

    te = _smooth_te(workdir / "TE.dat")
    pipe = l0_l1_l2(parse_te(te), 300.0)
    ok = _independent_recheck(te, {"junction": pipe}, 300.0)
    assert ok["passed"] is True and ok["severity"] == "info"
    bad = dict(pipe)
    bad["L0"] = pipe["L0"] * 1.05
    mismatch = _independent_recheck(te, {"junction": bad}, 300.0)
    assert mismatch["passed"] is False and mismatch["severity"] == "warn"
    assert "不一致" in mismatch["detail"]
    na = _independent_recheck(te, None, 300.0)
    assert na["applicable"] is False


def test_independent_recheck_skips_tiny_spectra(workdir):
    from smj.junction.validate import _independent_recheck

    te = workdir / "TE.dat"
    te.write_text("0 0 1\n0.01 0.01 1\n0.02 0.02 1\n", encoding="utf-8")
    row = _independent_recheck(te, None, 300.0)
    assert row["applicable"] is False
    assert "点过少" in row["detail"]


# --------------------------------------------------------------------------- #
# 结构检查点接线
# --------------------------------------------------------------------------- #
def _geometry(path: Path, dx: float = 0.0) -> Path:
    path.write_text("atom %.6f 0.000000 0.000000 S\n"
                    "atom 2.400000 0.000000 0.000000 Au\n"
                    % (0.0 + dx), encoding="utf-8")
    return path


def test_review_consumes_structure_checkpoint_records(workdir):
    run = workdir / "run"
    stage = run / "aims_tip_relax"
    stage.mkdir(parents=True)
    _geometry(stage / "geometry.in")
    write_review_xyz(stage, "tips", stage="aims_tip_relax")

    result = ReviewResult()
    _structure_checkpoint_issues(result, run)
    codes = [issue["code"] for issue in result.issues]
    assert "structure_check_unconfirmed" in codes

    record_confirmation(stage, "tips", who="tester", stage="aims_tip_relax")
    result = ReviewResult()
    _structure_checkpoint_issues(result, run)
    confirmed = [issue for issue in result.issues
                 if issue["code"] == "structure_checkpoint_confirmed"]
    assert confirmed and "tester" in confirmed[0]["message"]

    # 几何被改过 -> 旧确认失效，必须重新看
    _geometry(stage / "geometry.in", dx=0.1)
    result = ReviewResult()
    _structure_checkpoint_issues(result, run)
    stale = [issue for issue in result.issues if issue["code"] == "structure_check_stale"]
    assert stale and stale[0]["severity"] == "warning"


def test_review_stage_writes_checkpoint_xyz_on_prepare(workdir):
    """预置的 geometry 经 write_review_xyz 后 review 能看见——保证写/读同一契约。"""
    stage = workdir / "au60_build"
    stage.mkdir()
    _geometry(stage / "geometry.in")
    info = write_review_xyz(stage, "electrodes", stage="au60_build")
    assert Path(info["path"]).is_file()
    result = ReviewResult()
    _structure_checkpoint_issues(result, workdir)
    assert any(i["code"] == "structure_check_unconfirmed" for i in result.issues)


# --------------------------------------------------------------------------- #
# manifest 合并
# --------------------------------------------------------------------------- #
def test_junction_manifest_delegates_to_report(monkeypatch, workdir):
    run = workdir / "run"
    (run / "aitranss_te").mkdir(parents=True)
    (run / "aitranss_te" / "TE.dat").write_text("0 0 1\n0.1 0.1 1\n", encoding="utf-8")
    sentinel = {"run_id": "r", "root": str(run), "spec": {"a": 1},
                "stages": {"aitranss_te": {"n_files": 1}}, "checks": [{"name": "x"}],
                "verdict": {"passed": True, "physics_verdict": "passed"},
                "totals": {"bytes": 1}, "judgment": {"schema": "s"}}
    calls = []

    def fake(target, **kwargs):
        calls.append((target, kwargs))
        return sentinel

    monkeypatch.setattr("smj.report.build_manifest", fake)
    man = junction_manifest(run, run_id="r", spec={"a": 1}, system="BDT")
    assert calls and calls[0][1]["include_checks"] is True
    assert man["verdict"] == sentinel["verdict"]
    assert man["checks"] == sentinel["checks"]
    assert "judgment" not in man          # 投影回本模块的公开形状


def test_junction_manifest_without_te_does_not_invent_verdict(workdir):
    run = workdir / "run"
    (run / "gaussian_opt").mkdir(parents=True)
    (run / "gaussian_opt" / "mol.log").write_text("x\n", encoding="utf-8")
    man = junction_manifest(run)
    assert man["verdict"] is None
    assert man["totals"]["stages"] == 1
