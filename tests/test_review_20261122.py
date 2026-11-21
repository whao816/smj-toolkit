"""2026-11-22 复审修复批 1（判定链）的回归测试。

覆盖：
- manifest 记录判定指纹；判据变化/旧 manifest 必须重建，不得复用旧结论。
- 报告渲染 unregistered_checks 与 no_applicable_gate，不得再出现"无法判定 +
  全部判据通过"并存。
- spin_channel_consistency 在 <run>/junction_analyze/results.json 布局下生效。
- TE.dat 非有限值按证据损坏处理（te_present missing -> inconclusive），
  不得变成"否证 H0"。
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from smj.analysis import _write_json
from smj.junction.criteria import evaluate_run_criteria, physics_verdict
from smj.junction.landauer import analyze_te, read_te, transport_moments
from smj.report import generate_report, render_report
from smj.verdict import judgment_fingerprint


def _manifest(verdict: dict) -> dict:
    return {
        "schema": "smj.report.v1",
        "run_id": "t",
        "root": "C:/run/t",
        "generated_at": "2026-11-22T00:00:00+08:00",
        "state": {"stages": {}, "first_incomplete": None},
        "inputs": {},
        "totals": {"stages": 0, "gib": 0.0},
        "resources": {},
        "transport": {},
        "failures": [],
        "warnings": [],
        "verdict": verdict,
    }


# --------------------------------------------------------------------------- #
# P0-1：过期 manifest 不得复用
# --------------------------------------------------------------------------- #
def _write_manifest(path: Path, manifest: dict) -> None:
    path.write_text(json.dumps(manifest, ensure_ascii=False), encoding="utf-8")


def test_generate_report_rebuilds_manifest_when_fingerprint_changed(workdir):
    run = workdir / "run"
    run.mkdir()
    first = generate_report(run, include_checks=True, force=True)
    assert first["stale_manifest"] is None
    manifest_path = run / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["judgment"]["sha256"] == judgment_fingerprint()

    # 模拟"旧代码写下的 manifest"：把指纹改成别的值，报告文字也换成旧内容。
    manifest["judgment"]["sha256"] = "0" * 64
    _write_manifest(manifest_path, manifest)
    (run / "report.md").write_text("STALE REPORT", encoding="utf-8")

    second = generate_report(run, include_checks=True)
    assert second["stale_manifest"] == "判据实现已变化（judgment.sha256 不匹配）"
    assert second["warnings"]
    fresh = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert fresh["judgment"]["sha256"] == judgment_fingerprint()
    assert "STALE REPORT" not in (run / "report.md").read_text(encoding="utf-8")


def test_generate_report_rebuilds_legacy_manifest_without_fingerprint(workdir):
    run = workdir / "run"
    run.mkdir()
    generate_report(run, include_checks=True, force=True)
    manifest_path = run / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    del manifest["judgment"]
    _write_manifest(manifest_path, manifest)

    result = generate_report(run, include_checks=True)
    assert result["stale_manifest"] == "旧 manifest 未记录判据指纹"
    fresh = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert fresh["judgment"]["sha256"] == judgment_fingerprint()


def test_generate_report_reuses_fresh_manifest(workdir):
    run = workdir / "run"
    run.mkdir()
    generate_report(run, include_checks=True, force=True)
    before = (run / "manifest.json").read_text(encoding="utf-8")

    result = generate_report(run, include_checks=True)
    assert result["stale_manifest"] is None
    assert result["written"] == []
    assert (run / "manifest.json").read_text(encoding="utf-8") == before


# --------------------------------------------------------------------------- #
# P0-2：unregistered_checks / no_applicable_gate 渲染
# --------------------------------------------------------------------------- #
def test_report_renders_unregistered_checks():
    text = render_report(_manifest({
        "passed": False, "engineering_passed": True,
        "errors": 0, "warnings": 0,
        "physics_verdict": "inconclusive",
        "wording": "无法判定", "allowed_action": "停",
        "failed_criteria": [], "inconclusive_criteria": [],
        "missing_evidence_criteria": [],
        "unregistered_checks": ["unregistered_x"],
    }))
    assert "未登记判据" in text and "unregistered_x" in text
    # 有未登记判据时不得再写"全部判据通过"。
    assert "全部判据通过" not in text


def test_report_explains_no_applicable_gate_instead_of_warn_criteria():
    text = render_report(_manifest({
        "passed": False, "engineering_passed": True,
        "errors": 0, "warnings": 0,
        "physics_verdict": "inconclusive",
        "wording": "无法判定", "allowed_action": "停",
        "failed_criteria": [],
        "inconclusive_criteria": ["no_applicable_gate"],
        "missing_evidence_criteria": [],
        "unregistered_checks": [],
    }))
    assert "无适用族闸门" in text
    assert "未过的 warn 级判据" not in text


# --------------------------------------------------------------------------- #
# P0-3：spin_channel_consistency 路径回退
# --------------------------------------------------------------------------- #
def _write_te(path: Path, rows) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["# bias = 0.0000 V efermi = -0.188452434085 H"]
    for energy, alpha, beta in rows:
        lines.append("%.12f %.6f %.9f %.9f"
                     % (-0.188452434085 + energy / 27.211386245988,
                        energy, alpha, beta))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _te_rows():
    rows = []
    for i in range(121):
        e = -0.6 + 0.01 * i
        rows.append((e, 0.5 - 0.30 * math.tanh(e / 0.08),
                     0.5 + 0.20 * math.tanh(e / 0.12)))
    return rows


def test_spin_consistency_reads_junction_analyze_results(workdir):
    from smj.landauer_spin import analyze_te_spin

    run = workdir / "run"
    te = _write_te(run / "aitranss_te" / "TE.dat", _te_rows())
    result = dict(analyze_te(te, temperature=300.0))
    result["spin"] = analyze_te_spin(te, temperature=300.0)
    results_path = run / "junction_analyze" / "results.json"
    results_path.parent.mkdir(parents=True)
    results_path.write_text(json.dumps(result), encoding="utf-8")

    checks = evaluate_run_criteria(run / "aitranss_te", temperature=300.0,
                                   system="BDT", multiplicity=2)
    spin = {c["name"]: c for c in checks}["spin_channel_consistency"]
    assert spin["applicable"] is True
    assert spin["passed"] is True

    # 破坏内部恒等关系后，判据必须能读到并判失败。
    result["spin"] = dict(result["spin"])
    result["spin"]["g_over_g0_total"] = float(result["spin"]["g_over_g0_total"]) * 1.5
    results_path.write_text(json.dumps(result), encoding="utf-8")
    checks = evaluate_run_criteria(run / "aitranss_te", temperature=300.0,
                                   system="BDT", multiplicity=2)
    spin = {c["name"]: c for c in checks}["spin_channel_consistency"]
    assert spin["applicable"] is True
    assert spin["passed"] is False
    assert "不自洽" in spin["detail"]


# --------------------------------------------------------------------------- #
# P0-4：非有限透射值
# --------------------------------------------------------------------------- #
def test_te_nan_is_rejected_and_becomes_missing_evidence(workdir):
    run = workdir / "run"
    te = run / "aitranss_te" / "TE.dat"
    te.parent.mkdir(parents=True)
    te.write_text(
        "# efermi = -0.188 H\n"
        "-0.195 0.1 0.5 0.5\n"
        "-0.190 0.2 nan 0.5\n"
        "-0.185 0.3 0.5 0.5\n",
        encoding="utf-8")

    with pytest.raises(ValueError, match="非有限"):
        read_te(te)

    checks = evaluate_run_criteria(run)
    assert len(checks) == 1
    assert checks[0]["name"] == "te_present"
    assert checks[0]["missing"] is True

    verdict = physics_verdict(checks)
    assert verdict["physics_verdict"] == "inconclusive"
    assert "te_present" in verdict["missing_evidence_criteria"]
    assert verdict["failed_criteria"] == []


def test_transport_moments_never_returns_nan():
    out = transport_moments(float("nan"), 0.0, 0.0, 300.0)
    assert all(value is None for value in out.values())


def test_write_json_refuses_nan(workdir):
    target = workdir / "bad.json"
    with pytest.raises(ValueError):
        _write_json({"x": float("nan")}, target)
    assert not target.exists()


def test_negative_transmission_is_reported(workdir):
    te = workdir / "TE.dat"
    te.write_text(
        "# efermi = -0.188 H\n"
        "-0.190 -0.1 -0.2 -0.2\n"
        "-0.188 0.0 0.5 0.5\n"
        "-0.186 0.1 0.2 0.2\n",
        encoding="utf-8")
    result = analyze_te(te, temperature=300.0)
    assert result["negative_transmission_points"] == 1
    assert any("负透射" in w for w in result["warnings"])
