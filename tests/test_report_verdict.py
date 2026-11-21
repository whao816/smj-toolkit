"""报告判定行：工程判定与总体判定不是一个量。

2026-11-13 的验收 run 工程链 7/7 collected，但物理判据 t_ef_magnitude /
g_in_theory_band 未过；manifest 里 engineering_passed=true、passed=false。
报告曾把 passed 写成"工程判定 未通过"，把"链路是否跑通"和"物理是否支持"
混为一谈。
"""

from __future__ import annotations

from smj.report import render_report


def _manifest(**verdict_overrides):
    verdict = {
        "passed": False,
        "engineering_passed": True,
        "errors": 2,
        "warnings": 1,
        "physics_verdict": "failed",
        "wording": "否证 H0",
        "allowed_action": "停",
        "failed_criteria": ["t_ef_magnitude", "g_in_theory_band"],
        "inconclusive_criteria": ["t_ef_no_peak_near_ef"],
    }
    verdict.update(verdict_overrides)
    return {
        "schema": "smj.report.v1",
        "run_id": "t",
        "root": "C:/run/t",
        "generated_at": "2026-11-13T00:00:00+08:00",
        "state": {"stages": {}, "first_incomplete": None},
        "inputs": {"mode": "slurm-ssh"},
        "totals": {"stages": 7, "gib": 0.1},
        "resources": {},
        "transport": {},
        "failures": [],
        "warnings": [],
        "verdict": verdict,
    }


def test_report_engineering_row_uses_engineering_passed():
    text = render_report(_manifest())
    assert "| 工程判定 | 通过 |" in text


def test_report_renders_review_blocked_and_connection_state():
    """2026-11-22 复审 P1：报告必须标注复核结论，不能只写工程判定。"""
    man = _manifest()
    man["review"] = {
        "path": "C:/run/t/review.json",
        "status": "blocked",
        "ok": False,
        "accepted": True,
        "errors": 2,
        "warnings": 1,
        "connection_check": {
            "present": True,
            "valid": False,
            "reason": "geometry.in 缺失，无法按当前几何复核",
            "deviation_geom_A": 0.05,
        },
        "flags": ["connection_check_geometry_changed"],
    }
    text = render_report(man)
    assert "## 0.2 结构复核 / 连接核对" in text
    assert "记录已失效" in text
    assert "0.050 Å" in text
    assert "复核未通过（blocked）" in text
    assert "connection_check_geometry_changed" in text
    assert "人工豁免 | 是" in text


def test_report_review_record_unavailable_is_not_silent():
    man = _manifest()
    man["review"] = {"path": "C:/run/t/review.json",
                     "error": "review.json 不可解析: bad"}
    text = render_report(man)
    assert "复核记录不可用" in text
    assert "review.json 不可解析" in text
    assert "工程判定 | 未通过" not in text


def test_report_still_reports_physics_failure():
    text = render_report(_manifest())
    assert "否证 H0" in text
    assert "t_ef_magnitude" in text

# --- 2026-11-20 复核 P1-6 / R7：完工度闸门 -------------------------------- #

def test_completeness_check_fails_when_run_incomplete():
    from smj.verdict import apply_completeness

    verdict = {
        "checks": [], "passed": True, "engineering_passed": True,
        "errors": 0, "warnings": 0, "physics_verdict": "inconclusive",
        "failed_criteria": [], "inconclusive_criteria": [],
    }
    out = apply_completeness(verdict, is_done=False,
                             first_incomplete="molecule_opt",
                             missing=["molecule_opt", "aims_scf"])

    assert out["passed"] is False
    assert out["engineering_passed"] is False
    assert out["errors"] == 1
    assert "run_incomplete" in out["failed_criteria"]
    assert out["completeness"]["is_done"] is False
    assert out["completeness"]["stages_missing"] == ["molecule_opt", "aims_scf"]
    # 物理判定是独立归约，完工度不改它。
    assert out["physics_verdict"] == "inconclusive"


def test_completeness_check_passes_when_done():
    from smj.verdict import apply_completeness

    verdict = {"checks": [], "passed": True, "engineering_passed": True,
               "errors": 0, "warnings": 0, "physics_verdict": "passed",
               "failed_criteria": []}
    out = apply_completeness(verdict, is_done=True, first_incomplete=None,
                             missing=[])

    assert out["passed"] is True
    assert out["engineering_passed"] is True
    assert out["errors"] == 0
    assert out["completeness"]["is_done"] is True


def test_report_marks_te_only_run_as_incomplete(workdir):
    """复核现场：只放 aitranss_te/TE.dat 的 run 曾经 passed=true。"""
    from smj.report import build_manifest

    run_dir = workdir / "demo"
    (run_dir / "aitranss_te").mkdir(parents=True)
    manifest = build_manifest(run_dir, run_id="demo", include_checks=True)

    verdict = manifest["verdict"]
    assert verdict["passed"] is False
    assert verdict["engineering_passed"] is False
    assert verdict["errors"] >= 1
    assert verdict["completeness"]["is_done"] is False
    assert verdict["completeness"]["first_incomplete"] == "molecule_opt"
    names = [c["name"] for c in manifest["checks"]]
    assert "run_incomplete" in names

    text = render_report(manifest)
    assert "| 完工度 | 未完成" in text
    assert "molecule_opt" in text


# --- 2026-11-20 复审 P2-7 / R7：工程错误、物理失败与警告分账 ---------------- #
def _split_verdict():
    return {
        "checks": [
            {"name": "scf_geometry_consistency", "kind": "engineering",
             "severity": "error", "passed": False, "applicable": True},
            {"name": "g_in_theory_band", "kind": "physical",
             "severity": "error", "passed": False, "applicable": True},
            {"name": "k_e_zt_order", "kind": "physical",
             "severity": "warn", "passed": False, "applicable": True},
            {"name": "te_window_covers_efermi", "kind": "physical",
             "severity": "error", "passed": True, "applicable": True},
        ],
        "passed": False, "engineering_passed": False,
        "errors": 2, "warnings": 0, "physics_verdict": "failed",
        "failed_criteria": ["g_in_theory_band"],
    }


def test_engineering_and_physics_failures_are_counted_separately():
    from smj.verdict import apply_completeness

    out = apply_completeness(_split_verdict(), is_done=True,
                             first_incomplete=None, missing=[])

    assert out["engineering_errors"] == 1
    assert out["physics_failures"] == 1
    assert out["errors"] == 2          # 总体口径仍是"所有 error 级失败"
    assert out["engineering_passed"] is False
    assert out["warnings"] == 1        # 源头的 severity 是 "warn"


def test_physics_failure_alone_does_not_fail_engineering():
    from smj.verdict import apply_completeness

    verdict = _split_verdict()
    verdict["checks"] = [c for c in verdict["checks"]
                         if c["name"] != "scf_geometry_consistency"]
    out = apply_completeness(verdict, is_done=True,
                             first_incomplete=None, missing=[])

    assert out["engineering_errors"] == 0
    assert out["engineering_passed"] is True   # 物理失败不该拖垮工程判定
    assert out["passed"] is False              # 总体判定仍然是否证


# --- 2026-11-20 复审 P2-6 / R6：finished ≠ 已验收 --------------------------- #
def test_finished_stage_is_awaiting_not_accepted():
    from smj.verdict import apply_completeness, stage_acceptance

    statuses = {"molecule_opt": "collected", "aims_scf": "finished",
                "aitranss_te": "pending"}
    awaiting, missing = stage_acceptance(
        statuses, ["molecule_opt", "aims_scf", "aitranss_te"])
    assert awaiting == ["aims_scf"]
    assert missing == ["aitranss_te"]

    out = apply_completeness(
        {"checks": [], "passed": True, "engineering_passed": True,
         "errors": 0, "warnings": 0, "physics_verdict": "inconclusive",
         "failed_criteria": []},
        is_done=False, first_incomplete="aims_scf",
        missing=missing, awaiting=awaiting)

    assert out["completeness"]["stages_awaiting"] == ["aims_scf"]
    assert out["passed"] is False and out["engineering_passed"] is False
    detail = out["checks"][-1]["detail"]
    assert "待验收" in detail and "aims_scf" in detail


def test_skipped_stage_counts_as_accepted():
    from smj.verdict import stage_acceptance

    awaiting, missing = stage_acceptance(
        {"molecule_opt": "skipped", "aims_scf": "collected"},
        ["molecule_opt", "aims_scf"])
    assert awaiting == [] and missing == []


def test_manifest_does_not_call_finished_stage_done(workdir):
    import json

    from smj.report import build_manifest
    from smj.state import STAGE_ORDER

    run_dir = workdir / "run"
    run_dir.mkdir()
    stages = {s: {"status": "collected"} for s in STAGE_ORDER}
    stages["aims_scf"] = {"status": "finished"}
    (run_dir / "state.json").write_text(
        json.dumps({"schema": "smj.state.v1", "stages": stages}),
        encoding="utf-8")

    manifest = build_manifest(run_dir, run_id="t", include_checks=True)
    verdict = manifest["verdict"]
    assert verdict["completeness"]["is_done"] is False
    assert "aims_scf" in verdict["completeness"]["stages_awaiting"]
    assert verdict["engineering_passed"] is False
    text = render_report(manifest)
    assert "待验收" in text


def test_manifest_reports_skipped_stage_reason(workdir):
    """skipped 属于"已验收"（有替代证据），报告要说明为什么跳过。"""
    import json

    from smj.report import build_manifest
    from smj.state import STAGE_ORDER

    run_dir = workdir / "run"
    run_dir.mkdir()
    stages = {s: {"status": "collected"} for s in STAGE_ORDER}
    stages["molecule_opt"] = {
        "status": "skipped",
        "note": "直接使用库里的优化结构（library: azulene-26-dithiol）",
    }
    (run_dir / "state.json").write_text(
        json.dumps({"schema": "smj.state.v1", "stages": stages}),
        encoding="utf-8")

    manifest = build_manifest(run_dir, run_id="t", include_checks=True)
    assert manifest["verdict"]["completeness"]["is_done"] is True
    assert manifest["verdict"]["completeness"]["stages_awaiting"] == []
    text = render_report(manifest)
    assert "library: azulene-26-dithiol" in text
