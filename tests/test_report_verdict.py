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
