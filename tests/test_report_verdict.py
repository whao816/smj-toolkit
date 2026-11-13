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
