"""批 D / S5：开壳层协议自检（engineering, warn）与 no_applicable_gate 文案。

计划 docs/plan-radical-2026-11-21.md §3.5：
- 开壳层（M>1 或 BS 单重态）时，报告必须核对协议齐全：spin 块存在、channels==2、
  transport.extra_outputs 含 hirshfeld 或 spin_density；缺项记 warn，不升级为 error。
- 所有族闸门 not applicable 时，physics_verdict 的措辞要讲清"开壳层未校准"，
  而不是笼统的"无法判定"。
"""
from __future__ import annotations

import json

import pytest

from smj.junction.criteria import physics_verdict
from smj.junction.validate import validate_junction_run


def _write_te(directory, *, columns=4, points=101):
    directory.mkdir(parents=True, exist_ok=True)
    lines = ["#bias = 0.0000 V efermi = -0.188452434085 H"]
    for i in range(points):
        e = -0.5 + i * 0.01
        chans = [0.25, 0.25] if columns == 4 else [0.25]
        lines.append(
            "%.12f %.12f %s"
            % (-0.188452434085 + e / 27.211386245988, e,
               " ".join("%.6f" % t for t in chans))
        )
    (directory / "TE.dat").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return directory / "TE.dat"


def _run_dir(workdir, *, spin_block=True):
    run = workdir / "run"
    _write_te(run / "aitranss_te")
    junction = {"points": 101}
    if spin_block:
        junction["spin"] = {"channels": 2}
    (run / "junction_analyze").mkdir(parents=True, exist_ok=True)
    (run / "junction_analyze" / "results.json").write_text(
        json.dumps({"junction": junction}, ensure_ascii=False), encoding="utf-8")
    return run


def _check(verdict, name="open_shell_protocol"):
    return next(c for c in verdict["checks"] if c["name"] == name)


def _validate(run, *, multiplicity=2, broken_symmetry=False, extra_outputs=("hirshfeld",)):
    return validate_junction_run(
        run,
        results={"junction": {"points": 101, "spin": {"channels": 2}}},
        system="BDT", multiplicity=multiplicity, broken_symmetry=broken_symmetry,
        extra_outputs=extra_outputs,
    )


def test_open_shell_protocol_passes_when_complete(workdir):
    verdict = _validate(_run_dir(workdir))
    check = _check(verdict)
    assert check["passed"] is True
    assert check["severity"] == "info"
    assert verdict["engineering_passed"] is True


def test_open_shell_protocol_warns_without_spin_block(workdir):
    run = _run_dir(workdir, spin_block=False)
    verdict = validate_junction_run(
        run, results={"junction": {"points": 101}}, system="BDT",
        multiplicity=2, extra_outputs=("hirshfeld",),
    )
    check = _check(verdict)
    assert check["passed"] is False
    assert check["severity"] == "warn"
    assert "spin" in check["detail"]
    assert verdict["engineering_passed"] is True  # warn 不拉垮工程判定


def test_open_shell_protocol_warns_without_extra_outputs(workdir):
    verdict = _validate(_run_dir(workdir), extra_outputs=())
    check = _check(verdict)
    assert check["passed"] is False
    assert check["severity"] == "warn"
    assert "extra_outputs" in check["detail"]
    assert "hirshfeld" in check["detail"]


def test_open_shell_protocol_marks_closed_shell_not_applicable(workdir):
    verdict = _validate(_run_dir(workdir), multiplicity=1)
    check = _check(verdict)
    assert check["applicable"] is False
    assert "闭壳层" in check["detail"]


def _open_shell_gate_checks():
    """只有 applicable=False 的族闸门：归约到 no_applicable_gate。"""
    detail = "[适用域: 不适用(开壳层未校准, 仅记录); 体系族 BDT]"
    return [{"name": "g_in_theory_band", "value": 0.5, "passed": True,
             "severity": "info", "detail": detail, "applicable": False}]


def test_no_applicable_gate_wording_is_plain_by_default():
    verdict = physics_verdict(_open_shell_gate_checks())
    assert verdict["physics_verdict"] == "inconclusive"
    assert "no_applicable_gate" in verdict["inconclusive_criteria"]
    assert "开壳层判据未校准" not in verdict["wording"]


def test_no_applicable_gate_wording_mentions_open_shell_calibration():
    verdict = physics_verdict(_open_shell_gate_checks(), open_shell=True)
    assert verdict["physics_verdict"] == "inconclusive"
    assert "no_applicable_gate" in verdict["inconclusive_criteria"]
    assert "开壳层判据未校准" in verdict["wording"]
    assert "标定" in verdict["wording"]
    assert "g_in_theory_band" in verdict["wording"]
