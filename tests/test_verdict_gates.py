"""2026-11-20 复核 P1-2：能量窗越界与"无适用闸门"都不得判 passed。

复核证据：合成常数谱 T=0.25、200 点、E−E_F=+0.250…+0.449 eV、300 K、system=AZ
原来给出 ``physics_verdict=passed``（g_in_theory_band / seebeck_sign 记
applicable=False，没有任何闸门真正被评估过）。新增两条约束：

1. ``te_window_covers_efermi``（physical/error）：E−E_F 必须跨 0，否则 failed；
2. ``physics_verdict``：所有带 applies_to 的族闸门都 applicable=False → inconclusive。
"""
from __future__ import annotations

from pathlib import Path

import pytest

from smj._vendor.junction.criteria import evaluate_run_criteria, physics_verdict
from smj._vendor.junction.landauer import analyze_te

H_EV = 27.211386245988


def _write_te(path: Path, lo_ev: float, hi_ev: float, *, n: int = 200,
              t: float = 0.25, efermi_ha: float = -0.188452434085) -> Path:
    """写一份 E−E_F ∈ [lo, hi] eV 的常数透射谱（首列是 Ha 绝对能量）。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["#bias = 0.000000 V efermi = %.12f H" % efermi_ha]
    for i in range(n):
        e_ev = lo_ev + (hi_ev - lo_ev) * i / (n - 1)
        lines.append("%.10f %.10f %.10e" % (efermi_ha + e_ev / H_EV, e_ev, t))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _by_name(checks: list[dict]) -> dict[str, dict]:
    return {c["name"]: c for c in checks}


def test_window_above_ef_is_failed(workdir: Path) -> None:
    """整窗在 E_F 之上：te_window_covers_efermi 必须是 error 且未过 → failed。"""
    te = _write_te(workdir / "TE.dat", 0.250, 0.449)
    checks = evaluate_run_criteria(te.parent, temperature=300.0, system="AZ")
    cover = _by_name(checks)["te_window_covers_efermi"]
    assert cover["severity"] == "error"
    assert cover["passed"] is False
    assert "未跨越" in cover["detail"]
    assert physics_verdict(checks)["physics_verdict"] == "failed"


@pytest.mark.parametrize("lo,hi", [(0.25, 0.449), (-0.449, -0.25)])
def test_window_on_either_side_is_failed(workdir: Path, lo: float, hi: float) -> None:
    te = _write_te(workdir / "TE.dat", lo, hi)
    checks = evaluate_run_criteria(te.parent, temperature=300.0, system="BDT")
    assert _by_name(checks)["te_window_covers_efermi"]["passed"] is False
    assert physics_verdict(checks)["physics_verdict"] == "failed"


def test_window_crossing_ef_passes_coverage_check(workdir: Path) -> None:
    te = _write_te(workdir / "TE.dat", -2.0, 2.0)
    checks = evaluate_run_criteria(te.parent, temperature=300.0, system="AZ")
    assert _by_name(checks)["te_window_covers_efermi"]["passed"] is True


def test_az_constant_spectrum_is_inconclusive_not_passed(workdir: Path) -> None:
    """薁（AZ）没有可用的族闸门：常数谱必须 inconclusive，不能包装成"支持 H0"。"""
    te = _write_te(workdir / "TE.dat", -2.0, 2.0)
    checks = evaluate_run_criteria(te.parent, temperature=300.0, system="AZ")
    names = _by_name(checks)
    assert names["g_in_theory_band"]["applicable"] is False
    assert names["seebeck_sign"]["applicable"] is False
    verdict = physics_verdict(checks)
    assert verdict["physics_verdict"] == "inconclusive"
    assert "no_applicable_gate" in verdict["inconclusive_criteria"]
    assert "支持 H0" not in verdict["wording"]


def test_bdt_has_applicable_gate(workdir: Path) -> None:
    """BDT 是标定族：闸门适用，verdict 不能因为"无闸门"而 inconclusive。"""
    te = _write_te(workdir / "TE.dat", -2.0, 2.0)
    checks = evaluate_run_criteria(te.parent, temperature=300.0, system="BDT")
    verdict = physics_verdict(checks)
    assert "no_applicable_gate" not in verdict["inconclusive_criteria"]


def test_seebeck_zero_does_not_imply_pass(workdir: Path) -> None:
    """S≈0 的对称谱不能靠"非负即过"的 seebeck_sign 通过：AZ 下仍是 inconclusive。"""
    te = _write_te(workdir / "TE.dat", -2.0, 2.0)
    res = analyze_te(te, temperature=300.0)
    assert abs(res["S_uV_per_K"]) < 1e-6
    checks = evaluate_run_criteria(te.parent, temperature=300.0, system="AZ")
    assert physics_verdict(checks)["physics_verdict"] != "passed"
