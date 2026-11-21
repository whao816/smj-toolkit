"""2026-11-21 第三轮复审（master 6384ce0d）回归 —— 批 1（R1）。

对应附件 ``test_smj_review_6384ce0.py`` 第 8 条：

    粗网格常数谱（12 点，步长 0.2 eV）在 300 K 下只有 1 个采样点落在 ±kBT 里，
    梯形求积把权重积分低估 42%、一阶矩变成 -1.96 kBT —— 于是常数谱也报出
    S=+169 uV/K。旧判据只看端点覆盖（f(E_min)-f(E_max)=1.0），照样判 passed。

R1 用"求积充分性"（权重积分/精确覆盖、权重一阶矩）来区分"网格够密"与
"网格太粗"，缺失时按证据不足（inconclusive）归约，不判物理否证。
"""

from __future__ import annotations

import pytest

from smj.junction.criteria import evaluate_run_criteria, physics_verdict
from smj.junction.landauer import analyze_te

EF_HARTREE = -0.18
H_EV = 27.2114


def write_constant_te(directory, energies):
    """常数透射谱 T=0.25：首列 Hartree（efermi=-0.18），第二列 E-E_F (eV)。"""
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "TE.dat"
    lines = ["#bias = 0.0000 V efermi = -0.18 H"]
    lines.extend("%.12f %.12f 0.25" % (EF_HARTREE + e / H_EV, e) for e in energies)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_underresolved_constant_spectrum_cannot_pass_with_spurious_seebeck(workdir):
    coarse = write_constant_te(workdir / "coarse", (-1.065 + i * 0.2 for i in range(12)))
    fine = write_constant_te(workdir / "fine", (-1.065 + i * 0.002 for i in range(1101)))

    reference = analyze_te(fine, temperature=300)
    assert reference["L0"] == pytest.approx(0.25, abs=1e-6)
    assert abs(reference["S_uV_per_K"]) < 1e-6

    result = analyze_te(coarse, temperature=300)
    checks = evaluate_run_criteria(coarse.parent, temperature=300, system="BDT")
    verdict = physics_verdict(checks)["physics_verdict"]

    correct = abs(result["L0"] - 0.25) < 1e-6 and abs(result["S_uV_per_K"]) < 1e-6
    assert correct or verdict != "passed", (
        verdict, result["thermal_window"]["coverage"], result["L0"], result["S_uV_per_K"])

    # 判据层面：求积不足必须是 error 级、missing（证据不足），而不是静默放过。
    by_name = {c["name"]: c for c in checks}
    quad = by_name["te_fermi_quadrature"]
    assert quad["severity"] == "error"
    assert quad["passed"] is False
    assert quad["missing"] is True
    assert "求积" in quad["detail"] and "kBT" in quad["detail"]

    summary = physics_verdict(checks)
    assert summary["physics_verdict"] == "inconclusive"
    assert "te_fermi_quadrature" in summary["missing_evidence_criteria"]
