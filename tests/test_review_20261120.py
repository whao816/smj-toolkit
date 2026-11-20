"""2026-11-20 复审（0.2.2rc1 / 04739baa）回归：批 1（R2 → R3 → R1）。

三条断言，对应复核报告附件 ``test_smj_review_20261120.py`` 的第 1–3 条：

1. **R1 / P1-1** 窄能量窗：T=0.25 的常数谱、E−E_F∈[-0.025,+0.001] eV、300 K。
   窗口虽跨 0，但只覆盖约 23% 的费米导数权重 —— 截窗积分给出
   ``G/G0=0.0586``、``S=+38.4 µV/K``，真实值是 ``G=0.25``、``S=0``。
   修复后新判据 ``te_thermal_window_coverage`` 必须以 *证据不足*（missing）
   归约为 ``inconclusive``，不得判 passed。
2. **R3 / P1-2** 自旋通道：乱序 TE.dat 的 α/β 通道曾按文件行序积分，
   ``spin_sum`` 从 0.2498 变成 0.5413（主积分仍正确），同一份 results.json
   里出现两套数。修复后分通道与主积分都必须与行序无关。
3. **R2 / P1-3** 峰检测：``significant_peaks`` 直接把相邻行当能量相邻点，
   同一单调谱升序判 passed、乱序判 inconclusive。修复后积分/峰检测/绘图
   共用按 E 升序的行。

夹具统一用仓库 ``workdir``（``tmp_path`` 在 DSH 沙箱下不可写，见 conftest）。
"""
from __future__ import annotations

import math
import random
from pathlib import Path

import pytest

from smj.junction.criteria import evaluate_run_criteria, physics_verdict, significant_peaks
from smj.junction.landauer import (
    analyze_te,
    energy_order,
    read_te,
    read_te_ordered,
)
from smj.landauer_spin import analyze_te_spin

EF_HARTREE = -0.188452434085
H_EV = 27.211386245988


# --------------------------------------------------------------------------- #
# 夹具
# --------------------------------------------------------------------------- #
def _write_te(path: Path, rows, *, efermi_ha: float = EF_HARTREE) -> Path:
    """rows = [(E_Ha, E−E_F_eV, T)] 或 [(E_Ha, E−E_F_eV, Tα, Tβ)]。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["#bias = 0.000000 V efermi = %.12f H" % efermi_ha]
    lines += [" ".join("%.10f" % value if i < 2 else "%.10e" % value
                       for i, value in enumerate(row)) for row in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _constant(lo_ev: float, hi_ev: float, *, n: int = 200,
              t: float = 0.25) -> list[tuple[float, float, float]]:
    return [(EF_HARTREE + (lo_ev + (hi_ev - lo_ev) * i / (n - 1)) / H_EV,
             lo_ev + (hi_ev - lo_ev) * i / (n - 1), t) for i in range(n)]


def _tanh_spectrum(*, n: int = 201) -> list[tuple[float, float, float]]:
    """单调递减 tanh：T(E_F)=0.25（BDT 带内），S>0（p 型，seebeck_sign 过）。"""
    rows = []
    for i in range(n):
        e = -0.5 + 1.0 * i / (n - 1)
        rows.append((EF_HARTREE + e / H_EV, e, 0.25 - 0.2 * math.tanh(e / 0.05)))
    return rows


def _spin_rows(*, n: int = 201) -> list[tuple[float, float, float, float]]:
    """双通道谱：α 峰在 +0.12 eV、β 峰在 −0.05 eV，互不重叠。"""
    rows = []
    for i in range(n):
        e = -0.5 + 1.0 * i / (n - 1)
        alpha = 0.55 * math.exp(-((e - 0.12) / 0.08) ** 2) + 0.02
        beta = 0.45 * math.exp(-((e + 0.05) / 0.10) ** 2) + 0.015
        rows.append((EF_HARTREE + e / H_EV, e, alpha, beta))
    return rows


def _shuffled(rows, seed: int = 7):
    out = list(rows)
    random.Random(seed).shuffle(out)
    return out


def _run_dir(workdir: Path, name: str) -> Path:
    run = workdir / name
    (run / "aitranss_te").mkdir(parents=True)
    return run


# --------------------------------------------------------------------------- #
# R1 / P1-1：窄窗只覆盖 23% 费米权重
# --------------------------------------------------------------------------- #
def test_constant_spectrum_with_truncated_thermal_window(workdir: Path) -> None:
    run = _run_dir(workdir, "narrow")
    _write_te(run / "aitranss_te" / "TE.dat", _constant(-0.025, 0.001))

    checks = evaluate_run_criteria(run / "aitranss_te", temperature=300.0, system="BDT")
    by_name = {c["name"]: c for c in checks}
    cover = by_name["te_thermal_window_coverage"]
    assert cover["severity"] == "error", "覆盖不足必须是 error 级闸门"
    assert cover["passed"] is False
    assert cover.get("missing") is True, "采样不足是证据缺失，不是物理否证"
    assert "覆盖" in cover["detail"] and "kBT" in cover["detail"]

    verdict = physics_verdict(checks)
    assert verdict["physics_verdict"] == "inconclusive"
    assert verdict["missing_evidence_criteria"] == ["te_thermal_window_coverage"]
    assert verdict["failed_criteria"] == [], "证据不足不得写成'否证 H0'"


def test_wide_window_constant_spectrum_matches_analytic(workdir: Path) -> None:
    """宽窗（±0.5 eV ≈ ±19 kBT）覆盖≈1：G 回到解析值 0.25，S≈0。"""
    run = _run_dir(workdir, "wide")
    te = _write_te(run / "aitranss_te" / "TE.dat", _constant(-0.5, 0.5, n=501))

    result = analyze_te(te, temperature=300.0)
    assert result["L0"] == pytest.approx(0.25, rel=5e-3)
    assert abs(result["S_uV_per_K"]) < 2.0

    window = result["thermal_window"]
    assert window["coverage"] > 0.999
    assert window["e_min_ev"] == pytest.approx(-0.5, abs=1e-9)
    assert window["e_max_ev"] == pytest.approx(0.5, abs=1e-9)
    assert window["lower_kbt"] == pytest.approx(-0.5 / window["kbt_ev"], rel=1e-9)
    assert window["upper_kbt"] == pytest.approx(0.5 / window["kbt_ev"], rel=1e-9)

    checks = evaluate_run_criteria(run / "aitranss_te", temperature=300.0, system="BDT")
    # 只断言覆盖判据本身：seebeck_sign 是 error 闸门而这里 S≈0，符号无物理意义，
    # 拿整体 verdict 会把"覆盖正常"和"符号判据"两件事混在一起。
    cover = {c["name"]: c for c in checks}["te_thermal_window_coverage"]
    assert cover["passed"] is True and not cover.get("missing"), (
        "宽窗常数谱覆盖≈1，覆盖判据不得误判")


# --------------------------------------------------------------------------- #
# R3 / P1-2：自旋通道按能量排序后再积分
# --------------------------------------------------------------------------- #
def test_spin_sum_matches_total_for_unsorted_spectrum(workdir: Path) -> None:
    rows = _spin_rows()
    asc_te = _write_te(workdir / "spin_asc" / "TE.dat", rows)
    shuf_te = _write_te(workdir / "spin_shuf" / "TE.dat", _shuffled(rows))

    asc = analyze_te_spin(asc_te, temperature=300.0)
    shuf = analyze_te_spin(shuf_te, temperature=300.0)

    for key in ("g_over_g0_total", "g_over_g0_alpha", "g_over_g0_beta",
                "spin_polarization_L0", "spin_polarization_at_ef",
                "L0_alpha", "L0_beta"):
        assert shuf[key] == pytest.approx(asc[key], rel=1e-9, abs=1e-12), key

    # 分通道之和 == 主积分（同一 results.json 里两套数必须自洽）
    total = analyze_te(asc_te, temperature=300.0)["L0"]
    assert shuf["g_over_g0_total"] == pytest.approx(total, rel=1e-9)
    assert shuf["g_over_g0_alpha"] + shuf["g_over_g0_beta"] == pytest.approx(
        shuf["g_over_g0_total"], rel=1e-12)

    # 行序观测：让下游看得出"被排序纠正过"
    assert asc["input_energy_order"] == "ascending"
    assert shuf["input_energy_order"] == "unsorted"


# --------------------------------------------------------------------------- #
# R2 / P1-3：峰检测与积分共用同一套升序行
# --------------------------------------------------------------------------- #
def test_peak_check_is_invariant_to_row_order(workdir: Path) -> None:
    rows = _tanh_spectrum()
    asc_te = _write_te(workdir / "asc" / "TE.dat", rows)
    shuf_te = _write_te(workdir / "shuf" / "TE.dat", _shuffled(rows))

    asc_rows, asc_order = read_te_ordered(asc_te)
    shuf_rows, shuf_order = read_te_ordered(shuf_te)
    assert asc_order == "ascending"
    assert shuf_order == "unsorted"
    assert read_te(asc_te) == asc_rows, "升序文件不得被重排"
    assert shuf_rows == asc_rows, "乱序文件必须还原成同一套升序行"
    assert energy_order(shuf_rows) == "ascending"

    assert significant_peaks(asc_rows) == []
    assert significant_peaks(shuf_rows) == [], (
        "单调谱在**能量序**下没有显著峰；乱序造出的峰是行序伪影")

    va = physics_verdict(evaluate_run_criteria(asc_te.parent, temperature=300.0,
                                               system="BDT"))
    vb = physics_verdict(evaluate_run_criteria(shuf_te.parent, temperature=300.0,
                                               system="BDT"))
    assert va["physics_verdict"] == vb["physics_verdict"] == "passed"
    assert va["failed_criteria"] == vb["failed_criteria"] == []
    assert va["inconclusive_criteria"] == vb["inconclusive_criteria"] == []
