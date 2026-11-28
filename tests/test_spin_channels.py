"""自旋分辨热输运（批 B / S2）。

开壳层此前的 ``spin`` 块只有 L0 与 P：能说"两个通道电导不同"，但说不出
"Sα 与 Sβ 各是多少、总 S 是否等于按电导加权的平均"。这里把每通道的
L1/L2 也积出来（与 ``analyze_te`` 共用同一套公式），并用一条 engineering
判据 ``spin_channel_consistency`` 守住"spin 块与总输运自洽"。
"""

from __future__ import annotations

import json
import random
from pathlib import Path

import pytest

from smj.junction.criteria import evaluate_run_criteria
from smj.junction.landauer import analyze_te
from smj.landauer_spin import analyze_te_spin
from smj.report import _render_spin_section
from smj.spin_summary import normalize


def _write_te(path: Path, rows) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["#bias = 0.0000 V efermi = -0.188452434085 H"]
    for energy, alpha, beta in rows:
        lines.append("%.12f %.6f %.9f %.9f"
                     % (-0.188452434085 + energy / 27.211386245988,
                        energy, alpha, beta))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _grid(lo: float = -0.6, hi: float = 0.6, step: float = 0.01):
    n = int(round((hi - lo) / step)) + 1
    return [lo + i * step for i in range(n)]


def _asymmetric_rows():
    """α 为 p 型（随 E 下降 -> Sα>0），β 为 n 型（随 E 上升 -> Sβ<0）。

    两通道的幅度/宽度取得不同，避免 Sα+Sβ 恰好为 0（那时 P_S 无定义）。
    """
    import math
    return [(e, 0.5 - 0.30 * math.tanh(e / 0.08), 0.5 + 0.20 * math.tanh(e / 0.12))
            for e in _grid()]


def test_per_channel_seebeck_follows_each_spin_spectrum(workdir: Path):
    te = _write_te(workdir / "TE.dat", _asymmetric_rows())
    total = analyze_te(te, temperature=300.0)
    spin = analyze_te_spin(te, temperature=300.0)

    assert spin["S_alpha_uV_per_K"] > 0
    assert spin["S_beta_uV_per_K"] < 0
    weighted = (spin["L0_alpha"] * spin["S_alpha_uV_per_K"]
                + spin["L0_beta"] * spin["S_beta_uV_per_K"]) / (
                    spin["L0_alpha"] + spin["L0_beta"])
    assert weighted == pytest.approx(total["S_uV_per_K"], rel=1e-9)
    assert spin["S_total_uV_per_K"] == pytest.approx(total["S_uV_per_K"], rel=1e-9)
    assert spin["spin_polarization_S"] == pytest.approx(
        (spin["S_alpha_uV_per_K"] - spin["S_beta_uV_per_K"])
        / (spin["S_alpha_uV_per_K"] + spin["S_beta_uV_per_K"]), rel=1e-12)
    # 每通道 k_e / ZT 与总口径同形（非负；数值为常数谱的 0 或小量）
    assert spin["k_e_alpha_pW_per_K"] >= 0 and spin["k_e_beta_pW_per_K"] >= 0
    assert spin["ZT_alpha_electronic"] >= 0 and spin["ZT_beta_electronic"] >= 0
    assert spin["g_over_g0_total"] == pytest.approx(
        0.5 * (spin["L0_alpha"] + spin["L0_beta"]))


def test_identical_channels_have_no_seebeck_polarization(workdir: Path):
    rows = [(e, 0.25, 0.25) for e in _grid()]
    te = _write_te(workdir / "TE.dat", rows)
    total = analyze_te(te, temperature=300.0)
    spin = analyze_te_spin(te, temperature=300.0)

    assert spin["S_alpha_uV_per_K"] == pytest.approx(0.0, abs=1e-6)
    assert spin["S_beta_uV_per_K"] == pytest.approx(0.0, abs=1e-6)
    assert spin["S_alpha_uV_per_K"] == pytest.approx(spin["S_beta_uV_per_K"], abs=1e-12)
    assert spin["S_total_uV_per_K"] == pytest.approx(total["S_uV_per_K"], abs=1e-9)
    assert spin["spin_polarization_S"] in (None,) or abs(spin["spin_polarization_S"]) < 1e-6
    assert spin["spin_polarization_L0"] == pytest.approx(0.0, abs=1e-12)
    for channel in ("alpha", "beta"):
        assert spin[f"k_e_{channel}_pW_per_K"] == pytest.approx(total["k_e_pW_per_K"] / 2)


def test_unsorted_channels_give_the_same_per_channel_seebeck(workdir: Path):
    rows = _asymmetric_rows()
    sorted_te = _write_te(workdir / "sorted" / "TE.dat", rows)
    shuffled = list(rows)
    random.Random(7).shuffle(shuffled)
    shuffled_te = _write_te(workdir / "shuffled" / "TE.dat", shuffled)

    a = analyze_te_spin(sorted_te, temperature=300.0)
    b = analyze_te_spin(shuffled_te, temperature=300.0)
    for key in ("L1_alpha", "L1_beta", "L2_alpha", "L2_beta",
                "S_alpha_uV_per_K", "S_beta_uV_per_K", "spin_polarization_S"):
        assert b[key] == pytest.approx(a[key], rel=1e-10, abs=1e-12), key
    assert b["input_energy_order"] == "unsorted"


def _consistent_results(te: Path, *, corrupt: bool = False) -> dict:
    result = dict(analyze_te(te, temperature=300.0))
    spin = dict(analyze_te_spin(te, temperature=300.0))
    if corrupt:
        spin["g_over_g0_total"] = float(spin["g_over_g0_total"]) * 1.5
    result["spin"] = spin
    return result


def test_spin_channel_consistency_criterion(workdir: Path):
    te = _write_te(workdir / "TE.dat", _asymmetric_rows())
    run = workdir
    (run / "results.json").write_text(
        json.dumps(_consistent_results(te), ensure_ascii=False), encoding="utf-8")

    checks = evaluate_run_criteria(run, temperature=300.0, system="BDT", multiplicity=2)
    by_name = {c["name"]: c for c in checks}
    good = by_name["spin_channel_consistency"]
    assert good["passed"] is True
    assert good["applicable"] is True
    assert good["severity"] == "error"

    (run / "results.json").write_text(
        json.dumps(_consistent_results(te, corrupt=True), ensure_ascii=False),
        encoding="utf-8")
    checks = evaluate_run_criteria(run, temperature=300.0, system="BDT", multiplicity=2)
    bad = {c["name"]: c for c in checks}["spin_channel_consistency"]
    assert bad["passed"] is False
    assert "missing" not in bad
    assert "G_total" in bad["detail"] or "g_over_g0_total" in bad["detail"]


def test_spin_channel_consistency_not_applicable_for_closed_shell(workdir: Path):
    te = _write_te(workdir / "TE.dat", _asymmetric_rows())
    (workdir / "results.json").write_text(
        json.dumps(_consistent_results(te), ensure_ascii=False), encoding="utf-8")
    checks = evaluate_run_criteria(workdir, temperature=300.0, system="BDT",
                                   multiplicity=1)
    by_name = {c["name"]: c for c in checks}
    assert by_name["spin_channel_consistency"]["applicable"] is False


def test_report_and_normalize_keep_per_channel_seebeck():
    spin = normalize({
        "channels": 2,
        "S_alpha_uV_per_K": 12.5,
        "S_beta_uV_per_K": -7.5,
        "S_total_uV_per_K": 2.5,
        "spin_polarization_S": 2.0,
    })
    assert spin["S_alpha_uV_per_K"] == pytest.approx(12.5)
    assert spin["S_beta_uV_per_K"] == pytest.approx(-7.5)
    assert spin["spin_polarization_S"] == pytest.approx(2.0)

    lines: list[str] = []
    _render_spin_section(lines, spin)
    text = "\n".join(lines)
    assert "Sα" in text and "Sβ" in text and "P_S" in text
    assert "+12.5000" in text and "-7.5000" in text
