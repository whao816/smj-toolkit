"""M5-T2/T3/T4：热导的判定、报告节与 thermal.json 契约。"""
from __future__ import annotations

import json

import numpy as np
import pytest


def _spectrum(n: int = 245):
    """0-24.4 meV 的平滑谱：带顶（19.5 meV）以上严格为零。"""
    from smj.phonon.landauer import thermal_conductance

    energy = np.arange(n) * 0.1                     # meV
    t_ph = np.where(energy <= 19.5, 0.9 * (1.0 - (energy / 19.5) ** 2), 0.0)
    g_ph = thermal_conductance(energy * 1e-3, t_ph, 300.0)
    cumulative = np.linspace(0.0, g_ph, energy.size)
    return {
        "energy_meV": energy.tolist(),
        "transmission": t_ph.tolist(),
        "cumulative_g_ph_pw_per_k": cumulative.tolist(),
    }, g_ph


def _thermal(**over):
    spectrum, g_ph = _spectrum()
    from smj.phonon.units import gq_pw_per_k

    payload = {
        "schema": "smj.thermal.v1",
        "temperature_K": 300.0,
        "n_atoms": 90,
        "n_device_atoms": 54,
        "n_lead_atoms": 36,
        "g_ph_pw_per_k": g_ph,
        "g_quantum_pw_per_k": gq_pw_per_k(300.0),
        "n_channels_effective": g_ph / gq_pw_per_k(300.0),
        "acoustic_sum_residual": 1.0e-9,
        "acoustic_sum_relative_residual": 1.0e-11,
        "imaginary_modes": {"count": 0, "worst_cm1": None, "modes": [],
                            "tolerance_cm1": -10.0, "n_modes": 270},
        "n_k_points": 1,
        "band_top_meV": 19.5,
        "eta_meV": 0.5,
        "energy_step_meV": 0.1,
        "spectrum": spectrum,
        "zt": {"available": False, "reason": "no electron result"},
        "fc_asymmetry": 1.0e-3,
    }
    payload.update(over)
    return payload


def _write(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    return path


# -- M5-T4:契约 --------------------------------------------------------------- #

def test_read_thermal_accepts_valid_payload(workdir):
    from smj.analysis import read_thermal

    path = _write(workdir / "thermal.json", _thermal())
    data, warnings = read_thermal(path)
    assert data["schema"] == "smj.thermal.v1"
    assert warnings == []


def test_read_thermal_rejects_missing_key_and_bad_schema(workdir):
    from smj.analysis import read_thermal

    payload = _thermal()
    payload.pop("g_ph_pw_per_k")
    path = _write(workdir / "bad1.json", payload)
    with pytest.raises(ValueError, match="缺必需键"):
        read_thermal(path)

    path = _write(workdir / "bad2.json", _thermal(schema="smj.thermal.v0"))
    with pytest.raises(ValueError, match="schema"):
        read_thermal(path)


def test_read_thermal_rejects_ragged_spectrum(workdir):
    from smj.analysis import read_thermal

    payload = _thermal()
    payload["spectrum"]["transmission"] = payload["spectrum"]["transmission"][:-1]
    path = _write(workdir / "ragged.json", payload)
    with pytest.raises(ValueError, match="长度不一致"):
        read_thermal(path)


def test_read_thermal_warns_on_unknown_key(workdir):
    from smj.analysis import read_thermal

    path = _write(workdir / "extra.json", _thermal(mystery=1))
    _, warnings = read_thermal(path)
    assert warnings and "mystery" in warnings[0]


# -- M5-T2:判据评估 ------------------------------------------------------------ #

def test_evaluate_thermal_inconclusive_without_relax_evidence():
    from smj.phonon.checks import evaluate_thermal

    outcome = evaluate_thermal(_thermal())
    # 缺 max|F| 与 lead 胞证据 -> 不判 passed(方案 §11)
    assert outcome["verdict"] == "inconclusive"
    assert "ph_relax_converged" in outcome["inconclusive"]
    assert "ph_lead_cell_match" in outcome["inconclusive"]


def test_evaluate_thermal_passes_when_evidence_complete():
    from smj.phonon.checks import evaluate_thermal

    outcome = evaluate_thermal(
        _thermal(),
        relax_max_force=1.0e-3,
        relax_sc_accuracy_forces="5E-3",
        lead_cell_match=True,
    )
    assert outcome["verdict"] == "passed", outcome
    assert outcome["failed"] == [] and outcome["inconclusive"] == []


def test_grid_drift_threshold_follows_measured_noise_floor(monkeypatch):
    """D14: 门限按实测求积噪声(2.0%-2.5%)取 5%, 不再是不达的 1%。

    直接调门限比构造"恰好 3% 漂移"的谱更稳: 该检查的唯一可调量就是门限,
    而漂移本身由谱形状决定。两条边界都要钉住 —— 太小会误判, 太大等于不查。
    """
    from smj.phonon import checks

    assert checks.MAX_GRID_DRIFT == 5.0e-2

    monkeypatch.setattr(checks, "MAX_GRID_DRIFT", 1.0e-12)
    tight = checks.evaluate_thermal(
        _thermal(), relax_max_force=1.0e-3,
        relax_sc_accuracy_forces="5E-3", lead_cell_match=True)
    assert "ph_grid_independence" in tight["inconclusive"]

    monkeypatch.setattr(checks, "MAX_GRID_DRIFT", 1.0)
    loose = checks.evaluate_thermal(
        _thermal(), relax_max_force=1.0e-3,
        relax_sc_accuracy_forces="5E-3", lead_cell_match=True)
    assert "ph_grid_independence" not in loose["inconclusive"]


def test_evaluate_thermal_fails_above_band_top():
    from smj.phonon.checks import evaluate_thermal

    payload = _thermal()
    t_ph = list(payload["spectrum"]["transmission"])
    t_ph[-1] = 1.0e-3                      # 带顶之上冒出一个透射
    payload["spectrum"]["transmission"] = t_ph
    outcome = evaluate_thermal(
        payload, relax_max_force=1.0e-3, relax_sc_accuracy_forces="5E-3",
        lead_cell_match=True)
    assert outcome["verdict"] == "failed"
    assert "ph_band_top_zero" in outcome["failed"]


def test_evaluate_thermal_inconclusive_on_imaginary_modes():
    from smj.phonon.checks import evaluate_thermal

    payload = _thermal(imaginary_modes={"count": 2, "worst_cm1": -35.0, "modes": [],
                                        "tolerance_cm1": -10.0, "n_modes": 270})
    outcome = evaluate_thermal(
        payload, relax_max_force=1.0e-3, relax_sc_accuracy_forces="5E-3",
        lead_cell_match=True)
    assert outcome["verdict"] == "inconclusive"
    assert "ph_imaginary_modes" in outcome["inconclusive"]


# -- M5-T2:review 集成 --------------------------------------------------------- #

def test_review_reports_thermal_absent_without_json(workdir):
    from smj.review import review_run

    run = workdir / "run"
    run.mkdir(parents=True, exist_ok=True)
    result = review_run(run)
    data = dict(result)
    assert data["thermal"]["available"] is False
    assert "thermal_absent" in [i["code"] for i in data["issues"]]


def test_review_marks_inconclusive_and_never_passed(workdir):
    from smj.review import review_run

    run = workdir / "run"
    _write(run / "thermal_analyze" / "thermal.json", _thermal())
    data = dict(review_run(run))
    assert data["thermal"]["available"] is True
    assert data["thermal"]["verdict"] == "inconclusive"
    codes = [i["code"] for i in data["issues"]]
    assert "thermal_inconclusive" in codes
    # 判定不得是 passed:thermal_inconclusive 必须是 warning 级
    thermal = [i for i in data["issues"] if i["code"] == "thermal_inconclusive"]
    assert thermal and thermal[0]["severity"] == "warning"
    assert data["status"] != "passed"


def test_review_blocks_on_above_band_transmission(workdir):
    from smj.review import review_run

    run = workdir / "run"
    payload = _thermal()
    t_ph = list(payload["spectrum"]["transmission"])
    t_ph[-1] = 1.0e-2
    payload["spectrum"]["transmission"] = t_ph
    _write(run / "thermal_analyze" / "thermal.json", payload)
    data = dict(review_run(run))
    assert data["thermal"]["verdict"] == "failed"
    assert data["status"] == "blocked"


# -- M5-T3:报告节 -------------------------------------------------------------- #

def test_report_has_thermal_section_and_geometry_warning(workdir):
    from smj.report import render_report

    run = workdir / "run"
    _write(run / "thermal_analyze" / "thermal.json", _thermal())
    text = render_report({"run_id": "t", "root": str(run)})
    assert "## 5. 热导（声子）" in text
    assert "不可并列" in text
    assert "## 6. 失败原因" in text
    assert "## 7. 产物清单" in text


def test_report_without_thermal_does_not_crash(workdir):
    from smj.report import render_report

    run = workdir / "run"
    run.mkdir(parents=True, exist_ok=True)
    text = render_report({"run_id": "t", "root": str(run)})
    assert "未跑热导链" in text
    assert "## 6. 失败原因" in text


# -- D13:ASR 判据改用相对残差 ---------------------------------------------------- #

def _asr_check(payload):
    from smj.phonon.checks import evaluate_thermal
    report = evaluate_thermal(payload)
    for item in report["checks"]:
        if item["name"] == "ph_asr_residual":
            return item
    raise AssertionError("ph_asr_residual 没出现在检查结果里")


def test_asr_judged_on_relative_residual_not_absolute():
    """绝对残差 1e-2 eV/A^2 但相对只有 1e-5 —— 旧口径会误判失败。"""
    item = _asr_check(_thermal(acoustic_sum_residual=1.0e-2,
                               acoustic_sum_relative_residual=1.0e-5))
    assert item["passed"] is True
    assert "相对残差" in item["detail"]


def test_asr_fails_when_relative_residual_exceeds_threshold():
    item = _asr_check(_thermal(acoustic_sum_residual=1.0e-2,
                               acoustic_sum_relative_residual=5.0e-3))
    assert item["passed"] is False


def test_asr_falls_back_to_absolute_for_old_payloads():
    """旧 thermal.json 没有相对残差时,退回 1e-6 绝对口径而不是判缺证据。"""
    payload = _thermal(acoustic_sum_residual=1.0e-9)
    payload.pop("acoustic_sum_relative_residual")
    item = _asr_check(payload)
    assert item["passed"] is True
    assert "无相对残差记录" in item["detail"]



# -- M5-T3:谱图与分位 ---------------------------------------------------------- #

def test_cumulative_quantiles_interpolates_and_is_monotonic():
    from smj.plotting import cumulative_quantiles

    quantiles = cumulative_quantiles([0.0, 1.0, 2.0, 3.0], [0.0, 1.0, 2.0, 4.0])
    assert quantiles["25%"] == pytest.approx(1.0)
    assert quantiles["50%"] == pytest.approx(2.0)
    assert quantiles["75%"] == pytest.approx(2.5)
    assert quantiles["90%"] == pytest.approx(2.8)
    ordered = [quantiles[k] for k in ("25%", "50%", "75%", "90%")]
    assert ordered == sorted(ordered)


def test_cumulative_quantiles_empty_when_nothing_conducts():
    from smj.plotting import cumulative_quantiles

    assert cumulative_quantiles([0.0, 1.0], [0.0, 0.0]) == {}


def test_report_embeds_spectrum_figure_and_quantiles(workdir):
    """图与分位表要进报告——`G_ph` 是积分，看不出热流从哪段声子来。"""
    pytest.importorskip("matplotlib")
    from smj.report import render_report

    run = workdir / "run"
    _write(run / "thermal_analyze" / "thermal.json", _thermal())
    text = render_report({"run_id": "t", "root": str(run)})
    assert "![声子透射谱与累积热导](thermal_analyze/thermal_spectrum.png)" in text
    assert (run / "thermal_analyze" / "thermal_spectrum.png").is_file()
    assert "热流贡献分位" in text
    for label in ("25%", "50%", "75%", "90%"):
        assert "| %s |" % label in text
    # 分位必须在报告里给出波数,方便与振动模式指认对照
    assert "波数 (cm" in text
