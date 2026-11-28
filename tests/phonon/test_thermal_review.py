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
        # M5 收尾新增：温度曲线与软模归属（解释性字段，可选）
        "g_ph_by_temperature": {
            "temperature_K": [10.0, 50.0, 100.0, 300.0, 800.0],
            "g_ph_pw_per_k": [0.5 * g_ph / 6.0, 0.7 * g_ph, 0.92 * g_ph,
                              g_ph, 1.01 * g_ph],
        },
        "soft_modes": [
            {"index": 0, "cm1": -21.3, "participation_ratio": 14.0,
             "regions": {"molecule": 0.995, "device_au": 0.005, "lead": 0.0},
             "top_atoms": [{"atom": 3, "element": "C", "weight": 0.2}],
             "bonds": [{"pair": [3, 4], "label": "C-C", "stretch": 0.0105}],
             "rigid_body": {"fraction": 0.14}},
        ],
        "soft_mode_impact_relative_g_ph": 0.002,
        "g_ph_symmetric_pw_per_k": g_ph,
        "g_ph_soft_removed_pw_per_k": g_ph * 1.002,
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
        _thermal(transport_asr={"relative_residual": 1e-5}),
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


def test_review_needs_above_band_samples():
    from smj.phonon.checks import evaluate_thermal

    spectrum, g_ph = _spectrum(101)
    result = evaluate_thermal(_thermal(spectrum=spectrum, g_ph_pw_per_k=g_ph))
    assert "ph_band_top_zero" in result["inconclusive"]


def test_review_quantum_bound_is_inclusive_and_channel_specific():
    from smj.phonon.checks import evaluate_thermal
    from smj.phonon.units import gq_pw_per_k

    payload = _thermal(g_ph_pw_per_k=3 * gq_pw_per_k(300))
    assert "ph_quantum_bound" not in evaluate_thermal(payload)["failed"]
    payload["g_ph_pw_per_k"] = 3.99 * gq_pw_per_k(300)
    assert "ph_quantum_bound" in evaluate_thermal(payload)["failed"]
    assert "ph_quantum_bound" not in evaluate_thermal(payload, max_channels=4)["failed"]


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
                               acoustic_sum_relative_residual=1.0e-5,
                               transport_asr={"relative_residual": 1.0e-5}))
    assert item["passed"] is True
    assert "相对残差" in item["detail"]


def test_asr_fails_when_relative_residual_exceeds_threshold():
    item = _asr_check(_thermal(acoustic_sum_residual=1.0e-2,
                               acoustic_sum_relative_residual=5.0e-3))
    assert item["passed"] is False


def test_raw_closure_cannot_mask_a_truncated_transport_violation():
    payload = _thermal(acoustic_sum_relative_residual=0.,
                       transport_matrix_protocol="reciprocal_principal_layer_v1",
                       transport_asr={"relative_residual": .02,
                                      "residual_after_ev_per_A2": 0.})
    assert _asr_check(payload)["passed"] is False
    del payload["transport_asr"]
    assert _asr_check(payload)["passed"] is None


def test_asr_falls_back_to_absolute_for_old_payloads():
    """旧文件的完整行残差可评估，但缺少输运分区证据。"""
    payload = _thermal(acoustic_sum_residual=1.0e-9)
    payload.pop("acoustic_sum_relative_residual")
    item = _asr_check(payload)
    assert item["passed"] is None
    assert "无相对残差记录" in item["detail"]
    assert "输运分区 ASR 未记录" in item["detail"]


def test_archived_negative_spectra_are_rejected_on_read(workdir):
    from smj.analysis import read_thermal

    payload = _thermal()
    payload["spectrum"]["transmission"][6] = -.1
    with pytest.raises(ValueError, match="负"):
        read_thermal(_write(workdir / "thermal.json", payload))



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


def test_report_prints_environment_and_convergence_axes(workdir):
    """环境指纹与收敛轴要进报告：一个数配上"哪套工具算的、扫过什么"才可复核。"""
    from smj.report import render_report

    run = workdir / "run"
    _write(run / "thermal_analyze" / "thermal.json", _thermal())
    _write(run / "thermal_analyze" / "environment.json", {
        "schema": "smj.environment.v1",
        "generated_at": "2026-11-25T00:00:00+00:00",
        "smj": {"version": "0.2.2rc2", "root": "/x/src/smj"},
        "python": {"version": "3.12.10", "implementation": "CPython",
                   "executable": "/usr/bin/python3"},
        "platform": {"system": "Linux", "machine": "x86_64",
                     "release": "4.18.0", "processor": "x86_64"},
        "packages": {"numpy": "2.3.0"},
        "junction_provenance_sha256": "p" * 64,
        "thermal": {},
        "binaries": {"aims": {"path": "/opt/aims", "sha256": "a" * 64,
                              "exists": True}},
        "species": {"directory": "/s", "files": {}, "digest": "s" * 64,
                    "present": True},
        "lead_fc": {"path": "/l/lead_fc.json", "sha256": "l" * 64},
    })
    payload = _thermal()
    payload["environment"] = {"schema": "smj.environment.v1",
                              "digest": "0123456789abcdef" + "0" * 48,
                              "file": "environment.json"}
    payload["convergence"] = {
        "fc_cutoff": {
            "values_A": [8.0, 10.0, 12.0],
            "rows": [
                {"cutoff_A": 8.0, "g_ph_pw_per_k": 15.92639, "rejected": None},
                {"cutoff_A": 10.0, "g_ph_pw_per_k": 15.92639, "rejected": None},
                {"cutoff_A": 12.0, "g_ph_pw_per_k": 15.92639, "rejected": None},
            ],
            "relative_spread": 0.0983,
            "effective_span_A": 7.082,
            "production_A": 10.0,
            "nominal_shortfall": None,
            "nominal_covers_span": True,
            "beyond_nominal_checked": False,
            "converged": None,
        },
        "lead_layers": {"status": "unmeasured", "current": 4,
                        "requires": "新建 6 层版"},
        "device_layers": {"status": "unmeasured", "current": 2,
                          "requires": "device_layers=1"},
    }
    _write(run / "thermal_analyze" / "thermal.json", payload)
    text = render_report({"run_id": "t", "root": str(run)})

    assert "**运行环境**" in text
    assert "环境指纹" in text
    assert "**收敛轴**" in text
    # 名义截断覆盖了块内全部非零对时，不能报「未收敛」；但也不能把
    # 「后处理看不到更远耦合」说成收敛证据。
    assert "截断半径未收敛" not in text
    assert "名义截断已覆盖" in text
    assert "无法后处理验证" in text
    assert "9.8%" in text
    assert "lead 层数" in text and "未测" in text


def test_report_embeds_temperature_curve_and_soft_modes(workdir):
    """`G_ph(T)` 曲线与软模归属表都要进报告。

    这两块回答的是两类不同问题：曲线是"300 K 够不够收敛"，软模表是"那 2 个
    虚频到底动了什么、要不要紧"。仅有 `G_ph` 一个数时两者都无从判断。
    """
    pytest.importorskip("matplotlib")
    from smj.report import render_report

    run = workdir / "run"
    _write(run / "thermal_analyze" / "thermal.json", _thermal())
    text = render_report({"run_id": "t", "root": str(run)})
    assert "thermal_conductance_vs_temperature.png" in text
    assert (run / "thermal_analyze"
            / "thermal_conductance_vs_temperature.png").is_file()
    assert "G_ph(T)" in text
    assert "10 K" in text and "800 K" in text
    # 软模表：分子/器件 Au/lead 三列占比 + 参与比 + 刚体占比
    assert "软模归属" in text
    assert "| 1 | -21.30 |" in text
    assert "99.5%" in text
    assert "投影掉" in text
    assert "互易性" in text


def test_report_handles_zero_conductance_curve(workdir):
    from smj.report import render_report

    run = workdir / "run"
    payload = _thermal(g_ph_pw_per_k=0., g_ph_by_temperature={
        "temperature_K": [10., 800.], "g_ph_pw_per_k": [0., 0.]})
    _write(run / "thermal_analyze" / "thermal.json", payload)
    assert "0.00 pW/K" in render_report({"run_id": "zero", "root": str(run)})


def test_plot_thermal_vs_temperature_reports_reference_and_plateau(workdir):
    pytest.importorskip("matplotlib")
    from smj.plotting import plot_thermal_vs_temperature

    payload = _thermal()
    result = plot_thermal_vs_temperature(payload, workdir / "curve.png",
                                         formats=("png",))
    assert result["ok"] is True
    assert result["g_at_reference_pw_per_k"] == pytest.approx(
        payload["g_ph_pw_per_k"])
    assert result["g_high_temperature_pw_per_k"] == pytest.approx(
        payload["g_ph_by_temperature"]["g_ph_pw_per_k"][-1])
    assert (workdir / "curve.png").is_file()


def test_plot_thermal_vs_temperature_without_curve_is_not_ok(workdir):
    from smj.plotting import plot_thermal_vs_temperature

    payload = _thermal()
    payload.pop("g_ph_by_temperature")
    result = plot_thermal_vs_temperature(payload, workdir / "curve.png")
    assert result["ok"] is False
    assert "g_ph_by_temperature" in result["reason"]
