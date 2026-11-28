import json

from smj.config import RunConfig, RunInfo
from smj.report import generate_report
from smj.review import review_run
from smj.siesta import get_stage
from smj.siesta.io import dump, sha
from smj.siesta.report import _physics_verdict, _screening_text, _verdict_reason, render
from smj.siesta.workflow import stage_protocol
from smj.state import RunLayout, RunState


def completed(tmp_path):
    cfg = RunConfig(RunInfo("report", dir=str(tmp_path), mode="local"))
    cfg.transport.backend = "transiesta"
    layout = RunLayout(cfg.run_dir)
    state = RunState(cfg.name, backend="transiesta")
    contacts = {"passed": True, "contacts": {"left": {"motif": "flat", "final_site": "fcc"},
                                              "right": {"motif": "adatom", "final_site": "apex"}}}
    for name in state.order:
        stage = get_stage(name)
        directory = layout.ensure_stage(name)
        files = list(stage.outputs)
        if name == "tbtrans_te":
            files += ["device.TBT.nc", "tbtrans.out"]
        if name == "transiesta_scf":
            files.append("siesta.out")
        for f in files:
            (directory / f).write_bytes(b"offline fixture")
        if name == "siesta_junction_relax":
            dump(directory / "contacts.json", contacts)
        if name == "transiesta_scf":
            dump(directory / "scf.json", {"open_boundary_converged": True, "charge_deviation_e": .003,
                                          "screening": "unverified"})
        if name == "junction_analyze":
            dump(directory / "results.json", {"backend": "transiesta", "G_S": .00001,
                 "G_over_G0": .129, "temperature_K": 300., "smearing_eV": .025,
                 "spin": "unpolarized", "spin_convention": "degenerate", "fermi_eV": -4.7,
                 "physical_verdict": "inconclusive", "contacts": contacts, "warnings": []})
        state.set(name, status="collected", outputs_sha256={f: sha(directory / f) for f in files},
                  siesta_parents={p: state.stage(p).outputs_sha256 for p in stage.depends},
                  siesta_protocol=stage_protocol(cfg, name))
    state.save(layout.state_path)
    return cfg, state


def test_report_and_review_show_only_verified_siesta_results(tmp_path):
    cfg, state = completed(tmp_path)
    data = generate_report(cfg, include_checks=True)
    assert data["manifest"]["transport"]["G_over_G0"] == .129
    assert data["manifest"]["verdict"]["physics_verdict"] == "inconclusive"
    assert list(data["manifest"]["state"]["stages"]) == list(state.order)
    assert all(s in data["report"] for s in ("G/G₀", "fcc", "apex", "device.TBT.nc", "siesta.out"))
    assert review_run(cfg).ok
    (cfg.run_dir / "tbtrans_te/device.TBT.nc").write_bytes(b"changed")
    data = generate_report(cfg, include_checks=True)
    assert data["manifest"]["transport"] is None and "0.129" not in data["report"]
    assert not data["manifest"]["verdict"]["engineering_passed"]
    assert not review_run(cfg, force=False).ok
    assert json.loads((cfg.run_dir / "review.json").read_text(encoding="utf-8"))["status"] == "blocked"


def test_partial_and_reconfigured_runs_do_not_report_old_results(tmp_path):
    cfg, state = completed(tmp_path)
    cfg.transport.temperature = 100.
    assert generate_report(cfg)["manifest"]["transport"] is None
    state.set("siesta_junction_relax", status="prepared")
    state.save(cfg.run_dir / "state.json")
    assert generate_report(cfg.run_dir)["manifest"]["transport"] is None
    assert not review_run(cfg.run_dir).ok


def test_saved_route_without_state_is_incomplete(tmp_path):
    cfg, _ = completed(tmp_path)
    (cfg.run_dir / "state.json").unlink()
    dump(cfg.run_dir / "run_spec.json", {"config": {"transport": {"backend": "transiesta"}}})
    result = generate_report(cfg.run_dir)
    assert result["manifest"]["backend"] == "transiesta"
    assert result["manifest"]["transport"] is None
    assert not review_run(cfg.run_dir).ok


def test_physics_verdict_requires_both_numeric_checks():
    passed = {"status": "passed"}
    assert _physics_verdict({"screening": passed}, {"passed": True}) == "passed"
    assert _physics_verdict({"screening": passed}, None) == "inconclusive"
    assert _physics_verdict({"screening": "unverified"}, {"passed": True}) == "inconclusive"
    assert _physics_verdict({"screening": passed}, {"passed": False}) == "failed"
    assert _physics_verdict({"screening": {"status": "failed"}}, None) == "failed"
    assert _verdict_reason({"screening": passed}, None) == "横向 k 点收敛证据尚未生成"
    assert _verdict_reason({"screening": passed}, {"passed": False}) == "横向 k 点相对变化超过容差"
    assert _verdict_reason({"screening": {"status": "failed"}}, {"passed": True}) == "屏蔽判据未通过"
    assert _verdict_reason({"screening": "unverified"}, {"passed": True}) == "屏蔽证据缺失"


def test_render_shows_kconv_table_and_verdict():
    data = {
        "run_id": "demo", "stages": {}, "stage_issues": {}, "contacts": None,
        "scf": None, "transport": None,
        "kconv": {
            "base": {"run": "base", "G_over_G0": 1.8e-3, "fermi_eV": 0.19},
            "refined": {"run": "refined", "G_over_G0": 1.81e-3, "fermi_eV": 0.191},
            "metrics": {"G_over_G0": {"relative_change": 0.005}},
            "tolerance": 0.01, "passed": True,
        },
        "verdict": {"physics_verdict": "passed"},
    }
    text = render(data)
    assert "横向 k 点收敛" in text
    assert "passed" in text
    assert "不构成实验有效性证据" in text

    data["kconv"]["passed"] = False
    data["verdict"]["physics_verdict"] = "failed"
    failed = render(data)
    assert "物理判定为 **failed**：横向 k 点相对变化超过容差。" in failed


def test_screening_text_keeps_small_deviations_visible():
    text = _screening_text({
        "status": "passed",
        "electrode_layer_spread_e": 1.3e-05,
        "sides": {"left": {"outermost_delta_e": -7.4e-06,
                           "second_outermost_delta_e": 6.1e-04}},
        "reason": "bulk-like",
        "criterion": "outer two layers within 0.05 e",
    })
    assert "1.30e-05" in text
    assert "-7.40e-06" in text
    assert "6.10e-04" in text
