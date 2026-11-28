import json

import pytest

from smj.siesta.kconv import compare, load_run, render_markdown


def _write_run(root, name, g_over_g0, *, fermi=0.2, include_g=True):
    directory = root / name / "junction_analyze"
    directory.mkdir(parents=True)
    payload = {"G_S": 7.748091729863649e-05 * g_over_g0, "fermi_eV": fermi}
    if include_g:
        payload["G_over_G0"] = g_over_g0
    (directory / "results.json").write_text(json.dumps(payload), encoding="utf-8")
    return root / name


def test_compare_passes_within_tolerance(tmp_path):
    base = _write_run(tmp_path, "base", 0.0018)
    refined = _write_run(tmp_path, "refined", 0.00181)
    payload = compare(base, refined, tolerance=0.01)
    assert payload["passed"] is True
    assert payload["metrics"]["G_over_G0"]["relative_change"] == pytest.approx(1.0 / 180.0)


def test_compare_fails_outside_tolerance(tmp_path):
    base = _write_run(tmp_path, "base", 0.0018)
    refined = _write_run(tmp_path, "refined", 0.0020)
    payload = compare(base, refined, tolerance=0.01)
    assert payload["passed"] is False
    assert payload["metrics"]["G_over_G0"]["relative_change"] == pytest.approx(0.111111, rel=1e-5)
    assert "failed" in render_markdown(payload)


def test_load_run_rejects_missing_metric(tmp_path):
    run = _write_run(tmp_path, "bad", 0.0018, include_g=False)
    with pytest.raises(ValueError, match="G_over_G0"):
        load_run(run)


def test_render_markdown_keeps_scope_statement(tmp_path):
    base = _write_run(tmp_path, "base", 0.0018)
    text = render_markdown(compare(base, base))
    assert "passed" in text
    assert "agree within tolerance" in text
    assert "not evidence of experimental validity" in text

    refined = _write_run(tmp_path, "refined", 0.0020)
    failed = render_markdown(compare(base, refined))
    assert "failed" in failed
    assert "differ beyond tolerance" in failed
    assert "not evidence of experimental validity" in failed
