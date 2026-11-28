"""P0 展示层：report 的自旋小节与 review 的开壳层告警。

2026-11-20 复核后的 P0（``docs/自由基单分子结-实现对照与开发规划.md`` §6）：
此前 ``review.py``/``report.py`` 里 0 处提及自旋，P(E) 只躺在
``junction_analyze/results.json`` 与图里，评审者看不到"这是开壳层、自旋
通道占多少、解是否可信"。本文件锁住那两条展示路径。
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from smj.report import build_manifest, render_report
from smj.review import review_run

_SPIN = {
    "source": "aitranss_te/TE.dat",
    "channels": 2,
    "points": 200,
    "temperature_K": 300.0,
    "L0_alpha": 1.2e-5,
    "L0_beta": 0.4e-5,
    "g_over_g0_alpha": 1.2e-5,
    "g_over_g0_beta": 0.4e-5,
    "g_over_g0_total": 8.0e-6,
    "T_alpha_at_ef": 0.3486,
    "T_beta_at_ef": 0.1194,
    "spin_polarization_at_ef": 0.4897,
    "spin_polarization_L0": 0.5,
    "max_transmission_alpha": 0.40,
    "max_transmission_beta": 0.14,
}


def _results(run: Path, *, spin: dict | None) -> None:
    target = run / "junction_analyze"
    target.mkdir(parents=True, exist_ok=True)
    payload = {
        "points": 200,
        "temperature_K": 300.0,
        "L0": 1.6e-5,
        "log10_G_over_G0_like": -1.7959,
        "S_uV_per_K": 12.0,
        "k_e_pW_per_K": 0.0027,
        "ZT_electronic": 0.4,
        "fermi_nearest": {"T": 0.468},
    }
    if spin is not None:
        payload["spin"] = spin
    (target / "results.json").write_text(json.dumps(payload), encoding="utf-8")


def _te(run: Path, *, columns: int = 2, points: int = 40) -> None:
    target = run / "aitranss_te"
    target.mkdir(parents=True, exist_ok=True)
    lines = []
    for index in range(points):
        e_ev = -0.2 + 0.01 * index
        row = [-0.188452434085 + e_ev / 27.211386245988, e_ev]
        if columns >= 2:
            row += [0.30 + 0.0005 * index, 0.10 + 0.0002 * index]
        else:
            row.append(0.25)
        lines.append(" ".join("%.6f" % value for value in row))
    (target / "TE.dat").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _config(run: Path, *, charge: int, multiplicity: int,
            fixed_spin_moment: bool = True) -> None:
    (run / "mol.xyz").write_text("2\n\nS 0 0 0\nAu 0 0 3\n", encoding="utf-8")
    (run / "smj.yaml").write_text(
        "run:\n  name: spin-surface\n  mode: dry-run\n"
        "molecule:\n  xyz: mol.xyz\n  charge: %d\n  multiplicity: %d\n"
        "  fixed_spin_moment: %s\n" % (
            charge, multiplicity, "true" if fixed_spin_moment else "false"),
        encoding="utf-8")


def test_open_shell_report_has_spin_section(workdir: Path) -> None:
    run = workdir / "open"
    _results(run, spin=_SPIN)
    manifest = build_manifest(
        run, run_id="open",
        spec={"molecule": {"charge": 1, "multiplicity": 2,
                           "fixed_spin_moment": True}},
        include_checks=False)
    text = render_report(manifest)
    assert "| 电荷 / 多重度 | 1 / 2 |" in text
    assert "fixed_spin_moment（固定 Nα−Nβ = M−1）" in text
    assert "### 4.1 自旋通道（开壳层）" in text
    assert "| P(E_F) = (Tα−Tβ)/(Tα+Tβ) | +48.97% |" in text
    assert "| P(L0) 积分极化 | +50.00% |" in text
    assert "双通道均值" in text
    assert "EXPECTED_RESONANCE" in text
    # manifest 本身也带 spin 块，供下游程序读取
    assert manifest["transport"]["spin"]["multiplicity"] == 2


def test_closed_shell_report_has_no_spin_section(workdir: Path) -> None:
    run = workdir / "closed"
    _results(run, spin=None)
    manifest = build_manifest(
        run, run_id="closed",
        spec={"molecule": {"charge": 0, "multiplicity": 1}},
        include_checks=False)
    text = render_report(manifest)
    assert "| 电荷 / 多重度 | 0 / 1 |" in text
    assert "### 4.1 自旋通道" not in text
    assert "总自旋约束" not in text
    assert manifest["transport"].get("spin") is None


def test_review_warns_open_shell_uncalibrated(workdir: Path) -> None:
    run = workdir / "run"
    _te(run)
    _results(run, spin=_SPIN)
    _config(run, charge=1, multiplicity=2)
    result = review_run(run, force=True, review_config={"connection_check": False})
    codes = result.codes()
    assert "open_shell_uncalibrated" in codes
    assert "spin_fixed_constraint" in codes
    assert result.spin is not None
    assert result.spin["multiplicity"] == 2
    assert result.spin["spin_polarization_at_ef"] == pytest.approx(0.4897)
    # 判据族未校准只降级为 warning，不阻断
    assert result.status == "warning"
    assert result.ok is True


def test_review_warns_single_channel_te(workdir: Path) -> None:
    run = workdir / "run"
    _te(run, columns=1)
    _results(run, spin=None)
    _config(run, charge=1, multiplicity=2)
    result = review_run(run, force=True, review_config={"connection_check": False})
    assert "spin_channels_missing" in result.codes()
    assert result.spin is not None
    assert result.spin["channels"] == 1
    assert result.status == "warning"


def test_review_closed_shell_stays_quiet(workdir: Path) -> None:
    run = workdir / "run"
    _te(run, columns=1)
    _results(run, spin=None)
    _config(run, charge=0, multiplicity=1)
    result = review_run(run, force=True, review_config={"connection_check": False})
    codes = result.codes()
    assert "open_shell_uncalibrated" not in codes
    assert "spin_channels_missing" not in codes
    assert "spin_block_unexpected" not in codes


def test_review_two_channel_without_multiplicity_is_called_out(workdir: Path) -> None:
    """旧 run 没有 run_spec/multiplicity：两列谱不能被当成"确认闭壳层"。"""
    run = workdir / "legacy"
    _te(run, columns=2)
    _results(run, spin={**_SPIN, "multiplicity": None, "charge": None})
    result = review_run(run, force=True, review_config={"connection_check": False})
    codes = result.codes()
    assert "spin_channels_undeclared" in codes
    assert "spin_block_unexpected" not in codes
    assert "open_shell_uncalibrated" not in codes


@pytest.mark.parametrize("schema,scf_moment,expected", [
    (2, None, False), (3, None, True), (3, "released", False), (2, "0", True),
])
def test_bs_constraint_agrees_in_all_reports(workdir, schema, scf_moment, expected):
    from smj.spin_report import build_spin_report

    run = workdir / "bs"
    _results(run, spin=_SPIN)
    molecule = {"multiplicity": 1, "charge": 0, "broken_symmetry": True,
                "fixed_spin_moment": True}
    (run / "run_spec.json").write_text(json.dumps(
        {"protocol_schema": schema, "config": {"molecule": molecule}}), encoding="utf-8")
    if scf_moment is not None:
        scf = run / "aims_scf"
        scf.mkdir()
        control = "spin collinear\n"
        if scf_moment != "released":
            control += f"fixed_spin_moment {scf_moment}\n"
        (scf / "control.in").write_text(control, encoding="utf-8")
    manifest = build_manifest(run, include_checks=False)
    assert manifest["inputs"]["molecule"]["fixed_spin_moment"] is expected
    assert manifest["transport"]["spin"]["fixed_spin_moment"] is expected
    assert build_spin_report(run)["molecule"]["fixed_spin_moment"] is expected
    review = review_run(run, force=True)
    assert review.spin["fixed_spin_moment"] is expected
    assert ("spin_fixed_constraint" in review.codes()) is expected


def test_invalid_te_cannot_reinstate_cached_spin(workdir, capsys):
    from smj.cli import main, EXIT_FAILED

    run = workdir / "damaged"
    _te(run)
    (run / "aitranss_te" / "TE.dat").write_text("0 0 .1 .2\n.01 .1 .4\n", encoding="utf-8")
    _results(run, spin=_SPIN)
    assert main(["spin", str(run), "--json"]) == EXIT_FAILED
    assert "TE.dat" in capsys.readouterr().err
    review = review_run(run, force=True)
    assert "spin_te_unreadable" in review.codes()
    assert review.spin is None
    manifest = build_manifest(run, include_checks=False)
    assert manifest["transport"] is None
    assert any("TE.dat" in w for w in manifest["warnings"])
