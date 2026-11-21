"""S8：``smj spin`` 汇总子命令（批 D / 可选）。

汇总 run_spec 的分子自旋设置、TE.dat 的双通道积分（P(E)、Sα/Sβ）、
aims_scf/spin_partition.json 的 Hirshfeld 分区与 review 的自旋提示，
输出一张表或 ``--json``。
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from smj.cli import EXIT_FAILED, EXIT_OK, main

H_EV = 27.211386245988
EF_HA = -0.188452434085


def _write_te(path: Path, *, channels: int = 2, n: int = 201) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["#bias = 0.000000 V efermi = %.12f H" % EF_HA]
    for i in range(n):
        e = -1.0 + 2.0 * i / (n - 1)
        alpha = 0.20 + 0.10 * math.tanh(e / 0.10)
        beta = 0.20 - 0.05 * math.tanh(e / 0.15)
        values = [alpha, beta] if channels == 2 else [0.5 * (alpha + beta)]
        lines.append("%.10f %.10f %s" % (
            EF_HA + e / H_EV, e, " ".join("%.6f" % v for v in values)))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_run(root: Path, *, multiplicity: int = 2, charge: int = 0,
               broken_symmetry: bool = False, channels: int = 2,
               partition: bool = True, legacy_spin: bool = False) -> Path:
    run = root / "run"
    _write_te(run / "aitranss_te" / "TE.dat", channels=channels)
    spec = {"schema": "smj.run-spec.v1", "config": {
        "molecule": {"charge": charge, "multiplicity": multiplicity,
                     "fixed_spin_moment": True,
                     "broken_symmetry": broken_symmetry},
        "transport": {"temperature": 300.0}}}
    (run / "run_spec.json").write_text(json.dumps(spec), encoding="utf-8")
    if legacy_spin:
        # 批 B 之前的产物：spin 块只有 L0/T(E_F)/P，没有分通道 S。
        (run / "junction_analyze").mkdir(parents=True, exist_ok=True)
        (run / "junction_analyze" / "results.json").write_text(json.dumps({
            "source": str(run / "aitranss_te" / "TE.dat"),
            "temperature_K": 300.0, "points": 201, "L0": 0.2,
            "spin": {"channels": 2, "points": 201, "temperature_K": 300.0,
                     "L0_alpha": 0.2, "L0_beta": 0.2,
                     "spin_polarization_at_ef": 0.01,
                     "spin_polarization_L0": 0.0},
        }), encoding="utf-8")
    if partition:
        (run / "aims_scf").mkdir(parents=True, exist_ok=True)
        (run / "aims_scf" / "spin_partition.json").write_text(json.dumps({
            "schema": "smj.spin-partition.v1", "source": "aims.out",
            "sidecar": "spin_partition.json", "total": 1.0,
            "atoms": [{"index": 1, "element": "C", "spin_moment": 0.125},
                      {"index": 2, "element": "C", "spin_moment": 0.125},
                      {"index": 3, "element": "Au", "spin_moment": 0.375},
                      {"index": 4, "element": "Au", "spin_moment": 0.375}],
            "missing": 0,
            "groups": {
                "molecule": {"atoms": [1, 2], "spin_moment": 0.25, "fraction": 0.25,
                             "n_positive": 2, "n_negative": 0, "net": 0.25,
                             "max_abs": 0.125},
                "electrode": {"atoms": [3, 4], "spin_moment": 0.75, "fraction": 0.75,
                              "n_positive": 2, "n_negative": 0, "net": 0.75,
                              "max_abs": 0.375}},
            "n_positive": 4, "n_negative": 0, "localized": False,
            "net": 1.0, "max_abs": 0.375,
            "molecule_net": 0.25, "molecule_max_abs": 0.125,
            "molecule_n_positive": 2, "molecule_n_negative": 0,
            "localized_antiparallel": False,
        }), encoding="utf-8")
    return run


def _payload(capsys) -> dict:
    out = capsys.readouterr().out
    return json.loads(out)


def test_spin_json_summarizes_channels_and_partition(workdir: Path, capsys) -> None:
    run = _write_run(workdir)
    assert main(["spin", str(run), "--json"]) == EXIT_OK
    payload = _payload(capsys)
    assert payload["molecule"]["multiplicity"] == 2
    assert payload["open_shell"] is True
    assert payload["spin"]["channels"] == 2
    assert payload["spin"]["S_alpha_uV_per_K"] is not None
    assert payload["spin"]["S_beta_uV_per_K"] is not None
    assert payload["spin"]["spin_polarization_S"] is not None
    assert payload["spin_source"] == "TE.dat"
    assert payload["partition"]["groups"]["molecule"]["fraction"] == pytest.approx(0.25)
    codes = {item["code"] for item in payload["findings"]}
    assert "open_shell_uncalibrated" in codes


def test_spin_text_shows_table_and_caveats(workdir: Path, capsys) -> None:
    run = _write_run(workdir)
    assert main(["spin", str(run)]) == EXIT_OK
    text = capsys.readouterr().out
    for needle in ("自旋汇总", "Sα", "Sβ", "P(E_F)", "自旋分布", "open_shell_uncalibrated"):
        assert needle in text


def test_spin_recomputes_from_te_for_legacy_results(workdir: Path, capsys) -> None:
    run = _write_run(workdir, legacy_spin=True)
    assert main(["spin", str(run), "--json"]) == EXIT_OK
    payload = _payload(capsys)
    # 旧 spin 块没有分通道 S：必须从 TE.dat 重算，而不是显示 "—"。
    assert payload["spin_source"] == "TE.dat"
    assert payload["spin"]["S_alpha_uV_per_K"] is not None
    assert payload["spin"]["S_beta_uV_per_K"] is not None


def test_spin_closed_shell_prints_message(workdir: Path, capsys) -> None:
    run = _write_run(workdir, multiplicity=1, channels=1, partition=False)
    assert main(["spin", str(run)]) == EXIT_OK
    text = capsys.readouterr().out
    assert "闭壳层" in text
    assert "自旋汇总" in text


def test_spin_missing_run_fails_cleanly(workdir: Path, capsys) -> None:
    assert main(["spin", str(workdir / "nope")]) == EXIT_FAILED
    err = capsys.readouterr().err
    assert "找不到" in err
