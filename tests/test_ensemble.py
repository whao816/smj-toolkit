"""构象自动分支（smj.ensemble）：生成批处理脚本与汇总散布。"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

pytest.importorskip("rdkit")

from smj import ensemble  # noqa: E402  （rdkit 缺失时整个文件跳过）

BDT = "Sc1ccc(S)cc1"
CHAIN = "CCCCCCCCCCCCCCCCS"


def test_write_ensemble_writes_xyz_manifest_and_script(workdir):
    written = ensemble.write_ensemble(BDT, workdir / "ens", keep=1)
    data = written["data"]
    assert data["kept"] == 1
    assert data["smiles"] == BDT
    assert data["conformers"][0]["atom_count"] == 14
    assert (workdir / "ens" / "conf_01.xyz").exists()
    manifest = json.loads((workdir / "ens" / "ensemble.json").read_text(encoding="utf-8"))
    assert manifest["caveat"].startswith("MMFF94s")
    script = (workdir / "ens" / "run_ensemble.sh").read_text(encoding="utf-8")
    assert "set -eu" in script
    assert 'cd "$(dirname "$0")"' in script
    assert "--run-dir run/conf_01" in script
    assert "--xyz conf_01.xyz" in script


def test_report_marks_missing_runs(workdir):
    ensemble.write_ensemble(BDT, workdir / "ens", keep=1)
    payload = ensemble.report_ensemble(workdir / "ens")
    assert payload["rows"] == []
    assert len(payload["missing"]) == 1
    assert payload["spread"]["gap_ev"] is None


def _fake_results(run_dir: Path, gap: float, force: float = 0.004) -> None:
    stage = run_dir / "molecule_opt"
    stage.mkdir(parents=True, exist_ok=True)
    (stage / "results.json").write_text(json.dumps({
        "molecule_opt": {
            "atom_count": 51,
            "converged": True,
            "max_force_ev_per_A": force,
            "frontier_orbitals": {"gap_ev": gap},
            "electron_count": 200,
            "species": {"actual_tier": "light"},
        }
    }), encoding="utf-8")


def test_report_computes_gap_spread(workdir):
    written = ensemble.write_ensemble(CHAIN, workdir / "ens", keep=2)
    entries = written["conformers"]
    assert len(entries) == 2
    _fake_results(Path(entries[0]["run_dir"]), gap=1.20)
    _fake_results(Path(entries[1]["run_dir"]), gap=1.05, force=0.006)

    payload = ensemble.report_ensemble(workdir / "ens")
    assert [row["rank"] for row in payload["rows"]] == [1, 2]
    assert payload["spread"]["gap_ev"] == pytest.approx(0.15, abs=1e-9)
    assert payload["spread"]["max_force_ev_per_A"] == pytest.approx(0.002, abs=1e-9)
    assert payload["rows"][0]["tier"] == "light"


def test_report_requires_manifest(workdir):
    with pytest.raises(FileNotFoundError):
        ensemble.report_ensemble(workdir / "nothing")


def test_cli_ensemble_generate_and_report(workdir, capsys):
    from smj.cli import main

    assert main(["ensemble", "--smiles", BDT, "--dir", str(workdir / "ens"),
                 "--keep", "1"]) == 0
    out = capsys.readouterr().out
    assert "构象分支" in out
    assert "run_ensemble.sh" in out

    assert main(["ensemble", "--dir", str(workdir / "ens"), "--report"]) == 0
    report = capsys.readouterr().out
    assert "尚无结果" in report


def test_cli_ensemble_bad_smiles(workdir, capsys):
    from smj.cli import EXIT_FAILED, main

    assert main(["ensemble", "--smiles", "not-a-smiles", "--dir", str(workdir)]) == EXIT_FAILED
