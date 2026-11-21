"""S6：自由基支线端到端集成回归（真实 run 资产，默认不跑）。

默认测试集在 ``pyproject.toml`` 里用 ``addopts = ... -m 'not integration'``
排除了本目录；要跑：

    python -m pytest -m integration tests/integration -q

资产 ``runs/ttf-cation/ttf-cation`` 是 2026-11-19 的旧产物：``spin_partition``
（批次 A）与分通道 S（批次 B）在写入时还不存在，因此这里**离线重算**这些新
字段（只改临时副本，不动仓库里的 run），再验证报告与复核能把它们端到端带
到输出上。这样既守住真实数据的数值，又不要求重跑昂贵的 AIMS/TRANSIESTA。
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from smj.landauer_spin import analyze_te_spin
from smj.report import build_manifest, render_report
from smj.review import review_run
from smj.spin_partition import default_groups, parse_hirshfeld, summarize
from smj.state import RunState

pytestmark = pytest.mark.integration

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "runs" / "ttf-cation" / "ttf-cation"
STAGES = ("molecule_opt", "aims_tip_relax", "electrode_build", "aims_scf",
          "aitranss_te", "junction_analyze")


def _require_run() -> None:
    if not (RUN / "state.json").is_file():
        pytest.skip("真实 run 资产不存在: %s" % RUN)


def _partition() -> dict:
    """重算 ttf-cation 的 Hirshfeld 分区（分子 1–14 / 接触区 15–16 / 电极 17–88）。"""
    table = parse_hirshfeld(RUN / "aims_scf" / "aims.out")
    groups = default_groups(14, 88, contact=(15, 16))
    summary = summarize(table, groups)
    assert summary is not None, "aims.out 里应有 Hirshfeld 段"
    return summary


def test_ttf_cation_is_open_shell_and_spin_lives_on_electrodes():
    """真实开壳层 run：六阶段齐全、双通道、极化近零、自旋主要在电极上。"""
    _require_run()
    state = RunState.load(RUN / "state.json")
    assert state.statuses() == {name: "collected" for name in STAGES}

    spin = analyze_te_spin(RUN / "aitranss_te" / "TE.dat", temperature=300.0)
    assert spin["channels"] == 2
    assert abs(spin["spin_polarization_at_ef"]) < 0.05
    assert abs(spin["spin_polarization_L0"]) < 0.05

    summary = _partition()
    assert summary["total"] == pytest.approx(1.0, abs=1e-3)   # 实测 0.99997791（Hirshfeld 舍入）
    molecule = summary["groups"]["molecule"]
    assert molecule["fraction"] < 0.01          # 实测 0.0056：自旋几乎不在分子上
    assert summary["groups"]["electrode"]["fraction"] > 0.9
    assert summary["localized_antiparallel"] is False


def test_regenerated_report_and_review_expose_spin_evidence(workdir: Path):
    """把新字段注入旧 run 的副本，验证报告与复核端到端可见（不改原 run）。"""
    _require_run()
    copy = workdir / "run"
    shutil.copytree(RUN, copy)

    results_path = copy / "junction_analyze" / "results.json"
    payload = json.loads(results_path.read_text(encoding="utf-8"))
    payload["spin"] = analyze_te_spin(copy / "aitranss_te" / "TE.dat", temperature=300.0)
    results_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
                            encoding="utf-8")

    state = RunState.load(copy / "state.json")
    summary = {**_partition(), "source": "aims.out", "sidecar": "spin_partition.json"}
    aims = state.stage("aims_scf")
    aims.results["aims_scf"]["spin_partition"] = summary
    (copy / "aims_scf" / "spin_partition.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    text = render_report(build_manifest(copy, state=state))
    assert "4.1 自旋通道" in text
    assert "Sα" in text and "Sβ" in text
    assert "自旋分布" in text
    assert "0.0056" in text

    review = review_run(copy, force=True)
    codes = review.codes()
    assert "open_shell_uncalibrated" in codes
    assert "spin_localized_off_molecule" in codes
