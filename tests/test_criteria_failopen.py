"""F2 回归：缺产物不得判 passed（fail-closed），且"证据缺失"不等于"否证 H0"。

背景（2026-11-11 本地补丁）：旧实现里 ``evaluate_run_criteria`` 在没有 TE.dat 时返回
``[]``，``physics_verdict([])`` 落到 ``else: 'passed'``，于是只有空 ``aitranss_te/``
目录的 run 会被判 ``passed: true / 支持 H0``。本文件把修复后的语义钉死：

1. ``te_present`` 是 error 级工程闸门，缺 TE.dat 时 validate 必须 ``passed: False``；
2. 物理判定必须是 ``inconclusive``（无法判定），而不是 ``failed``（否证 H0）——
   没有观测，就没有可否证的东西；
3. ``smj report`` 即使没有任何产物也必须给出 verdict，除非显式 ``--no-checks``。

注：本机 %TEMP% 处在 DSH 沙箱的写限制之外（mkdir 报 WinError 5），pytest 自带的
``tmp_path`` 会因此在 ``tmp_path_factory.mktemp`` 处失败，故这里用自管理的 ``workdir``。
"""
from __future__ import annotations

import shutil
import uuid
from pathlib import Path

import pytest

from smj._vendor.junction.criteria import (
    CRITERIA,
    evaluate_run_criteria,
    physics_verdict,
)
from smj._vendor.junction.validate import validate_junction_run
from smj.report import build_manifest, render_report


@pytest.fixture()
def workdir() -> Path:
    base = Path(".tmp") / ("t-" + uuid.uuid4().hex[:8])
    base.mkdir(parents=True)
    try:
        yield base
    finally:
        shutil.rmtree(base, ignore_errors=True)


def _spec(name: str) -> dict | None:
    for crit in CRITERIA:
        if crit["name"] == name:
            return crit
    return None


def test_te_present_registered_as_error_gate():
    spec = _spec("te_present")
    assert spec is not None, "te_present 必须在 CRITERIA 登记，否则归约会把它当未知判据"
    assert spec["kind"] == "engineering"
    assert spec["gate"] == "error"
    assert spec["implemented"] is True


def test_empty_checks_is_inconclusive_not_passed():
    out = physics_verdict([])
    assert out["physics_verdict"] == "inconclusive"
    assert out["wording"] != "支持 H0"
    assert out["inconclusive_criteria"] == ["no_checks"]


def test_missing_te_is_missing_evidence_not_refutation(workdir: Path):
    (workdir / "aitranss_te").mkdir()
    checks = evaluate_run_criteria(workdir)
    assert len(checks) == 1
    check = checks[0]
    assert check["name"] == "te_present"
    assert check["passed"] is False
    assert check["severity"] == "error"
    assert check["missing"] is True

    out = physics_verdict(checks)
    assert out["physics_verdict"] == "inconclusive"
    assert out["failed_criteria"] == []
    assert out["missing_evidence_criteria"] == ["te_present"]


def test_validate_empty_te_dir_is_not_passed(workdir: Path):
    # 复刻原始 bug 的最小现场：只有空 aitranss_te/，没有任何产物文件。
    (workdir / "aitranss_te").mkdir()
    out = validate_junction_run(workdir)
    assert out["passed"] is False
    assert out["engineering_passed"] is False
    assert out["errors"] >= 1
    assert out["physics_verdict"] == "inconclusive"
    assert out["allowed_action"] != "可进入发布(A5)/预测卡(A4)"


def test_validate_honours_explicit_temperature(workdir: Path):
    # 温度不参与 te_present 判定，但接口必须接受显式温度（生产路径传 transport.temperature）。
    (workdir / "aitranss_te").mkdir()
    out = validate_junction_run(workdir, temperature=77.0)
    assert out["passed"] is False
    assert out["physics_verdict"] == "inconclusive"


def test_report_manifest_always_has_verdict(workdir: Path):
    run_dir = workdir / "demo"
    (run_dir / "aitranss_te").mkdir(parents=True)
    manifest = build_manifest(run_dir, run_id="demo", include_checks=True)
    assert manifest["verdict"] is not None
    assert manifest["verdict"]["passed"] is False
    assert manifest["verdict"]["physics_verdict"] == "inconclusive"

    text = render_report(manifest)
    assert "## 0.1 判定" in text
    assert "无法判定" in text
    assert "证据缺失" in text


def test_report_no_checks_leaves_verdict_null(workdir: Path):
    run_dir = workdir / "demo"
    run_dir.mkdir()
    manifest = build_manifest(run_dir, run_id="demo", include_checks=False)
    assert manifest["verdict"] is None
    assert "未计算判定" in render_report(manifest)
