"""CLI 契约回归（2026-11-22 复审批 3）。

守住四件用户能直接感知、此前全部失效的事：

* ``smj run --json`` 的 stdout 只能有 JSON（进度日志走 stderr）；
* 顶层未预期异常给一行可读错误 + 非 0 退出码，不把 traceback 当正常输出；
* ``smj report --strict`` 在判定失败/阻断时退出码为 1（默认仍是 0）；
* ``smj review --stage`` 真的透传到 ``review_run``（此前普通路径完全忽略）。
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from smj.cli import _make_runner, main
from smj.config import load_config


def _config(workdir):
    run = workdir / "run"
    run.mkdir()
    xyz = workdir / "mol.xyz"
    xyz.write_text("2\nprobe\nS 0.0 0.0 -1.0\nS 0.0 0.0 1.0\n",
                   encoding="utf-8")
    cfg = workdir / "smj.yaml"
    cfg.write_text(
        "run:\n  name: probe\n  dir: %s\n  mode: dry-run\n"
        "molecule:\n  xyz: %s\n"
        % (run.resolve().as_posix(), xyz.resolve().as_posix()),
        encoding="utf-8")
    return cfg


def test_runner_logger_goes_to_stderr(capsys, workdir):
    config = load_config(_config(workdir))
    runner = _make_runner(argparse.Namespace(quiet=False), config)
    runner._log("[molecule_opt] prepare")
    out, err = capsys.readouterr()
    assert out == ""
    assert "[molecule_opt] prepare" in err


def test_unexpected_exception_is_reported_without_traceback(
        monkeypatch, capsys, workdir):
    cfg = _config(workdir)
    monkeypatch.setattr("smj.cli.cmd_report",
                        lambda args: (_ for _ in ()).throw(RuntimeError("boom")))
    assert main(["report", "-c", str(cfg)]) == 1
    out, err = capsys.readouterr()
    assert "Traceback" not in err
    assert "RuntimeError: boom" in err


def test_unexpected_exception_json_mode_is_parseable(
        monkeypatch, capsys, workdir):
    cfg = _config(workdir)
    monkeypatch.setattr("smj.cli.cmd_report",
                        lambda args: (_ for _ in ()).throw(RuntimeError("boom")))
    assert main(["report", "-c", str(cfg), "--json"]) == 1
    out, err = capsys.readouterr()
    # stdout 必须整体可解析（_emit_json 是缩进 JSON，没有进度日志混入）。
    payload = json.loads(out)
    assert payload["ok"] is False
    assert "RuntimeError: boom" in payload["error"]
    assert "Traceback" not in err


def test_report_strict_exits_nonzero_when_judgment_failed(
        monkeypatch, capsys, workdir):
    cfg = _config(workdir)
    fake = {
        "ok": True,
        "manifest": {"verdict": {"physics_verdict": "failed",
                                 "wording": "否证 H0", "errors": 1}},
    }
    monkeypatch.setattr("smj.cli._load_config", lambda args: cfg)
    monkeypatch.setattr("smj.cli._optional_module",
                        lambda name, purpose: type("M", (), {"generate_report": staticmethod(lambda *a, **k: fake)}))
    assert main(["report", "-c", str(cfg)]) == 0
    assert main(["report", "-c", str(cfg), "--strict"]) == 1


def test_review_stage_is_forwarded(monkeypatch, capsys, workdir):
    cfg = _config(workdir)
    stage = workdir / "electrode_build"
    stage.mkdir()
    seen = {}

    import smj.review as review_mod
    monkeypatch.setattr(
        review_mod, "review_run",
        lambda target, **kwargs: (seen.update(kwargs), {"ok": True})[1])
    assert main(["review", "-c", str(cfg), "--stage", str(stage)]) == 0
    assert seen.get("stage") == str(stage)

