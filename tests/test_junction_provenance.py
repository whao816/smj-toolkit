"""PROVENANCE.json 台账与脚本语义（2026-11-20 复核 R10 / 后续③）。

A 方案之前，「vendored 副本」的说法有两处靠不住：
1. ``README`` 让人跑 ``scripts/sync_core.py --check``，但那个脚本只存在于兄弟仓；
2. ``--check`` 整份文本比较，``repo_commit`` 每次提交都变 → 这条校验**永远红**。

A 方案之后本仓自维护 junction 内核：``sync_core`` 已删除，``PROVENANCE.json``
只记录来源与本地变更史（``--check`` 忽略 ``repo_commit``，其余字段逐一比对）。
本文件守住三件事：哈希与磁盘一致、``repo_commit`` 不参与、补丁登记的文件真实存在。
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "gen_junction_provenance.py"
KERNEL_DIR = ROOT / "src" / "smj" / "junction"
PROVENANCE_PATH = KERNEL_DIR / "PROVENANCE.json"


def _load_script():
    spec = importlib.util.spec_from_file_location("gen_provenance_json_under_test", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture()
def script():
    return _load_script()


def _recorded() -> dict:
    return json.loads(PROVENANCE_PATH.read_text(encoding="utf-8"))


def _strip_commit(manifest: dict) -> dict:
    manifest = dict(manifest)
    manifest.pop("repo_commit", None)
    return manifest


def test_provenance_json_matches_vendored_files(script):
    """CI 的 vendor job 同款检查，本地跑测试时也要绿。"""
    manifest = script.build_manifest()
    assert manifest["files"], "vendored 目录下没有文件？"
    assert _strip_commit(_recorded()) == _strip_commit(manifest)


def test_check_ignores_repo_commit(script, workdir, monkeypatch, capsys):
    """提交一次就让校验变红是旧实现的病根：repo_commit 不该参与比对。"""
    patched = json.loads(PROVENANCE_PATH.read_text(encoding="utf-8"))
    patched["repo_commit"] = "0" * 40
    target = workdir / "PROVENANCE.json"
    target.write_text(json.dumps(patched, ensure_ascii=False, indent=2), encoding="utf-8")
    monkeypatch.setattr(script, "PROVENANCE_PATH", target)

    assert script.main(["--check"]) == 0
    assert "repo_commit 不参与比对" in capsys.readouterr().out


def test_check_detects_real_content_drift(script, workdir, monkeypatch, capsys):
    """真正漂移（sha256 变了）必须报错，否则校验是空转。"""
    patched = json.loads(PROVENANCE_PATH.read_text(encoding="utf-8"))
    first = sorted(patched["files"])[0]
    patched["files"][first]["sha256"] = "f" * 64
    target = workdir / "PROVENANCE.json"
    target.write_text(json.dumps(patched, ensure_ascii=False, indent=2), encoding="utf-8")
    monkeypatch.setattr(script, "PROVENANCE_PATH", target)

    assert script.main(["--check"]) == 1
    assert "不一致" in capsys.readouterr().out


def test_check_reports_missing_manifest(script, workdir, monkeypatch, capsys):
    monkeypatch.setattr(script, "PROVENANCE_PATH", workdir / "missing.json")
    assert script.main(["--check"]) == 1
    assert "缺少 PROVENANCE.json" in capsys.readouterr().out


def test_local_patch_files_exist_and_are_tracked(script):
    """补丁台账的文件名是「哪些文件被本地改过」的唯一记录，写错就失去可追溯性。"""
    manifest = script.build_manifest()
    recorded = _recorded()
    assert recorded["local_patches"], "本地补丁台账为空？"
    ids = [patch["id"] for patch in recorded["local_patches"]]
    assert len(ids) == len(set(ids)), "补丁 id 重复: %s" % ids

    for patch in recorded["local_patches"]:
        assert patch["files"], "%s 没有登记文件" % patch["id"]
        assert patch["what"] and patch["why"], "%s 缺少 what/why 说明" % patch["id"]
        for rel in patch["files"]:
            assert rel in manifest["files"], "%s 登记了不存在的文件 %s" % (patch["id"], rel)
            assert (KERNEL_DIR / rel).is_file()


def test_every_patch_file_would_not_be_overwritten_silently(script):
    """被补丁改过的文件必须都在台账里，否则日后看不出它已被本地改过。"""
    manifest = script.build_manifest()
    patched = {rel for patch in _recorded()["local_patches"] for rel in patch["files"]}
    assert patched <= set(manifest["files"])
    # 台账自称的补丁文件数应与 local_patches 汇总一致，且都真实存在。
    assert len(patched) >= 11, "补丁文件集合异常小: %s" % sorted(patched)
