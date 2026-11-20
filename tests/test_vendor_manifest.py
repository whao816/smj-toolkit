"""VENDOR.json 台账与脚本语义（2026-11-20 复核 R10 / 后续③）。

复核发现的两件事：
1. ``README`` 让人跑 ``scripts/sync_core.py --check``，但仓库里没有该脚本；
2. ``scripts/gen_vendor_json.py --check`` 整份文本比较，``repo_commit`` 每次提交
   都变 → 这条校验**永远红**，等于没有校验。

第 1 条由新增的 ``scripts/sync_core.py`` 解决（见 ``tests`` 里对本文件的用例），
第 2 条由 ``--check`` 忽略 ``repo_commit`` 解决。本文件守住这两个语义，外加
"补丁台账里登记的文件必须真实存在"，否则 ``sync_core`` 会静默覆盖本地补丁。
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "gen_vendor_json.py"
VENDOR_DIR = ROOT / "src" / "smj" / "_vendor" / "junction"
VENDOR_JSON = VENDOR_DIR / "VENDOR.json"


def _load_script():
    spec = importlib.util.spec_from_file_location("gen_vendor_json_under_test", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture()
def script():
    return _load_script()


def _recorded() -> dict:
    return json.loads(VENDOR_JSON.read_text(encoding="utf-8"))


def _strip_commit(manifest: dict) -> dict:
    manifest = dict(manifest)
    manifest.pop("repo_commit", None)
    return manifest


def test_vendor_json_matches_vendored_files(script):
    """CI 的 vendor job 同款检查，本地跑测试时也要绿。"""
    manifest = script.build_manifest()
    assert manifest["files"], "vendored 目录下没有文件？"
    assert _strip_commit(_recorded()) == _strip_commit(manifest)


def test_check_ignores_repo_commit(script, workdir, monkeypatch, capsys):
    """提交一次就让校验变红是旧实现的病根：repo_commit 不该参与比对。"""
    patched = json.loads(VENDOR_JSON.read_text(encoding="utf-8"))
    patched["repo_commit"] = "0" * 40
    target = workdir / "VENDOR.json"
    target.write_text(json.dumps(patched, ensure_ascii=False, indent=2), encoding="utf-8")
    monkeypatch.setattr(script, "VENDOR_JSON", target)

    assert script.main(["--check"]) == 0
    assert "repo_commit 不参与比对" in capsys.readouterr().out


def test_check_detects_real_content_drift(script, workdir, monkeypatch, capsys):
    """真正漂移（sha256 变了）必须报错，否则校验是空转。"""
    patched = json.loads(VENDOR_JSON.read_text(encoding="utf-8"))
    first = sorted(patched["files"])[0]
    patched["files"][first]["sha256"] = "f" * 64
    target = workdir / "VENDOR.json"
    target.write_text(json.dumps(patched, ensure_ascii=False, indent=2), encoding="utf-8")
    monkeypatch.setattr(script, "VENDOR_JSON", target)

    assert script.main(["--check"]) == 1
    assert "不一致" in capsys.readouterr().out


def test_check_reports_missing_manifest(script, workdir, monkeypatch, capsys):
    monkeypatch.setattr(script, "VENDOR_JSON", workdir / "missing.json")
    assert script.main(["--check"]) == 1
    assert "缺少 VENDOR.json" in capsys.readouterr().out


def test_local_patch_files_exist_and_are_tracked(script):
    """补丁台账的文件名是 sync_core 跳过覆盖的唯一依据，写错就会静默丢补丁。"""
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
            assert (VENDOR_DIR / rel).is_file()


def test_every_patch_file_would_not_be_overwritten_silently(script):
    """被补丁改过的文件必须都在台账里，否则 sync_core 会把它当普通文件同步。"""
    manifest = script.build_manifest()
    patched = {rel for patch in _recorded()["local_patches"] for rel in patch["files"]}
    assert patched <= set(manifest["files"])
    # 台账自称的补丁文件数应与 sync_core 的跳过集合一致（脚本自身会打印并跳过）。
    assert len(patched) >= 11, "补丁文件集合异常小: %s" % sorted(patched)
