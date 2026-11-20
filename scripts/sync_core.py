#!/usr/bin/env python3
"""把 canonical ``AgenticDFT`` junction 同步为本包 ``src/smj/_vendor/junction`` 副本。

为什么要副本
------------
本包不能依赖 AgenticDFT 主仓（100+ 模块的研究智能体），但 junction 的算法
（几何/电极/输入写出器/自能/Landauer/判据）必须与主仓**同源**，否则两条实现
会漂移。因此：

- 主仓 ``agenticdft/agenticdft/tools/junction/`` 是**唯一开发位置**（canonical）；
- 本包 ``src/smj/_vendor/junction/`` 是它的副本 + 仅 3 处 import 改写
  （``preflight.py`` 里对 ``...backends`` / ``..resources`` 的惰性引用改为
  ``smj.compat``，见 :mod:`smj.compat`）；
- ``VENDOR.json`` 记录来源、``repo_commit``、每个文件的 sha256 与本地补丁台账；
- ``python scripts/sync_core.py --check`` 验证"副本 = 主仓当前内容"。

与主仓同名脚本的两点差异（2026-11-20，复核 R10）
------------------------------------------------
1. **布局**：副本目标是 ``src/smj/_vendor/junction``（主仓的等价脚本写
   ``smj/smj/_vendor/junction``，那是旧布局）。
2. **补丁感知**：默认**不覆写** ``VENDOR.json`` 里 ``local_patches`` 登记过的
   文件（``criteria.py`` / ``validate.py`` 等 4 类本地补丁），避免全量 sync
   静默抹掉它们。要显式覆盖受管文件必须 ``--files`` 指定，脚本会给出警告。
3. ``--check`` 不比对 ``repo_commit``：该字段是"最近一次同步时主仓的 HEAD"，
   主仓每提交一次都会变，拿它当漂移信号会让检查每次必红。文件哈希才是漂移信号。

用法
----
    python scripts/sync_core.py                     # 同步非受管文件（幂等）
    python scripts/sync_core.py --check             # 只检查，漂移则退出码 1
    python scripts/sync_core.py --files landauer.py,independent_check.py
    python scripts/sync_core.py --source D:/path/to/junction
    python scripts/sync_core.py --list              # 列出受管文件与补丁状态
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

SKIP_DIRS = {"__pycache__", ".pytest_cache"}
SKIP_SUFFIXES = {".pyc", ".pyo"}

#: 只允许这 3 行跨包 import 被改写（多一处都说明主仓结构变了，必须重新审查）。
REWRITES: tuple[tuple[str, str], ...] = (
    ("from ...backends import get_backend",
     "from smj.compat import get_backend"),
    ("from ..resources import estimate_smj_resources, check_smj_budget",
     "from smj.compat import estimate_smj_resources, check_smj_budget"),
    ("from ..resources import walltime_seconds",
     "from smj.compat import walltime_seconds"),
)


def default_source() -> Path:
    """canonical junction 目录（``SMJ_CANONICAL_JUNCTION`` 可覆盖）。"""
    env = os.environ.get("SMJ_CANONICAL_JUNCTION")
    if env:
        return Path(env)
    return ROOT.parent / "AgenticDFT" / "agenticdft" / "agenticdft" / "tools" / "junction"


def default_dest() -> Path:
    return ROOT / "src" / "smj" / "_vendor" / "junction"


def _git_head(repo: Path) -> str | None:
    try:
        proc = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"],
                              capture_output=True, text=True, timeout=15)
    except (OSError, subprocess.TimeoutExpired):
        return None
    return proc.stdout.strip() if proc.returncode == 0 else None


def _files(root: Path) -> list[Path]:
    out: list[Path] = []
    if not root.is_dir():
        return out
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        if any(part in SKIP_DIRS for part in path.parts):
            continue
        if path.suffix in SKIP_SUFFIXES:
            continue
        out.append(path)
    return out


def _rel(path: Path, root: Path) -> str:
    return path.relative_to(root).as_posix()


def _rewrite(rel: str, text: str) -> str:
    """对 ``preflight.py`` 做 3 处受控改写；每处必须恰好命中一次。"""
    if rel != "preflight.py":
        return text
    for old, new in REWRITES:
        count = text.count(old)
        if count != 1:
            raise SystemExit(
                "错误: preflight.py 中预期的跨包 import %r 命中 %d 次（应为 1）。\n"
                "主仓结构可能已变化：请重新审查 junction 的跨包依赖，更新 REWRITES。"
                % (old, count))
        text = text.replace(old, new)
    return text


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_meta(dest: Path) -> dict:
    manifest = dest / "VENDOR.json"
    if not manifest.is_file():
        return {}
    try:
        return json.loads(manifest.read_text(encoding="utf-8"))
    except ValueError:
        return {}


def patched_files(meta: dict) -> dict[str, list[str]]:
    """``{相对路径: [补丁 id, ...]}`` —— 受本地补丁管理的文件。"""
    out: dict[str, list[str]] = {}
    for patch in meta.get("local_patches") or []:
        for rel in patch.get("files") or []:
            out.setdefault(rel, []).append(str(patch.get("id", "?")))
    return out


def managed_files(source: Path) -> list[str]:
    """源目录里会被同步的受管文件（``*.py`` 与 ``templates/*``）。"""
    keep: list[str] = []
    for path in _files(source):
        rel = _rel(path, source)
        if rel.endswith(".py") or rel.startswith("templates/"):
            keep.append(rel)
    return keep


def build_plan(source: Path, only: set[str] | None = None) -> dict[str, bytes]:
    """``{相对路径: 目标字节}``（含改写后的 preflight.py 与模板）。"""
    if not source.is_dir():
        raise SystemExit(
            "找不到 canonical junction 目录: %s\n"
            "请用 --source 指定，或设置环境变量 SMJ_CANONICAL_JUNCTION。" % source)
    allowed = set(managed_files(source))
    if only:
        unknown = sorted(only - allowed)
        if unknown:
            raise SystemExit("--files 里有不在 canonical 目录的条目: %s" % ", ".join(unknown))
    plan: dict[str, bytes] = {}
    for path in _files(source):
        rel = _rel(path, source)
        if rel not in allowed:
            continue
        if only and rel not in only:
            continue
        raw = path.read_text(encoding="utf-8")
        plan[rel] = _rewrite(rel, raw).encode("utf-8")
    if not plan:
        raise SystemExit("没有可同步的文件（--files 过滤后为空）")
    return plan


def build_manifest(dest: Path, meta: dict, repo_commit: str | None) -> dict:
    """以磁盘实际内容为准重算 ``files``；其余台账字段原样保留。"""
    files: dict[str, dict] = {}
    for path in _files(dest):
        rel = _rel(path, dest)
        if rel == "VENDOR.json":
            continue
        data = path.read_bytes()
        files[rel] = {"sha256": _sha(data), "bytes": len(data)}
    out = {
        "schema": meta.get("schema", "smj.vendor.v1"),
        "upstream": meta.get("upstream", "agenticdft/agenticdft/tools/junction"),
        "upstream_commit": meta.get("upstream_commit"),
        "rebuilt_from": meta.get("rebuilt_from", ""),
        "repo_commit": repo_commit,
        "files": files,
        "local_patches": meta.get("local_patches") or [],
        "note": meta.get("note", ""),
    }
    return out


def _manifest_text(manifest: dict) -> str:
    return json.dumps(manifest, ensure_ascii=False, indent=2) + "\n"


def check(dest: Path, source: Path, only: set[str] | None) -> list[str]:
    """比对副本与 canonical 目录；带本地补丁的文件不参与漂移判定。

    带补丁文件本来就应该与主仓不同（补丁尚未回迁），把它们算成漂移会让 CI 永久红，
    从而失去信号。``--files`` 显式点名时仍按字节比对。
    """
    problems: list[str] = []
    if not source.is_dir():
        return ["canonical 目录不存在: %s" % source]
    meta = load_meta(dest)
    patched = patched_files(meta)
    plan = build_plan(source, only)
    if only:
        pass
    else:
        expected_diff = {rel: data for rel, data in plan.items() if rel in patched}
        plan = {rel: data for rel, data in plan.items() if rel not in patched}
        if expected_diff:
            print("跳过本地补丁文件（%d 个，不参与漂移判定）: %s"
                  % (len(expected_diff), ", ".join(sorted(expected_diff))))
    for rel, data in sorted(plan.items()):
        target = dest / rel
        if not target.is_file():
            problems.append("缺失: %s" % rel)
        elif target.read_bytes() != data:
            problems.append("内容漂移: %s" % rel)
    # VENDOR.json 的 files 哈希必须与磁盘一致（覆盖"改了 vendored 文件没重算"的情况）
    recorded = meta.get("files") or {}
    on_disk: dict[str, str] = {}
    for path in _files(dest):
        rel = _rel(path, dest)
        if rel == "VENDOR.json":
            continue
        on_disk[rel] = _sha(path.read_bytes())
    for rel in sorted(set(recorded) | set(on_disk)):
        old = (recorded.get(rel) or {}).get("sha256")
        if old != on_disk.get(rel):
            problems.append("VENDOR.json 哈希过期: %s（请重跑 gen_vendor_json.py）" % rel)
    return problems


def sync(dest: Path, source: Path, only: set[str] | None = None) -> int:
    meta = load_meta(dest)
    patched = patched_files(meta)
    plan = build_plan(source, only)
    if not only:
        skipped = sorted(set(managed_files(source)) - set(plan))
        if skipped:
            print("跳过本地补丁受管文件（%d 个）: %s" % (len(skipped), ", ".join(skipped)))
    for rel in sorted(plan):
        if rel in patched and only:
            print("警告: %s 带本地补丁 %s，按 --files 显式要求覆写" % (rel, "/".join(patched[rel])))
    dest.mkdir(parents=True, exist_ok=True)
    for rel, data in sorted(plan.items()):
        target = dest / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
    manifest = build_manifest(dest, meta, _git_head(ROOT))
    (dest / "VENDOR.json").write_text(_manifest_text(manifest), encoding="utf-8", newline="\n")
    # 让 _vendor/ 与 _vendor/junction/ 成为可导入包（不覆盖来源文件）。
    for pkg in (dest.parent, dest):
        init = pkg / "__init__.py"
        if not init.exists():
            init.write_text('"""vendored 代码（由 scripts/sync_core.py 生成，勿手改）。"""\n',
                            encoding="utf-8", newline="\n")
    print("已同步 %d 个文件 -> %s（VENDOR.json 已重算）" % (len(plan), dest))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="同步 canonical junction 到本包 vendored 副本")
    parser.add_argument("--source", type=Path, default=None, help="canonical junction 目录")
    parser.add_argument("--dest", type=Path, default=None, help="副本目标目录")
    parser.add_argument("--files", default=None,
                        help="逗号分隔的相对路径；缺省=全部非受管文件")
    parser.add_argument("--check", action="store_true", help="只检查，不写文件")
    parser.add_argument("--list", action="store_true", help="列出受管文件与补丁状态后退出")
    args = parser.parse_args(argv)

    source = (args.source or default_source()).resolve()
    dest = (args.dest or default_dest()).resolve()
    only = None
    if args.files:
        only = {item.strip() for item in args.files.split(",") if item.strip()}

    if args.list:
        meta = load_meta(dest)
        patched = patched_files(meta)
        if not source.is_dir():
            print("canonical 目录不存在: %s" % source)
            return 2
        for rel in managed_files(source):
            mark = "补丁:" + "/".join(patched[rel]) if rel in patched else "同步"
            print("%-8s %s" % (mark, rel))
        return 0

    if args.check:
        problems = check(dest, source, only)
        if problems:
            print("smj vendor 副本与主仓不一致（%d 处）:" % len(problems))
            for item in problems:
                print("  - " + item)
            print("\n修复: python scripts/sync_core.py [--files ...] && python scripts/gen_vendor_json.py")
            return 1
        print("smj vendor 副本与主仓一致（%d 个文件）" % len(build_plan(source, only)))
        return 0

    return sync(dest, source, only)


if __name__ == "__main__":
    sys.exit(main())
