"""重新生成 ``src/smj/_vendor/junction/VENDOR.json``（本地补丁台账）。

为什么需要它：0.1.0 的 wheel **没有**带 VENDOR.json（README 却说它记录 sha256 与
主仓 commit），而 0.1.1 在 vendored 目录里做了 4 处本地补丁（F2 判据、F1 头回取、
F4 提示文案），必须留下可追溯的记录，否则下次从主仓 sync 会静默覆盖。

本脚本只用标准库，输出确定性内容（文件按路径排序；不写时间戳到 files 之外时
可选择 --no-timestamp 便于 diff）。用法::

    python scripts/gen_vendor_json.py           # 写入 VENDOR.json
    python scripts/gen_vendor_json.py --check   # 只比对，不一致返回 1
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VENDOR_DIR = ROOT / "src" / "smj" / "_vendor" / "junction"
VENDOR_JSON = VENDOR_DIR / "VENDOR.json"

# 0.1.1 的本地补丁台账：每条写明改了什么、为什么、被哪条测试守住。
LOCAL_PATCHES = [
    {
        "id": "F2-te-present",
        "commit": "bc7febf",
        "files": ["criteria.py", "validate.py"],
        "what": "无 TE.dat 时登记 te_present(engineering/error, missing=True)，"
                "physics_verdict([]) 返回 inconclusive 而非通过；validate_junction_run "
                "接受 temperature 参数。",
        "why": "原实现里空 aitranss_te/ 目录会得到 passed=True / 支持 H0，"
               "缺证据被当成通过（P0 fail-open）。",
        "tests": ["tests/test_criteria_failopen.py"],
    },
    {
        "id": "F1-mos-aims-head",
        "commit": "f7d3970",
        "files": ["artifacts.py"],
        "what": "新增 AIMS_HEAD_STAGES=('aims_scf',)，只有 aims_scf 的 retrieval_plan "
                "才带 mos.aims/omat.aims 的 head 条目。",
        "why": "原先所有阶段都收到 head 条目，gaussian_opt 等阶段凭空多出必然 missing "
               "的项；真正消费 head 的下载逻辑在 smj/remote.py（非 vendored）。",
        "tests": ["tests/test_remote_mos_head.py"],
    },
    {
        "id": "F4-pending-message",
        "commit": "fe48c29",
        "files": ["review.py"],
        "what": "pending_message 不再提示主仓的 `agenticdft smj-review`，改为说明独立包"
                "的真实入口 `smj review --approve --who <名字>`。",
        "why": "该提示指向独立包里不存在的命令（且该检查点模块本身未接线）。",
        "tests": ["tests/test_a2_approve.py"],
    },
]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git_head() -> str | None:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                             capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.SubprocessError):  # pragma: no cover
        return None
    return out.stdout.strip() or None if out.returncode == 0 else None


def build_manifest() -> dict:
    files = {}
    for path in sorted(VENDOR_DIR.rglob("*")):
        if not path.is_file() or path.name == "VENDOR.json":
            continue
        if "__pycache__" in path.parts:
            continue
        rel = path.relative_to(VENDOR_DIR).as_posix()
        files[rel] = {"sha256": _sha256(path), "bytes": path.stat().st_size}
    return {
        "schema": "smj.vendor.v1",
        "upstream": "agenticdft/agenticdft/tools/junction",
        "upstream_commit": None,
        "rebuilt_from": "smj_toolkit-0.1.0-py3-none-any.whl（原 wheel 未带 VENDOR.json）",
        "repo_commit": _git_head(),
        "files": files,
        "local_patches": LOCAL_PATCHES,
        "note": "0.1.1 本地补丁的权威记录；回迁主仓后应删除对应条目并重新 sync。"
                "本文件由 scripts/gen_vendor_json.py 生成，勿手改。",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true",
                        help="只比对现有 VENDOR.json，不一致返回 1（不写文件）")
    args = parser.parse_args(argv)

    manifest = build_manifest()
    text = json.dumps(manifest, ensure_ascii=False, indent=2) + "\n"
    if args.check:
        current = VENDOR_JSON.read_text(encoding="utf-8") if VENDOR_JSON.is_file() else ""
        if current != text:
            print("VENDOR.json 与磁盘内容不一致：请运行 python scripts/gen_vendor_json.py")
            return 1
        print("VENDOR.json 一致：%d 个文件" % len(manifest["files"]))
        return 0

    VENDOR_JSON.write_text(text, encoding="utf-8", newline="\n")
    print("写入 %s（%d 个文件，%d 条本地补丁）"
          % (VENDOR_JSON, len(manifest["files"]), len(manifest["local_patches"])))
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
