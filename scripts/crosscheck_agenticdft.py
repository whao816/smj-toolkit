"""把同一份 TE.dat 分别交给本仓内核与 AgenticDFT 主仓算，逐项对比。

A 方案之后两边的**代码**不再同步，唯一一致的是**结果口径**。这个脚本就是那条
契约的可执行检查：手动跑、不进 CI（CI 里没有主仓）。

用法::

    python scripts/crosscheck_agenticdft.py demo/azulene-26-dithiol/.../TE.dat
    python scripts/crosscheck_agenticdft.py <TE.dat> --upstream D:/dsh_projects/AgenticDFT

退出码：两边都有实现且差异在容差内 → 0；差异超容差或主仓不可用 → 1。
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_UPSTREAM = ROOT.parent / "AgenticDFT"

#: 两边必须一致的字段与容差（相对差或绝对差）。
TOLERANCE = {
    "L0": ("rel", 1e-9),
    "S_uV_per_K": ("abs", 1e-6),
    "k_e_pW_per_K": ("rel", 1e-9),
    "ZT_electronic": ("rel", 1e-9),
}


def _load(module_name: str, path: Path):
    """按文件路径加载一个 landauer 模块（两边同名，不能靠 import 区分）。"""
    import importlib.util

    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise SystemExit("无法加载 %s" % path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _own(te: Path) -> dict:
    from smj.junction.landauer import analyze_te

    return analyze_te(te)


def _upstream(te: Path, upstream: Path) -> dict:
    module_path = upstream / "agenticdft" / "agenticdft" / "tools" / "junction" / "landauer.py"
    if not module_path.is_file():
        raise FileNotFoundError(str(module_path))
    return _load("agenticdft_landauer_crosscheck", module_path).analyze_te(te)


def _diff(name: str, mine, theirs) -> tuple[bool, str]:
    if mine is None or theirs is None:
        return mine == theirs, "%s: %r vs %r" % (name, mine, theirs)
    kind, limit = TOLERANCE[name]
    delta = abs(float(mine) - float(theirs))
    if kind == "rel":
        base = max(abs(float(mine)), abs(float(theirs)), 1e-30)
        ok = delta / base <= limit
        return ok, "%s: %.12g vs %.12g（相对差 %.3g）" % (name, mine, theirs, delta / base)
    ok = delta <= limit
    return ok, "%s: %.12g vs %.12g（绝对差 %.3g）" % (name, mine, theirs, delta)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="本仓内核 vs AgenticDFT 主仓的 TE.dat 结果对账")
    parser.add_argument("te", type=Path, help="TE.dat 路径")
    parser.add_argument("--upstream", type=Path, default=DEFAULT_UPSTREAM, help="AgenticDFT 主仓根目录")
    args = parser.parse_args(argv)

    if not args.te.is_file():
        print("找不到 TE.dat: %s" % args.te)
        return 1

    mine = _own(args.te)
    print("[本仓] %s" % args.te)
    print("  S=%s uV/K  L0=%s  k_e=%s pW/K  ZT=%s"
          % (mine.get("S_uV_per_K"), mine.get("L0"), mine.get("k_e_pW_per_K"), mine.get("ZT_electronic")))

    try:
        theirs = _upstream(args.te, args.upstream)
    except Exception as exc:  # noqa: BLE001 - 主仓缺失是正常情况，不当作崩溃
        print("[主仓] 不可用：%s" % exc)
        print("结论：无法对账（本仓结果已打印）")
        return 1

    missing = [name for name in TOLERANCE if name not in theirs]
    if missing:
        print("[主仓] 缺字段：%s（版本过旧？）" % ", ".join(missing))
        return 1

    bad = []
    for name in TOLERANCE:
        ok, line = _diff(name, mine.get(name), theirs.get(name))
        print("  " + ("OK   " if ok else "DIFF ") + line)
        if not ok:
            bad.append(name)

    print("结论：%s" % ("两边一致" if not bad else "不一致字段: %s" % ", ".join(bad)))
    return 0 if not bad else 1


if __name__ == "__main__":
    raise SystemExit(main())
