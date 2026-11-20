#!/usr/bin/env python
"""生成/校验 demo 回归参考值（2026-11-20 复核 R9/后续①）。

demo/ 下的 TE.dat 是唯一随发布走的真实谱。以前它们只被"能不能读"覆盖，
没有任何数值回归：改了 landauer.py 的积分、单位或温度口径，CI 也不会红。
本脚本把每份谱的关键量算出来固化成 demo/demo-reference.json，测试
``tests/test_demo_regression.py`` 再用 ``analysis.compare_reference`` 对它做
容差比对。

用法::

    python scripts/gen_demo_reference.py            # 重写 demo/demo-reference.json
    python scripts/gen_demo_reference.py --check    # 只校验，不一致返回 1

``--check`` 守住的两件事是：谱文件被改动（te_sha256 变）与积分口径被改动
（S/L0/ZT 变）。它按字段比较而不是比较整份文本（2026-11-20 修复）：
``L0``/``ZT`` 这类量是浮点数的完整精度，``np.exp``/``np.cosh`` 在不同平台
（Windows 与 Linux 的 numpy 构建、不同 CPU 的 SIMD 分派）上可以有 1 ulp
级差异，逐字节比对会让 CI 在"什么都没改"的情况下变红。这里对数值字段用
``REL_TOL`` 比较，对哈希/点数/能量顺序等用精确比较。

故意改参考值时请跑不带 ``--check`` 的形式并提交新文件。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from smj.junction.landauer import analyze_te  # noqa: E402

DEMO_DIR = ROOT / "demo"
REFERENCE = DEMO_DIR / "demo-reference.json"
SCHEMA = "smj.demo-reference.v1"

#: 判据容差沿用 ``smj.analysis.DEFAULT_TOLERANCES`` 的口径：
#: log10(G/G0) 绝对容差 0.1；S / ZT / k_e 相对容差 10%。谱是同一份文件，
#: 容差放宽主要是给不同平台的浮点求和顺序留余量。
TOLERANCE = {
    "log10_G_over_G0_like": {"abs": 0.1},
    "S_uV_per_K": {"rel": 0.10},
    "ZT_electronic": {"rel": 0.10},
    "k_e_pW_per_K": {"rel": 0.10},
}

NOTE = (
    "demo 谱的数值回归基线。由 scripts/gen_demo_reference.py 生成；"
    "改动谱文件或积分口径后请重新生成并在评审里说明原因。"
)

#: ``--check`` 的数值比较口径。1e-9 足以抓住任何真实的积分/单位改动
#: （那类差异从 1e-3 量级起），同时容忍跨平台 numpy 的 1 ulp 级浮点差异。
REL_TOL = 1e-9
ABS_TOL = 1e-15

#: 必须精确相等的字段：含哈希与计数，任何变化都意味着谱被换掉。
_EXACT_KEYS = ("te_sha256", "points", "input_energy_order")
#: 按相对容差比较的浮点字段（写入时并不做四舍五入，因此逐字节比不可能稳）。
_FLOAT_KEYS = (
    "log10_G_over_G0_like",
    "L0",
    "S_uV_per_K",
    "k_e_pW_per_K",
    "ZT_electronic",
    "max_transmission_T",
    "max_transmission_at_eV",
)


def te_files() -> list[Path]:
    return sorted(DEMO_DIR.glob("**/aitranss_te/TE.dat"))


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_reference() -> dict:
    demos: dict[str, dict] = {}
    for te in te_files():
        result = analyze_te(te)
        demos[te.relative_to(DEMO_DIR).as_posix()] = {
            "te_sha256": _digest(te),
            "points": int(result["points"]),
            "temperature_K": float(result["temperature_K"]),
            "log10_G_over_G0_like": round(float(result["log10_G_over_G0_like"]), 9),
            "L0": float(result["L0"]),
            "S_uV_per_K": round(float(result["S_uV_per_K"]), 9),
            "k_e_pW_per_K": round(float(result["k_e_pW_per_K"]), 9),
            "ZT_electronic": float(result["ZT_electronic"]),
            "max_transmission_T": float(result["max_transmission"]["T"]),
            "max_transmission_at_eV": float(result["max_transmission"]["E_minus_Ef_eV"]),
            "input_energy_order": str(result["input_energy_order"]),
        }
    return {"schema": SCHEMA, "note": NOTE, "tolerance": TOLERANCE, "demos": demos}


def render(manifest: dict) -> str:
    return json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n"


def _compare_demo(expected: dict, current: dict) -> tuple[list[str], float]:
    """比较一份谱的参考条目，返回 (问题列表, 最大相对偏差)。"""
    problems: list[str] = []
    worst = 0.0
    for field in _EXACT_KEYS:
        if current.get(field) != expected[field]:
            problems.append("%s: 期望 %r 实际 %r" % (field, expected[field], current.get(field)))
    try:
        got_temperature = float(current["temperature_K"])
    except (KeyError, TypeError, ValueError):
        problems.append("temperature_K: 缺失或不是数值 (%r)" % (current.get("temperature_K"),))
    else:
        if abs(got_temperature - float(expected["temperature_K"])) > 1e-9:
            problems.append(
                "temperature_K: 期望 %r 实际 %r" % (expected["temperature_K"], got_temperature)
            )
    for field in _FLOAT_KEYS:
        try:
            value = float(current[field])
        except (KeyError, TypeError, ValueError):
            problems.append("%s: 缺失或不是数值 (%r)" % (field, current.get(field)))
            continue
        want = float(expected[field])
        diff = abs(value - want)
        scale = max(abs(want), abs(value), 1e-30)
        rel = diff / scale
        worst = max(worst, rel)
        if diff > ABS_TOL + REL_TOL * abs(want):
            problems.append(
                "%s: 期望 %.17g 实际 %.17g（相对差 %.2e）" % (field, want, value, rel)
            )
    return problems, worst


def compare_reference(expected: dict, current: dict) -> tuple[list[str], float]:
    """按字段比较参考值；返回 (问题列表, 最大相对偏差)。"""
    problems: list[str] = []
    worst = 0.0
    if not isinstance(current, dict):
        return ["参考值不是 JSON 对象"], worst
    for field in ("schema", "note", "tolerance"):
        if current.get(field) != expected[field]:
            problems.append("%s: 期望 %r 实际 %r" % (field, expected[field], current.get(field)))
    demos = current.get("demos")
    if not isinstance(demos, dict):
        return problems + ["demos: 缺失或不是对象"], worst
    wanted = expected["demos"]
    missing = sorted(set(wanted) - set(demos))
    extra = sorted(set(demos) - set(wanted))
    if missing or extra:
        problems.append("demo 谱集合不一致: 缺 %s；多 %s" % (missing, extra))
    for key in sorted(set(wanted) & set(demos)):
        entry = demos[key]
        if not isinstance(entry, dict):
            problems.append("%s: 条目不是对象" % key)
            continue
        found, rel = _compare_demo(wanted[key], entry)
        worst = max(worst, rel)
        problems.extend("%s → %s" % (key.split("/")[0], item) for item in found)
    return problems, worst


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="只校验，不一致返回 1")
    parser.add_argument("--output", default=str(REFERENCE), help="输出路径")
    args = parser.parse_args(argv)

    target = Path(args.output)
    expected = build_reference()
    text = render(expected)
    if args.check:
        if not target.is_file():
            print("缺少参考值文件: %s" % target, file=sys.stderr)
            return 1
        try:
            current = json.loads(target.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            print("参考值不是合法 JSON: %s (%s)" % (target, exc), file=sys.stderr)
            return 1
        problems, worst = compare_reference(expected, current)
        if problems:
            print(
                "demo 参考值与当前代码/谱文件不一致: %s\n  %s\n"
                "（谱被改动，或积分口径变了。确认是有意修改后运行 "
                "`python scripts/gen_demo_reference.py` 并提交新文件）"
                % (target, "\n  ".join(problems)),
                file=sys.stderr,
            )
            return 1
        print(
            "demo 参考值一致（%d 份谱；数值按相对 %.0e 比对，最大相对偏差 %.2e）"
            % (len(te_files()), REL_TOL, worst)
        )
        return 0

    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text, encoding="utf-8", newline="\n")
    print("写入 %s（%d 份谱）" % (target, len(te_files())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
