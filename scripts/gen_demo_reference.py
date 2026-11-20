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

``--check`` 比较的是"参考值文件内容 == 由当前代码重算出的内容"，
因此它同时守住两件事：谱文件被改动（te_sha256 变）与积分口径被改动
（S/L0/ZT 变）。故意改参考值时请跑不带 --check 的形式并提交新文件。
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

from smj._vendor.junction.landauer import analyze_te  # noqa: E402

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


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="只校验，不一致返回 1")
    parser.add_argument("--output", default=str(REFERENCE), help="输出路径")
    args = parser.parse_args(argv)

    target = Path(args.output)
    text = render(build_reference())
    if args.check:
        if not target.is_file():
            print("缺少参考值文件: %s" % target, file=sys.stderr)
            return 1
        current = target.read_text(encoding="utf-8")
        if current != text:
            print(
                "demo 参考值与当前代码/谱文件不一致: %s\n"
                "（谱被改动，或积分口径变了。确认是有意修改后运行 "
                "`python scripts/gen_demo_reference.py` 并提交新文件）" % target,
                file=sys.stderr,
            )
            return 1
        print("demo 参考值一致（%d 份谱）" % len(te_files()))
        return 0

    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text, encoding="utf-8")
    print("写入 %s（%d 份谱）" % (target, len(te_files())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
