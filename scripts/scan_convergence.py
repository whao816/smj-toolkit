#!/usr/bin/env python
"""收敛扫描：把「改一个参数重跑一遍」变成可复现的矩阵 + 对照表。

用法::

    # 1) 只出计划：按轴生成变体目录，每个目录一份独立 smj.yaml
    python scripts/scan_convergence.py plan smj.yaml \
        --axis electrode.nlayers=3,4,5 \
        --axis paths.species_dir=/data/species/light,/data/species/tight \
        --out scan/ --pyramid electrode.au_per_electrode=20,30

    # 2) 真跑（必须有 --yes；AGENTS.md：只跑用户点名的计算）
    python scripts/scan_convergence.py run scan/scan.json --yes

    # 3) 汇总：读各变体的 junction_analyze/analysis.json，给逐轴收敛判定
    python scripts/scan_convergence.py compare scan/scan.json

设计取舍：

- **不猜物理**。脚本只认 ``section.key=值1,值2`` 这种点路径覆盖，具体拿哪些轴
  收敛（层数/团簇尺寸/基组档位/能量窗）由调用者按体系决定。
- **基组档位就是路径**。FHI-aims 的 tier 不是配置键，而是 ``paths.species_dir``
  指向不同物种目录，所以轴值写成目录路径即可。
- **默认只写计划**。``plan`` 从不提交作业；``run`` 需要 ``--yes``，并且逐条打印
  实际执行的命令，方便先看一眼再决定。
- 变体的 ``smj.yaml`` 各自带 ``run.name`` 与 ``run.dir: run``，互不覆盖状态。
"""
from __future__ import annotations

import argparse
import itertools
import json
import subprocess
import sys
from pathlib import Path

import yaml

SCHEMA = "smj.convergence-scan.v1"
DEFAULT_ANALYSIS = "junction_analyze/analysis.json"
#: 判「收敛」的默认阈值：相邻两档相对变化 < 5%。
DEFAULT_TOL = 0.05
#: 对照表里关心的量（与 analysis.json / demo-reference.json 同名）。
METRICS = ("S_uV_per_K", "log10_G_over_G0_like", "ZT_electronic", "k_e_pW_per_K")


class ScanError(ValueError):
    """扫描规格不合法（轴解析失败、重名、缺基线）。"""


def typed(text: str) -> object:
    """把命令行里的字符串还原成 bool/int/float/None/str。"""
    raw = text.strip()
    low = raw.lower()
    if low in ("true", "false"):
        return low == "true"
    if low in ("none", "null"):
        return None
    for cast in (int, float):
        try:
            return cast(raw)
        except ValueError:
            continue
    return raw


def parse_axis(text: str) -> tuple[str, list[object]]:
    """``electrode.nlayers=3,4,5`` → ``("electrode.nlayers", [3, 4, 5])``。"""
    if "=" not in text:
        raise ScanError("轴必须写成 section.key=v1,v2,…，收到 %r" % text)
    key, values = text.split("=", 1)
    key = key.strip()
    if not key or "." not in key:
        raise ScanError("轴名必须是 section.key 形式，收到 %r" % key)
    parsed = [typed(part) for part in values.split(",") if part.strip()]
    if not parsed:
        raise ScanError("轴 %s 没有取值" % key)
    seen: list[object] = []
    for value in parsed:
        if value not in seen:
            seen.append(value)
    return key, seen


def set_path(data: dict, key: str, value: object) -> None:
    """按点路径写嵌套 dict（中间层缺失就建）。"""
    parts = key.split(".")
    node = data
    for part in parts[:-1]:
        child = node.get(part)
        if not isinstance(child, dict):
            child = {}
            node[part] = child
        node = child
    node[parts[-1]] = value


def get_path(data: dict, key: str) -> object:
    node: object = data
    for part in key.split("."):
        if not isinstance(node, dict) or part not in node:
            return None
        node = node[part]
    return node


def variant_name(index: int, values: dict[str, object]) -> str:
    """``v003__nlayers-5__species_dir-tight``：目录名自带轴值，肉眼可读。"""
    pieces = ["v%03d" % index]
    for key, value in values.items():
        leaf = key.split(".")[-1]
        shown = "None" if value is None else str(value)
        # 目录里塞完整路径太长，只留末段；完整值记在 scan.json 里。
        shown = shown.replace("\\", "/").rstrip("/").split("/")[-1] or "root"
        pieces.append("%s-%s" % (leaf, shown))
    return "__".join(pieces)


def build_plan(base_path: Path, axes: list[tuple[str, list[object]]],
               out_dir: Path, *, constants: dict[str, object] | None = None) -> dict:
    """把轴组合展开成变体清单（纯计算，不落盘）。"""
    if not axes:
        raise ScanError("至少要给一个 --axis")
    keys = [key for key, _ in axes]
    if len(set(keys)) != len(keys):
        raise ScanError("轴重复: %s" % keys)
    base = yaml.safe_load(base_path.read_text(encoding="utf-8")) or {}
    if not isinstance(base, dict):
        raise ScanError("%s 不是 mapping 配置" % base_path)
    variants: list[dict] = []
    for index, combo in enumerate(itertools.product(*[values for _, values in axes]), start=1):
        values = dict(zip(keys, combo))
        config = yaml.safe_load(yaml.safe_dump(base, allow_unicode=True)) or {}
        for key, value in (constants or {}).items():
            set_path(config, key, value)
        for key, value in values.items():
            set_path(config, key, value)
        config.setdefault("run", {})
        if not isinstance(config["run"], dict):
            raise ScanError("配置里的 run 段不是 mapping")
        config["run"]["name"] = variant_name(index, values)
        config["run"]["dir"] = "run"
        variants.append({"index": index, "name": config["run"]["name"],
                         "values": values, "config": config})
    return {"schema": SCHEMA, "base": str(base_path), "out_dir": str(out_dir),
            "axes": [{"key": key, "values": values} for key, values in axes],
            "constants": constants or {}, "variants": variants}


def write_plan(plan: dict) -> Path:
    out_dir = Path(plan["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    for variant in plan["variants"]:
        target = out_dir / variant["name"]
        target.mkdir(parents=True, exist_ok=True)
        (target / "smj.yaml").write_text(
            yaml.safe_dump(variant["config"], allow_unicode=True, sort_keys=False),
            encoding="utf-8", newline="\n")
    manifest = out_dir / "scan.json"
    public = dict(plan)
    public["variants"] = [
        {"index": v["index"], "name": v["name"], "values": v["values"],
         "dir": str(out_dir / v["name"]), "config": str(out_dir / v["name"] / "smj.yaml")}
        for v in plan["variants"]
    ]
    manifest.write_text(json.dumps(public, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8", newline="\n")
    return manifest


def load_manifest(path: Path) -> dict:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if data.get("schema") != SCHEMA:
        raise ScanError("%s 不是 %s 清单" % (path, SCHEMA))
    if not data.get("variants"):
        raise ScanError("%s 里没有变体" % path)
    return data


def variant_metrics(variant: dict) -> dict:
    """读一个变体的分析结果；缺文件就标 missing，不伪造 0。"""
    run_dir = Path(variant["dir"]) / "run"
    analysis = run_dir / DEFAULT_ANALYSIS
    if not analysis.is_file():
        return {"status": "missing", "run_dir": str(run_dir), "values": {}}
    payload = json.loads(analysis.read_text(encoding="utf-8"))
    values = {key: payload.get(key) for key in METRICS}
    return {"status": "ok", "run_dir": str(run_dir), "values": values}


def _relative(new: float, base: float) -> float:
    if base == 0:
        return float("inf") if new else 0.0
    return abs(new - base) / abs(base)


def compare_axis(axis: dict, rows: list[dict], tol: float) -> dict:
    """沿一条轴按取值顺序比较相邻两档，返回收敛判定。"""
    ordered = []
    for value in axis["values"]:
        for row in rows:
            if row["values"].get(axis["key"]) == value:
                ordered.append((value, row["metrics"]))
                break
    changes = []
    for (prev_value, prev), (next_value, nxt) in zip(ordered, ordered[1:]):
        if prev.get("status") != "ok" or nxt.get("status") != "ok":
            changes.append({"from": prev_value, "to": next_value, "status": "missing"})
            continue
        detail = {}
        worst = 0.0
        for metric in METRICS:
            before, after = prev["values"].get(metric), nxt["values"].get(metric)
            if before is None or after is None:
                continue
            rel = _relative(float(after), float(before))
            detail[metric] = rel
            worst = max(worst, rel)
        changes.append({"from": prev_value, "to": next_value,
                        "worst_relative_change": worst,
                        "converged": worst < tol, "detail": detail})
    if not changes:
        verdict = "insufficient"
    elif any(change.get("status") == "missing" for change in changes):
        verdict = "incomplete"
    elif all(change.get("converged") for change in changes):
        verdict = "converged"
    else:
        verdict = "not-converged"
    return {"axis": axis["key"], "tolerance": tol, "verdict": verdict,
            "changes": changes}


def table(plan: dict, rows: list[dict]) -> str:
    header = ["变体"] + [axis["key"] for axis in plan["axes"]] + ["状态"] + list(METRICS)
    lines = [" | ".join(header), " | ".join("---" for _ in header)]
    for row in rows:
        cells = [row["name"]]
        cells += ["%s" % row["values"].get(axis["key"]) for axis in plan["axes"]]
        cells.append(row["metrics"]["status"])
        for metric in METRICS:
            value = row["metrics"]["values"].get(metric)
            cells.append("-" if value is None else "%.6g" % value)
        lines.append(" | ".join(cells))
    return "\n".join(lines)


def cmd_plan(args: argparse.Namespace) -> int:
    axes = [parse_axis(text) for text in args.axis]
    constants = {}
    for item in args.set or []:
        key, _, raw = item.partition("=")
        constants[key.strip()] = typed(raw)
    plan = build_plan(Path(args.config), axes, Path(args.out), constants=constants)
    manifest = write_plan(plan)
    print("扫描计划: %s（%d 个变体）" % (manifest, len(plan["variants"])))
    for variant in plan["variants"]:
        print("  %-40s %s" % (variant["name"],
                              ", ".join("%s=%s" % kv for kv in variant["values"].items())))
    print("\n下一步: python scripts/scan_convergence.py run %s --yes" % manifest)
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    plan = load_manifest(Path(args.manifest))
    commands = []
    for variant in plan["variants"]:
        config_path = variant.get("config") or str(Path(variant["dir"]) / "smj.yaml")
        command = [sys.executable, "-m", "smj.cli", "run", "-c", config_path]
        if args.from_stage:
            command += ["--from", args.from_stage]
        if args.to_stage:
            command += ["--to", args.to_stage]
        commands.append(command)
    if not args.yes:
        print("以下命令尚未执行（加 --yes 才真跑；AGENTS.md：只跑用户点名的计算）：")
        for command in commands:
            print("  " + " ".join(command))
        return 0
    failed = 0
    for variant, command in zip(plan["variants"], commands):
        print("== %s ==" % variant["name"])
        result = subprocess.run(command)
        if result.returncode != 0:
            failed += 1
            print("  失败（退出码 %d），继续下一个变体" % result.returncode)
    print("\n完成: %d/%d 个变体成功" % (len(commands) - failed, len(commands)))
    return 0 if failed == 0 else 1


def cmd_compare(args: argparse.Namespace) -> int:
    plan = load_manifest(Path(args.manifest))
    rows = []
    for variant in plan["variants"]:
        rows.append({"name": variant["name"], "values": variant["values"],
                     "metrics": variant_metrics(variant)})
    print(table(plan, rows))
    print()
    for axis in plan["axes"]:
        result = compare_axis(axis, rows, args.tolerance)
        print("%s: %s（阈值 %.0f%%）" % (result["axis"], result["verdict"],
                                         result["tolerance"] * 100))
        for change in result["changes"]:
            if change.get("status") == "missing":
                print("  %s → %s: 缺分析结果" % (change["from"], change["to"]))
            else:
                print("  %s → %s: 最大相对变化 %.2e%s"
                      % (change["from"], change["to"],
                         change["worst_relative_change"],
                         "" if change["converged"] else "  ← 未收敛"))
    if args.json:
        print(json.dumps({"table": rows, "axes": [
            compare_axis(axis, rows, args.tolerance) for axis in plan["axes"]]},
            ensure_ascii=False, indent=2))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="scan_convergence", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    subs = parser.add_subparsers(dest="command", required=True)

    plan = subs.add_parser("plan", help="按轴展开变体并写 scan.json（不提交任何作业）")
    plan.add_argument("config", help="基线 smj.yaml")
    plan.add_argument("--axis", action="append", required=True, metavar="KEY=V1,V2",
                      help="扫描轴，可重复，例如 --axis electrode.nlayers=3,4,5")
    plan.add_argument("--set", action="append", metavar="KEY=VALUE",
                      help="所有变体共用的常量覆盖，可重复")
    plan.add_argument("-o", "--out", default="scan", metavar="DIR", help="输出目录（默认 scan/）")
    plan.set_defaults(handler=cmd_plan)

    run = subs.add_parser("run", help="按顺序执行变体（需要 --yes）")
    run.add_argument("manifest", help="scan.json")
    run.add_argument("--from", dest="from_stage", default=None, help="smj run --from")
    run.add_argument("--to", dest="to_stage", default=None, help="smj run --to")
    run.add_argument("--yes", action="store_true", help="确认真的提交/执行")
    run.set_defaults(handler=cmd_run)

    compare = subs.add_parser("compare", help="汇总各变体的 S/G/ZT 并给逐轴收敛判定")
    compare.add_argument("manifest", help="scan.json")
    compare.add_argument("--tolerance", type=float, default=DEFAULT_TOL,
                         help="相邻两档的判收敛阈值（相对，默认 0.05）")
    compare.add_argument("--json", action="store_true", help="额外打印 JSON")
    compare.set_defaults(handler=cmd_compare)

    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except ScanError as exc:
        print("错误: %s" % exc, file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
