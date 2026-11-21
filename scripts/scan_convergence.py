#!/usr/bin/env python
"""收敛扫描：把「改一个参数重跑一遍」变成可复现的矩阵 + 对照表。

用法::

    # 1) 只出计划：按轴生成变体目录，每个目录一份独立 smj.yaml
    python scripts/scan_convergence.py plan smj.yaml \
        --axis electrode.nlayers=3,4,5 \
        --axis cluster.species_dir=/data/species/light,/data/species/tight \
        --out scan/ --pyramid electrode.au_per_electrode=20,30

    # 2) 真跑（必须有 --yes；AGENTS.md：只跑用户点名的计算）
    python scripts/scan_convergence.py run scan/scan.json --yes

    # 3) 汇总：读各变体的 junction_analyze/results.json（旧 analysis.json 回退），
    #    给逐轴、逐切片的收敛判定
    python scripts/scan_convergence.py compare scan/scan.json

设计取舍：

- **不猜物理**。脚本只认 ``section.key=值1,值2`` 这种点路径覆盖，具体拿哪些轴
  收敛（层数/团簇尺寸/基组档位/能量窗）由调用者按体系决定。
- **基组档位就是路径**。FHI-aims 的 tier 不是配置键，而是 ``cluster.species_dir``
  （或 ``paths.species_defaults``）指向不同物种目录，所以轴值写成目录路径即可；
  ``paths.species_dir`` 只是 ``paths.species_defaults`` 的兼容别名。
- **变体配置必须能读回来**。轴键先按 ``load_config`` 的字段表校验；路径值在
  计划里绝对化，写盘后再回读自检，避免"轴写了但没生效"的假收敛。
- **默认只写计划**。``plan`` 从不提交作业；``run`` 需要 ``--yes``，并且逐条打印
  实际执行的命令，方便先看一眼再决定。
- 变体的 ``smj.yaml`` 各自带 ``run.name`` 与 ``run.dir: run``，互不覆盖状态；
  ``scan.json`` 记录每个变体的真实 ``run_dir``（由 run/dir 规则算出）。
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from smj.config import ConfigError, known_keys, load_config  # noqa: E402

SCHEMA = "smj.convergence-scan.v1"
#: 阶段产物的正式位置：junction_analyze 阶段写 results.json。
STAGE_RESULTS = "junction_analyze/results.json"
#: 旧数据里可能只有 analysis.json（report 兼容读取），汇总时回退。
LEGACY_ANALYSIS = "junction_analyze/analysis.json"
#: 判「收敛」的默认阈值：相邻两档相对变化 < 5%。
DEFAULT_TOL = 0.05
#: 对照表里关心的量（与 analysis.json / demo-reference.json 同名）。
METRICS = ("S_uV_per_K", "log10_G_over_G0_like", "ZT_electronic", "k_e_pW_per_K")
#: 需要绝对化的路径键：变体目录与基线目录不同，相对路径会解析到错误位置。
PATH_KEYS = ("molecule.xyz", "cluster.species_dir", "paths.species_defaults",
             "paths.species_dir", "paths.templates")


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


def _absolute(value: object, base_dir: Path) -> object:
    """把路径字符串扩成绝对路径；其余值原样返回。"""
    if not isinstance(value, str) or not value:
        return value
    text = os.path.expanduser(os.path.expandvars(value))
    path = Path(text)
    return str(path if path.is_absolute() else (base_dir / path).resolve())


def _absolutize(config: dict, base_dir: Path) -> None:
    for key in PATH_KEYS:
        value = get_path(config, key)
        if value is None:
            continue
        set_path(config, key, _absolute(value, base_dir))


def _canonicalize_paths(config: dict) -> None:
    """``paths.species_dir`` 是 ``paths.species_defaults`` 的兼容别名。"""
    paths = config.get("paths")
    if not isinstance(paths, dict) or "species_dir" not in paths:
        return
    value = paths.pop("species_dir")
    if value:
        paths["species_defaults"] = value


def axis_problem(key: str) -> str | None:
    """轴键必须是 ``smj.config`` 认得的 ``section.key``；否则返回原因。"""
    section, _, leaf = key.partition(".")
    if not leaf or "." in leaf:
        return "轴名必须是 section.key 形式"
    allowed = known_keys().get(section)
    if allowed is None:
        return "未知配置段 %r" % section
    if leaf not in allowed:
        return "配置段 %s 没有字段 %r" % (section, leaf)
    return None


def _resolved_value(config: object, key: str) -> object:
    """按点路径取 ``load_config`` 解析后的值（``paths.species_dir`` 走别名）。"""
    if key == "paths.species_dir":
        return config.paths.species_defaults
    section, _, leaf = key.partition(".")
    node = getattr(config, section, None)
    return getattr(node, leaf, None) if node is not None else None


def build_plan(base_path: Path, axes: list[tuple[str, list[object]]],
               out_dir: Path, *, constants: dict[str, object] | None = None) -> dict:
    """把轴组合展开成变体清单（纯计算，不落盘）。"""
    if not axes:
        raise ScanError("至少要给一个 --axis")
    keys = [key for key, _ in axes]
    if len(set(keys)) != len(keys):
        raise ScanError("轴重复: %s" % keys)
    for key in keys:
        problem = axis_problem(key)
        if problem is not None:
            raise ScanError("轴 %s 不可用：%s；基组档位请用 cluster.species_dir 或 "
                            "paths.species_defaults" % (key, problem))
    base_path = Path(base_path)
    base = yaml.safe_load(base_path.read_text(encoding="utf-8")) or {}
    if not isinstance(base, dict):
        raise ScanError("%s 不是 mapping 配置" % base_path)
    base_dir = base_path.parent.resolve()
    variants: list[dict] = []
    for index, combo in enumerate(itertools.product(*[values for _, values in axes]), start=1):
        values = dict(zip(keys, combo))
        config = yaml.safe_load(yaml.safe_dump(base, allow_unicode=True)) or {}
        for key, value in (constants or {}).items():
            set_path(config, key, value)
        for key, value in values.items():
            set_path(config, key, value)
        _absolutize(config, base_dir)
        _canonicalize_paths(config)
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


def _same_path(left: object, right: object) -> bool:
    if left is None or right is None:
        return left is right
    try:
        return Path(str(left)).resolve() == Path(str(right)).resolve()
    except OSError:
        return False


def _verify_variant(plan: dict, variant: dict, config: object, config_path: Path) -> None:
    """回读自检：轴真的写进了配置、路径真的存在、没有未知字段。"""
    unknown = [warning for warning in config.warnings if "未知字段" in warning]
    if unknown:
        raise ScanError("变体 %s 的 smj.yaml 含未知字段（%s）；轴键必须是真实配置键"
                        % (variant["name"], "；".join(unknown)))
    base_dir = Path(plan["base"]).parent.resolve()
    raw = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    for axis in plan["axes"]:
        key = axis["key"]
        planned = variant["values"].get(key)
        expected = _absolute(planned, base_dir) if key in PATH_KEYS else planned
        stored = get_path(raw, key)
        if stored is None and key == "paths.species_dir":
            stored = get_path(raw, "paths.species_defaults")
        parsed = _resolved_value(config, key)
        if isinstance(expected, str):
            same = _same_path(stored, expected) or _same_path(parsed, expected)
        else:
            same = stored == expected or parsed == expected
        if not same:
            raise ScanError("变体 %s 的轴 %s 实际解析为 %r，与计划值 %r 不一致"
                            % (variant["name"], key, stored, planned))
    base = yaml.safe_load(Path(plan["base"]).read_text(encoding="utf-8")) or {}
    base_xyz = get_path(base, "molecule.xyz")
    xyz = getattr(config, "molecule_xyz", lambda: None)()
    if xyz is not None and not Path(xyz).is_file() and base_xyz:
        if Path(str(_absolute(base_xyz, base_dir))).is_file():
            raise ScanError("变体 %s 的 molecule.xyz 不存在：%s" % (variant["name"], xyz))
    species = _resolved_value(config, "cluster.species_dir") or _resolved_value(config, "paths.species_defaults")
    base_species = (get_path(base, "cluster.species_dir")
                    or get_path(base, "paths.species_defaults")
                    or get_path(base, "paths.species_dir"))
    if species and not Path(str(species)).is_dir() and base_species:
        if Path(str(_absolute(base_species, base_dir))).is_dir():
            raise ScanError("变体 %s 的 species_dir 不存在：%s" % (variant["name"], species))


def write_plan(plan: dict) -> Path:
    out_dir = Path(plan["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    for variant in plan["variants"]:
        target = out_dir / variant["name"]
        target.mkdir(parents=True, exist_ok=True)
        config_path = target / "smj.yaml"
        config_path.write_text(
            yaml.safe_dump(variant["config"], allow_unicode=True, sort_keys=False),
            encoding="utf-8", newline="\n")
        try:
            config = load_config(config_path)
        except ConfigError as exc:
            raise ScanError("变体 %s 的 smj.yaml 读不回来：%s" % (variant["name"], exc)) from exc
        _verify_variant(plan, variant, config, config_path)
        variant["run_dir"] = str(config.run_dir)
    manifest = out_dir / "scan.json"
    public = dict(plan)
    public["variants"] = [
        {"index": v["index"], "name": v["name"], "values": v["values"],
         "dir": str(out_dir / v["name"]), "config": str(out_dir / v["name"] / "smj.yaml"),
         "run_dir": v["run_dir"]}
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


def _run_dir_for(variant: dict) -> Path:
    """变体的真实 run 目录：清单里记的值优先，否则按配置规则重算。"""
    stored = variant.get("run_dir")
    if stored:
        return Path(stored)
    config_path = variant.get("config")
    if config_path and Path(str(config_path)).is_file():
        try:
            return load_config(Path(str(config_path))).run_dir
        except ConfigError:
            pass
    return Path(variant["dir"]) / "run" / str(variant.get("name") or "run")


def variant_metrics(variant: dict) -> dict:
    """读一个变体的分析结果；缺文件就标 missing，不伪造 0。"""
    run_dir = _run_dir_for(variant)
    for relative in (STAGE_RESULTS, LEGACY_ANALYSIS):
        path = run_dir / relative
        if not path.is_file():
            continue
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            return {"status": "unreadable", "run_dir": str(run_dir),
                    "source": str(path), "values": {}, "error": str(exc)}
        if not isinstance(payload, dict):
            return {"status": "unreadable", "run_dir": str(run_dir),
                    "source": str(path), "values": {}}
        values = {key: payload.get(key) for key in METRICS}
        status = "ok" if any(value is not None for value in values.values()) else "empty"
        return {"status": status, "run_dir": str(run_dir), "source": str(path),
                "values": values}
    return {"status": "missing", "run_dir": str(run_dir), "values": {}}


def _relative(new: float, base: float) -> float:
    if base == 0:
        return float("inf") if new else 0.0
    return abs(new - base) / abs(base)


def _slice_change(prev: dict | None, nxt: dict | None) -> dict:
    """比较「其余轴取值相同」的一对切片。"""
    if prev is None or nxt is None:
        return {"status": "unmatched"}
    before, after = prev["metrics"], nxt["metrics"]
    if before.get("status") == "missing" or after.get("status") == "missing":
        return {"status": "missing"}
    if before.get("status") != "ok" or after.get("status") != "ok":
        return {"status": "incomplete"}
    detail: dict[str, float] = {}
    worst = 0.0
    for metric in METRICS:
        old, new = before["values"].get(metric), after["values"].get(metric)
        if old is None or new is None:
            continue
        rel = _relative(float(new), float(old))
        detail[metric] = rel
        worst = max(worst, rel)
    if not detail:
        return {"status": "empty"}
    return {"status": "ok", "worst_relative_change": worst, "detail": detail}


def _merge_changes(prev_value: object, next_value: object, slices: list[dict],
                   tol: float) -> dict:
    bad = [item for item in slices if item["status"] != "ok"]
    if bad:
        statuses = sorted({item["status"] for item in bad})
        return {"from": prev_value, "to": next_value,
                "status": statuses[0] if len(statuses) == 1 else "incomplete",
                "slices": len(slices), "problems": statuses}
    worst = max(item["worst_relative_change"] for item in slices)
    detail: dict[str, float] = {}
    for item in slices:
        for metric, rel in item["detail"].items():
            detail[metric] = max(detail.get(metric, 0.0), rel)
    return {"from": prev_value, "to": next_value, "worst_relative_change": worst,
            "converged": worst < tol, "detail": detail, "slices": len(slices)}


def compare_axis(axis: dict, rows: list[dict], tol: float) -> dict:
    """沿一条轴逐档比较，其余轴取值相同的行才配对（切片比较）。"""
    key = axis["key"]
    values = list(axis["values"])
    groups: dict[tuple, dict[object, dict]] = {}
    for row in rows:
        vals = row.get("values") or {}
        if key not in vals:
            continue
        signature = tuple(sorted((name, value) for name, value in vals.items()
                                 if name != key))
        groups.setdefault(signature, {})[vals[key]] = row
    changes: list[dict] = []
    for prev_value, next_value in zip(values, values[1:]):
        slices = [_slice_change(by_value.get(prev_value), by_value.get(next_value))
                  for _, by_value in sorted(groups.items(), key=lambda item: str(item[0]))]
        if not slices:
            changes.append({"from": prev_value, "to": next_value, "status": "unmatched"})
        else:
            changes.append(_merge_changes(prev_value, next_value, slices, tol))
    if not changes:
        verdict = "insufficient"
    elif any(change.get("status") in ("missing", "unmatched", "incomplete", "empty")
             for change in changes):
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
            status = change.get("status")
            if status is not None:
                reason = {"missing": "缺分析结果", "unmatched": "缺少对应的切片",
                          "incomplete": "结果不可用", "empty": "指标全空"}.get(status, status)
                print("  %s → %s: %s" % (change["from"], change["to"], reason))
            else:
                print("  %s → %s: 最大相对变化 %.2e%s（%d 个切片）"
                      % (change["from"], change["to"],
                         change["worst_relative_change"],
                         "" if change["converged"] else "  ← 未收敛",
                         change.get("slices", 1)))
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
