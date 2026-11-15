"""``smj`` 命令行入口。

子命令一览::

    smj init      写一份 smj.yaml 模板
    smj plan      打印 6 个阶段、依赖与当前状态
    smj preflight 提交前检查（调用 smj.preflight）
    smj stage     单阶段操作：--prepare / --submit / --wait / --collect / --run
    smj run       按 --from/--to 链式推进（prepare → submit → wait → collect）
    smj status    查看 state.json 与阶段状态
    smj report    生成人读报告（调用 smj.report）
    smj analyze   透射谱 / Landauer 分析（调用 smj.analysis）
    smj audit     本地一致性检查（自包含，不依赖远端）
    smj review    结构复核（调用 smj.review）
    smj verify    安装自检（依赖 / 模板 / 阶段模块）

约定：

- 退出码：``0`` 成功；``1`` 失败 / 配置或输入错误；``2`` 未完成（等待超时，
  作业可能仍在运行；argparse 的用法错误也是 2）。
- ``--json`` 把结果按 JSON 打到 stdout（进度信息仍走 stderr / 被 ``-q`` 抑制）。
- ``smj.stages`` 一律延迟导入：阶段模块缺失时 ``plan`` / ``status`` 仍可用，
  ``stage`` / ``run`` 给出可读错误。
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib
import inspect
import json
import re
import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Callable

from . import __version__
from .config import DEFAULT_CONFIG_NAME, MODES, RunConfig, load_config, write_default_config
from .errors import SmjError, TransportError
from .remote import SMJ_STAGE_SCRIPTS
from .resources import hhmmss
from .runner import StageOutcome, StageRunner, lazy_stage_getter
from .state import (
    REMOTE_STAGES,
    STAGE_ORDER,
    RunLayout,
    RunState,
    normalize_stage,
    stage_range,
)
from .transport import make_transport

# 退出码（README / 文档引用这里的常量语义）
EXIT_OK = 0
EXIT_FAILED = 1
EXIT_UNFINISHED = 2

STAGE_ALIASES = STAGE_ORDER  # 仅用于 help 文案

_AUX_MODULES = {
    "preflight": ("smj.preflight", ("run_preflight", "preflight", "check", "run")),
    "report": ("smj.report", ("generate_report", "build_report", "write_report", "run_report", "report")),
    "analyze": ("smj.analysis", ("analyze_run", "analyze", "run_analysis")),
    "review": ("smj.review", ("review_run", "review", "run_review")),
}


# --------------------------------------------------------------------------- #
# 小工具
# --------------------------------------------------------------------------- #
def _slug(text: str) -> str:
    """把目录名变成合法的非 VASP run 名（小写，允许 ``[a-z0-9._-]``）。"""
    slug = re.sub(r"[^A-Za-z0-9._-]+", "-", str(text or "")).strip("-._").lower()
    return slug or "my-run"


def _error(message: str, *, debug: bool = False) -> None:
    print("错误: %s" % message, file=sys.stderr)
    if debug:
        import traceback

        traceback.print_exc()


def _jsonable(value: Any) -> Any:
    """把 Path / dataclass / 带 to_dict() 的对象转成可 json.dumps 的结构。"""
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_jsonable(v) for v in value]
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return _jsonable(dataclasses.asdict(value))
    to_dict = getattr(value, "to_dict", None)
    if callable(to_dict):
        try:
            return _jsonable(to_dict())
        except Exception:  # noqa: BLE001 - 尽力而为
            pass
    return str(value)


def _as_mapping(result: Any) -> dict[str, Any]:
    """把 aux 模块的返回值（dict / dataclass / to_dict() / bool）统一成 dict。"""
    if result is None:
        return {}
    if isinstance(result, bool):
        return {"ok": result}
    if isinstance(result, Mapping):
        return dict(result)
    if dataclasses.is_dataclass(result) and not isinstance(result, type):
        return dataclasses.asdict(result)
    to_dict = getattr(result, "to_dict", None)
    if callable(to_dict):
        try:
            return dict(to_dict())
        except Exception:  # noqa: BLE001
            pass
    if isinstance(result, (list, tuple)):
        return {"issues": list(result)}
    return {"result": result}


def _truthy_ok(data: Mapping[str, Any]) -> bool:
    """从 preflight / review 的返回里推断整体是否通过。"""
    for key in ("ok", "passed", "success"):
        if key in data:
            return bool(data[key])
    if "blocked" in data:
        return not bool(data["blocked"])
    issues = data.get("issues") or data.get("problems") or []
    if isinstance(issues, (list, tuple)):
        for issue in issues:
            severity = ""
            if isinstance(issue, Mapping):
                severity = str(issue.get("severity") or issue.get("level") or "").lower()
            if severity in ("error", "fatal", "blocked"):
                return False
    return True


def _emit_json(payload: Any) -> None:
    print(json.dumps(_jsonable(payload), ensure_ascii=False, indent=2))


def _print_issues(issues: Any) -> None:
    if not issues:
        return
    for issue in issues:
        if isinstance(issue, Mapping):
            severity = issue.get("severity") or issue.get("level") or "info"
            stage = issue.get("stage")
            message = (issue.get("message") or issue.get("detail")
                       or issue.get("issue") or json.dumps(_jsonable(issue), ensure_ascii=False))
            prefix = "[%s]" % severity
            if stage:
                prefix += " %s:" % stage
            print("  %s %s" % (prefix, message))
        else:
            print("  - %s" % issue)


def _print_mapping(data: Mapping[str, Any], *, title: str | None = None) -> None:
    if title:
        print("%s:" % title)
    if not data:
        print("  （空）")
        return
    for key in sorted(data):
        value = data[key]
        if key in ("issues", "problems") and isinstance(value, (list, tuple)):
            _print_issues(value)
            continue
        if isinstance(value, (Mapping, list, tuple)):
            text = json.dumps(_jsonable(value), ensure_ascii=False)
            if len(text) > 400:
                text = text[:397] + "..."
            print("  %s: %s" % (key, text))
        else:
            print("  %s: %s" % (key, value))


def _print_table(headers: tuple[str, ...], rows: list[tuple[Any, ...]]) -> None:
    widths = [len(h) for h in headers]
    for row in rows:
        for idx, cell in enumerate(row):
            widths[idx] = max(widths[idx], len(str(cell)))
    fmt = "  ".join("%%-%ds" % width for width in widths)
    print(fmt % headers)
    print("  ".join("-" * width for width in widths))
    for row in rows:
        print(fmt % tuple(str(cell) for cell in row))


def _load_config(args: argparse.Namespace) -> RunConfig:
    config = load_config(Path(args.config))
    if not getattr(args, "quiet", False):
        for warning in config.warnings:
            print("警告: %s" % warning, file=sys.stderr)
    return config


def _overrides_from_args(args: argparse.Namespace) -> dict[str, Any]:
    overrides: dict[str, Any] = {}
    for key in ("xyz", "geometry", "scf_dir", "te", "metadata", "out"):
        value = getattr(args, key, None)
        if value:
            overrides[key] = value
    return overrides


def _make_runner(args: argparse.Namespace, config: RunConfig) -> StageRunner:
    logger = None if getattr(args, "quiet", False) else (lambda message: print(message))
    return StageRunner(
        config,
        logger=logger,
        overrides=_overrides_from_args(args),
        dry_run=True if getattr(args, "dry_run", False) else None,
        stage_getter=getattr(args, "stage_getter", None),
        force=bool(getattr(args, "force", False)),
    )


def _state_for(config: RunConfig) -> tuple[RunLayout, RunState]:
    layout = RunLayout(config.run_dir)
    state = RunState.load(layout.state_path, run=config.name)
    return layout, state


# --------------------------------------------------------------------------- #
# 阶段元数据（plan / audit 用；stages 缺失时退回静态表）
# --------------------------------------------------------------------------- #
def _fallback_meta(name: str) -> dict[str, Any]:
    index = STAGE_ORDER.index(name)
    return {
        "stage": name,
        "kind": "slurm" if name in REMOTE_STAGES else "local",
        "script": SMJ_STAGE_SCRIPTS.get(name),
        "depends": (STAGE_ORDER[index - 1],) if index > 0 else (),
        "input": "",
        "outputs": (),
        "needs_submit": name in REMOTE_STAGES,
        "available": False,
    }


def _stage_meta(name: str, stage_getter: Callable[[str], Any] | None = None) -> dict[str, Any]:
    getter = stage_getter or lazy_stage_getter
    try:
        stage = getter(name)
    except (SmjError, KeyError, ImportError, AttributeError):
        return _fallback_meta(name)
    needs_submit = bool(getattr(stage, "needs_submit", name in REMOTE_STAGES))
    return {
        "stage": name,
        "kind": "slurm" if needs_submit else "local",
        "script": getattr(stage, "script", None) or SMJ_STAGE_SCRIPTS.get(name),
        "depends": tuple(getattr(stage, "depends", ()) or ()),
        "input": str(getattr(stage, "input_hint", "") or ""),
        "outputs": tuple(getattr(stage, "outputs", ()) or ()),
        "needs_submit": needs_submit,
        "available": True,
    }


# --------------------------------------------------------------------------- #
# aux 模块延迟导入
# --------------------------------------------------------------------------- #
def _optional_module(module_name: str, label: str) -> Any:
    try:
        return importlib.import_module(module_name)
    except ImportError as exc:  # pragma: no cover - 正常安装时不会走到
        raise SmjError(
            "%s 功能暂不可用：导入 %s 失败（%s）。该功能是 smj 的可选组件，"
            "请确认安装完整或改用 `smj stage` / `smj status`。" % (label, module_name, exc)
        ) from exc


def _find_callable(module: Any, names: tuple[str, ...]) -> Callable[..., Any]:
    for name in names:
        candidate = getattr(module, name, None)
        if callable(candidate):
            return candidate
    raise SmjError(
        "模块 %s 没有预期的入口函数（找过: %s）" % (module.__name__, ", ".join(names))
    )


def _invoke(func: Callable[..., Any], config: RunConfig, **kwargs: Any) -> Any:
    """按签名过滤 kwargs 后调用；首个位置参数传 config（aux 模块的约定）。"""
    try:
        signature = inspect.signature(func)
    except (TypeError, ValueError):  # pragma: no cover - 内置函数
        return func(config, **kwargs)
    parameters = signature.parameters
    accepts_any = any(p.kind is inspect.Parameter.VAR_KEYWORD for p in parameters.values())
    if accepts_any:
        filtered = dict(kwargs)
    else:
        filtered = {key: value for key, value in kwargs.items() if key in parameters}
    positional = [
        p for p in parameters.values()
        if p.kind in (inspect.Parameter.POSITIONAL_ONLY,
                      inspect.Parameter.POSITIONAL_OR_KEYWORD)
    ]
    if positional:
        return func(config, **filtered)
    return func(**filtered)


# --------------------------------------------------------------------------- #
# 输出辅助
# --------------------------------------------------------------------------- #
def _print_prepared(payload: Mapping[str, Any], *, stage: str) -> None:
    print("[%s] 已生成输入（prepared）" % stage)
    for item in payload.get("inputs") or ():
        print("  输入: %s" % item)
    script = payload.get("script")
    if script:
        print("  作业脚本: %s" % script)
    for command in payload.get("commands") or ():
        print("  命令: %s" % command)
    for note in payload.get("notes") or ():
        print("  提示: %s" % note)


def _print_submission(payload: Mapping[str, Any], *, stage: str) -> None:
    job_id = payload.get("job_id")
    print("[%s] 已提交 job %s -> %s" % (stage, job_id, payload.get("remote_dir")))
    if payload.get("cleaned"):
        print("  清理旧产物: %s" % ", ".join(payload["cleaned"]))
    if payload.get("uploaded"):
        print("  上传: %s" % ", ".join(payload["uploaded"]))
    if payload.get("skipped"):
        print("  跳过: %s" % ", ".join(payload["skipped"]))


def _print_wait(payload: Mapping[str, Any], *, stage: str) -> None:
    if payload.get("failed"):
        print("[%s] 作业失败（exit_code=%s）：%s"
              % (stage, payload.get("exit_code"), payload.get("note") or "未知原因"))
    elif payload.get("finished"):
        print("[%s] 作业结束（exit_code=%s，耗时 %s）"
              % (stage, payload.get("exit_code"), hhmmss(payload.get("elapsed_s") or 0)))
    else:
        print("[%s] 等待超时：%s" % (stage, payload.get("note") or "作业可能仍在运行"))


def _print_collected(payload: Mapping[str, Any], *, stage: str) -> None:
    print("[%s] 已回取并解析（collected）" % stage)
    artifacts = payload.get("artifacts") or ()
    if artifacts:
        print("  产物: %s" % ", ".join(str(a) for a in artifacts))
    results = payload.get("results") or {}
    if results:
        _print_mapping(results, title="  结果")
    for note in payload.get("notes") or ():
        print("  提示: %s" % note)


def _print_outcome(outcome: StageOutcome) -> None:
    if outcome.ok:
        print("[%s] %s" % (outcome.stage, outcome.message or ("%s 完成" % outcome.status)))
        return
    print("[%s] %s（status=%s）" % (outcome.stage, outcome.message or "未完成", outcome.status))


# --------------------------------------------------------------------------- #
# 子命令实现
# --------------------------------------------------------------------------- #
def cmd_init(args: argparse.Namespace) -> int:
    target = Path(args.path)
    if target.suffix.lower() in (".yaml", ".yml"):
        config_path = target
        base = target.parent
    else:
        config_path = target / DEFAULT_CONFIG_NAME
        base = target
    name = args.name or _slug(base.resolve().name or "my-run")
    path = write_default_config(
        config_path, name=name, xyz=args.xyz, mode=args.mode, overwrite=args.force)
    print("已写入配置: %s" % path)
    print("下一步:")
    print("  1) 准备分子结构 %s（或改 molecule.xyz 指向真实文件）" % args.xyz)
    print("  2) smj plan -c %s" % path)
    print("  3) smj run -c %s --dry-run   # 只生成输入，不提交" % path)
    return EXIT_OK


def cmd_plan(args: argparse.Namespace) -> int:
    config = _load_config(args)
    layout, state = _state_for(config)
    rows = []
    for name in STAGE_ORDER:
        meta = _stage_meta(name, getattr(args, "stage_getter", None))
        script = meta["script"] or "-"
        if name == "molecule_opt":
            engine = str(getattr(config.molecule, "opt_engine", "aims") or "aims").lower()
            script = "sub_aims.sh" if engine == "aims" else "sub_g16.sh"
        rows.append({
            "stage": name,
            "kind": meta["kind"],
            "script": script,
            "depends": ", ".join(meta["depends"]) or "-",
            "input": meta["input"],
            "status": state.status(name),
        })
    next_stage = state.first_incomplete()
    payload = {
        "run": config.name,
        "mode": config.mode,
        "run_dir": str(config.run_dir),
        "state_path": str(layout.state_path),
        "stages": rows,
        "next": next_stage,
    }
    if args.json:
        _emit_json(payload)
        return EXIT_OK
    print("run: %s   mode: %s" % (config.name, config.mode))
    print("dir: %s" % config.run_dir)
    print()
    table = [
        (index, row["stage"], row["kind"], row["script"], row["depends"],
         row["status"], row["input"] or "-")
        for index, row in enumerate(rows, 1)
    ]
    _print_table(("#", "stage", "kind", "script", "depends", "status", "input"), table)
    if next_stage:
        print("\n下一步: smj run --from %s -c %s" % (next_stage, args.config))
    else:
        print("\n全部阶段已完成。")
    return EXIT_OK


def cmd_preflight(args: argparse.Namespace) -> int:
    config = _load_config(args)
    transport = None
    if config.mode == "slurm-ssh" and not args.no_remote:
        try:
            transport = make_transport(config)
        except TransportError as exc:
            if not args.quiet:
                print("警告: 无法建立远端传输（跳过远端检查）: %s" % exc, file=sys.stderr)
    stage = normalize_stage(args.stage) if args.stage else None
    module = _optional_module(_AUX_MODULES["preflight"][0], "preflight（提交前检查）")
    func = _find_callable(module, _AUX_MODULES["preflight"][1])
    result = _invoke(
        func, config, transport=transport, stage=stage,
        dry_run=config.mode == "dry-run")
    data = _as_mapping(result)
    ok = _truthy_ok(data)
    if args.json:
        _emit_json({"run": config.name, "ok": ok, **data})
        return EXIT_OK if ok else EXIT_FAILED
    print("preflight: %s（mode=%s）" % (config.name, config.mode))
    if "issues" in data or "problems" in data:
        _print_issues(data.get("issues") or data.get("problems"))
    _print_mapping(
        {k: v for k, v in data.items() if k not in ("issues", "problems", "ok")})
    print("结论: %s" % ("通过" if ok else "未通过"))
    return EXIT_OK if ok else EXIT_FAILED


def cmd_stage(args: argparse.Namespace) -> int:
    config = _load_config(args)
    stage = normalize_stage(args.stage)
    action = "run"
    for name in ("prepare", "submit", "wait", "collect", "run"):
        if getattr(args, name, False):
            action = name
            break

    dry_run = bool(args.dry_run) or config.mode == "dry-run"
    if action in ("submit", "wait") and dry_run:
        raise SmjError(
            "dry-run 模式不能 %s；把 run.mode 改成 local / slurm-ssh，"
            "或去掉 --dry-run" % action)

    runner = _make_runner(args, config)

    if action == "prepare":
        payload = runner.prepare(stage, force=args.force)
        if args.json:
            _emit_json(payload)
        else:
            _print_prepared(payload, stage=stage)
        return EXIT_OK

    if action == "submit":
        payload = runner.submit(stage, force=args.force)
        if args.json:
            _emit_json(payload)
        else:
            _print_submission(payload, stage=stage)
        return EXIT_OK

    if action == "wait":
        payload = runner.wait(stage)
        if args.json:
            _emit_json(payload)
        else:
            _print_wait(payload, stage=stage)
        if payload.get("failed"):
            return EXIT_FAILED
        if not payload.get("finished"):
            return EXIT_UNFINISHED
        return EXIT_OK

    if action == "collect":
        payload = runner.collect(stage, force=args.force)
        if args.json:
            _emit_json(payload)
        else:
            _print_collected(payload, stage=stage)
        return EXIT_OK

    outcome = runner.run_stage(stage, force=args.force)
    if args.json:
        _emit_json(outcome.to_dict())
    else:
        _print_outcome(outcome)
    if outcome.ok:
        return EXIT_OK
    return EXIT_UNFINISHED if outcome.status == "waiting" else EXIT_FAILED


def cmd_run(args: argparse.Namespace) -> int:
    config = _load_config(args)
    start = args.from_stage or args.stage
    end = args.to_stage or args.stage
    names = stage_range(start, end)  # 提前校验，避免跑一半才发现区间非法

    keep_going: bool | None = None
    if args.keep_going:
        keep_going = True
    elif args.fail_fast:
        keep_going = False

    runner = _make_runner(args, config)
    outcomes = runner.run(start, end, keep_going=keep_going, force=args.force)

    if args.json:
        _emit_json({
            "run": config.name,
            "mode": config.mode,
            "stages": names,
            "outcomes": [outcome.to_dict() for outcome in outcomes],
            "ok": all(outcome.ok for outcome in outcomes),
        })
    else:
        for outcome in outcomes:
            _print_outcome(outcome)

    if all(outcome.ok for outcome in outcomes):
        if not args.json:
            print("\n完成：%s" % ", ".join(names))
        return EXIT_OK
    unfinished = any(
        (not outcome.ok) and outcome.status == "waiting" for outcome in outcomes)
    return EXIT_UNFINISHED if unfinished else EXIT_FAILED


def cmd_status(args: argparse.Namespace) -> int:
    config = _load_config(args)
    layout, state = _state_for(config)
    if args.stage:
        only = normalize_stage(args.stage)
        rows = [row for row in _status_rows(state) if row["stage"] == only]
    else:
        rows = _status_rows(state)
    payload = {
        "run": config.name,
        "mode": config.mode,
        "run_dir": str(layout.run_dir),
        "state_path": str(layout.state_path),
        "state_exists": layout.state_path.is_file(),
        "stages": rows,
        "next": state.first_incomplete(),
    }
    if args.json:
        _emit_json(payload)
        return EXIT_OK
    print("run: %s   mode: %s" % (config.name, config.mode))
    print("dir: %s" % layout.run_dir)
    print("state: %s%s" % (
        layout.state_path, "" if layout.state_path.is_file() else "（尚未生成）"))
    print()
    table = [
        (row["stage"], row["status"], row["job_id"] or "-",
         hhmmss(row["elapsed_s"] or 0) if row["elapsed_s"] else "-",
         row["note"] or "-")
        for row in rows
    ]
    _print_table(("stage", "status", "job_id", "elapsed", "note"), table)
    if state.first_incomplete():
        print("\n下一步: smj run --from %s -c %s"
              % (state.first_incomplete(), args.config))
    return EXIT_OK


def _status_rows(state: RunState) -> list[dict[str, Any]]:
    rows = []
    for name in STAGE_ORDER:
        stage = state.stage(name)
        rows.append({
            "stage": name,
            "status": stage.status,
            "job_id": stage.job_id,
            "remote_dir": stage.remote_dir,
            "elapsed_s": stage.elapsed_s,
            "note": stage.note,
            "error": stage.error,
            "artifacts": list(stage.artifacts),
        })
    return rows


def cmd_report(args: argparse.Namespace) -> int:
    config = _load_config(args)
    module = _optional_module(_AUX_MODULES["report"][0], "report（报告生成）")
    func = _find_callable(module, _AUX_MODULES["report"][1])
    result = _invoke(func, config, force=args.force,
                     include_checks=not args.no_checks)
    data = _as_mapping(result)
    ok = _truthy_ok(data)
    if args.json:
        _emit_json(data)
    else:
        _print_mapping(data, title="report: %s" % config.name)
    return EXIT_OK if ok else EXIT_FAILED


def cmd_analyze(args: argparse.Namespace) -> int:
    config = _load_config(args)
    module = _optional_module(_AUX_MODULES["analyze"][0], "analyze（透射谱分析）")
    func = _find_callable(module, _AUX_MODULES["analyze"][1])
    result = _invoke(func, config, force=args.force)
    data = _as_mapping(result)
    ok = _truthy_ok(data)
    if args.json:
        _emit_json(data)
    else:
        _print_mapping(data, title="analyze: %s" % config.name)
    return EXIT_OK if ok else EXIT_FAILED


def cmd_review(args: argparse.Namespace) -> int:
    config = _load_config(args)

    # 2026-11-15：A2 人工确认更名为「连接核对」（connection check），默认自动登记。
    # --approve 是手工补签/重签路径（fail-closed），签完在同一条命令里看门过没过。
    if getattr(args, "approve", False):
        if not (args.who or "").strip():
            print("错误: --approve 需要 --who <签署人>（连接核对必须记录是谁看的）")
            return 2
        from .a2 import approve_connection_check  # 本地模块，非可选依赖

        approval = approve_connection_check(
            config, who=args.who, stage=args.stage, no_figure=args.no_figure,
            png=args.png, tcontrol=args.tcontrol, force=args.force)
        if args.json:
            _emit_json(approval.to_dict())
        else:
            print("连接核对: %s" % approval.status)
            print("  %s" % approval.reason)
            for key in ("build_dir", "who", "when", "fingerprint",
                        "audit_path", "png_path"):
                value = approval.get(key)
                if value:
                    print("  %s: %s" % (key, value))
            if approval.checks_failed:
                print("  未通过检查: %s" % ", ".join(approval.checks_failed))
        if not approval.ok:
            return approval.exit_code

    module = _optional_module(_AUX_MODULES["review"][0], "review（结构复核）")
    func = _find_callable(module, _AUX_MODULES["review"][1])
    kwargs: dict[str, Any] = {"force": args.force}
    if args.accept_deviation is not None:
        kwargs["accept_deviation"] = args.accept_deviation
    connection_check = getattr(args, "connection_check", None)
    if connection_check is None:
        # 旧 --surface-check/--no-surface-check、--a2/--no-a2 别名（隐藏，仍可用）
        connection_check = getattr(args, "legacy_surface_check", None)
    if connection_check is None:
        connection_check = getattr(args, "legacy_a2", None)
    if connection_check is not None:
        # 命令行开关只作用于本次 review 调用（不写回 smj.yaml、不进 run_spec）：
        # --no-connection-check 等价 connection_check=false，缺记录/指纹失效只降为提示。
        review_cfg = getattr(config, "review", None)
        if review_cfg is not None:
            kwargs["review_config"] = dataclasses.replace(
                review_cfg, connection_check=bool(connection_check))
    result = _invoke(func, config, **kwargs)
    data = _as_mapping(result)
    ok = _truthy_ok(data)
    if args.json:
        _emit_json({"run": config.name, "ok": ok, **data})
        return EXIT_OK if ok else EXIT_FAILED
    print("review: %s" % config.name)
    if "issues" in data or "problems" in data:
        _print_issues(data.get("issues") or data.get("problems"))
    _print_mapping(
        {k: v for k, v in data.items() if k not in ("issues", "problems", "ok")})
    print("结论: %s" % ("通过" if ok else "未通过"))
    return EXIT_OK if ok else EXIT_FAILED


def cmd_audit(args: argparse.Namespace) -> int:
    """本地一致性检查：不连远端，只看 state.json 与 run 目录。"""
    config = _load_config(args)
    layout, state = _state_for(config)
    issues: list[dict[str, Any]] = []

    def add(severity: str, code: str, message: str, stage: str | None = None) -> None:
        issues.append({"severity": severity, "code": code, "stage": stage,
                       "message": message})

    if not layout.state_path.is_file():
        # run 尚未开始：目录和 state 都不存在是正常状态，不算错误。
        add("warning", "state-missing",
            "没有 state.json，run 尚未开始: %s" % layout.state_path)
        if not layout.run_dir.is_dir():
            add("info", "run-dir-missing",
                "run 目录尚未创建: %s" % layout.run_dir)
    elif not layout.run_dir.is_dir():
        add("error", "run-dir-missing",
            "state.json 存在但 run 目录不存在: %s" % layout.run_dir)
    if state.run and state.run != config.name:
        add("warning", "run-name-mismatch",
            "state.json 记录的 run=%r 与配置 %r 不一致" % (state.run, config.name))

    # 2026-11 P0-1：运行参数快照必须能对账；漂移是 error，缺失在已开始的 run 里是 warning。
    if layout.state_path.is_file():
        if not layout.run_spec_path.is_file():
            add("warning", "run-spec-missing",
                "缺 run_spec.json：本次 run 的参数/软件版本没有随 run 落盘")
        else:
            from .provenance import load_run_spec, protocol_sha256
            spec = load_run_spec(layout.run_dir)
            if spec is None:
                add("error", "run-spec-unreadable",
                    "run_spec.json 存在但不可解析: %s" % layout.run_spec_path)
            else:
                if spec.get("reconstructed"):
                    add("warning", "run-spec-reconstructed",
                        "run_spec.json 是事后重建（reconstructed=true），不是原始参数快照")
                try:
                    current = protocol_sha256(config)
                except Exception as exc:  # noqa: BLE001 - 审计不应因配置异常崩掉
                    add("warning", "run-spec-check-failed",
                        "无法计算当前配置的协议指纹: %s" % exc)
                else:
                    stored = spec.get("protocol_sha256")
                    if stored and stored != current:
                        add("error", "run-spec-drift",
                            "run_spec.protocol_sha256=%s 与当前配置=%s 不一致"
                            % (str(stored)[:12], current[:12]))

    for name in STAGE_ORDER:
        stage = state.stage(name)
        meta = _stage_meta(name, getattr(args, "stage_getter", None))
        if (config.mode == "slurm-ssh" and stage.status in ("submitted", "waiting")
                and not stage.job_id):
            add("error", "missing-job-id",
                "status=%s 但没有 job_id" % stage.status, name)
        if stage.status in ("finished", "collected"):
            stage_dir = layout.stage_dir(name)
            if not stage_dir.is_dir():
                add("error", "stage-dir-missing",
                    "status=%s 但阶段目录不存在: %s" % (stage.status, stage_dir), name)
            else:
                # 优先用本阶段实际记录的 artifacts（引擎分派后 gaussian/aims
                # 的产物名不同）；只有旧 run 没记录时才回退模块声明的 outputs。
                expected = list(stage.artifacts) or list(meta.get("outputs") or ())
                for artifact in expected:
                    if not (stage_dir / artifact).exists():
                        add("warning", "artifact-missing",
                            "缺少产物 %s（可在远端日志确认）" % artifact, name)
        if stage.status == "collected" and not stage.artifacts:
            add("warning", "no-artifacts", "collected 但 artifacts 为空", name)
        if stage.status == "failed" and not stage.error:
            add("warning", "failed-without-error", "failed 但没有 error 信息", name)
        if stage.status == "pending" and layout.stage_dir(name).is_dir():
            add("info", "leftover-dir",
                "状态仍为 pending，但阶段目录已存在（可能是上次残留）", name)

    ok = not any(issue["severity"] == "error" for issue in issues)
    payload = {
        "run": config.name,
        "mode": config.mode,
        "run_dir": str(layout.run_dir),
        "ok": ok,
        "issues": issues,
        "stages": _status_rows(state),
    }
    if args.json:
        _emit_json(payload)
        return EXIT_OK if ok else EXIT_FAILED
    print("audit: %s（mode=%s）" % (config.name, config.mode))
    print("dir: %s" % layout.run_dir)
    print()
    _print_table(
        ("stage", "status", "job_id", "elapsed"),
        [(row["stage"], row["status"], row["job_id"] or "-",
          hhmmss(row["elapsed_s"] or 0) if row["elapsed_s"] else "-")
         for row in payload["stages"]])
    print()
    if issues:
        print("问题:")
        _print_issues(issues)
    else:
        print("问题: 无")
    errors = sum(1 for issue in issues if issue["severity"] == "error")
    warnings = sum(1 for issue in issues if issue["severity"] == "warning")
    print("结论: %s（error=%d, warning=%d）"
          % ("通过" if ok else "未通过", errors, warnings))
    return EXIT_OK if ok else EXIT_FAILED


def cmd_verify(args: argparse.Namespace) -> int:
    checks: list[dict[str, Any]] = []

    def check(name: str, ok: bool, detail: str = "", severity: str = "error") -> None:
        checks.append({"check": name, "ok": bool(ok), "detail": detail,
                       "severity": severity})

    version = "%d.%d.%d" % sys.version_info[:3]
    check("python", sys.version_info >= (3, 10), version)

    try:
        import yaml  # noqa: F401

        check("pyyaml", True, getattr(yaml, "__version__", "已安装"))
    except Exception as exc:  # noqa: BLE001
        check("pyyaml", False, "缺少 PyYAML（pip install PyYAML）: %s" % exc)

    try:
        import numpy  # noqa: F401

        check("numpy", True, getattr(numpy, "__version__", "已安装"))
    except Exception as exc:  # noqa: BLE001
        check("numpy", False, "缺少 numpy（pip install numpy）: %s" % exc)

    try:
        from ._vendor.junction import artifacts, jobwait  # noqa: F401

        check("vendored-junction", True, "jobwait / artifacts 可导入")
    except Exception as exc:  # noqa: BLE001
        check("vendored-junction", False, "vendored junction 不可导入: %s" % exc)

    templates = Path(__file__).resolve().parent / "_vendor" / "junction" / "templates"
    missing = [name for name in ("sub_g16.sh", "sub_aims.sh", "sub_te.sh")
               if not (templates / name).is_file()]
    check("templates", not missing,
          "缺少: %s" % ", ".join(missing) if missing else str(templates))

    try:
        import importlib as _importlib

        _stages = _importlib.import_module("smj.stages")
        from smj.stages.base import all_stages as _all_stages

        _registered = [s.name for s in _all_stages()]
        _missing = [name for name in STAGE_ORDER if name not in _registered]
        check("stages", not _missing,
              "%d 个阶段模块可导入%s" % (len(_registered),
                                       "" if not _missing else "；缺 %s" % ", ".join(_missing)))
    except Exception as exc:  # noqa: BLE001
        check("stages", False, "smj.stages 不可导入: %s" % exc)

    ok = all(item["ok"] for item in checks if item["severity"] == "error")
    if args.json:
        _emit_json({"ok": ok, "checks": checks})
        return EXIT_OK if ok else EXIT_FAILED
    print("smj verify（%s）" % ("通过" if ok else "未通过"))
    _print_table(
        ("check", "result", "detail"),
        [(item["check"], "ok" if item["ok"] else "FAIL", item["detail"])
         for item in checks])
    return EXIT_OK if ok else EXIT_FAILED


# --------------------------------------------------------------------------- #
# 参数解析
# --------------------------------------------------------------------------- #
def _add_common(parser: argparse.ArgumentParser, *, with_config: bool = True) -> None:
    if with_config:
        parser.add_argument(
            "-c", "--config", default=DEFAULT_CONFIG_NAME, metavar="PATH",
            help="配置文件路径（默认 %s）" % DEFAULT_CONFIG_NAME)
    parser.add_argument("-q", "--quiet", action="store_true", help="不打印进度信息")
    parser.add_argument("--json", action="store_true", help="以 JSON 输出结果")
    parser.add_argument("--debug", action="store_true", help="出错时打印完整 traceback")


def _add_overrides(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("输入覆盖")
    group.add_argument("--xyz", metavar="PATH", help="分子结构 xyz（覆盖 molecule.xyz）")
    group.add_argument("--geometry", metavar="PATH", help="电极结几何文件（覆盖 electrode_build 输入）")
    group.add_argument("--scf-dir", dest="scf_dir", metavar="DIR", help="aims_scf 输出目录覆盖")
    group.add_argument("--te", metavar="PATH", help="TE.dat / tcontrol 输入目录覆盖")
    group.add_argument("--metadata", metavar="PATH", help="电极元数据 JSON 覆盖")
    group.add_argument("--out", metavar="DIR", help="结果输出目录覆盖")


def _add_stage_group(parser: argparse.ArgumentParser) -> None:
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--prepare", action="store_true", help="只生成阶段输入")
    group.add_argument("--submit", action="store_true", help="提交已 prepare 的阶段")
    group.add_argument("--wait", action="store_true", help="等待已提交作业结束")
    group.add_argument("--collect", action="store_true", help="回取产物并解析")
    group.add_argument("--run", action="store_true", help="prepare→submit→wait→collect（默认）")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="smj",
        description="SMJ 单分子结计算工作流（Gaussian → FHI-aims → AITRANSS → Landauer）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "常用:\n"
            "  smj init .\n"
            "  smj plan\n"
            "  smj run --dry-run\n"
            "  smj run --from aims_scf --to aitranss_te\n"
            "  smj status\n"
        ),
    )
    parser.add_argument("--version", action="version", version="smj %s" % __version__)
    subparsers = parser.add_subparsers(dest="command", metavar="<command>")

    # init ------------------------------------------------------------------ #
    init = subparsers.add_parser(
        "init", help="写一份 smj.yaml 模板",
        description="在目录下写 smj.yaml（位置参数也可直接给 .yaml 路径）。")
    init.add_argument("path", nargs="?", default=".", metavar="PATH",
                      help="目标目录或 yaml 路径（默认当前目录）")
    init.add_argument("--name", help="run 名称（默认取目录名的小写 slug）")
    init.add_argument("--xyz", default="molecule.xyz", metavar="PATH",
                      help="molecule.xyz 路径（默认 molecule.xyz）")
    init.add_argument("--mode", choices=MODES, default="dry-run",
                      help="run.mode（默认 dry-run）")
    init.add_argument("--force", action="store_true", help="覆盖已存在的配置")
    init.add_argument("--debug", action="store_true", help="出错时打印完整 traceback")
    init.set_defaults(handler=cmd_init)

    # plan ------------------------------------------------------------------ #
    plan = subparsers.add_parser(
        "plan", help="打印 6 个阶段、依赖与当前状态",
        description="只读：不生成文件、不连远端。")
    _add_common(plan)
    plan.set_defaults(handler=cmd_plan)

    # preflight ------------------------------------------------------------- #
    preflight = subparsers.add_parser(
        "preflight", help="提交前检查（结构/输入/环境）",
        description="调用 smj.preflight；slurm-ssh 模式下会附带远端可达性检查。")
    _add_common(preflight)
    preflight.add_argument("--stage", metavar="STAGE", help="只检查某个阶段")
    preflight.add_argument("--no-remote", action="store_true",
                           help="跳过远端检查（离线时使用）")
    preflight.set_defaults(handler=cmd_preflight)

    # stage ----------------------------------------------------------------- #
    stage = subparsers.add_parser(
        "stage", help="单阶段操作：prepare / submit / wait / collect / run",
        description="不带动作 flag 时等价于 --run（全链路）。")
    stage.add_argument("stage", metavar="STAGE",
                       help="阶段名或别名: %s" % ", ".join(STAGE_ALIASES))
    _add_stage_group(stage)
    _add_overrides(stage)
    stage.add_argument("--dry-run", action="store_true", help="只生成输入，不提交")
    stage.add_argument("--force", action="store_true", help="允许重跑/覆盖已完成阶段")
    _add_common(stage)
    stage.set_defaults(handler=cmd_stage)

    # run ------------------------------------------------------------------- #
    run = subparsers.add_parser(
        "run", help="按 --from/--to 链式推进全部阶段",
        description="逐阶段 prepare→submit→wait→collect；失败默认停止。")
    run.add_argument("stage", nargs="?", metavar="STAGE",
                     help="只跑单个阶段（等价于 --from/--to 同一阶段）")
    run.add_argument("--from", dest="from_stage", metavar="STAGE", help="起始阶段（含）")
    run.add_argument("--to", dest="to_stage", metavar="STAGE", help="结束阶段（含）")
    keep = run.add_mutually_exclusive_group()
    keep.add_argument("--keep-going", action="store_true", help="阶段失败后继续后续阶段")
    keep.add_argument("--fail-fast", action="store_true", help="阶段失败立即停止（默认）")
    run.add_argument("--dry-run", action="store_true", help="只生成输入，不提交")
    run.add_argument("--force", action="store_true", help="重跑已完成阶段")
    _add_overrides(run)
    _add_common(run)
    run.set_defaults(handler=cmd_run)

    # status ---------------------------------------------------------------- #
    status = subparsers.add_parser(
        "status", help="查看 state.json 与阶段状态",
        description="只读；不连远端。")
    status.add_argument("stage", nargs="?", metavar="STAGE", help="只看某个阶段")
    _add_common(status)
    status.set_defaults(handler=cmd_status)

    # report / analyze / review -------------------------------------------- #
    report = subparsers.add_parser(
        "report", help="生成人读报告（report.md）",
        description="调用 smj.report；需要已 collect 的结果。")
    _add_common(report)
    report.add_argument("--force", action="store_true", help="覆盖已存在的报告")
    report.add_argument("--no-checks", action="store_true",
                        help="不计算物理判据（默认会算并写入 manifest.verdict）")
    report.set_defaults(handler=cmd_report)

    analyze = subparsers.add_parser(
        "analyze", help="透射谱 / Landauer 分析",
        description="调用 smj.analysis；读取 aitranss_te 的 TE.dat。")
    _add_common(analyze)
    analyze.add_argument("--force", action="store_true", help="覆盖已有分析结果")
    analyze.set_defaults(handler=cmd_analyze)

    review = subparsers.add_parser(
        "review", help="结构复核（含连接核对）",
        description="调用 smj.review；不通过时退出码为 1。"
                    "连接核对默认自动登记（review.auto_connection_check）；"
                    "手工补签用 --approve --who <名字>（审计/出图"
                    "不通过时退出码为 2）。")
    _add_common(review)
    review.add_argument("--force", action="store_true", help="覆盖已有复核结果")
    review.add_argument("--approve", action="store_true",
                        help="手工（重新）登记连接核对（须配合 --who；审计不通过则拒绝）")
    review.add_argument("--who", default=None, metavar="NAME",
                        help="连接核对签署人（写进 connection_check.json）")
    review.add_argument("--stage", default=None, metavar="DIR",
                        help="electrode_build 目录（缺省自动定位 run 下的 electrode_build）")
    review.add_argument("--png", default=None, metavar="PATH",
                        help="连接核对图输出路径（缺省 <electrode_build>/connection_check_viz.png）")
    review.add_argument("--tcontrol", default=None, metavar="PATH",
                        help="aitranss_te/tcontrol（用于核对表面索引一致性）")
    review.add_argument("--no-figure", action="store_true",
                        help="不出图直接登记（显式声明已看过图；无 matplotlib 时用）")
    review.add_argument("--accept-deviation", type=float, metavar="FRAC",
                        help="放宽连接核对相对偏差门槛（如 0.2）")
    review.add_argument(
        "--connection-check", dest="connection_check",
        action=argparse.BooleanOptionalAction, default=None,
        help="本次复核的连接核对开关：--no-connection-check 关闭（等价 "
             "review.connection_check: false，缺记录/指纹失效只降为提示），"
             "--connection-check 强制开启；缺省沿用 smj.yaml")
    # 更名前的旧开关（隐藏，仍可用）：先叫「表面核对」，更早叫 A2。
    review.add_argument("--surface-check", dest="legacy_surface_check",
                        action=argparse.BooleanOptionalAction, default=None,
                        help=argparse.SUPPRESS)
    review.add_argument("--a2", dest="legacy_a2", action="store_true",
                        default=None, help=argparse.SUPPRESS)
    review.add_argument("--no-a2", dest="legacy_a2", action="store_false",
                        default=None, help=argparse.SUPPRESS)
    review.set_defaults(handler=cmd_review)

    # audit / verify -------------------------------------------------------- #
    audit = subparsers.add_parser(
        "audit", help="本地一致性检查（state / 目录 / 产物）",
        description="只读；不连远端。有 error 时退出码为 1。")
    _add_common(audit)
    audit.set_defaults(handler=cmd_audit)

    verify = subparsers.add_parser(
        "verify", help="安装自检（依赖 / 模板 / 阶段模块）",
        description="检查 Python、PyYAML、numpy、vendored junction、模板与阶段模块。")
    _add_common(verify, with_config=False)
    verify.set_defaults(handler=cmd_verify)

    return parser


def main(argv: list[str] | None = None, *, stage_getter: Callable[[str], Any] | None = None) -> int:
    """CLI 入口；``stage_getter`` 仅测试注入用（默认延迟导入 smj.stages）。"""
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "command", None):
        parser.print_help()
        return EXIT_OK
    args.stage_getter = stage_getter
    handler = getattr(args, "handler", None)
    if handler is None:  # pragma: no cover - 每个子命令都 set_defaults 了
        parser.error("未知子命令: %s" % args.command)
    try:
        return int(handler(args))
    except SystemExit:
        raise
    except BrokenPipeError:  # pragma: no cover - 管道被关闭
        return EXIT_FAILED
    except SmjError as exc:
        _error(str(exc), debug=bool(getattr(args, "debug", False)))
        return EXIT_FAILED


__all__ = ["build_parser", "main", "EXIT_OK", "EXIT_FAILED", "EXIT_UNFINISHED"]


if __name__ == "__main__":  # pragma: no cover - 方便 `python -m smj.cli`
    raise SystemExit(main())
