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
    smj pdos      分组 PDOS（按 atom list 投影到 α/β，与 T(E) 同能量轴）
    smj plot      按统一版式画 T(E) 图（调用 smj.plotting）
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
from .engines import engine_for_config
from .errors import SmjError, TransportError
from .remote import SMJ_STAGE_SCRIPTS
from .resources import hhmmss
from .runner import StageOutcome, StageRunner, lazy_stage_getter
from .state import (
    REMOTE_STAGES,
    STAGE_ORDER,
    THERMAL_STAGE_ORDER,
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
    # 2026-11-22 复审 P1：进度日志必须走 stderr。旧实现把 logger 接成 stdout 的
    # print，`smj run --json` 的首行变成 "[molecule_opt] prepare ..."，
    # json.loads 直接失败——机器可读输出只能有一份 JSON。
    logger = (None if getattr(args, "quiet", False)
              else (lambda message: print(message, file=sys.stderr)))
    return StageRunner(
        config,
        run_dir=getattr(args, "run_dir", None),
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
def _style_version() -> str:
    """绘图版式版本（供运行时提示；plotting 模块顶层不 import matplotlib）。"""
    from .plotting import STYLE_VERSION

    return STYLE_VERSION


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


def cmd_conformers(args: argparse.Namespace) -> int:
    """SMILES → 构象粗筛（ETKDGv3 + MMFF94s + 重原子 RMSD 去重）。"""
    from . import conformers as conformers_mod

    try:
        kept = conformers_mod.select_conformers(
            args.smiles, n_conformers=args.n, seed=args.seed,
            rmsd_tol=args.rmsd_tol, keep=args.keep)
        written = conformers_mod.write_conformers(
            kept, args.dir, out=args.out, smiles=args.smiles,
            n_conformers=args.n, seed=args.seed, rmsd_tol=args.rmsd_tol)
    except (ValueError, RuntimeError) as exc:  # RuntimeError 兜住缺 RDKit
        _error(str(exc), debug=bool(getattr(args, "debug", False)))
        return EXIT_FAILED
    if args.json:
        _emit_json(written["data"])
        return EXIT_OK
    print("构象筛选: %s（ETKDGv3 %d 个候选 → MMFF94s 排序 → 去重保留 %d 个）"
          % (args.smiles, args.n, len(kept)))
    for item in kept:
        print("  %d  %9.3f kcal/mol  RMSD %.3f Å  %2d 原子"
              % (item.rank, item.energy_kcal, item.rmsd_to_first_A, item.atom_count))
    print("已写 %s（rank 1，可直接给 smj init --xyz）" % written["out"])
    print("已写 %s（清单；MMFF 只作粗筛，top-k 建议各跑一次下游链再比较）"
          % written["manifest"])
    return EXIT_OK


def cmd_variants(args: argparse.Namespace) -> int:
    """中性分子 → 硫醇盐 / 金帽对照臂（只生成输入，判定交给下游链）。"""
    from . import variants as variants_mod

    anchors = tuple(part.strip() for part in str(args.anchor).split(",") if part.strip())
    kinds = ("thiolate", "aucap") if args.kind == "both" else (args.kind,)
    try:
        atoms = variants_mod.load_source(args.xyz)
        built = variants_mod.build_variants(
            atoms, charge=args.charge, multiplicity=args.multiplicity,
            kinds=kinds, anchor_elements=anchors, au_distance=args.au_distance)
        written = variants_mod.write_variants(built, args.dir, source=args.xyz)
    except (ValueError, OSError) as exc:
        _error(str(exc), debug=bool(getattr(args, "debug", False)))
        return EXIT_FAILED
    if args.json:
        _emit_json(written["data"])
        return EXIT_OK
    print("对照臂: %s（%d 原子，charge %+d，multiplicity %d）"
          % (Path(args.xyz).name, len(atoms), args.charge, args.multiplicity))
    for entry in written["variants"]:
        print("  %-8s %3d 原子  charge %+d  电子数 %d  移除 H %s"
              % (entry["name"], entry["atom_count"], entry["charge"],
                 entry["electron_count"],
                 entry["removed_h_indices"] or "-"))
        print("           %s" % entry["xyz"])
    print("已写 %s（含 caveat 与 config_hint）" % written["manifest"])
    print("下一步：两条臂各跑一次 `smj stage molecule_opt --run --xyz <臂>.xyz`，"
          "按 variants.json 的 config_hint 设 charge/multiplicity，再比较 C–S 键长、"
          "gap 与 aims_tip_relax 后的锚–Au 接触几何")
    return EXIT_OK


def cmd_species(args: argparse.Namespace) -> int:
    """审计 species 目录：实际档位、基组规模、逐文件 sha256。"""
    from . import species as species_mod

    elements = [part.strip() for part in str(args.elements or "").split(",") if part.strip()]
    try:
        atoms = None
        if args.xyz:
            from .junction.geometry import read_xyz
            atoms = list(read_xyz(args.xyz))
        payload = species_mod.summary(
            args.dir, atoms, elements=elements or None, label=args.label)
        if args.out:
            written = species_mod.write_species_json(payload, args.out)
    except (OSError, ValueError) as exc:
        _error(str(exc), debug=bool(getattr(args, "debug", False)))
        return EXIT_FAILED
    if args.json:
        _emit_json(payload)
        return EXIT_OK
    print("species 目录: %s" % payload["species_dir"])
    print("实际档位: %s（配置标签 %s）"
          % (payload["actual_tier"], payload["declared_label"] or "未声明"))
    print("  %-4s %-8s %-4s %-4s %s" % ("元素", "档位", "基组", "最高l", "sha256[:12]"))
    for item in payload["elements"]:
        print("  %-4s %-8s %-4d %-4s %s"
              % (item["element"], item["tier"], item["basis_functions"],
                 item["max_angular"] or "-", item["sha256"][:12]))
    for warning in payload["warnings"]:
        print("警告: %s" % warning)
    if args.out:
        print("已写 %s" % written)
    return EXIT_OK


def cmd_ensemble(args: argparse.Namespace) -> int:
    """构象自动分支：生成批处理脚本，或汇总各构象 run 的散布。"""
    from . import ensemble as ensemble_mod

    if args.report:
        try:
            payload = ensemble_mod.report_ensemble(args.dir)
        except (OSError, ValueError) as exc:
            _error(str(exc), debug=bool(getattr(args, "debug", False)))
            return EXIT_FAILED
        if args.json:
            _emit_json(payload)
            return EXIT_OK
        print("构象散布: %s（%d 个构象，%d 个已有结果）"
              % (payload["smiles"], payload["kept"], len(payload["rows"])))
        print("  %-4s %-10s %-8s %-10s %-8s" % ("#", "gap/eV", "maxF", "E/kcal", "tier"))
        for row in payload["rows"]:
            print("  %-4d %-10s %-8s %-10s %-8s"
                  % (row["rank"],
                     "%.3f" % row["gap_ev"] if isinstance(row["gap_ev"], (int, float)) else "-",
                     ("%.4f" % row["max_force_ev_per_A"]
                      if isinstance(row["max_force_ev_per_A"], (int, float)) else "-"),
                     ("%.2f" % row["energy_kcal"]
                      if isinstance(row["energy_kcal"], (int, float)) else "-"),
                     row["tier"] or "-"))
        gap = payload["spread"]["gap_ev"]
        if gap is not None:
            print("gap 极差 %.3f eV（这就是构象贡献的量级；<0.05 eV 可认为构象不敏感）" % gap)
        else:
            print("至少需要两个构象的 results.json 才能给极差")
        if payload["missing"]:
            print("尚无结果: %s" % ", ".join(payload["missing"]))
        print("注意: %s" % payload["caveat"])
        return EXIT_OK

    try:
        written = ensemble_mod.write_ensemble(
            args.smiles, args.dir, n_conformers=args.n, keep=args.keep,
            seed=args.seed, rmsd_tol=args.rmsd_tol, run_root=args.run_root,
            config=args.config)
    except (ValueError, RuntimeError) as exc:  # RuntimeError 兜住缺 RDKit
        _error(str(exc), debug=bool(getattr(args, "debug", False)))
        return EXIT_FAILED
    if args.json:
        _emit_json(written["data"])
        return EXIT_OK
    print("构象分支: %s（去重保留 %d 个）" % (args.smiles, len(written["conformers"])))
    for entry in written["conformers"]:
        print("  #%d  %8.3f kcal/mol  %2d 原子  %s"
              % (entry["rank"], entry["energy_kcal"], entry["atom_count"], entry["xyz"]))
    print("已写 %s" % written["manifest"])
    print("已写 %s（每个构象一个隔离 --run-dir；跑完用 `smj ensemble --dir %s --report` 汇总）"
          % (written["script"], args.dir))
    return EXIT_OK


def cmd_plan(args: argparse.Namespace) -> int:
    config = _load_config(args)
    layout, state = _state_for(config)
    rows = []
    for name in STAGE_ORDER:
        meta = _stage_meta(name, getattr(args, "stage_getter", None))
        script = meta["script"] or "-"
        if name == "molecule_opt":
            # 脚本名由引擎声明（缺省引擎 aims），不再在这里写引擎分支。
            script = engine_for_config(config).script
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


def _resolve_retry_slots(spec: str, stage_dir: Path) -> list[int]:
    """把 ``--retry`` 的 SPEC 变成槽位下标列表。

    支持两种写法：显式下标（``17,203``），以及从远端收尾写下的
    ``fc_summary.json`` 里取 ``missing`` / ``not_converged``。后者是主路径：
    collect 的报错信息里给出的就是它可以接受的写法。
    """
    text = str(spec).strip()
    lowered = text.lower()
    if lowered in ("missing", "not-converged", "not_converged"):
        summary_path = Path(stage_dir) / "fc_summary.json"
        if not summary_path.is_file():
            raise SmjError(
                "--retry %s 需要 %s（远端数组作业收尾时写的收敛报告）。"
                "没有它就直接给下标：--retry 17,203"
                % (text, summary_path))
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise SmjError("读 %s 失败: %s" % (summary_path, exc)) from exc
        key = "missing" if lowered == "missing" else "not_converged"
        slots = [int(item) for item in (summary.get(key) or [])]
        if not slots:
            raise SmjError("%s 里 %s 为空：没有需要重跑的槽位" % (summary_path, key))
        return slots
    slots = []
    for part in text.replace(" ", "").split(","):
        if not part:
            continue
        try:
            slots.append(int(part))
        except ValueError:
            raise SmjError(
                "--retry 的下标必须是整数（收到 %r）；也可以写 missing "
                "或 not-converged" % part) from None
    if not slots:
        raise SmjError("--retry 没给下标")
    return slots


def cmd_thermal_compare(args: argparse.Namespace) -> int:
    """多 run 热导比较：先比协议，再比数。"""
    from .phonon.compare import compare, find_runs, load_run, render_markdown

    try:
        paths = find_runs([Path(p) for p in args.paths])
        records = [load_run(path) for path in paths]
    except (OSError, ValueError) as exc:
        print("错误: %s" % exc, file=sys.stderr)
        return EXIT_FAILED
    result = compare(records)
    markdown = render_markdown(result)
    if args.json:
        _emit_json({"consistent": result.consistent,
                    "mismatches": result.mismatches,
                    "rows": result.rows()})
    else:
        if args.out:
            target = Path(args.out)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(markdown, encoding="utf-8", newline="\n")
            print("已写 %s（%d 个 run）" % (target, len(records)))
        else:
            print(markdown)
    if args.strict and not result.consistent:
        return EXIT_FAILED
    return EXIT_OK


def _cmd_thermal_sensitivity(args, config, runner) -> int:
    """扫 η 与能量步长：回答"这个数有多依赖正则化参数"。

    不重跑 DFT —— 力常数已经在 `thermal_fc/force_constants.npy` 里，扫参只是
    把同一份力常数用不同 η 再算一遍透射。没有这一步，横向比较很容易变成
    "η 大的那次算低了"，而两次的差别可能比分子之间的真实差别还大。
    """
    from .stages.thermal_analyze import sensitivity

    fc_dir = runner.layout.ensure_stage("thermal_fc")
    if not (fc_dir / "force_constants.npy").is_file():
        raise SmjError("没找到 thermal_fc/force_constants.npy；先跑完 thermal_fc 再扫参")
    payload = sensitivity(
        fc_dir / "force_constants.npy",
        runner.layout.ensure_stage("thermal_device") / "device_metadata.json",
        fc_dir / "geometry.in",
        config.thermal,
        temperature=config.thermal.resolved_temperature(
            config.transport.temperature))
    target = runner.layout.ensure_stage("thermal_analyze") / "sensitivity.json"
    target.write_text(json.dumps(payload, indent=2, ensure_ascii=False,
                                 allow_nan=False) + "\n", encoding="utf-8")
    if args.json:
        _emit_json(payload)
        return EXIT_OK
    production = payload["production"]
    print("扫参（T = %.1f K，带顶 %.2f meV）：η 与步长都只改正则化，不改力常数"
          % (payload["temperature_K"], payload["band_top_meV"]))
    header = "  %-9s" % "η\\步长"
    for step in payload["energy_steps_meV"]:
        header += "%11s" % ("%.3f" % step)
    print(header)
    table: dict[float, dict[float, float | None]] = {}
    for row in payload["rows"]:
        table.setdefault(row["eta_meV"], {})[row["energy_step_meV"]] = \
            row["g_ph_pw_per_k"]
    for eta in payload["etas_meV"]:
        line = "  %-9.3f" % eta
        for step in payload["energy_steps_meV"]:
            value = table.get(eta, {}).get(step)
            line += "%11s" % ("—" if value is None else "%.3f" % value)
        if eta == production["eta_meV"]:
            line += "   ← 生产配置"
        print(line)
    extrapolated = payload.get("eta_zero_extrapolation")
    if extrapolated:
        print("  η→0 外推：%.3f pW/K（斜率 %.2f pW/K per meV，用 %d 个点）"
              % (extrapolated["g_ph_pw_per_k"],
                 extrapolated["slope_pw_per_k_per_meV"],
                 extrapolated["n_points"]))
    spread = payload.get("step_relative_spread") or {}
    if spread:
        print("  步长引起的相对变化最大 %.2f%%（梯形积分底噪，不是收敛趋势）"
              % (100.0 * max(spread.values())))
    print("  已写 %s" % target)
    return EXIT_OK


def _cmd_thermal_retry(args, config, runner, names) -> int:
    """只重跑 ``thermal_fc`` 数组作业里失败/未收敛的槽位。

    为什么单独一条路径：324 个位移点跑 9 小时，失败一两个点是常态。整批重跑
    既浪费机时也可能踩到同样的资源问题；而 Slurm 的命令行 ``--array=17,203``
    会覆盖脚本里的 ``#SBATCH --array`` 行，所以不必改脚本、不必重新 prepare。
    重跑任务的收尾逻辑按"已完成份数"触发，凑齐后会重写
    ``forces_all.json`` 与 ``fc_summary.json``。
    """
    if any(name != "thermal_fc" for name in names):
        raise SmjError("--retry 只能用于 thermal_fc；加 --stage thermal_fc 或去掉 --stage")
    if config.mode != "slurm-ssh":
        raise SmjError("--retry 只在 slurm-ssh 模式有意义（local 模式会顺序重跑整批）")
    stage_dir = Path(runner.layout.stage_dir("thermal_fc"))
    manifest = stage_dir / "displacements.json"
    if not manifest.is_file():
        raise SmjError("thermal_fc 还没 prepare（缺 %s）" % manifest)
    try:
        count = len(json.loads(manifest.read_text(encoding="utf-8")))
    except (OSError, ValueError) as exc:
        raise SmjError("读 %s 失败: %s" % (manifest, exc)) from exc
    slots = _resolve_retry_slots(args.retry, stage_dir)
    bad = sorted({slot for slot in slots if slot < 0 or slot >= count})
    if bad:
        raise SmjError("重跑下标越界: %s（有效范围 0–%d）" % (bad[:5], count - 1))
    parallel = max(1, int(getattr(config.thermal.fc, "array_parallel", 8) or 8))
    array_spec = "--array=%s%%%d" % (",".join(str(slot) for slot in slots), parallel)
    payload = runner.submit("thermal_fc", force=True, sbatch_args=(array_spec,))
    if args.json:
        _emit_json(payload)
    else:
        print("重跑 %d 个槽位: %s" % (len(slots), array_spec))
        _print_submission(payload, stage="thermal_fc")
        print("跑完后 `smj thermal --stage thermal_fc --collect` 会重拼"
              "forces_all.json / fc_summary.json。")
    return EXIT_OK


def cmd_thermal(args: argparse.Namespace) -> int:
    """声子热导链：``thermal_device`` → ``thermal_relax`` → ``thermal_fc`` →
    ``thermal_analyze``。

    独立于六阶段电子链（方案 D12）：``smj run`` / ``smj plan`` / ``smj status``
    的默认输出里看不到它，只有这一条命令能推进它。
    """
    config = _load_config(args)
    thermal = getattr(config, "thermal", None)
    if thermal is None or not thermal.enabled:
        raise SmjError(
            "thermal 段未启用：请在 smj.yaml 里加 `thermal:\\n  enabled: true`"
            "（并配好热导链要用的 geometry/fc/negf 参数）")

    if args.stage:
        name = normalize_stage(args.stage)
        if name not in THERMAL_STAGE_ORDER:
            raise SmjError("%r 不是热导链阶段；可选: %s"
                           % (args.stage, ", ".join(THERMAL_STAGE_ORDER)))
        names = (name,)
    else:
        names = THERMAL_STAGE_ORDER

    action = "status"
    for candidate in ("prepare", "submit", "wait", "collect", "run"):
        if getattr(args, candidate, False):
            action = candidate
            break
    if args.analyze_only:
        action = "analyze"

    dry_run = bool(args.dry_run) or config.mode == "dry-run"
    if action in ("submit", "wait") and dry_run:
        raise SmjError("dry-run 模式不能 %s；把 run.mode 改成 local / slurm-ssh，"
                       "或去掉 --dry-run" % action)

    runner = _make_runner(args, config)

    if args.sensitivity:
        return _cmd_thermal_sensitivity(args, config, runner)

    if args.retry is not None:
        return _cmd_thermal_retry(args, config, runner, names)

    if action == "status":
        rows = []
        for name in names:
            state = runner.state.stages.get(name)
            rows.append((name, state.status if state else "pending",
                         state.job_id if state else None))
        if args.json:
            _emit_json({"thermal": [
                {"stage": n, "status": s, "job_id": j} for n, s, j in rows]})
        else:
            for name, status, job_id in rows:
                print("  %-16s %-12s %s" % (name, status, job_id or ""))
        return EXIT_OK

    if action == "analyze":
        # 只重算 NEGF + Landauer：力常数已经在盘上，不必再碰集群
        payload = runner.collect("thermal_analyze", force=args.force)
        if args.json:
            _emit_json(payload)
        else:
            _print_collected(payload, stage="thermal_analyze")
        return EXIT_OK

    if action == "prepare":
        for name in names:
            payload = runner.prepare(name, force=args.force)
            if args.json:
                _emit_json(payload)
            else:
                _print_prepared(payload, stage=name)
        return EXIT_OK

    if action == "submit":
        for name in names:
            payload = runner.submit(name, force=args.force)
            if args.json:
                _emit_json(payload)
            else:
                _print_submission(payload, stage=name)
        return EXIT_OK

    if action == "wait":
        for name in names:
            payload = runner.wait(name)
            if args.json:
                _emit_json(payload)
            else:
                _print_wait(payload, stage=name)
            if payload.get("failed"):
                return EXIT_FAILED
            if not payload.get("finished"):
                return EXIT_UNFINISHED
        return EXIT_OK

    if action == "collect":
        for name in names:
            payload = runner.collect(name, force=args.force)
            if args.json:
                _emit_json(payload)
            else:
                _print_collected(payload, stage=name)
        return EXIT_OK

    # --run：逐阶段 prepare → (submit → wait) → collect，失败即停
    from .stages.base import get_stage

    for name in names:
        stage = get_stage(name)
        payload = runner.prepare(name, force=args.force)
        if not args.json:
            _print_prepared(payload, stage=name)
        if stage.needs_submit:
            if dry_run:
                print("  [dry-run] 跳过 %s 的提交" % name)
                continue
            runner.submit(name, force=args.force)
            waited = runner.wait(name)
            if waited.get("failed"):
                return EXIT_FAILED
            if not waited.get("finished"):
                return EXIT_UNFINISHED
        elif not args.json:
            print("  [本地] %s 无需提交，直接 collect" % name)
        collected = runner.collect(name, force=args.force)
        if args.json:
            _emit_json(collected)
        else:
            _print_collected(collected, stage=name)
    return EXIT_OK


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
            if "junction_analyze" in names:
                print("出图: junction_analyze/T_E.png + T_E.pdf（统一版式 style %s）；"
                      "换尺寸/叠加用 smj plot -c %s"
                      % (_style_version(), args.config))
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


def _report_blocked(data: Mapping[str, Any]) -> bool:
    """``smj report --strict``：判定失败/被阻断时不要让 shell 看到 0。"""
    manifest = data.get("manifest")
    verdict = manifest.get("verdict") if isinstance(manifest, Mapping) else None
    if not isinstance(verdict, Mapping):
        return bool(data.get("error"))
    if verdict.get("error"):
        return True
    if verdict.get("physics_verdict") == "failed":
        return True
    if verdict.get("engineering_passed") is False:
        return True
    if int(verdict.get("errors") or 0) > 0:
        return True
    return verdict.get("status") in ("blocked", "failed")


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
    if getattr(args, "strict", False) and _report_blocked(data):
        if not args.json:
            _error("判定未通过：报告已写出，但 --strict 要求退出码非 0")
        return EXIT_FAILED
    return EXIT_OK if ok else EXIT_FAILED


_PDOS_VECTOR_FILES = ("mos.aims", "alpha.aims", "beta.aims")


def _fetch_pdos_vectors(args: argparse.Namespace) -> None:
    """按显式清单整取波函数文件。

    默认回取策略把 ``mos.aims``/``omat.aims`` 归入 head（实测集群侧 0.96 GB），
    本地只留表头；PDOS 要的是完整系数，所以这里显式点名、按全量拉取。
    """
    from .config import load_config
    from .errors import SmjError
    from .runner import StageRunner

    config = load_config(args.config)
    run_dir = Path(args.run)
    runner = StageRunner(
        config, run_dir=run_dir,
        logger=None if args.quiet else (lambda message: print(message, file=sys.stderr)))
    session = runner.session()
    download = session.download("aims_scf", run_dir / "aims_scf",
                                names=_PDOS_VECTOR_FILES)
    fetched = list(download.get("fetched") or [])
    missing = list(download.get("missing") or [])
    if not fetched:
        raise SmjError(
            "没有从远端取回任何波函数文件（缺失: %s）。PDOS 需要完整的 "
            "mos.aims（闭壳层）或 alpha.aims/beta.aims（开壳层）；"
            "确认 run 目录与 remote 配置指向同一作业。" % (", ".join(missing) or "全部"))
    if not args.quiet:
        print("已取回: %s" % ", ".join(fetched))
        if missing:
            print("远端缺失/跳过: %s" % ", ".join(missing))


def cmd_pdos(args: argparse.Namespace) -> int:
    """分组 PDOS：数据不够就明确失败，不拿"只有本征值"的产物凑数。"""
    from .pdos import analyze_pdos

    plot = None if args.plot is None else (True if args.plot == "" else args.plot)
    if args.fetch:
        _fetch_pdos_vectors(args)
    result = analyze_pdos(args.run, groups=args.groups, sigma_ev=args.sigma,
                          te=args.te, output=args.output, plot=plot)
    if args.json:
        _emit_json(result)
        return EXIT_OK
    print("pdos: %s" % result["run_dir"])
    print("  E_F = %r Ha；σ = %.3g eV；网格 %d 点" % (
        result["efermi_Ha"], result["sigma_eV"], len(result["grid_eV"])))
    for group in result["groups"]:
        print("  分组 %-12s 原子 %d 个 / 基函数 %d 个"
              % (group["name"], len(group["atoms"]), group["functions"]))
    for name, channel in result["channels"].items():
        weights = ", ".join("%s=%.4f" % (key, value)
                            for key, value in channel["integrated_weight"].items())
        print("  通道 %-6s %d 个态（%s）；积分权重 %s"
              % (name, channel["states"], channel["spin"], weights))
    for note in result["notes"]:
        print("  注意: %s" % note)
    if result.get("json_path"):
        print("  结果: %s" % result["json_path"])
    if (result.get("plot") or {}).get("path"):
        print("  图: %s" % result["plot"]["path"])
    return EXIT_OK


def cmd_analyze(args: argparse.Namespace) -> int:
    config = _load_config(args)
    module = _optional_module(_AUX_MODULES["analyze"][0], "analyze（透射谱分析）")
    func = _find_callable(module, _AUX_MODULES["analyze"][1])
    plot = None if args.plot is None else (True if args.plot == "" else args.plot)
    result = _invoke(func, config, force=args.force, plot=plot)
    data = _as_mapping(result)
    ok = _truthy_ok(data)
    if args.json:
        _emit_json(data)
    else:
        _print_mapping(data, title="analyze: %s" % config.name)
    return EXIT_OK if ok else EXIT_FAILED


def cmd_plot(args: argparse.Namespace) -> int:
    """统一版式画 T(E)：单条 -> 一张图；多个 input -> 同坐标系叠加。"""
    from .plotting import plot_te, plot_te_multi, resolve_te

    sources = [str(item) for item in (args.inputs or [])]
    if not sources:
        config = _load_config(args)
        sources = [str(config.run_dir)]
    try:
        te_files = [resolve_te(item) for item in sources]
    except FileNotFoundError as exc:
        raise SmjError(str(exc)) from exc

    formats = tuple(f for f in (args.formats or "png,pdf").split(",") if f.strip())
    labels = [item.strip() for item in args.labels.split(",")] if args.labels else None
    if labels is not None and len(labels) != len(te_files):
        # 旧行为：zip 按较短的 labels 截断，剩余曲线被静默丢掉——图看起来"少了一条"，
        # 而命令仍退出 0。叠加图的标签必须与输入一一对应，否则拒绝画。
        raise SmjError("--labels 给了 %d 个标签，但有 %d 条输入；请按输入顺序逐一给出（或省略 --labels）"
                       % (len(labels), len(te_files)))
    out = args.out
    if out is None:
        base = Path(sources[0])
        if not base.is_dir():
            base = te_files[0].parent
        out = base / ("T_E.png" if len(te_files) == 1 else "T_E_overlay.png")
    common = dict(scale=args.scale, figsize=args.figsize, formats=formats, dpi=args.dpi,
                  xlim=args.xlim, ylim=args.ylim, title=args.title, panel=args.panel,
                  spin=args.spin)
    if len(te_files) == 1:
        meta = plot_te(te_files[0], out, label=args.label, annotate=not args.no_annotate,
                       **common)
    else:
        if args.spin and args.spin != "auto":
            # plot_te_multi 没有 spin 参数；旧代码把 --spin 从 common 里悄悄删掉，
            # 叠加图于是只画"总透射"，用户却以为看的是自旋分辨曲线。
            # 默认值 "auto" 等价于"没指定"，不应拦截普通叠加图。
            raise SmjError("叠加模式（多个 input）不支持 --spin：请对单个 TE.dat 单独出图，"
                           "或去掉 --spin 画总透射")
        multi = {k: v for k, v in common.items() if k != "spin"}
        meta = plot_te_multi(list(zip(te_files, labels or [None] * len(te_files))),
                             out, **multi)
    payload = {"paths": meta["paths"], "points": meta.get("points"),
               "scale": meta["scale"], "style_version": meta["style_version"],
               "spin": meta.get("spin"), "channels": meta.get("channels")}
    if args.json:
        _emit_json(payload)
    else:
        print("图已写入: %s" % ", ".join(meta["paths"]))
        print("  scale=%s  style=%s" % (meta["scale"], meta["style_version"]))
        if len(te_files) == 1:
            print("  spin=%s  channels=%s" % (meta.get("spin"), meta.get("channels")))
    return EXIT_OK


def cmd_review(args: argparse.Namespace) -> int:
    config = _load_config(args)

    # 2026-11-15：A2 人工确认更名为「连接核对」（connection check），默认自动登记。
    # --approve 是手工补签/重签路径（fail-closed），签完在同一条命令里看门过没过。
    if getattr(args, "approve", False):
        if not (args.who or "").strip():
            _error("--approve 需要 --who <签署人>（连接核对必须记录是谁看的）")
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
    if args.stage:
        # 2026-11-22 复审 P1：--stage 曾只用于 --approve，普通复核路径完全忽略，
        # 复核错了目录也会写 review.json。这里把显式的 electrode_build 目录透传。
        kwargs["stage"] = args.stage
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


def cmd_spin(args: argparse.Namespace) -> int:
    """S8：只读汇总一份 run 的自旋通道 / Sα-Sβ / Hirshfeld 分区。"""
    from pathlib import Path

    from .spin_report import build_spin_report, render_spin_report

    if getattr(args, "run", None):
        run = Path(args.run)
    else:
        run = _load_config(args).run_dir
    try:
        payload = build_spin_report(run, temperature=args.temperature)
    except FileNotFoundError as exc:
        _error(str(exc), debug=bool(getattr(args, "debug", False)))
        return EXIT_FAILED
    if args.json:
        _emit_json(payload)
        return EXIT_OK
    print(render_spin_report(payload))
    return EXIT_OK


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


def cmd_lead_fc(args: argparse.Namespace) -> int:
    """电极力常数资产：默认打印身份，``--export DIR`` 复制出去。"""
    from .phonon import assets as assets_mod

    try:
        if args.export:
            written = assets_mod.export_lead_fc(args.export)
            if args.json:
                print(json.dumps({"exported": [str(p) for p in written]},
                                 ensure_ascii=False, indent=2))
            else:
                for path in written:
                    print("已导出: %s" % path)
                print("在 smj.yaml 里用 thermal.lead_fc 指向该目录即可替换随包资产。")
            return 0
        info = assets_mod.lead_fc_info()
        if args.json:
            print(json.dumps(info, ensure_ascii=False, indent=2, default=str))
        else:
            for line in assets_mod.describe_lines():
                print(line)
        return 0
    except (FileNotFoundError, ValueError, OSError) as exc:
        print("错误: %s" % exc, file=sys.stderr)
        return 1


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
        from .junction import artifacts, jobwait  # noqa: F401

        check("junction-kernel", True, "jobwait / artifacts 可导入")
    except Exception as exc:  # noqa: BLE001
        check("junction-kernel", False, "junction 内核不可导入: %s" % exc)

    # PROVENANCE.json 台账：记录的 sha256 必须与磁盘一致（改了内核文件忘重算就失败）
    kernel_root = Path(__file__).resolve().parent / "junction"
    try:
        import hashlib as _hashlib
        import json as _json

        _meta = _json.loads((kernel_root / "PROVENANCE.json").read_text(encoding="utf-8"))
        _files = _meta.get("files") or {}
        _stale = []
        for _rel, _info in sorted(_files.items()):
            _path = kernel_root / _rel
            if not _path.is_file():
                _stale.append("%s(缺失)" % _rel)
            elif _hashlib.sha256(_path.read_bytes()).hexdigest() != _info.get("sha256"):
                _stale.append(_rel)
        check("kernel-provenance", not _stale,
              "%d 个文件、%d 条本地补丁" % (len(_files), len(_meta.get("local_patches") or []))
              if not _stale
              else "哈希过期: %s（重跑 scripts/gen_junction_provenance.py）" % ", ".join(_stale[:5]))
    except Exception as exc:  # noqa: BLE001
        check("kernel-provenance", False, "PROVENANCE.json 不可读: %s" % exc)

    templates = kernel_root / "templates"
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
        # 2026-11-22 复审 P2: verify 此前只查 STAGE_ORDER，退役阶段模块（如
        # au_relax）是否还活着完全没人看。它在 state.LEGACY_STAGES 里仍用于
        # 旧 run 的状态迁移，删模块会让旧 run 直接 respec/collect 失败。
        from smj.state import LEGACY_STAGES
        _legacy_missing = []
        for _name in LEGACY_STAGES:
            try:
                _importlib.import_module("smj.stages.%s" % _name)
            except Exception:  # noqa: BLE001
                _legacy_missing.append(_name)
        check("legacy-stages", not _legacy_missing,
              ("已退役阶段模块仍在（旧 run 兼容）: %s" % ", ".join(LEGACY_STAGES))
              if not _legacy_missing
              else "LEGACY_STAGES 列了但模块不可导入: %s" % ", ".join(_legacy_missing),
              severity="warning")
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
    group.add_argument("--run-dir", dest="run_dir", metavar="DIR",
                       help="覆盖 run 目录（构象/对照臂并行跑时各臂互不覆盖产物）")


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

    # conformers（分子前端） -------------------------------------------------- #
    from .conformers import DEFAULT_SEED as _CONF_SEED  # 与库默认值保持单一来源

    conformers = subparsers.add_parser(
        "conformers",
        help="SMILES → 构象粗筛 → molecule.xyz",
        description="ETKDGv3 生成候选 + MMFF94s 排序 + 重原子 RMSD 去重，写 top-k "
                    "构象与 conformers.json。需要可选依赖 RDKit"
                    "（pip install rdkit）。MMFF 只作粗筛，不是电子结构判据。")
    _add_common(conformers, with_config=False)
    conformers.add_argument("--smiles", required=True, metavar="SMILES",
                            help="分子 SMILES（写清 E/Z；只处理裸分子，不带金属）")
    conformers.add_argument("--dir", default=".", metavar="DIR",
                            help="输出目录（默认当前目录）")
    conformers.add_argument("--out", default="molecule.xyz", metavar="NAME",
                            help="rank 1 的文件名（默认 molecule.xyz）")
    conformers.add_argument("--n", type=int, default=20, metavar="N",
                            help="ETKDGv3 候选数（默认 20）")
    conformers.add_argument("--keep", type=int, default=3, metavar="K",
                            help="去重后保留数（默认 3；刚性分子自然少于 K）")
    conformers.add_argument("--seed", type=int, default=_CONF_SEED, metavar="SEED",
                            help="ETKDGv3 种子（默认 %d；0 会生成重复构象，被拒绝）"
                            % _CONF_SEED)
    conformers.add_argument("--rmsd-tol", type=float, default=0.25, metavar="A",
                            help="重原子 RMSD 去重阈值（Å，默认 0.25）")
    conformers.set_defaults(handler=cmd_conformers)

    # variants（化学态对照臂） ------------------------------------------------ #
    variants = subparsers.add_parser(
        "variants",
        help="分子化学态对照臂（硫醇盐 / 金帽）→ xyz + variants.json",
        description="把中性分子派生出两条对照臂：thiolate（去 S–H 的 H，电荷 -1/个）"
                    "与 aucap（S–H 的 H 换成沿原方向的 Au）。只生成输入与清单，"
                    "判定仍由 molecule_opt → aims_tip_relax 跑出来；气相/裸阴离子的"
                    "口径边界写进 variants.json 的 caveat。")
    _add_common(variants, with_config=False)
    variants.add_argument("--xyz", required=True, metavar="FILE",
                          help="中性分子几何（molecule.xyz / opt.xyz）")
    variants.add_argument("--dir", default="variants", metavar="DIR",
                          help="输出目录（默认 variants/）")
    variants.add_argument("--kind", default="both",
                          choices=("thiolate", "aucap", "both"),
                          help="生成哪条臂（默认 both）")
    variants.add_argument("--charge", type=int, default=0, metavar="Q",
                          help="输入分子的总电荷（默认 0）")
    variants.add_argument("--multiplicity", type=int, default=1, metavar="M",
                          help="输入分子的多重度（默认 1；电子数不变，两条臂沿用）")
    variants.add_argument("--anchor", default="S", metavar="ELEM",
                          help="锚定元素（默认 S；Se 等用逗号分隔）")
    variants.add_argument("--au-distance", type=float, default=2.40, metavar="A",
                          help="金帽沿原 S–H 方向的距离（Å，默认 2.40）")
    variants.set_defaults(handler=cmd_variants)

    # species（基组口径审计） ------------------------------------------------ #
    species = subparsers.add_parser(
        "species",
        help="审计 species 目录（实际档位 / 基组规模 / sha256）",
        description="molecule.opt_basis 只是配置标签，真正生效的是 species 文件。本命令读"
                    "目录里的 *_default 文件，报告文件头声明的档位、生效基组行数与"
                    "最高角动量、逐文件 sha256，并可导出 species.json 作溯源。")
    _add_common(species, with_config=False)
    species.add_argument("--dir", required=True, metavar="DIR",
                         help="species 目录（cluster.species_dir 或 paths.species_defaults）")
    species.add_argument("--xyz", metavar="FILE",
                         help="只审计这个分子用到的元素（省掉 --elements）")
    species.add_argument("--elements", metavar="LIST",
                         help="逗号分隔的元素（如 S,C,H,Au）")
    species.add_argument("--label", metavar="TIER",
                         help="配置里的标签（如 light）；与文件头不符时给警告")
    species.add_argument("--out", metavar="FILE",
                         help="把审计结果写成 JSON（建议 species.json）")
    species.set_defaults(handler=cmd_species)

    # ensemble（构象自动分支） ----------------------------------------------- #
    ensemble = subparsers.add_parser(
        "ensemble",
        help="构象自动分支：k 个构象各一个 run 目录，或汇总散布",
        description="生成 top-k 构象 + ensemble.json + run_ensemble.sh（每条命令用"
                    " --run-dir 隔离）；跑完用 `smj ensemble --dir <dir> --report` 汇总"
                    " gap / 最大力 / 构象能量的极差。MMFF94s 只作粗筛。")
    _add_common(ensemble)
    ensemble.add_argument("--dir", default="ensemble", metavar="DIR",
                          help="工作目录（默认 ensemble/）")
    ensemble.add_argument("--smiles", metavar="SMILES",
                          help="分子 SMILES（生成模式必填）")
    ensemble.add_argument("--n", type=int, default=20, metavar="N",
                          help="ETKDGv3 候选数（默认 20）")
    ensemble.add_argument("--keep", type=int, default=3, metavar="K",
                          help="去重后保留数（默认 3）")
    ensemble.add_argument("--seed", type=int, default=42, metavar="SEED",
                          help="ETKDGv3 种子（默认 42；0 会生成重复构象，被拒绝）")
    ensemble.add_argument("--rmsd-tol", type=float, default=0.25, metavar="A",
                          help="重原子 RMSD 去重阈值（Å，默认 0.25）")
    ensemble.add_argument("--run-root", default="run", metavar="SUBDIR",
                          help="各构象 run 目录的父目录（默认 run/）")
    ensemble.add_argument("--report", action="store_true",
                          help="只汇总已有 run 的 results.json（不生成）")
    ensemble.set_defaults(handler=cmd_ensemble)

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

    # thermal --------------------------------------------------------------- #
    thermal_cmd = subparsers.add_parser(
        "thermal", help="声子热导链：器件 → 弛豫 → 力常数 → 热导（独立于六阶段）",
        description=(
            "推进 thermal_device → thermal_relax → thermal_fc → thermal_analyze。"
            "不带动作 flag 时只打印这四个阶段的状态。"
            "需要 smj.yaml 里 thermal.enabled: true。"))
    thermal_cmd.add_argument("--stage", metavar="STAGE",
                             help="只操作一个阶段: %s"
                             % ", ".join(THERMAL_STAGE_ORDER))
    group = thermal_cmd.add_mutually_exclusive_group()
    group.add_argument("--prepare", action="store_true", help="只生成输入，不提交")
    group.add_argument("--submit", action="store_true", help="只提交")
    group.add_argument("--wait", action="store_true", help="只等待")
    group.add_argument("--collect", action="store_true", help="只解析产物")
    group.add_argument("--run", action="store_true", help="prepare→submit→wait→collect")
    group.add_argument("--analyze-only", dest="analyze_only", action="store_true",
                       help="只重算 NEGF/Landauer（力常数已在盘上，不碰集群）")
    group.add_argument("--retry", metavar="SPEC", default=None,
                       help=("只重跑 thermal_fc 数组作业里失败/未收敛的槽位。"
                             "SPEC 是逗号分隔的下标（如 17,203），"
                             "或 missing / not-converged（从 fc_summary.json 读）。"
                             "重跑后收尾步骤会自动重拼 forces_all.json"))
    group.add_argument("--sensitivity", action="store_true",
                       help=("用已有 force_constants.npy 扫 η 与能量步长"
                             "（不碰集群），写 sensitivity.json 并给出 η→0 外推"))
    thermal_cmd.add_argument("--dry-run", action="store_true", help="只生成输入")
    thermal_cmd.add_argument("--force", action="store_true", help="允许重跑/覆盖")
    _add_common(thermal_cmd)
    thermal_cmd.set_defaults(handler=cmd_thermal)

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
    report.add_argument("--strict", action="store_true",
                        help="判定失败/被阻断时退出码为 1（默认只要报告写成功就是 0）")
    report.set_defaults(handler=cmd_report)

    analyze = subparsers.add_parser(
        "analyze", help="透射谱 / Landauer 分析",
        description="调用 smj.analysis；读取 aitranss_te 的 TE.dat。")
    _add_common(analyze)
    analyze.add_argument("--force", action="store_true", help="覆盖已有分析结果")
    analyze.add_argument("--plot", nargs="?", const="", default=None, metavar="PATH",
                         help="同时出 T(E) 图（缺省写到 TE.dat 同目录 T_E_plot.png）")
    analyze.set_defaults(handler=cmd_analyze)

    pdos = subparsers.add_parser(
        "pdos", help="分组 PDOS（按 atom list 投影到 α/β，与 T(E) 同能量轴）",
        description="读 aims_scf 的 mos.aims/alpha.aims/beta.aims 与 basis-indices.out，"
                    "按原子分组投影到 TE.dat 的 E − E_F 网格；缺完整本征向量会显式失败。")
    _add_common(pdos, with_config=True)
    pdos.add_argument("--fetch", action="store_true",
                      help="先按显式清单从远端整取 mos.aims/alpha.aims/beta.aims"
                           "（默认回取策略只抓表头，PDOS 需要完整本征向量）")
    pdos.add_argument("run", metavar="RUN",
                      help="run 目录（含 aims_scf/ 与 aitranss_te/TE.dat）")
    pdos.add_argument("--groups", required=True, metavar="SPEC",
                      help="原子分组，例如 molecule=1-20,left=21-52,right=53-84（1-based）")
    pdos.add_argument("--sigma", type=float, default=0.1, metavar="EV",
                      help="高斯展宽 σ（eV，默认 0.1）")
    pdos.add_argument("--te", metavar="PATH", help="TE.dat 覆盖（默认 <run>/aitranss_te/TE.dat）")
    pdos.add_argument("-o", "--output", metavar="PATH",
                      help="把结果写成 JSON（默认只打印摘要）")
    pdos.add_argument("--plot", nargs="?", const="", default=None, metavar="PATH",
                      help="出 PDOS+T(E) 图（缺省写到 <run>/PDOS.png）")
    pdos.set_defaults(handler=cmd_pdos)

    plot = subparsers.add_parser(
        "plot", help="按统一版式画 T(E) 图",
        description="统一版式（文献式）：四边框、刻度朝内、E_F 虚线、log 纵轴、"
                    "panel 角标。可给 TE.dat / run 目录 / 多个输入叠加。")
    _add_common(plot)
    plot.add_argument("inputs", nargs="*", metavar="PATH",
                      help="TE.dat、run 目录或阶段目录（缺省用 -c 配置的 run）")
    plot.add_argument("-o", "--out", metavar="PATH",
                      help="输出路径（缺省 <run 目录>/T_E.png；输入为 TE.dat 时用它所在目录）")
    plot.add_argument("--scale", choices=("log", "linear"), default="log",
                      help="纵轴刻度（默认 log）")
    plot.add_argument("--figsize", choices=("single", "double", "wide"),
                      default="single", help="成图尺寸（默认 single，单栏）")
    plot.add_argument("--formats", default="png,pdf", metavar="LIST",
                      help="导出格式，逗号分隔（默认 png,pdf）")
    plot.add_argument("--dpi", type=float, default=600, help="位图分辨率（默认 600）")
    plot.add_argument("--xlim", nargs=2, type=float, metavar=("LO", "HI"), help="x 轴范围（eV）")
    plot.add_argument("--ylim", nargs=2, type=float, metavar=("LO", "HI"), help="y 轴范围")
    plot.add_argument("--title", default=None, help="图标题（默认无，符合期刊图注习惯）")
    plot.add_argument("--panel", default=None, metavar="LETTER",
                      help="panel 角标字母，如 a")
    plot.add_argument("--label", default=None,
                      help="单条曲线的图例标签（不给则无图例）")
    plot.add_argument("--labels", default=None, metavar="A,B,C",
                      help="多条曲线叠加时的图例标签（逗号分隔）")
    plot.add_argument("--spin", choices=("auto", "mean", "channels", "polarization"),
                      default="auto",
                      help="自旋通道画法（默认 auto：TE.dat 有 α/β 两列时出双 panel）")
    plot.add_argument("--no-annotate", action="store_true",
                      help="不标注 T(E_F) / G/G0")
    plot.set_defaults(handler=cmd_plot)

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

    # spin ---------------------------------------------------------------- #
    spin = subparsers.add_parser(
        "spin", help="开壳层自旋汇总（通道 / Sα-Sβ / Hirshfeld 分区）",
        description="只读：汇总 run_spec、TE.dat 的双通道积分、"
                    "aims_scf/spin_partition.json 与 review 的自旋提示；"
                    "不写任何文件（需要新结论请跑 smj review）。")
    _add_common(spin)
    spin.add_argument("run", nargs="?", default=None, metavar="RUN",
                      help="run 目录（缺省用配置里的 run_dir）")
    spin.add_argument("--temperature", type=float, default=None, metavar="K",
                      help="覆盖积分温度(K)；缺省取 run_spec/transport.temperature")
    spin.set_defaults(handler=cmd_spin)

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

    lead_fc = subparsers.add_parser(
        "lead-fc", help="查看/导出声子热导用的电极力常数资产",
        description=("声子 NEGF 需要一份 Au(111) 半无限电极的力常数。它由一次"
                     "体相有限位移得到（3x3x6 层、Gamma、54 个 SCF），与分子无关，"
                     "所以随包发布一份。--show 看身份，--export 导出后可用 "
                     "thermal.lead_fc 指向自己的版本。"))
    _add_common(lead_fc, with_config=False)
    lead_fc.add_argument("--show", action="store_true", help="打印资产身份（默认动作）")
    lead_fc.add_argument("--export", metavar="DIR", default=None,
                         help="把资产复制到 DIR（生成 lead_fc.npy + lead_fc.json）")
    lead_fc.set_defaults(handler=cmd_lead_fc)

    compare_cmd = subparsers.add_parser(
        "thermal-compare", help="把多个 run 的声子热导摆到一张表里",
        description=("读各 run 的 thermal_analyze/thermal.json，输出 markdown 表。"
                     "**先比协议再比数**：η、能量步长、温度、电极力常数任一不一致"
                     "都会显式列出（D15：跨协议比较等于换方法）；--strict 时直接失败。"))
    _add_common(compare_cmd, with_config=False)
    compare_cmd.add_argument("paths", nargs="+", metavar="PATH",
                             help="run 目录、其父目录（扫一层）或 thermal.json")
    compare_cmd.add_argument("--out", metavar="FILE", default=None,
                             help="把 markdown 表写到文件（默认打印到 stdout）")
    compare_cmd.add_argument("--strict", action="store_true",
                             help="协议不一致时返回非零退出码")
    compare_cmd.set_defaults(handler=cmd_thermal_compare)

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
    except Exception as exc:  # noqa: BLE001 - CLI 顶层兜底：给一行可读错误，不喷 traceback
        if getattr(args, "json", False):
            _emit_json({"ok": False, "error": "%s: %s" % (type(exc).__name__, exc)})
        _error("%s: %s（加 --debug 看完整 traceback）"
               % (type(exc).__name__, exc),
               debug=bool(getattr(args, "debug", False)))
        return EXIT_FAILED


__all__ = ["build_parser", "main", "EXIT_OK", "EXIT_FAILED", "EXIT_UNFINISHED"]


if __name__ == "__main__":  # pragma: no cover - 方便 `python -m smj.cli`
    raise SystemExit(main())
