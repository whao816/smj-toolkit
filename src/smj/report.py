"""run 报告：从 run 目录生成 ``manifest.json`` + 人读 ``report.md``。

输入是三样**已有事实**（不臆造）：

- ``<run>/state.json`` —— 阶段状态机（:class:`smj.state.RunState`）；
- ``<run>/<stage>/`` —— 各阶段产物目录（含 ``auN_build`` 变体目录）；
- 可选的既有 ``manifest.json`` / ``report.md`` —— 只读参考格式，**不覆盖**
  除非显式 ``force=True``。

报告固定包含：概览、阶段状态表、输入摘要、资源估算、输运结果、失败原因、产物
清单。与 vendored ``junction.manifest`` 的差异：本模块面向**独立 smj 包**的
``state.json`` 驱动流程，额外给出阶段状态与资源估算；产物清单/判据表沿用
既有格式，便于人读与机械解析（:func:`parse_report`）。

命名与 CLI 约定
---------------
``generate_report(config_or_run_dir, *, force=False)`` 是首选入口；``run_dir``
还可以直接是 :class:`smj.config.RunConfig`（取 ``config.run_dir``）。
``force=False`` 时**逐文件**跳过已存在的 ``manifest.json`` / ``report.md``，
并在返回值 ``skipped`` 里说明原因。
"""

from __future__ import annotations

import hashlib
import json
import os
import re
from pathlib import Path
from typing import Any, Iterable, Mapping

from .provenance import (
    load_run_spec,
    run_spec_payload,
    run_spec_sha256,
)
from .resources import check_smj_budget, estimate_smj_resources
from .state import STAGE_ORDER, RunLayout, RunState, now_iso

REPORT_SCHEMA = "smj.report.v1"

#: 单文件哈希上限（超过且未被策略标为 full 的只记大小，避免为 GB 级产物白算校验）。
_HASH_MAX_BYTES = 8 * 1024 * 1024

#: 已知阶段（用于阶段状态表；``auN_build`` 变体由目录发现另行补入）。
_BUILD_STAGE_RE = re.compile(r"^au\d+_build$")

#: 阶段默认几何产物（用于资源估算的原子数）。
_STAGE_GEOMETRY = {
    "gaussian_opt": ("opt.xyz", "molecule.xyz", "geometry.in"),
    "aims_tip_relax": ("geometry.in.next_step", "geometry.in"),
    "aims_scf": ("geometry.in",),
    "aitranss_te": ("geometry.in",),
}

_SECTION_ORDER = (
    "0. 概览",
    "1. 阶段状态",
    "2. 输入摘要",
    "3. 资源估算",
    "4. 输运结果",
    "5. 失败原因",
    "6. 产物清单",
)


# --------------------------------------------------------------------------- #
# 基础工具
# --------------------------------------------------------------------------- #
def _write_text(path: str | Path, text: str) -> Path:
    """原子写文本（``.tmp`` -> ``os.replace``），父目录自动创建。"""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_name(target.name + ".tmp")
    tmp.write_text(text, encoding="utf-8", newline="\n")
    os.replace(tmp, target)
    return target


def _write_json(payload: Mapping[str, Any], path: str | Path) -> Path:
    return _write_text(path, json.dumps(dict(payload), indent=2, ensure_ascii=False) + "\n")


def _run_dir_of(target: Any) -> tuple[Path, Any]:
    """接受 RunConfig 或路径 -> (run_dir, config|None)。"""
    if hasattr(target, "run_dir") and not isinstance(target, (str, Path)):
        return Path(target.run_dir), target
    return Path(target), None


def _sha256(path: Path, chunk: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(chunk)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def _retrieval_mode(stage: str, name: str) -> str:
    """vendored artifacts 的取回策略（缺失时 unknown，不臆断）。"""
    try:
        from smj._vendor.junction.artifacts import retrieval_mode  # lazy: vendored
    except Exception:  # noqa: BLE001
        return "unknown"
    try:
        return str(retrieval_mode(stage, name))
    except Exception:  # noqa: BLE001
        return "unknown"


def _atom_count(path: Path, *, xyz: bool = False) -> int | None:
    """读结构文件的原子数：失败返回 None（不猜）。"""
    if not path.is_file():
        return None
    try:
        if xyz or path.suffix.lower() == ".xyz":
            from smj._vendor.junction.geometry import read_xyz  # lazy: vendored
            return len(read_xyz(path))
        from smj._vendor.junction.geometry import read_aims_geometry  # lazy: vendored
        return len(read_aims_geometry(path))
    except Exception:  # noqa: BLE001
        return None


def _read_config_summary(run_dir: Path, warnings: list[str]) -> dict:
    """若 run 目录有 smj.yaml，读取模式/温度/电极摘要（失败只告警）。"""
    summary: dict[str, Any] = {}
    cfg_path = run_dir / "smj.yaml"
    if not cfg_path.is_file():
        return summary
    try:
        from .config import load_config
        cfg = load_config(cfg_path)
    except Exception as exc:  # noqa: BLE001 —— 旧 run 的配置可能不再合法
        warnings.append("smj.yaml 不可解析: %s" % exc)
        return summary
    summary.update({"mode": cfg.mode, "name": cfg.name,
                    "temperature_K": getattr(cfg.transport, "temperature", None),
                    "au_per_electrode": getattr(cfg.electrode, "au_per_electrode", None),
                    "source": str(cfg_path)})
    for warning in cfg.warnings:
        warnings.append("config: %s" % warning)
    return summary


# --------------------------------------------------------------------------- #
# 阶段发现 / 产物清单
# --------------------------------------------------------------------------- #
def _build_number(name: str) -> int:
    match = re.match(r"^au(\d+)_build$", name)
    return int(match.group(1)) if match else 10 ** 9


def discover_stages(run_dir: str | Path) -> list[str]:
    """按 run 目录**实际存在**的阶段子目录，按执行顺序返回阶段名。

    2026-11 P0-2：阶段清单的唯一事实源是 :data:`smj.state.STAGE_ORDER`。
    旧实现调用 vendored ``discover_smj_stages``，而那份清单缺 ``au_relax``，
    导致 demo 的 manifest/report 里阶段数少一个、A2 对象阶段显示为"—"。
    ``auN_build`` 变体目录仍作为扩展追加在 ``au60_build`` 之后。
    """
    root = Path(run_dir)
    if not root.is_dir():
        return []
    ordered = [stage for stage in STAGE_ORDER if (root / stage).is_dir()]
    extras = sorted(
        (child.name for child in root.iterdir()
         if child.is_dir() and _BUILD_STAGE_RE.match(child.name)
         and child.name not in STAGE_ORDER),
        key=lambda name: (_build_number(name), name))
    if not extras:
        return ordered
    if "au60_build" in ordered:
        index = ordered.index("au60_build") + 1
    elif "au_relax" in ordered:
        index = ordered.index("au_relax")
    elif "aims_tip_relax" in ordered:
        index = ordered.index("aims_tip_relax") + 1
    else:
        index = 0
    ordered[index:index] = extras
    return ordered


def _inventory(run_dir: Path, stages: Iterable[str]) -> dict[str, dict]:
    inventory: dict[str, dict] = {}
    for stage in stages:
        stage_dir = run_dir / stage
        if not stage_dir.is_dir():
            continue
        files = []
        total = 0
        for path in sorted(p for p in stage_dir.iterdir() if p.is_file()):
            size = path.stat().st_size
            total += size
            mode = _retrieval_mode(stage, path.name)
            entry: dict[str, Any] = {"name": path.name, "bytes": size, "retrieval": mode}
            if mode == "full" or size <= _HASH_MAX_BYTES:
                entry["sha256"] = _sha256(path)
            else:
                entry["sha256"] = None
                entry["note"] = ("超过 %d MiB 且 retrieval=%s，按策略免哈希"
                                 % (_HASH_MAX_BYTES // (1024 * 1024), mode))
            files.append(entry)
        inventory[stage] = {"dir": str(stage_dir), "n_files": len(files),
                            "bytes": total, "files": files}
    return inventory


# --------------------------------------------------------------------------- #
# 输入 / 资源 / 输运 / 失败
# --------------------------------------------------------------------------- #
def _find_build_dir(run_dir: Path) -> Path | None:
    # 0.1.2：报告里的"电极"几何优先用 au_relax（背层固定弛豫后、真正送 SCF 的那份）。
    relaxed = run_dir / "au_relax"
    if (relaxed / "geometry.in").is_file() and (relaxed / "electrode_metadata.json").is_file():
        return relaxed
    for child in sorted(run_dir.iterdir()):
        if child.is_dir() and _BUILD_STAGE_RE.match(child.name):
            return child
    return None


def _input_summary(run_dir: Path, spec: Mapping[str, Any],
                   warnings: list[str], *,
                   run_spec: Mapping[str, Any] | None = None) -> dict:
    spec_config = dict((run_spec or {}).get("config") or {})
    run_cfg = dict(spec_config.get("run") or {})
    transport_cfg = dict(spec_config.get("transport") or {})
    if run_spec is not None:
        # 2026-11 P0-1：mode/温度以 run_spec 为准；旧实现读 run_dir/smj.yaml，
        # 而 demo 的 yaml 在 run 的父目录，于是 mode/temperature 都是 null。
        mode = run_cfg.get("mode")
        temperature = transport_cfg.get("temperature")
    else:
        cfg = _read_config_summary(run_dir, warnings)
        mode = cfg.get("mode")
        temperature = cfg.get("temperature_K", spec.get("temperature_K"))
    effective_spec = dict(spec or {}) or spec_config
    inputs: dict[str, Any] = {
        "run_id": run_dir.name,
        "mode": mode,
        "temperature_K": temperature,
        "spec": effective_spec,
        "molecule": None,
        "electrode": None,
        "te": None,
    }
    # 分子结构：优先 gaussian_opt/opt.xyz（Gaussian 优化产物），否则 run 根的 xyz。
    molecule_xyz = run_dir / "gaussian_opt" / "opt.xyz"
    if not molecule_xyz.is_file():
        candidates = sorted(run_dir.glob("*.xyz"))
        molecule_xyz = candidates[0] if candidates else molecule_xyz
    if molecule_xyz.is_file():
        inputs["molecule"] = {
            "path": str(molecule_xyz),
            "atom_count": _atom_count(molecule_xyz, xyz=True),
        }
        # run_spec 记的是**输入**分子（优化前）的哈希，保留为独立字段，避免混淆。
        molecule_meta = (run_spec or {}).get("molecule") or {}
        if molecule_meta.get("sha256"):
            inputs["molecule"]["run_input_sha256"] = molecule_meta["sha256"]

    # 电极：auN_build/electrode_metadata.json 的既有键（不补造）。
    build = _find_build_dir(run_dir)
    if build is not None:
        metadata_path = build / "electrode_metadata.json"
        meta: dict[str, Any] = {}
        if metadata_path.is_file():
            try:
                meta = json.loads(metadata_path.read_text(encoding="utf-8"))
            except (OSError, ValueError) as exc:
                warnings.append("electrode_metadata.json 不可解析: %s" % exc)
        geometry = build / "geometry.in"
        inputs["electrode"] = {
            "stage": build.name,
            "geometry": str(geometry) if geometry.is_file() else None,
            "atom_count": _atom_count(geometry) if geometry.is_file() else None,
            "au_per_electrode": meta.get("au_per_electrode"),
            "left_surface_indices": meta.get("left_surface_indices"),
            "right_surface_indices": meta.get("right_surface_indices"),
            "anchor_s_au_distance_A": meta.get("anchor_s_au_distance_A"),
            "surface_fingerprint": meta.get("surface_fingerprint"),
        }

    te = run_dir / "aitranss_te" / "TE.dat"
    if te.is_file():
        try:
            from .analysis import te_header
            inputs["te"] = te_header(te)
        except Exception as exc:  # noqa: BLE001
            warnings.append("TE.dat 表头解析失败: %s" % exc)
            inputs["te"] = {"source": str(te)}
    return inputs


def _stage_atoms(run_dir: Path, stage: str) -> int | None:
    for name in _STAGE_GEOMETRY.get(stage, ()):
        path = run_dir / stage / name
        if path.is_file():
            return _atom_count(path, xyz=name.endswith(".xyz"))
    # aims_scf / aitranss_te 的几何落在上游 auN_build。
    if stage in ("aims_scf", "aitranss_te"):
        build = _find_build_dir(run_dir)
        if build is not None:
            return _atom_count(build / "geometry.in")
    return None


def _resource_summary(run_dir: Path, state: RunState,
                      warnings: list[str]) -> dict[str, dict]:
    summary: dict[str, dict] = {}
    for stage in ("gaussian_opt", "aims_tip_relax", "au_relax", "aims_scf", "aitranss_te"):
        recorded = dict(state.stage(stage).resources or {})
        if recorded:
            summary[stage] = {"source": "state", **recorded,
                              "budget": check_smj_budget(recorded)}
            continue
        atoms = _stage_atoms(run_dir, stage)
        if not atoms:
            continue
        estimate = estimate_smj_resources(stage, atoms)
        if estimate.get("unknown"):
            warnings.append("资源估算失败（%s）: %s" % (stage, estimate.get("basis")))
            continue
        summary[stage] = {"source": "estimate", "atoms": atoms, **estimate,
                          "budget": check_smj_budget(estimate)}
    return summary


_TRANSPORT_KEYS = (
    "points", "temperature_K", "L0", "L1", "L2", "log10_G_over_G0_like",
    "S_uV_per_K", "k_e_pW_per_K", "ZT_electronic", "fermi_nearest",
    "max_transmission",
)


def _transport_summary(run_dir: Path, state: RunState, temperature: float | None,
                       warnings: list[str]) -> dict | None:
    analysis_path = run_dir / "junction_analyze" / "analysis.json"
    if analysis_path.is_file():
        try:
            payload = json.loads(analysis_path.read_text(encoding="utf-8"))
            return {"source": str(analysis_path),
                    **{k: payload.get(k) for k in _TRANSPORT_KEYS}}
        except (OSError, ValueError) as exc:
            warnings.append("analysis.json 不可解析: %s" % exc)
    te = run_dir / "aitranss_te" / "TE.dat"
    if not te.is_file():
        return None
    try:
        from .analysis import analyze_te
        result = analyze_te(te, temperature=float(temperature or 300.0),
                            independent=False)
    except Exception as exc:  # noqa: BLE001 —— 报告不因输运分析失败而崩
        warnings.append("输运分析失败: %s" % exc)
        return None
    payload = {"source": str(te), **{k: result.get(k) for k in _TRANSPORT_KEYS}}
    recorded = dict(state.stage("aitranss_te").results or {})
    if recorded:
        payload["state_results"] = recorded
    return payload


def _failure_summary(state: RunState) -> list[dict]:
    failures = []
    for stage in STAGE_ORDER:
        item = state.stage(stage)
        if item.status != "failed":
            continue
        failures.append({
            "stage": stage,
            "status": item.status,
            "error": item.error,
            "note": item.note,
            "exit_code": item.exit_code,
            "job_id": item.job_id,
            "finished_at": item.finished_at,
        })
    return failures


# --------------------------------------------------------------------------- #
# manifest / report
# --------------------------------------------------------------------------- #
def build_manifest(target: Any, *, run_id: str | None = None,
                   spec: Mapping[str, Any] | None = None,
                   state: RunState | None = None,
                   stages: Iterable[str] | None = None,
                   include_checks: bool = False,
                   system: str | None = None) -> dict:
    """汇总 run 目录 -> manifest dict（不写文件）。"""
    run_dir, config = _run_dir_of(target)
    if system is None and config is not None:
        system = getattr(getattr(config, "molecule", None), "family", None)
    layout = RunLayout(run_dir)
    rid = run_id or (getattr(config, "name", None) if config is not None else None) or run_dir.name
    state = state or RunState.load(layout.state_path, run=rid)
    warnings: list[str] = list(state.warnings)

    # 2026-11 P0-1：manifest 绑定 run_spec（配置快照 + 文件哈希）。没有 run_spec
    # 的旧 run 允许用当前配置重建，但必须显式标 reconstructed 并告警，不能假装
    # 这份 manifest 能证明原始参数。
    run_spec = load_run_spec(run_dir)
    provenance: dict[str, Any] = {
        "path": None,
        "run_spec_sha256": None,
        "reconstructed": False,
        "protocol_sha256": None,
    }
    if run_spec is None:
        if config is not None:
            run_spec = run_spec_payload(config, run_dir, reconstructed=True)
            provenance["reconstructed"] = True
            warnings.append(
                "run_spec.json 缺失：manifest 用当前配置重建（reconstructed=true），"
                "不能证明与原始 run 的参数完全一致")
        else:
            warnings.append(
                "run_spec.json 缺失：run 目录里没有可加载的配置，manifest 无法绑定运行参数")
    else:
        provenance["path"] = str(run_dir / "run_spec.json")
        provenance["run_spec_sha256"] = run_spec_sha256(run_dir)
        provenance["reconstructed"] = bool(run_spec.get("reconstructed"))
        if provenance["reconstructed"]:
            warnings.append("run_spec.json 标记为 reconstructed：参数来自事后重建，不是原始快照")
    if run_spec is not None:
        provenance["protocol_sha256"] = run_spec.get("protocol_sha256")

    spec_data = dict(spec or {}) or dict((run_spec or {}).get("config") or {})
    input_summary = _input_summary(run_dir, spec_data, warnings, run_spec=run_spec)
    inventory = _inventory(run_dir, stages if stages is not None
                           else discover_stages(run_dir))
    resources = _resource_summary(run_dir, state, warnings)
    transport = _transport_summary(run_dir, state,
                                   input_summary.get("temperature_K"), warnings)
    failures = _failure_summary(state)

    totals = {
        "bytes": sum(item["bytes"] for item in inventory.values()),
        "gib": round(sum(item["bytes"] for item in inventory.values()) / (1024 ** 3), 3),
        "stages": len(inventory),
    }
    manifest: dict[str, Any] = {
        "schema": REPORT_SCHEMA,
        "run_id": rid,
        "root": str(run_dir),
        "generated_at": now_iso(),
        "spec": spec_data,
        "run_spec": run_spec,
        "provenance": provenance,
        "state": {
            "created_at": state.created_at,
            "updated_at": state.updated_at,
            "stages": state.statuses(),
            "details": {
                stage: {
                    "job_id": state.stage(stage).job_id,
                    "exit_code": state.stage(stage).exit_code,
                    "elapsed_s": state.stage(stage).elapsed_s,
                    "note": state.stage(stage).note,
                    "error": state.stage(stage).error,
                }
                for stage in STAGE_ORDER
            },
            "first_incomplete": state.first_incomplete(),
            "is_done": all(state.is_done(s) for s in STAGE_ORDER),
            "warnings": list(state.warnings),
        },
        "inputs": input_summary,
        "stages": inventory,
        "resources": resources,
        "transport": transport,
        "failures": failures,
        "warnings": warnings,
        "verdict": None,
    }

    if include_checks:
        # 2026-11-11 本地补丁 (smj F2): 不再用 `aitranss_te/` 目录是否存在当门。
        # 空目录曾经直接判 passed; 现在统一交给 validate_junction_run, 缺 TE.dat 时
        # criteria 会给出 error 级 te_present, 结果是 failed 而不是"支持 H0"。
        try:
            from smj._vendor.junction.validate import validate_junction_run  # lazy: vendored
            verdict = validate_junction_run(
                run_dir, results={"junction": transport or {}}, system=system,
                temperature=input_summary.get("temperature_K"))
            manifest["verdict"] = {
                key: verdict.get(key) for key in
                ("passed", "engineering_passed", "errors", "warnings",
                 "physics_verdict", "wording", "allowed_action",
                 "failed_criteria", "inconclusive_criteria",
                 "missing_evidence_criteria",
                 "not_applicable_criteria")
            }
            manifest["checks"] = verdict.get("checks", [])
        except Exception as exc:  # noqa: BLE001
            manifest["verdict"] = {"error": "validate_junction_run 失败: %s" % exc}
    manifest["totals"] = totals
    return manifest


def _fmt(value: Any, spec: str = "%s") -> str:
    if value is None:
        return "—"
    if isinstance(value, bool):
        return "true" if value else "false"
    try:
        return spec % value
    except (TypeError, ValueError):
        return str(value)


def _cell(value: Any) -> str:
    """Markdown 表格单元：压平换行与竖线，避免破坏表结构。"""
    text = "—" if value is None else str(value)
    return text.replace("|", "/").replace("\n", " ").strip()


def render_report(manifest: Mapping[str, Any]) -> str:
    """manifest -> Markdown 报告（机械可生成部分，不写文件）。"""
    man = dict(manifest or {})
    rid = man.get("run_id", "?")
    state = man.get("state") or {}
    totals = man.get("totals") or {}
    lines: list[str] = []
    lines.append("# smj run 报告: %s" % rid)
    lines.append("")
    lines.append("> 由 `smj.report` 从 run 目录（state.json / 阶段子目录 / manifest.json）机械生成；")
    lines.append("> 数值结论仍须人工复核，报告不含图与独立科学裁决。")
    lines.append("")

    lines.append("## 0. 概览")
    lines.append("")
    lines.append("| 项 | 值 |")
    lines.append("|---|---|")
    lines.append("| run_id | `%s` |" % _cell(rid))
    lines.append("| 根目录 | `%s` |" % _cell(man.get("root")))
    lines.append("| 模式 | %s |" % _cell((man.get("inputs") or {}).get("mode")))
    prov = man.get("provenance") or {}
    if prov:
        protocol = str(prov.get("protocol_sha256") or "")[:12]
        if prov.get("reconstructed"):
            snap = "重建（reconstructed，非原始快照）"
        elif prov.get("path"):
            snap = "run_spec.json（sha256 %s…）" % str(prov.get("run_spec_sha256") or "")[:12]
        else:
            snap = "缺失"
        lines.append("| 参数快照 | %s；protocol %s… |" % (_cell(snap), _cell(protocol)))
    lines.append("| 阶段数 / 文件总量 | %s / %.3f GiB |"
                 % (totals.get("stages"), float(totals.get("gib") or 0.0)))
    statuses = state.get("stages") or {}
    summary = ", ".join("%s=%s" % (k, statuses[k]) for k in STAGE_ORDER if k in statuses)
    lines.append("| 阶段状态 | %s |" % _cell(summary))
    lines.append("| 首个未完成阶段 | `%s` |" % _cell(state.get("first_incomplete")))
    lines.append("| 生成时间 | %s |" % _cell(man.get("generated_at")))
    lines.append("")

    lines.append("## 0.1 判定")
    lines.append("")
    verdict = man.get("verdict")
    if isinstance(verdict, dict) and verdict.get("physics_verdict"):
        lines.append("| 项 | 值 |")
        lines.append("|---|---|")
        lines.append("| 工程判定 | %s |"
                     % ("通过" if verdict.get("passed") else
                        "未通过（errors=%s）" % _cell(verdict.get("errors"))))
        lines.append("| 物理判定 | %s（%s） |"
                     % (_cell(verdict.get("wording")),
                        _cell(verdict.get("physics_verdict"))))
        lines.append("| 允许动作 | %s |" % _cell(verdict.get("allowed_action")))
        lines.append("")
        failed = list(verdict.get("failed_criteria") or [])
        inconclusive = list(verdict.get("inconclusive_criteria") or [])
        missing = list(verdict.get("missing_evidence_criteria") or [])
        if failed:
            lines.append("- 否证 H0 的 error 级判据: %s" % ", ".join("`%s`" % x for x in failed))
        if missing:
            lines.append("- 证据缺失（无法判定, 非否证）: %s"
                         % ", ".join("`%s`" % x for x in missing))
        if inconclusive:
            lines.append("- 未过的 warn 级判据: %s" % ", ".join("`%s`" % x for x in inconclusive))
        if not failed and not inconclusive and not missing:
            lines.append("- 全部判据通过（仅代表判据层面未被否证）。")
    elif isinstance(verdict, dict) and verdict.get("error"):
        lines.append("> 判定不可用: %s" % _cell(verdict.get("error")))
    else:
        lines.append("> 未计算判定（`--no-checks`，或 manifest 由旧版本生成）。")
    lines.append("")

    lines.append("## 1. 阶段状态")
    lines.append("")
    lines.append("| 阶段 | 状态 | 目录 | 作业 | 退出码 | 耗时 (s) | 说明 |")
    lines.append("|---|---|---|---|---|---|---|")
    inventory = man.get("stages") or {}
    details = state.get("details") or {}
    failures = {f["stage"]: f for f in (man.get("failures") or [])}
    for stage in STAGE_ORDER:
        status = statuses.get(stage, "pending")
        has_dir = stage in inventory
        detail = details.get(stage) or {}
        note = detail.get("note") or ""
        if stage in failures:
            note = failures[stage].get("error") or failures[stage].get("note") or note
        lines.append("| `%s` | `%s` | %s | %s | %s | %s | %s |" % (
            stage, status, "有" if has_dir else "—",
            _cell(detail.get("job_id")), _cell(detail.get("exit_code")),
            _cell(detail.get("elapsed_s")), _cell(note)))
    lines.append("")

    lines.append("## 2. 输入摘要")
    lines.append("")
    lines.append("| 项 | 值 |")
    lines.append("|---|---|")
    inputs = man.get("inputs") or {}
    molecule = inputs.get("molecule") or {}
    lines.append("| 分子结构 | %s |" % _cell(molecule.get("path")))
    lines.append("| 分子原子数 | %s |" % _cell(molecule.get("atom_count")))
    electrode = inputs.get("electrode") or {}
    lines.append("| 建结阶段 | %s |" % _cell(electrode.get("stage")))
    lines.append("| 结原子数 | %s |" % _cell(electrode.get("atom_count")))
    lines.append("| 每侧电极 Au | %s |" % _cell(electrode.get("au_per_electrode")))
    lines.append("| 输运温度 (K) | %s |" % _cell(inputs.get("temperature_K")))
    te = inputs.get("te") or {}
    lines.append("| TE.dat | %s |" % _cell(te.get("source")))
    lines.append("| TE 数据点 | %s |" % _cell(te.get("points")))
    if te.get("efermi_Ha") is not None:
        lines.append("| efermi (Ha) | %s |" % _cell(te.get("efermi_Ha")))
    lines.append("")

    lines.append("## 3. 资源估算")
    lines.append("")
    lines.append("| 阶段 | 墙钟 | --time | 内存 | rank/核 | 核时 | 来源 / 预算 |")
    lines.append("|---|---|---|---|---|---|---|")
    resources = man.get("resources") or {}
    for stage in ("gaussian_opt", "aims_tip_relax", "au_relax", "aims_scf", "aitranss_te"):
        item = resources.get(stage)
        if not item:
            lines.append("| `%s` | — | — | — | — | — | 未估算 |" % stage)
            continue
        budget = item.get("budget") or {}
        source = item.get("source", "estimate")
        tag = "%s; %s" % (source, "预算内" if budget.get("ok") else budget.get("reason"))
        if item.get("capped"):
            tag += "; capped=%s" % ",".join(item["capped"])
        lines.append("| `%s` | %s | %s | %s | %s | %s | %s |" % (
            stage, _cell(item.get("walltime")), _cell(item.get("time_limit")),
            _cell(item.get("mem")), _cell(item.get("cores") or item.get("ntasks")),
            _cell(item.get("core_h")), _cell(tag)))
    lines.append("")

    lines.append("## 4. 输运结果")
    lines.append("")
    transport = man.get("transport")
    if transport:
        lines.append("| 量 | 值 |")
        lines.append("|---|---|")
        nearest = transport.get("fermi_nearest") or {}
        lines.append("| T(E_F) [G0] | %s |" % _fmt(nearest.get("T"), "%.6f"))
        lines.append("| G/G0 = L0 | %s |" % _fmt(transport.get("L0"), "%.9f"))
        lines.append("| log10(G/G0) | %s |" % _fmt(transport.get("log10_G_over_G0_like"), "%.6f"))
        lines.append("| S (µV/K) | %s |" % _fmt(transport.get("S_uV_per_K"), "%.4f"))
        lines.append("| k_e (pW/K) | %s |" % _fmt(transport.get("k_e_pW_per_K"), "%.4f"))
        lines.append("| ZT_electronic | %s |" % _fmt(transport.get("ZT_electronic"), "%.6f"))
        lines.append("| 数据点 | %s |" % _cell(transport.get("points")))
        lines.append("| 来源 | `%s` |" % _cell(transport.get("source")))
    else:
        lines.append("- 无 TE.dat / analysis.json，输运结果不可得。")
    lines.append("")

    lines.append("## 5. 失败原因")
    lines.append("")
    failures = man.get("failures") or []
    if failures:
        for item in failures:
            detail = item.get("error") or item.get("note") or "（无错误信息）"
            lines.append("- `%s` (job=%s, exit=%s): %s" % (
                item.get("stage"), item.get("job_id"), item.get("exit_code"),
                _cell(detail)))
    else:
        lines.append("- 无失败阶段。")
    lines.append("")

    lines.append("## 6. 产物清单")
    lines.append("")
    lines.append("| 阶段 | 文件数 | 体积 (MiB) | 关键产物 sha256(前12) |")
    lines.append("|---|---|---|---|")
    for stage, info in inventory.items():
        key = [f for f in info.get("files", [])
               if f.get("name") in ("TE.dat", "geometry.in", "mol.log", "aims.out",
                                    "analysis.json")]
        shas = ", ".join("%s:%s" % (f["name"], (f.get("sha256") or "未哈希")[:12])
                         for f in key) or "—"
        lines.append("| `%s` | %d | %.2f | %s |" % (
            stage, info.get("n_files", 0), info.get("bytes", 0) / (1024 * 1024), shas))
    lines.append("")
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# 解析（报告往返 / 既有 report.md 复用）
# --------------------------------------------------------------------------- #
def _split_sections(text: str) -> dict[str, list[str]]:
    sections: dict[str, list[str]] = {}
    current: str | None = None
    for line in text.splitlines():
        if line.startswith("## "):
            current = line[3:].strip()
            sections[current] = []
        elif current is not None:
            sections[current].append(line)
    return sections


def _table_rows(lines: Iterable[str]) -> list[list[str]]:
    rows = []
    for line in lines:
        stripped = line.strip()
        if not stripped.startswith("|"):
            continue
        cells = [c.strip() for c in stripped.strip("|").split("|")]
        if not cells or all(set(c) <= set("-: ") for c in cells if c):
            continue
        if cells[0] == "项" or cells[0] == "阶段":
            continue
        rows.append(cells)
    return rows


def _two_col(lines: Iterable[str]) -> dict[str, str]:
    return {row[0]: row[1] for row in _table_rows(lines) if len(row) >= 2}


def _number(text: str) -> Any:
    cleaned = text.replace("`", "").strip()
    try:
        return float(cleaned)
    except ValueError:
        return cleaned


def parse_report(text_or_path: str | Path) -> dict:
    """解析本模块（或既有 AgenticDFT run）的 ``report.md``。

    返回 ``{run_id, overview, stages, transport, failures, sections}``；
    表格值做尽力数值化（无法转换则保留字符串）。
    """
    if isinstance(text_or_path, Path) or (
            isinstance(text_or_path, str) and "\n" not in text_or_path
            and Path(text_or_path).is_file()):
        text = Path(text_or_path).read_text(encoding="utf-8", errors="ignore")
    else:
        text = str(text_or_path)
    match = re.search(r"^#\s+.*?报告[:：]\s*(\S+)\s*$", text, re.MULTILINE)
    if not match:
        match = re.search(r"^#\s+(\S+)\s*$", text, re.MULTILINE)
    run_id = match.group(1).strip("`") if match else None
    sections = _split_sections(text)

    overview = _two_col(sections.get("0. 概览", []))
    stages: dict[str, dict] = {}
    for row in _table_rows(sections.get("1. 阶段状态", [])):
        if len(row) < 2:
            continue
        stages[row[0].strip("`")] = {
            "status": row[1].strip("`"),
            "dir": row[2] if len(row) > 2 else None,
            "job": row[3] if len(row) > 3 else None,
            "exit_code": row[4] if len(row) > 4 else None,
            "elapsed_s": row[5] if len(row) > 5 else None,
            "note": row[6] if len(row) > 6 else None,
        }
    transport = {key: _number(value)
                 for key, value in _two_col(sections.get("4. 输运结果", [])).items()}
    failures = []
    for line in sections.get("5. 失败原因", []):
        stripped = line.strip()
        if stripped.startswith("- ") and stripped not in ("- 无失败阶段。",):
            failures.append(stripped[2:].strip())
    return {"run_id": run_id, "overview": overview, "stages": stages,
            "transport": transport, "failures": failures, "sections": sections}


# --------------------------------------------------------------------------- #
# 写文件入口
# --------------------------------------------------------------------------- #
def read_manifest(run_dir: str | Path) -> dict | None:
    """读既有 manifest.json（不可解析返回 None，不抛错）。"""
    path = Path(run_dir) / "manifest.json"
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def write_manifest(manifest: Mapping[str, Any], path: str | Path) -> Path:
    return _write_json(manifest, path)


def generate_report(target: Any, *, run_id: str | None = None,
                    spec: Mapping[str, Any] | None = None,
                    force: bool = False,
                    include_checks: bool = False,
                    system: str | None = None,
                    state: RunState | None = None,
                    stages: Iterable[str] | None = None) -> dict:
    """生成 ``<run>/manifest.json`` + ``<run>/report.md``。

    ``force=False`` 时逐文件跳过已存在者（外部 run 的文件不动），原因记入
    ``skipped``。返回 ``{ok, run_id, manifest_path, report_path, manifest,
    report, written, skipped}``。
    """
    run_dir, _config = _run_dir_of(target)
    layout = RunLayout(run_dir).ensure()
    manifest = build_manifest(target, run_id=run_id, spec=spec, state=state,
                              stages=stages, include_checks=include_checks,
                              system=system)
    skipped: list[dict] = []
    written: list[str] = []

    if layout.manifest_path.exists() and not force:
        skipped.append({"path": str(layout.manifest_path),
                        "reason": "manifest.json 已存在（force=True 才覆盖）"})
        effective = read_manifest(run_dir) or manifest
    else:
        write_manifest(manifest, layout.manifest_path)
        written.append(str(layout.manifest_path))
        effective = manifest

    if layout.report_path.exists() and not force:
        skipped.append({"path": str(layout.report_path),
                        "reason": "report.md 已存在（force=True 才覆盖）"})
        report_text = layout.report_path.read_text(encoding="utf-8", errors="ignore")
    else:
        report_text = render_report(effective)
        _write_text(layout.report_path, report_text)
        written.append(str(layout.report_path))

    return {
        "ok": True,
        "run_id": effective.get("run_id", manifest.get("run_id")),
        "manifest_path": str(layout.manifest_path),
        "report_path": str(layout.report_path),
        "manifest": effective,
        "report": report_text,
        "written": written,
        "skipped": skipped,
    }


#: CLI/兼容别名（remote-builder 的首选名是 generate_report）。
build_report = generate_report
write_report = generate_report
generate = generate_report
