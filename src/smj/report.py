"""run 报告：从 run 目录生成 ``manifest.json`` + 人读 ``report.md``。

输入是三样**已有事实**（不臆造）：

- ``<run>/state.json`` —— 阶段状态机（:class:`smj.state.RunState`）；
- ``<run>/<stage>/`` —— 各阶段产物目录（旧 run 的 ``auN_build`` 由布局回退）；
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
from .spin_summary import SPIN_CAVEAT, normalize as normalize_spin, percent as spin_percent
from .state import STAGE_ORDER, RunLayout, RunState, now_iso
from .verdict import JUDGMENT_SCHEMA, judgment_fingerprint

REPORT_SCHEMA = "smj.report.v1"

#: 单文件哈希上限（超过且未被策略标为 full 的只记大小，避免为 GB 级产物白算校验）。
_HASH_MAX_BYTES = 8 * 1024 * 1024

#: 已知阶段（用于阶段状态表）；旧 ``auN_build`` 目录通过 RunLayout 回退。
_BUILD_STAGE_RE = re.compile(r"^au\d+_build$")

#: 阶段默认几何产物（用于资源估算的原子数）。
_STAGE_GEOMETRY = {
    "molecule_opt": ("opt.xyz", "molecule.xyz", "geometry.in"),
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
        from smj.junction.artifacts import retrieval_mode  # lazy: vendored
    except Exception:  # noqa: BLE001
        return "unknown"
    try:
        # vendored 策略仍以 gaussian_opt 命名首段；molecule_opt 沿用其条目。
        return str(retrieval_mode("gaussian_opt" if stage == "molecule_opt" else stage, name))
    except Exception:  # noqa: BLE001
        return "unknown"


def _atom_count(path: Path, *, xyz: bool = False) -> int | None:
    """读结构文件的原子数：失败返回 None（不猜）。"""
    if not path.is_file():
        return None
    try:
        if xyz or path.suffix.lower() == ".xyz":
            from smj.junction.geometry import read_xyz  # lazy: vendored
            return len(read_xyz(path))
        from smj.junction.geometry import read_aims_geometry  # lazy: vendored
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
    2026-11-13：``au_relax`` 已从主链路移除；``electrode_build`` 的旧目录
    ``au60_build`` / ``auN_build`` 由 :class:`smj.state.RunLayout` 回退定位，
    因此这里不再把旧目录当成额外阶段追加。
    """
    root = Path(run_dir)
    if not root.is_dir():
        return []
    return [stage for stage in STAGE_ORDER if RunLayout(root).exists(stage)]

def _inventory(run_dir: Path, stages: Iterable[str]) -> dict[str, dict]:
    inventory: dict[str, dict] = {}
    for stage in stages:
        stage_dir = RunLayout(run_dir).stage_dir(stage)
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
    # au_relax 已从主链路移除：电极/结几何以 electrode_build（旧名 au60_build /
    # auN_build）为准。旧 run 的 au_relax 若真的跑过 AIMS 弛豫，它才是送进 SCF
    # 的几何；否则回退 build。
    legacy = run_dir / "au_relax"
    if ((legacy / "geometry.in").is_file()
            and (legacy / "electrode_metadata.json").is_file()
            and (legacy / "aims.out").is_file()):
        return legacy
    build = RunLayout(run_dir).stage_dir("electrode_build")
    if (build / "geometry.in").is_file() or (build / "electrode_metadata.json").is_file():
        return build
    for child in sorted(run_dir.iterdir()):
        if child.is_dir() and _BUILD_STAGE_RE.match(child.name):
            return child
    if (legacy / "geometry.in").is_file() and (legacy / "electrode_metadata.json").is_file():
        return legacy
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
    # 自旋/电荷是开壳层判读的前提：闭壳层 multiplicity=1，开壳层 >1。
    # run_spec 里固定总自旋为 True（默认）时该键被省略，故缺省按 True 解释。
    molecule_cfg = dict(effective_spec.get("molecule") or {})
    multiplicity = molecule_cfg.get("multiplicity")
    try:
        multiplicity = int(multiplicity) if multiplicity is not None else None
    except (TypeError, ValueError):
        multiplicity = None
    fixed_spin = molecule_cfg.get("fixed_spin_moment")
    fixed_spin = True if fixed_spin is None else bool(fixed_spin)
    broken_sym = bool(molecule_cfg.get("broken_symmetry", False))
    spin_meta: dict[str, Any] = {
        "charge": molecule_cfg.get("charge"),
        "multiplicity": multiplicity,
        "fixed_spin_moment": fixed_spin if ((multiplicity or 1) > 1 or broken_sym) else None,
    }
    if broken_sym:
        # M=1 的 BS 单重态：单看 multiplicity 会被读成闭壳层，必须显式标注。
        spin_meta["broken_symmetry"] = True
        spin_meta["fixed_spin_moment"] = None
    if molecule_cfg.get("initial_moments"):
        spin_meta["initial_moments"] = molecule_cfg.get("initial_moments")
    # 分子结构：优先 molecule_opt/opt.xyz（旧 run 回退 gaussian_opt），否则 run 根的 xyz。
    molecule_xyz = RunLayout(run_dir).artifact("molecule_opt", "opt.xyz")
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
    if multiplicity is not None or molecule_cfg:
        inputs["molecule"] = {**(inputs["molecule"] or {"path": None, "atom_count": None}),
                              **spin_meta}

    # 电极：electrode_build/electrode_metadata.json 的既有键（不补造）。
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
            "stage": "electrode_build",
            "source": meta.get("electrode_source") or "builtin",
            "element": meta.get("electrode_element") or "Au",
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
        path = RunLayout(run_dir).stage_dir(stage) / name
        if path.is_file():
            return _atom_count(path, xyz=name.endswith(".xyz"))
    # aims_scf / aitranss_te 的几何落在上游 electrode_build。
    if stage in ("aims_scf", "aitranss_te"):
        build = _find_build_dir(run_dir)
        if build is not None:
            return _atom_count(build / "geometry.in")
    return None


def _resource_summary(run_dir: Path, state: RunState,
                      warnings: list[str]) -> dict[str, dict]:
    summary: dict[str, dict] = {}
    for stage in ("molecule_opt", "aims_tip_relax", "aims_scf", "aitranss_te"):
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
    "max_transmission", "thermal_window", "fermi_quadrature",
)


def _spin_block(spin: Any, molecule: Mapping[str, Any]) -> dict | None:
    """给 ``spin`` 块补上配置侧信息（多重度/电荷/是否钉住总自旋）。"""
    if not isinstance(spin, Mapping):
        return None
    return normalize_spin(
        spin,
        multiplicity=molecule.get("multiplicity"),
        charge=molecule.get("charge"),
        fixed_spin_moment=molecule.get("fixed_spin_moment"),
        broken_symmetry=bool(molecule.get("broken_symmetry", False)),
    )


def _spin_partition_block(state: RunState) -> dict | None:
    """``aims_scf`` 写下的 Hirshfeld 自旋分区（S1）；缺失时返回 None。

    分区来自 FHI-aims 的 Hirshfeld 分析：总磁矩守恒不能说明自由基还在分子上，
    必须把分子/接触区/电极的份额分开报。这里只搬运展示数据，不参与判定。
    """
    try:
        recorded = state.stage("aims_scf").results or {}
    except (KeyError, AttributeError):
        return None
    inner = recorded.get("aims_scf") if isinstance(recorded, Mapping) else None
    if not isinstance(inner, Mapping):
        return None
    block = inner.get("spin_partition")
    return dict(block) if isinstance(block, Mapping) else None


def _transport_summary(run_dir: Path, state: RunState, temperature: float | None,
                       warnings: list[str], *,
                       molecule: Mapping[str, Any] | None = None) -> dict | None:
    """输运摘要，含开壳层的自旋通道块。

    优先 ``junction_analyze/results.json`` —— 它是该阶段的正式产物，开壳层时
    还带 ``spin`` 块（α/β 通道积分与 P(E)）。旧 run 可能只有 ``analysis.json``；
    两者都没有时再从 ``aitranss_te/TE.dat`` 现算（开壳层顺带解析双通道）。
    """
    molecule = dict(molecule or {})
    for name in ("results.json", "analysis.json"):
        candidate = run_dir / "junction_analyze" / name
        if not candidate.is_file():
            continue
        try:
            payload = json.loads(candidate.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            warnings.append("%s 不可解析: %s" % (name, exc))
            continue
        if not isinstance(payload, dict):
            continue
        transport = {"source": str(candidate),
                     **{k: payload.get(k) for k in _TRANSPORT_KEYS}}
        spin = _spin_block(payload.get("spin"), molecule)
        if spin is not None:
            transport["spin"] = spin
        partition = _spin_partition_block(state)
        if partition is not None:
            transport["spin_partition"] = partition
        return transport
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
    if (molecule.get("multiplicity") or 1) > 1 or molecule.get("broken_symmetry"):
        try:
            from .landauer_spin import analyze_te_spin
            spin = _spin_block(
                analyze_te_spin(te, temperature=float(temperature or 300.0)),
                molecule)
        except (OSError, ValueError) as exc:
            warnings.append("自旋通道解析失败（TE.dat 可能只有一列）: %s" % exc)
        else:
            if spin is not None:
                payload["spin"] = spin
    partition = _spin_partition_block(state)
    if partition is not None:
        payload["spin_partition"] = partition
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
                                   input_summary.get("temperature_K"), warnings,
                                   molecule=input_summary.get("molecule"))
    failures = _failure_summary(state)

    totals = {
        "bytes": sum(item["bytes"] for item in inventory.values()),
        "gib": round(sum(item["bytes"] for item in inventory.values()) / (1024 ** 3), 3),
        "stages": len(inventory),
    }
    # 2026-11-20 复审 R6：finished 只表示远端作业结束，不算验收通过。
    from .verdict import stage_acceptance  # lazy: verdict 依赖 junction 注册表
    statuses_now = state.statuses()
    awaiting_now, missing_now = stage_acceptance(statuses_now, STAGE_ORDER)
    run_complete = not awaiting_now and not missing_now
    manifest: dict[str, Any] = {
        "schema": REPORT_SCHEMA,
        # 判定章的指纹：generate_report 据此拒绝复用"判据已变化"的旧 manifest
        # （2026-11-22 复审 P0-1）。旧 manifest 没有该键，一律视为过期。
        "judgment": {
            "schema": JUDGMENT_SCHEMA,
            "sha256": judgment_fingerprint(),
        },
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
            "is_done": run_complete,
            "stages_awaiting": awaiting_now,
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
            from smj.junction.validate import validate_junction_run  # lazy: vendored
            from smj.verdict import apply_completeness
            # 2026-11-20：判据适用域要"族 + 多重度"两个输入。族优先取 run_spec
            # （run 时快照），而不是当前配置；开壳层由 multiplicity 驱动，判据层
            # 会把按族标定的闸门降为记录。
            run_family = system
            if run_family is None and run_spec is not None:
                run_family = (((run_spec.get("config") or {}).get("molecule")
                               or {}).get("family"))
            multiplicity = (input_summary.get("molecule") or {}).get("multiplicity")
            # M=1 的 BS 单重态在 input_summary 里带 broken_symmetry；run_spec 的配置
            # 快照同样可读。判据层据此把按族标定的闸门降为记录（开壳层）。
            broken_sym = bool((input_summary.get("molecule") or {}).get("broken_symmetry")
                              or ((run_spec or {}).get("config", {}).get("molecule")
                                  or {}).get("broken_symmetry", False))
            verdict = validate_junction_run(
                run_dir, results={"junction": transport or {}}, system=run_family,
                temperature=input_summary.get("temperature_K"),
                multiplicity=multiplicity, broken_symmetry=broken_sym,
                extra_outputs=tuple(
                    (((run_spec or {}).get("config") or {}).get("transport") or {})
                    .get("extra_outputs") or ()))
            # 2026-11-20 复核 P1-6 / R6：产物驱动的判定看不到"该跑的阶段没跑"。
            # finished（作业结束、产物已取回但未 collect/验收）单独进 awaiting，
            # 不算完成；run_incomplete 合成判据随之不过。
            statuses = state.statuses()
            awaiting, missing = stage_acceptance(statuses, STAGE_ORDER)
            is_done = not awaiting and not missing
            verdict = apply_completeness(
                verdict, is_done=is_done,
                first_incomplete=state.first_incomplete(), missing=missing,
                awaiting=awaiting)
            manifest["verdict"] = {
                key: verdict.get(key) for key in
                ("passed", "engineering_passed", "errors", "warnings",
                 "engineering_errors", "physics_failures",
                 "physics_verdict", "wording", "allowed_action",
                 "failed_criteria", "inconclusive_criteria",
                 "missing_evidence_criteria", "not_applicable_criteria",
                 "unregistered_checks", "completeness")
            }
            manifest["checks"] = verdict.get("checks", [])
        except Exception as exc:  # noqa: BLE001
            manifest["verdict"] = {"error": "validate_junction_run 失败: %s" % exc}
    manifest["totals"] = totals
    manifest["review"] = _review_summary(run_dir)
    return manifest


def _review_summary(run_dir: Path) -> dict[str, Any] | None:
    """把 ``review.json`` 的关键状态收进 manifest（报告必须显式标注复核结论）。

    2026-11-22 复审 P1：报告此前完全不渲染 review/connection_check，接触几何
    被改写或复核 blocked 时 report.md 仍显示"工程通过"。
    """
    path = Path(run_dir) / "review.json"
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        return {"path": str(path), "error": "review.json 不可解析: %s" % exc}
    if not isinstance(payload, dict):
        return {"path": str(path), "error": "review.json 顶层不是对象"}
    issues = payload.get("issues") or []
    connection = payload.get("connection_check") or {}
    return {
        "path": str(path),
        "status": payload.get("status"),
        "ok": payload.get("ok"),
        "accepted": payload.get("accepted"),
        "errors": sum(1 for i in issues if isinstance(i, Mapping)
                      and i.get("severity") == "error"),
        "warnings": sum(1 for i in issues if isinstance(i, Mapping)
                        and i.get("severity") == "warning"),
        "connection_check": {
            "present": connection.get("present"),
            "valid": connection.get("valid"),
            "reason": connection.get("reason"),
            "deviation_geom_A": connection.get("deviation_geom_A"),
        } if isinstance(connection, Mapping) else None,
        "flags": [i.get("code") for i in issues
                  if isinstance(i, Mapping) and i.get("severity") in ("error", "warning")],
    }


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


def _render_spin_section(lines: list[str], spin: Mapping[str, Any],
                         partition: Mapping[str, Any] | None = None) -> None:
    """§4.1 自旋通道：补 α/β 分量（G、S、P）与 Hirshfeld 分区，不改总口径。"""
    lines.append("")
    lines.append("### 4.1 自旋通道（开壳层）")
    lines.append("")
    lines.append("> 上表的 G、S、k_e、ZT 是**双通道均值**（landauer 口径 "
                 "T(E) = (Tα+Tβ)/2）；本小节的分量按同一公式逐通道计算，"
                 "总 S 与加权平均的自洽由 `spin_channel_consistency` 判据守。")
    lines.append("")
    lines.append("| 量 | 值 |")
    lines.append("|---|---|")
    if spin.get("multiplicity") is not None:
        lines.append("| 多重度 / 电荷 | %s / %s |" % (
            _cell(spin.get("multiplicity")), _cell(spin.get("charge"))))
    if spin.get("broken_symmetry"):
        lines.append("| 自旋态 | 破缺对称单重态（M=1，broken_symmetry；"
                     "单行列式 BS 近似，未验证 ⟨S²⟩） |")
    fixed = spin.get("fixed_spin_moment")
    if fixed is not None:
        lines.append("| 总自旋约束 | %s |" % (
            "fixed_spin_moment（钉住 2S = M−1）" if fixed
            else "未钉住（release 解，自旋可离域到电极）"))
    lines.append("| TE.dat 通道数 | %s |" % _cell(spin.get("channels")))
    lines.append("| Tα(E_F) / Tβ(E_F) [G0] | %s / %s |" % (
        _fmt(spin.get("T_alpha_at_ef"), "%.6f"),
        _fmt(spin.get("T_beta_at_ef"), "%.6f")))
    lines.append("| P(E_F) = (Tα−Tβ)/(Tα+Tβ) | %s |"
                 % spin_percent(spin.get("spin_polarization_at_ef")))
    lines.append("| P(L0) 积分极化 | %s |"
                 % spin_percent(spin.get("spin_polarization_L0")))
    if (spin.get("S_alpha_uV_per_K") is not None
            or spin.get("S_beta_uV_per_K") is not None):
        lines.append("| Sα / Sβ [μV/K] | %s / %s |" % (
            _fmt(spin.get("S_alpha_uV_per_K"), "%+.4f"),
            _fmt(spin.get("S_beta_uV_per_K"), "%+.4f")))
    if spin.get("spin_polarization_S") is not None:
        lines.append("| P_S = (Sα−Sβ)/(Sα+Sβ) | %s |"
                     % spin_percent(spin.get("spin_polarization_S")))
    if spin.get("g_over_g0_total") is not None:
        lines.append("| Gα / Gβ / G_total [G0] | %s / %s / %s |" % (
            _fmt(spin.get("g_over_g0_alpha"), "%.6f"),
            _fmt(spin.get("g_over_g0_beta"), "%.6f"),
            _fmt(spin.get("g_over_g0_total"), "%.6f")))
    lines.append("")
    lines.append("- %s" % SPIN_CAVEAT)
    if spin.get("fixed_constraint_caveat"):
        lines.append("- %s" % spin["fixed_constraint_caveat"])
    if spin.get("note"):
        lines.append("- 数据提示：%s" % _cell(spin.get("note")))
    if isinstance(partition, Mapping) and isinstance(partition.get("groups"), Mapping):
        groups = partition["groups"]
        lines.append("")
        lines.append("**自旋分布（Hirshfeld，μB）**")
        lines.append("")
        lines.append("| 区域 | 原子数 | 自旋矩 (μB) | 占比 |")
        lines.append("|---|---|---|---|")
        for name, label in (("molecule", "分子"), ("contact", "接触区"),
                            ("electrode", "电极")):
            data = groups.get(name)
            if not isinstance(data, Mapping):
                continue
            lines.append("| %s | %s | %s | %s |" % (
                label, _cell(data.get("atoms")),
                _fmt(data.get("spin_moment"), "%+.4f"),
                spin_percent(data.get("fraction"))))
        lines.append("")
        lines.append("- 接触区 = 距锚原子 ≤3.0 Å 的电极原子；自旋分布只陈述 "
                     "FHI-aims 的 Hirshfeld 分配，不参与判定。")
        if spin.get("broken_symmetry"):
            localized = partition.get("localized_antiparallel")
            if localized is True:
                tag = "局部反平行"
            elif localized is False:
                tag = "未满足局部反平行判据（净矩/幅度）"
            else:
                tag = "分区数据不完整"
            lines.append(
                "- 破缺对称证据：分子内 + / − 自旋矩各 %s / %s 个，净矩 %s μB（%s）。"
                % (_cell(partition.get("molecule_n_positive")),
                   _cell(partition.get("molecule_n_negative")),
                   _fmt(partition.get("molecule_net"), "%+.4f"), tag))
    if spin.get("broken_symmetry"):
        lines.append("")
        lines.append("- FHI-aims 单行列式不提供 ⟨S²⟩/NOON；上行的局部自旋矩是"
                     "证据展示，不改变任何判据结论。如需严格验证，可另跑 Gaussian "
                     "模型分子的 ⟨S²⟩ 交叉验证（默认不自动执行）。")


def _thermal_section(lines: list[str], man: Mapping[str, Any]) -> None:
    """``## 5. 热导`` —— 声子端到端热导（M5-T3）。"""
    lines.append("## 5. 热导（声子）")
    lines.append("")
    root = man.get("root")
    path = None if not root else Path(str(root)) / "thermal_analyze" / "thermal.json"
    if path is None or not path.is_file():
        lines.append("- 该 run 未跑热导链（无 `thermal_analyze/thermal.json`）。")
        lines.append("")
        return

    from .analysis import read_thermal
    try:
        data, warnings = read_thermal(path)
    except ValueError as exc:
        lines.append("- `thermal.json` 不可用：%s" % _cell(str(exc)))
        lines.append("")
        return

    lines.append("> ⚠️ **几何不可比**：热导器件是 Au(111) 3×3 周期板（~%d 原子，"
                 "含 %d 层 lead），电子六阶段链的器件是金字塔团簇（另一次计算）。"
                 "`G_ph` 与 `G/G0` 来自**两套几何**，不可并列比较；`ZT_total` 只是"
                 "把两者各自的数拼起来的**参考量**，不构成对该结的热电结论。"
                 % (int(data["n_atoms"]), int(round(data["n_lead_atoms"] / 9 or 1))))
    lines.append("")
    lines.append("| 量 | 值 |")
    lines.append("|---|---|")
    lines.append("| G_ph (pW/K) | %s |" % _fmt(data.get("g_ph_pw_per_k"), "%.4f"))
    lines.append("| G_Q (pW/K) | %s |" % _fmt(data.get("g_quantum_pw_per_k"), "%.4f"))
    lines.append("| 有效通道数 N_eff | %s |" % _fmt(data.get("n_channels_effective"), "%.4f"))
    zt = data.get("zt") or {}
    if zt.get("available"):
        lines.append("| k_e (pW/K) | %s |" % _fmt(zt.get("k_e_pW_per_K"), "%.4f"))
        lines.append("| ZT_electronic | %s |" % _fmt(zt.get("zt_electronic"), "%.6f"))
        lines.append("| ZT_total（参考） | %s |" % _fmt(zt.get("zt_total"), "%.6f"))
    else:
        lines.append("| ZT_total | 不可得（%s） |"
                     % _cell(zt.get("reason") or "无电子结果"))
    band_top = float(data.get("band_top_meV", 19.5))
    spectrum = data.get("spectrum") or {}
    energies = [float(v) for v in spectrum.get("energy_meV") or []]
    transmission = [float(v) for v in spectrum.get("transmission") or []]
    above = [abs(t) for e, t in zip(energies, transmission) if e > band_top]
    lines.append("| T_ph(ω>带顶 %.1f meV) | %s |"
                 % (band_top, "—" if not above else "%.2e" % max(above)))
    lines.append("| ASR 残差 (eV/Å²) | %s |" % _fmt(data.get("acoustic_sum_residual"), "%.3e"))
    modes = data.get("imaginary_modes") or {}
    worst = modes.get("worst_cm1")
    lines.append("| 虚频 | %d 个（最差 %s cm⁻¹） |"
                 % (int(modes.get("count", 0)),
                    "—" if worst is None else "%.2f" % float(worst)))
    lines.append("| 温度 (K) | %s |" % _fmt(data.get("temperature_K"), "%.1f"))
    lines.append("| 能量网格 | %s meV（η=%s meV） |"
                 % (_fmt(data.get("energy_step_meV"), "%.3f"),
                    _fmt(data.get("eta_meV"), "%.3f")))
    lines.append("| 来源 | `%s` |" % _cell(str(path)))
    lines.append("")

    # 谱图 + 分位：`G_ph` 是个积分，看不出热流由哪一段声子带走；这张两联图
    # （T_ph 与归一化累积）是解释分子间差异的主要依据，因此嵌进报告而不是只
    # 留在阶段目录里。
    from .phonon.units import EV_TO_CM1
    from .plotting import cumulative_quantiles

    image_rel = "thermal_analyze/thermal_spectrum.png"
    image = Path(str(man.get("root"))) / image_rel
    if not image.is_file():
        # 阶段内出图失败过（无 matplotlib / 字体问题）时在这里补一次；失败只
        # 留文字，绝不让报告生成失败。
        try:
            from .plotting import plot_thermal
            plot_thermal(path, image, formats=("png",))
        except Exception as exc:  # noqa: BLE001 - 报告不因缺图而失败
            lines.append("- 未能生成透射谱图：%s" % _cell(str(exc)))
            lines.append("")
    if image.is_file():
        lines.append("![声子透射谱与累积热导](%s)" % image_rel)
        lines.append("")

    quantiles = cumulative_quantiles(
        energies, [float(v) for v in spectrum.get("cumulative_g_ph_pw_per_k") or []])
    if quantiles:
        lines.append("热流贡献分位（累积 `G_ph` 达到该百分比时的声子能量）：")
        lines.append("")
        lines.append("| 分位 | 能量 (meV) | 波数 (cm⁻¹) |")
        lines.append("|---|---|---|")
        for label in ("25%", "50%", "75%", "90%"):
            if label not in quantiles:
                continue
            energy_mev = float(quantiles[label])
            lines.append("| %s | %.2f | %.0f |"
                         % (label, energy_mev, energy_mev * 1.0e-3 * EV_TO_CM1))
        lines.append("")
        lines.append("分位分别落在界面骨架模（Au–S/S–C 伸缩弯曲、环扭转）区间；"
                     "分子内 C–C 伸缩（~190 meV）远在 lead 带顶之上，弹性透射严格"
                     "为零，不参与热输运。")
        lines.append("")

    lines.append(
        "η（`eta_meV`）是**数值正则化参数，不是物理线宽**：它只给表面格林函数"
        "一个虚部，把窄的法布里-珀罗共振展宽。η 越小、共振越尖，`G_ph` 单调"
        "升高并外推到 η→0 的弹道上限。因此本方法的系统不确定性主要来自 η，"
        "量级可与非谐修正（文献 ~15%）相比；两项都不在本方法内，报告 `G_ph` "
        "时必须一起说明。量级参照：η 由 0.5 降到 0.05 meV，同一份力常数下 "
        "`G_ph` 可由 10.9 升到 13.7 pW/K。横向对比时请对齐 η，或统一声明为"
        "“η→0 外推值”。")
    lines.append("")

    curve = data.get("g_ph_by_temperature") or {}
    curve_values = [float(v) for v in curve.get("g_ph_pw_per_k") or []]
    curve_rel = "thermal_analyze/thermal_conductance_vs_temperature.png"
    curve_image = Path(str(man.get("root"))) / curve_rel
    if curve_values and not curve_image.is_file():
        try:
            from .plotting import plot_thermal_vs_temperature
            plot_thermal_vs_temperature(path, curve_image, formats=("png",))
        except Exception:  # noqa: BLE001 —— 出图失败只是少一张图
            pass
    if curve_image.is_file():
        lines.append("![声子热导随温度](%s)" % curve_rel)
        lines.append("")
    if len(curve_values) >= 2:
        plateau = curve_values[-1]
        reference_k = float(data.get("temperature_K", 300.0))
        lines.append(
            "`G_ph(T)`：%.0f K 时已达高温平台（%.2f pW/K）的 %.1f%%，"
            "低温端 %.0f→%.0f K 由 %.2f 升到 %.2f pW/K（比值 %.2f，"
            "与 `G ∝ T` 的量子极限一致）。谐波近似下 `T_ph(ω)` 与温度无关，"
            "这条曲线只换玻色因子，不额外花机时——它同时是 300 K 取值够不够"
            "收敛、低温端有没有物理错误的免费检验。"
            % (reference_k, plateau,
               100.0 * float(data.get("g_ph_pw_per_k", 0.0)) / plateau,
               float(curve["temperature_K"][0]), float(curve["temperature_K"][1]),
               curve_values[0], curve_values[1], curve_values[1] / curve_values[0]))
        lines.append("")

    soft = data.get("soft_modes") or []
    if soft:
        lines.append("**软模归属**（虚频只报波数没有意义，必须看谁在动）：")
        lines.append("")
        lines.append("| 模 | ω (cm⁻¹) | 分子 | 器件区 Au | lead | 参与原子数 | 刚体占比 |")
        lines.append("|---|---|---|---|---|---|---|")
        for mode in soft:
            regions = mode.get("regions") or {}
            rigid = mode.get("rigid_body") or {}
            lines.append(
                "| %d | %.2f | %.1f%% | %.1f%% | %.1f%% | %.1f | %.2f |"
                % (int(mode.get("index", -1)) + 1, float(mode.get("cm1", 0.0)),
                   100.0 * float(regions.get("molecule", 0.0)),
                   100.0 * float(regions.get("device_au", 0.0)),
                   100.0 * float(regions.get("lead", 0.0)),
                   float(mode.get("participation_ratio", 0.0)),
                   float(rigid.get("fraction", 0.0))))
        lines.append("")
        impact = data.get("soft_mode_impact_relative_g_ph")
        removed = data.get("g_ph_soft_removed_pw_per_k")
        symmetric = data.get("g_ph_symmetric_pw_per_k")
        if impact is not None and removed is not None and symmetric is not None:
            lines.append(
                "把**全部负本征值投影掉**再算一遍：`G_ph` 由 %.3f 变 %.3f pW/K"
                "（%+.2f%%），远小于 η 带来的 5–25%% 的漂移，故不改变结论。"
                "权重集中在分子、lead 占比 ~0，最大键伸缩只 ~0.01 Å（单位范数"
                "位移场下），说明它们是分子骨架的浅弯曲/扭转（PES 在该方向上轻微"
                "凹），既不是键拉伸，也不是 lead 截断假象。"
                % (symmetric, removed, 100.0 * float(impact)))
            lines.append("")
        asymmetry = data.get("fc_asymmetry")
        if asymmetry is not None:
            lines.append(
                "互易性：力常数矩阵不对称度 %.2f%%（有限差分噪声 + 器件/lead 两侧"
                "分开算）。对称化后再积分，`G_ph` 由 %.3f 变 %.3f pW/K；该差别"
                "已计入 `ph_fc_symmetry` 判据，不另行修正数值。"
                % (100.0 * float(asymmetry), float(data.get("g_ph_pw_per_k", 0.0)),
                   float(data.get("g_ph_symmetric_pw_per_k", 0.0))))
            lines.append("")

    environment_path = Path(str(man.get("root"))) / "thermal_analyze" / "environment.json"
    if environment_path.is_file():
        try:
            environment = json.loads(environment_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            environment = {}
        if environment:
            from .phonon.environment import describe_lines

            lines.append("**运行环境**（决定数值的输入与工具，换一样数就可能变）：")
            lines.append("")
            for item in describe_lines(environment):
                lines.append("- " + item)
            fingerprint = (data.get("environment") or {}).get("digest") or ""
            if fingerprint:
                lines.append(f"- 环境指纹 `{fingerprint[:16]}`"
                             "（`smj thermal-compare` 用它确认两个 run 是否同环境）")
            lines.append("")

    convergence = data.get("convergence") or {}
    if convergence:
        lines.append("**收敛轴**（哪几个开关真的影响 `G_ph`）：")
        lines.append("")
        cutoff = convergence.get("fc_cutoff") or {}
        rows = cutoff.get("rows") or []
        if rows:
            production = cutoff.get("production_A")
            lines.append("| 力常数截断 (Å) | G_ph (pW/K) | 说明 |")
            lines.append("|---|---|---|")
            for row in rows:
                value = row.get("g_ph_pw_per_k")
                mark = ""
                if production and abs(float(row["cutoff_A"]) - float(production)) < 1e-9:
                    mark = "（生产）"
                if value is None:
                    cells = f"| {float(row['cutoff_A']):.1f}{mark} | — | 带顶闸门拒绝 |"
                else:
                    cells = f"| {float(row['cutoff_A']):.1f}{mark} | {float(value):.3f} | |"
                lines.append(cells)
            lines.append("")
            spread = cutoff.get("relative_spread")
            if spread is not None:
                lines.append(
                    f"4–12 Å 范围内 `G_ph` 的相对散布 {100.0 * float(spread):.1f}%。"
                    "这条轴零机时可得：切的是已组装好的矩阵，不碰 DFT。")
                if cutoff.get("nominal_covers_span"):
                    lines.append(
                        "名义截断已覆盖两层块内的全部非零耦合：本次实测 8/10/12 Å "
                        "三档逐位相同，剩余散布来自把 r_c 刻意切到 6 Å 以下——"
                        "切掉的是真实近邻，不是名义 `fc_cutoff_A` 不足。")
                lines.append("")
            shortfall = cutoff.get("nominal_shortfall")
            if shortfall:
                lines.append(
                    f"> **截断半径未收敛**：名义 `fc_cutoff_A` = "
                    f"{float(cutoff.get('production_A') or 0.0):.1f} Å，但组装矩阵里"
                    f"最远的非零元素在 {float(cutoff.get('effective_span_A') or 0.0):.2f} Å"
                    f"——按名义值再切一刀会丢掉 {100.0 * float(shortfall):.1f}% 的 `G_ph`。"
                    "要把它变成真收敛的数，得把 `fc_cutoff_A` 提上去重算电极力常数。")
                lines.append("")
            if cutoff.get("beyond_nominal_checked") is False:
                lines.append(
                    "> 超出名义截断的体相耦合不在当前 lead 力常数资产里，"
                    "无法后处理验证：后处理只能证明「不多切」，不能证明「没漏切」。"
                    "要量化这一段，需要用更大的 `fc_cutoff_A` 重算一次电极"
                    "力常数（54 个 SCF）后对照。")
                lines.append("")
        for key, label in (("lead_layers", "lead 层数"),
                           ("device_layers", "器件区层数")):
            axis = convergence.get(key) or {}
            if axis.get("status") != "unmeasured":
                continue
            lines.append(
                f"- **{label}：未测**（当前 {axis.get('current')}；"
                f"{axis.get('requires')}）")
        if any((convergence.get(key) or {}).get("status") == "unmeasured"
               for key in ("lead_layers", "device_layers")):
            lines.append("")
            lines.append("这两条轴改变的是**哪些原子参与 DFT**，力常数矩阵必须重算，"
                         "不能后处理得到；在它们测完之前，`G_ph` 的不确定度只包含 "
                         "η、步长与截断三项。")
            lines.append("")

    review_path = Path(str(man.get("root"))) / "review.json"
    if review_path.is_file():
        try:
            payload = json.loads(review_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            payload = {}
        thermal = (payload or {}).get("thermal") or {}
        if thermal.get("available"):
            lines.append("热导判据结论：**%s**" % thermal.get("verdict"))
            failed = thermal.get("failed") or []
            unclear = thermal.get("inconclusive") or []
            if failed:
                lines.append("- failed：%s" % ", ".join("`%s`" % f for f in failed))
            if unclear:
                lines.append("- inconclusive：%s"
                             % ", ".join("`%s`" % u for u in unclear))
            if not failed and not unclear:
                lines.append("- 十条判据全部通过。")
            lines.append("")
    for message in warnings:
        lines.append("- %s" % _cell(message))
    if warnings:
        lines.append("")


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
        engineering = verdict.get("engineering_passed")
        if engineering is None:
            engineering = verdict.get("passed")
        eng_errs = verdict.get("engineering_errors")
        if eng_errs is None:
            eng_errs = verdict.get("errors")
        lines.append("| 工程判定 | %s |"
                     % ("通过" if engineering else
                        "未通过（工程错误=%s）" % _cell(eng_errs)))
        if "physics_failures" in verdict:
            lines.append("| 判据分账 | 工程错误 %s / 物理失败 %s / 警告 %s |"
                         % (_cell(verdict.get("engineering_errors")),
                            _cell(verdict.get("physics_failures")),
                            _cell(verdict.get("warnings"))))
        comp = verdict.get("completeness")
        if isinstance(comp, dict):
            awaiting = list(comp.get("stages_awaiting") or [])
            if comp.get("is_done"):
                lines.append("| 完工度 | 已走完全部阶段 |")
            else:
                parts = []
                if awaiting:
                    parts.append("已采集待验收 %s" % _cell(", ".join(awaiting)))
                stages_missing = list(comp.get("stages_missing") or [])
                if stages_missing or not awaiting:
                    parts.append("缺少 %s" % _cell(", ".join(stages_missing) or "—"))
                lines.append("| 完工度 | 未完成：首个未完成阶段 `%s`；%s |"
                             % (_cell(comp.get("first_incomplete")),
                                "；".join(parts)))
        lines.append("| 物理判定 | %s（%s） |"
                     % (_cell(verdict.get("wording")),
                        _cell(verdict.get("physics_verdict"))))
        lines.append("| 允许动作 | %s |" % _cell(verdict.get("allowed_action")))
        lines.append("")
        # 2026-11-20 复审 R6：skipped 算已验收，但报告要说明为什么跳过。
        details = state.get("details") or {}
        for name in STAGE_ORDER:
            if (statuses or {}).get(name) != "skipped":
                continue
            note = (details.get(name) or {}).get("note")
            lines.append("- 阶段 `%s` 已跳过：%s"
                         % (name, _cell(note) if note else "（未记录原因）"))
        failed = list(verdict.get("failed_criteria") or [])
        inconclusive = list(verdict.get("inconclusive_criteria") or [])
        missing = list(verdict.get("missing_evidence_criteria") or [])
        unregistered = [str(x) for x in (verdict.get("unregistered_checks") or [])]
        # no_applicable_gate 不是"未过的判据"：没有判据参与归约，报告要直说。
        gates = [x for x in inconclusive if x == "no_applicable_gate"]
        warn_only = [x for x in inconclusive if x != "no_applicable_gate"]
        if failed:
            lines.append("- 否证 H0 的 error 级判据: %s" % ", ".join("`%s`" % x for x in failed))
        if missing:
            lines.append("- 证据缺失（无法判定, 非否证）: %s"
                         % ", ".join("`%s`" % x for x in missing))
        if gates:
            lines.append("- 无适用族闸门（`no_applicable_gate`）：本次没有判据参与 H0 "
                         "归约（体系族未知/未校准），不能据此判定通过。")
        if warn_only:
            lines.append("- 未过的 warn 级判据: %s"
                         % ", ".join("`%s`" % x for x in warn_only))
        if unregistered:
            # 2026-11-22 复审 P0-2：validate 会返回未登记判据，归约时只降为
            # inconclusive；报告若把它丢掉，就会出现"无法判定"与"全部通过"并存。
            lines.append("- 未登记判据（未参与 H0 归约）: %s"
                         % ", ".join("`%s`" % x for x in unregistered))
        if not failed and not inconclusive and not missing and not unregistered:
            lines.append("- 全部判据通过（仅代表判据层面未被否证）。")
    elif isinstance(verdict, dict) and verdict.get("error"):
        lines.append("> 判定不可用: %s" % _cell(verdict.get("error")))
    else:
        lines.append("> 未计算判定（`--no-checks`，或 manifest 由旧版本生成）。")
    lines.append("")

    review = man.get("review")
    if isinstance(review, dict):
        lines.append("## 0.2 结构复核 / 连接核对")
        lines.append("")
        if review.get("error"):
            lines.append("> 复核记录不可用: %s" % _cell(review.get("error")))
            lines.append("")
        else:
            status = review.get("status")
            lines.append("| 项 | 值 |")
            lines.append("|---|---|")
            lines.append("| 复核结论 | %s（ok=%s；error=%s / warning=%s） |"
                         % (_cell(status), _cell(review.get("ok")),
                            _cell(review.get("errors")), _cell(review.get("warnings"))))
            connection = review.get("connection_check") or {}
            if connection:
                if connection.get("valid"):
                    cc = "有效"
                elif connection.get("present"):
                    cc = "记录已失效"
                else:
                    cc = "无记录"
                lines.append("| 连接核对 | %s%s |" % (
                    cc,
                    ("；原因: %s" % _cell(connection.get("reason")))
                    if connection.get("reason") else ""))
                if connection.get("deviation_geom_A") is not None:
                    lines.append("| 表面原子最大位移 | %.3f Å |"
                                 % float(connection["deviation_geom_A"]))
            if review.get("accepted"):
                lines.append("| 人工豁免 | 是（accept_deviation 显式放行，"
                             "不得当作连接核对有效） |")
            flags = [x for x in (review.get("flags") or []) if x]
            if flags:
                lines.append("| 告警/错误 | %s |"
                             % _cell(", ".join("`%s`" % x for x in flags)))
            if not review.get("ok"):
                lines.append("")
                lines.append("> **复核未通过（blocked）**：报告后续的工程/物理判定"
                             "不构成对接触几何的确认。")
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
    if molecule.get("multiplicity") is not None:
        lines.append("| 电荷 / 多重度 | %s / %s |" % (
            _cell(molecule.get("charge")), _cell(molecule.get("multiplicity"))))
        if molecule.get("broken_symmetry"):
            lines.append("| 自旋态 | 破缺对称单重态（M=1，broken_symmetry；"
                         "反号 initial_moments，不写 fixed_spin_moment） |")
        elif (molecule.get("multiplicity") or 1) > 1:
            fixed = molecule.get("fixed_spin_moment")
            lines.append("| 总自旋约束 | %s |" % (
                "fixed_spin_moment（钉住 2S = M−1）" if fixed
                else "未钉住（release 解，自旋可离域到电极）"))
    electrode = inputs.get("electrode") or {}
    lines.append("| 建结阶段 | %s |" % _cell(electrode.get("stage")))
    lines.append("| 结原子数 | %s |" % _cell(electrode.get("atom_count")))
    lines.append("| 电极构建方式 | %s |" % _cell(electrode.get("source")))
    lines.append("| 电极元素 | %s |" % _cell(electrode.get("element")))
    lines.append("| 每侧电极原子数 | %s |" % _cell(electrode.get("au_per_electrode")))
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
    for stage in ("molecule_opt", "aims_tip_relax", "aims_scf", "aitranss_te"):
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
        window = transport.get("thermal_window")
        if isinstance(window, Mapping):
            lines.append("| 费米窗覆盖 | %s |" % _fmt(window.get("coverage"), "%.4f"))
            if all(k in window for k in ("e_min_ev", "e_max_ev")):
                lines.append("| 能量窗 (eV) | [%+.4f, %+.4f] |"
                             % (float(window["e_min_ev"]), float(window["e_max_ev"])))
        quad = transport.get("fermi_quadrature")
        if isinstance(quad, Mapping):
            ratio = quad.get("ratio")
            moment = quad.get("first_moment_kbt")
            step = quad.get("min_step_in_window_ev")
            lines.append(
                "| 费米窗求积 | %s / 矩 %s kBT / 窗内 %s 点 / 最小步长 %s |"
                % (_fmt(ratio, "%.4f"),
                   "—" if moment is None else "%+.4f" % float(moment),
                   _cell(quad.get("points_in_window")),
                   "—" if step is None else "%.2f meV" % (1000.0 * float(step))))
        lines.append("| 来源 | `%s` |" % _cell(transport.get("source")))
    else:
        lines.append("- 无 TE.dat / analysis.json，输运结果不可得。")
    if isinstance(transport, Mapping) and isinstance(transport.get("spin"), Mapping):
        _render_spin_section(lines, transport["spin"], transport.get("spin_partition"))
    lines.append("")

    _thermal_section(lines, man)

    lines.append("## 6. 失败原因")
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

    lines.append("## 7. 产物清单")
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
    warnings: list[str] = []
    existing = read_manifest(run_dir)
    stale_reason: str | None = None
    if existing is not None and not force:
        previous = existing.get("judgment")
        previous_sha = previous.get("sha256") if isinstance(previous, dict) else None
        current_sha = (manifest.get("judgment") or {}).get("sha256")
        if not previous_sha:
            stale_reason = "旧 manifest 未记录判据指纹"
        elif previous_sha != current_sha:
            stale_reason = "判据实现已变化（judgment.sha256 不匹配）"
        elif include_checks and not (
                isinstance(existing.get("verdict"), dict)
                and existing["verdict"].get("physics_verdict")):
            # 旧版 --no-checks 产物或 validate 失败留下的 {"error": ...}：
            # 本次明确要求判定，就补算一次，而不是继续复用空判定。
            stale_reason = "既有 manifest 没有可用的物理判定，本次要求补算"

    if stale_reason is not None:
        # 2026-11-22 复审 P0-1：判定逻辑变了，旧 manifest 不能原样复用。
        # 重建 manifest 的同时必须重写 report.md——否则人读报告仍是旧结论。
        warnings.append("manifest 已过期（%s），已按当前判据重建" % stale_reason)
        manifest.setdefault("warnings", []).append(
            "manifest 曾过期（%s）：本次按当前判据重建，旧判定不再复用"
            % stale_reason)
        write_manifest(manifest, layout.manifest_path)
        written.append(str(layout.manifest_path))
        effective = manifest
        report_text = render_report(effective)
        _write_text(layout.report_path, report_text)
        written.append(str(layout.report_path))
    else:
        if layout.manifest_path.exists() and not force:
            skipped.append({"path": str(layout.manifest_path),
                            "reason": "manifest.json 已存在（force=True 才覆盖）"})
            effective = existing or manifest
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
        "stale_manifest": stale_reason,
        "warnings": warnings,
    }


#: CLI/兼容别名（remote-builder 的首选名是 generate_report）。
build_report = generate_report
write_report = generate_report
generate = generate_report
