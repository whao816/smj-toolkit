"""结构 / 结果复核门（review）：连接核对、可接受偏差与 review.json。

复核分两层，**不混在一起**：

1. **结构几何** —— 有 ``electrode_build/geometry.in`` 时调 vendored
   ``validate_junction_geometry``（表面索引/输运轴/接触重叠/Au-S 锚定），
   error 级问题直接阻断；
2. **连接核对（connection check；2026-11-15 前叫 A2 人工确认）** —— 读
   ``electrode_build/connection_check.json``（旧 run 读 ``a2_confirmation.json``，
   vendored ``fingerprint.record_confirmation`` 的格式）与当前几何指纹：

   - ``ReviewConfig.connection_check=True`` 且**无核对记录** -> ``blocked``；
   - 确认有效（指纹一致）-> 通过；
   - 指纹失效（几何已变）时比较**表面索引偏差**
     ``1 - 交集 / 较大集合``（0..1，未匹配表面原子占比）：
     ``<= connection_check_deviation`` 视为可接受偏差；
     否则若显式给 ``accept_deviation`` 且 ``<= accept_deviation``，记为
     ``warning`` 并 ``accepted=True``（报告须标注）；再超 -> ``blocked``。
   - ``connection_check=False`` 时连接核对只是加分证据：缺失/失效记 info/warning，
     不阻断（但结构几何 error 仍然阻断）。

``ReviewConfig.auto_connection_check``（默认 True）让 ``smj review`` 在复核前
自动登记核对记录（审计通过才签；出图不可用时退化为只记指纹），
``runner`` 也在 electrode_build 采集后自动登记，所以默认路径不需要手工签署。

``review.json`` 是派生结果，默认覆盖（``force=False`` 才跳过）；与
``report.py`` 不覆盖外部 run 文件的策略不同，原因写在函数 docstring。

vendored 模块全部 lazy import。
"""

from __future__ import annotations

import json
import os
from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

from .config import ReviewConfig
from .provenance import (
    PROTOCOL_SCHEMA,
    compare_protocol_configs,
    load_run_spec,
    protocol_components,
    protocol_sha256,
)
from .errors import SmjError
from .spin_summary import findings as spin_findings, normalize as normalize_spin
from .state import RunLayout, normalize_stage, now_iso

REVIEW_SCHEMA = "smj.review.v1"

#: 连接核对相关文件（与 vendored fingerprint/a2_viz 同名；2026-11-15 更名）。
CONFIRM_FILE = "connection_check.json"
AUDIT_FILE = "connection_check_audit.json"
#: 更名前的兼容名（旧 run 目录只读）：先「表面核对」时代，再更早的 A2 时代。
LEGACY_CONFIRM_FILES = ("surface_check.json", "a2_confirmation.json")
LEGACY_AUDIT_FILES = ("surface_check_audit.json", "a2_surface_audit.json")
LEGACY_CONFIRM_FILE = LEGACY_CONFIRM_FILES[-1]
LEGACY_AUDIT_FILE = LEGACY_AUDIT_FILES[-1]
METADATA_FILE = "electrode_metadata.json"
GEOMETRY_FILE = "geometry.in"

SIDES = ("left", "right")


def _first_existing(directory: Path, names: tuple[str, ...]) -> Path | None:
    for name in names:
        candidate = directory / name
        if candidate.is_file():
            return candidate
    return None


# --------------------------------------------------------------------------- #
# 结果对象
# --------------------------------------------------------------------------- #
@dataclass
class ReviewResult(Mapping):
    status: str = "passed"
    issues: list[dict] = field(default_factory=list)
    connection_check: dict | None = None
    structure: dict | None = None
    spin: dict | None = None
    run_id: str | None = None
    accepted: bool = False
    review_path: str | None = None
    written: bool = False

    @property
    def ok(self) -> bool:
        """``blocked`` 才是不通过；warning/info 不阻断。"""
        return self.status != "blocked"

    @property
    def errors(self) -> int:
        return sum(1 for i in self.issues if i["severity"] == "error")

    @property
    def warnings(self) -> int:
        return sum(1 for i in self.issues if i["severity"] == "warning")

    def add(self, severity: str, code: str, message: str, **extra: Any) -> dict:
        if severity not in ("error", "warning", "info"):
            raise ValueError("非法 severity %r" % (severity,))
        issue = {"severity": severity, "code": code, "message": message}
        issue.update({k: v for k, v in extra.items() if v is not None})
        self.issues.append(issue)
        return issue

    def codes(self) -> list[str]:
        return [i.get("code", "") for i in self.issues]

    def finalize(self) -> "ReviewResult":
        """汇总状态：error 必阻断；结构几何的**方向性告警**只记录、不降级。

        连接核对门的 warning（缺记录但未要求、显式接受偏差、审计异常）才把整体
        状态降为 ``warning``；``geometry_missing`` 这类"复核无法执行"也算 warning。
        开壳层的 ``open_shell_uncalibrated``/``spin_channels_missing`` 同样降为
        ``warning``（判据族对该体系未校准，结论可信度必须显示出来），但不阻断。
        """
        if any(i["severity"] == "error" for i in self.issues):
            self.status = "blocked"
            return self
        gate_warning = any(
            i["severity"] == "warning" and (
                str(i.get("code", "")).startswith("connection_check_")
                or str(i.get("code", "")).startswith("a2_")
                or str(i.get("code", "")).startswith("provenance_")
                or i.get("code") in ("geometry_missing", "geometry_metadata_missing",
                                     "geometry_metadata_unreadable", "geometry_check_failed",
                                     "open_shell_uncalibrated", "spin_channels_missing"))
            for i in self.issues)
        self.status = "warning" if gate_warning else "passed"
        return self

    def to_dict(self, config: ReviewConfig | None = None) -> dict:
        if config is None:
            config = getattr(self, "config_used", None)
        return {
            "schema": REVIEW_SCHEMA,
            "generated_at": now_iso(),
            "run_id": self.run_id,
            "status": self.status,
            "ok": self.ok,
            "accepted": self.accepted,
            "config": {
                "connection_check": getattr(config, "connection_check", None),
                "connection_check_deviation": getattr(config, "connection_check_deviation", None),
                "auto_connection_check": getattr(config, "auto_connection_check", None),
                "accept_deviation": getattr(config, "accept_deviation", None),
            },
            "issues": list(self.issues),
            "errors": self.errors,
            "warnings": self.warnings,
            "connection_check": self.connection_check,
            "structure": self.structure,
            "spin": self.spin,
            "review_path": self.review_path,
            "written": self.written,
        }

    # -- dict 式访问 -------------------------------------------------------- #
    def __getitem__(self, key: str) -> Any:
        return self.to_dict()[key]

    def get(self, key: str, default: Any = None) -> Any:
        return self.to_dict().get(key, default)

    def __iter__(self):
        return iter(self.to_dict())

    def __len__(self) -> int:
        return len(self.to_dict())


# --------------------------------------------------------------------------- #
# 工具
# --------------------------------------------------------------------------- #
def _write_json(payload: Mapping[str, Any], path: str | Path) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_name(target.name + ".tmp")
    tmp.write_text(json.dumps(dict(payload), indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8", newline="\n")
    os.replace(tmp, target)
    return target


def _run_dir_of(target: Any) -> tuple[Path, Any]:
    if hasattr(target, "run_dir") and not isinstance(target, (str, Path)):
        return Path(target.run_dir), target
    return Path(target), None


def _has_confirmation(path: Path) -> bool:
    """该目录是否有连接核对记录（新名或更名前的旧名）。"""
    return _first_existing(path, (CONFIRM_FILE,) + LEGACY_CONFIRM_FILES) is not None


def find_build_dir(run_dir: str | Path) -> Path | None:
    """定位连接核对的对象：``electrode_build``（旧名 ``au60_build`` / ``auN_build``）。

    2026-11-13 起 ``au_relax`` 已从主链路移除，核对的对象就是贴装后的建结几何；
    旧 run（只有 ``au_relax`` 而没有 build 目录）作为兼容回退。
    """
    root = Path(run_dir)
    if not root.is_dir():
        return None
    preferred = RunLayout(root).stage_dir("electrode_build")
    legacy = root / "au_relax"
    if preferred.is_dir():
        # 旧 run 的核对记录签在 au_relax 几何上；若新口径的 electrode_build 还没有记录，
        # 仍以旧对象为准，避免显式兼容旧 run 时误报 connection_check_missing。
        if _has_confirmation(legacy) and not _has_confirmation(preferred):
            return legacy
        return preferred
    for child in sorted(root.iterdir()):
        if child.is_dir() and child.name.startswith("au") and child.name.endswith("_build"):
            return child
    if (legacy / GEOMETRY_FILE).is_file() and (legacy / METADATA_FILE).is_file():
        return legacy
    return None


def _indices(data: Mapping[str, Any] | None) -> dict[str, list[int]] | None:
    """从 confirmation/metadata 里取左右表面索引（兼容 indices 嵌套）。"""
    if not isinstance(data, Mapping):
        return None
    payload = data.get("indices") if isinstance(data.get("indices"), Mapping) else data
    surface = payload.get("surface_indices") if isinstance(
        payload.get("surface_indices"), Mapping) else payload
    out: dict[str, list[int]] = {}
    for side in SIDES:
        values = surface.get("%s_surface_indices" % side)
        if values is None:
            values = data.get("%s_surface_indices" % side)
        if values:
            try:
                out[side] = [int(v) for v in values]
            except (TypeError, ValueError):
                return None
    return out if len(out) == 2 else None


def surface_index_deviation(recorded: Mapping[str, Any] | None,
                            current: Mapping[str, Any] | None) -> float | None:
    """表面索引偏差：``1 - 交集 / 较大集合``（0=完全一致，1=完全换了一套原子）。

    两侧分别算后按基数加权汇总；任一侧缺索引 -> ``None``（无法量化，不臆断）。
    """
    rec = _indices(recorded)
    cur = _indices(current)
    if not rec or not cur:
        return None
    overlap = 0
    denom = 0
    for side in SIDES:
        a, b = set(rec[side]), set(cur[side])
        overlap += len(a & b)
        denom += max(len(a), len(b))
    return (1.0 - float(overlap) / denom) if denom else None


def _load_json(path: Path) -> tuple[dict | None, str | None]:
    if not path.is_file():
        return None, None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        return None, str(exc)
    return (payload if isinstance(payload, dict) else None), None


# --------------------------------------------------------------------------- #
# 结构几何
# --------------------------------------------------------------------------- #
def _structure_review(result: ReviewResult, build_dir: Path | None,
                      geometry: str | Path | None = None,
                      metadata: Mapping[str, Any] | None = None) -> None:
    if build_dir is None:
        result.add("warning", "geometry_missing",
                   "无 electrode_build 阶段目录，结构几何复核跳过（阶段未 prepare？）")
        return
    geo_path = Path(geometry) if geometry else build_dir / GEOMETRY_FILE
    if not geo_path.is_file():
        result.add("warning", "geometry_missing",
                   "无建结几何 %s，结构几何复核跳过" % geo_path)
        return
    meta = metadata
    if meta is None:
        meta, error = _load_json(build_dir / METADATA_FILE)
        if error:
            result.add("error", "geometry_metadata_unreadable",
                       "%s 不可解析: %s" % (METADATA_FILE, error))
            return
    if not meta:
        result.add("error", "geometry_metadata_missing",
                   "缺 %s：无法校验表面索引/输运轴" % METADATA_FILE)
        return
    try:
        from smj.junction.validate import validate_junction_geometry  # lazy: vendored
        verdict = validate_junction_geometry(build_dir, metadata=dict(meta),
                                             geometry=geo_path)
    except Exception as exc:  # noqa: BLE001
        result.add("error", "geometry_check_failed", "结构几何检查失败: %s" % exc)
        return
    checks = []
    for check in verdict.get("checks", []):
        # 连接核对由本模块单独复核（避免同一事实两条结论）。
        if check.get("name") in ("junction_connection_check", "junction_surface_check",
                                 "junction_a2_confirmation"):
            continue
        checks.append(check)
        if check.get("passed"):
            continue
        severity = {"error": "error", "warn": "warning"}.get(check.get("severity"), "info")
        result.add(severity, "structure:%s" % check.get("name", "unknown"),
                   str(check.get("detail", "")), value=check.get("value"))
    result.structure = {
        "passed": verdict.get("passed"),
        "errors": verdict.get("errors"),
        "warnings": verdict.get("warnings"),
        "geometry": str(geo_path),
        "metadata": str(build_dir / METADATA_FILE) if meta else None,
        "checks": checks,
    }
    if verdict.get("passed") and not verdict.get("warnings"):
        result.add("info", "structure_ok", "结构几何闸门通过")
    _anchor_report(result, geo_path, meta)


def _anchor_report(result: ReviewResult, geo_path: Path,
                   meta: Mapping[str, Any]) -> None:
    """锚定接触几何（M–X 距离 / C–X–M 角）——只报告，不判定。

    ``validate_junction_geometry`` 已核 Au–S 锚定是否在带内；这里补的是**角度**：
    S 若接近线形（~180°）通常意味着过度拉伸/插进金属，C–S–M 明显小于 90° 则可能
    被相邻原子挤压。合格区间取决于锚定化学与表面位点，所以只写 ``info`` + JSON。
    """
    try:
        from smj.junction.geometry import read_aims_geometry  # lazy: vendored

        from smj.checks import anchor_contacts

        atoms = read_aims_geometry(geo_path)
    except Exception:  # noexcept: 几何不可解析时 validate 已经报过
        return
    metal = "Au"
    for key in ("electrode", "electrode_element", "metal"):
        value = meta.get(key) if isinstance(meta, Mapping) else None
        if isinstance(value, str) and value:
            metal = value
        elif isinstance(value, Mapping) and value.get("element"):
            metal = str(value["element"])
    contacts = anchor_contacts(atoms, metal_element=metal)
    result.structure["anchors"] = contacts
    if not contacts:
        result.add("info", "structure:anchor_contacts",
                   "未找到贴着 %s 的锚定原子（S/N/Se）：接触几何无法报告" % metal)
        return
    parts = []
    for item in contacts:
        text = "%s%d–%s%d %.2f Å" % (
            item["anchor_element"], item["anchor"], item["metal_element"],
            item["metal"], item["metal_distance_A"])
        if item["angle_deg"] is not None:
            text += "，%s%d–%s%d–%s%d %.1f°" % (
                item["neighbor_element"], item["neighbor"], item["anchor_element"],
                item["anchor"], item["metal_element"], item["metal"], item["angle_deg"])
        parts.append(text)
    result.add("info", "structure:anchor_contacts",
               "锚定接触：%s（只报告；合格区间取决于锚定化学与位点）" % "；".join(parts))


# --------------------------------------------------------------------------- #
# 连接核对（connection check；2026-11-15 前叫 A2 人工确认）
# --------------------------------------------------------------------------- #
def _connection_check_state(build_dir: Path | None) -> dict:
    state: dict[str, Any] = {"present": False, "valid": False, "confirmation": None,
                             "audit": None, "recorded": None, "current": None,
                             "deviation": None, "reason": None}
    if build_dir is None:
        state["reason"] = "无 electrode_build 阶段目录"
        return state
    confirm_path = _first_existing(build_dir, (CONFIRM_FILE,) + LEGACY_CONFIRM_FILES)
    if confirm_path is None:
        confirm_path = build_dir / CONFIRM_FILE
    audit_path = _first_existing(build_dir, (AUDIT_FILE,) + LEGACY_AUDIT_FILES)
    if audit_path is None:
        audit_path = build_dir / AUDIT_FILE
    metadata, _meta_error = _load_json(build_dir / METADATA_FILE)
    state["current_indices"] = _indices(metadata) if metadata else None
    if audit_path.is_file():
        audit = None
        try:
            audit = json.loads(audit_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            audit = None
        state["audit"] = {"path": str(audit_path),
                          "ok": bool((audit or {}).get("ok")),
                          "checks": (audit or {}).get("checks")}
    if not confirm_path.is_file():
        state["reason"] = "尚无连接核对记录（%s）" % confirm_path
        return state
    state["present"] = True
    confirmation, error = _load_json(confirm_path)
    state["confirmation"] = str(confirm_path)
    if error or confirmation is None:
        state["reason"] = "连接核对记录不可解析: %s" % (error or "非对象")
        return state
    state["recorded"] = {key: confirmation.get(key) for key in
                         ("fingerprint", "who", "when", "indices")}
    state["recorded_indices"] = _indices(confirmation)
    state["deviation"] = surface_index_deviation(confirmation, metadata)
    try:
        from smj.junction.fingerprint import check_confirmation  # lazy: vendored
        checked = check_confirmation(build_dir)
    except Exception as exc:  # noqa: BLE001
        state["reason"] = "连接核对指纹检查失败: %s" % exc
        return state
    state["valid"] = bool(checked.get("valid"))
    state["present"] = bool(checked.get("present"))
    state["fingerprint_current"] = checked.get("current")
    state["fingerprint_recorded"] = checked.get("recorded")
    state["reason"] = checked.get("reason")
    return state


def _auto_connection_check(result: ReviewResult, run_dir: Path,
                        build_dir: Path | None, config: ReviewConfig,
                        *, force: bool = False) -> None:
    """默认自动登记连接核对（``ReviewConfig.auto_connection_check``）。

    审计通过才签；无 matplotlib 时退化为"只记指纹不出图"（图只是证据留存，
    自动模式本来就没有人工看图这一步）。拒绝时只记 warning —— 若闸门开着，
    随后的 connection-check 复核会把它升级为 error，保持 fail-closed。
    """
    if build_dir is None:
        return
    current = _connection_check_state(build_dir)
    if current["valid"] and not force:
        return
    try:
        from .a2 import approve_connection_check
    except Exception as exc:  # noqa: BLE001
        result.add("warning", "connection_check_auto_failed",
                   "自动连接核对不可用: %s" % exc)
        return
    try:
        approval = approve_connection_check(run_dir, who="auto", stage=build_dir,
                                         force=force)
        if not approval.ok and approval.code == "figure_failed":
            approval = approve_connection_check(run_dir, who="auto", stage=build_dir,
                                             no_figure=True, force=force)
    except Exception as exc:  # noqa: BLE001
        result.add("warning", "connection_check_auto_failed",
                   "自动连接核对失败: %s" % exc)
        return
    if approval.ok:
        result.add("info", "connection_check_auto_signed",
                   "已自动登记连接核对（指纹 %s）"
                   % str(approval.fingerprint or "")[:12])
    else:
        result.add("warning", "connection_check_auto_refused",
                   "自动连接核对未登记（%s）：%s" % (approval.code, approval.reason),
                   code_detail=approval.code)


def _connection_check_review(result: ReviewResult, build_dir: Path | None,
                          config: ReviewConfig) -> None:
    check = _connection_check_state(build_dir)
    result.connection_check = check
    if not check["present"]:
        if config.connection_check:
            result.add("error", "connection_check_missing",
                       "报错：缺连接核对记录（%s）；connection_check=true 要求先核对表面原子"
                       "（默认自动登记，可重跑 review 或设 review.auto_connection_check: true）"
                       % CONFIRM_FILE)
        else:
            result.add("info", "connection_check_disabled",
                       "连接核对未作为闸门（review.connection_check=false）；%s"
                       % (check.get("reason") or "无核对记录"))
        return
    if check.get("audit") and not check["audit"]["ok"]:
        severity = "error" if config.connection_check else "warning"
        result.add(severity, "connection_check_audit_failed",
                   "表面审计未通过（%s）：表面三元组/与 tcontrol 一致性有问题"
                   % check["audit"]["path"])
    if check["valid"]:
        result.add("info", "connection_check_confirmed",
                   "连接核对有效（指纹 %s）"
                   % str(check.get("fingerprint_current") or "")[:12])
        return

    deviation = check.get("deviation")
    tolerance = float(config.connection_check_deviation)
    accept = (float(config.accept_deviation)
              if config.accept_deviation is not None else None)
    if deviation is None:
        severity = "error" if config.connection_check else "warning"
        result.add(severity, "connection_check_stale_unquantified",
                   "连接核对记录已失效，且表面索引缺失、无法量化偏差：%s"
                   % (check.get("reason") or "指纹不符"))
        return
    if deviation <= tolerance:
        result.add("info", "connection_check_deviation_within",
                   "连接核对指纹不符，但表面索引偏差 %.3f ≤ connection_check_deviation %.3f，"
                   "按可接受偏差继续" % (deviation, tolerance), deviation=deviation)
        return
    if accept is not None and deviation <= accept:
        result.accepted = True
        result.add("warning", "connection_check_deviation_accepted",
                   "表面索引偏差 %.3f 超过 connection_check_deviation %.3f，但 ≤ "
                   "accept_deviation %.3f：已显式接受（报告必须标注这是人工豁免，"
                   "不是连接核对有效）" % (deviation, tolerance, accept),
                   deviation=deviation)
        return
    severity = "error" if config.connection_check else "warning"
    result.add(severity, "connection_check_deviation_exceeded",
               "表面索引偏差 %.3f > connection_check_deviation %.3f（accept_deviation=%s）——"
               " 几何/表面索引已变，须按新几何重新核对"
               % (deviation, tolerance, "未设" if accept is None else "%.3f" % accept),
               deviation=deviation)


# --------------------------------------------------------------------------- #
# 运行参数溯源（P0-1）
# --------------------------------------------------------------------------- #
_NOT_FAMILY_ONLY = object()


def _family_drift(config: Any, spec: Mapping[str, Any]) -> Any:
    """协议漂移是否只来自 ``molecule.family``；是则返回 ``(run 时, 当前)``。

    ``family`` 只是判据标定标注（决定哪些判据适用），不进入任何物理计算，
    事后补登记（如把 ``null`` 改成 ``BPY``）不应把整条 run 判成"参数变了"。
    验证方式：把 run_spec 快照里的 family 塞回当前配置重算协议指纹，若与
    run 开始时记录的一致，就说明除 family 外的物理参数逐个相同。

    任何异常（配置不是 dataclass、快照缺 config、哈希算不出）都返回哨兵，
    按"不是只有 family 漂移"处理 —— fail-closed。
    """
    run_family: Any = None
    cfg = spec.get("config")
    if isinstance(cfg, Mapping):
        molecule = cfg.get("molecule")
        if isinstance(molecule, Mapping):
            run_family = molecule.get("family")
    stored = spec.get("protocol_sha256")
    try:
        patched = replace(
            config, molecule=replace(config.molecule, family=run_family))
        if stored and protocol_sha256(patched) == stored:
            return run_family, getattr(config.molecule, "family", None)
    except Exception:  # noqa: BLE001 - 非 dataclass / 缺字段一律 fail-closed
        return _NOT_FAMILY_ONLY
    return _NOT_FAMILY_ONLY


def _schema_shape_drift(spec: dict, config: Any) -> dict | None:
    """协议指纹不同但只是 schema 形状演进时，返回比对明细；否则 None。

    代码升级会给协议子集新增/移除只影响记录方式的键（``stable_check``、
    ``xyz``→``xyz_sha256``、``library_dir`` 等），此时哈希必然变，但物理参数
    没变。判据：两侧都有的键取值全同 + :data:`PROTOCOL_CRITICAL` 一个不缺 +
    确实存在只在单侧的键（否则哈希不同就无法解释，仍按漂移 fail-closed）。
    """
    stored_config = spec.get("config")
    if not isinstance(stored_config, dict):
        return None
    stored = {key: stored_config.get(key) or {}
              for key in ("molecule", "electrode", "transport")}
    try:
        current = protocol_components(config)
    except Exception:  # noqa: BLE001 - 复核不应因配置对象异常崩掉
        return None
    comparison = compare_protocol_configs(stored, current)
    if comparison["differences"] or comparison["critical_missing"]:
        return None
    if not (comparison["only_run"] or comparison["only_now"]):
        return None
    comparison["run_schema"] = spec.get("protocol_schema")
    comparison["now_schema"] = PROTOCOL_SCHEMA
    return comparison


def _provenance_review(result: ReviewResult, run_dir: Path, config: Any) -> None:
    """检查 ``run_spec.json``：缺失/重建只告警，协议漂移是阻断项。"""
    spec = load_run_spec(run_dir)
    if spec is None:
        result.add(
            "warning", "provenance_missing",
            "run_spec.json 缺失：无法证明本次复核对应 run 开始时的参数快照")
        return
    if spec.get("reconstructed"):
        result.add(
            "warning", "provenance_reconstructed",
            "run_spec.json 是事后重建（reconstructed=true），不是原始参数快照")
    stored = spec.get("protocol_sha256")
    if config is not None:
        try:
            current = protocol_sha256(config)
        except Exception as exc:  # noqa: BLE001 - 复核不应因配置对象异常崩掉
            result.add("warning", "provenance_check_failed",
                       "无法计算当前配置的协议指纹: %s" % exc)
            return
        if stored and current != stored:
            drift = _family_drift(config, spec)
            if drift is not _NOT_FAMILY_ONLY:
                family_run, family_now = drift
                result.add(
                    "warning", "provenance_family_relabel",
                    "run_spec.protocol_sha256=%s，当前配置=%s；两者差异仅在 "
                    "molecule.family（run 时 %r → 当前 %r）。family 只决定判据"
                    "适用族、不改变物理协议（几何/电极/输运参数一致），本次复核"
                    "按当前标注重算判据"
                    % (str(stored)[:12], current[:12], family_run, family_now),
                    run_spec_protocol=stored, current_protocol=current,
                    family_run=family_run, family_current=family_now)
            else:
                shaped = _schema_shape_drift(spec, config)
                if shaped is not None:
                    result.add(
                        "warning", "provenance_schema_evolved",
                        "run_spec.protocol_sha256=%s，当前配置=%s；协议 schema 已演进"
                        "（run 时 %s → 当前 %s），两侧都存在的 %d 个键取值完全相同，"
                        "仅 %d 个键只在 run 侧、%d 个只在当前侧。属记录方式变化，"
                        "不是物理参数变化"
                        % (str(stored)[:12], current[:12],
                           shaped.get("run_schema"), shaped.get("now_schema"),
                           shaped["shared"], len(shaped["only_run"]),
                           len(shaped["only_now"])),
                        run_spec_protocol=stored, current_protocol=current,
                        **shaped)
                else:
                    result.add(
                        "error", "provenance_drift",
                        "run_spec.protocol_sha256=%s，当前配置=%s；配置在 run 开始后发生"
                        "物理相关变化，本次结果不能按原参数解读（需重建 run 或显式豁免）"
                        % (str(stored)[:12], current[:12]),
                        run_spec_protocol=stored, current_protocol=current)


def _spin_review(result: ReviewResult, run_dir: Path, config: Any) -> None:
    """开壳层的可读性：多重度、是否钉住总自旋、α/β 通道是否还在。

    只产出展示字段与 WARN/INFO，不改任何数值判定 —— 判据族对开壳层未校准，
    这一点本身要作为 warning 让评审者看到，而不是藏在 results.json 里。
    """
    molecule = getattr(config, "molecule", None)
    multiplicity = getattr(molecule, "multiplicity", None)
    charge = getattr(molecule, "charge", None)
    fixed = getattr(molecule, "fixed_spin_moment", None)
    broken = getattr(molecule, "broken_symmetry", None)
    if multiplicity is None or charge is None:
        spec_mol = ((load_run_spec(run_dir) or {}).get("config") or {}).get("molecule") or {}
        if multiplicity is None:
            multiplicity = spec_mol.get("multiplicity")
        if charge is None:
            charge = spec_mol.get("charge")
        if fixed is None:
            fixed = spec_mol.get("fixed_spin_moment")
        if broken is None:
            broken = spec_mol.get("broken_symmetry")
    try:
        multiplicity = int(multiplicity) if multiplicity is not None else None
    except (TypeError, ValueError):
        multiplicity = None
    broken = bool(broken)
    if broken:
        # 破缺对称单重态不写 fixed_spin_moment（M−1=0 会把态压成闭壳层），
        # 因此不能沿用 MoleculeConfig 的默认 True，否则报告会谎称"已固定总自旋"。
        fixed = None
    elif fixed is None:
        # provenance 把默认值 True 省掉了，缺省即按 True 解释。
        fixed = True
    te = run_dir / "aitranss_te" / "TE.dat"
    te_present = te.is_file()
    channels = None
    if te_present:
        try:
            from .landauer_spin import read_te_channels
            _, channels = read_te_channels(te)
        except (OSError, ValueError) as exc:
            result.add("info", "spin_te_unreadable",
                       "TE.dat 自旋通道不可读，自旋占比不可得: %s" % exc)
    raw: Any = None
    results_path = run_dir / "junction_analyze" / "results.json"
    if results_path.is_file():
        try:
            payload = json.loads(results_path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            result.add("info", "spin_results_unreadable",
                       "junction_analyze/results.json 不可解析，自旋块不可得: %s" % exc)
            payload = None
        if isinstance(payload, Mapping):
            raw = payload.get("spin")
    if raw is None and channels is not None:
        raw = {"channels": channels, "source": str(te)} if te_present else None
    summary = normalize_spin(raw, multiplicity=multiplicity, charge=charge,
                             fixed_spin_moment=fixed if ((multiplicity or 1) > 1 or broken)
                             else None,
                             broken_symmetry=broken)
    if summary is not None:
        result.spin = summary
    for item in spin_findings(summary, multiplicity=multiplicity, te_present=te_present,
                              broken_symmetry=broken):
        extra = {k: v for k, v in item.items()
                 if k not in ("severity", "code", "message")}
        result.add(item["severity"], item["code"], item["message"], **extra)


# --------------------------------------------------------------------------- #
# 主入口
# --------------------------------------------------------------------------- #
def review_run(target: Any, *, force: bool = True,
               accept_deviation: float | None = None,
               output: str | Path | None = None,
               geometry: str | Path | None = None,
               metadata: Mapping[str, Any] | None = None,
               stage: str | None = None,
               review_config: ReviewConfig | Mapping[str, Any] | None = None,
               config: Any = None) -> ReviewResult:
    """执行结构/A2 复核并写 ``review.json``。

    Parameters
    ----------
    target:
        run 目录，或 :class:`smj.config.RunConfig`（取其 ``review`` 配置与
        ``run_dir``）。
    force:
        ``review.json`` 已存在时，``False`` 跳过写入（返回已有路径）；默认覆盖。
    accept_deviation:
        显式人工豁免阈值（覆盖 ``ReviewConfig.accept_deviation``）。
    """
    run_dir, cfg = _run_dir_of(target)
    if cfg is None and config is not None:
        cfg = config
    if cfg is None:
        cfg_path = run_dir / "smj.yaml"
        if cfg_path.is_file():
            try:
                from .config import load_config
                cfg = load_config(cfg_path)
            except Exception:  # noqa: BLE001 —— 旧 run 配置可能不再合法
                cfg = None

    review_cfg = getattr(cfg, "review", None) if cfg is not None else None
    if review_config is not None:
        if isinstance(review_config, Mapping):
            review_cfg = ReviewConfig(**{k: v for k, v in review_config.items()
                                         if k in ReviewConfig.__dataclass_fields__})
        else:
            review_cfg = review_config
    review_cfg = review_cfg or ReviewConfig()
    if accept_deviation is not None:
        review_cfg = ReviewConfig(connection_check=review_cfg.connection_check,
                                  connection_check_deviation=review_cfg.connection_check_deviation,
                                  auto_connection_check=review_cfg.auto_connection_check,
                                  accept_deviation=float(accept_deviation))

    result = ReviewResult(run_id=getattr(cfg, "name", None) or run_dir.name)
    build_dir = None
    if stage:
        candidate = Path(stage)
        if not candidate.is_absolute():
            candidate = run_dir / candidate
        if candidate.is_dir():
            build_dir = candidate
        else:
            # 旧名 ``au60_build`` / 尺寸名 ``au20_build`` 走兼容回退。
            try:
                key = normalize_stage(str(stage))
            except SmjError:
                key = None
            if key is not None:
                fallback = RunLayout(run_dir).stage_dir(key)
                if fallback.is_dir():
                    build_dir = fallback
    if build_dir is None:
        build_dir = find_build_dir(run_dir)
    if review_cfg.auto_connection_check:
        _auto_connection_check(result, run_dir, build_dir, review_cfg)
    _structure_review(result, build_dir, geometry=geometry, metadata=metadata)
    _connection_check_review(result, build_dir, review_cfg)
    _provenance_review(result, run_dir, cfg)
    _spin_review(result, run_dir, cfg)
    result.finalize()
    # 供 Mapping 接口（CLI 的 _as_mapping 走 dict(result)）拿到本次实际配置。
    result.config_used = review_cfg

    target_path = Path(output) if output else (run_dir / "review.json")
    result.review_path = str(target_path)
    result.written = not (target_path.exists() and not force)
    payload = result.to_dict(review_cfg)
    if result.written:
        result.review_path = str(_write_json(payload, target_path))
    return result


#: 兼容别名（CLI 首选名是 review_run；备用 run_review/review）。
review = review_run
run_review = review_run
