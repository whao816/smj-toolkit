"""结构 / 结果复核门（review）：A2 确认、可接受偏差与 review.json。

复核分两层，**不混在一起**：

1. **结构几何** —— 有 ``auN_build/geometry.in`` 时调 vendored
   ``validate_junction_geometry``（表面索引/输运轴/接触重叠/Au-S 锚定），
   error 级问题直接阻断；
2. **A2 人工确认** —— 读 ``auN_build/a2_confirmation.json``（vendored
   ``fingerprint.record_confirmation`` 的格式）与当前几何指纹：

   - ``ReviewConfig.require_a2=True`` 且**无 A2 数据** -> ``blocked``；
   - 确认有效（指纹一致）-> 通过；
   - 指纹失效（几何已变）时比较**表面索引偏差**
     ``1 - 交集 / 较大集合``（0..1，未匹配表面原子占比）：
     ``<= a2_deviation`` 视为可接受偏差；
     否则若显式给 ``accept_deviation`` 且 ``<= accept_deviation``，记为
     ``warning`` 并 ``accepted=True``（报告须标注）；再超 -> ``blocked``。
   - ``require_a2=False`` 时 A2 只是加分证据：缺失/失效记 info/warning，
     不阻断（但结构几何 error 仍然阻断）。

``review.json`` 是派生结果，默认覆盖（``force=False`` 才跳过）；与
``report.py`` 不覆盖外部 run 文件的策略不同，原因写在函数 docstring。

vendored 模块全部 lazy import。
"""

from __future__ import annotations

import json
import os
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .config import ReviewConfig
from .provenance import load_run_spec, protocol_sha256
from .state import now_iso

REVIEW_SCHEMA = "smj.review.v1"

#: A2 相关文件（与 vendored fingerprint/a2_viz 同名）。
CONFIRM_FILE = "a2_confirmation.json"
AUDIT_FILE = "a2_surface_audit.json"
METADATA_FILE = "electrode_metadata.json"
GEOMETRY_FILE = "geometry.in"

SIDES = ("left", "right")


# --------------------------------------------------------------------------- #
# 结果对象
# --------------------------------------------------------------------------- #
@dataclass
class ReviewResult(Mapping):
    status: str = "passed"
    issues: list[dict] = field(default_factory=list)
    a2: dict | None = None
    structure: dict | None = None
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

        A2 门的 warning（缺确认但未要求、显式接受偏差、审计异常）才把整体状态
        降为 ``warning``；``geometry_missing`` 这类"复核无法执行"也算 warning。
        """
        if any(i["severity"] == "error" for i in self.issues):
            self.status = "blocked"
            return self
        gate_warning = any(
            i["severity"] == "warning" and (
                str(i.get("code", "")).startswith("a2_")
                or str(i.get("code", "")).startswith("provenance_")
                or i.get("code") in ("geometry_missing", "geometry_metadata_missing",
                                     "geometry_metadata_unreadable", "geometry_check_failed"))
            for i in self.issues)
        self.status = "warning" if gate_warning else "passed"
        return self

    def to_dict(self, config: ReviewConfig | None = None) -> dict:
        return {
            "schema": REVIEW_SCHEMA,
            "generated_at": now_iso(),
            "run_id": self.run_id,
            "status": self.status,
            "ok": self.ok,
            "accepted": self.accepted,
            "config": {
                "require_a2": getattr(config, "require_a2", None),
                "a2_deviation": getattr(config, "a2_deviation", None),
                "accept_deviation": getattr(config, "accept_deviation", None),
            },
            "issues": list(self.issues),
            "errors": self.errors,
            "warnings": self.warnings,
            "a2": self.a2,
            "structure": self.structure,
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


def find_build_dir(run_dir: str | Path) -> Path | None:
    """定位 A2 确认的对象：``auN_build``（默认 ``au60_build``）。

    2026-11-13 起 ``au_relax`` 已从主链路移除，A2 的对象就是贴装后的建结几何；
    旧 run（只有 ``au_relax`` 而没有 build 目录）作为兼容回退。
    """
    root = Path(run_dir)
    if not root.is_dir():
        return None
    preferred = root / "au60_build"
    legacy = root / "au_relax"
    if preferred.is_dir():
        # 旧 run 的 A2 签在 au_relax 几何上；若新口径的 au60 还没有签署，
        # 仍以旧签署对象为准，避免显式兼容旧 run 时误报 a2_missing。
        if ((legacy / "a2_confirmation.json").is_file()
                and not (preferred / "a2_confirmation.json").is_file()):
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
                   "无 auN_build 阶段目录，结构几何复核跳过（阶段未 prepare？）")
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
        from smj._vendor.junction.validate import validate_junction_geometry  # lazy: vendored
        verdict = validate_junction_geometry(build_dir, metadata=dict(meta),
                                             geometry=geo_path)
    except Exception as exc:  # noqa: BLE001
        result.add("error", "geometry_check_failed", "结构几何检查失败: %s" % exc)
        return
    checks = []
    for check in verdict.get("checks", []):
        # A2 由本模块单独复核（避免同一事实两条结论）。
        if check.get("name") == "junction_a2_confirmation":
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


# --------------------------------------------------------------------------- #
# A2 复核
# --------------------------------------------------------------------------- #
def _a2_state(build_dir: Path | None) -> dict:
    state: dict[str, Any] = {"present": False, "valid": False, "confirmation": None,
                             "audit": None, "recorded": None, "current": None,
                             "deviation": None, "reason": None}
    if build_dir is None:
        state["reason"] = "无 auN_build 阶段目录"
        return state
    confirm_path = build_dir / CONFIRM_FILE
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
        state["reason"] = "尚无 A2 确认记录（%s）" % confirm_path
        return state
    state["present"] = True
    confirmation, error = _load_json(confirm_path)
    state["confirmation"] = str(confirm_path)
    if error or confirmation is None:
        state["reason"] = "A2 确认记录不可解析: %s" % (error or "非对象")
        return state
    state["recorded"] = {key: confirmation.get(key) for key in
                         ("fingerprint", "who", "when", "indices")}
    state["recorded_indices"] = _indices(confirmation)
    state["deviation"] = surface_index_deviation(confirmation, metadata)
    try:
        from smj._vendor.junction.fingerprint import check_confirmation  # lazy: vendored
        checked = check_confirmation(build_dir)
    except Exception as exc:  # noqa: BLE001
        state["reason"] = "A2 指纹检查失败: %s" % exc
        return state
    state["valid"] = bool(checked.get("valid"))
    state["present"] = bool(checked.get("present"))
    state["fingerprint_current"] = checked.get("current")
    state["fingerprint_recorded"] = checked.get("recorded")
    state["reason"] = checked.get("reason")
    return state


def _a2_review(result: ReviewResult, build_dir: Path | None,
               config: ReviewConfig) -> None:
    a2 = _a2_state(build_dir)
    result.a2 = a2
    if not a2["present"]:
        if config.require_a2:
            result.add("error", "a2_missing",
                       "ReviewConfig.require_a2=True，但缺 A2 人工确认记录（%s）："
                       "须确认表面原子并把指纹写入 %s 后重新复核"
                       % (CONFIRM_FILE, CONFIRM_FILE))
        else:
            result.add("info", "a2_not_required",
                       "A2 未要求（ReviewConfig.require_a2=False）；%s"
                       % (a2.get("reason") or "无确认记录"))
        return
    if a2.get("audit") and not a2["audit"]["ok"]:
        severity = "error" if config.require_a2 else "warning"
        result.add(severity, "a2_audit_failed",
                   "A2 表面审计未通过（%s）：表面三元组/与 tcontrol 一致性有问题"
                   % a2["audit"]["path"])
    if a2["valid"]:
        result.add("info", "a2_confirmed",
                   "A2 确认有效（指纹 %s）"
                   % str(a2.get("fingerprint_current") or "")[:12])
        return

    deviation = a2.get("deviation")
    tolerance = float(config.a2_deviation)
    accept = (float(config.accept_deviation)
              if config.accept_deviation is not None else None)
    if deviation is None:
        severity = "error" if config.require_a2 else "warning"
        result.add(severity, "a2_stale_unquantified",
                   "A2 确认已失效，且表面索引缺失、无法量化偏差：%s"
                   % (a2.get("reason") or "指纹不符"))
        return
    if deviation <= tolerance:
        result.add("info", "a2_deviation_within",
                   "A2 指纹不符，但表面索引偏差 %.3f ≤ a2_deviation %.3f，按可接受偏差继续"
                   % (deviation, tolerance), deviation=deviation)
        return
    if accept is not None and deviation <= accept:
        result.accepted = True
        result.add("warning", "a2_deviation_accepted",
                   "表面索引偏差 %.3f 超过 a2_deviation %.3f，但 ≤ accept_deviation %.3f："
                   "已显式接受（报告必须标注这是人工豁免，不是 A2 有效）"
                   % (deviation, tolerance, accept), deviation=deviation)
        return
    severity = "error" if config.require_a2 else "warning"
    result.add(severity, "a2_deviation_exceeded",
               "表面索引偏差 %.3f > a2_deviation %.3f（accept_deviation=%s）——"
               " 几何/表面索引已变，须按新几何重新做 A2 确认"
               % (deviation, tolerance, "未设" if accept is None else "%.3f" % accept),
               deviation=deviation)


# --------------------------------------------------------------------------- #
# 运行参数溯源（P0-1）
# --------------------------------------------------------------------------- #
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
            result.add(
                "error", "provenance_drift",
                "run_spec.protocol_sha256=%s，当前配置=%s；配置在 run 开始后发生"
                "物理相关变化，本次结果不能按原参数解读（需重建 run 或显式豁免）"
                % (str(stored)[:12], current[:12]),
                run_spec_protocol=stored, current_protocol=current)


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
        review_cfg = ReviewConfig(require_a2=review_cfg.require_a2,
                                  a2_deviation=review_cfg.a2_deviation,
                                  accept_deviation=float(accept_deviation))

    result = ReviewResult(run_id=getattr(cfg, "name", None) or run_dir.name)
    build_dir = None
    if stage:
        candidate = Path(stage)
        if not candidate.is_absolute():
            candidate = run_dir / candidate
        if candidate.is_dir():
            build_dir = candidate
    if build_dir is None:
        build_dir = find_build_dir(run_dir)
    _structure_review(result, build_dir, geometry=geometry, metadata=metadata)
    _a2_review(result, build_dir, review_cfg)
    _provenance_review(result, run_dir, cfg)
    result.finalize()

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
