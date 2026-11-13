"""A2 人工确认（签署）的打包侧入口（F4，2026-11-11 本地补丁）。

**缺口**：独立包此前没有任何签署命令。``smj review`` 只**读取**
``electrode_build/a2_confirmation.json``；真正能登记的路径是 vendored 的
``python -m smj._vendor.junction.a2_viz --stage <electrode_build> --sign-who <名字>``，
既没写进 README，review 的提示还指向主仓的 ``agenticdft smj-review``
（独立包里不存在）。于是"红线"在默认配置下（``review.require_a2: false``）
根本不出现在任何路径上。

本模块把 vendored 的 ``surface_audit`` / ``render_figure`` /
``record_confirmation`` 组装成一条 fail-closed 的 CLI 路径
（``smj review --approve --who <名字>``）：

1. 审计不过 -> 拒绝签署（``refused``，CLI 退出码 2）：先修几何/表面索引；
2. 出图失败（如未装 matplotlib）-> 拒绝签署 —— 签署的语义是"人已经看过
   表面三元组图"，不能在没有图的情况下假装看过；``no_figure=True`` 是显式
   声明"已在别处看过图"的逃生口；
3. 已有一份**仍然有效**的确认（指纹与当前 geometry.in 一致）时默认不覆盖，
   ``force=True`` 才重签；
4. 记录绑定的指纹以**当前** ``geometry.in`` 现算，由 vendored
   ``fingerprint`` 负责（几何一变旧确认立刻失效）。

纯文件操作，不跑任何作业。
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .config import RunConfig
from .errors import SmjError
from .review import AUDIT_FILE, CONFIRM_FILE, GEOMETRY_FILE, METADATA_FILE, find_build_dir
from .state import RunLayout, normalize_stage, now_iso

A2_SCHEMA = "smj.a2.v1"

#: CLI 退出码：审计/出图不通过时用 2（与 vendored ``a2_viz`` 的 --sign-who 一致），
#: 与 review 的"复核不通过 = 1"区分开。
EXIT_REFUSED = 2


@dataclass
class A2Approval:
    """``approve_a2`` 的返回值；``ok`` 表示"确认记录已就位"（含本来就有）。"""

    status: str = "refused"          # signed / already_signed / refused
    ok: bool = False
    refused: bool = False
    code: str = ""
    reason: str = ""
    run_id: str | None = None
    build_dir: str | None = None
    who: str | None = None
    when: str | None = None
    fingerprint: str | None = None
    indices: dict[str, Any] | None = None
    audit_ok: bool | None = None
    checks_failed: list[str] = field(default_factory=list)
    audit_path: str | None = None
    png_path: str | None = None
    generated_at: str = field(default_factory=now_iso)

    @property
    def exit_code(self) -> int:
        if self.ok:
            return 0
        return EXIT_REFUSED if self.refused else 1

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": A2_SCHEMA,
            "status": self.status,
            "ok": self.ok,
            "refused": self.refused,
            "code": self.code,
            "reason": self.reason,
            "run_id": self.run_id,
            "build_dir": self.build_dir,
            "who": self.who,
            "when": self.when,
            "fingerprint": self.fingerprint,
            "indices": self.indices,
            "audit_ok": self.audit_ok,
            "checks_failed": list(self.checks_failed),
            "audit_path": self.audit_path,
            "png_path": self.png_path,
            "exit_code": self.exit_code,
            "generated_at": self.generated_at,
        }

    def __getitem__(self, key: str) -> Any:
        return self.to_dict()[key]

    def get(self, key: str, default: Any = None) -> Any:
        return self.to_dict().get(key, default)


# --------------------------------------------------------------------------- #
# 工具
# --------------------------------------------------------------------------- #
def _run_dir_of(target: Any) -> tuple[Path, RunConfig | None]:
    """接受 run 目录或 :class:`smj.config.RunConfig`。"""
    if isinstance(target, (str, Path)):
        return Path(target), None
    run_dir = getattr(target, "run_dir", None)
    if run_dir is None:
        raise SmjError("approve_a2 需要 run 目录或 RunConfig（缺 run_dir）")
    return Path(run_dir), target if isinstance(target, RunConfig) else None


def _resolve_build_dir(run_dir: Path, stage: str | Path | None) -> Path | None:
    if stage is not None:
        candidate = Path(stage)
        if not candidate.is_absolute():
            candidate = run_dir / candidate
        if candidate.is_dir():
            return candidate
        # ``--stage au60_build`` / ``au20_build`` 这类旧名或尺寸名，按
        # RunLayout 的兼容回退定位到规范目录。
        try:
            key = normalize_stage(str(stage))
        except SmjError:
            return None
        fallback = RunLayout(run_dir).stage_dir(key)
        return fallback if fallback.is_dir() else None
    return find_build_dir(run_dir)


def _refuse(code: str, reason: str, **fields: Any) -> A2Approval:
    return A2Approval(status="refused", ok=False, refused=True, code=code,
                      reason=reason, **fields)


def _write_json(payload: dict, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
                   encoding="utf-8", newline="\n")
    tmp.replace(path)
    return path


# --------------------------------------------------------------------------- #
# 主入口
# --------------------------------------------------------------------------- #
def approve_a2(target: Any, *, who: str, stage: str | Path | None = None,
               when: str | None = None, png: str | Path | None = None,
               no_figure: bool = False,
               tcontrol: str | Path | None = None,
               force: bool = False) -> A2Approval:
    """登记 electrode_build 的 A2 人工确认；审计/出图不通过时拒绝。

    Parameters
    ----------
    target:
        run 目录或 :class:`smj.config.RunConfig`（取其 ``run_dir``）。
    who:
        签署人（必填，写进 ``a2_confirmation.json``）。
    stage:
        ``electrode_build`` 目录（相对 run 目录或绝对路径）；缺省用
        :func:`smj.review.find_build_dir` 自动定位。
    no_figure:
        不出图直接签署（显式声明"已在别处看过图"；无 matplotlib 时用）。
    """
    who = str(who or "").strip()
    if not who:
        raise SmjError("A2 签署需要签署人：`smj review --approve --who <名字>`")

    run_dir, config = _run_dir_of(target)
    run_id = (getattr(config, "name", None) if config is not None else None) or run_dir.name
    base = {"run_id": run_id}

    build_dir = _resolve_build_dir(run_dir, stage)
    if build_dir is None:
        return _refuse(
            "no_build_dir",
            "找不到 electrode_build 阶段目录（A2 确认的对象是建结几何）：%s" % run_dir,
            **base)
    geometry = build_dir / GEOMETRY_FILE
    metadata_path = build_dir / METADATA_FILE
    base.update(build_dir=str(build_dir))
    if not geometry.is_file():
        return _refuse("no_geometry", "缺建结几何 %s" % geometry, **base)
    if not metadata_path.is_file():
        return _refuse("no_metadata",
                       "缺 %s：无法校验表面索引/输运轴，不能签署" % metadata_path, **base)
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        return _refuse("metadata_unreadable",
                       "%s 不可解析: %s" % (metadata_path, exc), **base)
    if not isinstance(metadata, dict) or not metadata:
        return _refuse("metadata_unreadable",
                       "%s 不是对象或为空" % metadata_path, **base)

    from smj._vendor.junction import a2_viz, fingerprint  # lazy: vendored

    if tcontrol is not None:
        tc_path: Path | None = Path(tcontrol)
    else:
        candidate = build_dir.parent / "aitranss_te" / "tcontrol"
        tc_path = candidate if candidate.is_file() else None

    try:
        audit = a2_viz.surface_audit(geometry, metadata,
                                     tcontrol=tc_path if tc_path and tc_path.is_file() else None)
    except Exception as exc:  # noqa: BLE001 —— 几何/metadata 不合法都算拒绝，不是崩溃
        return _refuse("audit_failed", "A2 表面审计执行失败: %s" % exc, **base)

    audit_path = _write_json(audit, build_dir / AUDIT_FILE)
    failed = [str(c.get("name")) for c in audit.get("checks", []) if not c.get("ok")]
    base.update(audit_ok=bool(audit.get("ok")), checks_failed=failed,
                audit_path=str(audit_path))
    if not audit.get("ok"):
        return _refuse(
            "audit_not_ok",
            "A2 表面审计未通过（%s）：%s；先修几何/表面索引，签署被拒绝"
            % (", ".join(failed) or "未命名检查",
               "%d 项检查未通过" % len(failed)),
            **base)

    # 已有仍有效的确认：默认不覆盖（幂等），--force 才重签。
    try:
        existing = fingerprint.check_confirmation(build_dir)
    except Exception as exc:  # noqa: BLE001
        return _refuse("check_failed", "读取现有 A2 确认失败: %s" % exc, **base)
    if existing.get("valid") and not force:
        return A2Approval(
            status="already_signed", ok=True, refused=False, code="already_signed",
            reason="已有与当前几何一致的 A2 确认（%s @ %s），未改动；重签请加 --force"
                   % (existing.get("who"), existing.get("when")),
            who=existing.get("who"), when=existing.get("when"),
            fingerprint=existing.get("current"), indices=None, **base)

    if no_figure:
        png_path = None
    else:
        png_path = Path(png) if png else (build_dir / "a2_surface_viz.png")
        try:
            a2_viz.render_figure(geometry, metadata, png_path, audit=audit,
                                 title="A2 surface-triple check: %s" % run_id)
        except RuntimeError as exc:
            return _refuse(
                "figure_failed",
                "出图失败，拒绝在无图情况下签署（人没有图 = 没有看过表面三元组）: %s；"
                "确认已在别处看过图可加 --no-figure" % exc, **base)
        except Exception as exc:  # noqa: BLE001
            return _refuse("figure_failed", "出图异常，拒绝签署: %s" % exc, **base)
        base["png_path"] = str(png_path)

    try:
        record = fingerprint.record_confirmation(build_dir, who=who, when=when)
    except (OSError, ValueError) as exc:
        return _refuse("record_failed", "写 A2 确认记录失败: %s" % exc, **base)

    return A2Approval(
        status="signed", ok=True, refused=False, code="signed",
        reason="已登记 A2 人工确认 -> %s" % (build_dir / CONFIRM_FILE),
        who=record.get("who"), when=record.get("when"),
        fingerprint=record.get("fingerprint"), indices=record.get("indices"),
        **base)


#: 兼容别名
approve = approve_a2
