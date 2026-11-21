"""完工度闸门（2026-11-20 复核 P1-6 / R7，P2-6 / R6）。

`validate_junction_run` 的检查项由"发现了哪些产物"驱动：只放一个
``aitranss_te/TE.dat`` 的空 run，缺失的阶段根本不产生 check，于是
``errors=0``、``engineering_passed=true``、甚至顶层 ``passed=true``。
判定只看"出现的产物是否合格"，不看"该跑的阶段是否跑完"。

这里补一条**合成判据** ``run_incomplete``（engineering/error）：只要 run 没有
走完 ``STAGE_ORDER``，它就不过，并让 ``passed``/``engineering_passed`` 同步为
false；同时把 ``completeness`` 摘要挂到 verdict 上，报告渲染可以直说
"缺哪些阶段 / 哪些阶段已采集待验收"。

2026-11-20 复审后的两条口径：

* R7：**工程错误、物理失败、警告分开计数**。``engineering_passed`` 只看工程类
  error 级失败；物理判据失败不再把工程判定拖成 false（它另有 ``physics_verdict``）。
  ``errors`` 仍保留"所有 error 级失败"的总体语义，``passed`` 仍等价于 ``errors==0``。
* R6：``finished`` 只表示远端作业结束、产物已取回，**不等于验收通过**；完工度只认
  ``ACCEPTED_STATUSES = ("collected", "skipped")``。``finished`` 进
  ``stages_awaiting``（待验收），``skipped`` 视为"有替代证据已验收"。

``physics_verdict`` 不在这里重算（它是物理闸门的归约，与完工度无关）。
"""
from __future__ import annotations

import hashlib
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from smj.junction.criteria import PHYSICAL_NAMES

CHECK_NAME = "run_incomplete"

#: 报告口径的"已验收"状态：``collected`` 是通过校验收下的，``skipped`` 有替代证据。
ACCEPTED_STATUSES: tuple[str, ...] = ("collected", "skipped")

#: 作业已结束但产物还没走 collect/验收——报告里显示为"待采集/待验收"。
AWAITING_STATUSES: tuple[str, ...] = ("finished",)

#: 判定口径的 schema 版本。语义变更时递增；源文件指纹负责捕获实现变更。
JUDGMENT_SCHEMA = "smj.judgment.v1"


@lru_cache(maxsize=1)
def judgment_fingerprint() -> str:
    """返回判定实现（criteria/validate/verdict 源码）的 sha256。

    manifest 里保存这个指纹，报告重跑时据此判断"旧判定是否仍由当前代码
    产生"。仅比较 ``generated_at`` 不够：同一份 manifest 可能由判据已变化的
    新版本复用，从而把旧结论（例如 passed）原样发布出去。

    指纹覆盖实际决定判定的三个模块；任一文件变化都会让旧 manifest 过期。
    读不到某个源文件时退化为已读文件之和，不抛错（打包/冻结环境）。
    """
    digest = hashlib.sha256()
    digest.update(JUDGMENT_SCHEMA.encode("ascii"))
    sources: dict[str, Path] = {str(Path(__file__).resolve()): Path(__file__)}
    for module_name in ("criteria", "validate"):
        try:
            module = __import__("smj.junction.%s" % module_name,
                                fromlist=[module_name])
            path = Path(getattr(module, "__file__", ""))
        except Exception:  # noqa: BLE001 - 指纹不能成为导入失败的替罪羊
            continue
        if path.is_file():
            sources[str(path.resolve())] = path
    for key in sorted(sources):
        source = sources[key]
        digest.update(b"\0")
        digest.update(source.name.encode("utf-8"))
        digest.update(b"\0")
        try:
            digest.update(source.read_bytes())
        except OSError:
            digest.update(b"<unreadable>")
    return digest.hexdigest()


def _is_physical(check: Mapping[str, Any]) -> bool:
    """判据是否属于物理闸门：优先看 check 自带的 ``kind``，退回判据注册表。"""
    if check.get("kind") == "physical":
        return True
    if check.get("kind") == "engineering":
        return False
    return str(check.get("name")) in PHYSICAL_NAMES


def completeness_check(*, is_done: bool, first_incomplete: str | None,
                       missing: Sequence[str],
                       awaiting: Sequence[str] = ()) -> dict[str, Any]:
    """生成 ``run_incomplete`` 合成判据。"""
    if is_done:
        detail = "全部阶段已验收"
    else:
        parts = ["run 未走完"]
        if first_incomplete:
            parts.append("首个未完成阶段 %s" % first_incomplete)
        if awaiting:
            parts.append("已采集待验收 %s" % ", ".join(awaiting))
        if missing:
            parts.append("缺少阶段 %s" % ", ".join(missing))
        detail = "；".join(parts)
    return {
        "name": CHECK_NAME,
        "kind": "engineering",
        "severity": "error",
        "passed": bool(is_done),
        "detail": detail,
        "applicable": True,
        "evidence": None,
    }


def _applicable(check: Mapping[str, Any]) -> bool:
    return check.get("applicable") is not False


def _failed(check: Mapping[str, Any]) -> bool:
    return _applicable(check) and not check.get("passed")


def apply_completeness(verdict: Mapping[str, Any], *,
                       is_done: bool, first_incomplete: str | None,
                       missing: Sequence[str],
                       awaiting: Sequence[str] = ()) -> dict[str, Any]:
    """把完工度判据并进 vendored verdict，并重算工程判定。

    返回**新** dict（不改入参）。``checks`` 里追加合成判据；
    ``errors``/``warnings``/``passed``/``engineering_passed`` 按同一口径重算，
    并额外给出 ``engineering_errors`` / ``physics_failures`` 两个分账计数。
    """
    checks = [dict(c) for c in (verdict.get("checks") or [])]
    if not any(c.get("name") == CHECK_NAME for c in checks):
        checks.append(completeness_check(is_done=is_done,
                                         first_incomplete=first_incomplete,
                                         missing=missing, awaiting=awaiting))

    error_checks = [c for c in checks
                    if c.get("severity") == "error" and _failed(c)]
    engineering_errors = [c for c in error_checks if not _is_physical(c)]
    physics_failures = [c for c in error_checks if _is_physical(c)]
    # 源头的 severity 是 "warn"（validate 的 _warn）；历史 verdict 里也有
    # "warning"，两种都认，否则工程的 warn 级判据永远不会出现在 warnings 里。
    warnings = [c for c in checks
                if c.get("severity") in ("warn", "warning") and _failed(c)]

    failed = [str(x) for x in (verdict.get("failed_criteria") or [])]
    if not is_done and CHECK_NAME not in failed:
        failed.append(CHECK_NAME)

    out = dict(verdict)
    out["checks"] = checks
    out["errors"] = len(error_checks)
    out["engineering_errors"] = len(engineering_errors)
    out["physics_failures"] = len(physics_failures)
    out["warnings"] = len(warnings)
    out["passed"] = not error_checks
    out["engineering_passed"] = not engineering_errors
    out["failed_criteria"] = failed
    out["completeness"] = {
        "is_done": bool(is_done),
        "first_incomplete": first_incomplete,
        "stages_missing": list(missing),
        "stages_awaiting": list(awaiting),
    }
    return out


def stage_acceptance(statuses: Mapping[str, str],
                     expected: Iterable[str]) -> tuple[list[str], list[str]]:
    """按 ``STAGE_ORDER`` 分出"已采集待验收"与"还没有产物"两类阶段。

    ``finished`` 只说明远端作业结束，产物可能还没取回/没跑 collect——放进
    ``awaiting``，不算完成；``skipped`` 有替代证据，算已验收。
    """
    awaiting: list[str] = []
    missing: list[str] = []
    for name in expected:
        status = statuses.get(name)
        if status in ACCEPTED_STATUSES:
            continue
        if status in AWAITING_STATUSES:
            awaiting.append(name)
        else:
            missing.append(name)
    return awaiting, missing


def missing_stages(statuses: Mapping[str, str],
                   expected: Iterable[str]) -> list[str]:
    """按 ``STAGE_ORDER`` 列出"没有验收"的阶段（missing/awaiting 都算）。"""
    awaiting, missing = stage_acceptance(statuses, expected)
    return [name for name in expected
            if name in set(awaiting) | set(missing)]
