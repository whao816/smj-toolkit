"""完工度闸门（2026-11-20 复核 P1-6 / R7）。

`validate_junction_run` 的检查项由"发现了哪些产物"驱动：只放一个
``aitranss_te/TE.dat`` 的空 run，缺失的阶段根本不产生 check，于是
``errors=0``、``engineering_passed=true``、甚至顶层 ``passed=true``。
判定只看"出现的产物是否合格"，不看"该跑的阶段是否跑完"。

这里补一条**合成判据** ``run_incomplete``（engineering/error）：只要 run 没有
走完 ``STAGE_ORDER``，它就不过，并让 ``passed``/``engineering_passed`` 同步为
false；同时把 ``completeness`` 摘要挂到 verdict 上，报告渲染可以直说
"缺哪些阶段"。

判定口径保持与 vendored 一致：``applicable is False`` 的 check 不作闸门；
``physics_verdict`` 不在这里重算（它是物理闸门的归约，与完工度无关）。
"""
from __future__ import annotations

from typing import Any, Iterable, Mapping, Sequence

from smj.junction.criteria import PHYSICAL_NAMES

CHECK_NAME = "run_incomplete"


def completeness_check(*, is_done: bool, first_incomplete: str | None,
                       missing: Sequence[str]) -> dict[str, Any]:
    """生成 ``run_incomplete`` 合成判据。"""
    if is_done:
        detail = "全部阶段已完成"
    elif missing:
        detail = "run 未走完：首个未完成阶段 %s；缺少阶段 %s" % (
            first_incomplete or "?", ", ".join(missing))
    else:
        detail = "run 未走完（首个未完成阶段 %s）" % (first_incomplete or "?")
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


def apply_completeness(verdict: Mapping[str, Any], *,
                       is_done: bool, first_incomplete: str | None,
                       missing: Sequence[str]) -> dict[str, Any]:
    """把完工度判据并进 vendored verdict，并重算工程判定。

    返回**新** dict（不改入参）。``checks`` 里追加合成判据；
    ``errors``/``warnings``/``passed``/``engineering_passed`` 按同一口径重算。
    """
    checks = [dict(c) for c in (verdict.get("checks") or [])]
    if not any(c.get("name") == CHECK_NAME for c in checks):
        checks.append(completeness_check(is_done=is_done,
                                         first_incomplete=first_incomplete,
                                         missing=missing))

    errors = sum(1 for c in checks
                 if c.get("severity") == "error" and not c.get("passed")
                 and _applicable(c))
    warnings = sum(1 for c in checks
                   if c.get("severity") == "warning" and not c.get("passed")
                   and _applicable(c))
    engineering = [c for c in checks
                   if c.get("name") not in PHYSICAL_NAMES and _applicable(c)]
    engineering_passed = not any(
        c.get("severity") == "error" and not c.get("passed") for c in engineering)

    failed = [str(x) for x in (verdict.get("failed_criteria") or [])]
    if not is_done and CHECK_NAME not in failed:
        failed.append(CHECK_NAME)

    out = dict(verdict)
    out["checks"] = checks
    out["errors"] = errors
    out["warnings"] = warnings
    out["passed"] = errors == 0
    out["engineering_passed"] = engineering_passed and errors == 0
    out["failed_criteria"] = failed
    out["completeness"] = {
        "is_done": bool(is_done),
        "first_incomplete": first_incomplete,
        "stages_missing": list(missing),
    }
    return out


def missing_stages(statuses: Mapping[str, str],
                   expected: Iterable[str]) -> list[str]:
    """按 ``STAGE_ORDER`` 列出"没有完成"的阶段（missing/skipped/failed 都算）。"""
    done = {"collected", "finished", "skipped"}
    return [name for name in expected if statuses.get(name) not in done]
