"""开壳层结的展示层汇总（P0）。

复核与报告要回答的是同一个问题：**这份结果是不是开壳层的、自旋通道各占多少、
这个解值不值得信**。原始数字散在两处——``run_spec`` 的 molecule 块（多重度/
电荷/是否钉住总自旋）与 ``junction_analyze/results.json`` 的 ``spin`` 块
（α/β 通道积分与 P(E)）。这里把两者压成一份稳定字段，供 ``report`` 与
``review`` 共用，避免两处各写一套口径。

2026-11-20 P0（``docs/自由基单分子结-实现对照与开发规划.md`` §6）：此前这些量
只躺在 results.json 与图里，评审者在 report/review 里看不到，容易把"总磁矩
保持"误读成"自由基在结中存活"。
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

#: 判据族未针对开壳层校准——这条提示必须跟着自旋数字一起出现。
SPIN_CAVEAT = (
    "判据族未针对开壳层校准：EXPECTED_RESONANCE 只登记闭壳层参照族，"
    "g_in_theory_band / seebeck_sign 对自由基体系只作参考，物理判定不应单独作为结论。"
)

#: 固定总自旋只证明约束被满足，不证明自旋属于分子。
FIXED_CONSTRAINT_CAVEAT = (
    "固定总自旋（fixed_spin_moment）只证明约束被满足，不证明自旋留在分子上；"
    "需与 fixed_spin_moment=false 的 release 解对照，并看分子分组的自旋投影。"
)

#: 开壳层但 TE.dat 只有一列 = 双通道信息在解析前就丢了。
SINGLE_CHANNEL_CAVEAT = (
    "开壳层体系的 TE.dat 只有 1 列：α/β 通道在写入时已合成均值，"
    "自旋极化无法从这份谱恢复；检查 FHI-aims 是否真的以 spin collinear 收敛。"
)


def _number(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if number != number or number in (float("inf"), float("-inf")):  # NaN/inf
        return None
    return number


def percent(value: Any) -> str:
    """极化率 0.0042 -> ``+0.42%``；不可用时返回 ``—``。"""
    number = _number(value)
    if number is None:
        return "—"
    return "%+.2f%%" % (100.0 * number)


def normalize(
    spin: Mapping[str, Any] | None,
    *,
    multiplicity: int | None = None,
    charge: int | None = None,
    fixed_spin_moment: bool | None = None,
) -> dict | None:
    """把 ``results.json`` 的 ``spin`` 块压成展示层字段；不是映射则返回 None。"""
    if not isinstance(spin, Mapping):
        return None
    channels = spin.get("channels")
    summary: dict[str, Any] = {
        "multiplicity": multiplicity,
        "charge": charge,
        "fixed_spin_moment": fixed_spin_moment,
        "channels": channels,
        "points": spin.get("points"),
        "temperature_K": spin.get("temperature_K"),
        "T_alpha_at_ef": _number(spin.get("T_alpha_at_ef")),
        "T_beta_at_ef": _number(spin.get("T_beta_at_ef")),
        "spin_polarization_at_ef": _number(spin.get("spin_polarization_at_ef")),
        "spin_polarization_L0": _number(spin.get("spin_polarization_L0")),
        "g_over_g0_alpha": _number(spin.get("g_over_g0_alpha")),
        "g_over_g0_beta": _number(spin.get("g_over_g0_beta")),
        "g_over_g0_total": _number(spin.get("g_over_g0_total")),
        "note": spin.get("note"),
        "caveat": SPIN_CAVEAT,
    }
    if summary["fixed_spin_moment"]:
        summary["fixed_constraint_caveat"] = FIXED_CONSTRAINT_CAVEAT
    return summary


def findings(
    summary: Mapping[str, Any] | None,
    *,
    multiplicity: int | None = None,
    open_shell: bool | None = None,
    te_present: bool = False,
) -> list[dict[str, Any]]:
    """按展示层字段给出复核提示（不改变任何数值判定）。

    ``open_shell`` 缺省时由 ``multiplicity`` 推断；返回的每项含
    ``severity``/``code``/``message``，与 :meth:`ReviewResult.add` 同形。
    """
    if open_shell is None:
        open_shell = bool(multiplicity and int(multiplicity) > 1)
    issues: list[dict[str, Any]] = []
    if open_shell:
        issues.append({
            "severity": "warning",
            "code": "open_shell_uncalibrated",
            "message": "开壳层结（multiplicity=%s）：%s" % (multiplicity, SPIN_CAVEAT),
            "multiplicity": multiplicity,
        })
        fixed = (summary or {}).get("fixed_spin_moment")
        if fixed:
            issues.append({
                "severity": "info",
                "code": "spin_fixed_constraint",
                "message": FIXED_CONSTRAINT_CAVEAT,
                "fixed_spin_moment": True,
            })
        channels = (summary or {}).get("channels")
        if summary is None:
            if te_present:
                issues.append({
                    "severity": "info",
                    "code": "spin_block_missing",
                    "message": "开壳层但 junction_analyze/results.json 没有 spin 块："
                               "报告里不会出现 α/β 通道数，自旋占比不可得",
                })
        elif channels is not None and int(channels) < 2:
            issues.append({
                "severity": "warning",
                "code": "spin_channels_missing",
                "message": SINGLE_CHANNEL_CAVEAT,
                "channels": channels,
            })
    elif summary is not None and summary.get("channels") is not None \
            and int(summary.get("channels") or 0) >= 2:
        if multiplicity is None:
            issues.append({
                "severity": "info",
                "code": "spin_channels_undeclared",
                "message": "TE.dat 有两列透射（自旋极化结果），但配置里没有 "
                           "multiplicity 记录（旧 run 缺 run_spec）：无法确认这是"
                           "开壳层还是误配的闭壳层谱，解读时需回看 control.in",
                "channels": summary.get("channels"),
            })
        else:
            issues.append({
                "severity": "info",
                "code": "spin_block_unexpected",
                "message": "配置为闭壳层（multiplicity=%s）但 TE.dat 有两列透射："
                           "输运是按自旋极化算的，检查 multiplicity/control.in 是否一致"
                           % multiplicity,
                "channels": summary.get("channels"),
            })
    return issues
