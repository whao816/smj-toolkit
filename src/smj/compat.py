"""独立包与 junction 内核之间的唯一接缝。

内核原先是 ``agenticdft/agenticdft/tools/junction/`` 的副本，其中 3 行跨包 import
指向 AgenticDFT 的 ``...backends`` 与 ``..resources``；本包在导入时把它们改为指向本模块
（自 0.2.2 起内核由本仓维护，这 3 处接缝保持原样即可，改动内核时不要把它们改回
AgenticDFT 的路径）：

- :func:`get_backend` —— 独立包不用 AgenticDFT 的调度后端体系（提交由
  :mod:`smj.transport` / :mod:`smj.remote` 负责）。这里返回一个 ``local`` 后端
  占位，让 vendored ``sbatch_available()`` 给出"不适用"而不是 ImportError。
- :func:`estimate_smj_resources` / :func:`check_smj_budget` / :func:`walltime_seconds`
  —— 直接转发到 :mod:`smj.resources`（数值与引擎同源）。
"""

from __future__ import annotations

from .resources import (  # noqa: F401  (re-export for the vendored code)
    check_smj_budget,
    estimate_smj_resources,
    walltime_seconds,
)

__all__ = [
    "get_backend",
    "estimate_smj_resources",
    "check_smj_budget",
    "walltime_seconds",
]


class _LocalBackend:
    """最小后端占位：本包自己管提交，vendored preflight 只需要它可探测。"""

    name = "local"
    remote = False

    def probe(self) -> dict:
        return {"ok": True, "slurm": "", "detail": "standalone（提交由 smj.remote 负责）"}

    def __repr__(self) -> str:  # pragma: no cover - 调试用
        return "<smj.compat._LocalBackend name='local' remote=False>"


def get_backend(settings=None) -> _LocalBackend:  # noqa: ARG001 - 保持 vendored 调用签名
    """返回占位后端（忽略 ``settings``，保留原函数签名的兼容性）。"""
    return _LocalBackend()
