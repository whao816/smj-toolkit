"""解析 FHI-aims 弛豫日志（``aims.out``）的收敛状态。

2026-11-20 复核 P1-3：``aims_tip_relax.collect`` 原来只判 ``geometry.in.next_step``
是否存在，从没打开 ``aims.out``；``au_relax.collect`` 更严但仍不解析收敛标记，却在
metadata 里硬写 ``converged: True``。两者都让"未收敛的弛豫"一路流到建电极。

这里只做一件事：把日志里的收敛事实读出来，供各阶段拒绝 ``collected``。

判定口径（与 FHI-aims 输出一致）：

* 每步几何更新末尾会打印 ``Present geometry is converged.`` 或
  ``Present geometry is not yet converged.``；**最后一次**出现的标记才是最终状态。
* 正常结束会打印 ``Have a nice day.``；没有它说明作业被截断/被杀死。
* 步数上限会打印 ``Reached the maximum number of relaxation steps`` 一类信息；
  注意启动时的 ``No maximum number of relaxation steps, defaulting to 1000`` 是
  "未设上限"的提示，**不是**命中上限，必须排除。

大日志（demo 里单个 ``aims.out`` 已有 20 MB，生产上常见数百 MB）只读尾部
``_TAIL_BYTES``；收敛标记总在最后几步里，尾部足够，且记录 ``truncated`` 以便
报错时说明"只看了尾部"。
"""

from __future__ import annotations

from pathlib import Path

CONVERGED_MARKER = "Present geometry is converged."
NOT_CONVERGED_MARKER = "Present geometry is not yet converged."
NORMAL_END_MARKER = "Have a nice day."

#: 命中步数上限的提示（大小写不敏感匹配）；必须排除 "No maximum number ... defaulting"。
_STEP_LIMIT_HINTS = (
    "reached the maximum number of relaxation",
    "reached maximum number of relaxation",
    "maximum number of relaxation steps reached",
    "exceeded the maximum number of relaxation",
)

#: 只读日志尾部的字节数；收敛标记就在最后几步内。
_TAIL_BYTES = 1 << 20


def _read_tail(path: Path, limit: int = _TAIL_BYTES) -> tuple[str, int, bool]:
    """读日志尾部（大文件只读最后 ``limit`` 字节，从行首对齐）。"""
    size = path.stat().st_size
    with path.open("rb") as fh:
        if size > limit:
            fh.seek(size - limit)
            chunk = fh.read()
            cut = chunk.find(b"\n")
            if cut >= 0:
                chunk = chunk[cut + 1:]
            return chunk.decode("utf-8", "replace"), size, True
        return fh.read().decode("utf-8", "replace"), size, False


def relax_status(aims_out: str | Path) -> dict | None:
    """返回 ``aims.out`` 的弛豫状态；文件不存在返回 ``None``。

    返回字段：``converged`` / ``terminated_normally`` / ``hit_step_limit`` /
    ``last_marker`` / ``steps`` / ``bytes`` / ``truncated`` / ``path``。
    ``steps`` 是尾部窗口里看到的几何更新次数（下界，日志被截断时只统计可见部分）。
    """
    path = Path(aims_out)
    if not path.is_file():
        return None
    text, size, truncated = _read_tail(path)

    last_marker = None
    steps = 0
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith(NOT_CONVERGED_MARKER):
            last_marker = "not_converged"
            steps += 1
        elif stripped.startswith(CONVERGED_MARKER):
            last_marker = "converged"
            steps += 1

    lowered = text.lower()
    hit_step_limit = any(hint in lowered for hint in _STEP_LIMIT_HINTS)
    terminated_normally = NORMAL_END_MARKER in text
    converged = (last_marker == "converged" and terminated_normally
                 and not hit_step_limit)
    return {
        "path": str(path),
        "converged": converged,
        "terminated_normally": terminated_normally,
        "hit_step_limit": hit_step_limit,
        "last_marker": last_marker,
        "steps": steps,
        "bytes": size,
        "truncated": truncated,
    }


def describe(status: dict | None) -> str:
    """把状态压成一句可放进报错信息的话。"""
    if status is None:
        return "未找到 aims.out，无法验收弛豫收敛"
    bits = ["最后标记=%s" % (status["last_marker"] or "无"),
            "尾部可见步数=%d" % status["steps"]]
    if not status["terminated_normally"]:
        bits.append("日志缺少 'Have a nice day.'（作业可能被截断或未正常结束）")
    if status["hit_step_limit"]:
        bits.append("命中弛豫步数上限")
    if status["truncated"]:
        bits.append("日志过大，只解析了尾部 %d 字节" % _TAIL_BYTES)
    return "，".join(bits)
