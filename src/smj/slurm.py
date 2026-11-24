"""Slurm 作业脚本的共用规则（全项目一份，不属于任何单条链）。

这里只放"渲染作业脚本时到处都要用、而且踩过坑"的纯函数，刻意不 import 本仓
其它模块，好让 ``junction.inputs``、``stages._thermal_common`` 这些位置都能安全
引用它（避免 ``smj.remote`` ↔ ``smj.junction`` 的循环导入）。

目前只有一条：GRES 行。背景见 :data:`DEFAULT_CPU_GRES`。
"""

from __future__ import annotations

#: 纯 CPU 作业默认显式申请的 GRES。
#:
#: 为什么要有这个默认值（2026-11-24 在 gpu1 上实测，属于**全项目**的调度经验）：
#: ``/etc/slurm/slurm.conf`` 里 ``NodeName=gpu1 ... Gres=gpu:4``，而分区 ``main``
#: 没有为 GRES 设默认值。Slurm 的默认行为是**作业按节点全部 GRES 分配**，于是
#: 任何不写 ``--gres`` 的作业都会把 4 张卡一起占下 —— 表现为「节点明明空着
#: （``CPUAlloc=8/CPUTot=64``），后面的作业却一直 PENDING(Resources)」，而
#: ``sbatch --test-only`` 给出的预计开始时间正好是前一个作业的时限终点。
#:
#: 判定实验：``sbatch --test-only -n 1 --mem=100M`` → 排队到 12 h 之后；
#: 加 ``--gres=gpu:0`` → 立刻可跑；``--gres=none`` **无效**（Slurm 不接受）。
#: 所以这里显式要 0 张卡，而不是省略这一行。
DEFAULT_CPU_GRES = "gpu:0"

#: 配置里用来表示"本集群没有 GRES，别渲染这一行"的取值。
_GRES_OFF = frozenset({"off", "none", "-", "no"})


def gres_line(configured: str | None = None, *, default: str = DEFAULT_CPU_GRES) -> str:
    """渲染 ``--gres=...`` 片段（不含 ``#SBATCH`` 前缀）。

    * ``configured`` 为空 → 用 ``default``（默认 ``gpu:0``，依据见
      :data:`DEFAULT_CPU_GRES`）；
    * ``configured`` 为 ``off`` / ``none`` / ``-`` / ``no`` → 返回空串，供没有
      GRES 的集群使用。注意 Slurm 不接受 ``--gres=none``，所以"关掉"只能是
      **不渲染**这一行，不能写成 ``none``；
    * 其余取值原样渲染。

    调用方负责把返回值拼进模板；返回空串时那一行应当整行消失（见
    ``sub_aims.sh`` / ``sub_thermal_*.sh`` 模板里 ``{gres_line}`` 的用法）。
    """
    value = "" if configured is None else str(configured).strip()
    if value.lower() in _GRES_OFF:
        return ""
    return "--gres=%s" % (value or default)


__all__ = ["DEFAULT_CPU_GRES", "gres_line"]
