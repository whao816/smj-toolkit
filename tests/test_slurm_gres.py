"""GRES 行的共用规则（``smj.slurm``）与它在 FHI-aims 作业脚本里的落地。

为什么要专门测这个：gpu1 的节点声明了 ``Gres=gpu:4`` 而分区没给默认值，
**不写 ``--gres`` 的作业会把整节点 4 张卡一起占下**，于是「节点明明空着、
后面的作业却一直排队」。这个坑吃过一次（324 位移数组排到 relax 时限终点之后），
所以把规则收进 ``smj.slurm`` 供全项目共用，并用测试钉住三条语义：

* 不配置 → 显式 ``gpu:0``（纯 CPU 作业不该占卡）；
* 配了值 → 原样渲染；
* ``off``/``none``/``-``/``no`` → **整行消失**（Slurm 不接受 ``--gres=none``，
  所以"关掉"只能是不渲染）。
"""

from __future__ import annotations

import pytest

from smj.junction.inputs import _render_aims_sub
from smj.slurm import DEFAULT_CPU_GRES, gres_line


def test_default_requests_zero_gpus():
    """没配置时必须显式要 0 张卡 —— 省略这一行等于要整节点。"""
    assert DEFAULT_CPU_GRES == "gpu:0"
    assert gres_line() == "--gres=gpu:0"
    assert gres_line(None) == "--gres=gpu:0"
    assert gres_line("") == "--gres=gpu:0"
    assert gres_line("   ") == "--gres=gpu:0"


def test_configured_value_is_rendered_verbatim():
    assert gres_line("gpu:1") == "--gres=gpu:1"
    assert gres_line(" gpu:2 ") == "--gres=gpu:2"
    assert gres_line("mps:1") == "--gres=mps:1"


@pytest.mark.parametrize("value", ["off", "OFF", "none", "-", "no"])
def test_off_values_drop_the_line_entirely(value):
    """Slurm 不接受 ``--gres=none``，所以"关掉"只能是返回空串、不渲染整行。"""
    assert gres_line(value) == ""


def test_custom_default_is_honoured():
    assert gres_line(None, default="gpu:4") == "--gres=gpu:4"


# --- 落到 sub_aims.sh：纯 CPU 分支必须出现 gpu:0 -------------------------

def _env(**over):
    base = {"partition": "main", "aims_bin": "/opt/aims/bin/aims.x"}
    base.update(over)
    return base


def test_aims_cpu_script_asks_for_zero_gpus():
    script = _render_aims_sub(_env(), ntasks=8)
    assert "#SBATCH --gres=gpu:0" in script
    # 纯 CPU：不该出现 GPU 版的 --gpus=N
    assert "--gpus=" not in script


def test_aims_gpu_mode_keeps_gpus_flag_not_gres():
    """aims_gpus 非空是真正的 GPU 模式（1 GPU ↔ 1 rank），走 --gpus=N。"""
    script = _render_aims_sub(_env(aims_gpus="4"), ntasks=8)
    assert "#SBATCH --gpus=4" in script
    assert "--gres=" not in script
    # rank 数跟随 GPU 数
    assert "mpirun -np 4" in script


def test_aims_sched_gres_override_is_rendered():
    """aims_gres_sched 只影响调度、不改 rank 数。"""
    script = _render_aims_sub(_env(aims_gres_sched="gpu:1"), ntasks=8)
    assert "#SBATCH --gres=gpu:1" in script
    assert "mpirun -np 8" in script


def test_aims_sched_gres_off_removes_the_line():
    """没有 GRES 的集群用 off 关掉整行，而不是渲染成 --gres=none。"""
    script = _render_aims_sub(_env(aims_gres_sched="off"), ntasks=8)
    assert "--gres" not in script
    assert "mpirun -np 8" in script
