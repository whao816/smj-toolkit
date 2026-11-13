"""回归：远端 prepare 必须有传输层；本地阶段 collect 不得回取。

两个 0.1.2 的现场 bug（薁 demo 首次真实跑通时暴露）：

1. ``StageRunner.prepare()`` 从不建传输层，``ctx.transport`` 一直是 ``None``，
   于是 ``resolve_cluster()`` 的"从集群 electrodes.library 取库团簇"分支永远
   走不到 —— 默认 ``electrode.source: library`` 在 slurm-ssh 下端直接失败，
   报错还是"未找到 electrodes.library"，把人指向错误的排查方向。
2. ``StageRunner.collect()`` 在 slurm-ssh 下对所有阶段都调用 ``download()``，
   而本地阶段（electrode_build / junction_analyze）不在 artifacts 台账里，于是
   ``smj stage electrode_build --collect`` 报"未知 smj 阶段"。

注：本机 %TEMP% 在 DSH 沙箱写限制之外，测试统一用 conftest 的 ``workdir``。
"""
from __future__ import annotations

from pathlib import Path

from smj.config import ClusterConfig, RunConfig, RunInfo
from smj.runner import StageRunner
from smj.testing import FakeTransport


def _config(workdir: Path) -> RunConfig:
    return RunConfig(
        run=RunInfo(name="t", dir=str(workdir / "run"), mode="slurm-ssh"),
        cluster=ClusterConfig(
            host="h", user="u", remote_root="/remote/runs",
            poll_seconds=1.0, timeout_hours=0.5),
    )


def _runner(workdir: Path, stage, transport=None) -> StageRunner:
    return StageRunner(
        _config(workdir), transport=transport, logger=lambda message: None,
        sleep=lambda seconds: None, stage_getter=lambda name: stage)


class _TransportSpyStage:
    """prepare 时记下 ctx.transport；collect 返回固定结果。"""

    name = "gaussian_opt"
    needs_submit = True
    script = "sub_g16.sh"
    outputs = ("mol.gjf",)

    def __init__(self) -> None:
        self.seen = {}

    def prepare(self, ctx, stage_dir):
        self.seen["transport"] = ctx.transport
        return {"artifacts": ["mol.gjf"], "inputs": [], "notes": [], "script": self.script}

    def collect(self, ctx, stage_dir):
        return {"artifacts": [], "results": {}, "notes": []}


class _LocalStage:
    """本地阶段：needs_submit=False，不应触发任何回取。"""

    name = "electrode_build"
    needs_submit = False
    script = None
    outputs = ("electrode_metadata.json",)

    def __init__(self) -> None:
        self.calls = 0

    def prepare(self, ctx, stage_dir):
        return {"artifacts": list(self.outputs), "inputs": [], "notes": []}

    def collect(self, ctx, stage_dir):
        self.calls += 1
        return {"artifacts": list(self.outputs), "results": {"ok": True}, "notes": []}


def test_prepare_builds_transport_for_slurm_ssh(workdir):
    stage = _TransportSpyStage()
    runner = _runner(workdir, stage)

    assert runner.context.transport is None, "前置：传输层是惰性的"
    runner.prepare("gaussian_opt")

    assert runner.transport is not None
    assert runner.context.transport is not None
    assert stage.seen["transport"] is runner.transport


def test_local_stage_collect_does_not_download(workdir):
    stage = _LocalStage()
    transport = FakeTransport()
    runner = _runner(workdir, stage, transport=transport)

    runner.prepare("electrode_build")
    runner.collect("electrode_build")

    assert stage.calls == 1
    downloads = [cmd for kind, cmd in transport.commands if kind in ("get", "run")]
    assert downloads == [], "本地阶段不应发起任何远端读写: %r" % (transport.commands,)
    assert runner.state.status("electrode_build") == "collected"
    assert runner.state.stage("electrode_build").results == {"ok": True}
