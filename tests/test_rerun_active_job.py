"""2026-11-20 复核 P1-5：正在远端运行的阶段不能被"重跑"打成 failed。

复核证据：``state=waiting`` + ``job_id=12345`` 时调 ``run_stage``，``prepare``
抛"正在运行"的保护异常，被 ``run_stage`` 的兜底 ``except SmjError`` 写成
``status=failed``；下一次 ``run`` 于是重新 prepare + sbatch，在集群上留下
第二个作业（实测新作业 424242 覆盖了 12345 的登记）。

这里锁住四条行为：

1. ``run_stage`` 对活跃作业返回 ``ok=False`` 但**保留** ``waiting``/``job_id``；
2. 紧接着 ``run`` 走"等待 + 补 collect"，**不重新提交**；
3. ``--force`` 仍然按显式语义重跑（覆盖保护）；
4. 等待超时保持 ``waiting``，不写 failed。
"""
from __future__ import annotations

from pathlib import Path

from smj.config import ClusterConfig, RunConfig, RunInfo
from smj.runner import StageRunner
from smj.testing import FakeTransport
from smj.transport import CommandResult

REMOTE_ROOT = "/remote/rerun"


class _FakeStage:
    """最小 aims_scf 阶段：prepare 落脚本，collect 返回一条结果。"""

    stage = "aims_scf"
    script = "sub_aims.sh"
    outputs = ("aims.out",)
    needs_submit = True

    def __init__(self) -> None:
        self.prepare_calls = 0
        self.collect_calls = 0

    def prepare(self, ctx, stage_dir):
        self.prepare_calls += 1
        Path(stage_dir, "sub_aims.sh").write_text(
            "#!/bin/sh\necho aims\n", encoding="utf-8")
        Path(stage_dir, "control.in").write_text(
            "# fake control\n", encoding="utf-8")
        return {"artifacts": ["aims.out"], "inputs": []}

    def collect(self, ctx, stage_dir):
        self.collect_calls += 1
        return {"artifacts": ["aims.out"],
                "results": {"band_gap_eV": 1.234}, "notes": []}


class _SbatchTransport(FakeTransport):
    """sbatch 返回 job id；轮询按场景返回"已结束"或"仍在队列"。"""

    def __init__(self, *, finish: bool = True) -> None:
        super().__init__()
        self.finish = finish
        self.sbatch_calls = 0
        self.probe_calls = 0

    def on_run(self, command, cwd, timeout):
        if "sbatch" in command:
            self.sbatch_calls += 1
            return CommandResult(0, "Submitted batch job 424242\n", "",
                                 command=command)
        if "@@SENT" in command:
            self.probe_calls += 1
            if self.finish:
                return CommandResult(
                    0,
                    "@@SENT 0 0 2026-11-20T00:00:00+08:00\n"
                    "@@QRC 0\n@@ID 424242\n@@END\n",
                    "", command=command)
            return CommandResult(
                0, "@@SENT 0 \n@@QRC 0\n@@ID 12345\n@@END\n", "",
                command=command)
        return None


def _config(run_parent: Path, name: str) -> RunConfig:
    return RunConfig(
        run=RunInfo(name=name, dir=str(run_parent), mode="slurm-ssh"),
        cluster=ClusterConfig(host="fake-host", user="u",
                              remote_root=REMOTE_ROOT, partition="p",
                              poll_seconds=0.01, timeout_hours=0.5),
    )


def _runner(workdir: Path, stage: _FakeStage, transport, *,
            name: str = "rerun", timeout_s: float | None = None) -> StageRunner:
    parent = workdir / "runs"
    parent.mkdir(parents=True, exist_ok=True)
    return StageRunner(_config(parent, name), transport=transport,
                       logger=lambda _m: None, sleep=lambda _s: None,
                       stage_getter=lambda _key: stage, timeout_s=timeout_s)


def _mark_waiting(runner: StageRunner, *, job_id: str = "12345") -> None:
    runner.state.set(
        "aims_scf", status="waiting", job_id=job_id,
        remote_dir="%s/%s/aims_scf" % (REMOTE_ROOT, runner.config.run.name),
        prepared_at="2026-11-20T00:00:00+08:00",
        submitted_at="2026-11-20T00:01:00+08:00",
        note="作业运行中")
    runner._save()


def _seed_remote(transport: _SbatchTransport, name: str = "rerun") -> None:
    """远端产物齐全：resume 分支只差"回取"，不该在缺件上失败。"""
    root = "%s/%s/aims_scf" % (REMOTE_ROOT, name)
    for rel, text in {
        "geometry.in": "lattice_vector 0 0 0\n",
        "control.in": "# fake control\n",
        "aims.out": ("Fermi level (Ha): -0.1884156757829439\n"
                     "Self-consistency cycle converged.\n"
                     "Have a nice day.\n"),
        "basis-indices.out": ("$basis.indices.aims\n"
                              "# amount of basis functions used:   4364\n"),
        "sub_aims.sh": "#!/bin/sh\necho aims\n",
        ".dsh_done": "0 2026-11-20T00:00:00+08:00\n",
        "mos.aims": "$scfmo.aims\nnsaos = 4364\n",
        "omat.aims": "$omat.aims   format(4d20.14)   nsaos=  4364\n",
    }.items():
        transport.add_text("%s/%s" % (root, rel), text)


def test_run_stage_keeps_waiting_state(workdir: Path):
    stage = _FakeStage()
    transport = _SbatchTransport(finish=True)
    runner = _runner(workdir, stage, transport)
    _mark_waiting(runner)
    _seed_remote(transport)

    outcome = runner.run_stage("aims_scf")

    assert outcome.ok is False
    assert outcome.status == "waiting"
    assert "正在远端运行" in outcome.message
    assert runner.state.status("aims_scf") == "waiting"
    assert runner.state.stage("aims_scf").job_id == "12345"
    assert stage.prepare_calls == 0
    assert transport.sbatch_calls == 0


def test_run_resumes_active_job_without_resubmit(workdir: Path):
    stage = _FakeStage()
    transport = _SbatchTransport(finish=True)
    runner = _runner(workdir, stage, transport)
    _mark_waiting(runner)
    _seed_remote(transport)

    outcomes = runner.run(start="aims_scf", end="aims_scf")

    assert len(outcomes) == 1
    assert outcomes[0].ok is True
    assert outcomes[0].status == "collected"
    assert outcomes[0].action == "resume"
    assert stage.collect_calls == 1
    assert stage.prepare_calls == 0
    assert transport.sbatch_calls == 0
    assert runner.state.status("aims_scf") == "collected"


def test_force_still_reruns_active_job(workdir: Path):
    stage = _FakeStage()
    transport = _SbatchTransport(finish=True)
    runner = _runner(workdir, stage, transport)
    _mark_waiting(runner)
    _seed_remote(transport)

    outcomes = runner.run(start="aims_scf", end="aims_scf", force=True)

    assert outcomes[0].ok is True
    assert stage.prepare_calls == 1
    assert transport.sbatch_calls == 1


def test_wait_timeout_keeps_waiting(workdir: Path):
    stage = _FakeStage()
    transport = _SbatchTransport(finish=False)
    runner = _runner(workdir, stage, transport, timeout_s=0.002)
    _mark_waiting(runner)
    _seed_remote(transport)

    outcomes = runner.run(start="aims_scf", end="aims_scf")

    assert outcomes[0].ok is False
    assert outcomes[0].status == "waiting", outcomes[0].to_dict()
    assert runner.state.status("aims_scf") == "waiting"
    assert runner.state.stage("aims_scf").job_id == "12345"
    assert transport.sbatch_calls == 0
