"""2026-11-22 复审批次 2：编排状态机不能把"未知"当"失败"。

回归的真实缺陷：
1. ``sbatch`` 超时/解析不到 job id 时 runner 统一标 ``failed``，下一次 ``run``
   会 prepare + 重新提交，而作业可能已经在集群排队；
2. 作业已 ``finished`` 后 collect 解析失败同样被标 ``failed``，下一次 run 直接
   重提交，丢掉远端已完成的结果；
3. ``hash_stage_outputs`` 对缺件/超限文件静默省略，``outputs_sha256`` 看起来完整。
"""
from __future__ import annotations

from pathlib import Path

from smj.config import ClusterConfig, RunConfig, RunInfo
from smj.errors import StageError, SubmitUncertain
from smj.provenance import hash_stage_outputs
from smj.runner import StageOutcome, StageRunner


class _Stage:
    name = "molecule_opt"
    needs_submit = True
    script = "sub_g16.sh"
    outputs = ("mol.log",)

    def __init__(self) -> None:
        self.prepare_calls = 0
        self.collect_calls = 0

    def prepare(self, ctx, stage_dir):
        self.prepare_calls += 1
        return {"artifacts": ["mol.log"], "inputs": [], "notes": [],
                "script": self.script}

    def collect(self, ctx, stage_dir):
        self.collect_calls += 1
        return {"artifacts": ["mol.log"], "results": {}, "notes": []}


class _UncertainSession:
    def stage_dir(self, stage: str) -> Path:
        return Path("/remote/runs/t") / stage

    def submit(self, stage, stage_dir, *, script=None, sbatch_args=()):
        raise SubmitUncertain("sbatch 超时（120 s），无法确认作业是否已创建",
                              remote_dir=str(self.stage_dir(stage)))


def _runner(workdir: Path, stage: _Stage) -> StageRunner:
    config = RunConfig(
        run=RunInfo(name="t", dir=str(workdir / "run"), mode="slurm-ssh"),
        cluster=ClusterConfig(host="h", user="u", remote_root="/remote/runs",
                              poll_seconds=1.0, timeout_hours=0.5),
    )
    return StageRunner(config, transport=object(), logger=lambda message: None,
                       sleep=lambda seconds: None,
                       stage_getter=lambda name: stage)


def test_sbatch_uncertain_keeps_unknown_status(workdir):
    stage = _Stage()
    runner = _runner(workdir, stage)
    runner._session = lambda: _UncertainSession()

    outcome = runner.run_stage("molecule_opt")

    assert outcome.ok is False
    assert outcome.status == "submit_unknown"
    assert runner.state.status("molecule_opt") == "submit_unknown"
    assert "超时" in (runner.state.stage("molecule_opt").error or "")

    # 第二次 --run 必须停下来要求人工核对，而不是 prepare 后重提。
    again = runner.run()
    assert again[-1].status == "submit_unknown"
    assert stage.prepare_calls == 1


def test_collect_failure_after_finished_does_not_require_resubmit(workdir):
    stage = _Stage()
    runner = _runner(workdir, stage)
    runner.state.set("molecule_opt", status="waiting", job_id="123")
    runner.state.save(runner.layout.state_path)

    def _broken_collect(name, **kw):
        raise StageError("TE.dat 解析失败（临时 I/O 错误）")

    runner.collect = _broken_collect
    outcome = runner._collect_after_job(
        "molecule_opt", StageOutcome(stage="molecule_opt"))
    runner._save()

    assert outcome.ok is False
    assert outcome.status == "finished"
    state = runner.state.stage("molecule_opt")
    assert state.status == "finished"
    assert state.extra.get("collect_failed") is True
    assert "不重提作业" in (state.note or "")
    assert stage.prepare_calls == 0


def test_hash_stage_outputs_reports_unhashed(workdir):
    missing: list[str] = []
    hashes = hash_stage_outputs(workdir, ["absent.dat"], missing=missing)

    assert hashes == {}
    assert missing == ["absent.dat"]


def test_unknown_status_is_accepted_by_state(workdir):
    runner = _runner(workdir, _Stage())
    runner.state.set("molecule_opt", status="submit_unknown")
    assert runner.state.status("molecule_opt") == "submit_unknown"
