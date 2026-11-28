"""F3 回归：状态机一致性（wait 幂等 / submit 白名单 / prepare 清旧结果 / skipped）。

背景（2026-11-11 本地补丁 F3）：报告里点出 5 处状态机偏差，本文件把它们逐一钉住：

1. ``wait()`` 在 slurm-ssh 下先写 ``status="waiting"`` 再轮询，对一个已
   collected 的阶段跑 ``--wait`` 会把状态机回退（README 推荐的
   ``--prepare --submit --wait`` 组合正好能触发）；
2. ``submit()`` 靠"没有分支拦住"让 failed 可提交，语义是偶然的；
3. ``prepare()`` 把 status 改回 prepared，却保留上一轮的 results 与
   ``download_*`` 台账，status/report 会把旧数值当成本轮产物；
4. ``skipped`` 在 ``is_done()`` 里算完成，``--run`` 却会静默重跑；
5. 脚本里仍是 ``UNCONFIGURED`` 占位时 sbatch 照提交，作业在远端才炸。

注：本机 %TEMP% 在 DSH 沙箱写限制之外，故用 conftest 里的 ``workdir``
（仓库内 ``.tmp/``）而不是 pytest 的 ``tmp_path``。
"""
from __future__ import annotations

from pathlib import Path

import pytest

from smj.config import ClusterConfig, RunConfig, RunInfo
from smj.errors import StageError
from smj.runner import StageRunner
from smj.state import DONE_STATUSES, RunState
from smj.testing import FakeTransport
from smj.transport import CommandResult


def _runner(workdir: Path, transport: FakeTransport | None = None,
            stage_getter=None) -> StageRunner:
    config = RunConfig(
        run=RunInfo(name="t", dir=str(workdir / "run"), mode="slurm-ssh"),
        cluster=ClusterConfig(
            host="h", user="u", remote_root="/home/u/smj-runs",
            poll_seconds=1.0, timeout_hours=0.5),
    )
    return StageRunner(
        config, transport=transport or FakeTransport(),
        logger=lambda message: None, sleep=lambda seconds: None,
        stage_getter=stage_getter)


class _FakeStage:
    """最小 Stage：prepare 落一个脚本，collect 返回一条 results。"""

    stage = "gaussian_opt"
    script = "sub_g16.sh"
    outputs = ("out.txt",)
    needs_submit = True

    def prepare(self, ctx, stage_dir):
        Path(stage_dir, "sub_g16.sh").write_text("#!/bin/sh\necho hi\n", encoding="utf-8")
        return {"artifacts": ["out.txt"], "notes": [], "inputs": []}

    def collect(self, ctx, stage_dir):
        return {"artifacts": ["out.txt"], "results": {"x": 1}, "notes": []}


class _SbatchTransport(FakeTransport):
    """让 sbatch 返回一个可解析的 job id。"""

    def on_run(self, command, cwd, timeout):
        if "sbatch" in command:
            return CommandResult(0, "Submitted batch job 424242\n", "", command=command)
        return None


def test_done_statuses_single_source():
    assert DONE_STATUSES == ("finished", "collected", "skipped")


def test_wait_refuses_done_status(workdir):
    runner = _runner(workdir)
    runner.state.set("gaussian_opt", status="collected", job_id="123")
    runner._save()

    with pytest.raises(StageError, match="已结束"):
        runner.wait("gaussian_opt")

    # 关键：状态不能被 --wait 打回 waiting
    assert runner.state.status("gaussian_opt") == "collected"


def test_submit_refuses_unconfigured_placeholder(workdir):
    runner = _runner(workdir)
    stage_dir = runner.layout.ensure_stage("gaussian_opt")
    (stage_dir / "sub_g16.sh").write_text(
        "#!/bin/bash\nG16_BIN=UNCONFIGURED\n", encoding="utf-8")
    (stage_dir / "mol.gjf").write_text("#p hf/sto-3g\n\n0 1\nH 0 0 0\n", encoding="utf-8")
    runner.state.set("gaussian_opt", status="prepared")
    runner._save()

    with pytest.raises(StageError, match="UNCONFIGURED"):
        runner.submit("gaussian_opt")

    assert runner.state.status("gaussian_opt") == "prepared"
    assert not any("sbatch" in command for _, command in runner.transport.commands)


def test_submit_allows_failed_for_retry(workdir):
    runner = _runner(workdir, transport=_SbatchTransport())
    stage_dir = runner.layout.ensure_stage("gaussian_opt")
    from smj.junction.inputs import SLURM_G16_SCRATCH
    (stage_dir / "sub_g16.sh").write_text(
        "#!/bin/sh\n" + SLURM_G16_SCRATCH + "\nsrun g16 \"$1\"\n", encoding="utf-8")
    (stage_dir / "mol.gjf").write_text("#p hf/sto-3g\n\n0 1\nH 0 0 0\n", encoding="utf-8")
    runner.state.set("gaussian_opt", status="failed", error="上一次 sbatch 失败")
    runner._save()

    submission = runner.submit("gaussian_opt")

    assert submission["job_id"] == "424242"
    assert runner.state.status("gaussian_opt") == "submitted"


@pytest.mark.parametrize("scratch", ["/data/tmp/g16", "/tmp/g16"])
def test_submit_rejects_prepared_external_gaussian_scratch(workdir, scratch):
    runner = _runner(workdir, transport=_SbatchTransport())
    stage_dir = runner.layout.ensure_stage("gaussian_opt")
    (stage_dir / "sub_g16.sh").write_text(
        f'#!/bin/bash\nmkdir -p "{scratch}"\nexport GAUSS_SCRDIR={scratch}\nsrun g16 "$1"\n',
        encoding="utf-8")
    runner.state.set("gaussian_opt", status="prepared")
    with pytest.raises(StageError, match="scratch"):
        runner.submit("gaussian_opt")
    assert runner.transport.commands == []


def test_submit_rejects_gaussian_scratch_override(workdir):
    from smj.junction.inputs import SLURM_G16_SCRATCH
    runner = _runner(workdir, transport=_SbatchTransport())
    stage_dir = runner.layout.ensure_stage("gaussian_opt")
    (stage_dir / "sub_g16.sh").write_text(
        SLURM_G16_SCRATCH + '\nexport GAUSS_SCRDIR=/tmp\nsrun g16 "$1"\n', encoding="utf-8")
    runner.state.set("gaussian_opt", status="prepared")
    with pytest.raises(StageError, match="scratch"):
        runner.submit("gaussian_opt")
    assert runner.transport.commands == []


def test_submit_still_refuses_collected(workdir):
    runner = _runner(workdir)
    runner.state.set("gaussian_opt", status="collected")
    runner._save()
    with pytest.raises(StageError, match="已完成"):
        runner.submit("gaussian_opt")


def test_prepare_clears_stale_results(workdir):
    runner = _runner(workdir, stage_getter=lambda name: _FakeStage())
    runner.state.set(
        "gaussian_opt", status="collected", results={"old": 123},
        download_fetched=["out.txt"], download_missing=["TE.dat"])
    runner._save()

    runner.prepare("gaussian_opt", force=True)

    stage = runner.state.stage("gaussian_opt")
    assert stage.results == {}
    assert runner.state.status("gaussian_opt") == "prepared"
    for stale in ("download_fetched", "download_missing"):
        assert stale not in stage.extra


def test_run_skips_skipped_stage(workdir):
    runner = _runner(workdir, stage_getter=lambda name: _FakeStage())
    runner.state.set("gaussian_opt", status="skipped")
    runner._save()

    outcomes = runner.run("gaussian_opt", "gaussian_opt")

    assert len(outcomes) == 1
    assert outcomes[0].status == "skipped"
    assert outcomes[0].ok is True
    assert "skipped" in outcomes[0].message


def test_skipped_is_done_and_blocks_nothing(workdir):
    runner = _runner(workdir)
    runner.state.set("electrode_build", status="skipped")
    assert runner.state.is_done("electrode_build")
    assert runner.state.first_incomplete() == "molecule_opt"
