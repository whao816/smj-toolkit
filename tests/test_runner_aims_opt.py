"""回归：AIMS 引擎的 gaussian_opt 回取文件名。

vendored ``retrieval_plan('gaussian_opt')`` 是 Gaussian 口径（mol.gjf/mol.log/
sub_g16.sh）；当 ``molecule.opt_engine: aims`` 时这些文件根本不存在。runner 必须
显式传入 AIMS 产物名，否则 collect 会在回取阶段错误地拒绝。
"""
from __future__ import annotations

from pathlib import Path

from smj.config import ClusterConfig, MoleculeConfig, RunConfig, RunInfo
from smj.runner import StageRunner
from smj.testing import FakeTransport

REMOTE = "/remote/runs/t/molecule_opt"


def _config(workdir: Path) -> RunConfig:
    xyz = workdir / "molecule.xyz"
    xyz.write_text("2\n\nC 0 0 0\nH 0 0 1\n", encoding="utf-8")
    return RunConfig(
        run=RunInfo(name="t", dir=str(workdir / "run"), mode="slurm-ssh"),
        molecule=MoleculeConfig(xyz=str(xyz), opt_engine="aims"),
        cluster=ClusterConfig(
            host="h", user="u", remote_root="/remote/runs",
            poll_seconds=1.0, timeout_hours=0.5),
    )


class _AimsStage:
    name = "molecule_opt"
    needs_submit = True
    script = "sub_aims.sh"
    outputs = ("opt.xyz", "aims.out")

    def prepare(self, ctx, stage_dir):
        return {
            "artifacts": ["geometry.in", "control.in", "sub_aims.sh", "input_summary.json"],
            "inputs": [],
            "notes": [],
            "script": "sub_aims.sh",
        }

    def collect(self, ctx, stage_dir):
        assert (stage_dir / "aims.out").exists()
        assert (stage_dir / "geometry.in.next_step").exists()
        return {"artifacts": ["opt.xyz", "aims.out"], "results": {"engine": "aims"}, "notes": []}


def test_aims_gaussian_opt_collect_downloads_aims_names(workdir):
    transport = FakeTransport()
    for name in ("geometry.in", "control.in", "sub_aims.sh", "input_summary.json",
                 "aims.out", "geometry.in.next_step", ".dsh_done"):
        transport.add_text(REMOTE + "/" + name, "x\n")
    stage = _AimsStage()
    runner = StageRunner(
        _config(workdir), transport=transport, logger=lambda message: None,
        sleep=lambda seconds: None, stage_getter=lambda name: stage)

    runner.prepare("gaussian_opt")
    runner.collect("gaussian_opt")

    gets = [command for kind, command in transport.commands if kind == "get"]
    assert any("aims.out" in command for command in gets)
    assert any("geometry.in.next_step" in command for command in gets)
    assert not any("mol.log" in command for command in gets)
    assert runner.state.status("molecule_opt") == "collected"
