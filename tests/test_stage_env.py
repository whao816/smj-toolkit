"""阶段环境变量：GPU 请求只给真正用 GPU 的阶段。

2026-11-13 的验收 run 里 aitranss_te 只排到 CPU，但 sub_te.sh 继承了
cluster.gres=gpu:1，集群 4 张卡被别的作业占满时它一直 PENDING(Resources)；
aitranss 是 OpenMP/CPU 程序，请求 GPU 只会无谓排队。
"""

from __future__ import annotations

from types import SimpleNamespace

from smj.config import ClusterConfig
from smj.context import StageContext


def _context(gres: str | None, env: dict | None = None) -> StageContext:
    config = SimpleNamespace(cluster=ClusterConfig(gres=gres, env=dict(env or {})))
    return SimpleNamespace(config=config)


def test_aitranss_does_not_inherit_gpu_request():
    ctx = _context("gpu:1")
    env = StageContext.stage_env(ctx, "aitranss_te", {"cpus": 4})
    assert "te_gres" not in env


def test_aitranss_does_not_inherit_bulk_memory_request():
    ctx = _context("gpu:1")
    ctx.config.cluster.mem = "24G"
    env = StageContext.stage_env(ctx, "aitranss_te", {"cpus": 4, "mem": "8G"})
    assert env["te_mem"] == "8G"


def test_aitranss_keeps_explicit_te_gres():
    ctx = _context("gpu:1", {"te_gres": "gpu:0"})
    env = StageContext.stage_env(ctx, "aitranss_te", {"cpus": 4})
    assert env["te_gres"] == "gpu:0"


def test_gaussian_still_inherits_gpu_request():
    ctx = _context("gpu:1")
    env = StageContext.stage_env(ctx, "gaussian_opt", {})
    assert env["g16_gres"] == "gpu:1"
