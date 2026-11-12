"""一次性：把已完成的 Gaussian 优化成果（mol.log）落成 gaussian_opt 的 collected 状态。

demo 说明：薁-2,6-二硫酚的 B3LYP/def2-TZVP opt+freq 已在主仓跑完（0 虚频），
本 demo 不重跑 Gaussian，只复用该日志；后续 6 个阶段全部真实执行。
"""
from pathlib import Path

from smj.config import load_config
from smj.context import StageContext
from smj.state import RunLayout, RunState
from smj.stages import get_stage

cfg = load_config(Path("smj.yaml"))
layout = RunLayout(cfg.run_dir)
state = RunState.load(layout.state_path, run=cfg.run.name)
ctx = StageContext(config=cfg, layout=layout, state=state)
stage = get_stage("gaussian_opt")
payload = stage.collect(ctx, layout.ensure_stage("gaussian_opt"))
state.set(
    "gaussian_opt",
    status="collected",
    artifacts=list(payload.artifacts),
    results=dict(payload.results),
    note="复用主仓 B3LYP/def2-TZVP opt+freq（0 虚频）日志",
)
state.save(layout.state_path)
print("collected:", payload.artifacts)
print("results:", payload.results)
