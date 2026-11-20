"""electrode_build 采集后自动登记连接核对（``review.auto_connection_check``）。

2026-11-15 用户要求「默认自动审核」：链式运行到 electrode_build collect 完成时，
直接在 build 目录写 ``connection_check.json``（who="auto"），不再需要人工
``smj review --approve --who``。本组测试锁住：

* 默认开启：collect 之后记录已落盘，可被 ``fingerprint.check_confirmation`` 验为有效；
* ``auto_connection_check: false``：不写任何记录（闸门仍在 review 侧 fail-closed）；
* 审计不过（表面索引指向非电极原子）时**不**写记录。
"""
from __future__ import annotations

import json
import math
from pathlib import Path

from smj.config import ReviewConfig, RunConfig, RunInfo
from smj.runner import StageRunner
from smj.junction.fingerprint import check_confirmation

TRI = ((0.0, 0.0), (2.88, 0.0), (1.44, 2.494))


def _geometry_text() -> str:
    atoms: list[tuple[str, float, float, float]] = []
    for k in range(6):
        ang = math.radians(60.0 * k)
        atoms.append(("C", 1.39 * math.cos(ang), 1.39 * math.sin(ang), 0.0))
    atoms.append(("S", 0.0, 0.0, -2.4))    # 7
    atoms.append(("S", 0.0, 0.0, 2.4))     # 8
    atoms.append(("Au", 0.0, 0.0, -4.8))   # 9
    atoms.append(("Au", 0.0, 0.0, 4.8))    # 10
    for z in (-7.15, -9.5):
        for dx, dy in TRI:
            atoms.append(("Au", dx, dy, z))  # 11-16
    for z in (7.15, 9.5):
        for dx, dy in TRI:
            atoms.append(("Au", dx, dy, z))  # 17-22
    return "".join("atom %.8f %.8f %.8f %s\n" % (x, y, z, e)
                   for e, x, y, z in atoms)


def _metadata(surface: tuple[int, int, int] = (14, 15, 16)) -> dict:
    return {
        "left_tip_index": 9,
        "right_tip_index": 10,
        "left_linker_index": 7,
        "right_linker_index": 8,
        "left_surface_indices": list(surface),
        "right_surface_indices": [20, 21, 22],
        "surface_rule": "det-v1",
    }


class _BuildStage:
    stage = "electrode_build"
    needs_submit = False
    outputs = ("geometry.in", "electrode_metadata.json")

    def __init__(self, surface: tuple[int, int, int] = (14, 15, 16)):
        self.surface = surface

    def collect(self, ctx, stage_dir):
        stage_dir = Path(stage_dir)
        (stage_dir / "geometry.in").write_text(_geometry_text(), encoding="utf-8")
        (stage_dir / "electrode_metadata.json").write_text(
            json.dumps(_metadata(self.surface), indent=2), encoding="utf-8")
        return {"artifacts": list(self.outputs), "results": {}, "notes": []}


def _runner(workdir: Path, *, auto: bool = True,
            surface: tuple[int, int, int] = (14, 15, 16)) -> StageRunner:
    config = RunConfig(
        run=RunInfo(name="run", dir=str(workdir), mode="dry-run"),
        review=ReviewConfig(auto_connection_check=auto),
    )
    stage = _BuildStage(surface)
    return StageRunner(config, dry_run=False, logger=lambda message: None,
                       stage_getter=lambda name: stage)


def _confirm_path(workdir: Path) -> Path:
    return workdir / "run" / "electrode_build" / "connection_check.json"


def test_collect_auto_registers_connection_check(workdir: Path):
    runner = _runner(workdir)
    runner.collect("electrode_build", force=True)

    assert _confirm_path(workdir).is_file()
    record = json.loads(_confirm_path(workdir).read_text(encoding="utf-8"))
    assert record["who"] == "auto"
    build = workdir / "run" / "electrode_build"
    assert check_confirmation(build)["valid"] is True
    # 审计文件同目录落盘，便于追溯
    assert (build / "connection_check_audit.json").is_file()


def test_collect_skips_when_auto_disabled(workdir: Path):
    runner = _runner(workdir, auto=False)
    runner.collect("electrode_build", force=True)
    assert not _confirm_path(workdir).exists()


def test_collect_refuses_when_audit_fails(workdir: Path):
    # 表面索引指向 C 原子（1-6）-> 审计必不过 -> 不得写记录
    runner = _runner(workdir, surface=(1, 2, 3))
    runner.collect("electrode_build", force=True)
    assert not _confirm_path(workdir).exists()
    assert (workdir / "run" / "electrode_build"
            / "connection_check_audit.json").is_file()
