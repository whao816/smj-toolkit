"""2026-11-20 复核 P1-3：弛豫未收敛不得被采集。

复核证据：``aims_tip_relax.collect`` 原来只判 ``geometry.in.next_step`` 存在，
从不打开 ``aims.out``；构造"步数上限 + ``Present geometry is not yet converged.``
+ ``Have a nice day.``"的日志仍被判 ``collected``，并给出 ``relaxed_au_s_A``，
``electrode_build`` 直接消费。``au_relax.collect`` 同样不解析收敛标记，却在
metadata 里硬写 ``converged: True``。

这里锁住四件事：

* 日志最后标记为未收敛 → 抛 ``StageError``（即使 ``next_step`` 在场）；
* 命中弛豫步数上限 → 同样拒绝；
* 缺 ``aims.out`` → 拒绝（"没日志"不等于"当作收敛"）；
* 日志收敛但缺 ``next_step`` → 拒绝（产物不完整）；
* 正向：真收敛的日志 + ``next_step`` → ``collected``；
* 下游双保险：``electrode_build.prepare`` 看到上游未收敛日志即 ``StageBlocked``。
"""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from smj.errors import StageBlocked, StageError
from smj.stages.aims_tip_relax import AimsTipRelaxStage
from smj.stages.electrode_build import ElectrodeBuildStage

_GEOMETRY = """\
atom 0.000000 0.000000 0.000000 S
atom 3.000000 0.000000 0.000000 S
atom 0.000000 0.000000 2.400000 Au
atom 3.000000 0.000000 2.400000 Au
"""

_CONVERGED_TAIL = """\
Maximum force component is  0.777213E-02 eV/A.
Present geometry is converged.
Self-consistency cycle converged.

Have a nice day.
"""

_UNCONVERGED_TAIL = """\
Present geometry is not yet converged.
Present geometry is not yet converged.
Present geometry is not yet converged.
Maximum force component is  0.910000E+00 eV/A.

Have a nice day.
"""

_STEP_LIMIT_TAIL = """\
Present geometry is not yet converged.
Reached the maximum number of relaxation steps (1000).

Have a nice day.
"""


class _Ctx:
    """aims_tip_relax.collect 需要的最小上下文。"""

    def __init__(self, run_dir: Path) -> None:
        self.layout = SimpleNamespace(run_dir=Path(run_dir))
        self.config = SimpleNamespace(
            electrode=SimpleNamespace(element="Au", left_anchor=None, right_anchor=None))

    def require_stage_file(self, stage: str, name: str) -> Path:
        path = self.layout.run_dir / stage / name
        if not path.is_file():
            raise StageBlocked(stage, "缺 %s" % path, [str(path)])
        return path


def _seed(stage_dir: Path, *, log_tail: str | None, next_step: bool = True) -> None:
    stage_dir.mkdir(parents=True, exist_ok=True)
    (stage_dir / "geometry.in").write_text(_GEOMETRY, encoding="utf-8")
    if next_step:
        (stage_dir / "geometry.in.next_step").write_text(_GEOMETRY, encoding="utf-8")
    if log_tail is not None:
        (stage_dir / "aims.out").write_text(
            "  FHI-aims version      : 221103\n"
            "  No maximum number of relaxation steps, defaulting to     1000\n"
            + log_tail,
            encoding="utf-8")
    (stage_dir / "input_summary.json").write_text(
        '{"left_anchor": 1, "right_anchor": 2, "tip_element": "Au",'
        ' "anchor_element": "S"}', encoding="utf-8")


def _collect(tmp_path: Path, *, log_tail: str | None, next_step: bool = True):
    run_dir = tmp_path / "run"
    stage_dir = run_dir / "aims_tip_relax"
    _seed(stage_dir, log_tail=log_tail, next_step=next_step)
    return AimsTipRelaxStage().collect(_Ctx(run_dir), stage_dir)


# --------------------------------------------------------------------------- #
# aims_tip_relax
# --------------------------------------------------------------------------- #
def test_tip_relax_rejects_unconverged_log(workdir: Path):
    """next_step 在场 + 日志最后标记未收敛：必须拒绝，不能给出"弛豫后"键长。"""
    with pytest.raises(StageError, match="未收敛"):
        _collect(workdir, log_tail=_UNCONVERGED_TAIL)


def test_tip_relax_rejects_step_limit(workdir: Path):
    with pytest.raises(StageError, match="步数上限"):
        _collect(workdir, log_tail=_STEP_LIMIT_TAIL)


def test_tip_relax_rejects_missing_log(workdir: Path):
    """没有 aims.out 时不能默认"当作收敛"：缺证据就是拒绝。"""
    with pytest.raises(StageError, match="缺少 aims.out"):
        _collect(workdir, log_tail=None)


def test_tip_relax_rejects_converged_without_next_step(workdir: Path):
    with pytest.raises(StageError, match="缺少 geometry.in.next_step"):
        _collect(workdir, log_tail=_CONVERGED_TAIL, next_step=False)


def test_tip_relax_accepts_converged(workdir: Path):
    out = _collect(workdir, log_tail=_CONVERGED_TAIL)
    payload = out.results["aims_tip_relax"]
    assert payload["relax_converged"] is True
    assert payload["relax_last_marker"] == "converged"
    assert payload["relax_steps"] >= 1
    assert payload["source"] == "geometry.in.next_step"
    assert set(payload["relaxed_au_s_A"]) == {"left", "right"}
    assert payload["relaxed_au_s_A"]["left"] == pytest.approx(2.4, abs=1e-6)


def test_tip_relax_default_step_limit_line_is_not_a_failure(workdir: Path):
    """"No maximum number of relaxation steps, defaulting to 1000" 不是命中上限。"""
    out = _collect(workdir, log_tail=_CONVERGED_TAIL)
    assert out.results["aims_tip_relax"]["relax_converged"] is True


# --------------------------------------------------------------------------- #
# electrode_build 双保险
# --------------------------------------------------------------------------- #
class _BuildCtx:
    def __init__(self, run_dir: Path) -> None:
        self.layout = SimpleNamespace(run_dir=Path(run_dir))
        self.config = SimpleNamespace(electrode=SimpleNamespace())

    def override_path(self, key: str):
        return None

    def artifact(self, stage: str, *parts: str) -> Path:
        return self.layout.run_dir / stage / parts[0]

    def upstream(self, stage: str, name: str, *, override=None) -> Path:
        path = self.artifact("aims_tip_relax", name)
        if not path.is_file():
            raise StageBlocked(stage, "缺 %s" % path, [str(path)])
        return path


def test_electrode_build_blocks_unconverged_upstream(workdir: Path):
    run_dir = workdir / "run"
    _seed(run_dir / "aims_tip_relax", log_tail=_UNCONVERGED_TAIL)
    with pytest.raises(StageBlocked, match="未收敛"):
        ElectrodeBuildStage().prepare(_BuildCtx(run_dir), run_dir / "electrode_build")
