"""S3：开壳层配单列 TE.dat 必须显式放行（2026-11-21 自由基支线改进计划）。

复核依据：开壳层体系的 ``TE.dat`` 只有 1 列时，α/β 通道在写入时已合成均值，
后续 ``analyze_te_spin`` 只能给出 P≡0；旧实现只把这件事记成 review warning，
阶段照样 ``collected``，下游无从分辨"真的没极化"与"自旋信息在解析前就丢了"。

这里锁住：
* 开壳层 + 单列 → ``StageError``（默认阻断）；
* ``transport.allow_single_channel: true`` → ``collected`` 且结果里留痕；
* 闭壳层 + 单列 → 照常 ``collected``；
* 开壳层 + 双列 → 照常 ``collected``（channels=2）。
"""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from smj.errors import StageError
from smj.stages.aitranss_te import AitranssTeStage

_4COL = """\
#bias = 0.0000 V efermi = -0.18 H
-0.1900000000 -0.20 0.2000 0.1900
-0.1800000000 0.00 0.2500 0.2400
-0.1700000000 0.20 0.2000 0.1900
"""
_3COL = """\
#bias = 0.0000 V efermi = -0.18 H
-0.1900000000 -0.20 0.1950
-0.1800000000 0.00 0.2450
-0.1700000000 0.20 0.1950
"""


class _Ctx:
    def __init__(self, run_dir: Path, *, multiplicity: int = 2,
                 allow_single_channel: bool = False) -> None:
        self.layout = SimpleNamespace(run_dir=Path(run_dir))
        self.config = SimpleNamespace(
            molecule=SimpleNamespace(multiplicity=multiplicity, charge=0,
                                     fixed_spin_moment=True,
                                     broken_symmetry=False),
            transport=SimpleNamespace(require_te_out=False,
                                      allow_single_channel=allow_single_channel),
        )

    def require_stage_file(self, stage: str, name: str) -> Path:
        path = self.layout.run_dir / stage / name
        if not path.is_file():
            raise StageError("缺 %s" % path)
        return path


def _run(workdir: Path, text: str) -> Path:
    run_dir = workdir / "run"
    stage_dir = run_dir / "aitranss_te"
    stage_dir.mkdir(parents=True)
    (stage_dir / "TE.dat").write_text(text, encoding="utf-8")
    return run_dir


def test_open_shell_with_single_column_te_is_rejected(workdir: Path):
    run_dir = _run(workdir, _3COL)
    with pytest.raises(StageError, match="自旋通道"):
        AitranssTeStage().collect(_Ctx(run_dir), run_dir / "aitranss_te")


def test_single_channel_is_allowed_only_when_configured(workdir: Path):
    run_dir = _run(workdir, _3COL)
    ctx = _Ctx(run_dir, allow_single_channel=True)
    out = AitranssTeStage().collect(ctx, run_dir / "aitranss_te")
    payload = out.results["aitranss_te"]
    assert payload["transmission_channels"] == 1
    assert payload["single_channel_allowed"] is True
    assert any("allow_single_channel" in note for note in out.notes)


def test_closed_shell_with_single_column_te_is_fine(workdir: Path):
    run_dir = _run(workdir, _3COL)
    out = AitranssTeStage().collect(_Ctx(run_dir, multiplicity=1), run_dir / "aitranss_te")
    payload = out.results["aitranss_te"]
    assert payload["transmission_channels"] == 1
    assert payload["single_channel_allowed"] is False


def test_open_shell_with_two_channels_is_fine(workdir: Path):
    run_dir = _run(workdir, _4COL)
    out = AitranssTeStage().collect(_Ctx(run_dir), run_dir / "aitranss_te")
    payload = out.results["aitranss_te"]
    assert payload["transmission_channels"] == 2
    assert payload["single_channel_allowed"] is False
