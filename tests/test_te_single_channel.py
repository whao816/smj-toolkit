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
from smj.analysis import analyze_te, te_header
from smj.junction.independent_check import parse_te
from smj.junction.validate import validate_junction_run
from smj.landauer_spin import analyze_te_spin
from smj.plotting import read_curve
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


@pytest.mark.parametrize("row", ["0 0 0.25", "0 0 0.25 ***", "0 0", "0 0 0.25 nan",
                                  "0 0 0.25 0.24 0.23"])
@pytest.mark.parametrize("reader", [analyze_te, analyze_te_spin, read_curve, parse_te, te_header])
def test_analysis_rejects_incomplete_or_invalid_spin_row(workdir, reader, row):
    te = workdir / "TE.dat"
    te.write_text(_4COL + row + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="第 5 行"):
        reader(te)


@pytest.mark.parametrize("allow_single", [False, True])
@pytest.mark.parametrize("spectrum", [
    _4COL + "0 0 0.25\n",
    _3COL + "0 0 0.25 0.24\n",
    "# E [Hartree] E-EF [eV] T(E).alpha T(E).beta\n" + _3COL,
    "#spin-polarized calculation\n" + _3COL,
], ids=["missing-beta", "added-beta", "declared-beta-column", "declared-spin-polarized"])
def test_collect_rejects_inconsistent_channels(workdir, allow_single, spectrum):
    run = _run(workdir, spectrum)
    ctx = _Ctx(run, allow_single_channel=allow_single)
    with pytest.raises(StageError, match="TE.dat.*第"):
        AitranssTeStage().collect(ctx, run / "aitranss_te")


def test_spin_reader_preserves_fortran_numbers_and_comments(workdir):
    te = workdir / "TE.dat"
    te.write_text(" # alpha and beta\n\n0 0 3D-1 2D-1 # first\n1 1 0.3 0.2\n",
                  encoding="utf-8")
    assert read_curve(te)["t"] == pytest.approx([0.25, 0.25])
    assert [row[2] for row in parse_te(te)] == pytest.approx([0.25, 0.25])


def test_validate_marks_corrupt_spectrum_unusable(workdir):
    run = _run(workdir, _4COL + "0 0 0.25\n")
    verdict = validate_junction_run(run)
    stage = next(check for check in verdict["checks"] if check["name"] == "aitranss_te")
    assert not stage["passed"]
    assert "第 5 行" in stage["detail"]
    assert verdict["physics_verdict"] == "inconclusive"
