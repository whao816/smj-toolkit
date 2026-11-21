"""2026-11-20 复核 P1-4：SCF 产物"存在"不等于"可用"。

复核证据：``aims_scf.collect`` 原来只 ``is_file()``；零字节的
``basis-indices.out`` / ``omat.aims`` / ``mos.aims`` 加上"未收敛但正常退出"的
日志仍判 ``collected``，而且 ``results`` 里的 ``normal_termination`` /
``has_basis_indices`` / ``has_omat`` 是**写死的 True**，不是观测值。

这里锁住：

* 零字节产物 → 拒绝；
* 截断产物（只有 10 字节、没有 FHI-aims 头部）→ 拒绝；
* 日志正常结束但缺 ``Self-consistency cycle converged`` → 拒绝；
* 正常 demo 形状的产物 → ``collected``，且三个布尔字段来自观测。
"""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from _aims_fixtures import SCF_CONVERGED, write_scf_artifacts
from smj.errors import StageBlocked, StageError
from smj.stages.aims_scf import AimsScfStage

_GEOMETRY = "atom 0.0 0.0 0.0 S\natom 3.0 0.0 0.0 S\n"


class _Ctx:
    def __init__(self, run_dir: Path) -> None:
        self.layout = SimpleNamespace(run_dir=Path(run_dir))
        self.config = SimpleNamespace()

    def require_stage_file(self, stage: str, name: str) -> Path:
        path = self.layout.run_dir / stage / name
        if not path.is_file():
            raise StageBlocked(stage, "缺 %s" % path, [str(path)])
        return path

    def artifact(self, stage: str, name: str) -> Path:
        path = self.layout.run_dir / stage / name
        if not path.is_file():
            raise StageBlocked(stage, "缺 %s" % path, [str(path)])
        return path


def _seed(stage_dir: Path, *, converged: bool = True, spin: bool = False) -> None:
    stage_dir.mkdir(parents=True, exist_ok=True)
    (stage_dir / "geometry.in").write_text(_GEOMETRY, encoding="utf-8")
    log = ("  FHI-aims version      : 221103\n"
           + (SCF_CONVERGED if converged else "Have a nice day.\n"))
    (stage_dir / "aims.out").write_text(log, encoding="utf-8")
    write_scf_artifacts(stage_dir, spin=spin)


def _collect(workdir: Path, *, converged: bool = True):
    run_dir = workdir / "run"
    stage_dir = run_dir / "aims_scf"
    _seed(stage_dir, converged=converged)
    return AimsScfStage().collect(_Ctx(run_dir), stage_dir)


def test_scf_collect_records_observed_fields(workdir: Path):
    out = _collect(workdir)
    payload = out.results["aims_scf"]
    assert payload["normal_termination"] is True
    assert payload["scf_converged"] is True
    assert payload["has_basis_indices"] is True
    assert payload["has_omat"] is True
    assert payload["basis_functions"] == 4364
    assert payload["omat_nsaos"] == 4364
    stage_dir = workdir / "run" / "aims_scf"
    assert payload["basis_indices_bytes"] == (stage_dir / "basis-indices.out").stat().st_size
    assert payload["omat_bytes"] == (stage_dir / "omat.aims").stat().st_size
    assert payload["wavefunction_bytes"] == (stage_dir / "mos.aims").stat().st_size
    assert payload["log_bytes"] > 0


@pytest.mark.parametrize("name", ["basis-indices.out", "omat.aims", "mos.aims"])
def test_scf_collect_rejects_zero_byte_artifacts(workdir: Path, name: str):
    run_dir = workdir / "run"
    stage_dir = run_dir / "aims_scf"
    _seed(stage_dir)
    (stage_dir / name).write_bytes(b"")

    with pytest.raises(StageError, match="零字节"):
        AimsScfStage().collect(_Ctx(run_dir), stage_dir)


@pytest.mark.parametrize("name", ["basis-indices.out", "omat.aims", "mos.aims"])
def test_scf_collect_rejects_truncated_artifacts(workdir: Path, name: str):
    run_dir = workdir / "run"
    stage_dir = run_dir / "aims_scf"
    _seed(stage_dir)
    (stage_dir / name).write_bytes(b"1234567890")

    with pytest.raises(StageError, match="头部缺少"):
        AimsScfStage().collect(_Ctx(run_dir), stage_dir)


def test_scf_collect_rejects_unconverged_scf(workdir: Path):
    """日志以 'Have a nice day.' 正常结束，但没有 SCF 收敛标记。"""
    with pytest.raises(StageError, match="SCF 未收敛"):
        _collect(workdir, converged=False)


def test_scf_collect_rejects_alpha_only_without_beta(workdir: Path):
    run_dir = workdir / "run"
    stage_dir = run_dir / "aims_scf"
    _seed(stage_dir)
    (stage_dir / "mos.aims").unlink()
    (stage_dir / "alpha.aims").write_text(
        "$scfmo.aims   alpha spin channel\n", encoding="utf-8")

    with pytest.raises(StageError, match="alpha.aims \\+ beta.aims"):
        AimsScfStage().collect(_Ctx(run_dir), stage_dir)

def test_stage_and_report_agree_on_unconverged_log(workdir: Path):
    """两级判定必须一致（复核 R13a）。

    R5 之前：stage 照收、validate 报 ``aims_scf`` error；同一份日志两套结论。
    R5 之后：stage 抛 StageError，validate 也报 error——任何一边漏掉，这个用例都会红。
    """
    from smj.junction.validate import validate_junction_run

    run_dir = workdir / "run"
    stage_dir = run_dir / "aims_scf"
    _seed(stage_dir, converged=False)
    (stage_dir / "aims.out").write_text(
        "SCF is not yet converged\nHave a nice day.\n", encoding="utf-8")

    with pytest.raises(StageError, match="SCF 未收敛"):
        AimsScfStage().collect(_Ctx(run_dir), stage_dir)

    checks = [c for c in validate_junction_run(run_dir)["checks"] if c["name"] == "aims_scf"]
    assert len(checks) == 1 and checks[0]["passed"] is False
    assert checks[0]["severity"] == "error"


# --------------------------------------------------------------------------- #
# 2026-11-20 复审 P1-4（自旋通道）：两份波函数都要验，不能只看 is_file()
# --------------------------------------------------------------------------- #
def _spin_run(workdir: Path) -> tuple[object, Path]:
    run_dir = workdir / "run"
    stage_dir = run_dir / "aims_scf"
    _seed(stage_dir, spin=True)  # 只写 alpha + beta，没有 mos
    ctx = _Ctx(run_dir)
    ctx.config = SimpleNamespace(molecule=SimpleNamespace(multiplicity=2))
    return ctx, stage_dir


def test_spin_scf_rejects_empty_beta_aims(workdir: Path):
    ctx, stage_dir = _spin_run(workdir)
    (stage_dir / "beta.aims").write_bytes(b"")
    with pytest.raises(StageError, match=r"beta\.aims"):
        AimsScfStage().collect(ctx, stage_dir)


def test_spin_scf_rejects_truncated_beta_aims(workdir: Path):
    ctx, stage_dir = _spin_run(workdir)
    (stage_dir / "beta.aims").write_bytes(b"1234567890")
    with pytest.raises(StageError, match=r"beta\.aims"):
        AimsScfStage().collect(ctx, stage_dir)


def test_spin_scf_collect_records_both_channels(workdir: Path):
    ctx, stage_dir = _spin_run(workdir)
    payload = AimsScfStage().collect(ctx, stage_dir).results["aims_scf"]
    assert payload["has_alpha"] is True and payload["has_beta"] is True
    assert payload["has_mos"] is False and payload["spin_unpolarized"] is False
    assert payload["wavefunction_files"] == ["alpha.aims", "beta.aims"]
    assert payload["wavefunction_bytes"] == sum(
        (stage_dir / name).stat().st_size for name in ("alpha.aims", "beta.aims"))


def test_open_shell_without_spin_wavefunctions_is_rejected(workdir: Path):
    """multiplicity>1 但只有闭壳层 mos.aims：不得静默降级成闭壳层。"""
    run_dir = workdir / "run"
    stage_dir = run_dir / "aims_scf"
    _seed(stage_dir)
    ctx = _Ctx(run_dir)
    ctx.config = SimpleNamespace(molecule=SimpleNamespace(multiplicity=3))
    with pytest.raises(StageError, match=r"alpha\.aims"):
        AimsScfStage().collect(ctx, stage_dir)


def test_aitranss_te_rejects_empty_beta_wavefunction(workdir: Path):
    """下游 tcontrol 自己也要验：不能指望上游 collect 永远先跑过。"""
    from _aims_fixtures import ALPHA
    from smj.stages.aitranss_te import AitranssTeStage

    scf_dir = workdir / "run" / "aims_scf"
    scf_dir.mkdir(parents=True)
    (scf_dir / "alpha.aims").write_text(ALPHA, encoding="utf-8")
    (scf_dir / "beta.aims").write_bytes(b"")

    with pytest.raises(StageBlocked, match=r"beta\.aims"):
        AitranssTeStage()._spin_nsaos(scf_dir)


# --------------------------------------------------------------------------- #
# 2026-11-21 自由基支线 S1：Hirshfeld 自旋矩 -> 分子/电极分区
# --------------------------------------------------------------------------- #
_HIRSHFELD = """\
  Performing Hirshfeld analysis of fragment charges and moments.
  ----------------------------------------------------------------------
  | Atom     1: S
  |   Hirshfeld charge        :     -0.06377126
  |   Hirshfeld spin moment   :      0.25000000
  | Atom     2: S
  |   Hirshfeld charge        :     -0.04902760
  |   Hirshfeld spin moment   :      0.25000000
  | Atom     3: Au
  |   Hirshfeld charge        :      0.01000000
  |   Hirshfeld spin moment   :      0.40000000
  | Atom     4: Au
  |   Hirshfeld charge        :      0.01000000
  |   Hirshfeld spin moment   :      0.10000000
  ----------------------------------------------------------------------
"""


def test_collect_records_hirshfeld_spin_partition(workdir: Path):
    run_dir = workdir / "run"
    stage_dir = run_dir / "aims_scf"
    stage_dir.mkdir(parents=True)
    (stage_dir / "geometry.in").write_text(
        "atom  0.0 0.0 0.0 S\n"
        "atom  3.0 0.0 0.0 S\n"
        "atom  5.4 0.0 0.0 Au\n"
        "atom -2.4 0.0 0.0 Au\n", encoding="utf-8")
    (stage_dir / "aims.out").write_text(
        "  FHI-aims version      : 221103\n" + SCF_CONVERGED + _HIRSHFELD,
        encoding="utf-8")
    write_scf_artifacts(stage_dir)
    tip = run_dir / "aims_tip_relax"
    tip.mkdir()
    (tip / "input_summary.json").write_text(
        json.dumps({"atom_count": 4, "au_tip_count": 2}), encoding="utf-8")

    ctx = _Ctx(run_dir)
    ctx.config = SimpleNamespace(
        molecule=SimpleNamespace(multiplicity=1),
        transport=SimpleNamespace(extra_outputs=("hirshfeld",)),
    )
    out = AimsScfStage().collect(ctx, stage_dir)
    partition = out.results["aims_scf"]["spin_partition"]
    assert partition["total"] == pytest.approx(1.0)
    assert partition["groups"]["molecule"]["spin_moment"] == pytest.approx(0.5)
    assert partition["groups"]["electrode"]["spin_moment"] == pytest.approx(0.5)
    sidecar = stage_dir / "spin_partition.json"
    assert sidecar.is_file()
    stored = json.loads(sidecar.read_text(encoding="utf-8"))
    assert stored["atoms"] == 4
    assert "spin_partition.json" in out.artifacts