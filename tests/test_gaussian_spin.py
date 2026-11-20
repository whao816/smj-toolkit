"""2026-11-20 规划 §6 P1-6：Gaussian 日志里的 ⟨S²⟩ 与 SCF 收敛证据。

开壳层单分子结用 Gaussian 预优化时，判"这个解是不是纯自旋态"靠 ``mol.log``：

* ``S**2 before annihilation X, after Y`` —— 湮灭前后的自旋平方期望值；
* ``SCF Done:`` 是否出现、有无 ``Convergence failure`` —— SCF 是否真的收敛；
* 可选 ``Stable=Opt`` —— 先做波函数稳定性检查再优化。

这里锁住解析与展示口径：多重度 M 的理想值 S(S+1)=((M−1)/2)((M+1)/2)，
偏差 > 0.1 要在 notes 里点名「自旋污染」。
"""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from smj.config import MoleculeConfig
from smj.engines.gaussian import GaussianEngine
from smj.errors import StageError
from smj.junction.inputs import make_gaussian_input
from smj.junction.spectra import parse_gaussian_log

_ORIENTATION = """\
                          Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          6            0        0.000000    0.000000    0.000000
      2          1            0        1.090000    0.000000    0.000000
 ---------------------------------------------------------------------
"""


def _log(*, s2: str = "", converged: bool = True,
         normal: bool = True, extra: str = "") -> str:
    parts = [
        " Entering Gaussian System\n",
        " SCF Done:  E(UB3LYP) =  -79.1234567890     A.U. after   12 cycles\n",
    ]
    if s2:
        parts.append(" Annihilation of the first spin contaminant:\n")
        parts.append(" S**2 before annihilation %s\n" % s2)
    if not converged:
        parts.append(" Convergence failure -- run terminated.\n")
    parts.append(_ORIENTATION)
    if extra:
        parts.append(extra)
    if normal:
        parts.append(" Normal termination of Gaussian 16\n")
    return "".join(parts)


def test_parse_gaussian_log_reads_s2_and_scf(workdir: Path):
    log = workdir / "mol.log"
    log.write_text(_log(s2="0.7601,   after     0.7500"), encoding="utf-8")
    summary = parse_gaussian_log(log)
    assert summary["scf_converged"] is True
    assert summary["scf_done_blocks"] == 1
    assert summary["s2_before_annihilation"] == pytest.approx(0.7601)
    assert summary["s2_after_annihilation"] == pytest.approx(0.7500)


def test_parse_gaussian_log_without_s2_lines(workdir: Path):
    log = workdir / "mol.log"
    log.write_text(_log(), encoding="utf-8")
    summary = parse_gaussian_log(log)
    assert summary["s2_before_annihilation"] is None
    assert summary["s2_after_annihilation"] is None
    assert summary["scf_converged"] is True


def test_parse_gaussian_log_flags_convergence_failure(workdir: Path):
    log = workdir / "mol.log"
    log.write_text(_log(converged=False), encoding="utf-8")
    summary = parse_gaussian_log(log)
    assert summary["scf_converged"] is False


def test_parse_gaussian_log_takes_last_s2_block(workdir: Path):
    """优化有多步，取最后一段（最终解）的自旋期望值。"""
    first = _log(s2="0.9000,   after     0.8900")
    second = _log(s2="0.7601,   after     0.7500")
    log = workdir / "mol.log"
    log.write_text(first + second, encoding="utf-8")
    summary = parse_gaussian_log(log)
    assert summary["s2_after_annihilation"] == pytest.approx(0.7500)
    assert summary["scf_done_blocks"] == 2


# --------------------------------------------------------------------------- #
# Stable=Opt 路由关键字
# --------------------------------------------------------------------------- #
def _env(root: Path) -> dict:
    return {"nproc": 4, "partition": "debug", "g16root": "/opt/g16",
            "g16_mem": "16G", "g16_walltime": "01:00:00"}


def test_stable_check_off_by_default(workdir: Path):
    xyz = workdir / "m.xyz"
    xyz.write_text("2\nm\nC 0.0 0.0 0.0\nH 1.09 0.0 0.0\n", encoding="utf-8")
    out = workdir / "mol"
    out.mkdir()
    make_gaussian_input(xyz, out, "mol", charge=0, multiplicity=2, env=_env(workdir))
    route = [line for line in (out / "mol.gjf").read_text(encoding="utf-8").splitlines()
             if line.startswith("#")][0]
    assert "Stable=Opt" not in route


def test_stable_check_appends_route_keyword(workdir: Path):
    xyz = workdir / "m.xyz"
    xyz.write_text("2\nm\nC 0.0 0.0 0.0\nH 1.09 0.0 0.0\n", encoding="utf-8")
    out = workdir / "mol"
    out.mkdir()
    make_gaussian_input(xyz, out, "mol", charge=0, multiplicity=2,
                        env=_env(workdir), stable_check=True)
    route = [line for line in (out / "mol.gjf").read_text(encoding="utf-8").splitlines()
             if line.startswith("#")][0]
    assert "Stable=Opt" in route
    assert "opt freq" in route


# --------------------------------------------------------------------------- #
# engine.collect：自旋污染进 notes，⟨S²⟩ 进 results
# --------------------------------------------------------------------------- #
class _Ctx:
    def __init__(self, run_dir: Path, *, multiplicity: int = 2) -> None:
        self.layout = SimpleNamespace(run_dir=Path(run_dir))
        self.config = SimpleNamespace(molecule=MoleculeConfig(
            charge=0, multiplicity=multiplicity))

    def require_stage_file(self, stage: str, name: str) -> Path:
        path = self.layout.run_dir / stage / name
        if not path.is_file():
            raise StageError("缺 %s" % path)
        return path

    def molecule_xyz(self):
        return None


def _collect(workdir: Path, *, multiplicity: int, s2: str = ""):
    run_dir = workdir / "run"
    stage_dir = run_dir / "molecule_opt"
    stage_dir.mkdir(parents=True)
    (stage_dir / "mol.log").write_text(_log(s2=s2), encoding="utf-8")
    return GaussianEngine().collect(_Ctx(run_dir, multiplicity=multiplicity), stage_dir, stage="molecule_opt")


def test_collect_reports_spin_contamination(workdir: Path):
    out = _collect(workdir, multiplicity=3, s2="2.5000,   after     2.1000")
    payload = out.results["gaussian"]
    assert payload["multiplicity"] == 3
    assert payload["s2_expected"] == pytest.approx(2.0)
    assert payload["s2_after_annihilation"] == pytest.approx(2.1)
    assert payload["spin_contamination"] == pytest.approx(0.1)
    assert payload["scf_converged"] is True


def test_collect_flags_large_contamination_in_notes(workdir: Path):
    out = _collect(workdir, multiplicity=3, s2="2.9000,   after     2.6000")
    assert any("自旋污染" in note for note in out.notes)


def test_collect_closed_shell_has_no_contamination(workdir: Path):
    out = _collect(workdir, multiplicity=1)
    payload = out.results["gaussian"]
    assert payload["s2_after_annihilation"] is None
    assert payload["spin_contamination"] is None
    assert payload["s2_expected"] == pytest.approx(0.0)


def test_collect_warns_when_scf_failed(workdir: Path):
    run_dir = workdir / "run"
    stage_dir = run_dir / "molecule_opt"
    stage_dir.mkdir(parents=True)
    (stage_dir / "mol.log").write_text(_log(converged=False), encoding="utf-8")
    out = GaussianEngine().collect(_Ctx(run_dir, multiplicity=1), stage_dir, stage="molecule_opt")
    assert out.results["gaussian"]["scf_converged"] is False
    assert any("SCF" in note for note in out.notes)
