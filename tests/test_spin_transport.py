"""开壳层结的输运口径：aitranss α/β 双通道 + junction_analyze 自旋块。

2026-11-19 Phase 2：tcontrol 在 spin collinear 下用 ``$uhfmo_alpha/$uhfmo_beta``
读两套波函数；``TE.dat`` 的两列分别是 α/β 通道，collect/analyze 不再把信息
在均值里丢光。
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from smj.config import load_config
from smj.errors import StageBlocked
from smj.landauer_spin import analyze_te_spin, read_te_channels
from smj.provenance import sha256_file
from smj.stages.aitranss_te import AitranssTeStage
from smj.stages.junction_analyze import JunctionAnalyzeStage

EFERMI = -0.180814064867
_AIMS_GEOMETRY = "\n".join(
    "atom %.3f 0.0 0.0 Au" % x
    for x in (-3.0, -2.0, -1.0, -0.5, 0.5, 1.0, 2.0, 3.0)) + "\n"


def _project(workdir: Path, *, molecule: list[str]) -> Path:
    project = (workdir / "proj").resolve()
    project.mkdir(parents=True)
    (project / "molecule.xyz").write_text("2\n\nS 0 0 0\nS 3 0 0\n", encoding="utf-8")
    lines = [
        "run:", "  name: t", "  mode: dry-run", "  dir: %s" % project.as_posix(),
        "molecule:", "  xyz: %s" % (project / "molecule.xyz").as_posix(),
    ]
    lines.extend("  " + row for row in molecule)
    cfg = project / "smj.yaml"
    cfg.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return cfg


class _Ctx:
    def __init__(self, scf: Path, build: Path, config) -> None:
        self.config = config
        self._scf = scf
        self._build = build

    def override_path(self, name: str):
        return None

    def upstream(self, stage: str, name: str, override=None) -> Path:
        if override is not None:
            return override
        return self._scf / name

    def artifact(self, stage: str, name: str) -> Path:
        return self._build / name

    def require_stage_file(self, stage: str, name: str) -> Path:
        return self._scf / name

    def estimate(self, stage: str, natoms: int) -> dict:
        return {"atoms": natoms}

    def stage_env(self, stage: str, estimate: dict) -> dict:
        return {}

    def log(self, message: str) -> None:
        pass


def _prepare_case(workdir: Path, *, molecule: list[str]):
    run_dir = (workdir / "run").resolve()
    scf = run_dir / "aims_scf"
    build = run_dir / "electrode_build"
    stage_dir = run_dir / "aitranss_te"
    for path in (scf, build, stage_dir):
        path.mkdir(parents=True)
    geometry = scf / "geometry.in"
    geometry.write_text(_AIMS_GEOMETRY, encoding="utf-8")
    (build / "geometry.in").write_text(_AIMS_GEOMETRY, encoding="utf-8")
    (scf / "aims.out").write_text(
        "| Chemical potential (Fermi level):   %.12f Ha\n" % EFERMI, encoding="utf-8")
    (scf / "alpha.aims").write_text("nsaos = 4210\n", encoding="utf-8")
    (scf / "beta.aims").write_text("nsaos = 4210\n", encoding="utf-8")
    metadata = {
        "left_surface_indices": [1, 2, 3],
        "right_surface_indices": [5, 6, 7],
        "left_tip_index": 1,
        "right_tip_index": 5,
        "left_electrode_indices": [1, 2, 3, 4],
        "right_electrode_indices": [5, 6, 7, 8],
        "surface_fingerprint": sha256_file(geometry),
    }
    (build / "electrode_metadata.json").write_text(
        json.dumps(metadata), encoding="utf-8")
    return _project(workdir, molecule=molecule), scf, build, stage_dir


# --------------------------------------------------------------------------- #
# TE.dat 双通道解析
# --------------------------------------------------------------------------- #
_TE_SINGLE = "0.0 0.0 0.25\n0.1 0.1 0.5\n"
_TE_DOUBLE = "0.0 0.0 0.30 0.10\n0.1 0.1 0.60 0.20\n"


def test_read_te_channels_single_column(workdir):
    te = workdir / "TE.dat"
    te.write_text(_TE_SINGLE, encoding="utf-8")
    rows, channels = read_te_channels(te)
    assert channels == 1
    assert rows[0][2] == pytest.approx(0.25)
    assert rows[0][3] == pytest.approx(0.25)


def test_read_te_channels_two_columns(workdir):
    te = workdir / "TE.dat"
    te.write_text(_TE_DOUBLE, encoding="utf-8")
    rows, channels = read_te_channels(te)
    assert channels == 2
    assert [row[2] for row in rows] == pytest.approx([0.30, 0.60])
    assert [row[3] for row in rows] == pytest.approx([0.10, 0.20])


def test_analyze_te_spin_reports_polarization(workdir):
    te = workdir / "TE.dat"
    te.write_text(_TE_DOUBLE, encoding="utf-8")
    spin = analyze_te_spin(te)
    assert spin["channels"] == 2
    assert spin["T_alpha_at_ef"] == pytest.approx(0.30)
    assert spin["T_beta_at_ef"] == pytest.approx(0.10)
    assert spin["spin_polarization_at_ef"] == pytest.approx(0.5)
    assert spin["g_over_g0_total"] == pytest.approx(
        0.5 * (spin["L0_alpha"] + spin["L0_beta"]))


def test_analyze_te_spin_single_column_zero_polarization(workdir):
    te = workdir / "TE.dat"
    te.write_text(_TE_SINGLE, encoding="utf-8")
    spin = analyze_te_spin(te)
    assert spin["channels"] == 1
    assert spin["spin_polarization_at_ef"] == pytest.approx(0.0)
    assert "note" in spin


# --------------------------------------------------------------------------- #
# aitranss_te.prepare：uhfmo 双波函数 + 价电子扣电荷
# --------------------------------------------------------------------------- #
def test_open_shell_prepare_writes_uhfmo(workdir):
    cfg, scf, build, stage_dir = _prepare_case(
        workdir, molecule=["multiplicity: 2"])
    config = load_config(cfg)
    ctx = _Ctx(scf, build, config)
    result = AitranssTeStage().prepare(ctx, stage_dir)
    text = (stage_dir / "tcontrol").read_text(encoding="utf-8")
    assert "$uhfmo_alpha file=alpha.aims" in text
    assert "$uhfmo_beta file=beta.aims" in text
    assert "$scfmo" not in text
    summary = json.loads((stage_dir / "input_summary.json").read_text(encoding="utf-8"))
    assert summary["spin_polarized"] is True
    assert summary["wavefunction_files"] == ["alpha.aims", "beta.aims"]
    assert any("uhfmo" in note for note in result.notes)


def test_charged_junction_adjusts_valence_electrons(workdir):
    cfg, scf, build, stage_dir = _prepare_case(
        workdir, molecule=["multiplicity: 2", "charge: -1"])
    config = load_config(cfg)
    ctx = _Ctx(scf, build, config)
    AitranssTeStage().prepare(ctx, stage_dir)
    text = (stage_dir / "tcontrol").read_text(encoding="utf-8")
    # 8 个 Au × 31 = 248 价电子；电荷 −1 → 249
    assert "$valence_electrons  249" in text


def test_open_shell_prepare_blocks_without_alpha_beta(workdir):
    cfg, scf, build, stage_dir = _prepare_case(
        workdir, molecule=["multiplicity: 2"])
    (scf / "alpha.aims").unlink()
    (scf / "beta.aims").unlink()
    config = load_config(cfg)
    ctx = _Ctx(scf, build, config)
    with pytest.raises(StageBlocked, match="alpha.aims"):
        AitranssTeStage().prepare(ctx, stage_dir)


# --------------------------------------------------------------------------- #
# junction_analyze：开壳层结果里带 spin 块
# --------------------------------------------------------------------------- #
def test_junction_analyze_writes_spin_block(workdir):
    cfg, _, _, stage_dir = _prepare_case(workdir, molecule=["multiplicity: 2"])
    te = stage_dir / "TE.dat"
    te.write_text(_TE_DOUBLE, encoding="utf-8")
    config = load_config(cfg)
    ctx = _Ctx(stage_dir, stage_dir, config)
    result = JunctionAnalyzeStage().prepare(ctx, stage_dir)
    data = json.loads((stage_dir / "results.json").read_text(encoding="utf-8"))
    assert data["spin"]["channels"] == 2
    assert data["spin"]["spin_polarization_at_ef"] == pytest.approx(0.5)
    assert data["points"] == 2
    assert any("自旋通道" in note for note in result.notes)
