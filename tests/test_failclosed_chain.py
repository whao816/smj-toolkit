"""第 1 步：堵住主链上仍然"绿色但不可信"的 fail-open。

锁四件事：

* ``aims_scf.collect`` 不再接受未正常结束 / 缺 ``basis-indices.out`` /
  ``omat.aims`` / ``mos.aims`` 的结果（这些是 tcontrol 的直接输入）；
* ``aitranss_te.prepare`` 要求 SCF 几何、A2/metadata 几何与 ``self.energy.in``
  依据的几何是同一份；
* ``aitranss_te.collect`` 默认要求 ``te.out`` 在场，只有显式
  ``transport.require_te_out: false`` 才豁免；
* collect 时从 ``aims.out`` / ``te.out`` 解析到的远端软件版本追加进
  ``run_spec.software.observed``（legacy run 落到 ``software_observed.json``）。
"""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from smj.config import load_config
from smj.errors import StageBlocked, StageError
from smj.provenance import sha256_file
from smj.stages.aims_scf import AimsScfStage
from smj.stages.aitranss_te import AitranssTeStage
from _aims_fixtures import write_scf_artifacts, write_scf_log


class _Ctx:
    """只实现 collect 需要的两个接口，避免把整个 runner 拖进单元测试。"""

    def __init__(self, run_dir: Path, stage_dir: Path, *,
                 require_te_out: bool = True) -> None:
        self.layout = SimpleNamespace(run_dir=Path(run_dir))
        self.config = SimpleNamespace(
            transport=SimpleNamespace(require_te_out=require_te_out))
        self._stage_dir = Path(stage_dir)

    def require_stage_file(self, stage: str, name: str) -> Path:
        path = self._stage_dir / name
        if not path.is_file():
            raise StageBlocked(stage, "缺 %s" % path, [str(path)])
        return path


def _seed_aims_scf(stage_dir: Path) -> None:
    stage_dir.mkdir(parents=True)
    (stage_dir / "geometry.in").write_text(
        "atom 0.0 0.0 0.0 Au\natom 1.0 0.0 0.0 Au\n", encoding="utf-8")
    write_scf_log(stage_dir)
    write_scf_artifacts(stage_dir)


def _seed_te(stage_dir: Path, *, te_out: bool = True) -> None:
    stage_dir.mkdir(parents=True)
    (stage_dir / "TE.dat").write_text(
        "# E(Ha)  E-EF(eV)  T\n"
        "-0.200000 0.000000 0.50\n"
        "-0.199800 0.010000 0.60\n",
        encoding="utf-8",
    )
    if te_out:
        (stage_dir / "te.out").write_text(
            "\n 2026-11-12 15:20:42.231\n\n"
            " invoking aitranss (version: 071813) ... \n"
            " compiled on 2026/09/23 at 15:38:33 on host gpu1\n",
            encoding="utf-8",
        )


def _write_run_spec(run_dir: Path) -> Path:
    path = run_dir / "run_spec.json"
    path.write_text(json.dumps({
        "schema": "smj.run_spec.v1",
        "software": {"executables": {}},
    }), encoding="utf-8")
    return path


# --------------------------------------------------------------------------- #
# aims_scf fail-closed
# --------------------------------------------------------------------------- #
def test_aims_scf_collect_records_versions(workdir: Path):
    run_dir = workdir / "run"
    stage_dir = run_dir / "aims_scf"
    _seed_aims_scf(stage_dir)
    spec_path = _write_run_spec(run_dir)

    out = AimsScfStage().collect(_Ctx(run_dir, stage_dir), stage_dir)

    assert out.results["aims_scf"]["normal_termination"] is True
    assert out.results["aims_scf"]["has_omat"] is True
    assert out.results["aims_scf"]["software"]["version"] == "221103"
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    observed = spec["software"]["observed"]["aims_scf"]
    assert observed["version"] == "221103"
    assert observed["aims_uuid"] == "TEST-UUID-1234"
    assert "te.out" not in spec["software"]["observed"]


def test_aims_scf_collect_rejects_abnormal_termination(workdir: Path):
    run_dir = workdir / "run"
    stage_dir = run_dir / "aims_scf"
    _seed_aims_scf(stage_dir)
    (stage_dir / "aims.out").write_text(
        "FHI-aims version      : 221103\nno goodbye here\n", encoding="utf-8")

    with pytest.raises(StageError, match="未正常结束"):
        AimsScfStage().collect(_Ctx(run_dir, stage_dir), stage_dir)


@pytest.mark.parametrize("missing", ["basis-indices.out", "omat.aims", "mos.aims"])
def test_aims_scf_collect_rejects_missing_tcontrol_inputs(workdir: Path, missing: str):
    run_dir = workdir / "run"
    stage_dir = run_dir / "aims_scf"
    _seed_aims_scf(stage_dir)
    (stage_dir / missing).unlink()

    with pytest.raises(StageError, match=missing):
        AimsScfStage().collect(_Ctx(run_dir, stage_dir), stage_dir)


def test_aims_scf_collect_rejects_empty_output(workdir: Path):
    run_dir = workdir / "run"
    stage_dir = run_dir / "aims_scf"
    _seed_aims_scf(stage_dir)
    (stage_dir / "aims.out").write_text("", encoding="utf-8")

    with pytest.raises(StageError, match="为空"):
        AimsScfStage().collect(_Ctx(run_dir, stage_dir), stage_dir)


# --------------------------------------------------------------------------- #
# aitranss_te collect：te.out 默认必须存在
# --------------------------------------------------------------------------- #
def test_aitranss_collect_requires_te_out_by_default(workdir: Path):
    run_dir = workdir / "run"
    stage_dir = run_dir / "aitranss_te"
    _seed_te(stage_dir, te_out=False)

    with pytest.raises(StageError, match="te.out"):
        AitranssTeStage().collect(
            _Ctx(run_dir, stage_dir, require_te_out=True), stage_dir)


def test_aitranss_collect_explicitly_allows_missing_te_out(workdir: Path):
    run_dir = workdir / "run"
    stage_dir = run_dir / "aitranss_te"
    _seed_te(stage_dir, te_out=False)

    out = AitranssTeStage().collect(
        _Ctx(run_dir, stage_dir, require_te_out=False), stage_dir)

    assert out.results["aitranss_te"]["has_te_out"] is False
    assert out.results["aitranss_te"]["require_te_out"] is False
    assert any("豁免" in note for note in out.notes)


def test_aitranss_collect_records_versions(workdir: Path):
    run_dir = workdir / "run"
    stage_dir = run_dir / "aitranss_te"
    _seed_te(stage_dir)
    spec_path = _write_run_spec(run_dir)

    out = AitranssTeStage().collect(_Ctx(run_dir, stage_dir), stage_dir)

    assert out.results["aitranss_te"]["software"]["version"] == "071813"
    assert "te.out" in out.artifacts
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    observed = spec["software"]["observed"]["aitranss_te"]
    assert observed["version"] == "071813"
    assert "gpu1" in observed["compiled"]


def test_observed_software_falls_back_to_file_for_legacy_run(workdir: Path):
    run_dir = workdir / "run"
    stage_dir = run_dir / "aitranss_te"
    _seed_te(stage_dir)

    AitranssTeStage().collect(_Ctx(run_dir, stage_dir), stage_dir)

    fallback = json.loads((run_dir / "software_observed.json").read_text(
        encoding="utf-8"))
    assert fallback["aitranss_te"]["version"] == "071813"
    assert not (run_dir / "run_spec.json").exists()


# --------------------------------------------------------------------------- #
# aitranss_te prepare：SCF 几何必须与上游/A2 几何一致
# --------------------------------------------------------------------------- #
def test_verify_geometry_chain_accepts_same_geometry(workdir: Path):
    geometry = workdir / "au_relax_geometry.in"
    scf_geometry = workdir / "aims_scf_geometry.in"
    text = "atom 0.0 0.0 0.0 Au\n"
    geometry.write_text(text, encoding="utf-8")
    scf_geometry.write_text(text, encoding="utf-8")
    metadata = workdir / "electrode_metadata.json"
    metadata.write_text(json.dumps({
        "surface_fingerprint": sha256_file(scf_geometry)}), encoding="utf-8")

    AitranssTeStage()._verify_geometry_chain(geometry, scf_geometry, metadata)


def test_verify_geometry_chain_rejects_scf_geometry_mismatch(workdir: Path):
    geometry = workdir / "au_relax_geometry.in"
    scf_geometry = workdir / "aims_scf_geometry.in"
    geometry.write_text("atom 0.0 0.0 0.0 Au\n", encoding="utf-8")
    scf_geometry.write_text("atom 1.0 0.0 0.0 Au\n", encoding="utf-8")
    metadata = workdir / "electrode_metadata.json"
    # 指纹与 SCF 几何一致，确保失败来自"上游几何与 SCF 几何不一致"这一条。
    metadata.write_text(json.dumps({
        "surface_fingerprint": sha256_file(scf_geometry)}), encoding="utf-8")

    with pytest.raises(StageBlocked, match="SCF 几何"):
        AitranssTeStage()._verify_geometry_chain(geometry, scf_geometry, metadata)


# --------------------------------------------------------------------------- #
# 配置入口
# --------------------------------------------------------------------------- #
def test_transport_require_te_out_parses(workdir: Path):
    (workdir / "molecule.xyz").write_text(
        "2\nH2\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")
    config = workdir / "smj.yaml"
    config.write_text(
        "run:\n  name: t\n  mode: dry-run\n  dir: .\n"
        "molecule:\n  xyz: molecule.xyz\n"
        "transport:\n  require_te_out: false\n",
        encoding="utf-8",
    )

    cfg = load_config(config)

    assert cfg.transport.require_te_out is False

