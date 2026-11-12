"""P0 溯源/阶段清单回归：run_spec.json、哈希闭环、report 阶段发现。

背景：架构评审指出三处可复现性缺口——
1. 运行参数（含默认值）、软件版本、输入哈希没有随 run 落盘；
2. report 用 vendored ``discover_smj_stages``（缺 au_relax），与
   ``state.STAGE_ORDER`` 形成双源；
3. 阶段之间只校验"文件存在"，没有父/子哈希，几何被替换也无法发现。

本文件只测行为，不测内部实现；所有路径都在 conftest 的 ``workdir`` 下。
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from smj.config import (
    ClusterConfig,
    ElectrodeConfig,
    MoleculeConfig,
    RunConfig,
    RunInfo,
    TransportConfig,
)
from smj.errors import ConfigError, StageBlocked
from smj.runner import StageRunner
from smj.stages.base import Collected, Prepared
from smj.state import STAGE_ORDER


class _FakeStage:
    """最小阶段：prepare 写一个输入，collect 写一个输出。"""

    name = "gaussian_opt"
    script = "sub_g16.sh"
    outputs = ("out.txt",)
    needs_submit = True

    def prepare(self, ctx, stage_dir):
        stage_dir = Path(stage_dir)
        (stage_dir / "input.txt").write_text("hello\n", encoding="utf-8")
        return Prepared(
            stage=self.name,
            inputs=[stage_dir / "input.txt", ctx.molecule_xyz()],
            artifacts=["out.txt"],
            resources={},
            script=self.script,
        )

    def collect(self, ctx, stage_dir):
        stage_dir = Path(stage_dir)
        (stage_dir / "out.txt").write_text("done\n", encoding="utf-8")
        return Collected(stage=self.name, artifacts=["out.txt"],
                         results={"ok": True}, output_dir=stage_dir)


def _config(workdir: Path, *, temperature: float = 300.0) -> RunConfig:
    xyz = (workdir / "molecule.xyz").resolve()
    xyz.write_text("4\nsmoke\nH 0 0 0\nH 0 0 1\nH 0 1 0\nH 1 0 0\n", encoding="utf-8")
    return RunConfig(
        run=RunInfo(name="t", dir=".", mode="dry-run"),
        molecule=MoleculeConfig(xyz=str(xyz)),
        electrode=ElectrodeConfig(au_per_electrode=8, source="builtin"),
        transport=TransportConfig(temperature=temperature),
        cluster=ClusterConfig(host="h", user="u", remote_root="/remote",
                              partition="p"),
        source=workdir / "smj.yaml",
    )


def _runner(config: RunConfig) -> StageRunner:
    return StageRunner(
        config, logger=lambda message: None, sleep=lambda seconds: None,
        stage_getter=lambda name: _FakeStage())


def _spec(run_dir: Path) -> dict:
    return json.loads((run_dir / "run_spec.json").read_text(encoding="utf-8"))


def _sha256_file(path: Path) -> str:
    import hashlib
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_prepare_writes_run_spec_with_resolved_config_and_protocol(workdir):
    config = _config(workdir)
    runner = _runner(config)
    runner.prepare("gaussian_opt")

    spec = _spec(runner.run_dir)
    assert spec["schema"] == "smj.run_spec.v1"
    assert spec["reconstructed"] is False
    # 完整解析配置：默认值也在里面，不能只记用户显式写过的键。
    assert spec["config"]["transport"]["temperature"] == 300.0
    assert spec["config"]["molecule"]["method"] == "b3lyp"
    assert spec["config"]["electrode"]["au_per_electrode"] == 8
    assert spec["config"]["cluster"]["partition"] == "p"
    # 输入与协议指纹
    assert spec["molecule"]["sha256"] == _sha256_file(workdir / "molecule.xyz")
    from smj.provenance import protocol_sha256
    assert spec["protocol_sha256"] == protocol_sha256(config)
    assert spec["software"]["smj"]
    assert spec["software"]["python"]


def test_prepare_refuses_config_drift_and_force_rebinds(workdir):
    first = _runner(_config(workdir, temperature=300.0))
    first.prepare("gaussian_opt")
    first_protocol = _spec(first.run_dir)["protocol_sha256"]

    changed = _runner(_config(workdir, temperature=250.0))
    with pytest.raises(ConfigError, match="run_spec"):
        changed.prepare("gaussian_opt")

    changed.prepare("gaussian_opt", force=True)
    rebound = _spec(changed.run_dir)
    assert rebound["config"]["transport"]["temperature"] == 250.0
    assert rebound.get("rebound") is True
    assert (changed.run_dir / "run_spec.previous.json").is_file()
    assert rebound["protocol_sha256"] != first_protocol


def test_manifest_binds_run_spec(workdir):
    config = _config(workdir)
    runner = _runner(config)
    runner.prepare("gaussian_opt")

    from smj.report import build_manifest
    manifest = build_manifest(config)

    assert manifest["spec"]["transport"]["temperature"] == 300.0
    assert manifest["inputs"]["mode"] == "dry-run"
    assert manifest["inputs"]["temperature_K"] == 300.0
    assert manifest["run_spec"]["schema"] == "smj.run_spec.v1"
    assert manifest["provenance"]["reconstructed"] is False
    assert manifest["provenance"]["run_spec_sha256"] == _sha256_file(
        runner.run_dir / "run_spec.json")


def test_manifest_marks_reconstructed_spec_when_run_spec_missing(workdir):
    config = _config(workdir)
    runner = _runner(config)
    runner.prepare("gaussian_opt")
    (runner.run_dir / "run_spec.json").unlink()

    from smj.report import build_manifest
    manifest = build_manifest(config)

    assert manifest["run_spec"]["reconstructed"] is True
    assert manifest["provenance"]["reconstructed"] is True
    assert manifest["inputs"]["mode"] == "dry-run"
    assert manifest["inputs"]["temperature_K"] == 300.0
    assert any("run_spec.json 缺失" in warning for warning in manifest["warnings"])


def test_prepare_and_collect_record_input_output_hashes(workdir):
    config = _config(workdir)
    runner = _runner(config)
    runner.prepare("gaussian_opt")

    stage = runner.state.stage("gaussian_opt")
    stage_dir = runner.run_dir / "gaussian_opt"
    # 断言按磁盘字节哈希（Windows 的 write_text 会把 \n 写成 \r\n）。
    assert _sha256_file(stage_dir / "input.txt") in set(stage.inputs_sha256.values())
    assert _sha256_file(workdir / "molecule.xyz") in set(stage.inputs_sha256.values())

    runner.collect("gaussian_opt")
    stage = runner.state.stage("gaussian_opt")
    assert stage.outputs_sha256["out.txt"] == _sha256_file(stage_dir / "out.txt")
    saved = json.loads((runner.run_dir / "state.json").read_text(encoding="utf-8"))
    assert saved["stages"]["gaussian_opt"]["outputs_sha256"]["out.txt"] == \
        _sha256_file(stage_dir / "out.txt")


def test_run_spec_hashes_electrode_library(workdir):
    """库电极的内容哈希进入 protocol，防止同名团簇换文件被当成同一输入。"""
    from smj.provenance import protocol_sha256, sha256_file

    rel = "fcc111.au.clusters.37_atoms/fcc111.cluster.natoms37.x1.aims"
    lib = workdir / rel
    lib.parent.mkdir(parents=True)
    lib.write_text("atom 0 0 0 Au\n", encoding="utf-8")
    config = _config(workdir)
    config.electrode = ElectrodeConfig(source="library", library_dir=str(workdir),
                                        library_cluster=rel)
    runner = _runner(config)
    runner.prepare("gaussian_opt")

    spec = _spec(runner.run_dir)
    assert spec["electrode_library"]["sha256"] == sha256_file(lib)
    assert spec["protocol_sha256"] == protocol_sha256(config)

    lib.write_text("atom 0 0 1 Au\n", encoding="utf-8")
    assert protocol_sha256(config) != spec["protocol_sha256"]


def test_guard_geometry_fingerprint_rejects_mismatch(workdir):
    from smj.provenance import verify_parent_geometry

    geometry = workdir / "geometry.in"
    geometry.write_text("atom 0 0 0 Au\n", encoding="utf-8")
    metadata = workdir / "electrode_metadata.json"
    metadata.write_text(json.dumps({"surface_fingerprint": "deadbeef"}),
                        encoding="utf-8")

    with pytest.raises(StageBlocked, match="指纹"):
        verify_parent_geometry("aims_scf", geometry, metadata)

    metadata.write_text(json.dumps(
        {"surface_fingerprint": _sha256_file(geometry)}), encoding="utf-8")
    assert verify_parent_geometry("aims_scf", geometry, metadata)


def test_discover_stages_includes_au_relax_even_if_vendored_lists_are_stale(
        workdir, monkeypatch):
    from smj import report as report_mod
    import smj._vendor.junction.manifest as vendored

    run_dir = workdir / "run"
    for stage in STAGE_ORDER:
        (run_dir / stage).mkdir(parents=True)
    monkeypatch.setattr(vendored, "discover_smj_stages",
                        lambda root: ["gaussian_opt", "au60_build"], raising=False)

    assert report_mod.discover_stages(run_dir) == list(STAGE_ORDER)


def test_discover_stages_keeps_extra_auN_build_dirs(workdir):
    from smj.report import discover_stages

    run_dir = workdir / "run"
    for stage in STAGE_ORDER:
        (run_dir / stage).mkdir(parents=True)
    (run_dir / "au20_build").mkdir()

    expected = list(STAGE_ORDER)
    expected.insert(expected.index("au60_build") + 1, "au20_build")
    assert discover_stages(run_dir) == expected


def test_review_flags_protocol_drift(workdir):
    config = _config(workdir)
    runner = _runner(config)
    runner.prepare("gaussian_opt")

    from smj.review import review_run
    changed = _config(workdir, temperature=250.0)
    result = review_run(changed, force=True)

    assert "provenance_drift" in result.codes()
    assert result.status == "blocked"
