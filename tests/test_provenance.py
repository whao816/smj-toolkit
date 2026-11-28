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


@pytest.mark.parametrize("legacy_setting", [True, False])
@pytest.mark.parametrize("old_transport", [False, True])
def test_legacy_bs_hash_resumes_only_released_run(workdir, legacy_setting, old_transport):
    from smj.provenance import (canonical_json, ensure_run_spec, protocol_components,
                                run_spec_payload, sha256_text)

    cfg = _config(workdir)
    mol = cfg.molecule
    mol.broken_symmetry = True
    mol.initial_moments = {1: 1.0, 2: -1.0}
    mol.fixed_spin_moment = False
    snapshot = run_spec_payload(cfg, workdir)
    snapshot["protocol_schema"] = 2
    snapshot["config"]["molecule"]["fixed_spin_moment"] = legacy_setting
    protocol = protocol_components(cfg)
    if legacy_setting:
        protocol["molecule"].pop("fixed_spin_moment", None)
    else:
        protocol["molecule"]["fixed_spin_moment"] = False
    if old_transport:
        snapshot["config"]["transport"].pop("allow_single_channel")
        protocol["transport"].pop("allow_single_channel")
    snapshot["protocol_sha256"] = sha256_text(canonical_json(protocol))
    (workdir / "run_spec.json").write_text(json.dumps(snapshot), encoding="utf-8")
    assert ensure_run_spec(cfg, workdir) == snapshot
    mol.fixed_spin_moment = True
    with pytest.raises(ConfigError, match="protocol_sha256"):
        ensure_run_spec(cfg, workdir)


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
    stage_dir = runner.layout.stage_dir("gaussian_opt")
    # 断言按磁盘字节哈希（Windows 的 write_text 会把 \n 写成 \r\n）。
    assert _sha256_file(stage_dir / "input.txt") in set(stage.inputs_sha256.values())
    assert _sha256_file(workdir / "molecule.xyz") in set(stage.inputs_sha256.values())

    runner.collect("gaussian_opt")
    stage = runner.state.stage("gaussian_opt")
    assert stage.outputs_sha256["out.txt"] == _sha256_file(stage_dir / "out.txt")
    saved = json.loads((runner.run_dir / "state.json").read_text(encoding="utf-8"))
    assert saved["stages"]["molecule_opt"]["outputs_sha256"]["out.txt"] == \
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


def test_run_spec_keeps_remote_posix_paths_verbatim(workdir):
    """slurm-ssh 的 POSIX 路径不能被 Windows Path() 改写成反斜杠。

    2026-11-13 的验收 run 里 run_spec 把 /data/apps/... 记成 \\data\\apps\\... 且
    exists=false，读起来像"二进制缺失"。远端路径应原样记录，exists 记 None
    （本机不知道），而不是谎报 False。
    """
    from smj.provenance import run_spec_payload

    config = _config(workdir)
    config.run.mode = "slurm-ssh"
    config.cluster.env = {"AIMS_BIN": "/data/apps/aims/bin/aims.x"}
    payload = run_spec_payload(config, workdir / "run")

    entry = payload["software"]["executables"]["AIMS_BIN"]
    assert entry["path"] == "/data/apps/aims/bin/aims.x"
    assert entry["exists"] is None


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


def test_guard_geometry_fingerprint_accepts_structural_when_tags_added(workdir):
    """开壳层/带电 SCF 往 geometry.in 追加 initial_* 子 tag 后，
    字节 sha256 变了，但结构没变 —— 结构指纹仍应通过。"""
    from smj.provenance import structural_fingerprint, verify_parent_geometry

    base = "atom 0.0 0.0 0.0 Au\natom 1.0 0.0 0.0 S\n"
    geometry = workdir / "geometry.in"
    geometry.write_text(base, encoding="utf-8")
    tagged = base + "  initial_charge  0.5\n  initial_moment  0.1\n"
    scf_geometry = workdir / "scf_geometry.in"
    scf_geometry.write_text(tagged, encoding="utf-8")

    meta = {"surface_fingerprint": _sha256_file(geometry),
            "surface_fingerprint_structural": structural_fingerprint(geometry)}
    metadata = workdir / "electrode_metadata.json"
    metadata.write_text(json.dumps(meta), encoding="utf-8")

    assert structural_fingerprint(geometry) == structural_fingerprint(scf_geometry)
    assert verify_parent_geometry("aims_scf", scf_geometry, metadata)

    meta_bad = dict(meta, surface_fingerprint_structural="0" * 64)
    metadata.write_text(json.dumps(meta_bad), encoding="utf-8")
    with pytest.raises(StageBlocked, match="结构指纹"):
        verify_parent_geometry("aims_scf", scf_geometry, metadata)


def test_structural_fingerprint_ignores_initial_tags_only(workdir):
    from smj.provenance import structural_fingerprint

    a = workdir / "a.in"
    a.write_text("atom 0 0 0 Au\n", encoding="utf-8")
    b = workdir / "b.in"
    b.write_text("atom 0 0 0 Au\n  initial_charge 1.0\n", encoding="utf-8")
    c = workdir / "c.in"
    c.write_text("atom 0 0 0 Ag\n", encoding="utf-8")

    assert structural_fingerprint(a) == structural_fingerprint(b)
    assert structural_fingerprint(a) != structural_fingerprint(c)
    assert structural_fingerprint(workdir / "missing.in") is None


def test_discover_stages_uses_state_order_not_vendored(
        workdir, monkeypatch):
    from smj import report as report_mod
    import smj.junction.manifest as vendored

    run_dir = workdir / "run"
    for stage in STAGE_ORDER:
        (run_dir / stage).mkdir(parents=True)
    monkeypatch.setattr(vendored, "discover_smj_stages",
                        lambda root: ["gaussian_opt", "au60_build"], raising=False)

    assert report_mod.discover_stages(run_dir) == list(STAGE_ORDER)


def test_discover_stages_collapses_legacy_auN_build_dirs(workdir):
    from smj.report import discover_stages

    run_dir = workdir / "run"
    for stage in STAGE_ORDER:
        (run_dir / stage).mkdir(parents=True)
    # 旧 run 的 au20_build 目录由 RunLayout 回退到 electrode_build，不应再作为
    # 额外阶段出现在阶段序列里。
    (run_dir / "au20_build").mkdir()

    assert discover_stages(run_dir) == list(STAGE_ORDER)

def test_review_flags_protocol_drift(workdir):
    config = _config(workdir)
    runner = _runner(config)
    runner.prepare("gaussian_opt")

    from smj.review import review_run
    changed = _config(workdir, temperature=250.0)
    result = review_run(changed, force=True)

    assert "provenance_drift" in result.codes()
    assert result.status == "blocked"

# --- 2026-11-20 复核 P2-1 / R8：采样覆盖必须进协议指纹 ------------------ #

def _protocol_config_of(config: RunConfig) -> dict:
    from smj.provenance import _protocol_config
    return _protocol_config(config)


def test_protocol_hash_ignores_explicit_estep_override(workdir):
    """复核实测：eV 相对窗下改 transport.estep，哈希必须变。"""
    from smj.provenance import protocol_sha256

    base = _config(workdir)
    overridden = _config(workdir)
    overridden.transport.estep = 0.001

    assert protocol_sha256(base) != protocol_sha256(overridden)
    assert _protocol_config_of(overridden)["transport"]["estep"] == 0.001


def test_protocol_hash_ignores_explicit_window_edges(workdir):
    from smj.provenance import protocol_sha256

    base = _config(workdir)
    changed = _config(workdir)
    changed.transport.ener = -0.2
    changed.transport.eend = -0.1

    assert protocol_sha256(base) != protocol_sha256(changed)
    transport = _protocol_config_of(changed)["transport"]
    assert transport["ener"] == -0.2
    assert transport["eend"] == -0.1


def test_protocol_default_config_keeps_old_shape(workdir):
    """默认配置（ener/estep/eend 皆 None）不得新增键——否则历史 run_spec 集体漂移。"""
    transport = _protocol_config_of(_config(workdir))["transport"]
    for key in ("ener", "estep", "eend"):
        assert key not in transport
    assert transport["energy_window_unit"] == "eV"


def test_run_spec_records_resolved_window_unit(workdir):
    config = _config(workdir)
    runner = _runner(config)
    runner.prepare("gaussian_opt")
    transport = _spec(runner.run_dir)["config"]["transport"]
    # eV 口径需要 efermi 才能解析；resolve 不了时该键可以缺席，出现就必须带单位。
    if "resolved_window" in transport:
        assert transport["resolved_window_unit"] == "Hartree"


# --- 2026-11-22 P0-5：执行面（species/templates/二进制）必须进 run_spec ------ #

def _execution_config(workdir: Path):
    templates = workdir / "templates"
    templates.mkdir(exist_ok=True)
    script = templates / "sub.sh"
    script.write_text("v1\n", encoding="utf-8")
    aims = workdir / "aims.bin"
    aims.write_text("aims-v1", encoding="utf-8")
    config = _config(workdir)
    config.paths.templates = str(templates)
    config.cluster.env["aims_bin"] = str(aims)
    return config, templates, aims


def test_run_spec_records_execution_face(workdir):
    """species defaults/模板/二进制的内容哈希必须随 run_spec 落盘。"""
    from smj.provenance import execution_sha256

    config, templates, aims = _execution_config(workdir)
    runner = _runner(config)
    runner.prepare("gaussian_opt")

    spec = _spec(runner.run_dir)
    assert spec["execution_sha256"] == execution_sha256(config)
    face = spec["execution"]
    assert face["templates"]["exists"] is True
    assert face["templates"]["files"] == 1
    assert face["executables"]["aims_bin"]["sha256"] == _sha256_file(aims)

    before = spec["execution_sha256"]
    (templates / "sub.sh").write_text("v2\n", encoding="utf-8")
    assert execution_sha256(config) != before


def test_prepare_blocks_execution_drift_and_force_rebinds(workdir):
    """改模板/二进制不改变 protocol，但必须触发执行面漂移门。"""
    config, templates, _ = _execution_config(workdir)
    _runner(config).prepare("gaussian_opt")
    first = _spec(_runner(config).run_dir)["execution_sha256"]

    (templates / "sub.sh").write_text("v2\n", encoding="utf-8")
    with pytest.raises(ConfigError, match="执行面"):
        _runner(config).prepare("gaussian_opt")

    _runner(config).prepare("gaussian_opt", force=True)
    rebound = _spec(_runner(config).run_dir)
    assert rebound["rebound"] is True
    assert rebound["execution_sha256"] != first
    assert rebind_path_exists(_runner(config).run_dir)


def rebind_path_exists(run_dir: Path) -> bool:
    return (run_dir / "run_spec.previous.json").is_file()


def test_legacy_run_spec_backfills_execution_face(workdir):
    """旧 run_spec 只有 protocol：补齐执行面字段，但不改写已知的物理协议。"""
    from smj.provenance import ensure_run_spec, execution_sha256

    config, _, _ = _execution_config(workdir)
    runner = _runner(config)
    runner.prepare("gaussian_opt")
    path = runner.run_dir / "run_spec.json"
    legacy = json.loads(path.read_text(encoding="utf-8"))
    protocol = legacy["protocol_sha256"]
    legacy.pop("execution_sha256")
    legacy.pop("execution")
    path.write_text(json.dumps(legacy, indent=2) + "\n", encoding="utf-8")

    restored = ensure_run_spec(config, runner.run_dir)
    assert restored["protocol_sha256"] == protocol
    assert restored["execution_sha256"] == execution_sha256(config)
    assert restored["execution_backfilled_at"]
    assert "无法追溯" in restored["execution_backfill_note"]


def test_execution_face_ignores_scheduler_env(workdir):
    """account/partition 等调度键不改变执行面，不该阻塞续跑。"""
    from smj.provenance import ensure_run_spec

    config, _, _ = _execution_config(workdir)
    config.cluster.env["account"] = "acct-a"
    _runner(config).prepare("gaussian_opt")

    config.cluster.env["account"] = "acct-b"
    again = ensure_run_spec(config, _runner(config).run_dir)
    assert again.get("rebound") is not True

