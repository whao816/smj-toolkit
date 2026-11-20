"""输运能量窗口径：默认 eV 且相对费米能级（E_F ±2 eV），旧 Ha 绝对口径兼容。

锁四件事：

* 默认 ``energy_window`` 是 eV 相对口径（[-2, 0.005, +2] eV），没有 efermi 时
  换算成 aitranss 的 Ha 绝对窗口必须 fail-closed，不猜；
* 0.2.x 的 Ha 绝对三元组按量级识别（跨度 ≤0.1 Ha、量级 <2 Ha），原样生效并提示；
* ``energy_window_unit`` 显式声明优先于量级探测；
* ``aitranss_te.prepare`` 真的把 eV 口径按 ``aims.out`` 的费米能级平移到
  tcontrol 的 Ha 绝对值。
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from smj import provenance
from smj.config import HARTREE_EV, load_config
from smj.errors import ConfigError, StageBlocked
from smj.provenance import sha256_file
from smj.stages.aitranss_te import AitranssTeStage

EFERMI = -0.180814064867


def _project(workdir: Path, *, transport: list[str] | None = None) -> Path:
    project = (workdir / "proj").resolve()
    project.mkdir(parents=True)
    (project / "molecule.xyz").write_text("2\n\nC 0 0 0\nH 1 0 0\n", encoding="utf-8")
    lines = [
        "run:", "  name: t", "  mode: dry-run", "  dir: %s" % project.as_posix(),
        "molecule:", "  xyz: %s" % (project / "molecule.xyz").as_posix(),
    ]
    if transport is not None:
        lines.append("transport:")
        lines.extend("  " + row for row in transport)
    cfg = project / "smj.yaml"
    cfg.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return cfg


# --------------------------------------------------------------------------- #
# 配置语义
# --------------------------------------------------------------------------- #
def test_default_window_is_ev_relative_to_fermi(workdir: Path):
    config = load_config(_project(workdir))
    transport = config.transport
    assert transport.energy_window == (-2.0, 0.005, 2.0)
    assert transport.energy_window_unit == "eV"
    assert transport.relative_window() == (-2.0, 0.005, 2.0)
    # 没有 efermi 不得静默换算出 Ha 绝对窗口
    with pytest.raises(ValueError, match="efermi"):
        transport.resolved_window()


def test_ev_window_is_shifted_by_efermi(workdir: Path):
    config = load_config(_project(
        workdir, transport=["energy_window: [-1.0, 0.01, 1.0]"])).transport
    assert config.energy_window_unit == "eV"
    ener, estep, eend = config.resolved_window(EFERMI)
    assert ener == pytest.approx(EFERMI - 1.0 / HARTREE_EV, abs=1e-12)
    assert eend == pytest.approx(EFERMI + 1.0 / HARTREE_EV, abs=1e-12)
    assert estep == pytest.approx(0.01 / HARTREE_EV, abs=1e-15)


def test_legacy_hartree_window_detected_and_kept(workdir: Path):
    config = load_config(_project(
        workdir, transport=["energy_window: [-0.21, 0.0002, -0.17]"])).transport
    assert config.energy_window_unit == "Hartree"
    assert config.resolved_window() == pytest.approx((-0.21, 0.0002, -0.17))
    assert any("旧口径" in w for w in load_config(
        _project(workdir / "again", transport=["energy_window: [-0.21, 0.0002, -0.17]"])
    ).warnings)


def test_explicit_hartree_unit_beats_heuristic(workdir: Path):
    config = load_config(_project(workdir, transport=[
        "energy_window: [-0.30, 0.001, -0.10]",  # 跨度 0.2 Ha：量级探测会误判成 eV
        "energy_window_unit: Hartree",
    ])).transport
    assert config.energy_window_unit == "Hartree"
    assert config.resolved_window() == pytest.approx((-0.30, 0.001, -0.10))


def test_explicit_ev_unit_on_narrow_window(workdir: Path):
    config = load_config(_project(workdir, transport=[
        "energy_window: [-0.1, 0.001, 0.1]",
        "energy_window_unit: eV",
    ]))
    assert config.transport.energy_window_unit == "eV"
    assert config.transport.resolved_window(EFERMI)[0] == pytest.approx(
        EFERMI - 0.1 / HARTREE_EV, abs=1e-12)


def test_invalid_unit_rejected(workdir: Path):
    with pytest.raises(ConfigError, match="energy_window_unit"):
        load_config(_project(workdir, transport=["energy_window_unit: hartreee"]))


def test_partial_override_stays_ha_absolute(workdir: Path):
    config = load_config(_project(workdir, transport=[
        "energy_window: [-2.0, 0.005, 2.0]",
        "ener: -0.19",  # Ha 绝对：只覆盖起点
    ])).transport
    ener, estep, eend = config.resolved_window(EFERMI)
    assert ener == pytest.approx(-0.19)                    # 覆盖值不作单位换算
    assert estep == pytest.approx(0.005 / HARTREE_EV)      # 步长仍按 eV 换算
    assert eend == pytest.approx(EFERMI + 2.0 / HARTREE_EV)


# --------------------------------------------------------------------------- #
# protocol 指纹：旧 Ha 口径保持历史形状，eV 口径显式标单位
# --------------------------------------------------------------------------- #
def test_protocol_keeps_legacy_shape_for_hartree(workdir: Path):
    legacy = load_config(_project(
        workdir, transport=["energy_window: [-0.21, 0.0002, -0.17]"]))
    transport = provenance._protocol_config(legacy)["transport"]
    assert transport["energy_window"] == pytest.approx([-0.21, 0.0002, -0.17])
    assert "energy_window_unit" not in transport


def test_protocol_records_ev_triplet(workdir: Path):
    modern = load_config(_project(
        workdir, transport=["energy_window: [-1.0, 0.01, 1.0]"]))
    transport = provenance._protocol_config(modern)["transport"]
    assert transport["energy_window"] == [-1.0, 0.01, 1.0]
    assert transport["energy_window_unit"] == "eV"


# --------------------------------------------------------------------------- #
# aitranss_te.prepare：eV -> Ha 绝对值真的写进 tcontrol
# --------------------------------------------------------------------------- #
_AIMS_GEOMETRY = "\n".join(
    "atom %.3f 0.0 0.0 Au" % x
    for x in (-3.0, -2.0, -1.0, -0.5, 0.5, 1.0, 2.0, 3.0)) + "\n"


class _Ctx:
    """只实现 prepare 需要的接口，避免把整个 runner 拖进单元测试。"""

    def __init__(self, scf: Path, build: Path, config) -> None:
        self.config = config
        self._scf = scf
        self._build = build

    def override_path(self, name: str):
        return None

    def upstream(self, stage: str, name: str) -> Path:
        return self._scf / name

    def artifact(self, stage: str, name: str) -> Path:
        return self._build / name

    def estimate(self, stage: str, natoms: int) -> dict:
        return {"atoms": natoms}

    def stage_env(self, stage: str, estimate: dict) -> dict:
        return {}

    def log(self, message: str) -> None:
        pass


def _prepare_case(workdir: Path, *, transport: list[str] | None = None,
                  ferry: bool = True) -> tuple[Path, Path]:
    run_dir = (workdir / "run").resolve()
    scf = run_dir / "aims_scf"
    build = run_dir / "electrode_build"
    stage_dir = run_dir / "aitranss_te"
    scf.mkdir(parents=True)
    build.mkdir(parents=True)
    stage_dir.mkdir(parents=True)
    geometry = scf / "geometry.in"
    geometry.write_text(_AIMS_GEOMETRY, encoding="utf-8")
    (build / "geometry.in").write_text(_AIMS_GEOMETRY, encoding="utf-8")
    if ferry:
        (scf / "aims.out").write_text(
            "| Chemical potential (Fermi level):   %.12f Ha\n" % EFERMI,
            encoding="utf-8")
    (scf / "mos.aims").write_text("nsaos = 4210\n", encoding="utf-8")
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
    cfg = _project(workdir, transport=transport)
    return cfg, stage_dir


def test_prepare_shifts_ev_window_by_scf_fermi(workdir: Path):
    cfg, stage_dir = _prepare_case(workdir)
    config = load_config(cfg)
    ctx = _Ctx(stage_dir.parent / "aims_scf", stage_dir.parent / "electrode_build",
               config)

    result = AitranssTeStage().prepare(ctx, stage_dir)

    ener = EFERMI - 2.0 / HARTREE_EV
    eend = EFERMI + 2.0 / HARTREE_EV
    text = (stage_dir / "tcontrol").read_text(encoding="utf-8")
    assert "%.10f" % ener in text
    assert "%.10f" % eend in text
    assert result.stage == "aitranss_te"
    assert any("E−E_F" in note and "eV" in note for note in result.notes)


def test_prepare_blocks_without_fermi_level(workdir: Path):
    cfg, stage_dir = _prepare_case(workdir, ferry=False)
    config = load_config(cfg)
    ctx = _Ctx(stage_dir.parent / "aims_scf", stage_dir.parent / "electrode_build",
               config)
    with pytest.raises(StageBlocked, match="efermi"):
        AitranssTeStage().prepare(ctx, stage_dir)


def test_prepare_blocks_window_not_covering_fermi(workdir: Path):
    """2026-11-20 复核 P1-2：整窗在 E_F 一侧时必须在生成输入前拦住。"""
    cfg, stage_dir = _prepare_case(
        workdir, transport=["energy_window: [0.25, 0.0005, 0.45]"])
    config = load_config(cfg)
    assert config.transport.energy_window_unit == "eV"
    ctx = _Ctx(stage_dir.parent / "aims_scf", stage_dir.parent / "electrode_build",
               config)
    with pytest.raises(StageBlocked, match="未覆盖 E_F"):
        AitranssTeStage().prepare(ctx, stage_dir)
    assert not (stage_dir / "tcontrol").exists()


def test_prepare_keeps_legacy_hartree_window(workdir: Path):
    cfg, stage_dir = _prepare_case(
        workdir, transport=["energy_window: [-0.20, 0.0002, -0.16]"], ferry=False)
    config = load_config(cfg)
    assert config.transport.energy_window_unit == "Hartree"
    ctx = _Ctx(stage_dir.parent / "aims_scf", stage_dir.parent / "electrode_build",
               config)
    # 旧 Ha 口径不需要 efermi 也能写 tcontrol（aims.out 照样被 make_tcontrol 解析）
    (stage_dir.parent / "aims_scf" / "aims.out").write_text(
        "Fermi level: -0.2000 Ha\n", encoding="utf-8")
    AitranssTeStage().prepare(ctx, stage_dir)
    text = (stage_dir / "tcontrol").read_text(encoding="utf-8")
    assert "-0.2000000000" in text
    assert "-0.1600000000" in text
