"""2026-11-20 规划 §6 P1-5：``transport.extra_outputs`` 的自旋密度/Hirshfeld 证据链。

开壳层单分子结要回答"自旋是否留在分子上"，需要 SCF 给出
``output hirshfeld``（分原子磁矩）与 ``output cube spin_density``（共线自旋密度）。
这里锁住四件事：

* 白名单校验：写错的关键字在**校验期**报错，而不是让 FHI-aims 远端静默忽略；
* 默认配置逐字节不变（老 run 的 control.in 与 protocol 指纹不受影响）；
* 开壳层才允许 ``spin_density``（闭壳层没有共线自旋密度可写）；
* 请求了就必须真的产出——``spin_density.cube`` 缺失/为空、aims.out 里没有
  Hirshfeld 段，都拒绝 ``collected``。
"""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from _aims_fixtures import SCF_CONVERGED, write_scf_artifacts
from smj.config import (ConfigError, EXTRA_OUTPUTS, MoleculeConfig, TransportConfig,
                        _parse_extra_outputs)
from smj.errors import StageBlocked, StageError
from smj.junction.inputs import make_aims_scf_input
from smj.provenance import protocol_sha256, sha256_file
from smj.stages.aims_scf import AimsScfStage

_GEOMETRY = "atom 0.0 0.0 0.0 S\natom 3.0 0.0 0.0 S\n"


def _species_dir(root: Path) -> Path:
    species = root / "species"
    species.mkdir(parents=True, exist_ok=True)
    for element in ("Au", "C", "H", "S"):
        (species / ("%s_default" % element)).write_text(
            "  species %s\n    # dummy\n" % element, encoding="utf-8")
    return species


# --------------------------------------------------------------------------- #
# 配置
# --------------------------------------------------------------------------- #
def test_extra_outputs_whitelist_matches_tcontrol():
    """config 的白名单与 tcontrol 的翻译表必须同步（两份清单，一个口径）。"""
    from smj.junction.tcontrol import _EXTRA_OUTPUT_LINES

    assert set(EXTRA_OUTPUTS) == set(_EXTRA_OUTPUT_LINES)


def test_parse_extra_outputs_accepts_known_names_and_dedupes():
    assert _parse_extra_outputs(None) == ()
    assert _parse_extra_outputs(["hirshfeld"]) == ("hirshfeld",)
    assert _parse_extra_outputs(["spin_density", "hirshfeld", "hirshfeld"]) == (
        "spin_density", "hirshfeld")
    assert _parse_extra_outputs("hirshfeld") == ("hirshfeld",)


@pytest.mark.parametrize("bad", [["hirshfeld", "charge_density"], [1], {"a": 1}])
def test_parse_extra_outputs_rejects_unknown(bad):
    with pytest.raises(ConfigError):
        _parse_extra_outputs(bad)


def test_transport_config_defaults_to_no_extra_outputs():
    assert TransportConfig().extra_outputs == ()


def test_default_control_in_is_byte_identical():
    """默认配置不写任何额外的 ``output`` 行（老 run 的 control.in 保持不变）。"""
    cfg = TransportConfig()
    assert cfg.extra_outputs == ()


# --------------------------------------------------------------------------- #
# control.in 渲染
# --------------------------------------------------------------------------- #
def test_make_aims_scf_input_emits_extra_output_lines(workdir: Path):
    geometry = workdir / "geometry.in"
    geometry.write_text(_GEOMETRY, encoding="utf-8")
    out = workdir / "aims_scf"
    out.mkdir()
    make_aims_scf_input(
        geometry, out, species_dir=_species_dir(workdir),
        extra_outputs=("hirshfeld", "spin_density"),
    )
    control = (out / "control.in").read_text(encoding="utf-8")
    assert "  output hirshfeld\n" in control
    assert "  output cube spin_density\n" in control
    # 仍保留 aitranss/mulliken 与 KS_method serial（aitranss 依赖的布局）
    assert "  output aitranss\n" in control
    assert "  output mulliken\n" in control
    assert "  KS_method serial\n" in control


def test_make_aims_scf_input_default_has_no_extra_output_lines(workdir: Path):
    geometry = workdir / "geometry.in"
    geometry.write_text(_GEOMETRY, encoding="utf-8")
    out = workdir / "aims_scf"
    out.mkdir()
    make_aims_scf_input(geometry, out, species_dir=_species_dir(workdir))
    control = (out / "control.in").read_text(encoding="utf-8")
    assert "output hirshfeld" not in control
    assert "output cube" not in control


def test_extra_outputs_do_not_change_protocol_hash(workdir: Path):
    """额外输出只决定多存哪些文件，不改求解参数 -> 指纹不变。"""
    from smj.config import RunConfig, RunInfo

    config = RunConfig(run=RunInfo(name="t"))
    before = protocol_sha256(config)
    config.transport = TransportConfig(extra_outputs=("hirshfeld", "spin_density"))
    assert protocol_sha256(config) == before


# --------------------------------------------------------------------------- #
# prepare：闭壳层不能要自旋密度
# --------------------------------------------------------------------------- #
class _PrepareCtx:
    def __init__(self, run_dir: Path, *, multiplicity: int = 1,
                 extra_outputs: tuple[str, ...] = ()) -> None:
        self.layout = SimpleNamespace(run_dir=Path(run_dir))
        self.config = SimpleNamespace(
            molecule=MoleculeConfig(multiplicity=multiplicity),
            transport=TransportConfig(extra_outputs=extra_outputs),
            electrode=SimpleNamespace(vdw=False),
            cluster=SimpleNamespace(species_dir=None),
        )
        self._species = _species_dir(Path(run_dir).parent)

    def override_path(self, key: str):
        return None

    def upstream(self, stage: str, name: str, override=None) -> Path:
        return Path(self.layout.run_dir) / "electrode_build" / "geometry.in"

    def artifact(self, stage: str, name: str) -> Path:
        return Path(self.layout.run_dir) / stage / name

    def species_dir(self):
        return self._species

    def estimate(self, name: str, natoms: int) -> dict:
        return {"atoms": natoms}

    def stage_env(self, name: str, estimate: dict) -> dict:
        return {}

    def log(self, message: str) -> None:
        pass


def _seed_electrode_geometry(run_dir: Path) -> None:
    stage = run_dir / "electrode_build"
    stage.mkdir(parents=True, exist_ok=True)
    geometry = stage / "geometry.in"
    geometry.write_text(_GEOMETRY, encoding="utf-8")
    # 2026-11-22 起 aims_scf.prepare 要求父几何与电极 metadata 的指纹一致
    # （缺 metadata 直接 StageBlocked）；测试也要按真实链路给出确认过的指纹。
    (stage / "electrode_metadata.json").write_text(
        json.dumps({"surface_fingerprint": sha256_file(geometry)}), encoding="utf-8")


def test_prepare_blocks_when_parent_metadata_missing(workdir: Path):
    """没有 electrode_metadata 就不知道 SCF 几何是否还是 A2 确认的那份。"""
    run_dir = workdir / "run"
    stage = run_dir / "electrode_build"
    stage.mkdir(parents=True)
    (stage / "geometry.in").write_text(_GEOMETRY, encoding="utf-8")
    stage_dir = run_dir / "aims_scf"
    stage_dir.mkdir()

    with pytest.raises(StageBlocked, match="electrode_metadata"):
        AimsScfStage().prepare(_PrepareCtx(run_dir), stage_dir)


def test_prepare_explicit_geometry_keeps_standalone_mode_visible(
        workdir: Path, monkeypatch):
    """--geometry 是唯一允许跳过父指纹的模式，但必须在 notes 里写清楚。"""
    run_dir = workdir / "run"
    stage = run_dir / "electrode_build"
    stage.mkdir(parents=True)
    geometry = stage / "geometry.in"
    geometry.write_text(_GEOMETRY, encoding="utf-8")
    stage_dir = run_dir / "aims_scf"
    stage_dir.mkdir()

    class _OverrideCtx(_PrepareCtx):
        def override_path(self, key: str):
            return geometry if key == "geometry" else None

    def _fake_scf_input(geometry_in, outdir, **kwargs):
        outdir = Path(outdir)
        control = outdir / "control.in"
        control.write_text("xc pbe\n", encoding="utf-8")
        script = outdir / "sub_aims.sh"
        script.write_text("#!/bin/bash\n", encoding="utf-8")
        summary = outdir / "input_summary.json"
        summary.write_text(json.dumps({"spin": kwargs.get("spin")}),
                           encoding="utf-8")
        return {"geometry": str(geometry_in), "control": str(control),
                "submit_script": str(script), "summary": json.loads(
                    summary.read_text(encoding="utf-8"))}

    monkeypatch.setattr("smj.stages.aims_scf.make_aims_scf_input",
                        _fake_scf_input)
    prepared = AimsScfStage().prepare(_OverrideCtx(run_dir), stage_dir)

    assert any("独立运行模式" in note for note in prepared.notes)


def test_prepare_rejects_spin_density_for_closed_shell(workdir: Path):
    run_dir = workdir / "run"
    _seed_electrode_geometry(run_dir)
    stage_dir = run_dir / "aims_scf"
    stage_dir.mkdir()
    ctx = _PrepareCtx(run_dir, multiplicity=1, extra_outputs=("spin_density",))
    with pytest.raises(StageError, match="闭壳层"):
        AimsScfStage().prepare(ctx, stage_dir)


def test_prepare_writes_extra_outputs_for_open_shell(workdir: Path):
    run_dir = workdir / "run"
    _seed_electrode_geometry(run_dir)
    stage_dir = run_dir / "aims_scf"
    stage_dir.mkdir()
    ctx = _PrepareCtx(run_dir, multiplicity=2, extra_outputs=("spin_density",))
    prepared = AimsScfStage().prepare(ctx, stage_dir)
    control = (stage_dir / "control.in").read_text(encoding="utf-8")
    assert "  output cube spin_density\n" in control
    assert "spin_density" in " ".join(prepared.notes)


# --------------------------------------------------------------------------- #
# collect：请求了就必须真的产出
# --------------------------------------------------------------------------- #
class _CollectCtx:
    def __init__(self, run_dir: Path, *, extra_outputs: tuple[str, ...] = ()) -> None:
        self.layout = SimpleNamespace(run_dir=Path(run_dir))
        self.config = SimpleNamespace(transport=TransportConfig(extra_outputs=extra_outputs))

    def require_stage_file(self, stage: str, name: str) -> Path:
        path = self.layout.run_dir / stage / name
        if not path.is_file():
            raise StageError("缺 %s" % path)
        return path


def _seed_stage(stage_dir: Path, *, log_text: str | None = None) -> None:
    stage_dir.mkdir(parents=True, exist_ok=True)
    (stage_dir / "geometry.in").write_text(_GEOMETRY, encoding="utf-8")
    (stage_dir / "aims.out").write_text(
        log_text if log_text is not None else
        "  FHI-aims version      : 221103\n" + SCF_CONVERGED,
        encoding="utf-8")
    write_scf_artifacts(stage_dir)


def test_collect_rejects_missing_spin_density_cube(workdir: Path):
    run_dir = workdir / "run"
    stage_dir = run_dir / "aims_scf"
    _seed_stage(stage_dir)
    ctx = _CollectCtx(run_dir, extra_outputs=("spin_density",))
    with pytest.raises(StageError, match="spin_density.cube"):
        AimsScfStage().collect(ctx, stage_dir)


def test_collect_rejects_empty_spin_density_cube(workdir: Path):
    run_dir = workdir / "run"
    stage_dir = run_dir / "aims_scf"
    _seed_stage(stage_dir)
    (stage_dir / "spin_density.cube").write_text("", encoding="utf-8")
    ctx = _CollectCtx(run_dir, extra_outputs=("spin_density",))
    with pytest.raises(StageError, match="spin_density.cube"):
        AimsScfStage().collect(ctx, stage_dir)


def test_collect_records_spin_density_cube(workdir: Path):
    run_dir = workdir / "run"
    stage_dir = run_dir / "aims_scf"
    _seed_stage(stage_dir)
    (stage_dir / "spin_density.cube").write_text(
        "cube header\ndata\n", encoding="utf-8")
    ctx = _CollectCtx(run_dir, extra_outputs=("spin_density",))
    out = AimsScfStage().collect(ctx, stage_dir)
    payload = out.results["aims_scf"]
    assert payload["extra_outputs"] == ["spin_density"]
    assert payload["spin_density_cube_bytes"] == \
        (stage_dir / "spin_density.cube").stat().st_size
    assert "spin_density.cube" in out.artifacts


def test_collect_rejects_hirshfeld_without_section(workdir: Path):
    """aims.out 只有 vdw 关键字回显、没有 Hirshfeld 分析段 -> 关键字没被接受。"""
    run_dir = workdir / "run"
    stage_dir = run_dir / "aims_scf"
    _seed_stage(stage_dir, log_text=(
        "  FHI-aims version      : 221103\n"
        "  vdw_correction_hirshfeld\n" + SCF_CONVERGED))
    ctx = _CollectCtx(run_dir, extra_outputs=("hirshfeld",))
    with pytest.raises(StageError, match="Hirshfeld"):
        AimsScfStage().collect(ctx, stage_dir)


def test_collect_accepts_hirshfeld_section(workdir: Path):
    run_dir = workdir / "run"
    stage_dir = run_dir / "aims_scf"
    _seed_stage(stage_dir, log_text=(
        "  FHI-aims version      : 221103\n"
        "  Performing Hirshfeld analysis of fragment charges and moments.\n"
        "  ----------------------------------------------------------------------\n"
        "  | Atom     1: C\n"
        "  |   Hirshfeld charge        :     -0.06377126\n"
        "  |   Hirshfeld spin moment   :      0.00010109\n" + SCF_CONVERGED))
    ctx = _CollectCtx(run_dir, extra_outputs=("hirshfeld",))
    out = AimsScfStage().collect(ctx, stage_dir)
    payload = out.results["aims_scf"]
    assert payload["hirshfeld_mentions"] == 3
    # 2026-11-21 自由基支线 S1：Hirshfeld 段被解析成自旋分区侧车。
    assert payload["spin_partition"]["atoms"] == 1
    assert payload["spin_partition"]["source"] == "aims.out"
    assert "spin_partition.json" in out.artifacts


def test_collect_rejects_hirshfeld_section_without_atoms(workdir: Path):
    """段名在、但解析不出 ``| Atom N:`` 行 -> 没有证据，不标 collected。"""
    run_dir = workdir / "run"
    stage_dir = run_dir / "aims_scf"
    _seed_stage(stage_dir, log_text=(
        "  FHI-aims version      : 221103\n"
        "  Performing Hirshfeld analysis of fragment charges and moments.\n"
        "  Hirshfeld charge      0.1234\n" + SCF_CONVERGED))
    ctx = _CollectCtx(run_dir, extra_outputs=("hirshfeld",))
    with pytest.raises(StageError, match="Atom"):
        AimsScfStage().collect(ctx, stage_dir)


def test_collect_without_extra_outputs_is_unchanged(workdir: Path):
    run_dir = workdir / "run"
    stage_dir = run_dir / "aims_scf"
    _seed_stage(stage_dir)
    out = AimsScfStage().collect(_CollectCtx(run_dir), stage_dir)
    payload = out.results["aims_scf"]
    assert payload["extra_outputs"] == []
    assert "spin_density_cube_bytes" not in payload
    assert "hirshfeld_mentions" not in payload
