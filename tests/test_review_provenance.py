"""协议指纹：schema 形状演进不能算物理漂移，物理参数变化仍必须阻断。

2026-11-20 ttf-release 复核实测：老 run 的 ``protocol_sha256`` 与当前配置不同，
但两侧都有的键（温度、能量窗、锚定、tip 距离、电荷、多重度）取值完全相同，差异
全来自后来新增/改名的记录键。把这种差异判成 ``provenance_drift`` error，会让每次
代码升级都把历史 run 打成"不可按原参数解读"。
"""
from __future__ import annotations

import json
from pathlib import Path

from smj.config import (
    ClusterConfig,
    ElectrodeConfig,
    MoleculeConfig,
    RunConfig,
    RunInfo,
    TransportConfig,
)
from smj.provenance import compare_protocol_configs, protocol_components
from smj.review import review_run
from smj.runner import StageRunner
from smj.stages.base import Collected, Prepared


class _FakeStage:
    name = "gaussian_opt"
    script = "sub_g16.sh"
    outputs = ("out.txt",)
    needs_submit = True

    def prepare(self, ctx, stage_dir):
        stage_dir = Path(stage_dir)
        (stage_dir / "input.txt").write_text("hello\n", encoding="utf-8")
        return Prepared(stage=self.name, inputs=[stage_dir / "input.txt"],
                        artifacts=["out.txt"], resources={},
                        script=self.script)

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


def _prepared_run(workdir: Path):
    config = _config(workdir)
    runner = StageRunner(
        config, logger=lambda message: None, sleep=lambda seconds: None,
        stage_getter=lambda name: _FakeStage())
    runner.prepare("gaussian_opt")
    return config, runner.run_dir


def _rewrite_spec(run_dir: Path, mutate) -> None:
    path = run_dir / "run_spec.json"
    spec = json.loads(path.read_text(encoding="utf-8"))
    # 伪造"老 run 的哈希"：内容一律对不上，差异只能由逐键比较解释。
    spec["protocol_sha256"] = "0" * 64
    mutate(spec)
    path.write_text(json.dumps(spec, ensure_ascii=False, indent=2),
                    encoding="utf-8")


def _codes(config, run_dir):
    return review_run(config, force=True,
                      review_config={"connection_check": False}).codes()


def test_schema_shape_evolution_is_warning_not_drift(workdir):
    config, run_dir = _prepared_run(workdir)

    def mutate(spec):
        spec.pop("protocol_schema", None)                    # 老 run 没有这个键
        spec["config"]["molecule"]["xyz"] = "/old/path.xyz"  # 旧键名
        spec["config"]["electrode"]["tip_direction"] = "axis"

    _rewrite_spec(run_dir, mutate)
    codes = _codes(config, run_dir)

    assert "provenance_schema_evolved" in codes
    assert "provenance_drift" not in codes


def test_changed_physics_parameter_is_still_drift(workdir):
    config, run_dir = _prepared_run(workdir)

    def mutate(spec):
        spec["config"]["transport"]["temperature"] = 250.0

    _rewrite_spec(run_dir, mutate)

    assert "provenance_drift" in _codes(config, run_dir)


def test_missing_critical_key_is_still_drift(workdir):
    config, run_dir = _prepared_run(workdir)

    def mutate(spec):
        spec["config"]["transport"].pop("temperature")

    _rewrite_spec(run_dir, mutate)

    assert "provenance_drift" in _codes(config, run_dir)


def test_compare_protocol_configs_separates_shape_from_physics(workdir):
    config = _config(workdir)
    current = protocol_components(config)

    same = {key: dict(value) for key, value in current.items()}
    same["molecule"]["legacy_only"] = 1
    compared = compare_protocol_configs(same, current)
    assert compared["differences"] == []
    assert compared["critical_missing"] == []
    assert compared["only_run"] == ["molecule.legacy_only"]

    changed = {key: dict(value) for key, value in current.items()}
    changed["transport"]["temperature"] = 42.0
    compared = compare_protocol_configs(changed, current)
    assert [item["path"] for item in compared["differences"]] == [
        "transport.temperature"]
