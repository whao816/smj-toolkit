"""M=1 的破缺对称（BS）单重态：配置/自旋状态/判据/协议指纹。

以两个中心的反号初始磁矩寻找 BS 解。FHI-aims 的共线自旋可固定 N↑−N↓=0，
也可释放净磁矩约束；这两种输入都允许非零局域自旋密度。

这些测试锁住四条硬前提，防止 BS 被静默当成闭壳层：
1. 只允许 multiplicity=1；
2. 必须有 initial_moments，合计 0；
3. 必须至少一正一负；
4. 判据层把 BS 当开壳层（按族标定的闸门降为记录），协议指纹也随 BS 变。
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from smj.config import MoleculeConfig, load_config
from smj.context import StageContext
from smj.engines.aims import AimsEngine
from smj.errors import ConfigError
from smj.junction.criteria import evaluate_run_criteria, is_open_shell
from smj.provenance import protocol_components
from smj.spin import SpinState, record_junction_spin, spin_state
from smj.stages import get_stage
from smj.state import RunLayout, RunState

H_EV = 27.211386245988


def _config(root: Path, molecule_yaml: str):
    cfg = root / "run.yaml"
    cfg.write_text(
        "run:\n  name: t\n  mode: dry-run\n"
        "molecule:\n  xyz: molecule.xyz\n" + molecule_yaml,
        encoding="utf-8")
    return load_config(cfg)


BS_YAML = ("  multiplicity: 1\n"
           "  broken_symmetry: true\n"
           "  initial_moments:\n"
           "    7: 1.0\n"
           "    18: -1.0\n")


# --------------------------------------------------------------------------- #
# 配置解析
# --------------------------------------------------------------------------- #
def test_broken_symmetry_parses(workdir):
    config = _config(workdir, BS_YAML)
    assert config.molecule.broken_symmetry is True
    assert config.molecule.multiplicity == 1
    assert config.molecule.initial_moments == {7: 1.0, 18: -1.0}


def test_broken_symmetry_default_off(workdir):
    assert _config(workdir, "  multiplicity: 3\n").molecule.broken_symmetry is False


def test_broken_symmetry_requires_multiplicity_one(workdir):
    yaml = BS_YAML.replace("multiplicity: 1", "multiplicity: 3")
    with pytest.raises(ConfigError, match="broken_symmetry.*multiplicity=1"):
        _config(workdir, yaml)


def test_broken_symmetry_requires_initial_moments(workdir):
    yaml = ("  multiplicity: 1\n  broken_symmetry: true\n")
    with pytest.raises(ConfigError, match="需要 molecule.initial_moments"):
        _config(workdir, yaml)


def test_broken_symmetry_requires_zero_sum(workdir):
    yaml = BS_YAML.replace("18: -1.0", "18: -0.5")
    with pytest.raises(ConfigError, match="合计应为 0"):
        _config(workdir, yaml)


def test_broken_symmetry_requires_both_signs(workdir):
    # 合计为 0 但没有任何正分量：必须命中"至少一正一负"这条
    yaml = ("  multiplicity: 1\n  broken_symmetry: true\n"
            "  initial_moments:\n    7: 0.0\n    18: 0.0\n")
    with pytest.raises(ConfigError, match="至少要有一个正、一个负"):
        _config(workdir, yaml)


# --------------------------------------------------------------------------- #
# SpinState
# --------------------------------------------------------------------------- #
def test_spin_state_broken_symmetry_is_open_shell_without_constraint(workdir):
    state = spin_state(_config(workdir, BS_YAML))
    assert state.broken_symmetry is True
    assert state.open_shell is True          # 两套自旋通道
    assert state.spin == "collinear"          # FHI-aims 关键字
    assert state.fix_spin_moment is False     # BS 省略开关时沿用历史释放行为
    assert state.fixed_spin_moment == 0.0


def test_spin_state_honors_explicit_fixed_zero_for_broken_symmetry():
    state = SpinState(multiplicity=1, charge=0, fix_spin_moment=True,
                      broken_symmetry=True)
    assert state.fix_spin_moment is True
    assert state.fixed_spin_moment == 0.0


@pytest.mark.parametrize("setting,fixed", [("", False), ("null", False),
                                           ("true", True), ("false", False)])
def test_bs_config_resolves_spin_constraint(workdir, setting, fixed):
    option = "  fixed_spin_moment: %s\n" % setting if setting else ""
    config = _config(workdir, BS_YAML + option)
    assert config.molecule.fixed_spin_moment is fixed
    assert spin_state(config).fix_spin_moment is fixed


def test_spin_state_rejects_broken_symmetry_with_triplet():
    with pytest.raises(ValueError, match="multiplicity=1"):
        SpinState(multiplicity=3, charge=0, broken_symmetry=True)


def test_record_junction_spin_marks_broken_symmetry(workdir):
    summary = workdir / "input_summary.json"
    summary.write_text(json.dumps({"stage": "aims_scf"}), encoding="utf-8")
    state = SpinState(multiplicity=1, charge=0, broken_symmetry=True)
    record_junction_spin(summary, state, initial_moments={7: 1.0, 18: -1.0})
    data = json.loads(summary.read_text(encoding="utf-8"))
    assert data["broken_symmetry"] is True
    assert data["spin"] == "collinear"
    assert data["fixed_spin_moment"] is None   # 解除约束，不是 0
    assert data["initial_moments"] == {"7": 1.0, "18": -1.0}


# --------------------------------------------------------------------------- #
# aims 引擎
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("stage", ["molecule_opt", "aims_tip_relax", "aims_scf"])
@pytest.mark.parametrize("setting,fixed", [("", False), ("false", False), ("true", True)])
def test_bs_stage_inputs_preserve_local_moments_and_net_constraint(
        workdir, stage, setting, fixed):
    option = "  fixed_spin_moment: %s\n" % setting if setting else ""
    config = _config(workdir, "  broken_symmetry: true\n"
                     "  initial_moments: {1: 1.0, 2: -1.0}\n"
                     + option + "cluster:\n  species_dir: species\n")
    species = workdir / "species"
    species.mkdir()
    for element in ("S", "Au"):
        (species / (element + "_default")).write_text(
            "  species %s\n" % element, encoding="utf-8")
    xyz = workdir / "molecule.xyz"
    xyz.write_text("2\nS2\nS 0 0 0\nS 0 0 3\n", encoding="utf-8")
    geometry = workdir / "geometry.in"
    geometry.write_text("atom 0 0 0 S\natom 0 0 3 S\n", encoding="utf-8")
    ctx = StageContext(config, RunLayout(workdir / "run"), RunState("t"), dry_run=True,
                       overrides={"geometry": geometry if stage == "aims_scf" else xyz})
    stage_dir = ctx.ensure_stage_dir(stage)
    get_stage(stage).prepare(ctx, stage_dir)

    lines = (stage_dir / "control.in").read_text(encoding="utf-8").splitlines()
    control = [line.split() for line in lines]
    assert ["spin", "collinear"] in control
    assert [float(parts[1]) for parts in control if parts[:1] == ["fixed_spin_moment"]] == (
        [0.0] if fixed else [])
    geom = (stage_dir / "geometry.in").read_text(encoding="utf-8").splitlines()
    assert [float(line.split()[1]) for line in geom
            if line.split()[:1] == ["initial_moment"]] == [1.0, -1.0]
    summary = json.loads((stage_dir / "input_summary.json").read_text(encoding="utf-8"))
    key = "spin_moment" if stage == "molecule_opt" else "fixed_spin_moment"
    assert summary[key] == (0.0 if fixed else None)


def test_aims_engine_spin_and_moments_accept_broken_symmetry():
    molecule = MoleculeConfig(multiplicity=1, broken_symmetry=True,
                              initial_moments={7: 1.0, 18: -1.0})
    assert AimsEngine._spin(molecule) == "collinear"
    assert AimsEngine._initial_moments(molecule, 78) == {7: 1.0, 18: -1.0}


def test_aims_engine_still_rejects_moments_for_closed_shell():
    molecule = MoleculeConfig(multiplicity=1, initial_moments={7: 1.0})
    with pytest.raises(Exception, match="initial_moments"):
        AimsEngine._initial_moments(molecule, 78)


# --------------------------------------------------------------------------- #
# 判据 / 协议
# --------------------------------------------------------------------------- #
def test_is_open_shell_accepts_broken_symmetry():
    assert is_open_shell(1) is False
    assert is_open_shell(1, broken_symmetry=True) is True
    assert is_open_shell(None, broken_symmetry=True) is True


def _write_te(path: Path, n: int = 200, t: float = 0.25) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    efermi = -0.188452434085
    lines = ["#bias = 0.000000 V efermi = %.12f H" % efermi]
    for i in range(n):
        e_ev = -2.0 + 4.0 * i / (n - 1)
        lines.append("%.10f %.10f %.10e" % (efermi + e_ev / H_EV, e_ev, t))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_broken_symmetry_disables_family_gates(workdir):
    te = _write_te(workdir / "TE.dat")
    # "BDT" 是 SYSTEM_FAMILIES 的 system 键，映射到判据族 "Au-BDT-Au"。
    checks = {c["name"]: c for c in evaluate_run_criteria(
        te.parent, system="BDT", multiplicity=1, broken_symmetry=True)}
    for gate in ("g_in_theory_band", "seebeck_sign"):
        assert checks[gate]["applicable"] is False
    # 同一个 TE.dat、族匹配，但按闭壳层读 -> 闸门仍然适用（对照）
    closed = {c["name"]: c for c in evaluate_run_criteria(
        te.parent, system="BDT", multiplicity=1)}
    assert closed["g_in_theory_band"]["applicable"] is not False


def test_protocol_includes_broken_symmetry_only_when_enabled(workdir):
    on = protocol_components(_config(workdir, BS_YAML))
    assert on["molecule"]["broken_symmetry"] is True
    (workdir / "off").mkdir(exist_ok=True)
    off = protocol_components(_config(workdir / "off", "  multiplicity: 3\n"))
    assert "broken_symmetry" not in off["molecule"]
