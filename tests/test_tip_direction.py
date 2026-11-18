"""``electrode.tip_direction``：axis（默认，旧行为）/ lone_pair（沿锚原子孤对外推）。

``lone_pair`` 解决的是间位取代 SMe 分子里 S···S 连线与 S–C(Me) 近共线的问题：
沿 axis 放 tip 会把金属压进甲基，此时应沿 C–S–C 角平分线的反向放置。

本组测试只锁几何方向与配置校验，不跑 prepare/collect。
"""
from __future__ import annotations

from pathlib import Path

import pytest

from smj.config import load_config
from smj.errors import ConfigError, StageError
from smj.stages.aims_tip_relax import AimsTipRelaxStage

# 两个合成的 S(C)(C) 中心，相距 10 Å；每个 S 恰有 2 个 1.5 Å 近邻，夹角 90°。
# 孤对方向 = -(C1-S + C2-S) 的单位向量 = (-0.7071, -0.7071, 0)。
_SYM_XYZ = """6
synthetic C-S-C x2
S   0.000  0.000  0.000
C   1.500  0.000  0.000
C   0.000  1.500  0.000
S  10.000  0.000  0.000
C  11.500  0.000  0.000
C  10.000  1.500  0.000
"""

# 左侧 S 只有 1 个近邻：lone_pair 无法定义角平分线，必须 fail-closed。
_ONE_NEIGHBOR_XYZ = """2
one-neighbor S
S  0.000  0.000  0.000
C  1.500  0.000  0.000
"""


def _config_yaml(root: Path, value: str) -> Path:
    cfg = root / "run.yaml"
    cfg.write_text(
        "run:\n  name: t\n  mode: dry-run\n"
        "molecule:\n  xyz: molecule.xyz\n"
        "electrode:\n  tip_direction: %s\n" % value,
        encoding="utf-8")
    return cfg


def _xyz(root: Path, text: str) -> Path:
    path = root / "molecule.xyz"
    path.write_text(text, encoding="utf-8")
    return path


def test_tip_direction_defaults_to_axis(workdir):
    cfg = workdir / "run.yaml"
    cfg.write_text("run:\n  name: t\n  mode: dry-run\n"
                   "molecule:\n  xyz: molecule.xyz\n", encoding="utf-8")
    assert load_config(cfg).electrode.tip_direction == "axis"


def test_tip_direction_accepts_lone_pair_and_normalizes_case(workdir):
    assert load_config(_config_yaml(workdir, "LONE_PAIR")).electrode.tip_direction == "lone_pair"


def test_tip_direction_rejects_unknown_value(workdir):
    with pytest.raises(ConfigError, match="electrode.tip_direction"):
        load_config(_config_yaml(workdir, "sideways"))


def test_axis_returns_none_directions(workdir):
    xyz = _xyz(workdir, _SYM_XYZ)
    assert AimsTipRelaxStage._tip_directions(xyz, 1, 4, "axis") == (None, None)


def test_lone_pair_points_opposite_the_two_bonds(workdir):
    xyz = _xyz(workdir, _SYM_XYZ)
    left, right = AimsTipRelaxStage._tip_directions(xyz, 1, 4, "lone_pair")
    for vec in (left, right):
        assert vec == pytest.approx((-0.70710678, -0.70710678, 0.0), abs=1e-6)
        assert sum(component * component for component in vec) == pytest.approx(1.0)


def test_lone_pair_requires_exactly_two_neighbors(workdir):
    xyz = _xyz(workdir, _ONE_NEIGHBOR_XYZ)
    with pytest.raises(StageError, match="恰有 2 个近邻"):
        AimsTipRelaxStage._tip_directions(xyz, 1, 1, "lone_pair")


def test_unknown_mode_fails_closed(workdir):
    xyz = _xyz(workdir, _SYM_XYZ)
    with pytest.raises(StageError, match="tip_direction"):
        AimsTipRelaxStage._tip_directions(xyz, 1, 4, "sideways")
