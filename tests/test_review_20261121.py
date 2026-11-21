"""2026-11-21 第三轮复审（master 6384ce0d）—— 附件之外的补充断言。

附件 8 条已逐条入库到 ``tests/test_review_20261121_regressions.py``（文件头有
"附件测试 → 修复项"映射）。这里放评审之外、但同一条修复链上必须钉住的补充
行为：

- R2 的轴键校验：拼错/不存在的 ``section.key`` 必须报错，而不是静默忽略
  （架空的轴就是"看起来扫了、其实没扫"）；
- R2 的 ``paths.species_dir`` 兼容别名：单给时并入 ``species_defaults`` 并告警，
  与 ``species_defaults`` 同时给时报错；
- R6 的开壳层配对：末两张表分别是 α/β（行数可以不同），占据表行数与对应
  通道的波函数态数不符时抛 ``PdosError``，不按索引静默截断。
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

from smj.config import ConfigError, load_config
from smj.pdos import PdosError, analyze_pdos
from smj.runner import StageRunner

from test_pdos import _aims_out, _run as _pdos_run, _uhfmo

ROOT = Path(__file__).resolve().parents[1]


def _scan():
    spec = importlib.util.spec_from_file_location(
        "scan_convergence", ROOT / "scripts" / "scan_convergence.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _base_config(workdir: Path) -> Path:
    target = workdir / "smj.yaml"
    target.write_text(
        "run:\n  name: base\n  mode: local\n"
        "molecule:\n  xyz: mol.xyz\n"
        "electrode:\n  nlayers: 3\n  au_per_electrode: 20\n",
        encoding="utf-8")
    return target


# --------------------------------------------------------------------------- #
# R2：轴键校验与 paths.species_dir 兼容别名
# --------------------------------------------------------------------------- #
def test_plan_rejects_unknown_axis_key(workdir):
    scan = _scan()
    with pytest.raises(scan.ScanError, match="electrode.typox"):
        scan.build_plan(_base_config(workdir), [("electrode.typox", [1])], workdir / "scan")


def test_paths_species_dir_alias_and_conflict(workdir):
    light = (workdir / "light").resolve()
    tight = (workdir / "tight").resolve()
    alias = workdir / "alias.yaml"
    alias.write_text(
        "run:\n  name: r\n  mode: local\nmolecule:\n  xyz: mol.xyz\n"
        "paths:\n  species_dir: " + light.as_posix() + "\n", encoding="utf-8")
    config = load_config(alias)
    assert config.paths.species_defaults is not None
    assert Path(config.paths.species_defaults) == light
    assert any("species_dir 已并入" in warning for warning in config.warnings), config.warnings
    assert not any("未知字段" in warning for warning in config.warnings), config.warnings

    clash = workdir / "clash.yaml"
    clash.write_text(
        "run:\n  name: r\n  mode: local\nmolecule:\n  xyz: mol.xyz\n"
        "paths:\n  species_dir: " + light.as_posix() + "\n"
        "  species_defaults: " + tight.as_posix() + "\n", encoding="utf-8")
    with pytest.raises(ConfigError, match="species_dir"):
        load_config(clash)


# --------------------------------------------------------------------------- #
# R6：开壳层配对与态数校验（附件第 7 条只管闭壳层）
# --------------------------------------------------------------------------- #
def test_open_shell_pairs_last_two_tables_per_channel(workdir):
    """开壳层：末两张表（α、β 行数可以不同）各自对本通道，早期迭代表忽略。"""
    run = _pdos_run(workdir, natoms=2, mos=_uhfmo("alpha", (-0.2, 0.2), natoms=2))
    scf = run / "aims_scf"
    (scf / "alpha.aims").write_text(_uhfmo("alpha", (-0.2, 0.2), natoms=2), encoding="utf-8")
    (scf / "beta.aims").write_text(_uhfmo("beta", (-0.1,), natoms=2), encoding="utf-8")
    (scf / "aims.out").write_text(
        _aims_out([1.0, 0.0], [1.0]) + _aims_out([2.0, 0.0], [0.0]),
        encoding="utf-8")

    result = analyze_pdos(run, groups="mol=1-2", sigma_ev=0.01)
    alpha = result["channels"]["alpha"]
    beta = result["channels"]["beta"]
    assert alpha["states"] == 2 and beta["states"] == 1
    assert sum(alpha["integrated_weight"].values()) == pytest.approx(2.0, abs=1e-9)
    assert sum(beta["integrated_weight"].values()) == pytest.approx(0.0, abs=1e-9)


def test_occupation_table_must_match_state_count(workdir):
    """闭壳层占据表行数与波函数态数不符 -> PdosError（不按索引截断）。"""
    run = _pdos_run(workdir, natoms=2, offsets=(-0.2, 0.2), grid=(-0.2, 0.0, 0.2))
    (run / "aims_scf" / "aims.out").write_text(
        _aims_out([2.0, 0.0, 0.0]), encoding="utf-8")
    with pytest.raises(PdosError, match="占据"):
        analyze_pdos(run, groups="mol=1-2", sigma_ev=0.01)


def test_open_shell_table_count_mismatch_raises(workdir):
    """开壳层占据表行数与通道态数不符 -> PdosError（不按索引截断）。"""
    run = _pdos_run(workdir, natoms=2, mos=_uhfmo("alpha", (-0.2, 0.2), natoms=2))
    scf = run / "aims_scf"
    (scf / "alpha.aims").write_text(_uhfmo("alpha", (-0.2, 0.2), natoms=2), encoding="utf-8")
    (scf / "beta.aims").write_text(_uhfmo("beta", (-0.1, 0.1), natoms=2), encoding="utf-8")
    (scf / "aims.out").write_text(
        _aims_out([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]), encoding="utf-8")
    with pytest.raises(PdosError, match="占据"):
        analyze_pdos(run, groups="mol=1-2", sigma_ev=0.01)
