"""0.2.0rc1 发布默认值回归：energy_window / au_per_electrode / preflight 探测。"""
from __future__ import annotations

import types
from pathlib import Path

from smj.config import load_config
from smj.preflight import _needs_gaussian_check


def _project(workdir: Path, *, electrode: list[str] | None = None,
             transport: list[str] | None = None) -> Path:
    project = (workdir / "proj").resolve()
    project.mkdir(parents=True)
    (project / "molecule.xyz").write_text("2\n\nC 0 0 0\nH 1 0 0\n", encoding="utf-8")
    lines = [
        "run:", "  name: t", "  mode: dry-run", "  dir: %s" % project.as_posix(),
        "molecule:", "  xyz: %s" % (project / "molecule.xyz").as_posix(),
    ]
    if electrode is not None:
        lines.append("electrode:")
        lines.extend("  " + row for row in electrode)
    if transport is not None:
        lines.append("transport:")
        lines.extend("  " + row for row in transport)
    cfg = project / "smj.yaml"
    cfg.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return cfg


def test_energy_window_default_around_fermi(workdir):
    config = load_config(_project(workdir))
    # 默认：eV 且相对 E_F（E_F ±2 eV），步长 5 meV。
    assert config.transport.energy_window == (-2.0, 0.005, 2.0)
    assert config.transport.energy_window_unit == "eV"
    assert config.transport.relative_window() == (-2.0, 0.005, 2.0)
    # 换算成 aitranss 的 Ha 绝对值必须显式给 efermi（不猜）。
    try:
        config.transport.resolved_window()
    except ValueError as exc:
        assert "efermi" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("eV 口径缺 efermi 时应报错")
    efermi = -0.180814064867
    ener, estep, eend = config.transport.resolved_window(efermi)
    assert abs(ener - (efermi - 2.0 / 27.211386245988)) < 1e-12
    assert abs(eend - (efermi + 2.0 / 27.211386245988)) < 1e-12
    assert abs(estep - 0.005 / 27.211386245988) < 1e-15


def test_energy_window_explicit_override_kept(workdir):
    # 旧口径（Ha 绝对）三元组：按量级识别并原样生效，附迁移 warning。
    config = load_config(_project(
        workdir, transport=["energy_window: [-0.22, 0.0001, -0.16]"]))
    assert config.transport.energy_window_unit == "Hartree"
    assert config.transport.resolved_window() == (-0.22, 0.0001, -0.16)
    assert any("旧口径" in w for w in config.warnings)


def test_library_ignores_au_per_electrode(workdir):
    config = load_config(_project(
        workdir, electrode=["source: library", "au_per_electrode: 12"]))
    assert config.electrode.au_per_electrode is None
    assert any("au_per_electrode" in w for w in config.warnings)


def test_builtin_keeps_au_per_electrode(workdir):
    config = load_config(_project(
        workdir, electrode=["source: builtin", "au_per_electrode: 12"]))
    assert config.electrode.au_per_electrode == 12
    assert not any("au_per_electrode" in w for w in config.warnings)


def test_preflight_gaussian_probe_only_for_gaussian_engine():
    aims = types.SimpleNamespace(molecule=types.SimpleNamespace(opt_engine="aims"))
    gauss = types.SimpleNamespace(molecule=types.SimpleNamespace(opt_engine="gaussian"))
    assert _needs_gaussian_check(aims) is False
    assert _needs_gaussian_check(gauss) is True
