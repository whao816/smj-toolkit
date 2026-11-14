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
    assert config.transport.energy_window == (-0.21, 0.0002, -0.17)
    assert config.transport.resolved_window() == (-0.21, 0.0002, -0.17)


def test_energy_window_explicit_override_kept(workdir):
    config = load_config(_project(
        workdir, transport=["energy_window: [-0.22, 0.0001, -0.16]"]))
    assert config.transport.resolved_window() == (-0.22, 0.0001, -0.16)


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
