"""打包契约回归（2026-11-22 复审批 4）。

* 构建脚本在 Python 3.10（没有 tomllib）上也要能读到版本号；
* wheel METADATA 的 ``Provides-Extra: all`` 必须与 pyproject 的 extras 对得上；
* sdist 要带 README 链接到的 docs/ 与 demo 数值基线。
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import smj


def _build_wheel_module():
    path = Path(__file__).resolve().parents[1] / "scripts" / "build_wheel.py"
    spec = importlib.util.spec_from_file_location("smj_build_wheel", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_version_reads_pyproject_and_matches_package():
    bw = _build_wheel_module()
    assert bw._version() == smj.__version__


def test_metadata_extra_all_matches_pyproject():
    bw = _build_wheel_module()
    text = bw._metadata(bw._version())
    extras = bw._extras()
    assert set(extras) == {"plot", "siesta", "dev", "all"}
    for name in extras:
        assert "Provides-Extra: %s" % name in text
    assert 'Requires-Dist: netCDF4>=1.6; extra == "siesta"' in text
    assert 'Requires-Dist: matplotlib>=3.7; extra == "dev"' in text
    pyproject = (Path(__file__).resolve().parents[1] / "pyproject.toml").read_text(
        encoding="utf-8")
    assert "siesta = [" in pyproject and "dev = [" in pyproject


def test_sdist_includes_docs_and_demo_reference():
    bw = _build_wheel_module()
    names = {name for name, _ in bw._sdist_files(bw._version())}
    assert "docs/plot_style.md" in names
    assert "demo/demo-reference.json" in names
