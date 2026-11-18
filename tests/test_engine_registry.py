"""引擎注册表：config 白名单、资源锚点、预检探测都由插件声明派生。"""

from types import SimpleNamespace

import pytest

from smj.config import _parse_molecule
from smj.engines import (
    DEFAULT_ENGINE,
    all_engines,
    engine_for_config,
    engine_names,
    get_engine,
)
from smj.errors import ConfigError
from smj.preflight import _needs_gaussian_check


def test_engine_names_keep_registration_order():
    assert engine_names() == ("aims", "gaussian")


def test_every_engine_declares_its_contract():
    for engine in all_engines():
        assert engine.name and engine.script
        assert engine.outputs
        assert engine.resource_anchor and engine.env_anchor


def test_get_engine_is_case_insensitive():
    assert get_engine("AIMS").name == "aims"
    assert get_engine(" Gaussian ").name == "gaussian"


def test_unknown_engine_raises_config_error():
    with pytest.raises(ConfigError, match="aims 或 gaussian"):
        get_engine("qe")


def test_engine_for_config_defaults_to_default_engine():
    assert engine_for_config(
        SimpleNamespace(molecule=SimpleNamespace())).name == DEFAULT_ENGINE
    assert engine_for_config(
        SimpleNamespace(molecule=SimpleNamespace(opt_engine="gaussian"))).name == "gaussian"


def test_config_whitelist_comes_from_registry():
    # 报错文案里的可选值来自 engine_names()，不是 config 里的硬编码元组。
    with pytest.raises(ConfigError, match="aims 或 gaussian"):
        _parse_molecule({"opt_engine": "qe"}, [])


def test_gaussian_probe_follows_engine_declaration():
    aims = SimpleNamespace(molecule=SimpleNamespace(opt_engine="aims"))
    gaussian = SimpleNamespace(molecule=SimpleNamespace(opt_engine="gaussian"))
    assert _needs_gaussian_check(aims) is False
    assert _needs_gaussian_check(gaussian) is True


def test_remote_files_are_declared_per_engine():
    assert "geometry.in.next_step" in get_engine("aims").remote_files
    # None = 沿用 vendored 的 gaussian_opt 回取计划（mol.gjf/mol.log/sub_g16.sh）。
    assert get_engine("gaussian").remote_files is None
