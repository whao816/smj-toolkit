"""molecule_opt 的引擎插件层与注册表。

``molecule_opt`` 只保留与引擎无关的公共流程：电子数/多重度校验、把插件的
prepare 结果补齐为 :class:`smj.stages.base.Prepared`、按注册表 dispatch。
引擎插件自己声明：

* ``script`` / ``outputs`` —— 作业脚本名与其产物契约；
* ``resource_anchor`` / ``env_anchor`` —— 资源估算、作业 env 用哪个阶段的标定
  （Gaussian 用 ``molecule_opt``，aims 借 ``aims_tip_relax``）；
* ``probes`` —— preflight 除 FHI-aims / aitranss 之外还要探测的软件 id；
* ``skip_probe_note`` —— 不需要额外探测时给用户的说明；
* ``remote_files`` —— slurm-ssh 回取清单（``None`` = 沿用 vendored 计划）。

新增引擎 = 在本包里注册一个 :class:`Engine` 子类，必要时在
:mod:`smj.preflight` 的 ``_ENGINE_PROBES`` 登记探测函数；不再改
``config`` / ``cli`` / ``preflight`` / ``runner`` 里的 if-else。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, ClassVar

from ..errors import ConfigError

#: 缺省引擎（与 ``MoleculeConfig.opt_engine`` 的默认值一致）
DEFAULT_ENGINE = "aims"


@dataclass
class EnginePrepared:
    """引擎 prepare 的产物（阶段再补资源/命令等公共字段）。"""

    inputs: list[Path] = field(default_factory=list)
    artifacts: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    script: str = ""


@dataclass
class EngineCollected:
    """引擎 collect 的产物（阶段再补 stage/output_dir）。"""

    artifacts: list[str] = field(default_factory=list)
    results: dict[str, Any] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)


class Engine:
    """molecule_opt 引擎插件协议。"""

    name: ClassVar[str] = ""
    script: ClassVar[str] = ""
    outputs: ClassVar[tuple[str, ...]] = ()
    resource_anchor: ClassVar[str] = ""
    env_anchor: ClassVar[str] = ""
    probes: ClassVar[tuple[str, ...]] = ()
    skip_probe_note: ClassVar[str | None] = None
    remote_files: ClassVar[tuple[str, ...] | None] = None

    def prepare(self, ctx, stage_dir: Path, *, stage: str, molecule: Any,
                xyz: Path, natoms: int, electron_note: str,
                estimate: dict, env: dict) -> EnginePrepared:
        raise NotImplementedError

    def collect(self, ctx, stage_dir: Path, *, stage: str) -> EngineCollected:
        raise NotImplementedError


_REGISTRY: dict[str, Engine] = {}


def register(engine_cls: type[Engine]) -> type[Engine]:
    """类装饰器：按小写引擎名登记插件（重名直接报错）。"""
    if not engine_cls.name:
        raise ValueError("引擎类 %s 未声明 name" % engine_cls.__name__)
    key = engine_cls.name.lower()
    if key in _REGISTRY:
        raise ValueError("引擎名重复注册: %s" % key)
    _REGISTRY[key] = engine_cls()
    return engine_cls


def engine_names() -> tuple[str, ...]:
    """已注册引擎名（注册顺序，供 config 白名单与报错文案使用）。"""
    return tuple(_REGISTRY)


def all_engines() -> list[Engine]:
    return list(_REGISTRY.values())


def get_engine(name: str) -> Engine:
    key = str(name or "").strip().lower()
    if key not in _REGISTRY:
        raise ConfigError(
            "molecule.opt_engine 应为 %s，实际是 %r"
            % (" 或 ".join(engine_names()), name))
    return _REGISTRY[key]


def engine_for_config(config: Any) -> Engine:
    """从 RunConfig（或同形对象）的 ``molecule.opt_engine`` 取插件。"""
    molecule = getattr(config, "molecule", None)
    name = getattr(molecule, "opt_engine", DEFAULT_ENGINE) or DEFAULT_ENGINE
    return get_engine(name)


# 注册内置引擎（顺序即 engine_names 的顺序）。
from . import aims  # noqa: E402,F401
from . import gaussian  # noqa: E402,F401

__all__ = [
    "DEFAULT_ENGINE",
    "Engine",
    "EngineCollected",
    "EnginePrepared",
    "all_engines",
    "engine_for_config",
    "engine_names",
    "get_engine",
    "register",
]
