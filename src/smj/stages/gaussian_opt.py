"""兼容旧名字：``gaussian_opt`` -> :mod:`smj.stages.molecule_opt`。

阶段规范名在 2026-11-13 改为 ``molecule_opt``（引擎由 ``molecule.opt_engine``
选择）。本模块只保留旧 import 路径，不再单独注册阶段。
"""

from __future__ import annotations

from .molecule_opt import MoleculeOptStage

#: 旧类名别名（第三方脚本可能直接 import）。
GaussianOptStage = MoleculeOptStage

__all__ = ["GaussianOptStage", "MoleculeOptStage"]
