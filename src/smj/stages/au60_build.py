"""旧阶段名 ``au60_build`` 的兼容 shim。

2026-11-13 起规范名是 :class:`smj.stages.electrode_build.ElectrodeBuildStage`
（``electrode_build``）。本模块不注册阶段，只保留旧的 import/类名，避免第三方
脚本在导入时断掉；CLI 的 ``smj stage au60_build ...`` 由
:func:`smj.state.normalize_stage` 映射到 ``electrode_build``。
"""

from __future__ import annotations

from .electrode_build import ElectrodeBuildStage

# 旧类名保留为别名。
Au60BuildStage = ElectrodeBuildStage

__all__ = ["ElectrodeBuildStage", "Au60BuildStage"]
