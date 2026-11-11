"""六个计算阶段的实现（按 :data:`smj.state.STAGE_ORDER` 顺序导入并注册）。

导入本模块即完成注册；:func:`smj.stages.base.get_stage` 之后可按名字取阶段实例。
"""

from __future__ import annotations

from .base import Collected, Prepared, Stage, all_stages, get_stage, register

# 具体阶段（每个模块用 @register 装饰器登记；导入顺序无依赖）
from . import gaussian_opt  # noqa: F401,E402
from . import aims_tip_relax  # noqa: F401,E402
from . import au60_build  # noqa: F401,E402
from . import aims_scf  # noqa: F401,E402
from . import aitranss_te  # noqa: F401,E402
from . import junction_analyze  # noqa: F401,E402

__all__ = ["Stage", "Prepared", "Collected", "register", "get_stage", "all_stages"]
