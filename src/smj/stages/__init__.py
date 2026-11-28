"""计算阶段的实现（按 :data:`smj.state.STAGE_ORDER` 顺序导入并注册）。

2026-11-13 起首段规范名为 ``molecule_opt``（旧名 ``gaussian_opt`` 保留为
兼容 shim），第三段规范名为 ``electrode_build``（旧名 ``au60_build`` 保留为
兼容 shim）。两个旧模块都不重复注册；导入本模块即完成注册，
:func:`smj.stages.base.get_stage` 之后可按名字取阶段实例。
"""

from __future__ import annotations

from .base import Collected, Prepared, Stage, all_stages, get_stage, register

# 具体阶段（每个模块用 @register 装饰器登记；导入顺序无依赖）
from . import molecule_opt  # noqa: F401,E402
from . import aims_tip_relax  # noqa: F401,E402
from . import electrode_build  # noqa: F401,E402
from . import au_relax  # noqa: F401,E402
from . import aims_scf  # noqa: F401,E402
from . import aitranss_te  # noqa: F401,E402
from . import junction_analyze  # noqa: F401,E402

# 热导链（D12：独立于六阶段电子链，只按名字显式调用）
from . import thermal_device  # noqa: F401,E402
from . import thermal_relax  # noqa: F401,E402
from . import thermal_fc  # noqa: F401,E402
from . import thermal_analyze  # noqa: F401,E402

__all__ = ["Stage", "Prepared", "Collected", "register", "get_stage", "all_stages"]
