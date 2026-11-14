"""smj-toolkit —— 单分子结 (single-molecule junction) 端到端计算工具包。

一条完整的计算链（每一步都可以单独运行）：

    molecule_opt     分子 XYZ -> FHI-aims PBE 优化（默认）或 Gaussian -> opt.xyz（阶段名 molecule_opt，旧名 gaussian_opt 兼容）
    aims_tip_relax   opt.xyz -> 加 Au 尖端 -> FHI-aims 弛豫 -> geometry.in.next_step
    electrode_build       geometry.in.next_step -> 拼接 Au(111) 电极 -> geometry.in（本地）
    aims_scf         geometry.in -> FHI-aims SCF -> omat.aims/mos.aims/aims.out
    aitranss_te      SCF 产物 -> Aitranss 透射谱 -> TE.dat
    junction_analyze TE.dat -> Landauer 电导/Seebeck/ZT -> results.json（本地）

设计原则：
- 每个阶段都是一个 :class:`smj.stages.base.Stage` 子类，``prepare()`` 只生成输入，
  ``collect()`` 只解析产物；阶段之间仅通过 ``<run>/<stage>/`` 目录里的文件耦合。
- 独立运行：``smj stage <name> --prepare/--submit/--wait/--collect/--run``，
  上游产物可以用 ``--xyz/--geometry/--scf-dir/--te`` 显式指定，无需重跑整条链。
- 真实计算在这台机器之外（Slurm 集群）；本地只做输入装配、提交、等待、解析与分析。
"""

__version__ = "0.2.0"

from .errors import (  # noqa: F401
    ConfigError,
    DependencyMissing,
    SmjError,
    StageBlocked,
    StageError,
    TransportError,
)

__all__ = [
    "__version__",
    "ConfigError",
    "DependencyMissing",
    "SmjError",
    "StageBlocked",
    "StageError",
    "TransportError",
]
