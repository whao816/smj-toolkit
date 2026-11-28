"""``smj.yaml`` 的加载、校验与路径解析。

一个配置文件描述**一条 run**（分子 + 电极 + 输运 + 集群）。相对路径一律相对
配置文件所在目录解析，所以 ``smj.yaml`` 可以和 ``molecule.xyz`` 放在同一个目录里。

与 AgenticDFT 的 ``config.yaml`` 不共享代码：独立包只保留计算所需的字段，字段语义
对齐 ``tools/junction`` 输入写出器的 ``env`` 约定（``g16root`` / ``aims_bin`` /
``aitranss_bin`` … 都可以在 ``cluster.env`` 里显式给出）。
"""

from __future__ import annotations

import math
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import yaml

from .errors import ConfigError
from .cluster_paths import validate_run_name
from .electrode_library import DEFAULT_CLUSTER
from .junction.electrode import FCC_LATTICE_A

DEFAULT_CONFIG_NAME = "smj.yaml"
MODES = ("dry-run", "local", "slurm-ssh")
_RUN_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_ELEMENT_RE = re.compile(r"^[A-Za-z]{1,3}$")


# --------------------------------------------------------------------------- #
# 基础工具
# --------------------------------------------------------------------------- #
def _expand(value: str) -> str:
    return os.path.expanduser(os.path.expandvars(str(value)))


def _fmt(where: str, err: Exception | str) -> str:
    return "%s: %s" % (where, err)


def _as_str(value: Any, where: str, *, allow_empty: bool = False) -> str:
    if value is None:
        return ""
    if not isinstance(value, (str, int, float)):
        raise ConfigError(_fmt(where, "应为字符串，实际是 %r" % (value,)))
    out = _expand(str(value)).strip()
    if not out and not allow_empty:
        raise ConfigError(_fmt(where, "不能为空"))
    return out


def _as_int(value: Any, where: str, *, minimum: int | None = None) -> int:
    if value is None:
        raise ConfigError(_fmt(where, "缺少整数值"))
    if isinstance(value, bool) or not isinstance(value, (int, str)):
        raise ConfigError(_fmt(where, "应为整数，实际是 %r" % (value,)))
    try:
        out = int(str(value).strip())
    except ValueError:
        raise ConfigError(_fmt(where, "应为整数，实际是 %r" % (value,))) from None
    if minimum is not None and out < minimum:
        raise ConfigError(_fmt(where, "应 >= %d，实际是 %d" % (minimum, out)))
    return out


def _as_float(value: Any, where: str) -> float:
    if value is None:
        raise ConfigError(_fmt(where, "缺少数值"))
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ConfigError(_fmt(where, "应为数值，实际是 %r" % (value,)))
    try:
        return float(str(value).strip())
    except ValueError:
        raise ConfigError(_fmt(where, "应为数值，实际是 %r" % (value,))) from None


def _as_bool(value: Any, where: str) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, int) and value in (0, 1):
        return bool(value)
    text = str(value).strip().lower()
    if text in ("1", "true", "yes", "on"):
        return True
    if text in ("0", "false", "no", "off", ""):
        return False
    raise ConfigError(_fmt(where, "应为布尔值 (true/false)，实际是 %r" % (value,)))


def _as_choice(value: Any, choices: tuple[str, ...], where: str) -> str:
    text = str(value).strip().lower()
    if text not in choices:
        raise ConfigError(_fmt(where, "应为 %s 之一，实际是 %r"
                              % ("/".join(choices), value)))
    return text


def _section(data: Mapping[str, Any], key: str, where: str, warnings: list[str]) -> Mapping[str, Any]:
    value = data.get(key)
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise ConfigError(_fmt("%s.%s" % (where, key), "应为映射 (yaml 缩进块)"))
    return value


def _warn_unknown(where: str, data: Mapping[str, Any], known: set[str], warnings: list[str]) -> None:
    unknown = sorted(set(data) - known - {"$schema"})
    if unknown:
        warnings.append("%s 含未知字段，已忽略: %s" % (where, ", ".join(unknown)))


# --------------------------------------------------------------------------- #
# 配置段
# --------------------------------------------------------------------------- #
@dataclass
class RunInfo:
    """运行标识与执行模式。"""

    name: str
    dir: str = "."                    # 相对 smj.yaml 的 run 目录
    mode: str = "slurm-ssh"           # dry-run | local | slurm-ssh
    clean: bool = False               # stage --prepare 时先清空该阶段目录
    keep_going: bool = False          # smj run：某阶段失败后是否继续（默认停）


@dataclass
class MoleculeConfig:
    """分子优化阶段参数。

    ``opt_engine`` 选择分子预优化引擎：``aims``（默认，与后续 tip/SCF 同为
    FHI-aims PBE/NAO）或 ``gaussian``（兼容旧 run）。Gaussian 分支继续用
    ``method``/``basis``；FHI-aims 分支用 ``opt_functional``/``opt_basis``。
    """

    xyz: str | None = None
    charge: int = 0
    multiplicity: int = 1
    #: 开壳层初始磁矩（**1-based 原子序号 -> μB**）。``spin collinear`` 下 FHI-aims
    #: 要求给出初始自旋密度，否则默认初猜可能收敛到错误自旋分布（甚至更慢）。
    #: 只影响收敛路径/落到哪个解，不改变 ``fixed_spin_moment``（=M−1）约束的总自旋；
    #: 各分量之和不必等于 M−1，但引擎会在 notes 里提示。闭壳层设了会被拒绝。
    initial_moments: dict[int, float] | None = None
    #: 是否写 ``fixed_spin_moment``（=M−1）把总自旋钉住。False = 解除约束的对照解：
    #: 只写 ``spin collinear`` + 初始磁矩，让 SCF 自己决定自旋态（与固定解比较）。
    fixed_spin_moment: bool = True
    #: 破缺对称（broken-symmetry, BS）开壳层单重态：M=1，但写 ``spin collinear`` +
    #: 两个自由基中心**反号**的 ``initial_moments``，且**不写** ``fixed_spin_moment``
    #: （2S=M−1=0 会把它直接压回闭壳层）。要求 ``multiplicity: 1``、初矩合计 0、
    #: 至少一正一负。FHI-aims 的单行列式只能给出 BS 近似解，且不打印 ⟨S²⟩：
    #: 是否真的落到开壳层单重态，必须由自旋密度/自然轨道占据数另行确认。
    broken_symmetry: bool = False
    method: str = "b3lyp"          # 仅 opt_engine=gaussian
    basis: str = "6-31g(d)"        # 仅 opt_engine=gaussian
    family: str | None = None  # 判据标定族（如 AZ/BDT），见 criteria.SYSTEM_FAMILIES
    opt_engine: str = "aims"       # aims | gaussian
    opt_functional: str = "pbe"    # 仅 aims；与 aims_tip_relax/aims_scf 一致
    opt_basis: str = "light"       # 仅 aims，记录用；实际 species 由 species_dir 决定
    # 仅 aims：写进 `relax_geometry bfgs <值>`；FHI-aims 口径是**最大力分量**，
    # 单位 eV/Å（不是 Hartree/Bohr）。默认 1.e-2 eV/Å ≈ 1.9e-4 Ha/Bohr，已经偏紧。
    opt_force_threshold: str = "1.e-2"
    #: 仅 Gaussian 分支：开壳层优化加 ``Stable=Opt``（先做波函数稳定性检查再优化），
    #: 用于怀疑当前解不是最低能量自旋态时。默认关闭——它让每次优化多一遍稳定性分析。
    stable_check: bool = False


@dataclass
class ElectrodeConfig:
    """电极 / 结几何参数（electrode_build 与 aims_tip_relax、aitranss_te 共用）。"""

    au_per_electrode: int | None = 30     # 仅 builtin 生效；None = 完整 60/侧；library 忽略
    left_anchor: str = "S"
    right_anchor: str = "S"
    au_tip_distance: float = 2.35
    #: tip 放置目标（Å）。None = 只用 au_tip_distance 放初始几何，不做任何锚–金属
    #: 距离约束；最终键长/角度完全由 aims_tip_relax 的全原子弛豫决定。
    #: 2.40 Å 只对**去质子硫醇盐**有参考意义；-SCH3 硫醚或 N/Se 锚不应把它当目标键长。
    anchor_au_s_target: float | None = None
    #: True = 把两个锚定原子和两个 tip 一起冻结（旧 W0-1 距离扫描口径；
    #: pinned 对照/敏感性专用，只有用户明确要求时才用）。默认 False = 全原子自由弛豫。
    pin_contact: bool = False
    #: tip 放置方向：``axis``（默认）= 沿两锚原子连线外推（旧行为逐字节不变）；
    #: ``lone_pair`` = 沿锚原子的孤对方向外推（S 为例：C–S–C 角平分线的反向）。
    #: 间位取代的 SMe 分子里 S···S 连线可能与 S–C(Me) 键几乎同向，沿 axis 放 tip
    #: 会把金属压进甲基（Au–C(Me) ~1.2 Å），此时用 lone_pair 才是物理接触方向。
    tip_direction: str = "axis"
    #: 去质子口径：None = 仅当锚定原子是 S 时自动移除 S–H（N/Se 等其它锚默认保留 H）；
    #: True = 强制移除锚定原子上的 H（每个锚原子恰好 1 个）；False = 完全保留。
    strip_thiol_h: bool | None = None
    vdw: bool = True
    nlayers: int | None = None            # Aitranss 中心区层数；None = 按库文件头部建议
    source: str = "library"               # library = FHI-aims 预定义团簇；builtin = 原三角岛
    library_dir: str | None = None        # electrodes.library 目录；None = 常见路径/远端探测
    library_cluster: str = DEFAULT_CLUSTER
    #: 电极金属元素（fcc(111) 单晶表面）。tip 也用该元素；builtin 直接按它生成三角岛，
    #: library 模式下必须与库团簇元素一致（否则 electrode_build 拒绝装配）。
    element: str = "Au"
    #: fcc 晶格常数（Å）。None = 用 FCC_LATTICE_A 的登记值（Au 4.07823、Ag 4.08530 …）；
    #: 未登记元素必须显式给出。
    a_lattice: float | None = None
    #: 旧口径兼容：仅为显式重跑 legacy ``au_relax`` 阶段保留。主链路在
    #: ``electrode_build`` 后即冻结几何，不再消费此字段。
    relax: bool = False
    relax_force_threshold: str = "1.e-2"  # BFGS 力收敛阈值（eV/Å，最大力分量）
    #: 旧口径兼容：``au_relax`` 冻结谁（backplane / electrodes）；主链路不使用
    relax_freeze: str = "backplane"


#: 1 Ha = 27.211386245988 eV（CODATA 2018）。
HARTREE_EV = 27.211386245988

#: ``transport.energy_window_unit`` 的合法取值。
WINDOW_UNITS = ("eV", "Hartree")


#: ``transport.extra_outputs`` 的白名单（键名 -> 由 ``junction.tcontrol`` 翻成
#: ``control.in`` 里的 ``output`` 行）。``spin_density`` 只对共线自旋（开壳层）有意义，
#: 由 ``aims_scf`` 阶段再校验。两份清单必须一致，测试 ``test_extra_outputs_*`` 守着。
EXTRA_OUTPUTS = ("hirshfeld", "spin_density")


@dataclass
class TransportConfig:
    """透射谱 / Landauer 分析参数。

    ``energy_window`` 默认是 **相对费米能级的 eV**（``E − E_F``，默认 ±2 eV）；
    旧配置里的 Ha 绝对能量三元组仍可读（``energy_window_unit: Hartree`` 或按量级
    自动判定），两种口径最终都在 ``aitranss_te`` 里化成 aitranss 需要的 Ha 绝对值。
    """

    backend: str = field(default="aitranss", kw_only=True)
    temperature: float = 300.0
    #: [起点, 步长, 终点]，单位见 ``energy_window_unit``。默认 eV 且以 E_F 为参考：
    #: E − E_F ∈ [−2, +2] eV、步长 5 meV（801 点）。换金属后 E_F 会移动，但相对
    #: 窗口不跟着漂；只算电导可收窄，想看峰再放宽。
    energy_window: tuple[float, float, float] = (-2.0, 0.005, 2.0)
    #: "eV"（默认）= 相对 E_F 的 eV；"Hartree" = 旧口径，Ha 绝对能量。
    energy_window_unit: str = "eV"
    #: 显式覆盖：分量按 **Ha 绝对能量** 解释（优先级最高，可单独给一个）。
    ener: float | None = None
    estep: float | None = None
    eend: float | None = None
    #: collect 时是否要求 aitranss 主输出 ``te.out`` 存在（结果可复核的一部分）。
    #: 仅外部/历史 run 才应显式设为 false 豁免；不改变物理结果，故不进 protocol。
    require_te_out: bool = True
    #: 开壳层 TE.dat 只有 1 个自旋通道时是否放行。默认阻断：tcontrol 若仍用
    #: ``$scfmo`` 或 beta 波函数缺失，aitranss 会把 α/β 合并成均值，自旋极化不可
    #: 恢复；只有确知该谱只能单列（外部/历史数据）才显式打开。只影响 collect 闸门，
    #: 不改变求解与积分口径，故与 ``require_te_out`` 一样**不进** protocol。
    allow_single_channel: bool = False
    #: 额外请求的 FHI-aims 输出（静态 SCF，``aims_scf``）：``hirshfeld``、``spin_density``。
    #: 开壳层判据要"自旋是否留在分子上"的第一手证据时才打开；默认关闭，避免所有 run
    #: 产物变重。只影响产物、不影响求解参数，故与 ``require_te_out`` 一样**不进** protocol。
    extra_outputs: tuple[str, ...] = ()

    # -- 静态 SCF 的收敛参数（aims_scf）--------------------------------------- #
    #: ``occupation_type gaussian`` 的展宽宽度（eV）。默认 0.01 是分子口径；有限 Au
    #: 团簇的近 E_F 态很密，0.01 eV 下占据数逐迭代跳变，配合 fixed_spin_moment 常
    #: 落进 limit cycle（实测 152 原子 IFS 结 40 次迭代 Δρ 卡在 0.02-0.08）。金属
    #: 电极结给 0.05（ttf-fixed 先例：88 原子结 85 次迭代收敛）。
    occupation_width: float = 0.01
    #: ``charge_mix_param``：Pulay 混合步长。金属 + 固定自旋时 0.2 会震荡，0.05 稳定。
    charge_mix_param: float = 0.2
    #: ``n_max_pulay``：Pulay 历史长度。加到 20-30 通常能救回缓慢收敛的金属结。
    n_max_pulay: int = 10

    def relative_window(self) -> tuple[float, float, float]:
        """相对 E_F 的 eV 三元组（仅 ``energy_window_unit: eV`` 口径）。"""
        if self.energy_window_unit != "eV":
            raise ValueError(
                "energy_window_unit=%r 是 Ha 绝对口径，没有相对 E_F 的三元组"
                % (self.energy_window_unit,))
        return (float(self.energy_window[0]), float(self.energy_window[1]),
                float(self.energy_window[2]))

    def resolved_window(self, efermi_hartree: float | None = None
                        ) -> tuple[float, float, float]:
        """实际写进 tcontrol 的 Ha 绝对窗口 ``(ener, estep, eend)``。

        eV 口径（默认）需要 SCF 的 ``efermi`` 才能把相对能量平移到 Ha 绝对坐标；
        缺 ``efermi`` 时**不猜**，直接报错。``ener/estep/eend`` 显式覆盖恒按 Ha
        绝对能量解释，逐分量替换。
        """
        ener, estep, eend = (float(v) for v in self.energy_window)
        if self.energy_window_unit == "eV":
            estep = estep / HARTREE_EV
            if self.ener is None or self.eend is None:
                if efermi_hartree is None:
                    raise ValueError(
                        "energy_window 为 eV 相对口径（默认 E_F ±2 eV），需要 SCF 的 "
                        "efermi(Ha) 才能换成 aitranss 的 Ha 绝对窗口；"
                        "请先跑 aims_scf，或显式给 transport.ener/eend（Ha 绝对）")
                if self.ener is None:
                    ener = float(efermi_hartree) + ener / HARTREE_EV
                if self.eend is None:
                    eend = float(efermi_hartree) + eend / HARTREE_EV
        if self.ener is not None:
            ener = float(self.ener)
        if self.estep is not None:
            estep = float(self.estep)
        if self.eend is not None:
            eend = float(self.eend)
        return ener, estep, eend


@dataclass
class ThermalElectrodeConfig:
    """声子热导的电极几何。

    **刻意不使用电子链的金字塔团簇**（``electrode`` 段继续管电子链，两者互不影响）：
    声子 NEGF 需要半无限 lead 与横向元胞，只有二维周期板才给得出来（方案 §4.1）。
    """

    #: Au fcc 晶格常数(Å)。最近邻 = a/√2 = 2.8845 Å。
    a_lattice: float = 4.078230
    #: 面内超胞 (nx, ny)：3×3 给 9 原子/层，总原子数压到 ~90。
    in_plane: tuple[int, int] = (3, 3)
    #: 每侧 lead 层数。
    layers: int = 4
    #: 每侧参与位移的层数（靠近分子的那些层）；其余层用冻结 lead 力常数。
    device_layers: int = 2
    #: 沿 c 轴的真空厚度(Å)。二维板沿 c 是真空，所以 kz 恒为 1。
    vacuum: float = 20.0
    #: 表面 Au 平面到 S 的距离(Å)。
    au_s_distance: float = 2.35
    #: 两个 S 之间的距离(Å)。苯/萘二硫醇 8.0–8.7 Å。
    gap: float = 8.4
    element: str = "Au"
    anchor_element: str = "S"
    #: 建结前去 S–H。**默认必须为真**：molecule_opt 产出的是中性二硫醇，
    #: 直接贴 Au 就是报告 §6 陷阱⑥。关掉只在明确想复现"带 H 的错模型"时用。
    strip_thiol_h: bool = True
    #: 弛豫时冻结的最外层数（每侧）。
    freeze_outer_layers: int = 1


@dataclass
class ThermalFCConfig:
    """有限位移力常数扫描的数值参数。"""

    #: 中心差分位移幅度(Å)。0.01 Å 是 Sadeghi–Lambert 组的标准做法；
    #: 本仓更小会让力噪声(1e-4 eV/Å)主导，更大则非谐截断误差上升。
    displacement: float = 0.01
    #: 力常数截断半径(Å)。<7 Å 直接拒绝：4 层 lead 沿 z 的跨度就是 7.06 Å。
    cutoff: float = 10.0
    #: Slurm 数组并发数。gpu1 是 64 核、每作业 8 核 → 上限 8（不是 16）。
    array_parallel: int = 8
    array_time_limit: str = "06:00:00"
    #: ``--mem``。实测 92 原子峰值 871.9 MB/task；12G×8=96G < 124.7G 物理内存。
    memory: str = "12G"
    #: ``sc_accuracy_forces``(eV/Å)：力的收敛阈值，必须与位移幅度配套。
    sc_accuracy_forces: str = "5E-3"
    #: 弛豫的 ``relax_geometry bfgs`` 阈值(eV/Å)。
    relax_force_threshold: str = "1.e-2"
    #: 面内 k 网格。``kz`` 必须是 1（沿 c 是真空）。
    k_grid: tuple[int, int, int] = (3, 3, 1)


@dataclass
class ThermalNEGFConfig:
    """声子 NEGF 的谱参数。

    这两个数是**冻结的生产协议**（方案 D15），不是可随手调的收敛参数：

    * η 是数值展宽，不是物理线宽。同一份力常数上 η 从 0.5 降到 0.05 meV，
      ``G_ph(300 K)`` 从 10.9 单调升到 13.72 pW/K（η→0 外推 ≈13.8）；
      0.3 meV 给 12.37。**跨分子比较必须用同一个 η**，否则 25% 的差别
      全部来自这个参数。
    * 步长再加密也压不下 1–2% 的梯形求积底噪（η=0.05/0.1/0.2/0.3 下分别
      1.99%/1.58%/2.46%/2.21%），所以判据门限重标到 5%（D14），步长取 0.05
      只是为了把谱画得准、不是为了收敛。
    """

    #: 虚部展宽 η(meV)。冻结协议值：0.05 meV。
    eta_meV: float = 0.05
    #: 能量步长(meV)。冻结协议值：0.05 meV（0–19.5 meV 共 391 点）。
    energy_step_meV: float = 0.05
    #: Au 声子带顶(meV)。实测 fcc-Au 最高支在 L 点 4.71 THz = 19.5 meV
    #: (Lynn, Smith & Nicklow, PRB 8, 3493 (1973))。超过它的弹性透射必须为零。
    band_top_meV: float = 19.5


@dataclass
class ThermalConfig:
    """声子热导链的配置（``thermal`` 段）。

    整段默认关闭：不写 ``thermal`` 段时 ``enabled=False``，电子链行为零变化。
    """

    enabled: bool = False
    #: 计算温度(K)。``None`` = 继承 ``transport.temperature``（默认 300 K）。
    temperature_K: float | None = None
    #: 电极力常数（lead FC）的位置。``"auto"`` = 用随包发布的 Au(111) 资产
    #: （``smj.phonon.assets``），写路径（目录或 JSON）则覆盖它。
    #:
    #: 为什么默认给资产而不是要求先算一次：电极是金、与分子无关，那份力常数
    #: 只需要算一次（3×3×6 层、Γ、54 个 SCF）。要重算或换成自己的版本，
    #: 用 ``smj lead-fc --export DIR`` 导出后按本字段指过去。
    lead_fc: str | None = "auto"
    electrode: ThermalElectrodeConfig = field(default_factory=ThermalElectrodeConfig)
    fc: ThermalFCConfig = field(default_factory=ThermalFCConfig)
    negf: ThermalNEGFConfig = field(default_factory=ThermalNEGFConfig)

    def resolved_temperature(self, fallback: float = 300.0) -> float:
        """``temperature_K`` 为 None 时回落到 ``transport.temperature``。"""
        if self.temperature_K is not None:
            return float(self.temperature_K)
        return float(fallback)


@dataclass
class ClusterConfig:
    """Slurm / SSH 集群与外部程序部署路径。"""

    host: str | None = None               # ~/.ssh/config 别名或 hostname
    user: str | None = None
    port: int = 22
    key: str | None = None                # 私钥路径（None = ssh 默认）
    remote_root: str | None = None        # 默认远端登录账户的 ~/smj-runs
    partition: str | None = None
    ntasks: int | None = None             # 覆盖 Gaussian/aims 的并行度（None = 自动）
    gres: str | None = None               # 如 gpu:1
    account: str | None = None
    qos: str | None = None
    species_dir: str | None = None        # FHI-aims species_defaults 目录
    time_limit: str | None = None         # 覆盖 --time（默认按资源估算）
    mem: str | None = None                # 覆盖 --mem（默认按资源估算）
    poll_seconds: float = 60.0
    timeout_hours: float = 6.0
    env: dict[str, str] = field(default_factory=dict)  # 原始 smj env 键，如 g16root/aims_bin

    @property
    def ssh_alias(self) -> str:
        return self.host or ""

    @property
    def ssh_host(self) -> str:
        return self.host or ""

    @property
    def ssh_target(self) -> str:
        """``user@host``（user 为空时就是 host）。"""
        host = self.host or ""
        if self.user and "@" not in host:
            return "%s@%s" % (self.user, host)
        return host


@dataclass
class PathsConfig:
    """可选路径覆盖；默认使用包内 vendored 模板。"""

    templates: str | None = None
    species_defaults: str | None = None


@dataclass
class ReviewConfig:
    """结构复核 / 报告（附带功能）的开关。

    2026-11-15 更名：原「A2 人工确认」→「表面核对」→「连接核对」（connection check）
    ——核对每侧三个表面原子与 tcontrol/输运轴是否一致，并把几何指纹绑进记录。
    默认**自动登记**（``auto_connection_check``）：electrode_build 采集后自动审计
    并写记录，不再需要手工 ``smj review --approve --who``。

    旧键 ``surface_check`` / ``surface_check_deviation`` / ``auto_surface_check`` 与更早的
    ``require_a2`` / ``a2_deviation`` 仍接受（等价新键）；旧 run 目录里的
    ``surface_check.json`` / ``a2_confirmation.json`` 也仍被读取（见 :mod:`smj.review`）。
    """

    connection_check: bool = True            # 连接核对是否作为复核闸门（默认开）
    connection_check_deviation: float = 0.15 # 表面索引偏差（1-交集/较大集）的门槛
    auto_connection_check: bool = True       # 默认自动登记连接核对（无需人工签署）
    accept_deviation: float | None = None # 显式人工豁免的放宽门槛

    # -- 旧名兼容（先叫 A2，后叫「表面核对」） ------------------------------ #
    @property
    def surface_check(self) -> bool:
        return self.connection_check

    @surface_check.setter
    def surface_check(self, value: bool) -> None:
        self.connection_check = bool(value)

    @property
    def surface_check_deviation(self) -> float:
        return self.connection_check_deviation

    @surface_check_deviation.setter
    def surface_check_deviation(self, value: float) -> None:
        self.connection_check_deviation = float(value)

    @property
    def auto_surface_check(self) -> bool:
        return self.auto_connection_check

    @auto_surface_check.setter
    def auto_surface_check(self, value: bool) -> None:
        self.auto_connection_check = bool(value)

    # -- 旧名兼容（更名前叫 A2） ------------------------------------------- #
    @property
    def require_a2(self) -> bool:
        return self.connection_check

    @require_a2.setter
    def require_a2(self, value: bool) -> None:
        self.connection_check = bool(value)

    @property
    def a2_deviation(self) -> float:
        return self.connection_check_deviation

    @a2_deviation.setter
    def a2_deviation(self, value: float) -> None:
        self.connection_check_deviation = float(value)


@dataclass
class RunConfig:
    """一条 smj run 的完整配置。"""

    run: RunInfo
    molecule: MoleculeConfig = field(default_factory=MoleculeConfig)
    electrode: ElectrodeConfig = field(default_factory=ElectrodeConfig)
    transport: TransportConfig = field(default_factory=TransportConfig)
    thermal: ThermalConfig = field(default_factory=ThermalConfig)
    siesta: dict[str, Any] = field(default_factory=dict, kw_only=True)
    cluster: ClusterConfig = field(default_factory=ClusterConfig)
    paths: PathsConfig = field(default_factory=PathsConfig)
    review: ReviewConfig = field(default_factory=ReviewConfig)
    source: Path | None = None
    warnings: list[str] = field(default_factory=list)

    # -- 基本属性 ---------------------------------------------------------- #
    @property
    def name(self) -> str:
        return self.run.name

    @property
    def mode(self) -> str:
        return self.run.mode

    @property
    def base_dir(self) -> Path:
        return self.source.parent if self.source else Path.cwd()

    @property
    def run_dir(self) -> Path:
        """``<run.dir>/<run.name>``：yaml 里的 dir 是 run 根，不是某个 run 本身。"""
        base = self.resolve(self.run.dir) or self.base_dir
        return base / self.run.name

    def resolve(self, value: str | Path | None) -> Path | None:
        """相对路径按配置文件目录解析（``~`` 与 ``${VAR}`` 先展开）。"""
        if value is None:
            return None
        text = _expand(str(value))
        if not text:
            return None
        path = Path(text)
        if path.is_absolute():
            return path
        return (self.base_dir / path).resolve()

    # -- 便捷方法 ---------------------------------------------------------- #
    def molecule_xyz(self, override: str | None = None) -> Path | None:
        raw = override or self.molecule.xyz
        return self.resolve(raw)


# --------------------------------------------------------------------------- #
# 解析
# --------------------------------------------------------------------------- #
_KNOWN_SECTIONS = {"siesta", "run", "molecule", "electrode", "transport", "thermal",
                   "cluster", "paths", "review"}

_RUN_KEYS = {"name", "dir", "mode", "clean", "keep_going"}
_MOLECULE_KEYS = {
    "xyz", "charge", "multiplicity", "initial_moments", "fixed_spin_moment",
    "broken_symmetry",
    "method", "basis", "family",
    "opt_engine", "opt_functional", "opt_basis", "opt_force_threshold",
    "stable_check",
}
_ELECTRODE_KEYS = {
    "au_per_electrode", "left_anchor", "right_anchor", "au_tip_distance",
    "anchor_au_s_target", "pin_contact", "tip_direction", "strip_thiol_h", "vdw", "nlayers",
    "source", "library_dir", "library_cluster", "element", "a_lattice",
    "relax", "relax_force_threshold", "relax_freeze",
}
_TRANSPORT_KEYS = {"backend", "temperature", "energy_window", "energy_window_unit",
                   "ener", "estep", "eend", "require_te_out", "allow_single_channel",
                   "extra_outputs", "occupation_width", "charge_mix_param",
                   "n_max_pulay"}
_THERMAL_KEYS = {"enabled", "temperature", "lead_fc", "electrode", "fc", "negf"}
_THERMAL_ELECTRODE_KEYS = {
    "a_lattice", "in_plane", "layers", "device_layers", "vacuum",
    "au_s_distance", "gap", "element", "anchor_element", "strip_thiol_h",
    "freeze_outer_layers",
}
_THERMAL_FC_KEYS = {
    "displacement", "cutoff", "array_parallel", "array_time_limit", "memory",
    "sc_accuracy_forces", "relax_force_threshold", "k_grid",
}
_THERMAL_NEGF_KEYS = {"eta_mev", "energy_step_mev", "band_top_mev"}
_CLUSTER_KEYS = {
    "host", "user", "port", "key", "remote_root", "partition", "ntasks", "gres",
    "account", "qos", "species_dir", "poll_seconds", "timeout_hours", "env",
    "time_limit", "mem",
}
_PATHS_KEYS = {"templates", "species_defaults", "species_dir"}
_REVIEW_KEYS = {"connection_check", "connection_check_deviation", "auto_connection_check",
                "accept_deviation",
                # 2026-11-15 更名前的旧键（仍接受，等价于 connection_check*）
                "surface_check", "surface_check_deviation", "auto_surface_check",
                # 更早的 A2 时代旧键
                "require_a2", "a2_deviation"}


def known_keys() -> dict[str, set[str]]:
    """各 section 允许的字段名（供 `scan_convergence` 之类工具校验轴键）。"""
    return {
        "run": set(_RUN_KEYS), "molecule": set(_MOLECULE_KEYS),
        "electrode": set(_ELECTRODE_KEYS), "transport": set(_TRANSPORT_KEYS),
        "thermal": set(_THERMAL_KEYS),
        "cluster": set(_CLUSTER_KEYS), "paths": set(_PATHS_KEYS),
        "review": set(_REVIEW_KEYS),
    }


def _parse_run(data: Mapping[str, Any], warnings: list[str]) -> RunInfo:
    _warn_unknown("run", data, _RUN_KEYS, warnings)
    name = data.get("name")
    if name is None:
        raise ConfigError("run.name 不能为空（每条 run 需要唯一名字，用作远端/本地目录名）")
    name = _as_str(name, "run.name")
    if not _RUN_NAME_RE.match(name):
        raise ConfigError("run.name 只允许字母/数字/._-，且不能以符号开头: %r" % name)
    mode = _as_str(data.get("mode") or "slurm-ssh", "run.mode")
    if mode not in MODES:
        raise ConfigError("run.mode 应为 %s 之一，实际是 %r" % ("/".join(MODES), mode))
    run_dir = _as_str(data.get("dir") or ".", "run.dir", allow_empty=True) or "."
    return RunInfo(
        name=name,
        dir=run_dir,
        mode=mode,
        clean=_as_bool(data.get("clean", False), "run.clean"),
        keep_going=_as_bool(data.get("keep_going", False), "run.keep_going"),
    )


def _parse_initial_moments(value: Any) -> dict[int, float] | None:
    """``molecule.initial_moments``：{1-based 原子序号: μB}。

    YAML 的键可能是整数也可能是字符串，这里统一成 ``int``；空映射等价于未设置。
    只做结构校验（序号 >=1、数值有限、|μB| 合理），原子序号上界要等读到结构
    才知道，留给引擎在 prepare 时 fail-closed。
    """
    if value is None:
        return None
    if not isinstance(value, Mapping):
        raise ConfigError(
            "molecule.initial_moments 应为映射（1-based 原子序号 -> μB），实际是 %r"
            % (value,))
    parsed: dict[int, float] = {}
    for key, raw in value.items():
        try:
            index = int(str(key).strip())
        except (TypeError, ValueError):
            raise ConfigError(
                "molecule.initial_moments 的键应为 1-based 原子序号（整数），实际是 %r"
                % (key,)) from None
        if index < 1:
            raise ConfigError(
                "molecule.initial_moments 的原子序号应 >= 1（1-based），实际是 %d" % index)
        moment = _as_float(raw, "molecule.initial_moments[%d]" % index)
        if not math.isfinite(moment) or abs(moment) > 20.0:
            raise ConfigError(
                "molecule.initial_moments[%d] 应为有限且 |μB| <= 20 的数值，实际是 %r"
                % (index, raw))
        parsed[index] = moment
    return dict(sorted(parsed.items())) or None


def _parse_molecule(data: Mapping[str, Any], warnings: list[str]) -> MoleculeConfig:
    _warn_unknown("molecule", data, _MOLECULE_KEYS, warnings)
    xyz = data.get("xyz")
    cfg = MoleculeConfig(
        xyz=_as_str(xyz, "molecule.xyz", allow_empty=True) or None,
        charge=_as_int(data.get("charge", 0), "molecule.charge"),
        multiplicity=_as_int(data.get("multiplicity", 1), "molecule.multiplicity", minimum=1),
        initial_moments=_parse_initial_moments(data.get("initial_moments")),
        fixed_spin_moment=_as_bool(
            data.get("fixed_spin_moment", True), "molecule.fixed_spin_moment"),
        broken_symmetry=_as_bool(
            data.get("broken_symmetry", False), "molecule.broken_symmetry"),
        method=_as_str(data.get("method") or "b3lyp", "molecule.method"),
        basis=_as_str(data.get("basis") or "6-31g(d)", "molecule.basis"),
        family=_as_str(data.get("family"), "molecule.family", allow_empty=True) or None,
        opt_engine=_as_str(data.get("opt_engine") or "aims", "molecule.opt_engine").lower(),
        opt_functional=_as_str(
            data.get("opt_functional") or "pbe", "molecule.opt_functional"),
        opt_basis=_as_str(data.get("opt_basis") or "light", "molecule.opt_basis"),
        opt_force_threshold=_as_str(
            data.get("opt_force_threshold") or "1.e-2", "molecule.opt_force_threshold"),
        stable_check=_as_bool(data.get("stable_check", False), "molecule.stable_check"),
    )
    from .engines import engine_names  # lazy: 避免配置模块顶层拉进产物解析依赖

    valid = engine_names()
    if cfg.opt_engine not in valid:
        raise ConfigError(
            "molecule.opt_engine 应为 %s，实际是 %r"
            % (" 或 ".join(valid), cfg.opt_engine))
    try:
        threshold = float(cfg.opt_force_threshold)
    except ValueError:
        raise ConfigError(
            "molecule.opt_force_threshold 应为数字（单位 eV/Å），实际是 %r"
            % cfg.opt_force_threshold) from None
    # 上限 0.1 eV/Å：再松就不是"优化过的分子"了；下限只是防写错（0 或负数会让
    # relax_geometry 永远不收敛）。
    if not 0.0 < threshold <= 0.1:
        raise ConfigError(
            "molecule.opt_force_threshold 应在 (0, 0.1] eV/Å（FHI-aims 的 relax_geometry "
            "判据是最大力分量，不是 Hartree/Bohr），实际是 %r" % cfg.opt_force_threshold)
    if cfg.broken_symmetry:
        # BS 单重态的四条硬前提：M=1、有初矩、合计 0（N↑=N↓）、正负都要有
        # （否则 SCF 从处处同号/为零的初猜出发只会落回闭壳层解）。
        if cfg.multiplicity != 1:
            raise ConfigError(
                "molecule.broken_symmetry=true 只适用于 multiplicity=1 的开壳层单重态，"
                "当前 multiplicity=%d；三重态（M>1）用 multiplicity 本身即可"
                % cfg.multiplicity)
        moments = cfg.initial_moments or {}
        if not moments:
            raise ConfigError(
                "molecule.broken_symmetry=true 需要 molecule.initial_moments："
                "给两个自由基中心各一个反号初猜（如 {7: 1.0, 18: -1.0}）；"
                "否则 SCF 没有破缺对称的起点")
        total = sum(moments.values())
        if abs(total) > 1e-6:
            raise ConfigError(
                "molecule.broken_symmetry=true 时 molecule.initial_moments 的合计应为 0"
                "（BS 单重态 N↑=N↓），当前合计 %g μB" % total)
        if not (any(v > 0 for v in moments.values())
                and any(v < 0 for v in moments.values())):
            raise ConfigError(
                "molecule.broken_symmetry=true 时 molecule.initial_moments 至少要有一个"
                "正、一个负分量；全同号的初猜收敛到的是固定自旋的另一种解")
    return cfg


def _parse_electrode(data: Mapping[str, Any], warnings: list[str]) -> ElectrodeConfig:
    _warn_unknown("electrode", data, _ELECTRODE_KEYS, warnings)
    source = _as_str(data.get("source") or "library", "electrode.source")
    if source not in ("library", "builtin"):
        raise ConfigError("electrode.source 应为 library/builtin 之一，实际是 %r" % source)
    au = data.get("au_per_electrode", 30)
    au_value = None if au is None else _as_int(au, "electrode.au_per_electrode", minimum=1)
    if au_value is not None and au_value > 30:
        raise ConfigError(
            "electrode.au_per_electrode 最多 30（builtin 按完整原子层截断），实际是 %d" % au_value)
    if source == "library":
        # 0.2.0rc1：library 模式完全忽略该键 —— 团簇大小由 electrodes.library 文件决定。
        if "au_per_electrode" in data:
            warnings.append(
                "electrode.au_per_electrode 仅对 source: builtin 生效；"
                "当前 source: library 按库团簇整簇贴装，该键已忽略")
        au_value = None
    target_raw = data.get("anchor_au_s_target", None)
    target = None if target_raw is None else _as_float(
        target_raw, "electrode.anchor_au_s_target")
    if target is not None and not (1.5 <= target <= 4.0):
        raise ConfigError("electrode.anchor_au_s_target 应在 1.5–4.0 Å，实际是 %s" % target)
    nlayers_raw = data.get("nlayers", None)
    nlayers = None if nlayers_raw is None else _as_int(nlayers_raw, "electrode.nlayers", minimum=1)
    if nlayers is not None and not (1 <= nlayers <= 8):
        raise ConfigError("electrode.nlayers 应在 1–8，实际是 %d（null = 按库文件建议）" % nlayers)
    library_dir = data.get("library_dir")
    library_dir = _as_str(library_dir, "electrode.library_dir") if library_dir else None
    cluster_rel = _as_str(data.get("library_cluster") or DEFAULT_CLUSTER, "electrode.library_cluster")
    if not cluster_rel.strip():
        raise ConfigError("electrode.library_cluster 不能为空（应为 library_dir 下的相对路径）")
    element = _as_str(data.get("element") or "Au", "electrode.element")
    if not _ELEMENT_RE.match(element):
        raise ConfigError(
            "electrode.element 应为元素符号（如 Au/Ag/Cu/Pt），实际是 %r" % element)
    element = element[0].upper() + element[1:].lower()
    a_lattice_raw = data.get("a_lattice", None)
    a_lattice = None if a_lattice_raw is None else _as_float(
        a_lattice_raw, "electrode.a_lattice")
    if a_lattice is not None and not (2.0 <= a_lattice <= 6.0):
        raise ConfigError("electrode.a_lattice 应在 2.0–6.0 Å，实际是 %s" % a_lattice)
    if a_lattice is None and element not in FCC_LATTICE_A:
        raise ConfigError(
            "electrode.element %r 没有内置 fcc 晶格常数：请显式给出 electrode.a_lattice"
            "（已登记: %s）" % (element, "/".join(sorted(FCC_LATTICE_A))))
    threshold = _as_str(data.get("relax_force_threshold") or "1.e-2", "electrode.relax_force_threshold")
    freeze = _as_str(data.get("relax_freeze") or "backplane", "electrode.relax_freeze")
    if freeze not in ("backplane", "electrodes"):
        raise ConfigError(
            "electrode.relax_freeze 应为 backplane 或 electrodes，实际是 %r" % freeze)
    cfg = ElectrodeConfig(
        au_per_electrode=au_value,
        left_anchor=_as_str(data.get("left_anchor") or "S", "electrode.left_anchor"),
        right_anchor=_as_str(data.get("right_anchor") or "S", "electrode.right_anchor"),
        au_tip_distance=_as_float(data.get("au_tip_distance", 2.35), "electrode.au_tip_distance"),
        anchor_au_s_target=target,
        pin_contact=_as_bool(data.get("pin_contact", False), "electrode.pin_contact"),
        tip_direction=_as_choice(data.get("tip_direction", "axis"),
                                 ("axis", "lone_pair"), "electrode.tip_direction"),
        strip_thiol_h=(None if data.get("strip_thiol_h") is None
                       else _as_bool(data.get("strip_thiol_h"), "electrode.strip_thiol_h")),
        vdw=_as_bool(data.get("vdw", True), "electrode.vdw"),
        nlayers=nlayers,
        source=source,
        library_dir=library_dir,
        library_cluster=cluster_rel,
        element=element,
        a_lattice=a_lattice,
        relax=_as_bool(data.get("relax", False), "electrode.relax"),
        relax_force_threshold=threshold,
        relax_freeze=freeze,
    )
    if cfg.pin_contact:
        warnings.append(
            "electrode.pin_contact=true 是 pinned 对照/距离扫描口径，不是生产默认；"
            "只有用户明确要求做 pinned 对照时才应启用（生产默认 false）")
    if not (1.0 <= cfg.au_tip_distance <= 5.0):
        raise ConfigError("electrode.au_tip_distance 应在 1.0–5.0 Å，实际是 %s" % cfg.au_tip_distance)
    return cfg


def _parse_transport(data: Mapping[str, Any], warnings: list[str]) -> TransportConfig:
    _warn_unknown("transport", data, _TRANSPORT_KEYS, warnings)
    window = data.get("energy_window", [-2.0, 0.005, 2.0])
    if not isinstance(window, (list, tuple)) or len(window) != 3:
        raise ConfigError("transport.energy_window 应为 [起点, 步长, 终点] 三个数，实际是 %r" % (window,))
    ener = _as_float(window[0], "transport.energy_window[0]")
    estep = _as_float(window[1], "transport.energy_window[1]")
    eend = _as_float(window[2], "transport.energy_window[2]")
    if estep <= 0 or eend <= ener:
        raise ConfigError("transport.energy_window 要求 步长 > 0 且 终点 > 起点（当前 %r）" % (window,))
    unit = _parse_window_unit(data.get("energy_window_unit"), ener, eend, warnings)
    if unit == "eV" and (eend - ener) > 20.0:
        warnings.append(
            "transport.energy_window 跨度 %.1f eV 偏大（默认 ±2 eV）；太宽会让 aitranss "
            "点数暴涨且远离 E_F 的透射没有意义" % (eend - ener))
    cfg = TransportConfig(
        backend=_as_choice(data.get("backend", "aitranss"), ("aitranss", "transiesta"), "transport.backend"),
        temperature=_as_float(data.get("temperature", 300.0), "transport.temperature"),
        energy_window=(ener, estep, eend),
        energy_window_unit=unit,
        ener=_as_float(data["ener"], "transport.ener") if data.get("ener") is not None else None,
        estep=_as_float(data["estep"], "transport.estep") if data.get("estep") is not None else None,
        eend=_as_float(data["eend"], "transport.eend") if data.get("eend") is not None else None,
        require_te_out=_as_bool(data.get("require_te_out", True),
                                "transport.require_te_out"),
        allow_single_channel=_as_bool(data.get("allow_single_channel", False),
                                      "transport.allow_single_channel"),
        extra_outputs=_parse_extra_outputs(data.get("extra_outputs")),
        occupation_width=_as_float(data.get("occupation_width", 0.01),
                                   "transport.occupation_width"),
        charge_mix_param=_as_float(data.get("charge_mix_param", 0.2),
                                   "transport.charge_mix_param"),
        n_max_pulay=_as_int(data.get("n_max_pulay", 10), "transport.n_max_pulay"),
    )
    if not math.isfinite(cfg.temperature) or cfg.temperature < 0 or (cfg.temperature == 0 and cfg.backend != "transiesta"):
        raise ConfigError("transport.temperature 应 > 0 K，实际是 %s" % cfg.temperature)
    if not 0.0 < cfg.occupation_width <= 0.5:
        raise ConfigError("transport.occupation_width 应在 (0, 0.5] eV，实际是 %s"
                          % cfg.occupation_width)
    if not 0.0 < cfg.charge_mix_param <= 1.0:
        raise ConfigError("transport.charge_mix_param 应在 (0, 1]，实际是 %s"
                          % cfg.charge_mix_param)
    if cfg.n_max_pulay < 0:
        raise ConfigError("transport.n_max_pulay 应 >= 0，实际是 %s" % cfg.n_max_pulay)
    return cfg


def _parse_thermal(data: Mapping[str, Any], warnings: list[str]) -> ThermalConfig:
    """``thermal`` 段。整段缺省 = 关闭，电子链行为零变化。"""
    _warn_unknown("thermal", data, _THERMAL_KEYS, warnings)

    electrode_raw = _section(data, "electrode", "thermal", warnings)
    fc_raw = _section(data, "fc", "thermal", warnings)
    negf_raw = _section(data, "negf", "thermal", warnings)
    # 2026-11-23：``_meV`` 这类键在文档里写成了 ``band_top_meV``（保留物理量写法），
    # 而本仓其它段的键一律小写。两种写法都接受：先统一小写再校验、再取值。
    negf_raw = {str(k).lower(): v for k, v in dict(negf_raw).items()}
    _warn_unknown("thermal.electrode", electrode_raw, _THERMAL_ELECTRODE_KEYS, warnings)
    _warn_unknown("thermal.fc", fc_raw, _THERMAL_FC_KEYS, warnings)
    _warn_unknown("thermal.negf", negf_raw, _THERMAL_NEGF_KEYS, warnings)

    in_plane = _as_pair(electrode_raw.get("in_plane", [3, 3]),
                        "thermal.electrode.in_plane")
    k_grid = _as_triple(fc_raw.get("k_grid", [3, 3, 1]), "thermal.fc.k_grid")
    if k_grid[2] != 1:
        raise ConfigError(
            "thermal.fc.k_grid 的 kz 必须为 1（得到 %d）：器件是带真空的二维板，"
            "沿 c 轴的 k 采样没有物理意义" % k_grid[2])

    electrode = ThermalElectrodeConfig(
        a_lattice=_as_float(electrode_raw.get("a_lattice", 4.078230),
                            "thermal.electrode.a_lattice"),
        in_plane=in_plane,
        layers=_as_int(electrode_raw.get("layers", 4), "thermal.electrode.layers",
                       minimum=2),
        device_layers=_as_int(electrode_raw.get("device_layers", 2),
                              "thermal.electrode.device_layers", minimum=1),
        vacuum=_as_float(electrode_raw.get("vacuum", 20.0),
                         "thermal.electrode.vacuum"),
        au_s_distance=_as_float(electrode_raw.get("au_s_distance", 2.35),
                                "thermal.electrode.au_s_distance"),
        gap=_as_float(electrode_raw.get("gap", 8.4), "thermal.electrode.gap"),
        element=_as_str(electrode_raw.get("element", "Au"),
                        "thermal.electrode.element"),
        anchor_element=_as_str(electrode_raw.get("anchor_element", "S"),
                               "thermal.electrode.anchor_element"),
        strip_thiol_h=_as_bool(electrode_raw.get("strip_thiol_h", True),
                               "thermal.electrode.strip_thiol_h"),
        freeze_outer_layers=_as_int(electrode_raw.get("freeze_outer_layers", 1),
                                    "thermal.electrode.freeze_outer_layers",
                                    minimum=0),
    )
    if electrode.device_layers >= electrode.layers:
        raise ConfigError(
            "thermal.electrode.device_layers(%d) 必须小于 layers(%d)："
            "至少留一层给冻结 lead" % (electrode.device_layers, electrode.layers))

    fc = ThermalFCConfig(
        displacement=_as_float(fc_raw.get("displacement", 0.01),
                               "thermal.fc.displacement"),
        cutoff=_as_float(fc_raw.get("cutoff", 10.0), "thermal.fc.cutoff"),
        array_parallel=_as_int(fc_raw.get("array_parallel", 8),
                               "thermal.fc.array_parallel", minimum=1),
        array_time_limit=_as_str(fc_raw.get("array_time_limit", "06:00:00"),
                                 "thermal.fc.array_time_limit"),
        memory=_as_str(fc_raw.get("memory", "12G"), "thermal.fc.memory"),
        sc_accuracy_forces=_as_str(fc_raw.get("sc_accuracy_forces", "5E-3"),
                                   "thermal.fc.sc_accuracy_forces"),
        relax_force_threshold=_as_str(fc_raw.get("relax_force_threshold", "1.e-2"),
                                      "thermal.fc.relax_force_threshold"),
        k_grid=k_grid,
    )
    if not 0.0 < fc.displacement <= 0.05:
        raise ConfigError(
            "thermal.fc.displacement 应在 (0, 0.05] Å，实际是 %s（0.01 Å 是标准值）"
            % fc.displacement)
    if fc.cutoff < 7.0:
        raise ConfigError(
            "thermal.fc.cutoff 应 >= 7.0 Å，实际是 %s：4 层 lead 沿 z 的跨度就是 7.06 Å，"
            "更小的截断会把相邻层的耦合切掉" % fc.cutoff)

    negf = ThermalNEGFConfig(
        eta_meV=_as_float(negf_raw.get("eta_mev", 0.05), "thermal.negf.eta_mev"),
        energy_step_meV=_as_float(negf_raw.get("energy_step_mev", 0.05),
                                  "thermal.negf.energy_step_mev"),
        band_top_meV=_as_float(negf_raw.get("band_top_mev", 19.5),
                               "thermal.negf.band_top_mev"),
    )
    if negf.eta_meV <= 0:
        raise ConfigError("thermal.negf.eta_mev 应 > 0，实际是 %s" % negf.eta_meV)
    if negf.energy_step_meV <= 0:
        raise ConfigError("thermal.negf.energy_step_mev 应 > 0，实际是 %s"
                          % negf.energy_step_meV)

    temperature = data.get("temperature")
    temperature_K = None if temperature is None else _as_float(
        temperature, "thermal.temperature")
    if temperature_K is not None and temperature_K <= 0:
        raise ConfigError("thermal.temperature 应 > 0 K，实际是 %s" % temperature_K)

    lead_fc = data.get("lead_fc")
    return ThermalConfig(
        enabled=_as_bool(data.get("enabled", False), "thermal.enabled"),
        temperature_K=temperature_K,
        lead_fc=_as_str(lead_fc, "thermal.lead_fc") if lead_fc else "auto",
        electrode=electrode,
        fc=fc,
        negf=negf,
    )


def _as_pair(raw: Any, where: str) -> tuple[int, int]:
    if not isinstance(raw, (list, tuple)) or len(raw) != 2:
        raise ConfigError("%s 应为两个整数，实际是 %r" % (where, raw))
    values = tuple(_as_int(v, "%s[%d]" % (where, i), minimum=1)
                   for i, v in enumerate(raw))
    return (values[0], values[1])


def _as_triple(raw: Any, where: str) -> tuple[int, int, int]:
    if not isinstance(raw, (list, tuple)) or len(raw) != 3:
        raise ConfigError("%s 应为三个整数，实际是 %r" % (where, raw))
    values = tuple(_as_int(v, "%s[%d]" % (where, i), minimum=1)
                   for i, v in enumerate(raw))
    return (values[0], values[1], values[2])


def _parse_extra_outputs(raw: Any) -> tuple[str, ...]:
    """``transport.extra_outputs`` -> 去重、保序的白名单元组。

    只接受已知关键字：写错的键必须在校验期报错，而不是让 FHI-aims 在远端静默
    忽略一行 ``output``（那会让"我有自旋密度"变成一个假假设）。
    """
    if raw is None:
        return ()
    if isinstance(raw, str):
        items: list[Any] = [raw]
    elif isinstance(raw, (list, tuple)):
        items = list(raw)
    else:
        raise ConfigError(
            "transport.extra_outputs 应为字符串列表（%s），实际是 %r"
            % (" / ".join(EXTRA_OUTPUTS), raw))
    result: list[str] = []
    for item in items:
        name = _as_str(item, "transport.extra_outputs[]").strip().lower()
        if name not in EXTRA_OUTPUTS:
            raise ConfigError(
                "transport.extra_outputs 只支持 %s（FHI-aims 静态 SCF 的额外输出），"
                "实际是 %r" % (" / ".join(EXTRA_OUTPUTS), item))
        if name not in result:
            result.append(name)
    return tuple(result)


def _parse_window_unit(raw: Any, ener: float, eend: float,
                       warnings: list[str]) -> str:
    """能量窗单位：显式键优先；缺省时按量级识别 0.2.x 的 Ha 绝对口径。"""
    if raw is not None:
        text = _as_str(raw, "transport.energy_window_unit")
        normalized = {"ev": "eV", "hartree": "Hartree", "ha": "Hartree"}.get(text.strip().lower())
        if normalized is None:
            raise ConfigError(
                "transport.energy_window_unit 只能是 eV（相对 E_F）或 Hartree（Ha 绝对），"
                "实际是 %r" % (text,))
        return normalized
    # 旧口径探测：0.2.x 模板写 Ha 绝对三元组，跨度 ≤0.1 Ha（≈2.7 eV）、量级 <2 Ha。
    # 新口径默认跨度 4 eV；窄到 0.1 eV 以下的 eV 窗极少见，且显式 unit 可覆盖。
    span = abs(eend - ener)
    if span <= 0.1 and max(abs(ener), abs(eend)) < 2.0:
        warnings.append(
            "transport.energy_window=%s 形如旧口径（Ha 绝对能量，跨度 %.3g ≤ 0.1）：按 "
            "Hartree 解释。新口径是 eV 且相对费米能级（默认 E_F ±2 eV），"
            "建议写 energy_window_unit: eV 显式声明" % ([ener, eend], span))
        return "Hartree"
    return "eV"


def _parse_cluster(data: Mapping[str, Any], warnings: list[str]) -> ClusterConfig:
    _warn_unknown("cluster", data, _CLUSTER_KEYS, warnings)
    env_raw = data.get("env") or {}
    if not isinstance(env_raw, Mapping):
        raise ConfigError("cluster.env 应为 {键: 值} 映射")
    env: dict[str, str] = {}
    for key, value in env_raw.items():
        if value is None:
            continue
        env[str(key)] = _expand(str(value))
    cfg = ClusterConfig(
        host=(_as_str(data["host"], "cluster.host") if data.get("host") else None),
        user=(_as_str(data["user"], "cluster.user") if data.get("user") else None),
        port=_as_int(data.get("port", 22), "cluster.port", minimum=1),
        key=(_as_str(data["key"], "cluster.key") if data.get("key") else None),
        remote_root=data.get("remote_root"),  # 在 SSH 登录节点解析 ~ 和 $HOME
        partition=(_as_str(data["partition"], "cluster.partition") if data.get("partition") else None),
        ntasks=(_as_int(data["ntasks"], "cluster.ntasks", minimum=1) if data.get("ntasks") is not None else None),
        gres=(_as_str(data["gres"], "cluster.gres") if data.get("gres") else None),
        account=(_as_str(data["account"], "cluster.account") if data.get("account") else None),
        qos=(_as_str(data["qos"], "cluster.qos") if data.get("qos") else None),
        species_dir=(_as_str(data["species_dir"], "cluster.species_dir") if data.get("species_dir") else None),
        time_limit=(_as_str(data["time_limit"], "cluster.time_limit") if data.get("time_limit") else None),
        mem=(_as_str(data["mem"], "cluster.mem") if data.get("mem") else None),
        poll_seconds=_as_float(data.get("poll_seconds", 60.0), "cluster.poll_seconds"),
        timeout_hours=_as_float(data.get("timeout_hours", 6.0), "cluster.timeout_hours"),
        env=env,
    )
    if cfg.port < 1 or cfg.port > 65535:
        raise ConfigError("cluster.port 超出范围: %d" % cfg.port)
    return cfg


def _parse_paths(data: Mapping[str, Any], warnings: list[str]) -> PathsConfig:
    _warn_unknown("paths", data, _PATHS_KEYS, warnings)
    species_dir = data.get("species_dir")
    species_defaults = data.get("species_defaults")
    if species_dir is not None and species_defaults is not None:
        raise ConfigError("paths.species_dir 与 paths.species_defaults 不能同时给（前者是兼容别名）")
    if species_dir is not None:
        warnings.append("paths.species_dir 已并入 paths.species_defaults（旧键仍生效）")
        species_defaults = species_dir
    return PathsConfig(
        templates=(_as_str(data["templates"], "paths.templates") if data.get("templates") else None),
        species_defaults=(_as_str(species_defaults, "paths.species_defaults") if species_defaults else None),
    )


def _parse_review(data: Mapping[str, Any], warnings: list[str]) -> ReviewConfig:
    _warn_unknown("review", data, _REVIEW_KEYS, warnings)
    accept = data.get("accept_deviation")
    if data.get("surface_check") is not None and data.get("connection_check") is None:
        warnings.append("review.surface_check 已更名为 review.connection_check（旧键仍生效）")
    if data.get("surface_check_deviation") is not None and data.get("connection_check_deviation") is None:
        warnings.append("review.surface_check_deviation 已更名为 review.connection_check_deviation（旧键仍生效）")
    if data.get("auto_surface_check") is not None and data.get("auto_connection_check") is None:
        warnings.append("review.auto_surface_check 已更名为 review.auto_connection_check（旧键仍生效）")
    if data.get("require_a2") is not None and data.get("connection_check") is None:
        warnings.append("review.require_a2 已更名为 review.connection_check（旧键仍生效）")
    if data.get("a2_deviation") is not None and data.get("connection_check_deviation") is None:
        warnings.append("review.a2_deviation 已更名为 review.connection_check_deviation（旧键仍生效）")
    check = data.get("connection_check", data.get("surface_check", data.get("require_a2", True)))
    deviation = data.get(
        "connection_check_deviation",
        data.get("surface_check_deviation", data.get("a2_deviation", 0.15)),
    )
    auto = data.get("auto_connection_check", data.get("auto_surface_check", True))
    return ReviewConfig(
        connection_check=_as_bool(check, "review.connection_check"),
        connection_check_deviation=_as_float(deviation, "review.connection_check_deviation"),
        auto_connection_check=_as_bool(auto, "review.auto_connection_check"),
        accept_deviation=_as_float(accept, "review.accept_deviation") if accept is not None else None,
    )


def load_config(path: str | Path | None = None) -> RunConfig:
    """读取并校验 ``smj.yaml``。

    未知字段只产生 warning（向前兼容），非法值抛 :class:`ConfigError`。
    """
    cfg_path = Path(path or DEFAULT_CONFIG_NAME)
    if not cfg_path.is_absolute():
        cfg_path = (Path.cwd() / cfg_path).resolve()
    if not cfg_path.exists():
        raise ConfigError("找不到配置文件: %s\n提示: 先运行 `smj init <目录>` 生成 smj.yaml" % cfg_path)
    try:
        raw = yaml.safe_load(cfg_path.read_text(encoding="utf-8")) or {}
    except yaml.YAMLError as exc:
        raise ConfigError("YAML 解析失败 (%s): %s" % (cfg_path, exc)) from exc
    if not isinstance(raw, Mapping):
        raise ConfigError("配置文件顶层应为映射: %s" % cfg_path)

    warnings: list[str] = []
    unknown_sections = sorted(set(raw) - _KNOWN_SECTIONS)
    if unknown_sections:
        warnings.append("未知配置段，已忽略: %s" % ", ".join(unknown_sections))

    run = _parse_run(_section(raw, "run", "config", warnings), warnings)
    config = RunConfig(
        run=run,
        molecule=_parse_molecule(_section(raw, "molecule", "config", warnings), warnings),
        electrode=_parse_electrode(_section(raw, "electrode", "config", warnings), warnings),
        transport=_parse_transport(_section(raw, "transport", "config", warnings), warnings),
        siesta=dict(_section(raw, "siesta", "config", warnings)),
        thermal=_parse_thermal(_section(raw, "thermal", "config", warnings), warnings),
        cluster=_parse_cluster(_section(raw, "cluster", "config", warnings), warnings),
        paths=_parse_paths(_section(raw, "paths", "config", warnings), warnings),
        review=_parse_review(_section(raw, "review", "config", warnings), warnings),
        source=cfg_path,
        warnings=warnings,
    )

    # 跨段校验
    if config.mode == "slurm-ssh":
        if not config.cluster.host:
            raise ConfigError("run.mode=slurm-ssh 时必须设置 cluster.host（SSH 别名/主机名）")
        if not config.cluster.partition:
            raise ConfigError("run.mode=slurm-ssh 时必须设置 cluster.partition（Slurm 分区）")
    if config.transport.backend == "transiesta":
        from .siesta.stages import validate_settings
        from .siesta.stages import DEFAULTS
        known = set(DEFAULTS) | {"electrode", "contact", "pseudopotentials", "basis",
                                 "tbtrans_kgrid", "electrode_reuse", "device_mixing_weight",
                                 "reuse_relax_from", "reuse_scf_from"}
        unknown = config.siesta.keys() - known
        if unknown:
            raise ConfigError(f"未知 siesta 配置字段: {', '.join(sorted(unknown))}")
        validate_settings(config.siesta)
    if config.siesta and config.transport.backend != "transiesta":
        raise ConfigError("siesta 配置需要 transport.backend: transiesta")
    if config.transport.backend == "transiesta" and any(getattr(config.transport, n) is not None for n in ("ener", "eend", "estep")):
        raise ConfigError("SIESTA 能窗请用相对化学势的 transport.energy_window")
    if config.transport.backend == "transiesta" and config.transport.energy_window_unit != "eV":
        raise ConfigError("TBtrans 能窗采用相对共同化学势的 eV")
    if config.mode == "local" and config.transport.backend == "aitranss" and not config.cluster.species_dir and not config.paths.species_defaults:
        # 只有本地真跑 AIMS 时才需要 species；这里只给 warning，不阻断 plan/prepare。
        warnings.append("未配置 cluster.species_dir / paths.species_defaults：本地真跑 AIMS 阶段会失败")

    # 本地路径统一按 yaml 所在目录解析（remote_root 是远端 POSIX 路径，保持原样）
    if config.cluster.key:
        config.cluster.key = str(config.resolve(config.cluster.key))
    if config.cluster.species_dir:
        config.cluster.species_dir = str(config.resolve(config.cluster.species_dir))
    if config.paths.templates:
        config.paths.templates = str(config.resolve(config.paths.templates))
    if config.paths.species_defaults:
        config.paths.species_defaults = str(config.resolve(config.paths.species_defaults))
    if config.electrode.library_dir:
        config.electrode.library_dir = str(config.resolve(config.electrode.library_dir))
    # "auto" 是关键字，不能当相对路径解析（否则会被拼成 yaml 目录下的
    # 一个名叫 auto 的文件）。
    if config.thermal.lead_fc and config.thermal.lead_fc != "auto":
        config.thermal.lead_fc = str(config.resolve(config.thermal.lead_fc))
    return config


# --------------------------------------------------------------------------- #
# smj init 模板
# --------------------------------------------------------------------------- #
DEFAULT_CONFIG_TEMPLATE = """\
# smj-toolkit 运行配置（smj init 生成）
# 路径均相对本文件所在目录；${{VAR}} / ~ 会先展开。
run:
  name: {name}
  dir: .
  mode: {mode}            # dry-run | local | slurm-ssh

molecule:
  xyz: {xyz}
  charge: 0
  multiplicity: 1
  initial_moments: null     # 开壳层初始磁矩 {{1-based 原子序号: μB}}，如 {{7: 1.0}}；null = 不写 initial_moment
  fixed_spin_moment: true   # 开壳层结是否钉住总自旋（M−1）；false = 解除约束的对照解
  broken_symmetry: false    # M=1 的 BS 单重态：spin collinear + 反号 initial_moments，不写 fixed_spin_moment（需 initial_moments 合计 0）
  opt_engine: aims          # aims = 与 tip/SCF 同为 FHI-aims PBE/NAO（推荐）；gaussian = 旧口径
  opt_functional: pbe       # 仅 opt_engine=aims
  opt_basis: light          # 仅 opt_engine=aims；实际 species 由 cluster.species_dir 决定
  opt_force_threshold: 1.e-2
  method: b3lyp             # 仅 opt_engine=gaussian
  basis: 6-31g(d)           # 仅 opt_engine=gaussian
  stable_check: false       # 仅 gaussian：加 Stable=Opt（先做波函数稳定性检查再优化）
  family: null              # 判据标定族（AZ/BDT/…）；未登记体系按"全部判据适用"处理

electrode:
  source: library           # 构建方式：library = FHI-aims electrodes.library 预定义团簇；builtin = 本地三角岛
  element: Au               # 电极金属（fcc(111) 单晶表面）：Au / Ag / Cu / Ni / Pt / Pd / Al；tip 同元素
  a_lattice: null           # fcc 晶格常数(Å)；null = 用内置登记值，未登记元素必须显式给出
  library_dir: null         # electrodes.library 目录；null = 本机常见路径 / slurm-ssh 远端探测
  library_cluster: fcc111.au.clusters.37_atoms/fcc111.cluster.natoms37.x1.aims
  # au_per_electrode: 30    # 仅 builtin：每侧电极金属原子数（null = 完整 60/侧）；library 模式忽略该键；按需取消注释
  left_anchor: S
  right_anchor: S
  au_tip_distance: 2.35
  anchor_au_s_target: null  # null = 只按 au_tip_distance 放置，锚–金属距离由全原子弛豫决定
  pin_contact: false        # false = aims_tip_relax 全原子自由弛豫（推荐）；true 仅用于用户明确要求的 pinned 对照/距离扫描
  tip_direction: axis       # axis = 沿两锚原子连线外推（默认）；lone_pair = 沿锚原子孤对外推
                            # （S···S 连线与 S–C(Me) 近共线的间位取代 SMe 分子才需要）
  strip_thiol_h: null       # null = 仅当锚定原子为 S 时自动去质子（其它锚保留 H）；true = 强制去；false = 不去
  vdw: true
  relax: false              # false = 电极扩展后几何冻结，SCF/输运不再弛豫（推荐）
  relax_force_threshold: 1.e-2
  relax_freeze: backplane   # 仅 relax: true 时生效；backplane = 冻背层，electrodes = 冻全部电极金属
  nlayers: null             # Aitranss 中心区层数；null = 按库文件头部 -nlayers 建议

transport:
  temperature: 300.0
  energy_window: [-2.0, 0.005, 2.0]  # [起点, 步长, 终点] eV，相对 E_F（默认 E_F ±2 eV，步长 5 meV = 801 点）
  energy_window_unit: eV    # eV（默认）= 相对费米能级；Hartree = 旧口径（Ha 绝对能量三元组）
  require_te_out: true      # false = 允许缺 te.out 的外部/历史 run（显式豁免，不建议）
  # extra_outputs: [hirshfeld, spin_density]   # 开壳层证据：Hirshfeld 磁矩 + 自旋密度 cube
  occupation_width: 0.01    # 静态 SCF 的 gaussian 展宽(eV)；金属电极结建议 0.05（0.01 常令 fixed_spin_moment 震荡）
  charge_mix_param: 0.2     # Pulay 混合步长；金属 + 固定自旋建议 0.05
  n_max_pulay: 10           # Pulay 历史长度；难收敛的金属结可加到 20-30

cluster:
  host: null                # ~/.ssh/config 别名，如 cluster
  user: null
  port: 22
  key: null
  remote_root: ~/smj-runs   # 远端登录用户的工作根，可省略
  partition: null           # 如 normal
  ntasks: null              # 覆盖并行度；null = 按体系/阶段自动估算
  gres: null
  account: null
  species_dir: null         # FHI-aims species_defaults 目录
  poll_seconds: 60
  timeout_hours: 6
  env:                      # 集群部署路径（对应 AgenticDFT 的 SMJ_* 环境变量）
    # g16root: /opt/gaussian/g16
    # g16_module: Gaussian/g16.c02
    # aims_bin: /opt/fhi-aims/bin/aims.x
    # aims_ld_module: oneapi/2024
    # aitranss_bin: /opt/aitranss/aitranss.x
    # mkl_lib: /opt/intel/oneapi/mkl/latest/lib/intel64

paths:
  templates: null
  species_defaults: null

review:
  connection_check: true          # 连接核对（核对每侧三个表面原子与 tcontrol 是否一致）
  connection_check_deviation: 0.15
  auto_connection_check: true     # 默认自动登记（electrode_build 采集后自动审计并写记录）
"""


def default_config_text(name: str = "my-run", xyz: str = "molecule.xyz", mode: str = "dry-run") -> str:
    if mode not in MODES:
        raise ConfigError("未知 mode: %r（可选 %s）" % (mode, "/".join(MODES)))
    if not _RUN_NAME_RE.match(name):
        raise ConfigError("run 名只允许字母/数字/._-: %r" % name)
    if mode == "slurm-ssh":
        validate_run_name(name)
    return DEFAULT_CONFIG_TEMPLATE.format(name=name, xyz=xyz, mode=mode)


def write_default_config(path: str | Path, *, name: str = "my-run", xyz: str = "molecule.xyz",
                         mode: str = "dry-run", overwrite: bool = False) -> Path:
    """写一份带注释的模板配置（``smj init`` 使用）。"""
    target = Path(path)
    if target.is_dir():
        target = target / DEFAULT_CONFIG_NAME
    if target.exists() and not overwrite:
        raise ConfigError("配置文件已存在: %s（加 --force 覆盖）" % target)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(default_config_text(name=name, xyz=xyz, mode=mode), encoding="utf-8", newline="\n")
    return target
