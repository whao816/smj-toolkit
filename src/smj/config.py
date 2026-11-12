"""``smj.yaml`` 的加载、校验与路径解析。

一个配置文件描述**一条 run**（分子 + 电极 + 输运 + 集群）。相对路径一律相对
配置文件所在目录解析，所以 ``smj.yaml`` 可以和 ``molecule.xyz`` 放在同一个目录里。

与 AgenticDFT 的 ``config.yaml`` 不共享代码：独立包只保留计算所需的字段，字段语义
对齐 ``tools/junction`` 输入写出器的 ``env`` 约定（``g16root`` / ``aims_bin`` /
``aitranss_bin`` … 都可以在 ``cluster.env`` 里显式给出）。
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import yaml

from .errors import ConfigError
from .electrode_library import DEFAULT_CLUSTER

DEFAULT_CONFIG_NAME = "smj.yaml"
MODES = ("dry-run", "local", "slurm-ssh")
_RUN_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


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
    """分子（Gaussian 优化阶段）参数。"""

    xyz: str | None = None
    charge: int = 0
    multiplicity: int = 1
    method: str = "b3lyp"
    basis: str = "6-31g(d)"


@dataclass
class ElectrodeConfig:
    """电极 / 结几何参数（au60_build 与 aims_tip_relax、aitranss_te 共用）。"""

    au_per_electrode: int | None = 30     # 仅 builtin：每侧 Au 数；None/30 = 完整 Au60
    left_anchor: str = "S"
    right_anchor: str = "S"
    au_tip_distance: float = 2.35
    anchor_au_s_target: float = 2.40
    strip_thiol_h: bool = True
    vdw: bool = True
    nlayers: int | None = None            # Aitranss 中心区层数；None = 按库文件头部建议
    source: str = "library"               # library = FHI-aims 预定义团簇；builtin = 原三角岛
    library_dir: str | None = None        # electrodes.library 目录；None = 常见路径/远端探测
    library_cluster: str = DEFAULT_CLUSTER
    relax: bool = True                    # au_relax：固定背层的几何弛豫（手册 §5.3.1 要求）
    relax_force_threshold: str = "1.e-2"  # BFGS 力收敛阈值（Hartree/Bohr）
    #: au_relax 冻结谁：backplane = 只冻结两侧背层（手册口径）；electrodes = 冻结全部 Au，只让分子弛豫
    relax_freeze: str = "backplane"


@dataclass
class TransportConfig:
    """透射谱 / Landauer 分析参数。"""

    temperature: float = 300.0
    energy_window: tuple[float, float, float] = (-2.5, 0.02, 2.5)  # (ener, estep, eend)
    ener: float | None = None             # 显式覆盖 energy_window 的分量
    estep: float | None = None
    eend: float | None = None

    def resolved_window(self) -> tuple[float, float, float]:
        ener, estep, eend = self.energy_window
        if self.ener is not None:
            ener = self.ener
        if self.estep is not None:
            estep = self.estep
        if self.eend is not None:
            eend = self.eend
        return float(ener), float(estep), float(eend)


@dataclass
class ClusterConfig:
    """Slurm / SSH 集群与外部程序部署路径。"""

    host: str | None = None               # ~/.ssh/config 别名或 hostname
    user: str | None = None
    port: int = 22
    key: str | None = None                # 私钥路径（None = ssh 默认）
    remote_root: str | None = None        # 远端 run 根，如 /scratch/$USER/smj
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
    """结构复核 / 报告（附带功能）的开关。"""

    require_a2: bool = True               # 默认强制 A2 表面确认（本地补丁 F4：红线默认生效）
    a2_deviation: float = 0.15            # A2 对比的相对偏差门槛
    accept_deviation: float | None = None # 覆盖 A2 默认门槛


@dataclass
class RunConfig:
    """一条 smj run 的完整配置。"""

    run: RunInfo
    molecule: MoleculeConfig = field(default_factory=MoleculeConfig)
    electrode: ElectrodeConfig = field(default_factory=ElectrodeConfig)
    transport: TransportConfig = field(default_factory=TransportConfig)
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
_KNOWN_SECTIONS = {"run", "molecule", "electrode", "transport", "cluster", "paths", "review"}

_RUN_KEYS = {"name", "dir", "mode", "clean", "keep_going"}
_MOLECULE_KEYS = {"xyz", "charge", "multiplicity", "method", "basis"}
_ELECTRODE_KEYS = {
    "au_per_electrode", "left_anchor", "right_anchor", "au_tip_distance",
    "anchor_au_s_target", "strip_thiol_h", "vdw", "nlayers",
    "source", "library_dir", "library_cluster", "relax", "relax_force_threshold",
    "relax_freeze",
}
_TRANSPORT_KEYS = {"temperature", "energy_window", "ener", "estep", "eend"}
_CLUSTER_KEYS = {
    "host", "user", "port", "key", "remote_root", "partition", "ntasks", "gres",
    "account", "qos", "species_dir", "poll_seconds", "timeout_hours", "env",
    "time_limit", "mem",
}
_PATHS_KEYS = {"templates", "species_defaults"}
_REVIEW_KEYS = {"require_a2", "a2_deviation", "accept_deviation"}


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


def _parse_molecule(data: Mapping[str, Any], warnings: list[str]) -> MoleculeConfig:
    _warn_unknown("molecule", data, _MOLECULE_KEYS, warnings)
    xyz = data.get("xyz")
    cfg = MoleculeConfig(
        xyz=_as_str(xyz, "molecule.xyz", allow_empty=True) or None,
        charge=_as_int(data.get("charge", 0), "molecule.charge"),
        multiplicity=_as_int(data.get("multiplicity", 1), "molecule.multiplicity", minimum=1),
        method=_as_str(data.get("method") or "b3lyp", "molecule.method"),
        basis=_as_str(data.get("basis") or "6-31g(d)", "molecule.basis"),
    )
    return cfg


def _parse_electrode(data: Mapping[str, Any], warnings: list[str]) -> ElectrodeConfig:
    _warn_unknown("electrode", data, _ELECTRODE_KEYS, warnings)
    source = _as_str(data.get("source") or "library", "electrode.source")
    if source not in ("library", "builtin"):
        raise ConfigError("electrode.source 应为 library/builtin 之一，实际是 %r" % source)
    au = data.get("au_per_electrode", 30)
    au_value = None if au is None else _as_int(au, "electrode.au_per_electrode", minimum=1)
    if au_value is not None and au_value > 30:
        raise ConfigError("electrode.au_per_electrode 最多 30（Au60 每侧 30 个原子），实际是 %d" % au_value)
    target = _as_float(data.get("anchor_au_s_target", 2.40), "electrode.anchor_au_s_target")
    if not (1.5 <= target <= 4.0):
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
        strip_thiol_h=_as_bool(data.get("strip_thiol_h", True), "electrode.strip_thiol_h"),
        vdw=_as_bool(data.get("vdw", True), "electrode.vdw"),
        nlayers=nlayers,
        source=source,
        library_dir=library_dir,
        library_cluster=cluster_rel,
        relax=_as_bool(data.get("relax", True), "electrode.relax"),
        relax_force_threshold=threshold,
        relax_freeze=freeze,
    )
    if not (1.0 <= cfg.au_tip_distance <= 5.0):
        raise ConfigError("electrode.au_tip_distance 应在 1.0–5.0 Å，实际是 %s" % cfg.au_tip_distance)
    return cfg


def _parse_transport(data: Mapping[str, Any], warnings: list[str]) -> TransportConfig:
    _warn_unknown("transport", data, _TRANSPORT_KEYS, warnings)
    window = data.get("energy_window", [-2.5, 0.02, 2.5])
    if not isinstance(window, (list, tuple)) or len(window) != 3:
        raise ConfigError("transport.energy_window 应为 [ener, estep, eend] 三个数，实际是 %r" % (window,))
    ener = _as_float(window[0], "transport.energy_window[0]")
    estep = _as_float(window[1], "transport.energy_window[1]")
    eend = _as_float(window[2], "transport.energy_window[2]")
    if estep <= 0 or eend <= ener:
        raise ConfigError("transport.energy_window 要求 estep > 0 且 eend > ener（当前 %r）" % (window,))
    cfg = TransportConfig(
        temperature=_as_float(data.get("temperature", 300.0), "transport.temperature"),
        energy_window=(ener, estep, eend),
        ener=_as_float(data["ener"], "transport.ener") if data.get("ener") is not None else None,
        estep=_as_float(data["estep"], "transport.estep") if data.get("estep") is not None else None,
        eend=_as_float(data["eend"], "transport.eend") if data.get("eend") is not None else None,
    )
    if cfg.temperature <= 0:
        raise ConfigError("transport.temperature 应 > 0 K，实际是 %s" % cfg.temperature)
    return cfg


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
        remote_root=(_as_str(data["remote_root"], "cluster.remote_root") if data.get("remote_root") else None),
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
    return PathsConfig(
        templates=(_as_str(data["templates"], "paths.templates") if data.get("templates") else None),
        species_defaults=(_as_str(data["species_defaults"], "paths.species_defaults") if data.get("species_defaults") else None),
    )


def _parse_review(data: Mapping[str, Any], warnings: list[str]) -> ReviewConfig:
    _warn_unknown("review", data, _REVIEW_KEYS, warnings)
    accept = data.get("accept_deviation")
    return ReviewConfig(
        require_a2=_as_bool(data.get("require_a2", True), "review.require_a2"),
        a2_deviation=_as_float(data.get("a2_deviation", 0.15), "review.a2_deviation"),
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
        if not config.cluster.remote_root:
            raise ConfigError("run.mode=slurm-ssh 时必须设置 cluster.remote_root（远端 run 根目录）")
        if not config.cluster.partition:
            raise ConfigError("run.mode=slurm-ssh 时必须设置 cluster.partition（Slurm 分区）")
    if config.mode == "local" and not config.cluster.species_dir and not config.paths.species_defaults:
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
  method: b3lyp
  basis: 6-31g(d)

electrode:
  source: library           # library = FHI-aims electrodes.library 预定义团簇；builtin = 自建三角岛
  library_dir: null         # electrodes.library 目录；null = 本机常见路径 / slurm-ssh 远端探测
  library_cluster: fcc111.au.clusters.37_atoms/fcc111.cluster.natoms37.x1.aims
  au_per_electrode: 30      # 仅 builtin：每侧 Au 原子数（最大 30 = 完整 Au60）
  left_anchor: S
  right_anchor: S
  au_tip_distance: 2.35
  anchor_au_s_target: 2.40
  strip_thiol_h: true
  vdw: true
  relax: true                     # au_relax：固定两侧背层做几何弛豫（手册 §5.3.1）
  relax_force_threshold: 1.e-2
  relax_freeze: backplane   # backplane = 只冻结两侧背层；electrodes = 冻结全部 Au，只让分子弛豫
  nlayers: null             # Aitranss 中心区层数；null = 按库文件头部 -nlayers 建议

transport:
  temperature: 300.0
  energy_window: [-2.5, 0.02, 2.5]   # [ener, estep, eend] (eV)

cluster:
  host: null                # ~/.ssh/config 别名，如 cluster
  user: null
  port: 22
  key: null
  remote_root: null         # 如 /scratch/$USER/smj
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
  require_a2: true
  a2_deviation: 0.15
"""


def default_config_text(name: str = "my-run", xyz: str = "molecule.xyz", mode: str = "dry-run") -> str:
    if mode not in MODES:
        raise ConfigError("未知 mode: %r（可选 %s）" % (mode, "/".join(MODES)))
    if not _RUN_NAME_RE.match(name):
        raise ConfigError("run 名只允许字母/数字/._-: %r" % name)
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
