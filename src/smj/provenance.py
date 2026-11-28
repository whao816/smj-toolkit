"""运行溯源：``run_spec.json``、协议指纹与阶段间哈希闭环。

为什么需要这个模块（2026-11 架构评审的 P0）：

- ``state.json`` 只记"哪个阶段到了什么状态"，不记**这次 run 到底用了什么参数**。
  模板里的默认值、用户没显式写的键、软件版本都不在档案里，事后无法复核；
- A2 只绑定 ``geometry.in`` 的 sha256，**配置**（方法/基组/电极距离/输运窗/温度）
  在 run 开始后被改掉，没有任何东西会发现；
- 阶段之间只检查"文件存在"，上游几何被替换、或 SCF 用的不是 A2 确认的那份
  几何，链路本身无法证明。

本模块给出三样东西：

1. :func:`ensure_run_spec` —— 在第一次 ``prepare`` 时把**完整解析后的配置**、
   分子/电极库哈希、代码版本、软件路径与本地可得的二进制哈希写成
   ``<run>/run_spec.json``（schema ``smj.run_spec.v1``）。之后配置漂移默认拒绝，
   ``--force`` 才重绑，并把旧 spec 留成 ``run_spec.previous.json``；
2. :func:`protocol_sha256` —— 只覆盖**会改变物理结果**的配置子集（分子/电极/
   输运），用于 review 的漂移门；
3. 输入/输出哈希与 :func:`verify_parent_geometry` —— 让"这一阶段的几何是不是
   上游/A2 确认的那一份"变成可验证的断言，而不是默认假设。

设计约束：本模块不导入阶段实现，不发起任何远端命令；远端二进制只记录路径，
本地可读的才记 sha256。所有失败都是可解释的 :class:`ConfigError` /
:class:`StageBlocked`，不静默降级。
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

from .errors import ConfigError, StageBlocked
from .state import now_iso

RUN_SPEC_SCHEMA = "smj.run_spec.v1"
RUN_SPEC_FILE = "run_spec.json"
RUN_SPEC_PREVIOUS = "run_spec.previous.json"

#: 协议指纹的"形状版本"。新增/移除只影响记录方式的配置键时递增；老 run 的
#: ``protocol_sha256`` 会因此不同，但物理参数没变，不该判成漂移。
PROTOCOL_SCHEMA = 2

#: 单文件哈希上限。超过这个大小只记 path/size，不把 GB 级 ``mos.aims`` 读进内存。
MAX_HASH_BYTES = 128 * 1024 * 1024

#: legacy run 没有 run_spec.json 时，collect 观测到的远端软件版本落到这里。
SOFTWARE_OBSERVED_FILE = "software_observed.json"

_AIMS_VERSION_RE = re.compile(r"FHI-aims version\s*:\s*(\S+)")
_AIMS_UUID_RE = re.compile(r"aims_uuid\s*:\s*(\S+)")
_AIMS_COMMIT_RE = re.compile(r"Commit number\s*:\s*(\S+)")
_AITRANSS_VERSION_RE = re.compile(r"invoking aitranss \(version:\s*([^)]+)\)")
_AITRANSS_COMPILED_RE = re.compile(r"compiled on ([^\n]+)")


# --------------------------------------------------------------------------- #
# 基础哈希
# --------------------------------------------------------------------------- #
def sha256_file(path: str | Path, *, max_bytes: int = MAX_HASH_BYTES) -> str | None:
    """文件 sha256（十六进制）；缺失/不可读/超限返回 ``None``。"""
    p = Path(path)
    try:
        if not p.is_file():
            return None
        if p.stat().st_size > max_bytes:
            return None
        digest = hashlib.sha256()
        with p.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()
    except OSError:
        return None


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def structural_fingerprint(path: str | Path) -> str | None:
    """几何的**结构**指纹：只哈希原子行（元素 + 坐标），忽略 initial_charge /
    initial_moment 等子 tag。

    为什么要单独定义一个：FHI-aims 要求带电/自旋极化体系把 per-atom
    ``initial_charge``/``initial_moment`` 写进 geometry.in，同一个结构的字节
    sha256 会随这些子 tag 变化。阶段间核对「是不是同一份几何」应当看结构本身，
    而不是看这些初始猜测量。缺失或读不到返回 ``None``。
    """
    p = Path(path)
    try:
        if not p.is_file():
            return None
        atoms: list[str] = []
        for raw in p.read_text(encoding="utf-8", errors="replace").splitlines():
            fields = raw.split()
            if not fields:
                continue
            head = fields[0]
            if head in ("atom", "lattice_vector"):
                # 坐标与元素做规范化（空白折叠），子 tag 行（initial_*）不参与。
                atoms.append(" ".join(fields))
        if not atoms:
            return None
        return sha256_text("\n".join(atoms))
    except OSError:
        return None


def canonical_json(payload: Any) -> str:
    """稳定 JSON：排序键、紧凑分隔符、不允许 NaN。"""
    return json.dumps(payload, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False, allow_nan=False)


def hash_stage_inputs(paths: Iterable[str | Path]) -> dict[str, str]:
    """把 prepare 的输入文件哈希成 ``{名字: sha256}``（缺失/超限的跳过）。

    键用文件名；同名不同内容时改用 ``父目录名/文件名``，保证不会互相覆盖。
    """
    out: dict[str, str] = {}
    for raw in paths or ():
        path = Path(raw)
        digest = sha256_file(path)
        if digest is None:
            continue
        key = path.name
        if key in out and out[key] != digest:
            key = "%s/%s" % (path.parent.name, path.name)
        out[key] = digest
    return out


def hash_stage_outputs(stage_dir: str | Path, artifacts: Iterable[str], *,
                       missing: list[str] | None = None) -> dict[str, str]:
    """把 collect 声明的产物哈希成 ``{相对名: sha256}``。

    缺件或超过 ``MAX_HASH_BYTES`` 的文件无法哈希：旧实现静默省略，``outputs_sha256``
    看起来"全都核过"，实际有洞。现在把名字追加到 ``missing``（由调用方写进
    state 的 ``output_hash_missing`` 与 note），让缺失可见。
    """
    root = Path(stage_dir)
    out: dict[str, str] = {}
    for raw in artifacts or ():
        name = str(raw)
        path = Path(name)
        if not path.is_absolute():
            path = root / path
        digest = sha256_file(path)
        if digest is not None:
            out[name] = digest
        elif missing is not None:
            missing.append(name)
    return out


# --------------------------------------------------------------------------- #
# 配置 <-> JSON
# --------------------------------------------------------------------------- #
def _jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if is_dataclass(value) and not isinstance(value, type):
        return {key: _jsonable(item) for key, item in asdict(value).items()}
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(item) for item in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def config_to_dict(config: Any) -> dict[str, Any]:
    """完整解析配置（含 dataclass 默认值）-> JSON 可写 dict。

    ``source`` / ``warnings`` 属于加载过程元数据，不进入配置快照。
    """
    data = {
        "run": _jsonable(config.run),
        "molecule": _jsonable(config.molecule),
        "electrode": _jsonable(config.electrode),
        "transport": _jsonable(config.transport),
        "cluster": _jsonable(config.cluster),
        "paths": _jsonable(config.paths),
        "review": _jsonable(config.review),
    }
    if config.transport.backend == "transiesta":
        data["siesta"] = _jsonable(config.siesta)
    else:
        data["transport"].pop("backend", None)
    # 输运窗同时保留"用户写的三元组"和"实际生效的窗口"（ener/estep/eend 覆盖）。
    try:
        data["transport"]["energy_window"] = list(config.transport.energy_window)
        data["transport"]["resolved_window"] = list(config.transport.resolved_window())
        # 2026-11-20 复核 P2-1：resolved_window 恒为 Ha 绝对能量，标注单位避免误读。
        data["transport"]["resolved_window_unit"] = "Hartree"
        if config.transport.backend == "transiesta":
            data["transport"]["resolved_window"] = list(config.transport.relative_window())
            data["transport"]["resolved_window_unit"] = "eV relative to common chemical potential"
    except Exception:  # noqa: BLE001 - 兼容手写/旧配置对象
        pass
    # 分子路径按绝对路径入档；实际文件哈希在 spec["molecule"] 里。
    try:
        xyz = config.molecule_xyz()
        if xyz is not None:
            data["molecule"]["xyz"] = str(Path(xyz))
    except Exception:  # noqa: BLE001
        pass
    return data


def _electrode_library_info(config: Any) -> dict[str, Any]:
    """电极库文件（FHI-aims ``electrodes.library`` 的预定义团簇）的定位与哈希。

    远端库在本机不可读时只记 ``exists=false``；不假装知道内容。该哈希进入
    protocol 指纹，避免"同名团簇换了库文件"被当成同一物理输入。
    """
    electrode = getattr(config, "electrode", None)
    if electrode is None or getattr(electrode, "source", "builtin") != "library":
        return {"path": None, "exists": False, "sha256": None}
    raw_dir = (getattr(electrode, "library_dir", None)
               or os.environ.get("SMJ_ELECTRODE_LIBRARY"))
    rel = getattr(electrode, "library_cluster", None)
    if not raw_dir or not rel:
        return {"path": None, "exists": False, "sha256": None}
    path = Path(os.path.expandvars(os.path.expanduser(str(raw_dir)))) / str(rel)
    return {
        "path": str(path),
        "exists": path.is_file(),
        "sha256": sha256_file(path),
    }


#: 这些 ``cluster.env`` 键决定"实际用哪套基组/哪个可执行文件"，会改变算出来的数；
#: 其余环境键（account/qos/partition 等）只影响调度，不进执行面指纹。
_EXECUTION_ENV_KEYS = (
    "aims_bin", "aitranss_bin", "aitransss_bin",
    "g16root", "G16ROOT", "MKL_LIB", "AIMS_SPECIES_DEFAULTS",
)

#: 目录内容哈希的读取上限：避免把整个 ``species_defaults`` 树读进内存/拖慢 prepare。
_EXECUTION_TREE_LIMIT = 128
_EXECUTION_LIST_LIMIT = 5000


def _resolve_user_path(raw: Any) -> Path | None:
    if raw is None or not str(raw).strip():
        return None
    return Path(os.path.expandvars(os.path.expanduser(str(raw))))


def _executable_probe(raw: Any) -> dict[str, Any]:
    """可执行文件：路径 + 存在性 + 本地可读时的 sha256（远端路径只记路径）。"""
    path = _resolve_user_path(raw)
    if path is None:
        return {"path": None, "exists": False, "sha256": None}
    exists = path.is_file()
    return {
        "path": str(path),
        "exists": exists,
        "sha256": sha256_file(path) if exists else None,
    }


def _tree_probe(raw: Any) -> dict[str, Any]:
    """目录（species defaults / templates）的路径、文件清单与内容指纹。

    清单哈希覆盖所有文件的相对路径与大小（最多 ``_EXECUTION_LIST_LIMIT`` 个），
    内容哈希只读前 ``_EXECUTION_TREE_LIMIT`` 个文件；``truncated`` 明确标注未
    覆盖的范围——不把"只哈希了一部分"包装成"整个目录已绑定"。
    """
    path = _resolve_user_path(raw)
    empty = {"path": None, "exists": False, "files": 0,
             "listing_sha256": None, "content_sha256": None, "truncated": False}
    if path is None:
        return empty
    if path.is_file():
        try:
            size = path.stat().st_size
        except OSError:
            size = None
        return {"path": str(path), "exists": True, "files": 1,
                "listing_sha256": sha256_text(canonical_json([[path.name, size]])),
                "content_sha256": sha256_file(path), "truncated": False}
    if not path.is_dir():
        return {**empty, "path": str(path)}
    listing: list[list[Any]] = []
    content: list[list[Any]] = []
    truncated = False
    for child in sorted(path.rglob("*")):
        if not child.is_file():
            continue
        if len(listing) >= _EXECUTION_LIST_LIMIT:
            truncated = True
            break
        try:
            size = child.stat().st_size
        except OSError:
            size = None
        rel = child.relative_to(path).as_posix()
        listing.append([rel, size])
        if len(content) < _EXECUTION_TREE_LIMIT:
            content.append([rel, size, sha256_file(child)])
        else:
            truncated = True
    return {
        "path": str(path),
        "exists": True,
        "files": len(listing),
        "listing_sha256": sha256_text(canonical_json(listing)),
        "content_sha256": sha256_text(canonical_json(content)),
        "truncated": truncated,
    }


def _execution_config(config: Any) -> dict[str, Any]:
    """会改变"算出来的数"、却不在物理协议里的执行面：基组、模板、二进制。

    物理协议（:func:`_protocol_config`）管分子/电极/输运；species defaults、
    提交模板和 aims/aitranss 可执行文件同样决定结果。2026-11 复审的 P0-5 是：
    改这些配置不会触发 run_spec 重绑，快照与实际计算脱节。这里单独成一份
    指纹，而不是塞进 protocol——否则换集群路径/排队参数也会让历史 run 漂移。
    """
    cluster = getattr(config, "cluster", None)
    paths = getattr(config, "paths", None)
    env = getattr(cluster, "env", None)
    env = env if isinstance(env, Mapping) else {}
    if getattr(config.transport, "backend", "aitranss") == "transiesta":
        return {"module": "siesta/5.4.2-gpu-netcdf", "version": "5.4.2"}
    species_raw = (getattr(cluster, "species_dir", None)
                   or getattr(paths, "species_defaults", None))
    templates_raw = getattr(paths, "templates", None)
    executables: dict[str, Any] = {}
    for key in ("aims_bin", "aitranss_bin", "aitransss_bin", "g16root"):
        value = env.get(key)
        if value:
            executables[key] = _executable_probe(value)
    selected_env = {
        str(key): str(env[key])
        for key in _EXECUTION_ENV_KEYS
        if env.get(key) is not None
    }
    return {
        "species_defaults": _tree_probe(species_raw),
        "templates": _tree_probe(templates_raw),
        "executables": executables,
        "cluster_env": selected_env,
    }


def _protocol_config(config: Any) -> dict[str, Any]:
    """只取会改变物理结果的配置子集（路径/资源/报告设置不算）。"""
    molecule = _jsonable(config.molecule)
    xyz = None
    try:
        xyz = config.molecule_xyz()
    except Exception:  # noqa: BLE001
        pass
    molecule["xyz_sha256"] = sha256_file(xyz) if xyz is not None else None
    molecule.pop("xyz", None)
    # 0.2.x 时代的 protocol 里没有 initial_moments；未设时与旧行为逐字节等价，
    # 所以只在非空时才入协议，避免历史 run_spec 集体判为漂移。
    if not molecule.get("initial_moments"):
        molecule.pop("initial_moments", None)
    # fixed_spin_moment=True 是旧行为（模板一直写 M−1）；默认值不入协议，
    # 只有显式解除约束（false）才需要触发漂移检查。
    if molecule.get("fixed_spin_moment") is True:
        molecule.pop("fixed_spin_moment", None)
    # broken_symmetry=False 是旧行为（没有 BS 途径）；同 fixed_spin_moment，只在
    # 显式开启时入协议，避免历史 run_spec 集体判为漂移。
    if molecule.get("broken_symmetry") is False:
        molecule.pop("broken_symmetry", None)

    electrode = _jsonable(config.electrode)
    # 库文件位置不是物理参数；库文件内容哈希才是。
    electrode.pop("library_dir", None)
    electrode["library_file_sha256"] = _electrode_library_info(config).get("sha256")
    # 0.2.x 时代的 protocol 里没有 tip_direction 这个键；默认值 axis 与旧行为逐字节
    # 等价，所以只在非默认（lone_pair）时才入协议，避免历史 run_spec 集体判为漂移。
    if electrode.get("tip_direction") == "axis":
        electrode.pop("tip_direction", None)

    transport = _jsonable(config.transport)
    # 口径标记：eV 相对口径把"用户写的三元组"入协议（不随 SCF 的 efermi 变）；
    # 旧 Hartree 绝对口径保持 0.2.x 的哈希形状（不新增键），避免历史 run_spec
    # 仅仅因为新增了 unit 字段就集体判为漂移。
    unit = transport.pop("energy_window_unit", "eV")
    try:
        if unit == "eV":
            transport["energy_window"] = list(config.transport.relative_window())
            transport["energy_window_unit"] = "eV"
        else:
            transport["energy_window"] = list(config.transport.resolved_window())
    except Exception:  # noqa: BLE001 - 兼容手写/旧配置对象
        pass
    # 2026-11-20 复核 P2-1：ener/estep/eend 是"采样覆盖"——它们确实会改变
    # 积分网格（config.resolved_window 里 estep 覆盖生效），所以**显式写了的
    # 覆盖值必须进协议**。旧实现无条件 pop 掉三个键，于是 eV 相对窗下改
    # transport.estep 不改变 protocol_sha256（复核实测）。
    # 默认配置（三者皆 None）不加键，保持 0.2.x 的哈希形状逐字节不变。
    for key in ("ener", "estep", "eend"):
        transport.pop(key, None)
        value = getattr(config.transport, key, None)
        if value is not None:
            transport[key] = float(value)
    # 这是"缺 te.out 是否拒绝"的验证策略，不改变物理结果，不进 protocol。
    transport.pop("require_te_out", None)
    # 额外输出只决定"多存哪些证据文件"，不改 SCF/输运的任何求解参数；与
    # require_te_out 同理不进 protocol，否则给老 run 打开开关就会集体判为漂移。
    transport.pop("extra_outputs", None)
    if transport.get("backend") == "transiesta":
        settings = _jsonable(config.siesta)
        pseudos = settings.get("pseudopotentials", {})
        settings["pseudopotential_sha256"] = {element: sha256_file(config.resolve(path))
                                              for element, path in pseudos.items()}
        return {"molecule": molecule, "electrode": electrode, "transport": transport, "siesta": settings}
    transport.pop("backend", None)
    return {"molecule": molecule, "electrode": electrode, "transport": transport}


def protocol_sha256(config: Any) -> str:
    """物理协议指纹：分子/电极/输运配置的稳定哈希。"""
    return sha256_text(canonical_json(_protocol_config(config)))


def protocol_components(config: Any) -> dict[str, Any]:
    """``protocol_sha256`` 的输入字典（供 drift 明细/比对使用）。"""
    return _protocol_config(config)


def execution_components(config: Any) -> dict[str, Any]:
    """执行面指纹的输入字典（species/templates/二进制；见 P0-5）。"""
    return _execution_config(config)


def execution_sha256(config: Any) -> str:
    """执行面指纹：species defaults、提交模板与可执行文件的内容哈希。"""
    return sha256_text(canonical_json(_execution_config(config)))


#: 会直接改变物理结果的键：schema 演进时也必须两边都在且相等，否则仍判漂移。
PROTOCOL_CRITICAL = (
    "molecule.charge",
    "molecule.multiplicity",
    "molecule.fixed_spin_moment",
    "molecule.broken_symmetry",
    "transport.temperature",
    "transport.energy_window",
    "transport.energy_window_unit",
    "electrode.left_anchor",
    "electrode.right_anchor",
    "electrode.au_tip_distance",
)

#: 这些键"缺失"等价于某个默认值（``True``/``False`` 是旧行为，协议里不落键）。
PROTOCOL_CRITICAL_DEFAULTS: dict[str, Any] = {
    "molecule.fixed_spin_moment": True,
    "molecule.broken_symmetry": False,
}


def _flatten_protocol(data: Any, prefix: str = "") -> dict[str, str]:
    """把协议 dict 摊平成 ``路径 -> 规范化 JSON 标量``。"""
    flat: dict[str, str] = {}
    if isinstance(data, dict):
        for key, value in data.items():
            path = "%s.%s" % (prefix, key) if prefix else str(key)
            flat.update(_flatten_protocol(value, path))
    elif isinstance(data, (list, tuple)):
        flat[prefix] = canonical_json(list(data))
    else:
        flat[prefix] = canonical_json(data)
    return flat


def compare_protocol_configs(stored: Any, current: Any) -> dict[str, Any]:
    """比较 run 时与当前的协议子集，区分"形状变了"和"物理参数变了"。

    协议 schema 会随版本演进（新增/移除只影响记录方式的键），此时
    ``protocol_sha256`` 必然不同，但只要**两边都有的键**全部相同、且
    :data:`PROTOCOL_CRITICAL` 一个不缺，就不能判成物理漂移。
    """
    a = _flatten_protocol(stored or {})
    b = _flatten_protocol(current or {})
    shared = sorted(set(a) & set(b))
    differences = [
        {"path": path, "run": a[path], "now": b[path]}
        for path in shared if a[path] != b[path]
    ]
    only_run = sorted(set(a) - set(b))
    only_now = sorted(set(b) - set(a))
    critical_missing: list[str] = []
    for path in PROTOCOL_CRITICAL:
        if path in a and path in b:
            continue                      # 已落在 shared 的相等/差异比较里
        if path in PROTOCOL_CRITICAL_DEFAULTS:
            # 缺键等价于默认值（旧行为不落键），按默认值补上再比。
            default = canonical_json(PROTOCOL_CRITICAL_DEFAULTS[path])
            run_value = a.get(path, default)
            now_value = b.get(path, default)
            if run_value != now_value:
                differences.append({"path": path, "run": run_value,
                                    "now": now_value})
            continue
        critical_missing.append(path)
    return {
        "shared": len(shared),
        "differences": differences,
        "only_run": only_run,
        "only_now": only_now,
        "critical_missing": critical_missing,
    }


# --------------------------------------------------------------------------- #
# 软件/代码信息
# --------------------------------------------------------------------------- #
def _package_version() -> str:
    try:
        from importlib.metadata import version
        return version("smj-toolkit")
    except Exception:  # noqa: BLE001
        try:
            from . import __version__
            return __version__
        except Exception:  # noqa: BLE001
            return "unknown"


def _git_info() -> dict[str, Any]:
    root: Path | None = None
    for parent in Path(__file__).resolve().parents:
        if (parent / ".git").exists():
            root = parent
            break
    info: dict[str, Any] = {"root": str(root) if root else None,
                            "commit": None, "dirty": None}
    if root is None:
        return info

    def _run(*args: str) -> str | None:
        try:
            proc = subprocess.run(
                ["git", "-C", str(root), *args], capture_output=True,
                text=True, timeout=3, check=False)
        except (OSError, subprocess.SubprocessError):
            return None
        return proc.stdout.strip() if proc.returncode == 0 else None

    info["commit"] = _run("rev-parse", "HEAD")
    status = _run("status", "--porcelain")
    info["dirty"] = None if status is None else bool(status)
    return info


def _software_info(config: Any) -> dict[str, Any]:
    env = {str(key).upper(): str(value) for key, value in
           (getattr(getattr(config, "cluster", None), "env", {}) or {}).items()}
    executables: dict[str, Any] = {}
    for key in ("G16ROOT", "G16_BIN", "AIMS_BIN", "AITRANSS_BIN", "MKL_LIB"):
        raw = env.get(key) or os.environ.get(key)
        if not raw:
            continue
        raw = os.path.expandvars(os.path.expanduser(raw))
        remote_posix = (raw.startswith("/")
                        and getattr(getattr(config, "run", None), "mode", None) == "slurm-ssh")
        if remote_posix:
            # 集群上的 POSIX 路径在 Windows 上会被 Path() 改写成 \data\...；
            # 原样记录，exists 用 None 表示“本机不知道”，不谎报 False。
            entry: dict[str, Any] = {
                "path": raw,
                "exists": None,
                "note": "slurm-ssh 远端 POSIX 路径；本机未检查",
            }
        else:
            path = Path(raw)
            entry = {"path": str(path), "exists": path.exists()}
            if path.is_file():
                digest = sha256_file(path)
                if digest is None:
                    entry["note"] = "文件超过 %d MiB，未哈希" % (MAX_HASH_BYTES // (1024 * 1024))
                else:
                    entry["sha256"] = digest
        executables[key] = entry
    return {
        "smj": _package_version(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "git": _git_info(),
        "executables": executables,
        "cluster_env_keys": sorted(env),
    }


# --------------------------------------------------------------------------- #
# 远端软件版本（collect 时从输出文件解析；run_spec 的 append-only 观测）
# --------------------------------------------------------------------------- #
def parse_aims_versions(text: str) -> dict[str, Any]:
    """从 ``aims.out`` 头部解析 FHI-aims 版本/uuid/commit。

    只读头几十行也够用；解析不到时返回空 dict（调用方记 note，不伪造版本）。
    """
    head = "\n".join(text.splitlines()[:120])
    version = _AIMS_VERSION_RE.search(head)
    uuid = _AIMS_UUID_RE.search(head)
    commit = _AIMS_COMMIT_RE.search(head)
    if not (version or uuid):
        return {}
    result: dict[str, Any] = {
        "program": "FHI-aims",
        "source": "aims.out",
        "version": version.group(1) if version else None,
    }
    if uuid:
        result["aims_uuid"] = uuid.group(1)
    if commit:
        result["commit"] = commit.group(1)
    return result


def parse_aitranss_versions(text: str) -> dict[str, Any]:
    """从 ``te.out`` 头部解析 aitranss 版本与编译信息。"""
    head = "\n".join(text.splitlines()[:80])
    version = _AITRANSS_VERSION_RE.search(head)
    compiled = _AITRANSS_COMPILED_RE.search(head)
    if not (version or compiled):
        return {}
    return {
        "program": "aitranss",
        "source": "te.out",
        "version": version.group(1).strip() if version else None,
        "compiled": compiled.group(1).strip() if compiled else None,
    }


def record_observed_software(run_dir: str | Path, stage: str,
                             software: Mapping[str, Any] | None) -> dict[str, Any] | None:
    """把 collect 时观测到的远端软件版本追加进 run_spec。

    设计边界：``config`` / ``protocol_sha256`` 仍是 run 级不可变快照；这里只追加
    ``software.observed[stage]``，不改变物理协议，也不触发 config 漂移门。legacy
    run 没有 ``run_spec.json`` 时写 ``<run>/software_observed.json``，不假装重建 spec。
    """
    if not software:
        return None
    payload = {key: value for key, value in software.items() if value is not None}
    if not payload:
        return None
    record = {"observed_at": now_iso(), **payload}
    root = Path(run_dir)
    spec = load_run_spec(root)
    if spec is not None:
        observed = spec.setdefault("software", {}).setdefault("observed", {})
        observed[stage] = record
        _write_json(spec, root / RUN_SPEC_FILE)
        return record

    path = root / SOFTWARE_OBSERVED_FILE
    data: dict[str, Any] = {}
    if path.is_file():
        try:
            loaded = json.loads(path.read_text(encoding="utf-8"))
            if isinstance(loaded, dict):
                data = loaded
        except (OSError, ValueError):
            data = {}
    data[stage] = record
    _write_json(data, path)
    return record


# --------------------------------------------------------------------------- #
# run_spec
# --------------------------------------------------------------------------- #
def run_spec_payload(config: Any, run_dir: str | Path, *,
                     reconstructed: bool = False,
                     rebound: bool = False,
                     reason: str | None = None) -> dict[str, Any]:
    """构造 run_spec dict（不写文件）。"""
    root = Path(run_dir).resolve()
    config_data = config_to_dict(config)
    spec: dict[str, Any] = {
        "schema": RUN_SPEC_SCHEMA,
        "created_at": now_iso(),
        "reconstructed": bool(reconstructed),
        "rebound": bool(rebound),
        "run_id": str(getattr(config, "name", "") or root.name),
        "run_dir": str(root),
        "config_sha256": sha256_text(canonical_json(config_data)),
        "protocol_sha256": protocol_sha256(config),
        "protocol_schema": PROTOCOL_SCHEMA,
        "execution_sha256": execution_sha256(config),
        "execution": execution_components(config),
        "config": config_data,
        "electrode_library": _electrode_library_info(config),
        "software": _software_info(config),
    }
    if getattr(config.transport, "backend", "aitranss") == "transiesta":
        spec.pop("electrode_library", None)
        spec["software"] = {"backend": "transiesta", "target_version": "5.4.2"}
    if reason:
        spec["reason"] = reason
    if reconstructed:
        spec["reconstructed_at"] = spec["created_at"]
        spec["reconstruction_reason"] = reason or (
            "legacy run 没有 run_spec.json；用当前可加载的配置重建，"
            "不能证明与原始 run 完全一致")

    xyz = None
    try:
        xyz = config.molecule_xyz()
    except Exception:  # noqa: BLE001
        pass
    spec["molecule"] = {
        "path": str(Path(xyz)) if xyz is not None else None,
        "exists": bool(xyz is not None and Path(xyz).is_file()),
        "sha256": sha256_file(xyz) if xyz is not None else None,
    }
    return spec


def _write_json(payload: Mapping[str, Any], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8", newline="\n")
    os.replace(tmp, path)
    return path


def load_run_spec(run_dir: str | Path) -> dict[str, Any] | None:
    """读 ``<run>/run_spec.json``；不存在/不可解析返回 ``None``。"""
    path = Path(run_dir) / RUN_SPEC_FILE
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return payload if isinstance(payload, dict) else None


def run_spec_sha256(run_dir: str | Path) -> str | None:
    """``run_spec.json`` 文件本身的 sha256（manifest 绑定用）。"""
    return sha256_file(Path(run_dir) / RUN_SPEC_FILE)


def ensure_run_spec(config: Any, run_dir: str | Path, *,
                    force: bool = False) -> dict[str, Any]:
    """确保 run_spec 存在且与当前配置的**物理协议**一致。

    - 首次 prepare：写 ``run_spec.json``；
    - 已有且协议一致：原样返回（config/protocol 是 run 级不可变快照；
      ``software.observed`` 由 :func:`record_observed_software` 在 collect 后追加）；
    - 已有但协议漂移：默认 :class:`ConfigError`，``force=True`` 时把旧文件
      复制成 ``run_spec.previous.json`` 后重绑（新 spec 带 ``rebound: true``）。
    """
    root = Path(run_dir)
    path = root / RUN_SPEC_FILE
    current = protocol_sha256(config)
    current_exec = execution_sha256(config)

    def _rebind(reason: str) -> dict[str, Any]:
        previous = root / RUN_SPEC_PREVIOUS
        try:
            shutil.copy2(path, previous)
        except OSError as exc:
            raise ConfigError("无法保留旧 run_spec（%s）: %s" % (previous, exc)) from exc
        spec = run_spec_payload(config, root, rebound=True, reason=reason)
        _write_json(spec, path)
        return spec

    if not path.is_file():
        spec = run_spec_payload(config, root)
        _write_json(spec, path)
        return spec

    try:
        existing = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ConfigError(
            "run_spec.json 不可解析（%s）；请检查文件或删除后重跑 --prepare，"
            "不要让它静默失效" % (path,)) from exc
    if not isinstance(existing, dict):
        raise ConfigError("run_spec.json 顶层必须是对象: %s" % path)

    stored = existing.get("protocol_sha256")
    stored_exec = existing.get("execution_sha256")
    if stored != current:
        if not force:
            raise ConfigError(
                "run_spec 与当前配置不一致（protocol_sha256: run=%s, 当前=%s）。"
                "run_spec 是 run 级不可变快照：请恢复原配置，或确认要重绑后加 --force"
                % (str(stored)[:12] or "缺失", current[:12]))
        return _rebind("物理协议漂移后 --force 重绑；旧快照见 %s" % RUN_SPEC_PREVIOUS)

    if stored_exec is None:
        # 旧 run_spec 只有物理协议、没有执行面。这里只补齐字段，不改写
        # config/protocol——不伪装成重绑；同时明确标注"基线是本次才建立的"。
        patched = dict(existing)
        patched["execution"] = execution_components(config)
        patched["execution_sha256"] = current_exec
        patched["execution_backfilled_at"] = now_iso()
        patched["execution_backfill_note"] = (
            "旧 run_spec 未记录执行面（species defaults/模板/二进制）；"
            "本次 prepare 补齐基线，无法追溯此前是否已改过这些输入")
        _write_json(patched, path)
        return patched

    if stored_exec == current_exec:
        return existing
    if not force:
        raise ConfigError(
            "run_spec 与当前执行面不一致（execution_sha256: run=%s, 当前=%s）："
            "species defaults/templates/可执行文件会改变算出来的数，"
            "但不会改变 protocol_sha256。请恢复原环境，或确认要重绑后加 --force"
            % (str(stored_exec)[:12], current_exec[:12]))
    return _rebind("执行面（species defaults/模板/二进制）漂移后 --force 重绑；"
                   "旧快照见 %s" % RUN_SPEC_PREVIOUS)


# --------------------------------------------------------------------------- #
# 阶段间指纹校验
# --------------------------------------------------------------------------- #
def verify_parent_geometry(stage: str, geometry: str | Path,
                           metadata: str | Path | Mapping[str, Any] | None) -> bool:
    """校验 ``geometry`` 的 sha256 与电极 metadata 的 ``surface_fingerprint``。

    这是"SCF/输运用的几何 == A2 确认/上游产出的几何"的硬断言；不一致时
    :class:`StageBlocked`，绝不继续。
    """
    geo = Path(geometry)
    if not geo.is_file():
        raise StageBlocked(stage, "父几何文件不存在: %s" % geo, [str(geo)])
    expected = None
    source: str | None = None
    if isinstance(metadata, Mapping):
        expected = metadata.get("surface_fingerprint")
        source = str(metadata.get("_path") or "metadata")
    elif metadata is not None:
        meta_path = Path(metadata)
        source = str(meta_path)
        if not meta_path.is_file():
            raise StageBlocked(stage, "电极 metadata 不存在: %s" % meta_path,
                              [str(meta_path)])
        try:
            payload = json.loads(meta_path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise StageBlocked(stage, "电极 metadata 不可解析: %s" % meta_path,
                              [str(meta_path)]) from exc
        expected = payload.get("surface_fingerprint") if isinstance(payload, dict) else None
        if isinstance(payload, dict) and payload.get("surface_fingerprint_structural"):
            # 结构指纹优先：initial_charge/initial_moment 会改变字节 sha256，
            # 但不应改变「这是同一份几何」这一事实。
            expected = payload["surface_fingerprint_structural"]
            actual_struct = structural_fingerprint(geo)
            if actual_struct is None:
                raise StageBlocked(stage, "几何文件没有可识别的原子行: %s" % geo,
                                  [str(geo)])
            if actual_struct != expected:
                raise StageBlocked(
                    stage, "几何结构指纹与电极 metadata 不一致（几何=%s，记录=%s，来源=%s）"
                    % (str(actual_struct)[:12], str(expected)[:12], source or "未提供"),
                    [str(geo)])
            return True

    actual = sha256_file(geo)
    if not expected:
        raise StageBlocked(
            stage, "电极 metadata 缺 surface_fingerprint，无法校验几何指纹（%s）"
            % (source or "未提供"), [str(geo)])
    if actual != expected:
        raise StageBlocked(
            stage, "几何指纹与电极 metadata 不一致（几何 sha256=%s，记录=%s，来源=%s）"
            % (str(actual)[:12], str(expected)[:12], source or "未提供"),
            [str(geo)])
    return True


__all__ = [
    "RUN_SPEC_SCHEMA",
    "RUN_SPEC_FILE",
    "RUN_SPEC_PREVIOUS",
    "PROTOCOL_SCHEMA",
    "MAX_HASH_BYTES",
    "sha256_file",
    "sha256_text",
    "canonical_json",
    "hash_stage_inputs",
    "hash_stage_outputs",
    "config_to_dict",
    "protocol_sha256",
    "protocol_components",
    "PROTOCOL_CRITICAL",
    "PROTOCOL_CRITICAL_DEFAULTS",
    "compare_protocol_configs",
    "run_spec_payload",
    "load_run_spec",
    "run_spec_sha256",
    "ensure_run_spec",
    "verify_parent_geometry",
    "SOFTWARE_OBSERVED_FILE",
    "parse_aims_versions",
    "parse_aitranss_versions",
    "record_observed_software",
]
