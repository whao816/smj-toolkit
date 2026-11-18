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


def hash_stage_outputs(stage_dir: str | Path, artifacts: Iterable[str]) -> dict[str, str]:
    """把 collect 声明的产物哈希成 ``{相对名: sha256}``。"""
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
    # 输运窗同时保留"用户写的三元组"和"实际生效的窗口"（ener/estep/eend 覆盖）。
    try:
        data["transport"]["energy_window"] = list(config.transport.energy_window)
        data["transport"]["resolved_window"] = list(config.transport.resolved_window())
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
    transport.pop("ener", None)
    transport.pop("estep", None)
    transport.pop("eend", None)
    # 这是"缺 te.out 是否拒绝"的验证策略，不改变物理结果，不进 protocol。
    transport.pop("require_te_out", None)
    return {"molecule": molecule, "electrode": electrode, "transport": transport}


def protocol_sha256(config: Any) -> str:
    """物理协议指纹：分子/电极/输运配置的稳定哈希。"""
    return sha256_text(canonical_json(_protocol_config(config)))


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
        "config": config_data,
        "electrode_library": _electrode_library_info(config),
        "software": _software_info(config),
    }
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
    if stored == current:
        return existing
    if not force:
        raise ConfigError(
            "run_spec 与当前配置不一致（protocol_sha256: run=%s, 当前=%s）。"
            "run_spec 是 run 级不可变快照：请恢复原配置，或确认要重绑后加 --force"
            % (str(stored)[:12] or "缺失", current[:12]))
    previous = root / RUN_SPEC_PREVIOUS
    try:
        shutil.copy2(path, previous)
    except OSError as exc:
        raise ConfigError("无法保留旧 run_spec（%s）: %s" % (previous, exc)) from exc
    spec = run_spec_payload(
        config, root, rebound=True,
        reason="配置漂移后 --force 重绑；旧快照见 %s" % RUN_SPEC_PREVIOUS)
    _write_json(spec, path)
    return spec


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
    "MAX_HASH_BYTES",
    "sha256_file",
    "sha256_text",
    "canonical_json",
    "hash_stage_inputs",
    "hash_stage_outputs",
    "config_to_dict",
    "protocol_sha256",
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
