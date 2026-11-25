"""运行环境清单：决定数值的那几样东西，单独落盘。

为什么需要它：``thermal.json`` 里的 ``G_ph`` 只有在"同一套数值环境"下才谈得上
可复现。真正决定这个数的是四样东西 —— Python 与 numpy 的版本、FHI-aims 的
二进制、species 目录里每一份基组文件、以及随包发布的电极力常数。任何一处变了，
13.7 pW/K 就可能变成 13.2 或 14.1，而产物里**看不出差别**。

所以这里把它们收成一份 ``environment.json``，算一个 digest 写进 ``thermal.json``；
``smj thermal-compare`` 用同一套键比对各 run 的环境，不一致时显式警告。这与
D15 的协议比对是同一个思路：先让"不同"变得看得见，再谈谁对。

为什么不是 ``pip freeze`` 全量：这份清单要能回答"为什么这个数变了"，不是
"我装了哪些包"。所以只钉**参与计算的输入与工具**，其余一概不记。

``generated_at`` 只用于人看，**不进 digest** —— 否则同一环境跑两次也会不等。
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

SCHEMA = "smj.environment.v1"

#: 比对环境时要看的字段（``identity`` 的输出）。顺序即报告里的行序。
IDENTITY_KEYS = (
    "smj_version",
    "python_version",
    "numpy_version",
    "aims_sha256",
    "species_digest",
    "lead_fc_sha256",
)


def sha256_file(path: str | Path, *, chunk: int = 1 << 20) -> str:
    """文件的 sha256 十六进制串。"""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def _package_version(name: str) -> str | None:
    try:
        from importlib.metadata import PackageNotFoundError, version
    except ImportError:  # pragma: no cover - Python < 3.8 不在支持范围
        return None
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def _smj_version() -> str:
    version = _package_version("smj-toolkit")
    if version:
        return version
    try:
        import smj

        return str(getattr(smj, "__version__", "unknown"))
    except Exception:  # noqa: BLE001 - 版本号取不到不该让计算失败
        return "unknown"


def _package_root() -> Path | None:
    try:
        import smj

        return Path(smj.__file__).resolve().parent
    except Exception:  # noqa: BLE001
        return None


def _provenance_sha(root: Path | None) -> str | None:
    if root is None:
        return None
    path = root / "junction" / "PROVENANCE.json"
    if not path.is_file():
        return None
    return sha256_file(path)


def _thermal_snapshot(config) -> dict:
    thermal = getattr(config, "thermal", None)
    if thermal is None:
        return {}
    try:
        return dataclasses.asdict(thermal)
    except TypeError:  # pragma: no cover - 非 dataclass 时退回空表
        return {}


def _species_snapshot(config) -> dict:
    """species 目录里每份 ``*.in`` 基组文件的哈希。

    只收 ``.in``：FHI-aims 的 species 目录里也可能放 README 之类的东西，
    那些不影响数值，进清单只会让 digest 无谓地变。
    """
    cluster = getattr(config, "cluster", None)
    raw = getattr(cluster, "species_dir", None) if cluster is not None else None
    if raw is None:
        return {"directory": None, "files": {}, "digest": None}
    directory = Path(raw)
    if not directory.is_dir():
        return {"directory": str(directory), "files": {}, "digest": None,
                "present": False}
    files = {}
    for path in sorted(directory.glob("*.in")):
        files[path.name] = sha256_file(path)
    combined = hashlib.sha256()
    for name in sorted(files):
        combined.update(name.encode("utf-8"))
        combined.update(files[name].encode("ascii"))
    return {"directory": str(directory), "files": files,
            "digest": combined.hexdigest() if files else None, "present": True}


def _aims_snapshot(config) -> dict:
    cluster = getattr(config, "cluster", None)
    env = getattr(cluster, "env", None) or {}
    configured = env.get("aims_bin") if isinstance(env, dict) else None
    if not configured or configured == "UNCONFIGURED":
        return {"configured": None, "exists": False, "sha256": None,
                "module": (env or {}).get("aims_ld_module") if isinstance(env, dict) else None}
    path = Path(str(configured))
    # 远程 run（slurm-ssh）时本地通常没有这个二进制：记路径，哈希留空，
    # 由 submit 脚本所在的集群负责版本一致性，不在这里假装知道。
    exists = path.is_file()
    return {
        "configured": str(path),
        "exists": exists,
        "sha256": sha256_file(path) if exists else None,
        "module": (env or {}).get("aims_ld_module") if isinstance(env, dict) else None,
    }


def _lead_fc_snapshot(config) -> dict:
    thermal = getattr(config, "thermal", None)
    value = getattr(thermal, "lead_fc", None) if thermal is not None else None
    if value in (None, "auto"):
        try:
            from .assets import lead_fc_info

            info = lead_fc_info()
        except Exception:  # noqa: BLE001 - 资产缺失由 thermal_analyze.prepare 报
            return {"source": "auto", "sha256": None, "in_plane": None}
        return {"source": "auto", "sha256": info["sha256"],
                "in_plane": info.get("in_plane"),
                "layers": info.get("layers"),
                "a_lattice_A": info.get("a_lattice_A")}
    path = Path(str(value))
    if not path.is_file():
        return {"source": str(path), "sha256": None, "in_plane": None}
    payload = json.loads(path.read_text(encoding="utf-8"))
    return {"source": str(path), "sha256": payload.get("sha256"),
            "in_plane": payload.get("in_plane"),
            "layers": payload.get("layers"),
            "a_lattice_A": payload.get("a_lattice_A")}


def collect(config=None, *, extra: dict | None = None) -> dict:
    """收一份环境清单。

    ``config`` 是可选的：纯 API 调用（不给 RunConfig）时只有解释器与依赖
    版本，配置相关字段为空表 —— 这样单测与脚本也能用同一份产出。
    """
    root = _package_root()
    payload = {
        "schema": SCHEMA,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "smj": {"version": _smj_version(),
                "root": str(root) if root else None},
        "python": {
            "version": platform.python_version(),
            "implementation": platform.python_implementation(),
            "executable": sys.executable,
        },
        "platform": {
            "system": platform.system(),
            "machine": platform.machine(),
            "release": platform.release(),
            "processor": platform.processor() or None,
        },
        "packages": {
            "numpy": _package_version("numpy"),
            "scipy": _package_version("scipy"),
            "matplotlib": _package_version("matplotlib"),
            "PyYAML": _package_version("PyYAML"),
        },
        # 电子链内核的 sha256 清单：热导链复用 junction.geometry 等模块，
        # 那边一动这个数就变，提醒"基线是否还可比"。
        "junction_provenance_sha256": _provenance_sha(root),
        "thermal": _thermal_snapshot(config),
        "binaries": {"aims": _aims_snapshot(config) if config is not None else {}},
        "species": _species_snapshot(config) if config is not None else {},
        "lead_fc": _lead_fc_snapshot(config) if config is not None else {},
    }
    if extra:
        payload["extra"] = dict(extra)
    return payload


def canonical(payload: dict) -> str:
    """稳定的 JSON 文本：键排序、无多余空白、``generated_at`` 除外。"""
    stripped = {key: value for key, value in payload.items() if key != "generated_at"}
    return json.dumps(stripped, sort_keys=True, ensure_ascii=False,
                      separators=(",", ":"))


def digest(payload: dict) -> str:
    """环境清单的指纹（不含 ``generated_at``）。"""
    return hashlib.sha256(canonical(payload).encode("utf-8")).hexdigest()


def identity(payload: dict | None) -> dict:
    """把清单压成可跨 run 比对的几个键。

    只留"换了它就可能换取数"的东西。路径不进比对：两台机器上同一个 FHI-aims
    装在两个目录里，仍然是同一个二进制。
    """
    payload = payload or {}
    packages = payload.get("packages") or {}
    python = payload.get("python") or {}
    aims = (payload.get("binaries") or {}).get("aims") or {}
    species = payload.get("species") or {}
    lead = payload.get("lead_fc") or {}
    return {
        "smj_version": (payload.get("smj") or {}).get("version"),
        "python_version": python.get("version"),
        "numpy_version": packages.get("numpy"),
        "aims_sha256": aims.get("sha256"),
        "species_digest": species.get("digest"),
        "lead_fc_sha256": lead.get("sha256"),
    }


def describe_lines(payload: dict) -> list[str]:
    """给人看的三五行事。"""
    ident = identity(payload)
    python = payload.get("python") or {}
    species = payload.get("species") or {}
    aims = (payload.get("binaries") or {}).get("aims") or {}

    def short(value, width: int = 12) -> str:
        return str(value)[:width] if value else "—"

    return [
        "smj %s；Python %s（%s）" % (ident["smj_version"] or "unknown",
                                    ident["python_version"] or "unknown",
                                    python.get("implementation") or "?"),
        "numpy %s" % (ident["numpy_version"] or "—"),
        "FHI-aims：%s（sha256 %s）"
        % (aims.get("configured") or "未配置", short(ident["aims_sha256"])),
        "species：%s（%d 份，digest %s）"
        % (species.get("directory") or "未配置",
           len(species.get("files") or {}), short(ident["species_digest"])),
        "电极力常数 sha256 %s" % short(ident["lead_fc_sha256"]),
    ]


__all__ = ["SCHEMA", "IDENTITY_KEYS", "collect", "identity", "digest",
           "canonical", "describe_lines", "sha256_file"]
