"""提交前预检（preflight）：把 vendored 检查包成**结构化 issues**。

检查项（顺序即报告顺序）：

1. 分子结构可读 / 原子重叠；
2. ``species_dir`` 存在且覆盖分子元素（缺目录或缺元素 = error）；
3. 建结几何闸门（有 ``auN_build/geometry.in`` 时调 vendored
   ``validate_junction_geometry``：表面索引、输运轴、接触重叠、Au-S 锚定、A2 确认）；
4. 资源估算（每阶段墙钟/内存，超过 72 h / 128 GiB 上限 = error）；
5. 集群与 SSH 连通（``slurm-ssh`` 且**没有传输层**时明确"未执行远端检查"，
   不静默当作通过；给了 transport 才做一次真实往返）；
6. 远端目录（``remote_root`` 存在且可写）；
7. 二进制/许可（dry-run 明确跳过；local/slurm-ssh 才核对）。

返回 :class:`PreflightResult`：``.ok`` 属性、``.issues``（每项
``{severity: error|warning|info, code, message, ...}``）、``.errors/.warnings``。
它同时支持 ``result["ok"]`` / ``result.get(...)`` 这类 dict 访问，便于 CLI 与
测试两种风格。

vendored 模块全部 lazy import；本模块 import 不触网、不做真实计算。
"""

from __future__ import annotations

import os
import shlex
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .resources import MEM_CAP_GIB, WALLTIME_CAP_S, check_smj_budget, estimate_smj_resources
from .state import RunLayout

SEVERITIES = ("error", "warning", "info")

#: 资源硬上限（与 smj.resources 同源）。
CAPS = {"walltime_s": WALLTIME_CAP_S, "mem_gib": MEM_CAP_GIB}

#: 配置里小写 env 键 -> vendored preflight 认识的大写环境变量。
_ENV_ALIASES = {
    "aims_bin": "AIMS_BIN",
    "aims_module": "AIMS_MODULE",
    "aims_ld_module": "AIMS_LD_MODULE",
    "aitranss_bin": "AITRANSS_BIN",
    "g16root": "G16ROOT",
    "g16_module": "G16_MODULE",
    "gauss_lic": "GAUSS_LIC",
    "g16_scratch": "G16_SCRATCH",
}

#: 资源检查的阶段（本地阶段 au60_build 不计集群预算）。
_RESOURCE_STAGES = ("gaussian_opt", "aims_tip_relax", "au_relax", "aims_scf", "aitranss_te")

#: 分子重叠检查的规模上限（超过只提示，不做 O(N^2)）。
_OVERLAP_MAX_ATOMS = 2000


# --------------------------------------------------------------------------- #
# 结果对象
# --------------------------------------------------------------------------- #
@dataclass
class PreflightResult(Mapping):
    """预检结果。``ok`` 仅由 error 级 issue 决定（warning/info 不阻断）。"""

    issues: list[dict] = field(default_factory=list)
    remote_checked: bool = False
    notes: list[str] = field(default_factory=list)
    estimates: dict[str, dict] = field(default_factory=dict)
    elements: list[str] = field(default_factory=list)
    species_dir: str | None = None
    mode: str | None = None
    run: str | None = None

    def add(self, severity: str, code: str, message: str, **extra: Any) -> dict:
        if severity not in SEVERITIES:
            raise ValueError("非法 severity %r（可选 %s）" % (severity, "/".join(SEVERITIES)))
        issue = {"severity": severity, "code": code, "message": message}
        issue.update({k: v for k, v in extra.items() if v is not None})
        self.issues.append(issue)
        return issue

    @property
    def ok(self) -> bool:
        return not any(i["severity"] == "error" for i in self.issues)

    @property
    def errors(self) -> int:
        return sum(1 for i in self.issues if i["severity"] == "error")

    @property
    def warnings(self) -> int:
        return sum(1 for i in self.issues if i["severity"] == "warning")

    @property
    def infos(self) -> int:
        return sum(1 for i in self.issues if i["severity"] == "info")

    def codes(self) -> list[str]:
        return [i.get("code", "") for i in self.issues]

    def to_dict(self) -> dict:
        return {
            "ok": self.ok,
            "errors": self.errors,
            "warnings": self.warnings,
            "issues": list(self.issues),
            "remote_checked": self.remote_checked,
            "notes": list(self.notes),
            "estimates": dict(self.estimates),
            "elements": list(self.elements),
            "species_dir": self.species_dir,
            "mode": self.mode,
            "run": self.run,
        }

    # -- dict 式访问（CLI 两种风格都能用） ---------------------------------- #
    def __getitem__(self, key: str) -> Any:
        return self.to_dict()[key]

    def get(self, key: str, default: Any = None) -> Any:
        return self.to_dict().get(key, default)

    def __contains__(self, key: object) -> bool:
        return key in self.to_dict()

    def keys(self):
        return self.to_dict().keys()

    def items(self):
        return self.to_dict().items()

    def __iter__(self):
        return iter(self.to_dict())

    def __len__(self) -> int:
        return len(self.to_dict())


# --------------------------------------------------------------------------- #
# 环境 / species / 结构
# --------------------------------------------------------------------------- #
def build_env(config: Any, env: Mapping[str, Any] | None = None) -> dict[str, str]:
    """os.environ + ``config.cluster.env`` + 显式 ``env``，并补大写别名。"""
    merged: dict[str, str] = dict(os.environ)
    layers = [getattr(getattr(config, "cluster", None), "env", None), env]
    for layer in layers:
        for key, value in dict(layer or {}).items():
            if value is None:
                continue
            merged[str(key)] = str(value)
            alias = _ENV_ALIASES.get(str(key).lower())
            if alias:
                merged[alias] = str(value)
    return merged


def resolve_species_dir(config: Any, species_dir: str | Path | None = None) -> Path | None:
    """``显式参数 > cluster.species_dir > paths.species_defaults``，解析为绝对路径。"""
    raw = species_dir
    if raw is None:
        raw = (getattr(getattr(config, "cluster", None), "species_dir", None)
               or getattr(getattr(config, "paths", None), "species_defaults", None))
    if not raw:
        return None
    resolver = getattr(config, "resolve", None)
    if callable(resolver):
        resolved = resolver(raw)
        return Path(resolved) if resolved is not None else Path(raw)
    path = Path(raw)
    return path if path.is_absolute() else path.resolve()


def molecule_atoms(path: str | Path) -> list[Any]:
    """读分子结构（xyz 或 FHI-aims geometry.in），返回 vendored ``Atom`` 列表。"""
    p = Path(path)
    from smj._vendor.junction.geometry import read_aims_geometry, read_xyz  # lazy: vendored
    if p.suffix.lower() == ".xyz" or p.name.lower().endswith(".xyz"):
        return read_xyz(p)
    try:
        return read_aims_geometry(p)
    except ValueError:
        return read_xyz(p)


def _molecule_check(result: PreflightResult, config: Any,
                    elements: list[str] | None) -> list[Any] | None:
    """分子结构存在/可读 + 原子重叠；返回 atoms（供元素/规模复用）。"""
    path = None
    try:
        path = config.molecule_xyz()
    except Exception as exc:  # noqa: BLE001
        result.add("error", "molecule_unresolved", "无法解析分子结构路径: %s" % exc)
        return None
    if path is None:
        result.add("error", "molecule_missing",
                   "未指定分子结构（molecule.xyz / --xyz），gaussian_opt 无法 prepare")
        return None
    path = Path(path)
    if not path.is_file():
        result.add("error", "molecule_not_found", "分子结构文件不存在: %s" % path)
        return None
    try:
        atoms = molecule_atoms(path)
    except Exception as exc:  # noqa: BLE001
        result.add("error", "molecule_unreadable", "分子结构不可解析 (%s): %s" % (path, exc))
        return None
    labels = list(elements or [])
    if not labels:
        try:
            from smj._vendor.junction.geometry import elements_in_geometry  # lazy: vendored
            labels = elements_in_geometry(atoms)
        except Exception:  # noqa: BLE001
            labels = []
    result.add("info", "molecule_ok",
               "分子结构 %s：%d 原子，元素 %s"
               % (path.name, len(atoms), ", ".join(labels) or "未解析"))
    if len(atoms) < _OVERLAP_MAX_ATOMS:
        _overlap_check(result, atoms)
    else:
        result.notes.append("分子 %d 原子，跳过 O(N^2) 重叠预检（由结构门在 prepare 阶段兜底）"
                            % len(atoms))
    return atoms


def _overlap_check(result: PreflightResult, atoms: list[Any]) -> None:
    """元素对 + 共价半径下限判原子重叠（与 vendored validate 同一口径）。"""
    try:
        from smj._vendor.junction.validate import (  # lazy: vendored
            COVALENT_RADII, MIN_BOND_FACTOR, PAIR_CUTOFF)
    except Exception:  # noqa: BLE001
        result.notes.append("vendored validate 不可用，跳过分子重叠预检")
        return
    import math

    def floor(e1: str, e2: str):
        r1, r2 = COVALENT_RADII.get(e1), COVALENT_RADII.get(e2)
        return None if r1 is None or r2 is None else MIN_BOND_FACTOR * (r1 + r2)

    overlaps = []
    coords = [a.coord() for a in atoms]
    elems = [a.element for a in atoms]
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            limit = floor(elems[i], elems[j])
            if limit is None:
                continue
            distance = math.dist(coords[i], coords[j])
            if distance < PAIR_CUTOFF and distance < limit:
                overlaps.append({"i": i + 1, "j": j + 1,
                                 "elements": "%s-%s" % (elems[i], elems[j]),
                                 "d": round(distance, 3), "floor": round(limit, 3)})
    if overlaps:
        result.add("error", "molecule_overlap",
                   "%d 对原子短于元素对键长下限（0.7×共价半径和）: %s"
                   % (len(overlaps), overlaps[:3]), overlaps=overlaps[:8])
    else:
        result.add("info", "molecule_no_overlap",
                   "%d 原子无重叠（元素对下限 + %.1f A 预筛）"
                   % (len(atoms), PAIR_CUTOFF))


def _species_check(result: PreflightResult, config: Any,
                   species_dir: str | Path | None, elements: list[str]) -> None:
    path = resolve_species_dir(config, species_dir)
    result.species_dir = str(path) if path else None
    if path is None:
        result.add("error", "species_dir_missing",
                   "未配置 species_dir（cluster.species_dir / paths.species_defaults）——"
                   " FHI-aims 会在远端才因缺 species 失败，浪费一次提交+排队")
        return
    if not path.is_dir():
        result.add("error", "species_dir_not_found", "species_dir 不存在: %s" % path)
        return
    if not elements:
        result.add("warning", "species_elements_unknown",
                   "species_dir 已给（%s），但分子元素未知，未逐元素核对" % path)
        return
    from smj._vendor.junction.preflight import species_available  # lazy: vendored
    report = species_available(path, elements)
    if report.get("ok"):
        result.add("info", "species_ok", report.get("detail", "species 齐备"),
                   species=report.get("found"))
    else:
        result.add("error", "species_missing", report.get("detail", "species 缺失"),
                   missing=report.get("missing"), species_dir=str(path))


def _find_build_dir(run_dir: Path) -> Path | None:
    try:
        for child in sorted(run_dir.iterdir()):
            if child.is_dir() and child.name.startswith("au") and child.name.endswith("_build"):
                return child
    except OSError:
        return None
    return None


def _geometry_check(result: PreflightResult, run_dir: Path, layout: RunLayout) -> None:
    build = layout.stage_dir("au60_build")
    if not build.is_dir():
        build = _find_build_dir(run_dir)
    geometry = build / "geometry.in" if build is not None else None
    if geometry is None or not geometry.is_file():
        result.add("info", "geometry_not_ready",
                   "尚无建结 geometry.in（将由 auN_build 阶段生成），结构几何闸门未执行")
        return
    try:
        from smj._vendor.junction.validate import validate_junction_geometry  # lazy: vendored
        verdict = validate_junction_geometry(build, geometry=geometry)
    except Exception as exc:  # noqa: BLE001
        result.add("error", "geometry_check_failed", "建结几何检查失败: %s" % exc)
        return
    for check in verdict.get("checks", []):
        if check.get("passed"):
            continue
        severity = {"error": "error", "warn": "warning"}.get(check.get("severity"), "info")
        result.add(severity, "geometry:%s" % check.get("name", "unknown"),
                   str(check.get("detail", "")), stage=build.name,
                   value=check.get("value"))
    if verdict.get("passed") and not verdict.get("warnings"):
        result.add("info", "geometry_ok",
                   "建结几何闸门通过（%s）" % verdict.get("checks", [{}])[0].get("detail", ""))
    elif verdict.get("passed"):
        result.notes.append("建结几何闸门通过，但有 %s 条方向性告警"
                            % verdict.get("warnings"))


# --------------------------------------------------------------------------- #
# 资源 / 集群 / 环境
# --------------------------------------------------------------------------- #
def _molecule_count(run_dir: Path, config: Any) -> int | None:
    path = None
    try:
        path = config.molecule_xyz()
    except Exception:  # noqa: BLE001
        path = None
    if path and Path(path).is_file():
        try:
            return len(molecule_atoms(path))
        except Exception:  # noqa: BLE001
            return None
    return None


def _geometry_count(path: Path) -> int | None:
    if not path.is_file():
        return None
    try:
        from smj._vendor.junction.geometry import read_aims_geometry  # lazy: vendored
        return len(read_aims_geometry(path))
    except Exception:  # noqa: BLE001
        return None


def _resource_check(result: PreflightResult, config: Any, run_dir: Path,
                    layout: RunLayout, stage: str | None,
                    natoms: int | None) -> None:
    molecule_n = natoms or _molecule_count(run_dir, config)
    build = _find_build_dir(run_dir)
    build_n = _geometry_count(build / "geometry.in") if build is not None else None
    tip_path = layout.artifact("aims_tip_relax", "geometry.in.next_step")
    tip_n = _geometry_count(tip_path)
    counts = {
        "gaussian_opt": molecule_n,
        "aims_tip_relax": tip_n or ((molecule_n + 2) if molecule_n else None),
        "aims_scf": build_n or natoms,
        "aitranss_te": build_n or natoms,
    }
    stages = [stage] if stage else list(_RESOURCE_STAGES)
    for name in stages:
        if name == "au60_build":
            n = build_n or molecule_n
            if not n:
                result.add("warning", "resource_skipped",
                           "au60_build: 无原子数，资源估算跳过")
                continue
            estimate = estimate_smj_resources(name, int(n))
            result.estimates[name] = {**estimate,
                                      "budget": {"ok": True, "reason": "本地阶段"}}
            result.add("info", "resource_local",
                       "au60_build 为本地阶段，不占集群预算（%s 原子）" % n, stage=name)
            continue
        n = counts.get(name)
        if not n:
            result.add("warning", "resource_skipped",
                       "%s: 无原子数（结构产物未就绪），资源估算跳过" % name)
            continue
        estimate = estimate_smj_resources(name, int(n))
        budget = check_smj_budget(estimate)
        result.estimates[name] = {**estimate, "budget": budget}
        if estimate.get("capped"):
            result.add("error", "resource_capped",
                       "%s 预估超过资源上限（%s），已封顶 72h/128GiB: %s"
                       % (name, ", ".join(estimate["capped"]), estimate.get("basis")),
                       stage=name, capped=estimate["capped"])
        elif not budget.get("ok"):
            result.add("error", "resource_budget",
                       "%s 资源超预算: %s" % (name, budget.get("reason")), stage=name)
        else:
            result.add("info", "resource_ok",
                       "%s: walltime %s / --time %s / mem %s / cores %s"
                       % (name, estimate.get("walltime"), estimate.get("time_limit"),
                          estimate.get("mem"), estimate.get("cores")),
                       stage=name)


def _cluster_check(result: PreflightResult, config: Any, transport: Any,
                   *, require_remote: bool, timeout: float) -> None:
    mode = getattr(config, "mode", None)
    cluster = getattr(config, "cluster", None)
    if mode != "slurm-ssh":
        result.add("info", "remote_skipped",
                   "run.mode=%s：无集群提交，未执行 SSH/远端目录检查" % mode)
        return
    missing = []
    if not getattr(cluster, "remote_root", None):
        missing.append("remote_root")
    if not getattr(cluster, "host", None):
        missing.append("host")
    if not getattr(cluster, "partition", None):
        missing.append("partition")
    if missing:
        result.add("error", "cluster_incomplete",
                   "slurm-ssh 模式缺 cluster 配置: %s" % ", ".join(missing))
        return
    target = getattr(cluster, "ssh_target", None) or getattr(cluster, "host", "")
    result.notes.append("远端目标: %s；工作根: %s" % (target, cluster.remote_root))
    if transport is None:
        severity = "error" if require_remote else "warning"
        result.add(severity, "ssh_skipped",
                   "未提供传输层：**未执行远端检查**（不触网）。提交前请手工 "
                   "`ssh %s true` 确认免密登录" % target)
        return
    result.remote_checked = True
    probe = transport.run("true", timeout=timeout)
    if not probe.ok:
        result.add("error", "ssh_unreachable",
                   "SSH 连通性检查失败 (rc=%s): %s"
                   % (probe.code, (probe.err or probe.out).strip()[:300]))
        return
    root = shlex.quote(str(cluster.remote_root))
    exists = transport.run("test -d %s" % root, timeout=timeout)
    if not exists.ok:
        result.add("error", "remote_dir_missing",
                   "远端工作根不存在: %s" % cluster.remote_root)
        return
    writable = transport.run("test -w %s" % root, timeout=timeout)
    if not writable.ok:
        result.add("error", "remote_dir_not_writable",
                   "远端工作根不可写: %s" % cluster.remote_root)
        return
    result.add("info", "remote_dir_ok", "SSH 连通、远端工作根存在且可写: %s"
               % cluster.remote_root)


def _env_check(result: PreflightResult, config: Any, mode: str,
               env: Mapping[str, Any] | None) -> None:
    merged = build_env(config, env)
    if mode == "dry-run":
        result.add("info", "env_not_checked",
                   "dry-run 只生成输入：未核对 Gaussian / FHI-aims / aitranss 二进制与许可"
                   "（提交前在目标机器上再跑一次非 dry-run 预检）")
        return
    from smj._vendor.junction import preflight as vendored  # lazy: vendored
    for code, name, func in (
            ("aims_unavailable", "FHI-aims", vendored.aims_executable),
            ("gaussian_unavailable", "Gaussian", vendored.gaussian_license),
            ("aitranss_unavailable", "aitranss", vendored.aitranss_binary)):
        try:
            report = func(merged)
        except Exception as exc:  # noqa: BLE001
            result.add("error", code, "%s 检查失败: %s" % (name, exc))
            continue
        severity = "info" if report.get("ok") else "error"
        result.add(severity, code, str(report.get("detail", name)))
    result.notes.append("sbatch/调度由 smj.remote 负责，vendored sbatch 探测不适用")


# --------------------------------------------------------------------------- #
# 主入口
# --------------------------------------------------------------------------- #
def run_preflight(config: Any, *, transport: Any = None, stage: str | None = None,
                  dry_run: bool = False, layout: RunLayout | None = None,
                  elements: list[str] | None = None,
                  species_dir: str | Path | None = None,
                  natoms: int | None = None,
                  env: Mapping[str, Any] | None = None,
                  require_remote: bool = False,
                  check_structure: bool = True,
                  check_resources: bool = True,
                  check_env: bool = True,
                  timeout: float = 15.0) -> PreflightResult:
    """提交前预检。``config`` 为 :class:`smj.config.RunConfig`。

    - ``transport=None`` 时**不做**任何网络往返；``slurm-ssh`` 模式下明确记为
      "未执行远端检查"（warning；``require_remote=True`` 则升为 error）；
    - ``dry_run`` 只会**加强** dry-run 语义（配置本身就是 dry-run 时始终按 dry-run 处理）；
    - ``stage`` 可只检查单个阶段（资源/环境过滤；结构检查仍按 run 目录进行）。
    """
    run_dir = Path(getattr(config, "run_dir", "."))
    layout = layout or RunLayout(run_dir)
    mode = str(getattr(config, "mode", "dry-run") or "dry-run")
    effective_dry = mode == "dry-run" or bool(dry_run)
    result = PreflightResult(mode=mode, run=getattr(config, "name", run_dir.name))
    requested_stage = None
    if stage:
        from .state import normalize_stage
        requested_stage = normalize_stage(stage)
        result.notes.append("仅检查阶段: %s" % requested_stage)

    atoms = _molecule_check(result, config, elements)
    if elements is None and atoms:
        from smj._vendor.junction.geometry import elements_in_geometry  # lazy: vendored
        elements = elements_in_geometry(atoms)
    result.elements = list(elements or [])

    _species_check(result, config, species_dir, result.elements)
    if check_structure:
        _geometry_check(result, run_dir, layout)
    if check_resources:
        _resource_check(result, config, run_dir, layout, requested_stage, natoms)
    _cluster_check(result, config, transport, require_remote=require_remote,
                   timeout=timeout)
    if check_env:
        _env_check(result, config, mode if not effective_dry else "dry-run", env)
    result.notes.append("预检不触网、不提交；error 级 issue 阻断提交")
    return result


#: 兼容别名（CLI 的首选名是 run_preflight）。
preflight = run_preflight
check = run_preflight
