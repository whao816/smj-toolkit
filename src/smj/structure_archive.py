"""结构归档：把每次 prepare/collect 看到或产出的 xyz 与计算结果一起存盘。

设计口径（2026-11-18）：

* 归档目录固定在 ``<run>/saved_structures/``——与 run 同级，跨阶段可查，且
  ``.gitignore`` 已排除，不会把运行期结构误提交进仓库；
* 每个条目在 ``manifest.json`` 里登记 ``stage`` / ``kind`` / ``engine`` /
  ``sha256`` / ``natoms`` / ``source`` / ``created_at``，并可选挂一份
  ``<seq>_<stage>_<kind>.results.json`` 结果快照；
* 结构文件按字节原样复制（不改注释、不改坐标），因此 ``sha256`` 可与源文件
  逐字节核对；"与计算结果挂钩"通过 manifest 的 ``results_file`` +
  ``results_sha256`` + ``summary`` 以及结果快照本身保证。

归档是审计动作，不是科学判据：写入失败只降级为阶段 note，不阻断
prepare/collect（与 species 审计同一口径）。
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

ARCHIVE_DIRNAME = "saved_structures"
MANIFEST_NAME = "manifest.json"
SCHEMA = 1

#: 结果快照里保留的标量键（其余结构如 geometry_check/scf 仍留在快照文件里）。
_SUMMARY_KEYS = (
    "engine",
    "normal_termination",
    "converged",
    "relax_steps",
    "max_force_ev_per_A",
    "force_threshold_ev_per_A",
    "final_scf_energy_hartree",
    "imaginary_frequency_count",
    "electron_count",
    "atom_count",
    "final_atom_count",
    "surface_fingerprint",
    "surface_indices_verified",
)

#: 归档格式 → 落盘后缀。``aims`` 是 FHI-aims 的 geometry.in（electrode_build /
#: aims_tip_relax / au_relax 的产物），``xyz`` 是优化阶段的 opt.xyz。
_FORMAT_SUFFIX = {"xyz": ".xyz", "aims": ".in"}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    """先写临时文件再 replace，避免中断留下半个 manifest。"""
    text = json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n"
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)


def _load_manifest(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {"schema": SCHEMA, "entries": []}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        # 损坏的 manifest 不阻塞归档：另存一份 .bad 再重建，保留现场。
        path.replace(path.with_name(path.name + ".bad"))
        return {"schema": SCHEMA, "entries": []}
    if not isinstance(data, dict) or not isinstance(data.get("entries"), list):
        path.replace(path.with_name(path.name + ".bad"))
        return {"schema": SCHEMA, "entries": []}
    data.setdefault("schema", SCHEMA)
    return data


def _atom_count(structure: Path, fmt: str = "xyz") -> int:
    """数原子。优先用解析器（xyz / aims 各有语义），失败再退回行数。"""
    try:
        if fmt == "aims":
            from ._vendor.junction.geometry import read_aims_geometry
            return len(read_aims_geometry(structure))
        from ._vendor.junction.geometry import read_xyz
        return len(read_xyz(structure))
    except Exception:  # noqa: BLE001 - 归档不因解析失败而中断
        pass
    try:
        lines = structure.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return 0
    if fmt == "aims":
        return sum(1 for line in lines
                   if line.split()[:1] in (["atom"], ["atom_frac"]))
    try:
        return int(lines[0].strip()) if lines else 0
    except (TypeError, ValueError):
        return max(len(lines) - 2, 0)


def _source_label(source: Path, run_dir: Path) -> str:
    try:
        return source.resolve().relative_to(run_dir.resolve()).as_posix()
    except ValueError:
        return str(source)


def summarize_results(results: Mapping[str, Any], stage: str | None = None) -> dict[str, Any]:
    """从阶段 results 里摘出与结构对应的标量结果（供 manifest 快速检索）。

    结果块的键名随阶段/引擎不同（``molecule_opt`` / ``gaussian`` /
    ``electrode_build`` / ``aims_tip_relax`` / ``au_relax``）：优先按 ``stage``
    取，其次按两套分子优化引擎取，最后取任意第一个映射块；非标量结构
    （几何检查、SCF 明细、表面索引列表）留在结果快照里。
    """
    candidates: list[str] = []
    if stage:
        candidates.append(stage)
    candidates.extend(("molecule_opt", "gaussian"))
    candidates.extend(str(key) for key in results)
    seen: set[str] = set()
    for key in candidates:
        if key in seen:
            continue
        seen.add(key)
        block = results.get(key)
        if isinstance(block, Mapping):
            summary = {name: block[name] for name in _SUMMARY_KEYS if name in block}
            summary["results_key"] = key
            return summary
    return {}


def archive_structure(
    run_dir: str | Path,
    source: str | Path,
    *,
    stage: str,
    kind: str,
    engine: str | None = None,
    results: Mapping[str, Any] | None = None,
    fmt: str = "xyz",
    natoms: int | None = None,
) -> dict[str, Any]:
    """把 ``source`` 原样复制进 ``<run>/saved_structures/`` 并登记 manifest。

    参数
    ----
    run_dir: run 根目录（``ctx.layout.run_dir``）。
    source: 要归档的结构文件；必须存在且可读。
    stage / kind: 例如 ``molecule_opt`` / ``input``、``electrode_build`` / ``junction``。
    engine: ``aims`` / ``gaussian``；仅为可读性记录。
    results: 阶段 results 字典；给出时额外写一份结果快照并记 sha256。
    fmt: ``xyz``（优化产物）或 ``aims``（FHI-aims geometry.in 原生格式）。
    natoms: 已知原子数时直接传入（阶段刚校验过），否则按 ``fmt`` 解析。

    返回新增的 manifest 条目。序号只增不减：重复 prepare/collect 也是一次
    独立的审计记录，便于区分"哪一次优化/装配"产出了哪个结构。
    """
    run_path = Path(run_dir)
    src = Path(source)
    if not src.is_file():
        raise FileNotFoundError("结构文件不存在: %s" % src)
    if fmt not in _FORMAT_SUFFIX:
        raise ValueError("未知的结构归档格式: %r（支持 %s）"
                         % (fmt, "/".join(_FORMAT_SUFFIX)))

    root = run_path / ARCHIVE_DIRNAME
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / MANIFEST_NAME
    manifest = _load_manifest(manifest_path)
    entries = manifest["entries"]
    seq = len(entries) + 1
    stem = "%04d_%s_%s" % (seq, stage, kind)
    dest = root / (stem + _FORMAT_SUFFIX[fmt])

    data = src.read_bytes()
    dest.write_bytes(data)
    atom_count = natoms if natoms is not None else _atom_count(dest, fmt)

    entry: dict[str, Any] = {
        "seq": seq,
        "stage": stage,
        "kind": kind,
        "engine": engine,
        "format": fmt,
        "file": dest.name,
        "sha256": _sha256_bytes(data),
        "natoms": atom_count,
        "source": _source_label(src, run_path),
        "created_at": _now(),
        "results_file": None,
        "results_sha256": None,
        "summary": {},
    }
    if results is not None:
        payload = (json.dumps(dict(results), indent=2, ensure_ascii=False,
                              sort_keys=True) + "\n").encode("utf-8")
        snapshot = root / (stem + ".results.json")
        snapshot.write_bytes(payload)
        entry["results_file"] = snapshot.name
        entry["results_sha256"] = _sha256_bytes(payload)
        entry["summary"] = summarize_results(results, stage=stage)

    entries.append(entry)
    manifest["run"] = run_path.name
    manifest["updated_at"] = _now()
    _write_json(manifest_path, manifest)
    return entry


def safe_archive(
    run_dir: str | Path,
    source: str | Path,
    *,
    stage: str,
    kind: str,
    engine: str | None = None,
    results: Mapping[str, Any] | None = None,
    fmt: str = "xyz",
    natoms: int | None = None,
) -> list[str]:
    """归档并转成阶段 notes；任何异常只降级为一行 note，不阻断阶段。

    归档是审计动作，不是科学判据：写盘失败（磁盘满、权限、manifest 损坏）
    都不应让已经算完的结构与结果丢掉。
    """
    try:
        entry = archive_structure(
            run_dir, source, stage=stage, kind=kind, engine=engine,
            results=results, fmt=fmt, natoms=natoms)
    except Exception as exc:  # noqa: BLE001 - 归档失败不能吞掉科学结果
        return ["结构归档失败（不影响本阶段）: %s" % exc]
    return [archive_note(entry)]


def archive_note(entry: Mapping[str, Any]) -> str:
    """给阶段 note 用的一行摘要。"""
    text = "结构已归档: %s/%s" % (ARCHIVE_DIRNAME, entry.get("file"))
    if entry.get("results_file"):
        text += "（结果快照 %s）" % entry["results_file"]
    return text
