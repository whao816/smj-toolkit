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
)


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


def _atom_count(xyz: Path) -> int:
    try:
        return len(xyz.read_text(encoding="utf-8").strip().splitlines()) - 2
    except OSError:
        return 0


def _source_label(source: Path, run_dir: Path) -> str:
    try:
        return source.resolve().relative_to(run_dir.resolve()).as_posix()
    except ValueError:
        return str(source)


def summarize_results(results: Mapping[str, Any]) -> dict[str, Any]:
    """从阶段 results 里摘出与结构对应的标量结果（供 manifest 快速检索）。

    两种引擎的结果块键名不同（``molecule_opt`` / ``gaussian``），这里按出现
    顺序取第一个映射块；其余非标量结构（几何检查、SCF 明细）留在结果快照里。
    """
    for key in ("molecule_opt", "gaussian"):
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
) -> dict[str, Any]:
    """把 ``source`` 原样复制进 ``<run>/saved_structures/`` 并登记 manifest。

    参数
    ----
    run_dir: run 根目录（``ctx.layout.run_dir``）。
    source: 要归档的 xyz；必须存在且可读。
    stage / kind: 例如 ``molecule_opt`` / ``input``、``molecule_opt`` / ``optimized``。
    engine: ``aims`` / ``gaussian``；仅为可读性记录。
    results: 阶段 results 字典；给出时额外写一份结果快照并记 sha256。

    返回新增的 manifest 条目。序号只增不减：重复 prepare/collect 也是一次
    独立的审计记录，便于区分"哪一次优化"产出了哪个结构。
    """
    run_path = Path(run_dir)
    src = Path(source)
    if not src.is_file():
        raise FileNotFoundError("结构文件不存在: %s" % src)

    root = run_path / ARCHIVE_DIRNAME
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / MANIFEST_NAME
    manifest = _load_manifest(manifest_path)
    entries = manifest["entries"]
    seq = len(entries) + 1
    stem = "%04d_%s_%s" % (seq, stage, kind)
    dest = root / (stem + ".xyz")

    data = src.read_bytes()
    dest.write_bytes(data)

    entry: dict[str, Any] = {
        "seq": seq,
        "stage": stage,
        "kind": kind,
        "engine": engine,
        "file": dest.name,
        "sha256": _sha256_bytes(data),
        "natoms": _atom_count(dest),
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
        entry["summary"] = summarize_results(results)

    entries.append(entry)
    manifest["run"] = run_path.name
    manifest["updated_at"] = _now()
    _write_json(manifest_path, manifest)
    return entry


def archive_note(entry: Mapping[str, Any]) -> str:
    """给阶段 note 用的一行摘要。"""
    text = "结构已归档: %s/%s" % (ARCHIVE_DIRNAME, entry.get("file"))
    if entry.get("results_file"):
        text += "（结果快照 %s）" % entry["results_file"]
    return text
