# -*- coding: utf-8 -*-
"""Export the small, git-tracked slice of the 2026-11-22 validation runs.

``runs/`` is gitignored, so the evidence behind
``docs/validation-light-vs-tight-2026-11-22.md`` and
``docs/validation-vdw-onoff-2026-11-22.md`` otherwise lives only on the
machine that produced it.  This script copies the reviewable artifacts of the
documented arms into ``docs/validation-data/`` and writes a sha256 manifest:

    python scripts/export_validation_data.py           # (re)export
    python scripts/export_validation_data.py --check   # verify repo copy

Copied per arm: ``state.json``, ``run_spec.json``, stage ``input_summary.json``
/ ``control.in`` / geometries, ``saved_structures/``, ``TE.dat`` and
``junction_analyze/results.json``.  ``aims.out`` (tens of MB), PDF/PNG renders
and comparison plots are intentionally excluded.  The abandoned az-thiolate
arms are not exported; they live under
``runs/light-tight-2026-11-22/_abandoned/``.
"""
from __future__ import annotations

import argparse
import fnmatch
import hashlib
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs" / "validation-data"
SCHEMA = "smj.validation-data.v1"

#: Documented arms only.  The abandoned az-thiolate arms are deliberately absent.
SETS: dict[str, dict[str, object]] = {
    "light-tight-2026-11-22": {
        "arms": ["2cn-light", "2cn-tight", "az-dithiol-light", "az-dithiol-tight"],
        "run_prefix": "lt-",
    },
    "vdw-2026-11-22": {
        "arms": ["2cn-novdw", "az-dithiol-novdw"],
        "run_prefix": "vdw-",
    },
}

ARM_FILES = ("molecule.xyz", "smj.yaml")

RUN_FILES = frozenset(
    {
        "run_spec.json",
        "state.json",
        "molecule_opt/control.in",
        "molecule_opt/geometry.in",
        "molecule_opt/geometry.in.next_step",
        "molecule_opt/opt.xyz",
        "molecule_opt/input_summary.json",
        "molecule_opt/species.json",
        "aims_tip_relax/control.in",
        "aims_tip_relax/geometry.in",
        "aims_tip_relax/geometry.in.next_step",
        "aims_tip_relax/review_tips.xyz",
        "aims_tip_relax/input_summary.json",
        "electrode_build/geometry.in",
        "electrode_build/review_junction.xyz",
        "electrode_build/connection_check.json",
        "electrode_build/connection_check_audit.json",
        "electrode_build/electrode_metadata.json",
        "aims_scf/control.in",
        "aims_scf/geometry.in",
        "aims_scf/input_summary.json",
        "aitranss_te/input_summary.json",
        "aitranss_te/TE.dat",
        "aitranss_te/tcontrol",
        "aitranss_te/self.energy.in",
        "junction_analyze/results.json",
    }
)
RUN_GLOBS = ("saved_structures/*", "electrodes/*.aims")

README = """# 验证数据（2026-11-22）

本目录是下列两份验证记录所引用 run 的**精简导出**，用于在没有 `runs/` 的
机器上核对文档里的数字，也是这些数据进入版本库的唯一途径（`runs/` 被
`.gitignore` 忽略）。

- [light vs tight 整链对照](../validation-light-vs-tight-2026-11-22.md)
- [vdW 开关整链对照](../validation-vdw-onoff-2026-11-22.md)

## 内容与范围

每个臂保留：`state.json`、`run_spec.json`、各阶段的 `input_summary.json`、
`control.in`、输入/优化几何、`saved_structures/`、`aitranss_te/TE.dat` 与
`junction_analyze/results.json`。文档里的 L0 / S / κe / ZT 可直接由
`TE.dat` + `results.json` 复算。

**不包含**：`aims.out`（单臂 10–25 MB）、PDF/PNG 渲染和对照图（图在
[`../figures/`](../figures/)）。本目录也不含未跑完的
`az-thiolate-*` 臂——它们已移到 `runs/light-tight-2026-11-22/_abandoned/`，
状态是 SCF 不收敛、无 `results.json`，不属于对照数据。

## 导出与校验

```bash
python scripts/export_validation_data.py           # 需要本机有 runs/，重新导出
python scripts/export_validation_data.py --check   # 校验 MANIFEST 与入库副本
```

`--check` 先核对入库副本的 sha256/字节数，再在 `runs/` 存在时逐字节比对
源文件；两边一致才算通过。`MANIFEST.json` 记录每个文件的来源路径与哈希，
`docs/` 其余文件请见 [文档索引](../README.md)。
"""


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def _selected(rel: str) -> bool:
    if rel in RUN_FILES:
        return True
    return any(fnmatch.fnmatchcase(rel, pattern) for pattern in RUN_GLOBS)


def _find_run_dir(arm_dir: Path, prefix: str) -> Path:
    candidates = [
        child
        for child in arm_dir.iterdir()
        if child.is_dir() and child.name.startswith(prefix)
    ]
    if len(candidates) != 1:
        names = ", ".join(sorted(c.name for c in candidates)) or "(none)"
        raise SystemExit(
            f"error: expected exactly one {prefix}* run dir under "
            f"{arm_dir.relative_to(ROOT)}, found: {names}"
        )
    return candidates[0]


def _collect() -> tuple[Path, list[dict[str, object]]]:
    records: list[dict[str, object]] = []
    for set_name, spec in SETS.items():
        set_dir = ROOT / "runs" / set_name
        if not set_dir.is_dir():
            raise SystemExit(f"error: missing source set {_rel(set_dir)}")
        for arm in spec["arms"]:  # type: ignore[union-attr]
            arm_dir = set_dir / arm
            if not arm_dir.is_dir():
                raise SystemExit(f"error: missing arm {_rel(arm_dir)}")
            run_dir = _find_run_dir(arm_dir, str(spec["run_prefix"]))
            pairs: list[tuple[Path, Path]] = []
            for name in ARM_FILES:
                src = arm_dir / name
                if src.is_file():
                    pairs.append((src, DATA / set_name / arm / name))
            for src in sorted(run_dir.rglob("*")):
                if not src.is_file():
                    continue
                rel = src.relative_to(run_dir).as_posix()
                if _selected(rel):
                    pairs.append(
                        (src, DATA / set_name / arm / run_dir.name / rel)
                    )
            if not pairs:
                raise SystemExit(f"error: nothing selected under {_rel(run_dir)}")
            for src, dst in pairs:
                data = src.read_bytes()
                records.append(
                    {
                        "path": _rel(dst),
                        "bytes": len(data),
                        "sha256": _sha256(data),
                        "source": _rel(src),
                    }
                )
    records.sort(key=lambda rec: str(rec["path"]))
    return DATA, records


def export() -> int:
    if DATA.exists():
        manifest = DATA / "MANIFEST.json"
        if not manifest.is_file():
            raise SystemExit(
                f"error: {_rel(DATA)} exists without MANIFEST.json; "
                "refusing to overwrite an unknown tree"
            )
        old = json.loads(manifest.read_text(encoding="utf-8"))
        if old.get("schema") != SCHEMA:
            raise SystemExit(
                f"error: {_rel(manifest)} has schema {old.get('schema')!r}, "
                f"expected {SCHEMA!r}"
            )
        shutil.rmtree(DATA)

    data_dir, records = _collect()
    for record in records:
        src = ROOT / str(record["source"])
        dst = ROOT / str(record["path"])
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)

    manifest = {
        "schema": SCHEMA,
        "sets": {
            name: {
                "arms": list(spec["arms"]),  # type: ignore[call-overload]
                "source": f"runs/{name}",
            }
            for name, spec in SETS.items()
        },
        "files": records,
    }
    (DATA / "MANIFEST.json").write_bytes(
        (json.dumps(manifest, ensure_ascii=False, indent=2) + "\n").encode("utf-8")
    )
    (DATA / "README.md").write_bytes(README.encode("utf-8"))

    total = sum(int(record["bytes"]) for record in records)
    print(
        f"exported {len(records)} files, {total} bytes "
        f"from {len(SETS)} set(s) to {_rel(DATA)}"
    )
    for name in SETS:
        count = sum(
            1 for record in records if str(record["path"]).startswith(f"docs/validation-data/{name}/")
        )
        print(f"  {name}: {count} files")
    return 0


def check() -> int:
    manifest_path = DATA / "MANIFEST.json"
    if not manifest_path.is_file():
        print(f"[fail] missing {_rel(manifest_path)}; run the exporter first")
        return 1
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != SCHEMA:
        print(f"[fail] unsupported schema {manifest.get('schema')!r}")
        return 1

    errors: list[str] = []
    listed: set[str] = set()
    checked_sources = 0
    missing_sources = 0
    for record in manifest["files"]:
        rel = str(record["path"])
        listed.add(rel)
        path = ROOT / rel
        if not path.is_file():
            errors.append(f"missing file: {rel}")
            continue
        data = path.read_bytes()
        if len(data) != record["bytes"]:
            errors.append(
                f"size mismatch: {rel} ({len(data)} != {record['bytes']})"
            )
        if _sha256(data) != record["sha256"]:
            errors.append(f"sha256 mismatch: {rel}")
        source = ROOT / str(record["source"])
        if source.is_file():
            checked_sources += 1
            if source.read_bytes() != data:
                errors.append(f"repo copy differs from source run: {rel}")
        else:
            missing_sources += 1

    for path in sorted(DATA.rglob("*")):
        if path.is_file():
            rel = _rel(path)
            if rel not in listed and path.name not in {"MANIFEST.json", "README.md"}:
                errors.append(f"unlisted file: {rel}")

    for issue in errors:
        print("[fail]", issue)
    if errors:
        return 1

    source_note = (
        f"all {checked_sources} source files match"
        if missing_sources == 0
        else f"{checked_sources} source files match, {missing_sources} not present on this machine"
    )
    print(
        f"[ok] {len(listed)} files match MANIFEST.json ({_rel(DATA)}); {source_note}"
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="verify docs/validation-data against its manifest (and runs/ if present)",
    )
    args = parser.parse_args(argv)
    return check() if args.check else export()


if __name__ == "__main__":
    sys.exit(main())
