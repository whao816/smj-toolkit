"""协议漂移只来自 ``molecule.family`` 时不得阻断复核。

``family`` 只是判据标定标注（决定哪些判据适用），不改变物理协议。
典型场景：run 做完才登记族（``null -> BPY``），run_spec 记录的
``protocol_sha256`` 必然变；若照旧报 ``provenance_drift``，整条 run 会被
迫重建。这里固定两条底线：

* 仅 family 不同 -> ``provenance_family_relabel`` warning，复核通过；
* 任何物理参数不同 -> 仍是 ``provenance_drift`` error，复核 blocked。
"""
from __future__ import annotations

import json
import math
from pathlib import Path

from smj.cli import main
from smj.config import load_config
from smj.provenance import run_spec_payload
from smj.review import _NOT_FAMILY_ONLY, _family_drift

TRI = ((0.0, 0.0), (2.88, 0.0), (1.44, 2.494))


def _junction_atoms() -> list[tuple[str, float, float, float]]:
    atoms: list[tuple[str, float, float, float]] = []
    for k in range(6):
        ang = math.radians(60.0 * k)
        atoms.append(("C", 1.39 * math.cos(ang), 1.39 * math.sin(ang), 0.0))
    atoms.append(("S", 0.0, 0.0, -2.4))    # 7  左 linker
    atoms.append(("S", 0.0, 0.0, 2.4))     # 8  右 linker
    atoms.append(("Au", 0.0, 0.0, -4.8))   # 9  左 tip
    atoms.append(("Au", 0.0, 0.0, 4.8))    # 10 右 tip
    for z in (-7.15, -9.5):
        for dx, dy in TRI:
            atoms.append(("Au", dx, dy, z))  # 11-16
    for z in (7.15, 9.5):
        for dx, dy in TRI:
            atoms.append(("Au", dx, dy, z))  # 17-22
    return atoms


def _metadata() -> dict:
    return {
        "left_tip_index": 9,
        "right_tip_index": 10,
        "left_linker_index": 7,
        "right_linker_index": 8,
        "left_surface_indices": [14, 15, 16],
        "right_surface_indices": [20, 21, 22],
        "surface_rule": "det-v1",
    }


def _write_case(workdir: Path, *, family: str | None = None,
                tip_distance: float = 2.35) -> Path:
    """最小 run：<workdir>/run/{electrode_build,smj.yaml}。"""
    build = workdir / "run" / "electrode_build"
    build.mkdir(parents=True, exist_ok=True)
    (build / "geometry.in").write_text(
        "".join("atom %.8f %.8f %.8f %s\n" % (x, y, z, e)
                for e, x, y, z in _junction_atoms()),
        encoding="utf-8")
    (build / "electrode_metadata.json").write_text(
        json.dumps(_metadata(), indent=2), encoding="utf-8")

    xyz = workdir / "molecule.xyz"
    xyz.write_text("2\nprobe\nS 0.0 0.0 -1.0\nS 0.0 0.0 1.0\n",
                   encoding="utf-8")
    family_line = ""
    if family is not None:
        family_line = "  family: %s\n" % family
    cfg = workdir / "smj.yaml"
    cfg.write_text(
        "run:\n  name: run\n  dir: %s\n  mode: dry-run\n"
        "molecule:\n  xyz: %s\n%s"
        "electrode:\n  au_tip_distance: %s\n"
        % (workdir.resolve().as_posix(), xyz.resolve().as_posix(),
           family_line, tip_distance),
        encoding="utf-8")
    return cfg


def _write_run_spec(cfg: Path, run_dir: Path) -> dict:
    spec = run_spec_payload(load_config(cfg), run_dir)
    (run_dir / "run_spec.json").write_text(
        json.dumps(spec, ensure_ascii=False, indent=2), encoding="utf-8")
    return spec


def _review(workdir: Path, cfg: Path) -> tuple[int, dict]:
    rc = main(["review", "-c", str(cfg), "--force"])
    payload = json.loads(
        (workdir / "run" / "review.json").read_text(encoding="utf-8"))
    return rc, payload


def _codes(payload: dict) -> list[str]:
    return [str(issue.get("code", "")) for issue in payload["issues"]]


# --------------------------------------------------------------------------- #
# 判定函数（单元）
# --------------------------------------------------------------------------- #
def test_family_only_drift_detected(workdir: Path):
    cfg = _write_case(workdir, family=None)
    spec = _write_run_spec(cfg, workdir / "run")

    relabeled = load_config(_write_case(workdir, family="BPY"))
    assert _family_drift(relabeled, spec) == (None, "BPY")


def test_physical_drift_not_tolerated(workdir: Path):
    cfg = _write_case(workdir, family=None)
    spec = _write_run_spec(cfg, workdir / "run")

    changed = load_config(_write_case(workdir, family="BPY", tip_distance=2.60))
    assert _family_drift(changed, spec) is _NOT_FAMILY_ONLY


def test_spec_without_config_snapshot_not_tolerated(workdir: Path):
    cfg = _write_case(workdir, family=None)
    _write_run_spec(cfg, workdir / "run")
    relabeled = load_config(_write_case(workdir, family="BPY"))
    assert _family_drift(relabeled, {"protocol_sha256": "deadbeef"}) is _NOT_FAMILY_ONLY


# --------------------------------------------------------------------------- #
# 端到端（CLI review 闸门）
# --------------------------------------------------------------------------- #
def test_review_passes_when_only_family_changed(workdir: Path):
    cfg = _write_case(workdir, family=None)
    _write_run_spec(cfg, workdir / "run")

    relabeled = _write_case(workdir, family="BPY")
    rc, payload = _review(workdir, relabeled)
    codes = _codes(payload)
    assert rc == 0
    assert payload["ok"] is True
    assert payload["status"] in ("passed", "warning")
    assert "provenance_family_relabel" in codes
    assert "provenance_drift" not in codes


def test_review_blocks_physical_drift(workdir: Path):
    cfg = _write_case(workdir, family=None)
    _write_run_spec(cfg, workdir / "run")

    changed = _write_case(workdir, family="BPY", tip_distance=2.60)
    rc, payload = _review(workdir, changed)
    codes = _codes(payload)
    assert rc == 1
    assert payload["ok"] is False
    assert "provenance_drift" in codes
    assert "provenance_family_relabel" not in codes
