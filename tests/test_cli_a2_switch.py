"""CLI 开关 ``smj review --a2 / --no-a2``：只覆盖本次调用的 review.require_a2。

底线：

* 不加开关时沿用 ``smj.yaml``（默认 true）——缺 A2 记录必须 blocked；
* ``--no-a2`` 只关闭"A2 人工确认"这一项，不绕过结构几何闸门；
* ``--a2`` 能在 yaml 写 false 时强制开启；
* 开关不写回 ``smj.yaml``，也不会凭空生成 ``a2_confirmation.json``。
"""
from __future__ import annotations

import json
import math
from pathlib import Path

from smj.cli import main

# 表面三角的三个 (x, y) 偏移（约 2.88 A 边长，面积 > 0.5 A^2 下限）
TRI = ((0.0, 0.0), (2.88, 0.0), (1.44, 2.494))


def _junction_atoms() -> list[tuple[str, float, float, float]]:
    """BDT 式结：C6 环 (z=0)、两 S (z=±2.4)、两 Au tip (z=±4.8)，

    每侧再两层 Au（层距 2.35 A），最外层各 3 个原子（索引 14-16 / 20-22）。
    S-Au=2.4 A、Au-Au≥2.35 A，结构闸门（接触/表面侧别）应全部通过。
    """
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


def _write_case(workdir: Path, *, require_a2: bool | None = None,
                metadata: dict | None = None) -> Path:
    """建一个最小 run：<workdir>/run/{electrode_build,review.json}。"""
    build = workdir / "run" / "electrode_build"
    build.mkdir(parents=True)
    (build / "geometry.in").write_text(
        "".join("atom %.8f %.8f %.8f %s\n" % (x, y, z, e)
                for e, x, y, z in _junction_atoms()),
        encoding="utf-8")
    (build / "electrode_metadata.json").write_text(
        json.dumps(metadata if metadata is not None else _metadata(), indent=2),
        encoding="utf-8")

    xyz = workdir / "molecule.xyz"
    xyz.write_text("2\nprobe\nS 0.0 0.0 -1.0\nS 0.0 0.0 1.0\n", encoding="utf-8")
    review_block = ""
    if require_a2 is not None:
        review_block = "review:\n  require_a2: %s\n" % (
            "true" if require_a2 else "false")
    cfg = workdir / "smj.yaml"
    cfg.write_text(
        "run:\n  name: run\n  dir: %s\n  mode: dry-run\n"
        "molecule:\n  xyz: %s\n%s"
        % (workdir.resolve().as_posix(), xyz.resolve().as_posix(), review_block),
        encoding="utf-8")
    return cfg


def _payload(workdir: Path) -> dict:
    return json.loads(
        (workdir / "run" / "review.json").read_text(encoding="utf-8"))


def test_default_uses_config_and_blocks(workdir: Path):
    cfg = _write_case(workdir)
    assert main(["review", "-c", str(cfg)]) == 1
    payload = _payload(workdir)
    assert payload["config"]["require_a2"] is True
    assert payload["status"] == "blocked"
    assert any(issue["code"] == "a2_missing" for issue in payload["issues"])


def test_no_a2_turns_off_for_this_call(workdir: Path):
    cfg = _write_case(workdir)
    assert main(["review", "--no-a2", "-c", str(cfg)]) == 0
    payload = _payload(workdir)
    assert payload["config"]["require_a2"] is False
    # 无 A2 记录不应再产生 error；run_spec 缺失只会把结论降为 warning。
    assert payload["ok"] is True
    assert payload["status"] in ("passed", "warning")
    assert not any(issue["severity"] == "error" for issue in payload["issues"])
    assert any(issue["code"] == "a2_not_required" for issue in payload["issues"])
    # 只影响本次调用：不写回配置，也不生成签署记录
    assert "require_a2" not in cfg.read_text(encoding="utf-8")
    assert not (workdir / "run" / "electrode_build"
                / "a2_confirmation.json").exists()


def test_a2_flag_forces_on_when_config_off(workdir: Path):
    cfg = _write_case(workdir, require_a2=False)
    assert main(["review", "--a2", "-c", str(cfg)]) == 1
    payload = _payload(workdir)
    assert payload["config"]["require_a2"] is True
    assert payload["status"] == "blocked"


def test_no_a2_does_not_bypass_structure_gate(workdir: Path):
    bad = _metadata()
    bad["left_surface_indices"] = [1, 2, 3]   # C 原子冒充表面 Au -> 结构闸门必须报错
    cfg = _write_case(workdir, metadata=bad)
    assert main(["review", "--no-a2", "-c", str(cfg)]) == 1
    payload = _payload(workdir)
    assert payload["status"] == "blocked"
    assert any(str(issue.get("code", "")).startswith("structure:")
               for issue in payload["issues"])
