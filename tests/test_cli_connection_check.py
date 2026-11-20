"""CLI 连接核对（connection check，原 A2）：默认自动登记 + 开关只覆盖本次调用。

底线（2026-11-15 用户要求：改个易懂的名字、默认自动审核）：

* 默认（``auto_connection_check: true``）不带 ``--approve`` 也能自动登记，
  记录 ``who="auto"``，review 直接通过；
* ``auto_connection_check: false`` + ``connection_check: true`` -> 缺记录必须 blocked；
* ``--no-connection-check`` 只关闭这一项闸门，不绕过结构几何闸门；
* ``--connection-check`` 能在 yaml 写 false 时强制开启；
* 旧开关 ``--surface-check / --no-surface-check``、``--a2 / --no-a2`` 与旧配置键
  ``surface_check*`` / ``require_a2`` / ``a2_deviation`` 仍生效（带迁移 warning）；
* 旧 run 的 ``surface_check.json`` / ``a2_confirmation.json`` 仍被识别；
* 开关/自动登记都不写回 ``smj.yaml``。
"""
from __future__ import annotations

import json
import math
from pathlib import Path

from smj.cli import main
from smj.config import load_config

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


def _write_case(workdir: Path, *, review_lines: list[str] | None = None,
                metadata: dict | None = None) -> Path:
    """建一个最小 run：<workdir>/run/{electrode_build,smj.yaml}。"""
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
    if review_lines:
        review_block = "review:\n" + "".join("  %s\n" % line for line in review_lines)
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


def _codes(payload: dict) -> list[str]:
    return [str(issue.get("code", "")) for issue in payload["issues"]]


def _confirm(workdir: Path) -> Path:
    return workdir / "run" / "electrode_build" / "connection_check.json"


# --------------------------------------------------------------------------- #
# 默认自动
# --------------------------------------------------------------------------- #
def test_default_auto_registers_and_passes(workdir: Path):
    cfg = _write_case(workdir)
    assert main(["review", "-c", str(cfg)]) == 0
    payload = _payload(workdir)
    assert payload["config"]["connection_check"] is True
    assert payload["config"]["auto_connection_check"] is True
    assert payload["ok"] is True
    assert payload["status"] in ("passed", "warning")
    assert "connection_check_auto_signed" in _codes(payload)
    assert payload["connection_check"]["present"] is True
    assert payload["connection_check"]["valid"] is True
    record = json.loads(_confirm(workdir).read_text(encoding="utf-8"))
    assert record["who"] == "auto"
    # 自动登记不写回 smj.yaml
    assert "auto_connection_check" not in cfg.read_text(encoding="utf-8")


def test_auto_off_blocks_when_gate_on(workdir: Path):
    cfg = _write_case(workdir, review_lines=[
        "connection_check: true", "auto_connection_check: false"])
    assert main(["review", "-c", str(cfg)]) == 1
    payload = _payload(workdir)
    assert payload["status"] == "blocked"
    assert "connection_check_missing" in _codes(payload)
    assert not _confirm(workdir).exists()


# --------------------------------------------------------------------------- #
# 本次调用开关
# --------------------------------------------------------------------------- #
def test_no_connection_check_turns_off_for_this_call(workdir: Path):
    cfg = _write_case(workdir, review_lines=["auto_connection_check: false"])
    assert main(["review", "--no-connection-check", "-c", str(cfg)]) == 0
    payload = _payload(workdir)
    assert payload["config"]["connection_check"] is False
    assert payload["ok"] is True
    assert "connection_check_disabled" in _codes(payload)
    assert not any(issue["severity"] == "error" for issue in payload["issues"])
    assert not _confirm(workdir).exists()
    # 本次调用的开关不写回配置（yaml 里只有原本的 auto_connection_check）
    assert "\n  connection_check:" not in cfg.read_text(encoding="utf-8")


def test_connection_check_flag_forces_on_when_config_off(workdir: Path):
    cfg = _write_case(workdir, review_lines=["connection_check: false"])
    assert main(["review", "--connection-check", "-c", str(cfg)]) == 0
    payload = _payload(workdir)
    assert payload["config"]["connection_check"] is True
    assert payload["connection_check"]["valid"] is True


def test_no_connection_check_does_not_bypass_structure_gate(workdir: Path):
    bad = _metadata()
    bad["left_surface_indices"] = [1, 2, 3]   # C 原子冒充表面 Au -> 结构闸门必须报错
    cfg = _write_case(workdir, review_lines=["auto_connection_check: false"],
                      metadata=bad)
    assert main(["review", "--no-connection-check", "-c", str(cfg)]) == 1
    payload = _payload(workdir)
    assert payload["status"] == "blocked"
    assert any(code.startswith("structure:") for code in _codes(payload))


# --------------------------------------------------------------------------- #
# 兼容：旧开关 / 旧配置键 / 旧记录文件
# --------------------------------------------------------------------------- #
def test_legacy_no_a2_flag_still_works(workdir: Path):
    cfg = _write_case(workdir, review_lines=["auto_connection_check: false"])
    assert main(["review", "--no-a2", "-c", str(cfg)]) == 0
    assert _payload(workdir)["config"]["connection_check"] is False


def test_legacy_a2_flag_still_works(workdir: Path):
    cfg = _write_case(workdir, review_lines=["connection_check: false"])
    assert main(["review", "--a2", "-c", str(cfg)]) == 0
    assert _payload(workdir)["config"]["connection_check"] is True


def test_legacy_config_keys_parse_with_warning(workdir: Path):
    cfg = _write_case(workdir, review_lines=[
        "require_a2: false", "a2_deviation: 0.2"])
    config = load_config(cfg)
    assert config.review.connection_check is False
    assert config.review.connection_check_deviation == 0.2
    assert any("require_a2" in w and "connection_check" in w
               for w in config.warnings)
    assert any("a2_deviation" in w for w in config.warnings)


def test_legacy_confirmation_file_is_honored(workdir: Path):
    from smj.junction.fingerprint import record_confirmation

    cfg = _write_case(workdir, review_lines=["auto_connection_check: false"])
    build = workdir / "run" / "electrode_build"
    record_confirmation(build, who="王昊")
    (build / "connection_check.json").rename(build / "a2_confirmation.json")

    assert main(["review", "-c", str(cfg)]) == 0
    payload = _payload(workdir)
    assert payload["connection_check"]["present"] is True
    assert payload["connection_check"]["valid"] is True


def test_legacy_surface_check_flag_still_works(workdir: Path):
    cfg = _write_case(workdir, review_lines=["auto_connection_check: false"])
    assert main(["review", "--no-surface-check", "-c", str(cfg)]) == 0
    assert _payload(workdir)["config"]["connection_check"] is False


def test_legacy_surface_config_keys_parse_with_warning(workdir: Path):
    cfg = _write_case(workdir, review_lines=[
        "surface_check: false", "surface_check_deviation: 0.3",
        "auto_surface_check: false"])
    config = load_config(cfg)
    assert config.review.connection_check is False
    assert config.review.connection_check_deviation == 0.3
    assert config.review.auto_connection_check is False
    joined = " ".join(config.warnings)
    for old, new in (("surface_check", "connection_check"),
                     ("surface_check_deviation", "connection_check_deviation"),
                     ("auto_surface_check", "auto_connection_check")):
        assert any(old in w and new in w for w in config.warnings), joined


def test_legacy_surface_confirmation_file_is_honored(workdir: Path):
    from smj.junction.fingerprint import record_confirmation

    cfg = _write_case(workdir, review_lines=["auto_connection_check: false"])
    build = workdir / "run" / "electrode_build"
    record_confirmation(build, who="王昊")
    (build / "connection_check.json").rename(build / "surface_check.json")

    assert main(["review", "-c", str(cfg)]) == 0
    payload = _payload(workdir)
    assert payload["connection_check"]["present"] is True
    assert payload["connection_check"]["valid"] is True
