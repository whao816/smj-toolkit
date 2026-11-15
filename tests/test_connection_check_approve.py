"""F4 回归：连接核对（connection check；原 A2 人工确认）的独立登记入口，fail-closed。

背景（2026-11-11 补丁 F4；2026-11-15 更名并默认自动登记）：独立包此前没有任何
签署命令，``smj review`` 只读 ``a2_confirmation.json``。本组测试锁住登记入口的
几条底线：

* 没有建结目录 / 没有 metadata / 审计不过 / 出图失败 -> 拒绝登记，且**不写**
  ``connection_check.json``（旧名 ``a2_confirmation.json`` 仅用于读旧 run）；
* 审计不过时把 ``connection_check_audit.json`` 落盘，便于人看是哪几项不过；
* 几何不可解析时返回 refused 而不是抛异常；
* ``--no-figure`` 是显式逃生口，签出来的记录绑定当前 geometry.in 的 sha256；
* 二次登记幂等（already_signed），``force`` 才覆盖；
* 旧公开名 ``approve_a2`` 仍可用（别名）。
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import pytest

from smj.a2 import EXIT_REFUSED, approve_a2, approve_connection_check
from smj.cli import main
from smj.errors import SmjError

# 表面三角的三个 (x, y) 偏移（约 2.88 A 边长，面积 ~3.6 A^2 > 0.5 A^2 下限）
TRI = ((0.0, 0.0), (2.88, 0.0), (1.44, 2.494))


def _atoms() -> list[tuple[str, float, float, float]]:
    """1-6 C 环 / 7-8 S 锚 / 9-10 Au tip / 11-16 左电极两层 / 17-22 右电极两层。

    左外层 = 14,15,16（z=-4.0），右外层 = 20,21,22（z=+4.0）；两侧最外层各恰好
    3 个 Au，故 det-v1 重算必然选中同一组（audit 的 `deterministic_indices` 硬检查）。
    """
    atoms: list[tuple[str, float, float, float]] = []
    for k in range(6):
        ang = math.radians(60.0 * k)
        atoms.append(("C", 1.39 * math.cos(ang), 1.39 * math.sin(ang), 0.0))
    atoms.append(("S", 0.0, 0.0, -2.0))      # 7  左 linker（分子锚原子）
    atoms.append(("S", 0.0, 0.0, 2.0))       # 8  右 linker
    atoms.append(("Au", 0.0, 0.0, -2.6))     # 9  左 tip
    atoms.append(("Au", 0.0, 0.0, 2.6))      # 10 右 tip
    for z in (-3.4, -4.0):
        for dx, dy in TRI:
            atoms.append(("Au", dx, dy, z))
    for z in (3.4, 4.0):
        for dx, dy in TRI:
            atoms.append(("Au", dx, dy, z))
    return atoms


def _geometry_text(atoms) -> str:
    return "".join("atom %.8f %.8f %.8f %s\n" % (x, y, z, e) for e, x, y, z in atoms)


def _metadata(**overrides) -> dict:
    payload = {
        "left_tip_index": 9,
        "right_tip_index": 10,
        "left_linker_index": 7,
        "right_linker_index": 8,
        "left_surface_indices": [14, 15, 16],
        "right_surface_indices": [20, 21, 22],
        "surface_rule": "det-v1",
    }
    payload.update(overrides)
    return payload


def _case(workdir: Path, *, geometry_text: str | None = None,
          metadata: dict | None = None, with_metadata: bool = True,
          build_name: str = "electrode_build") -> tuple[Path, Path]:
    run_dir = workdir / "run"
    build = run_dir / build_name
    build.mkdir(parents=True)
    (build / "geometry.in").write_text(
        geometry_text if geometry_text is not None else _geometry_text(_atoms()),
        encoding="utf-8")
    if with_metadata:
        (build / "electrode_metadata.json").write_text(
            json.dumps(metadata if metadata is not None else _metadata(), indent=2),
            encoding="utf-8")
    return run_dir, build


def _fingerprint(build: Path) -> str:
    return hashlib.sha256((build / "geometry.in").read_bytes()).hexdigest()


# --------------------------------------------------------------------------- #
# 拒绝路径
# --------------------------------------------------------------------------- #
def test_refuses_without_build_dir(workdir: Path):
    result = approve_connection_check(workdir / "run", who="tester")
    assert result.status == "refused"
    assert result.code == "no_build_dir"
    assert result.exit_code == EXIT_REFUSED
    assert not result.ok


def test_refuses_without_metadata(workdir: Path):
    run_dir, build = _case(workdir, with_metadata=False)
    result = approve_connection_check(run_dir, who="tester")
    assert result.code == "no_metadata"
    assert result.exit_code == EXIT_REFUSED
    assert not (build / "connection_check.json").exists()


def test_refuses_when_audit_not_ok(workdir: Path):
    # 内层 11,12,13 不是最外层 -> same_surface_layer 等检查必不过
    run_dir, build = _case(workdir, metadata=_metadata(
        left_surface_indices=[11, 12, 13], right_surface_indices=[11, 12, 13]))
    result = approve_connection_check(run_dir, who="tester", no_figure=True)
    assert result.status == "refused"
    assert result.code == "audit_not_ok"
    assert result.exit_code == EXIT_REFUSED
    assert result.checks_failed, "拒绝时必须列出未通过的检查项"
    assert (build / "connection_check_audit.json").is_file()
    assert not (build / "connection_check.json").exists()


def test_refuses_when_figure_fails(workdir: Path, monkeypatch):
    from smj._vendor.junction import a2_viz

    def boom(*args, **kwargs):
        raise RuntimeError("出图需要 matplotlib(可选依赖)")

    monkeypatch.setattr(a2_viz, "render_figure", boom)
    run_dir, build = _case(workdir)
    result = approve_connection_check(run_dir, who="tester")          # 默认要求出图
    assert result.code == "figure_failed"
    assert result.exit_code == EXIT_REFUSED
    assert not (build / "connection_check.json").exists()


def test_unparsable_geometry_is_refused_not_crash(workdir: Path):
    run_dir, build = _case(workdir, geometry_text="this is not an aims geometry\n")
    result = approve_connection_check(run_dir, who="tester", no_figure=True)
    assert result.status == "refused"
    assert result.code == "audit_failed"
    assert result.exit_code == EXIT_REFUSED
    assert not (build / "connection_check.json").exists()


def test_empty_who_raises(workdir: Path):
    run_dir, _build = _case(workdir)
    with pytest.raises(SmjError, match="签署人"):
        approve_connection_check(run_dir, who="   ")


# --------------------------------------------------------------------------- #
# 签署路径
# --------------------------------------------------------------------------- #
def test_signs_with_no_figure(workdir: Path):
    run_dir, build = _case(workdir)
    result = approve_connection_check(run_dir, who="王昊", no_figure=True)
    assert result.ok, result.reason
    assert result.status == "signed"
    assert result.exit_code == 0
    confirm = json.loads((build / "connection_check.json").read_text(encoding="utf-8"))
    assert confirm["who"] == "王昊"
    assert confirm["fingerprint"] == _fingerprint(build)
    assert (build / "connection_check_audit.json").is_file()


def test_second_approval_is_idempotent(workdir: Path):
    run_dir, build = _case(workdir)
    first = approve_connection_check(run_dir, who="tester", no_figure=True)
    assert first.ok
    recorded = json.loads((build / "connection_check.json").read_text(encoding="utf-8"))

    again = approve_connection_check(run_dir, who="other", no_figure=True)
    assert again.status == "already_signed"
    assert again.ok and again.exit_code == 0
    # 未加 force：记录必须原样不动
    assert json.loads((build / "connection_check.json").read_text(
        encoding="utf-8"))["who"] == recorded["who"]


def test_force_resigns(workdir: Path):
    run_dir, build = _case(workdir)
    assert approve_connection_check(run_dir, who="tester", no_figure=True).ok
    again = approve_connection_check(run_dir, who="王昊", no_figure=True, force=True)
    assert again.status == "signed"
    assert json.loads((build / "connection_check.json").read_text(
        encoding="utf-8"))["who"] == "王昊"


def test_signature_follows_geometry(workdir: Path):
    """改了几何 -> 旧签名失效（指纹绑字节），重签才回到有效。"""
    from smj._vendor.junction.fingerprint import check_confirmation

    run_dir, build = _case(workdir)
    assert approve_connection_check(run_dir, who="tester", no_figure=True).ok
    assert check_confirmation(build)["valid"] is True

    atoms = _atoms()
    atoms[0] = ("C", 1.5, 0.0, 0.0)
    (build / "geometry.in").write_text(_geometry_text(atoms), encoding="utf-8")
    assert check_confirmation(build)["valid"] is False

    assert approve_connection_check(run_dir, who="tester", no_figure=True, force=True).ok
    assert check_confirmation(build)["valid"] is True


# --------------------------------------------------------------------------- #
# CLI 接线
# --------------------------------------------------------------------------- #
def _config_file(workdir: Path, run_dir: Path) -> Path:
    """run_dir = <run.dir>/<run.name>（yaml 的 dir 是 run 根），故 dir=workdir、name=run。"""
    xyz = workdir / "molecule.xyz"
    xyz.write_text("2\nprobe\nS 0.0 0.0 -1.0\nS 0.0 0.0 1.0\n", encoding="utf-8")
    cfg = workdir / "smj.yaml"
    cfg.write_text(
        "run:\n  name: %s\n  dir: %s\n  mode: dry-run\n"
        "molecule:\n  xyz: %s\n"
        % (run_dir.name, run_dir.parent.resolve().as_posix(), xyz.resolve().as_posix()),
        encoding="utf-8")
    return cfg


def test_cli_requires_who(workdir: Path):
    run_dir, _build = _case(workdir)
    cfg = _config_file(workdir, run_dir)
    assert main(["review", "--approve", "-c", str(cfg)]) == 2


def test_cli_refusal_exit_code_is_2(workdir: Path):
    run_dir, build = _case(workdir, metadata=_metadata(
        left_surface_indices=[11, 12, 13]))
    cfg = _config_file(workdir, run_dir)
    code = main(["review", "--approve", "--who", "tester", "--no-figure",
                 "--json", "-c", str(cfg)])
    assert code == 2
    assert not (build / "connection_check.json").exists()


def test_cli_signs(workdir: Path):
    run_dir, build = _case(workdir)
    cfg = _config_file(workdir, run_dir)
    code = main(["review", "--approve", "--who", "tester", "--no-figure",
                 "-c", str(cfg)])
    assert (build / "connection_check.json").is_file(), "CLI 必须真的把确认写进 build 目录"
    # 签署后继续走 review；这里只断言签署结果已落盘，不约束 review 自身的结论
    assert code in (0, 1)
