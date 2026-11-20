"""2026-11-20：开壳层 / 未登记族不得借用 BDT 的判据包络当闸门。

背景（用户提问触发的现场核实）：``criteria.py`` 里只有两条判据声明了
``applies_to = {'families': ('Au-BDT-Au',)}`` —— ``g_in_theory_band``
（派生包络 G/G0 ∈ [0.03, 0.40]）与 ``seebeck_sign``（"S 为正"铁律）。
旧实现 ``applies()`` 在 ``family is None`` 时 fail-open（"宁可多判一条"），
于是 ``runs/ttf-release``（family 未登记、多重度 2 的自由基结）因为
G=0.2337 G0 落进 BDT 包络、S=+10.72 µV/K 为正，被判成 ``passed / 支持 H0``
——两条数字碰巧落在别人的带里，不是这个体系的证据。

现在的口径：

1. 开壳层（multiplicity > 1）：判据族只针对闭壳层标定过，一律降为记录；
2. 族未知（family is None）：没有标定依据，同样降为记录；
3. 族已知且匹配（Au-BDT-Au 闭壳层）：仍然作闸门（回归）。
"""
from __future__ import annotations

from pathlib import Path

import pytest

from smj.junction.criteria import (applies, evaluate_run_criteria, is_open_shell,
                                   physics_verdict)

H_EV = 27.211386245988
GATES = ("g_in_theory_band", "seebeck_sign")


def _write_te(path: Path, *, n: int = 200, t: float = 0.25,
              efermi_ha: float = -0.188452434085) -> Path:
    """写一份 E−E_F ∈ [-2, 2] eV 的常数透射谱（与 ttf-release 的量级接近）。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["#bias = 0.000000 V efermi = %.12f H" % efermi_ha]
    for i in range(n):
        e_ev = -2.0 + 4.0 * i / (n - 1)
        lines.append("%.10f %.10f %.10e" % (efermi_ha + e_ev / H_EV, e_ev, t))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _by_name(checks: list[dict]) -> dict[str, dict]:
    return {c["name"]: c for c in checks}


@pytest.mark.parametrize("multiplicity,expected", [
    (None, False), (1, False), (2, True), (3, True), ("2", True), ("abc", False),
])
def test_is_open_shell_edges(multiplicity, expected: bool) -> None:
    assert is_open_shell(multiplicity) is expected


def test_applies_requires_known_family() -> None:
    criterion = {"applies_to": {"families": ("Au-BDT-Au",)}}
    assert applies(criterion, "Au-BDT-Au") is True
    # 族未知：不能用别的族的标定当闸门（旧实现返回 True -> 误判 passed）
    assert applies(criterion, None) is False
    # 开壳层：即使族匹配也降级
    assert applies(criterion, "Au-BDT-Au", open_shell=True) is False
    # 未声明适用域的判据仍普适
    assert applies({"name": "te_window_covers_efermi"}, None) is True


def test_open_shell_unknown_family_is_inconclusive(workdir: Path) -> None:
    """ttf-release 的场景：族未登记 + 多重度 2 -> 两条族闸门都不适用。"""
    te = _write_te(workdir / "TE.dat")
    checks = evaluate_run_criteria(te.parent, temperature=300.0, system=None,
                                   multiplicity=2)
    names = _by_name(checks)
    for gate in GATES:
        assert names[gate]["applicable"] is False
        assert "开壳层未校准" in names[gate]["detail"]
    verdict = physics_verdict(checks)
    assert verdict["physics_verdict"] == "inconclusive"
    assert "no_applicable_gate" in verdict["inconclusive_criteria"]
    assert "支持 H0" not in verdict["wording"]


def test_closed_shell_unknown_family_is_inconclusive(workdir: Path) -> None:
    """闭壳层但族未知：同样没有标定依据（旧实现 fail-open 会当闸门）。"""
    te = _write_te(workdir / "TE.dat")
    checks = evaluate_run_criteria(te.parent, temperature=300.0, system="NOPE",
                                   multiplicity=1)
    names = _by_name(checks)
    assert names["g_in_theory_band"]["applicable"] is False
    assert "体系族未知" in names["g_in_theory_band"]["detail"]
    assert physics_verdict(checks)["physics_verdict"] == "inconclusive"


def test_bdt_closed_shell_keeps_gates(workdir: Path) -> None:
    """回归：Au-BDT-Au 闭壳层仍按标定作闸门。"""
    te = _write_te(workdir / "TE.dat")
    checks = evaluate_run_criteria(te.parent, temperature=300.0, system="BDT",
                                   multiplicity=1)
    names = _by_name(checks)
    for gate in GATES:
        assert names[gate]["applicable"] is True
    assert "no_applicable_gate" not in physics_verdict(checks)["inconclusive_criteria"]


def test_bdt_open_shell_drops_gates(workdir: Path) -> None:
    """BDT 骨架上带自由基（多重度 2）也不能套用闭壳层标定。"""
    te = _write_te(workdir / "TE.dat")
    checks = evaluate_run_criteria(te.parent, temperature=300.0, system="BDT",
                                   multiplicity=2)
    for gate in GATES:
        assert _by_name(checks)[gate]["applicable"] is False


def test_report_manifest_does_not_hail_passed_for_open_shell(workdir: Path) -> None:
    """端到端：manifest 的 physics_verdict 必须是 inconclusive，不能写"支持 H0"。"""
    from smj.report import build_manifest

    run = workdir / "run"
    _write_te(run / "aitranss_te" / "TE.dat")
    manifest = build_manifest(
        run, include_checks=True,
        spec={"molecule": {"charge": 1, "multiplicity": 2, "family": None}},
    )
    verdict = manifest["verdict"]
    assert verdict["physics_verdict"] == "inconclusive"
    assert "no_applicable_gate" in verdict["inconclusive_criteria"]
    molecules = manifest["inputs"]["molecule"]
    assert molecules["multiplicity"] == 2
