"""TE 谱行序无关性（复核报告 R1 / 原报告 P1-1）。

TE.dat 的**行序不是物理量**，但旧实现按行序梯形积分、且只对偶数阶矩取 ``abs``，
于是倒序文件会让唯一带符号的 L1 反号 —— Seebeck 静默变号而 G/κe/ZT 不变。
本文件锁四件事：

* 倒序、乱序谱与升序谱给出同一组 ``L0/L1/L2/S/κe/ZT``；
* ``analyze_te`` 显式报告 ``input_energy_order`` / ``reordered`` / ``warnings``；
* 修好之后的 Seebeck 符号与"排序后再积分"的参考值一致（不是被 abs 掩盖的巧合）；
* 独立复核实现 ``independent_check.l0_l1_l2`` 同样顺序无关。
"""
from __future__ import annotations

import math
import random
from pathlib import Path

import pytest

from smj._vendor.junction import independent_check, landauer

EF_HARTREE = -0.188452434085
S_KEYS = ("L0", "L1", "L2", "S_uV_per_K", "k_e_pW_per_K", "ZT_electronic")


def _spectrum(n: int = 201, center: float = 0.30, width: float = 0.15,
              baseline: float = 0.05) -> list[tuple[float, float, float]]:
    """(E_Ha, E-E_F eV, T)：一个偏离 E_F 的洛伦兹状峰 + 本底。"""
    rows = []
    for i in range(n):
        e = -1.0 + 2.0 * i / (n - 1)
        t = 0.6 * math.exp(-((e - center) / width) ** 2) + baseline
        rows.append((EF_HARTREE + e, e, t))
    return rows


def _write(path: Path, rows: list[tuple[float, float, float]]) -> Path:
    lines = ["# bias = 0.000000 V efermi = %.12f H" % EF_HARTREE]
    lines += ["%.10f %.10f %.10e" % row for row in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _reference_integral(rows, temperature: float = 300.0) -> dict:
    """顺序无关参考积分：先按 E-E_F 升序排序，再独立梯形求积。"""
    ordered = sorted(rows, key=lambda r: r[1])
    kb = landauer.KB_EV_PER_K
    e = [r[1] for r in ordered]
    t = [max(r[2], 0.0) for r in ordered]
    w = [1.0 / (4.0 * math.cosh(x / (2.0 * kb * temperature)) ** 2) for x in e]
    l0 = landauer.trapz([ti * wi for ti, wi in zip(t, w)], e)
    l1 = landauer.trapz([ei * ti * wi for ei, ti, wi in zip(e, t, w)], e)
    l2 = landauer.trapz([ei * ei * ti * wi for ei, ti, wi in zip(e, t, w)], e)
    return {"L0": abs(l0), "L1": l1, "L2": abs(l2),
            "S_uV_per_K": 1e6 * (-l1) / (abs(l0) * temperature)}


def test_descending_te_matches_ascending(tmp_path: Path) -> None:
    rows = _spectrum()
    asc = landauer.analyze_te(_write(tmp_path / "asc.dat", rows))
    desc = landauer.analyze_te(_write(tmp_path / "desc.dat", list(reversed(rows))))

    for key in S_KEYS:
        assert desc[key] == pytest.approx(asc[key], rel=1e-9, abs=1e-12), key
    assert asc["input_energy_order"] == "ascending"
    assert asc["reordered"] is False
    assert asc["warnings"] == []
    assert desc["input_energy_order"] == "descending"
    assert desc["reordered"] is True
    assert any("升序" in w for w in desc["warnings"])


def test_unsorted_te_is_sorted_before_integration(tmp_path: Path) -> None:
    rows = _spectrum()
    shuffled = list(rows)
    random.Random(20261120).shuffle(shuffled)
    base = landauer.analyze_te(_write(tmp_path / "asc.dat", rows))
    messy = landauer.analyze_te(_write(tmp_path / "shuffled.dat", shuffled))

    for key in S_KEYS:
        assert messy[key] == pytest.approx(base[key], rel=1e-9, abs=1e-12), key
    assert messy["input_energy_order"] == "unsorted"
    assert messy["reordered"] is True


def test_seebeck_matches_order_invariant_reference(tmp_path: Path) -> None:
    """符号由谱本身决定：峰在 E_F 上方 → S<0；倒序文件不得翻转它。"""
    rows = _spectrum(center=0.30)
    reference = _reference_integral(rows)
    assert reference["S_uV_per_K"] < 0.0
    asc = landauer.analyze_te(_write(tmp_path / "asc.dat", rows))
    desc = landauer.analyze_te(_write(tmp_path / "desc.dat", list(reversed(rows))))
    assert asc["S_uV_per_K"] == pytest.approx(reference["S_uV_per_K"], rel=1e-10)
    assert desc["S_uV_per_K"] == pytest.approx(reference["S_uV_per_K"], rel=1e-10)


def test_independent_check_is_order_invariant() -> None:
    rows = _spectrum()
    asc = independent_check.l0_l1_l2(list(rows))
    desc = independent_check.l0_l1_l2(list(reversed(rows)))
    for key in ("L0", "L1", "L2"):
        assert desc[key] == pytest.approx(asc[key], rel=1e-9, abs=1e-15), key


def test_analysis_surfaces_reorder_warning(tmp_path: Path) -> None:
    """wrapper 不得把 vendored 的排序告警清空（旧实现 ``result["warnings"]=[]``）。"""
    from smj.analysis import analyze_te

    rows = _spectrum()
    result = analyze_te(_write(tmp_path / "TE.dat", list(reversed(rows))),
                        independent=False)
    assert result["input_energy_order"] == "descending"
    assert result["reordered"] is True
    assert any("升序" in w for w in result["warnings"])
