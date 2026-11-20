"""t_ef_resonance_attribution: E_F 近旁有峰 -> 按族归因, 不再是"错误信号"。

2026-11-15 判据更新: 原 t_ef_no_peak_near_ef 只看"1.0 eV 内有没有峰", 峰在预期通道侧
也可能被判 warn。新口径把它改为**归因**判据: 无峰 -> pass; 峰落在族的预期通道侧且
|ΔE|∈[0.15,1.0] eV -> pass(物理信号已归因); 侧别不符/距 E_F 过近/通道未知 -> warn(需解释)。
未登记族按实测 Seebeck 符号自洽归因(T(E) 与 S 是同一组数据的不同矩)。
"""
from __future__ import annotations

import math
from pathlib import Path

from smj.junction.criteria import (
    EXPECTED_RESONANCE,
    evaluate_run_criteria,
    physics_verdict,
)

EV_PER_HA = 27.211386245988


def _te(path: Path, *, slope: float, peak_eV: float | None = None,
        peak_T: float = 0.6, width_eV: float = 0.06) -> None:
    """写合成 TE.dat: 线性背景 + 可选高斯峰。

    slope<0 -> T 随 E 下降 -> S>0(p 型/HOMO); slope>0 -> T 随 E 上升 -> S<0(n 型/LUMO)。
    """
    lines = ["#non-spin-polarized calculation",
             "#   E [Hartree]  E-EF [eV]  T(E) per spin"]
    for i in range(601):
        e = -1.5 + 0.005 * i
        t = 0.3 + slope * e
        if peak_eV is not None:
            t += peak_T * math.exp(-0.5 * ((e - peak_eV) / width_eV) ** 2)
        lines.append("  %.10f  %.8f  %.10e" % (e / EV_PER_HA, e, max(t, 0.0)))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _by_name(workdir: Path, **kw) -> dict:
    return {c["name"]: c for c in evaluate_run_criteria(workdir, **kw)}


def test_expected_resonance_table_carries_sources():
    for fam in ("Au-BDT-Au", "Au-(4,4'-bipyridine)-Au", "Au-(2,6-azulenedithiol)-Au",
                "Au-alkane-dithiol-Au"):
        entry = EXPECTED_RESONANCE[fam]
        assert entry["channel"] in ("HOMO", "LUMO", None), fam
        src = entry.get("source") or {}
        assert src.get("doi") and src.get("locator"), fam


def test_no_peak_passes_even_for_unknown_family(workdir: Path):
    _te(workdir / "TE.dat", slope=-0.02)
    check = _by_name(workdir)["t_ef_resonance_attribution"]
    assert check["passed"] is True


def test_bdt_homo_peak_below_ef_is_attributed(workdir: Path):
    _te(workdir / "TE.dat", slope=-0.05, peak_eV=-0.5)
    check = _by_name(workdir, system="BDT")["t_ef_resonance_attribution"]
    assert check["passed"] is True
    assert "已归因" in check["detail"]


def test_bpy_lumo_peak_above_ef_is_attributed(workdir: Path):
    _te(workdir / "TE.dat", slope=0.05, peak_eV=0.5)
    check = _by_name(workdir, system="BPY")["t_ef_resonance_attribution"]
    assert check["passed"] is True
    assert "已归因" in check["detail"]


def test_peak_too_close_to_ef_warns(workdir: Path):
    _te(workdir / "TE.dat", slope=0.05, peak_eV=0.05)
    check = _by_name(workdir, system="BPY")["t_ef_resonance_attribution"]
    assert check["passed"] is False
    assert "过近" in check["detail"]


def test_wrong_side_warns(workdir: Path):
    _te(workdir / "TE.dat", slope=-0.05, peak_eV=0.5)
    check = _by_name(workdir, system="BDT")["t_ef_resonance_attribution"]
    assert check["passed"] is False
    assert "侧别不符" in check["detail"]


def test_unregistered_family_falls_back_to_seebeck_sign(workdir: Path):
    # slope<0 -> S>0 -> HOMO; 峰在 E_F 下方 0.4 eV -> 自洽归因 pass
    _te(workdir / "TE.dat", slope=-0.05, peak_eV=-0.4)
    check = _by_name(workdir, system="2CN")["t_ef_resonance_attribution"]
    assert check["passed"] is True
    assert "自洽" in check["detail"]


def test_alkane_none_channel_any_peak_warns(workdir: Path):
    _te(workdir / "TE.dat", slope=-0.05, peak_eV=-0.4)
    check = _by_name(workdir, system="HDT")["t_ef_resonance_attribution"]
    assert check["passed"] is False


def test_resonance_warn_makes_verdict_inconclusive(workdir: Path):
    _te(workdir / "TE.dat", slope=0.05, peak_eV=0.05)
    out = physics_verdict(list(_by_name(workdir, system="BPY").values()))
    assert out["physics_verdict"] == "inconclusive"
    assert "t_ef_resonance_attribution" in out["inconclusive_criteria"]
