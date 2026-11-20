"""T(E) 绘图：纵轴一律对数（2026-11-16 需求），且 T<=0 的点不得让 log 轴炸掉。

背景：``plot_transmission`` 早期是线性纵轴，T(E) 的深谷（反共振、~1e-3）
在线性轴上几乎看不见。改成 log 之后需要保证：

1. ``_log_floor`` 永远返回正数（T=0 / 负值 / 全零都不炸）；
2. ``plot_transmission`` 写出的图元数据带 ``log_scale: True`` 与正 ``floor``；
3. ``analyze_te(plot=True)`` 走通并把结果塞进 ``plot`` 键。
"""
from __future__ import annotations

import math

import pytest

pytest.importorskip("matplotlib")

from smj.analysis import T_LOG_FLOOR, _log_floor, analyze_te, plot_transmission  # noqa: E402

TE_TEMPLATE = """#bias = 0.000000 V efermi = -0.200000000000 H
{rows}
#end
"""


def _write_te(path, energies_ev, transmission):
    rows = []
    for e_ev, t in zip(energies_ev, transmission):
        rows.append("%.10f %.10f %.10f" % (-0.2 + e_ev / 27.211386, e_ev, t))
    path.write_text(TE_TEMPLATE.format(rows="\n".join(rows)), encoding="utf-8")
    return path


def test_log_floor_is_positive_for_degenerate_input():
    assert _log_floor([]) > 0.0
    assert _log_floor([0.0, 0.0]) >= T_LOG_FLOOR
    assert _log_floor([-1.0, -2.0]) >= T_LOG_FLOOR
    floor = _log_floor([0.0, -1.0, 4.0e-3])
    assert floor == pytest.approx(2.0e-3)


def test_plot_transmission_uses_log_axis(workdir):
    energies = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6]
    transmission = [3.0e-4, 1.0e-2, 0.5, 0.12, 0.4, 3.0e-3, 1.0e-3]
    te = _write_te(workdir / "TE.dat", energies, transmission)
    out = workdir / "T_E_plot.png"
    meta = plot_transmission(te, {"fermi_nearest": {"E_minus_Ef_eV": 0.0, "T": 0.12}}, out)

    assert meta["log_scale"] is True
    assert meta["points"] == len(energies)
    assert meta["floor"] > 0.0
    assert meta["floor"] <= min(transmission)
    assert out.is_file() and out.stat().st_size > 0
    assert meta["path"] == str(out)


def test_plot_transmission_survives_zero_transmission(workdir):
    te = _write_te(workdir / "TE.dat", [-0.2, 0.0, 0.2], [0.0, 0.0, 0.0])
    meta = plot_transmission(te, {}, workdir / "zero.png")
    assert meta["floor"] >= T_LOG_FLOOR
    assert math.isfinite(meta["floor"])


def test_analyze_te_plot_metadata(workdir):
    te = _write_te(workdir / "TE.dat", [-0.5, -0.25, 0.0, 0.25, 0.5],
                   [0.01, 0.6, 0.30, 0.55, 0.02])
    result = analyze_te(te, plot=True, independent=False)
    assert result["plot"]["log_scale"] is True
    assert result["plot"]["floor"] > 0.0
    assert (workdir / "T_E_plot.png").is_file()
