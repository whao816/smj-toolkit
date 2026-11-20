"""统一 T(E) 绘图风格（``smj.plotting``）：版式作用域、导出格式、叠加与元数据。

风格目的是让**每次 SMJ 跑完后的 T(E) 图长得一样**：四边框、刻度朝内、
E_F 虚线、log 纵轴、固定调色板、PNG+PDF 同源导出。这里只测"骨架"是否
真的固定下来，以及元数据是否够下游引用（style_version / floor / points）。
"""
from __future__ import annotations

import json
import math

import pytest

pytest.importorskip("matplotlib")

from smj.plotting import (  # noqa: E402
    FIG_SIZES,
    LOG_FLOOR,
    PALETTE,
    STYLE,
    STYLE_VERSION,
    log_floor,
    plot_te,
    plot_te_multi,
    read_curve,
    resolve_te,
    te_style,
)

TE_TEMPLATE = """#bias = 0.000000 V efermi = -0.200000000000 H
{rows}
#end
"""


def _write_te(path, energies_ev, transmission):
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for e_ev, t in zip(energies_ev, transmission):
        rows.append("%.10f %.10f %.10f" % (-0.2 + e_ev / 27.211386, e_ev, t))
    path.write_text(TE_TEMPLATE.format(rows="\n".join(rows)), encoding="utf-8")
    return path


def _write_te_spin(path, energies_ev, alpha, beta):
    """写自旋极化 TE.dat（四列：E(Ha) E−EF(eV) Tα Tβ）。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for e_ev, t_a, t_b in zip(energies_ev, alpha, beta):
        rows.append("%.10f %.10f %.10f %.10f"
                    % (-0.2 + e_ev / 27.211386, e_ev, t_a, t_b))
    body = "#spin-polarized calculation\n" + TE_TEMPLATE.format(rows="\n".join(rows))
    path.write_text(body, encoding="utf-8")
    return path


def _spectrum():
    energies = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6]
    transmission = [3.0e-4, 1.0e-2, 0.5, 0.12, 0.4, 3.0e-3, 1.0e-3]
    return energies, transmission


def test_style_is_scoped_and_does_not_leak(workdir):
    import matplotlib

    before = matplotlib.rcParams["axes.linewidth"]
    with te_style():
        assert matplotlib.rcParams["axes.linewidth"] == STYLE["axes.linewidth"]
        assert matplotlib.rcParams["xtick.direction"] == "in"
        assert matplotlib.rcParams["xtick.top"] is True
        assert matplotlib.rcParams["legend.frameon"] is False
    assert matplotlib.rcParams["axes.linewidth"] == before


def test_plot_te_writes_all_formats_with_metadata(workdir):
    energies, transmission = _spectrum()
    te = _write_te(workdir / "TE.dat", energies, transmission)

    meta = plot_te(te, workdir / "T_E.png", formats=("png", "pdf"))

    assert meta["points"] == len(energies)
    assert meta["log_scale"] is True
    assert meta["style_version"] == STYLE_VERSION
    assert meta["floor"] == pytest.approx(min(transmission) * 0.5)
    assert meta["paths"] == [str(workdir / "T_E.png"), str(workdir / "T_E.pdf")]
    assert (workdir / "T_E.png").is_file()
    assert (workdir / "T_E.pdf").is_file()


def test_plot_te_annotates_fermi_values(workdir):
    energies, transmission = _spectrum()
    te = _write_te(workdir / "TE.dat", energies, transmission)

    meta = plot_te(te, workdir / "T_E.png", result={
        "fermi_nearest": {"E_minus_Ef_eV": 0.0, "T": 0.12},
        "L0": 0.25,
    }, return_figure=True)
    _fig, ax = meta["figure"]
    texts = [item.get_text() for item in ax.texts]
    assert any("0.12" in text or "1.20" in text for text in texts)
    assert any("0.25" in text for text in texts)
    # 统一版式：不画网格
    assert not any(line.get_visible() for line in ax.get_xgridlines())


def test_plot_te_linear_scale_and_xlim(workdir):
    energies, transmission = _spectrum()
    te = _write_te(workdir / "TE.dat", energies, transmission)

    meta = plot_te(te, workdir / "T_E_lin.png", scale="linear",
                   xlim=(-0.5, 0.5), formats=("png",), return_figure=True)
    _fig, ax = meta["figure"]
    assert meta["log_scale"] is False
    assert ax.get_xlim() == pytest.approx((-0.5, 0.5))
    assert ax.get_yscale() == "linear"
    assert (workdir / "T_E_lin.png").is_file()


def test_plot_te_multi_uses_fixed_palette(workdir):
    energies, transmission = _spectrum()
    first = _write_te(workdir / "run_a" / "aitranss_te" / "TE.dat", energies, transmission)
    second = _write_te(workdir / "run_b" / "aitranss_te" / "TE.dat", energies,
                       [value * 0.5 for value in transmission])

    meta = plot_te_multi([first, second], workdir / "overlay.png",
                         labels=["run-a", "run-b"], figsize="double",
                         formats=("png",), return_figure=True)

    assert len(meta["curves"]) == 2
    assert [curve["label"] for curve in meta["curves"]] == ["run-a", "run-b"]
    _fig, ax = meta["figure"]
    colors = [line.get_color() for line in ax.get_lines()]
    assert colors[0] == PALETTE[0] and colors[1] == PALETTE[1]
    assert ax.get_legend() is not None
    assert (workdir / "overlay.png").is_file()


def test_log_floor_and_read_curve_round_trip(workdir):
    assert log_floor([0.0, 0.0]) == LOG_FLOOR
    assert log_floor([4.0e-3]) == pytest.approx(2.0e-3)

    energies, transmission = _spectrum()
    te = _write_te(workdir / "aitranss_te" / "TE.dat", energies, transmission)
    curve = read_curve(workdir)
    assert curve["e"] == pytest.approx(energies)
    assert curve["t"] == pytest.approx(transmission)
    assert resolve_te(workdir) == te


def test_figure_sizes_are_fixed():
    assert FIG_SIZES["single"][0] < FIG_SIZES["double"][0]
    assert all(math.isfinite(value) and value > 0 for size in FIG_SIZES.values() for value in size)


def test_cli_plot_renders_single_and_overlay(workdir):
    from smj.cli import main

    energies, transmission = _spectrum()
    _write_te(workdir / "run_a" / "aitranss_te" / "TE.dat", energies, transmission)
    _write_te(workdir / "run_b" / "aitranss_te" / "TE.dat", energies,
              [value * 0.5 for value in transmission])

    single = workdir / "fig" / "te.png"
    assert main(["plot", str(workdir / "run_a"), "-o", str(single),
                 "--formats", "png", "--dpi", "120"]) == 0
    assert single.is_file()

    overlay = workdir / "fig" / "cmp.png"
    assert main(["plot", str(workdir / "run_a"), str(workdir / "run_b"),
                 "--labels", "a,b", "-o", str(overlay), "--formats", "png",
                 "--dpi", "120", "--figsize", "double"]) == 0
    assert overlay.is_file()


def test_cli_plot_missing_te_returns_failed(workdir, capsys):
    from smj.cli import main

    rc = main(["plot", str(workdir / "nowhere"), "--formats", "png"])
    assert rc == 1
    assert "TE.dat" in capsys.readouterr().err


def test_junction_analyze_stage_writes_unified_plot(workdir):
    from types import SimpleNamespace

    from smj.stages.junction_analyze import JunctionAnalyzeStage

    energies, transmission = _spectrum()
    te = _write_te(workdir / "aitranss_te" / "TE.dat", energies, transmission)

    class Ctx:
        config = SimpleNamespace(transport=SimpleNamespace(temperature=300.0))

        def override_path(self, key):
            return None

        def upstream(self, stage, name, override=None):
            return override or te

        def log(self, message):
            pass

    prepared = JunctionAnalyzeStage().prepare(Ctx(), workdir)
    assert (workdir / "results.json").is_file()
    assert (workdir / "T_E.png").is_file()
    assert (workdir / "T_E.pdf").is_file()
    assert "T_E.png" in prepared.artifacts
    assert any("统一版式图" in note for note in prepared.notes)


# --------------------------------------------------------------------------- #
# 自旋自适应：单列退化为均值，双列出双 panel
# --------------------------------------------------------------------------- #
def test_read_curve_separates_spin_channels(workdir):
    energies, transmission = _spectrum()
    beta = [value * 0.8 for value in transmission]
    _write_te_spin(workdir / "aitranss_te" / "TE.dat", energies, transmission, beta)

    curve = read_curve(workdir)
    assert curve["channels"] == 2
    assert curve["t_alpha"] == pytest.approx(transmission)
    assert curve["t_beta"] == pytest.approx(beta)
    assert curve["t"] == pytest.approx([(a + b) / 2 for a, b in zip(transmission, beta)])


def test_plot_te_auto_spin_polarized_draws_two_panels(workdir):
    energies = _spectrum()[0]
    alpha = [0.10, 0.20, 0.50, 0.30, 0.45, 0.25, 0.15]
    beta = [0.08, 0.16, 0.40, 0.24, 0.36, 0.20, 0.12]
    te = _write_te_spin(workdir / "aitranss_te" / "TE.dat", energies, alpha, beta)

    meta = plot_te(te, workdir / "T_E_spin.png", formats=("png",), dpi=120,
                   panel="a", return_figure=True)

    assert meta["spin"] == "polarization"
    assert meta["channels"] == 2
    fig, axes = meta["figure"]
    assert len(axes) == 2
    # 单栏尺寸自动升到双栏
    assert fig.get_size_inches()[0] == pytest.approx(FIG_SIZES["double"][0])

    labels = [line.get_label() for line in axes[0].get_lines()]
    assert r"$T_\alpha$" in labels and r"$T_\beta$" in labels
    left_texts = [item.get_text() for item in axes[0].texts]
    assert "(a)" in left_texts
    assert any(r"T_\alpha(E_\mathrm{F})" in text for text in left_texts)

    right_texts = [item.get_text() for item in axes[1].texts]
    assert "(b)" in right_texts
    assert any("P(E_\\mathrm{F})" in text.replace(" ", "") for text in right_texts)
    assert axes[1].get_ylabel() == r"$P(E)$ (%)"
    assert (workdir / "T_E_spin.png").is_file()


def test_plot_te_spin_mean_forces_single_axes(workdir):
    energies, transmission = _spectrum()
    te = _write_te_spin(workdir / "aitranss_te" / "TE.dat", energies,
                        transmission, [value * 0.8 for value in transmission])

    meta = plot_te(te, workdir / "T_E_mean.png", spin="mean", formats=("png",),
                   return_figure=True)
    assert meta["spin"] == "mean"
    _fig, ax = meta["figure"]
    assert not isinstance(ax, (list, tuple))


def test_plot_te_single_channel_polarization_falls_back_to_mean(workdir):
    energies, transmission = _spectrum()
    te = _write_te(workdir / "aitranss_te" / "TE.dat", energies, transmission)

    meta = plot_te(te, workdir / "T_E_fallback.png", spin="polarization",
                   formats=("png",), return_figure=True)
    assert meta["channels"] == 1
    assert meta["spin"] == "mean"


def test_plot_te_rejects_unknown_spin_mode(workdir):
    energies, transmission = _spectrum()
    te = _write_te(workdir / "aitranss_te" / "TE.dat", energies, transmission)

    with pytest.raises(ValueError, match="spin"):
        plot_te(te, workdir / "T_E_bad.png", spin="bogus")


def test_cli_plot_spin_option(workdir, capsys):
    from smj.cli import main

    energies, transmission = _spectrum()
    _write_te_spin(workdir / "run_s" / "aitranss_te" / "TE.dat", energies,
                   transmission, [value * 0.5 for value in transmission])

    out = workdir / "fig" / "spin_mean.png"
    assert main(["plot", str(workdir / "run_s"), "-o", str(out),
                 "--formats", "png", "--dpi", "120", "--spin", "mean",
                 "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["spin"] == "mean"
    assert payload["channels"] == 2
    assert out.is_file()

    auto = workdir / "fig" / "spin_auto.png"
    assert main(["plot", str(workdir / "run_s"), "-o", str(auto),
                 "--formats", "png", "--dpi", "120", "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["spin"] == "polarization"
    assert auto.is_file()
