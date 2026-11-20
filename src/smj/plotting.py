"""统一的 T(E) 绘图风格（发表级版式）。

风格依据
--------
* Yi et al., *Chin. Phys. B* **25**, 128503 (2016), Fig. 6：四边封闭边框、
  刻度朝内、E_F 处竖直虚线、log 透射纵轴、panel 内 (a)(b)(c) 角标、
  图例在 panel 内且无边框。
* Guan et al., *Nanotechnology* **34**, 095401 (2023), Figs. 2-3：四边框、
  刻度朝内、panel 内小图例、E_F 参考虚线、紧凑的单栏/双栏尺寸。

固定下来的是**排版骨架**（字体、字号、线宽、边框、刻度、调色板、导出格式），
物理口径不变：

* x 轴恒为 ``E - E_F``（eV）；E_F 位置恒画竖直虚线；
* y 轴恒为 ``T(E)``，单位 ``G0 = 2e^2/h``（per spin，与 aitranss 口径一致）；
* log 模式下低于 ``LOG_FLOOR`` 的点按 floor 截断（不静默丢点）。

matplotlib 属可选依赖（``pip install -e ".[plot]"``）。本模块**顶层不导入
matplotlib**：所有绘图函数内部才 import，因此 ``smj analyze`` 在没有
matplotlib 的环境里仍可优雅跳过出图。
"""

from __future__ import annotations

import math
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

#: 风格版本号；写进每次出图的元数据，便于日后追溯版式变更。
STYLE_VERSION = "2026.11"

#: log 轴的下限（与 :mod:`smj.analysis` 的 ``T_LOG_FLOOR`` 同源）。
LOG_FLOOR = 1e-12

#: 单条曲线的主色（深蓝），与文献中常见的黑/蓝主曲线一致。
DEFAULT_COLOR = "#1f4e79"

#: 多曲线叠加时的固定配色顺序（不随调用次数漂移）。
PALETTE: tuple[str, ...] = (
    "#1f4e79",  # 深蓝
    "#c0392b",  # 砖红
    "#1e8449",  # 深绿
    "#7d3c98",  # 紫
    "#b9770e",  # 赭黄
    "#117a8b",  # 青
)

#: 成图尺寸（英寸）：单栏 / 双栏 / 宽版。
FIG_SIZES: dict[str, tuple[float, float]] = {
    "single": (3.35, 2.55),
    "double": (7.10, 3.00),
    "wide": (6.20, 2.80),
}

#: E_F 参考线颜色。
EF_COLOR = "#7f7f7f"

#: 固定 rcParams（不修改全局 rcParams，只在 :func:`te_style` 里生效）。
STYLE: dict[str, Any] = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8.0,
    "mathtext.fontset": "dejavusans",
    "axes.linewidth": 0.8,
    "axes.labelsize": 9.0,
    "axes.labelpad": 3.0,
    "axes.titlesize": 9.0,
    "axes.titlepad": 6.0,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "xtick.minor.size": 1.8,
    "ytick.minor.size": 1.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.minor.width": 0.6,
    "ytick.minor.width": 0.6,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.labelsize": 8.0,
    "ytick.labelsize": 8.0,
    "lines.linewidth": 1.1,
    "lines.solid_capstyle": "round",
    "legend.frameon": False,
    "legend.fontsize": 7.0,
    "legend.handlelength": 1.6,
    "legend.handletextpad": 0.5,
    "legend.labelspacing": 0.25,
    "legend.borderaxespad": 0.4,
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}

_SUFFIXES = {"png", "pdf", "svg", "eps", "jpg", "jpeg", "tif", "tiff"}


def _pyplot():
    """拿 pyplot，并尽量落到无界面 Agg 后端（CLI 只往文件里写图）。"""
    import matplotlib

    try:
        matplotlib.use("Agg", force=True)
    except Exception:  # pragma: no cover - 交互式后端已在用时不强切
        pass
    import matplotlib.pyplot as plt

    return plt


# --------------------------------------------------------------------------- #
# 风格与工具
# --------------------------------------------------------------------------- #
@contextmanager
def te_style():
    """在固定版式下作图；退出后恢复调用前的 rcParams（不改全局状态）。"""
    import matplotlib

    with matplotlib.rc_context(STYLE):
        yield


def log_floor(values: Iterable[float]) -> float:
    """log 轴下限：最小正值的 1/2，但不低于 :data:`LOG_FLOOR`。"""
    finite = [float(v) for v in values if float(v) > 0.0 and math.isfinite(float(v))]
    if not finite:
        return LOG_FLOOR
    return max(min(finite) * 0.5, LOG_FLOOR)


def resolve_te(path: str | Path) -> Path:
    """定位 TE.dat；接受文件、run 根目录或阶段目录（复用 analysis 的规则）。"""
    from .analysis import resolve_te as _resolve

    return Path(_resolve(path))


def read_curve(path: str | Path) -> dict[str, Any]:
    """读 TE.dat 为 ``{"e", "t", "t_alpha", "t_beta", "channels", "source"}``。

    ``t`` 恒为双通道均值 (Tα+Tβ)/2（与 vendored ``read_te`` 口径一致，保证旧
    调用不变）；``t_alpha``/``t_beta`` 是分通道曲线，单列文件里两者相等且
    ``channels == 1``。支持 D 指数、跳过注释。
    """
    from .junction.landauer import energy_order, sort_by_energy
    from .landauer_spin import read_te_channels

    te = resolve_te(path)
    rows, channels = read_te_channels(te)
    order = energy_order(rows)
    if order != "ascending":
        # 2026-11-20 复审 P1-3: 曲线与峰检测/积分共用同一能量序, 乱序文件
        # 不再把"文件相邻行"画成"能量相邻点"。
        rows = sort_by_energy(rows)
    alpha = [max(float(r[2]), 0.0) for r in rows]
    beta = [max(float(r[3]), 0.0) for r in rows]
    return {
        "source": str(te),
        "e": [float(r[1]) for r in rows],
        "t": [(a + b) / 2.0 for a, b in zip(alpha, beta)],
        "t_alpha": alpha,
        "t_beta": beta,
        "channels": channels,
        "energy_order": order,
    }


def _fmt_num(value: float | None) -> str:
    """把数值写成 3 位有效数字；极小/极大用 mathtext 科学计数法。"""
    if value is None or not math.isfinite(float(value)):
        return "n/a"
    x = float(value)
    if x == 0.0:
        return "0"
    exp = int(math.floor(math.log10(abs(x))))
    if -2 <= exp <= 3:
        return f"{x:.3g}"
    mantissa = x / 10 ** exp
    return rf"{mantissa:.2f}\times10^{{{exp}}}"


def _summary(te: Path, temperature: float, result: Mapping[str, Any] | None) -> Mapping[str, Any]:
    """优先用调用方给的 Landauer 结果，否则现算一份（不写文件）。"""
    if isinstance(result, Mapping) and "fermi_nearest" in result:
        return result
    from .junction.landauer import analyze_te as _analyze

    return _analyze(te, temperature=float(temperature))


def _prepare_axes(ax: Any, scale: str) -> None:
    from matplotlib.ticker import (AutoMinorLocator, LogFormatterMathtext,
                                   LogLocator, MaxNLocator, NullFormatter)

    for spine in ax.spines.values():
        spine.set_linewidth(STYLE["axes.linewidth"])
        spine.set_visible(True)
    ax.tick_params(which="both", direction="in", top=True, right=True,
                   width=STYLE["xtick.major.width"])
    ax.tick_params(which="major", size=STYLE["xtick.major.size"])
    ax.tick_params(which="minor", size=STYLE["xtick.minor.size"])
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    if scale == "log":
        ax.set_yscale("log")
        ax.yaxis.set_major_locator(LogLocator(base=10.0))
        ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10.0))
        ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=tuple(range(2, 10))))
        ax.yaxis.set_minor_formatter(NullFormatter())
    else:
        ax.set_yscale("linear")
        ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))


def _fermi_line(ax: Any, label: bool = True) -> None:
    ax.axvline(0.0, color=EF_COLOR, linewidth=0.8, linestyle=(0, (4, 2)),
               zorder=1)
    if label:
        ax.annotate(r"$E_\mathrm{F}$", xy=(0.0, 0.945),
                    xycoords=("data", "axes fraction"), xytext=(3, 0),
                    textcoords="offset points", ha="left", va="top",
                    fontsize=7.0, color=EF_COLOR)


def _panel_label(ax: Any, panel: str | None) -> None:
    if not panel:
        return
    ax.text(0.025, 0.965, "(%s)" % panel, transform=ax.transAxes,
            ha="left", va="top", fontsize=9.0, fontweight="bold")


def _annotate_stats(ax: Any, summary: Mapping[str, Any], *, gap: tuple[float, float] | None) -> None:
    nearest = summary.get("fermi_nearest") or {}
    lines = []
    if nearest.get("T") is not None:
        lines.append(r"$T(E_\mathrm{F}) = %s$" % _fmt_num(nearest.get("T")))
    l0 = summary.get("L0")
    if l0 is not None:
        lines.append(r"$G/G_0 = %s$" % _fmt_num(l0))
    if gap:
        lines.append(r"$E_\mathrm{gap} = %.2f\ \mathrm{eV}$" % (gap[1] - gap[0]))
    if not lines:
        return
    ax.text(0.03, 0.04, "\n".join(lines), transform=ax.transAxes,
            ha="left", va="bottom", fontsize=7.0, linespacing=1.35)


def _log_top(max_value: float) -> float:
    """log 纵轴上限：数据贴近整十次幂时留 50% 余量，避免曲线顶到边框。"""
    decade = 10.0 ** math.ceil(math.log10(max_value))
    return max_value * 1.5 if max_value >= 0.9 * decade else decade


def _log_ylim(values: Sequence[float], floor: float) -> tuple[float, float]:
    return floor * 0.6, max(_log_top(max(max(values), floor)), floor * 10.0)


def _save(fig: Any, out: str | Path, formats: Sequence[str] | None, dpi: float) -> list[Path]:
    base = Path(out)
    if formats is None:
        formats = (base.suffix.lstrip(".").lower() or "png",)
    formats = tuple(str(f).strip().lower().lstrip(".") for f in formats if str(f).strip())
    bad = [f for f in formats if f not in _SUFFIXES]
    if bad:
        raise ValueError("不支持的图片格式: %s（可选 %s）" % (", ".join(bad), ", ".join(sorted(_SUFFIXES))))
    stem = base.with_suffix("") if base.suffix.lower().lstrip(".") in _SUFFIXES else base
    paths = []
    for ext in formats:
        target = stem.with_suffix("." + ext)
        target.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(target, dpi=dpi)
        paths.append(target)
    return paths


def _default_label(source: str | Path) -> str:
    p = Path(source)
    te = resolve_te(p) if not p.is_file() else p
    return te.parent.parent.name if te.parent.name == "aitranss_te" else te.stem


# --------------------------------------------------------------------------- #
# 自旋通道（α/β）
# --------------------------------------------------------------------------- #
SPIN_MODES: tuple[str, ...] = ("auto", "mean", "channels", "polarization")

#: α、β 两条通道的固定配色（深蓝 / 砖红，与 :data:`PALETTE` 同一套）。
ALPHA_COLOR = PALETTE[0]
BETA_COLOR = PALETTE[1]
POLARIZATION_COLOR = PALETTE[2]


def _spin_mode(requested: str, channels: int) -> str:
    """把 ``--spin`` 解析成实际画法；单通道数据上 channels/polarization 退化为 mean。"""
    if requested not in SPIN_MODES:
        raise ValueError("spin 只能是 %s，收到 %r" % ("/".join(SPIN_MODES), requested))
    if requested == "auto":
        requested = "polarization" if channels >= 2 else "mean"
    if requested in ("channels", "polarization") and channels < 2:
        return "mean"
    return requested


def _spin_polarization(rows: Mapping[str, Any]) -> list[float]:
    """P(E) = (Tα − Tβ) / (Tα + Tβ)；两通道都为 0 时记 0。"""
    return [(a - b) / (a + b) if (a + b) > 0.0 else 0.0
            for a, b in zip(rows["t_alpha"], rows["t_beta"])]


def _nearest_index(e: Sequence[float]) -> int:
    return min(range(len(e)), key=lambda i: abs(float(e[i])))


def _annotate_spin_channels(ax: Any, rows: Mapping[str, Any],
                            summary: Mapping[str, Any] | None) -> None:
    i = _nearest_index(rows["e"])
    lines = [r"$T_\alpha(E_\mathrm{F}) = %s$" % _fmt_num(rows["t_alpha"][i]),
             r"$T_\beta(E_\mathrm{F}) = %s$" % _fmt_num(rows["t_beta"][i])]
    l0 = summary.get("L0") if summary else None
    if l0 is not None:
        lines.append(r"$G/G_0 = %s$" % _fmt_num(l0))
    ax.text(0.03, 0.04, "\n".join(lines), transform=ax.transAxes,
            ha="left", va="bottom", fontsize=7.0, linespacing=1.35)


def _annotate_polarization(ax: Any, rows: Mapping[str, Any]) -> None:
    p = _spin_polarization(rows)
    i = _nearest_index(rows["e"])
    lines = [r"$P(E_\mathrm{F}) = %.2f\,\%%$" % (100.0 * p[i]),
             r"$\max|P| = %.2f\,\%%$" % (100.0 * max(abs(v) for v in p))]
    ax.text(0.03, 0.04, "\n".join(lines), transform=ax.transAxes,
            ha="left", va="bottom", fontsize=7.0, linespacing=1.35)


def _next_panel(panel: str | None) -> str | None:
    """panel 角标递增：'a' -> 'b'，其余原样返回 None。"""
    if not panel or len(panel) != 1 or not panel.isalpha():
        return None
    return chr(ord(panel) + 1)


def _spin_figsize(figsize: str | Sequence[float]) -> tuple[float, float]:
    """双 panel 成图尺寸：默认从单栏自适应到双栏。"""
    if isinstance(figsize, str):
        return FIG_SIZES["double"] if figsize == "single" else FIG_SIZES[figsize]
    return tuple(float(v) for v in figsize)


# --------------------------------------------------------------------------- #
# 单曲线
# --------------------------------------------------------------------------- #
def plot_te(
    source: str | Path,
    out: str | Path | None = None,
    *,
    scale: str = "log",
    figsize: str | Sequence[float] = "single",
    formats: Sequence[str] | None = None,
    dpi: float = 600,
    xlim: Sequence[float] | None = None,
    ylim: Sequence[float] | None = None,
    title: str | None = None,
    panel: str | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    annotate: bool = True,
    gap: tuple[float, float] | None = None,
    color: str | None = None,
    label: str | None = None,
    result: Mapping[str, Any] | None = None,
    temperature: float = 300.0,
    spin: str = "auto",
    return_figure: bool = False,
) -> dict[str, Any]:
    """画 T(E) 曲线（统一风格）；自旋极化数据自适应成双 panel。

    Parameters
    ----------
    source:
        TE.dat 文件、run 根目录，或含 ``aitranss_te/TE.dat`` 的阶段目录。
    out:
        输出路径；后缀决定单格式。给了 ``formats`` 时按格式列表**同时**导出。
    scale:
        ``"log"``（文献默认）或 ``"linear"``。
    figsize:
        ``"single"`` / ``"double"`` / ``"wide"``，或 ``(宽, 高)`` 英寸。
        自旋双 panel 在 ``"single"`` 下自动改用 ``"double"`` 宽度。
    annotate:
        是否在 panel 内标注 ``T(E_F)`` 与 ``G/G0``（由 ``result`` 或现算得到）。
    gap:
        ``(E_HOMO, E_LUMO)``（eV），给出时以浅灰底色标出能隙（仅均值面板）。
    result:
        已算好的 Landauer 结果（``smj.analysis`` 返回的 dict）；为空时自动算。
    spin:
        ``"auto"``（默认）：TE.dat 有 α/β 两列时画双 panel——(a) 两通道
        T(E)，(b) 自旋极化率 ``P(E)=(Tα−Tβ)/(Tα+Tβ)``；单列时与旧行为一致。
        ``"mean"`` 只画双通道均值；``"channels"`` 单 panel 叠画 α/β；
        ``"polarization"`` 强制双 panel（单通道数据退化为 mean）。
    """
    if scale not in ("log", "linear"):
        raise ValueError("scale 只能是 'log' 或 'linear'，收到 %r" % (scale,))
    rows = read_curve(source)
    e, t = rows["e"], rows["t"]
    if len(e) < 2:
        raise ValueError("TE.dat 数据点不足（%d 个）：%s" % (len(e), rows["source"]))

    channels = int(rows.get("channels", 1))
    mode = _spin_mode(spin, channels)
    floor = log_floor(t if mode == "mean" else list(rows["t_alpha"]) + list(rows["t_beta"]))
    size = FIG_SIZES[figsize] if isinstance(figsize, str) else tuple(float(v) for v in figsize)

    plt = _pyplot()
    out_path = Path(out) if out is not None else (
        Path(rows["source"]).parent / ("T_E.%s" % (formats[0] if formats else "png")))
    summary = _summary(Path(rows["source"]), temperature, result) if annotate else None

    axes: Any
    with te_style():
        if mode == "mean":
            fig, ax = plt.subplots(figsize=size)
            try:
                if gap:
                    ax.axvspan(float(gap[0]), float(gap[1]), color="#ececec", zorder=0, lw=0.0)
                y_plot = [max(v, floor) for v in t] if scale == "log" else t
                ax.plot(e, y_plot, color=color or DEFAULT_COLOR, lw=STYLE["lines.linewidth"],
                        label=label or None, zorder=3)
                _prepare_axes(ax, scale)
                _fermi_line(ax)

                if xlim:
                    ax.set_xlim(float(xlim[0]), float(xlim[1]))
                else:
                    ax.set_xlim(min(e), max(e))
                if ylim:
                    ax.set_ylim(float(ylim[0]), float(ylim[1]))
                elif scale == "log":
                    ax.set_ylim(*_log_ylim(t, floor))
                else:
                    ax.set_ylim(0.0, max(max(t) * 1.08, 1e-3))

                ax.set_xlabel(xlabel or r"$E - E_\mathrm{F}$ (eV)")
                ax.set_ylabel(ylabel or r"$T(E)$ ($G_0$)")
                if title:
                    ax.set_title(title, fontsize=STYLE["axes.titlesize"])
                _panel_label(ax, panel)
                if annotate and summary is not None:
                    _annotate_stats(ax, summary, gap=gap)
                if label:
                    ax.legend(loc="upper right", frameon=True, facecolor="white",
                              edgecolor="#b0b0b0", framealpha=0.92, borderpad=0.3)

                axes = ax
                paths = _save(fig, out_path, formats, dpi)
            finally:
                if not return_figure:
                    plt.close(fig)
        else:
            ncols = 2 if mode == "polarization" else 1
            figsize_eff = _spin_figsize(figsize) if ncols == 2 else size
            fig, axes_raw = plt.subplots(
                1, ncols, figsize=figsize_eff, sharex=True,
                gridspec_kw=({"width_ratios": [1.45, 1.0]} if ncols == 2 else None))
            ax_list = list(axes_raw) if ncols > 1 else [axes_raw]
            try:
                ax_t = ax_list[0]
                for values, curve_color, curve_label in (
                        (rows["t_alpha"], color or ALPHA_COLOR, r"$T_\alpha$"),
                        (rows["t_beta"], BETA_COLOR, r"$T_\beta$")):
                    y_plot = [max(v, floor) for v in values] if scale == "log" else values
                    ax_t.plot(e, y_plot, color=curve_color, lw=STYLE["lines.linewidth"],
                              label=curve_label, zorder=3)
                _prepare_axes(ax_t, scale)
                _fermi_line(ax_t)

                if xlim:
                    ax_t.set_xlim(float(xlim[0]), float(xlim[1]))
                else:
                    ax_t.set_xlim(min(e), max(e))
                if ylim:
                    ax_t.set_ylim(float(ylim[0]), float(ylim[1]))
                elif scale == "log":
                    ax_t.set_ylim(*_log_ylim(
                        list(rows["t_alpha"]) + list(rows["t_beta"]), floor))
                else:
                    ax_t.set_ylim(0.0, max(max(rows["t_alpha"] + rows["t_beta"]) * 1.08, 1e-3))

                ax_t.set_xlabel(xlabel or r"$E - E_\mathrm{F}$ (eV)")
                ax_t.set_ylabel(ylabel or r"$T(E)$ ($G_0$)")
                if title:
                    ax_t.set_title(title, fontsize=STYLE["axes.titlesize"])
                _panel_label(ax_t, panel)
                if annotate and summary is not None:
                    _annotate_spin_channels(ax_t, rows, summary)
                ax_t.legend(loc="upper right", frameon=True, facecolor="white",
                            edgecolor="#b0b0b0", framealpha=0.92, borderpad=0.3)

                if ncols == 2:
                    ax_p = ax_list[1]
                    p = _spin_polarization(rows)
                    p_pct = [100.0 * value for value in p]
                    ax_p.axhline(0.0, color="#b0b0b0", lw=0.8,
                                 linestyle=(0, (4, 2)), zorder=1)
                    ax_p.plot(e, p_pct, color=POLARIZATION_COLOR,
                              lw=STYLE["lines.linewidth"], zorder=3)
                    _prepare_axes(ax_p, "linear")
                    _fermi_line(ax_p, label=False)
                    span = max(max(abs(v) for v in p_pct), 1e-4) * 1.18
                    ax_p.set_ylim(-span, span)
                    ax_p.set_xlabel(xlabel or r"$E - E_\mathrm{F}$ (eV)")
                    ax_p.set_ylabel(r"$P(E)$ (%)")
                    _panel_label(ax_p, _next_panel(panel))
                    if annotate:
                        _annotate_polarization(ax_p, rows)

                axes = ax_list if ncols > 1 else ax_t
                paths = _save(fig, out_path, formats, dpi)
            finally:
                if not return_figure:
                    plt.close(fig)

    meta: dict[str, Any] = {
        "path": str(paths[0]),
        "paths": [str(p) for p in paths],
        "points": len(e),
        "log_scale": scale == "log",
        "scale": scale,
        "floor": floor,
        "style_version": STYLE_VERSION,
        "spin": mode,
        "channels": channels,
    }
    if return_figure:
        meta["figure"] = (fig, axes)
    return meta


# --------------------------------------------------------------------------- #
# 多曲线叠加（同一坐标系；固定配色与图例位置）
# --------------------------------------------------------------------------- #
def plot_te_multi(
    curves: Sequence[Any],
    out: str | Path | None = None,
    *,
    labels: Sequence[str | None] | None = None,
    colors: Sequence[str | None] | None = None,
    scale: str = "log",
    figsize: str | Sequence[float] = "double",
    formats: Sequence[str] | None = None,
    dpi: float = 600,
    xlim: Sequence[float] | None = None,
    ylim: Sequence[float] | None = None,
    title: str | None = None,
    panel: str | None = None,
    legend_loc: str = "upper right",
    return_figure: bool = False,
) -> dict[str, Any]:
    """把多条 T(E) 画进同一坐标系（对比不同 run / 构象 / 偏压）。

    ``curves`` 的每一项可以是路径，或 ``(路径, 标签)``，或
    ``(路径, 标签, 颜色)``。
    """
    if scale not in ("log", "linear"):
        raise ValueError("scale 只能是 'log' 或 'linear'，收到 %r" % (scale,))
    entries = []
    for item in curves:
        if isinstance(item, (str, Path)):
            entries.append((item, None, None))
        else:
            parts = list(item)
            entries.append((parts[0], parts[1] if len(parts) > 1 else None,
                            parts[2] if len(parts) > 2 else None))
    if not entries:
        raise ValueError("curves 为空")
    if labels is not None and len(labels) != len(entries):
        raise ValueError("labels 与 curves 数量不一致")

    rows = [read_curve(src) for src, _, _ in entries]
    floor = log_floor(v for row in rows for v in row["t"])
    size = FIG_SIZES[figsize] if isinstance(figsize, str) else tuple(float(v) for v in figsize)

    plt = _pyplot()

    with te_style():
        fig, ax = plt.subplots(figsize=size)
        try:
            all_t: list[float] = []
            for index, ((src, label0, color0), row) in enumerate(zip(entries, rows)):
                label = (labels[index] if labels is not None else label0) or _default_label(src)
                color = color0 or PALETTE[index % len(PALETTE)]
                t = row["t"]
                all_t.extend(t)
                y_plot = [max(v, floor) for v in t] if scale == "log" else t
                ax.plot(row["e"], y_plot, color=color, lw=STYLE["lines.linewidth"],
                        label=label, zorder=3 + index)
            _prepare_axes(ax, scale)
            _fermi_line(ax)
            xs = [x for row in rows for x in row["e"]]
            if xlim:
                ax.set_xlim(float(xlim[0]), float(xlim[1]))
            else:
                ax.set_xlim(min(xs), max(xs))
            if ylim:
                ax.set_ylim(float(ylim[0]), float(ylim[1]))
            elif scale == "log":
                ax.set_ylim(*_log_ylim(all_t, floor))
            else:
                ax.set_ylim(0.0, max(max(all_t) * 1.08, 1e-3))
            ax.set_xlabel(r"$E - E_\mathrm{F}$ (eV)")
            ax.set_ylabel(r"$T(E)$ ($G_0$)")
            if title:
                ax.set_title(title, fontsize=STYLE["axes.titlesize"])
            _panel_label(ax, panel)
            if len(entries) > 1 or any(entry[1] for entry in entries):
                ax.legend(loc=legend_loc, frameon=True, facecolor="white",
                          edgecolor="#b0b0b0", framealpha=0.92, borderpad=0.3)

            if out is None:
                out = Path(rows[0]["source"]).parent / ("T_E_overlay.%s" % (formats[0] if formats else "png"))
            paths = _save(fig, out, formats, dpi)
        finally:
            if not return_figure:
                plt.close(fig)

    meta: dict[str, Any] = {
        "path": str(paths[0]),
        "paths": [str(p) for p in paths],
        "curves": [
            {"source": row["source"], "label": (labels[i] if labels is not None else entries[i][1])
             or _default_label(entries[i][0]), "points": len(row["e"])}
            for i, row in enumerate(rows)
        ],
        "log_scale": scale == "log",
        "scale": scale,
        "floor": floor,
        "style_version": STYLE_VERSION,
    }
    if return_figure:
        meta["figure"] = (fig, ax)
    return meta
