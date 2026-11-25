"""画单个 run 的声子透射谱与累积热导，用于论文/报告图.

输入 ``thermal.json``（``smj.thermal.v1``），输出两联图：

* (a) ``T_ph(E)`` 本身，标出 lead 声子带顶（Au 4.71 THz = 19.5 meV）并把带顶
  以上阴影化——那里弹性透射**必须**为零（方案 §3.6 的 fail-closed 判据）。
* (b) 归一化累积热导，标出 25/50/75/90% 的贡献分位。单个 G_ph 数字看不出
  热流从哪来；这条曲线直接回答"哪一段声子把热带过去"，是解释分子间差异
  的主要依据（低频 vs 骨架模）。

用法示例::

    python scripts/plot_thermal_spectrum.py \\
        --thermal runs/nap-sh-thermal/nap-sh-thermal/thermal_analyze/thermal.json \\
        --output runs/nap-sh-thermal/nap-sh-thermal/thermal_analyze \\
        --title "naphthalene-2,6-dithiol / Au(111) 3x3"
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)

#: lead 声子带顶参考值（meV）。fcc-Au 最高支在 L 点 4.71 THz，见
#: Lynn/Smith/Nicklow, PRB 8, 3493 (1973)。只作为**参考线**用；真正的判据
#: 用的是 lead 模型自己的带顶（thermal.json 里的 ``band_top_meV``）。
AU_BAND_TOP_MEV = 19.5

#: 累积曲线标注的分位（占 G_ph 总量的百分比）。
QUANTILES = (25.0, 50.0, 75.0, 90.0)

#: 论文图统一样式。
STYLE = {
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "legend.fontsize": 9.5,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.4,
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--thermal", type=Path, required=True,
                   help="thermal.json 路径（smj.thermal.v1）")
    p.add_argument("--output", type=Path, default=None,
                   help="输出目录（默认与 thermal.json 同目录）")
    p.add_argument("--stem", type=str, default="thermal_transmission",
                   help="输出文件名主干（默认 thermal_transmission）")
    p.add_argument("--title", type=str, default=None,
                   help="图题；缺省用 thermal.json 里的体系信息拼一个")
    p.add_argument("--xmax-meV", type=float, default=22.0,
                   help="横轴上限（meV，默认 22，刚好把带顶与带外都框进来）")
    p.add_argument("--formats", type=str, default="png,pdf",
                   help="输出格式，逗号分隔（默认 png,pdf）")
    p.add_argument("--verbose", action="store_true", help="输出调试日志")
    return p.parse_args()


def load_spectrum(path: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    """读 thermal.json，返回 ``(能量 meV, T_ph, 顶层元数据)``。"""
    payload = json.loads(path.read_text(encoding="utf-8"))
    schema = payload.get("schema")
    if schema != "smj.thermal.v1":
        raise ValueError("schema 应为 smj.thermal.v1，得到 %r（%s）" % (schema, path))
    spectrum = payload.get("spectrum") or {}
    try:
        energy = np.asarray(spectrum["energy_meV"], dtype=float)
        t_ph = np.asarray(spectrum["transmission"], dtype=float)
    except KeyError as exc:
        raise ValueError("thermal.json 缺 spectrum.%s" % exc.args[0]) from exc
    if energy.shape != t_ph.shape or energy.size < 2:
        raise ValueError("能量与透射谱长度不一致或点数不足: %r vs %r（%s）"
                         % (energy.shape, t_ph.shape, path))
    if not (np.all(np.isfinite(energy)) and np.all(np.isfinite(t_ph))):
        raise ValueError("谱里含非有限数值（%s）" % path)
    return energy, t_ph, payload


def cumulative_percent(energy: np.ndarray, t_ph: np.ndarray,
                       temperature_k: float, band_top_mev: float) -> np.ndarray:
    """用与阶段内完全相同的积分器算累积热导，再归一化成百分比.

    复用 ``smj.phonon.landauer.cumulative_thermal_conductance``，避免图里的
    曲线与 ``thermal.json`` 里的数字出自两套积分。
    """
    from smj.phonon.landauer import cumulative_thermal_conductance

    energy_ev = energy * 1.0e-3
    cumulative = cumulative_thermal_conductance(
        energy_ev, t_ph, temperature_k) * 1.0e12
    total = float(cumulative[-1])
    if total <= 0.0:
        raise ValueError("累积热导为 0，无法归一化（G_ph=%r pW/K）" % total)
    return cumulative / total * 100.0


def crossing(energy: np.ndarray, cumulative_pct: np.ndarray,
             level: float) -> float:
    """累积曲线的 ``level``% 分位对应的能量（meV），线性插值."""
    index = int(np.searchsorted(cumulative_pct, level))
    if index <= 0:
        return float(energy[0])
    if index >= cumulative_pct.size:
        return float(energy[-1])
    y0, y1 = cumulative_pct[index - 1], cumulative_pct[index]
    if y1 == y0:
        return float(energy[index])
    frac = (level - y0) / (y1 - y0)
    return float(energy[index - 1] + frac * (energy[index] - energy[index - 1]))


def make_figure(energy: np.ndarray, t_ph: np.ndarray, payload: dict,
                title: str | None, xmax_mev: float):
    """画 (a) 透射谱 + (b) 归一化累积热导两联图."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with matplotlib.rc_context(STYLE):
        temperature_k = float(payload.get("temperature_K", 300.0))
        band_top = float(payload.get("band_top_meV", AU_BAND_TOP_MEV))
        g_ph = float(payload.get("g_ph_pw_per_k", float("nan")))
        g_q = float(payload.get("g_quantum_pw_per_k", float("nan")))
        n_eff = float(payload.get("n_channels_effective", float("nan")))
        eta = float(payload.get("eta_meV", float("nan")))

        keep = energy <= xmax_mev
        energy, t_ph = energy[keep], t_ph[keep]
        cumulative_pct = cumulative_percent(
            energy, t_ph, temperature_k, band_top)

        fig, (ax_top, ax_bottom) = plt.subplots(
            2, 1, figsize=(7.2, 6.6), sharex=True,
            gridspec_kw={"height_ratios": [1.35, 1.0], "hspace": 0.08})

        # ---- (a) 透射谱 ----
        ax_top.plot(energy, t_ph, color="#1f4e79", zorder=3)
        ax_top.axvspan(band_top, xmax_mev, color="#b0b0b0", alpha=0.22,
                       lw=0, zorder=1)
        ax_top.axvline(band_top, color="#7f7f7f", ls="--", lw=1.0, zorder=2)
        ax_top.set_ylabel(r"$T_{\mathrm{ph}}(\hbar\omega)$")
        ax_top.set_ylim(bottom=0.0)
        top = max(float(np.max(t_ph)), 1.0)
        ax_top.set_ylim(0.0, top * 1.12)
        ax_top.text(band_top - 0.35, top * 1.06,
                    "lead band top %.1f meV" % band_top,
                    fontsize=8.5, color="#4d4d4d", va="top", ha="right")
        ax_top.text(0.02, 0.96,
                    r"$G_{\mathrm{ph}}$(%d K) = %.2f pW/K" % (int(temperature_k), g_ph)
                    + "\n" + r"$G_Q$ = %.1f, $N_{\mathrm{eff}}$ = %.3f"
                    % (g_q, n_eff)
                    + "\n" + r"$\eta$ = %.2f meV" % eta,
                    transform=ax_top.transAxes, va="top", ha="left",
                    fontsize=9,
                    bbox={"facecolor": "white", "alpha": 0.85,
                          "edgecolor": "#cccccc", "boxstyle": "round,pad=0.35"})

        # ---- (b) 归一化累积热导 ----
        ax_bottom.plot(energy, cumulative_pct, color="#a83232", zorder=3)
        for level in QUANTILES:
            e_level = crossing(energy, cumulative_pct, level)
            ax_bottom.axhline(level, color="#a83232", ls=":", lw=0.8, alpha=0.45)
            ax_bottom.plot([e_level], [level], marker="o", ms=4.2,
                           color="#a83232", zorder=4)
            ax_bottom.annotate("%.0f%% @ %.1f meV" % (level, e_level),
                               xy=(e_level, level), xytext=(4, -11),
                               textcoords="offset points", fontsize=8.2,
                               color="#7a2020")
        ax_bottom.set_xlabel(r"phonon energy $\hbar\omega$ (meV)")
        ax_bottom.set_ylabel(r"cumulative $G_{\mathrm{ph}}$ (%)")
        ax_bottom.set_ylim(0.0, 102.0)
        ax_bottom.set_xlim(0.0, xmax_mev)

        if title:
            fig.suptitle(title, fontsize=11.5, y=0.965)
        fig.subplots_adjust(left=0.13, right=0.97, top=0.90, bottom=0.10)
        return fig


def main() -> None:
    args = parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    energy, t_ph, payload = load_spectrum(args.thermal)
    logger.info("读到 %d 点谱，E ∈ [%.2f, %.2f] meV，G_ph = %.3f pW/K",
                energy.size, energy.min(), energy.max(),
                float(payload.get("g_ph_pw_per_k", float("nan"))))

    output_dir = args.output or args.thermal.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    title = args.title
    if title is None:
        title = "%s  (N=%s atoms, %d K)" % (
            payload.get("system") or "single-molecule junction",
            payload.get("n_atoms", "?"), int(payload.get("temperature_K", 300)))

    fig = make_figure(energy, t_ph, payload, title, args.xmax_meV)
    written = []
    for suffix in [s.strip() for s in args.formats.split(",") if s.strip()]:
        path = output_dir / ("%s.%s" % (args.stem, suffix))
        fig.savefig(path, bbox_inches="tight")
        written.append(str(path))
        logger.info("写入 %s", path)

    # 谱数据落 CSV（表头带单位），便于外部复算/重画
    temperature_k = float(payload.get("temperature_K", 300.0))
    band_top = float(payload.get("band_top_meV", AU_BAND_TOP_MEV))
    keep = energy <= args.xmax_meV
    pct = cumulative_percent(energy[keep], t_ph[keep], temperature_k, band_top)
    csv_path = output_dir / ("%s.csv" % args.stem)
    lines = ["energy_meV,transmission,cumulative_g_ph_pct"]
    for e, t, c in zip(energy[keep], t_ph[keep], pct):
        lines.append("%.6f,%.8e,%.4f" % (e, t, c))
    csv_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    written.append(str(csv_path))
    logger.info("写入 %s", csv_path)

    # 参数快照，便于图能追回到具体数值文件
    snapshot = {
        "source": str(args.thermal.resolve()),
        "n_atoms": payload.get("n_atoms"),
        "temperature_K": payload.get("temperature_K"),
        "g_ph_pw_per_k": payload.get("g_ph_pw_per_k"),
        "g_quantum_pw_per_k": payload.get("g_quantum_pw_per_k"),
        "n_channels_effective": payload.get("n_channels_effective"),
        "band_top_meV": payload.get("band_top_meV"),
        "eta_meV": payload.get("eta_meV"),
        "energy_step_meV": payload.get("energy_step_meV"),
        "quantile_energies_meV": {
            "%d%%" % q: crossing(energy[energy <= args.xmax_meV],
                                 cumulative_percent(
                                     energy[energy <= args.xmax_meV],
                                     t_ph[energy <= args.xmax_meV],
                                     float(payload.get("temperature_K", 300.0)),
                                     float(payload.get("band_top_meV",
                                                       AU_BAND_TOP_MEV))),
                                 q)
            for q in QUANTILES},
        "written": written,
    }
    snapshot_path = output_dir / ("%s.meta.json" % args.stem)
    snapshot_path.write_text(
        json.dumps(snapshot, indent=2, ensure_ascii=False), encoding="utf-8")
    logger.info("写入 %s", snapshot_path)
    logger.info("分位: %s", snapshot["quantile_energies_meV"])


if __name__ == "__main__":
    main()
