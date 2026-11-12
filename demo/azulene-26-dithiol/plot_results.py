r"""把一次 SMJ run 的 TE.dat / analysis.json 画成一张概览图。

用法（在 demo 工程目录下）：
  ..\..\.venv\Scripts\python.exe plot_results.py <run目录>

输出：<run目录>/junction_analyze/TE_overview.png（与原 TE.pdf 并存，便于快查）
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def read_te(path: Path) -> tuple[np.ndarray, np.ndarray]:
    energies, transmissions = [], []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        text = line.split("#", 1)[0].strip()
        if not text:
            continue
        parts = text.replace(",", " ").split()
        if len(parts) < 2:
            continue
        try:
            # aitranss 的 TE.dat 三列：E[Hartree]  E-EF[eV]  T(E)；两列时第一列为横轴
            energies.append(float(parts[1]) if len(parts) >= 3 else float(parts[0]))
            transmissions.append(float(parts[-1]))
        except ValueError:
            continue
    return np.asarray(energies), np.asarray(transmissions)


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 2
    run_dir = Path(argv[1]).resolve()
    te = run_dir / "aitranss_te" / "TE.dat"
    analysis_path = run_dir / "junction_analyze" / "analysis.json"
    if not te.is_file():
        raise SystemExit("找不到 %s" % te)
    energy, transmission = read_te(te)
    if energy.size == 0:
        raise SystemExit("TE.dat 没有可解析的数据行: %s" % te)

    analysis = json.loads(analysis_path.read_text(encoding="utf-8")) if analysis_path.is_file() else {}
    nearest = analysis.get("fermi_nearest") or {}
    e_ef = nearest.get("E_minus_Ef_eV") if isinstance(nearest, dict) else None
    fermi = 0.0 if e_ef is None else float(e_ef)

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.0), sharex=True,
                             gridspec_kw={"height_ratios": [3, 1.2]})
    ax = axes[0]
    ax.plot(energy, transmission, color="#1f4e79", lw=1.2)
    ax.axvline(fermi, color="#c00000", lw=0.9, ls="--")
    ax.set_ylabel("T(E)")
    ax.set_title("azulene-26-dithiol | Au(111) library electrode  (E - E_F, eV)")
    ax.grid(alpha=0.25)

    text = []
    for key, label, fmt in (
            ("L0", "L0 (1/eV)", "%.4f"),
            ("S_uV_per_K", "S (µV/K)", "%+.2f"),
            ("k_e_pW_per_K", "k_e (pW/K)", "%.1f"),
            ("ZT_electronic", "ZT_e", "%.4f"),
            ("max_transmission", "max T", "%.3f"),
    ):
        if key in analysis:
            text.append("%s = %s" % (label, fmt % float(analysis[key])))
    if "log10_G_over_G0_like" in analysis:
        g_over_g0 = 10.0 ** float(analysis["log10_G_over_G0_like"])
        text.insert(1, "G / G0 = %.4f" % g_over_g0)
        text.insert(2, "G (µS) = %.2f" % (g_over_g0 * 77.4809))
    if text:
        ax.text(0.02, 0.97, "\n".join(text), transform=ax.transAxes, va="top",
                fontsize=8.5, family="monospace",
                bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#999999", alpha=0.9))

    ax2 = axes[1]
    dx = np.gradient(energy) if energy.size > 1 else np.ones_like(energy)
    dedx = np.gradient(transmission) / np.where(dx == 0, np.nan, dx)
    ax2.plot(energy, -dedx, color="#7f6000", lw=1.0)
    ax2.axvline(fermi, color="#c00000", lw=0.9, ls="--")
    ax2.set_xlabel("E - E_F (eV)")
    ax2.set_ylabel("-dT/dE")
    ax2.grid(alpha=0.25)
    fig.tight_layout()

    out = run_dir / "junction_analyze" / "TE_overview.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200)
    print("已写 %s（%d 个点，E_F=%.4f eV）" % (out, energy.size, fermi))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
