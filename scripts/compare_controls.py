"""Compare SMJ control runs: baseline vs dense grid vs constrained Au-S.

Reads each run's TE.dat / results.json / geometry and prints a markdown table,
writes JSON + a T(E) overlay figure under runs/_compare/.
"""
import json
import math
import sys
from pathlib import Path

ROOT = Path(r"D:/dsh_projects/AutoSMJ")
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

from smj._vendor.junction.criteria import significant_peaks
from smj._vendor.junction.geometry import read_aims_geometry

RUNS = {
    "baseline": ROOT / "runs/bdt/bdt",
    "dense": ROOT / "runs/bdt-dense/bdt-dense",
    "anchor2.42": ROOT / "runs/bdt-anchor24/bdt-anchor24",
    "anchor2.42-dense": ROOT / "runs/bdt-anchor24-dense/bdt-anchor24-dense",
}
OUT = ROOT / "runs/_compare"


def load_te(path):
    rows = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        f = s.split()
        if len(f) < 3:
            continue
        t = [float(x) for x in f[2:4]] if len(f) >= 4 else [float(f[2])]
        rows.append((float(f[0]), float(f[1]), sum(t) / len(t)))
    return rows


def t_at_fermi(rows):
    rows = sorted(rows, key=lambda r: r[1])
    e = np.array([r[1] for r in rows])
    t = np.array([r[2] for r in rows])
    i = int(np.argmin(np.abs(e)))
    if i == 0 or i == len(e) - 1:
        return float(t[i]), float(e[i])
    # linear interpolation between neighbours around E_F=0
    if e[i] <= 0 <= e[i + 1] or e[i - 1] <= 0 <= e[i]:
        j = i if e[i] <= 0 else i - 1
        x0, x1 = e[j], e[j + 1]
        y0, y1 = t[j], t[j + 1]
        return float(y0 + (y1 - y0) * (0 - x0) / (x1 - x0)), 0.0
    return float(t[i]), float(e[i])


def au_s(run):
    if not (run / "electrode_build/geometry.in").exists():
        fb = {"dense": "baseline", "anchor2.42-dense": "anchor2.42"}
        run = RUNS.get(fb.get(name, "baseline"), RUNS["baseline"])
    geo = read_aims_geometry(run / "electrode_build/geometry.in")
    meta = json.loads((run / "electrode_build/electrode_metadata.json").read_text(encoding="utf-8"))
    out = {}
    for side in ("left", "right"):
        a = geo[meta[f"{side}_tip_index"] - 1]
        b = geo[meta[f"{side}_linker_index"] - 1]
        out[side] = round(math.dist((a.x, a.y, a.z), (b.x, b.y, b.z)), 4)
    return out


def stage_times(run):
    st = json.loads((run / "state.json").read_text(encoding="utf-8"))
    return {k: v.get("elapsed_s") for k, v in st.get("stages", {}).items() if v.get("elapsed_s")}


summary = {}
for name, run in RUNS.items():
    if not (run / "aitranss_te/TE.dat").exists():
        summary[name] = {"status": "missing TE.dat"}
        continue
    rows = load_te(run / "aitranss_te/TE.dat")
    res = json.loads((run / "junction_analyze/results.json").read_text(encoding="utf-8"))
    tf, dn = t_at_fermi(rows)
    peaks = significant_peaks([(r[0], r[1], r[2]) for r in rows])
    summary[name] = {
        "points": len(rows),
        "estep_meV": round((rows[1][1] - rows[0][1]) * 1000.0, 4),
        "EF_Ha": round(rows[0][0] - rows[0][1] / 27.211386, 6),
        "T_fermi_interp": round(tf, 6),
        "T_fermi_nearest": res["fermi_nearest"]["T"],
        "nearest_dE_meV": round(res["fermi_nearest"]["E_minus_Ef_eV"] * 1000.0, 4),
        "L0": res["L0"],
        "log10": res["log10_G_over_G0_like"],
        "S_uV": res["S_uV_per_K"],
        "k_e_pW": res["k_e_pW_per_K"],
        "ZT_e": res["ZT_electronic"],
        "max_T": res["max_transmission"]["T"],
        "max_E_eV": res["max_transmission"]["E_minus_Ef_eV"],
        "peaks": [(round(p[0], 4), round(p[1], 5)) for p in peaks],
        "AuS": au_s(run),
        "times": stage_times(run),
    }

OUT.mkdir(parents=True, exist_ok=True)
(OUT / "compare.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

cols = ["L0", "S_uV", "k_e_pW", "ZT_e", "T_fermi_interp", "max_T", "points"]
avail = [n for n in RUNS if "L0" in summary[n]]
lines = ["| 指标 | " + " | ".join(avail) + " |",
         "|---|" + "---|" * len(avail)]
for c in cols:
    lines.append(f"| {c} | " + " | ".join(str(summary[n][c]) for n in avail) + " |")
for n in avail:
    s = summary[n]
    lines.append(f"| Au-S ({n}) / Å | left {s['AuS']['left']}, right {s['AuS']['right']} |" + " |" * (len(avail) - 1))
print("\n".join(lines))

# overlay T(E)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

styles = {
    "baseline": dict(color="#1f77b4", lw=0.9, marker="o", ms=2.2, label="baseline (301 pts)"),
    "dense": dict(color="#ff7f0e", lw=1.3, label="dense grid (1201 pts)"),
    "anchor2.42": dict(color="#2ca02c", lw=1.3, label="Au-S = 2.42 A (anchor24)"),
    "anchor2.42-dense": dict(color="#9467bd", lw=1.0, ls="--", label="Au-S = 2.42 A dense (1201 pts)"),
}
fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.6))
curves = {}
for n in avail:
    rows = load_te(RUNS[n] / "aitranss_te/TE.dat")
    e = np.array([r[1] for r in rows])
    t = np.array([r[2] for r in rows])
    curves[n] = (e, t)
    axes[0].plot(e, t, **styles[n])
    m = (e >= -0.35) & (e <= 0.35)
    axes[1].plot(e[m], t[m], **styles[n])

for ax in axes[:2]:
    ax.axvline(0, color="k", lw=0.8, ls="--")
    ax.axhline(0.27885397, color="gray", lw=0.7, ls=":")
    ax.set_xlabel("E - E_F (eV)")
    ax.set_ylabel("T(E)  [G0 per spin]")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
axes[0].set_title("full window")
axes[1].set_title("zoom |E-EF| <= 0.35 eV")

if "dense" in curves and "baseline" in curves:
    eb, tb = curves["baseline"]
    ed, td = curves["dense"]
    ti = np.interp(ed, eb, tb)
    axes[2].plot(ed, td - ti, color="#d62728", lw=0.9)
    axes[2].axvline(0, color="k", lw=0.8, ls="--")
    axes[2].set_xlabel("E - E_F (eV)")
    axes[2].set_ylabel("T_dense - T_baseline(interp)  [G0]")
    axes[2].set_title(f"grid convergence: max|dT| = {np.max(np.abs(td - ti)):.4f} G0")
    axes[2].grid(alpha=0.3)
axes[2].set_xlim(-0.35, 0.35)
fig.suptitle("BDT/Au(111): T(E) comparison (baseline vs dense grid vs Au-S=2.42 A)")
fig.tight_layout()
fig.savefig(OUT / "TE_compare.png", dpi=160)
print("wrote", OUT / "TE_compare.png", OUT / "compare.json")


