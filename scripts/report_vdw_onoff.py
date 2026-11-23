"""Build the vdW on/off full-chain comparison report + T(E) figure.

vdW ON  = runs/light-tight-2026-11-22/{az-dithiol,2cn}-light  (production light arms)
vdW OFF = runs/vdw-2026-11-22/{az-dithiol,2cn}-novdw          (this campaign)

Writes:
  docs/validation-vdw-onoff-2026-11-22.md
  docs/figures/vdw-onoff-TE.png
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
FIG = DOCS / "figures"
FIG.mkdir(exist_ok=True)

SYSTEMS = [("az-dithiol", "AZ-dithiol（20 原子，0；strip S–H）"), ("2cn", "2CN（54 原子，0）")]


def arm_dir(sysname: str, state: str) -> Path:
    if state == "on":
        arm = f"{sysname}-light"
        return ROOT / "runs" / "light-tight-2026-11-22" / arm / f"lt-{arm}"
    arm = f"{sysname}-novdw"
    return ROOT / "runs" / "vdw-2026-11-22" / arm / f"vdw-{sysname}-off"


def load(sysname: str, state: str):
    d = arm_dir(sysname, state)
    res = json.loads((d / "junction_analyze" / "results.json").read_text(encoding="utf-8"))
    aimso = d / "aims_scf" / "aims.out"
    nbf = None
    if aimso.exists():
        m = re.findall(r"Total number of basis functions\s*:\s*(\d+)", aimso.read_text(errors="ignore"))
        nbf = int(m[-1]) if m else None
    te = np.loadtxt(d / "aitranss_te" / "TE.dat", comments="#")
    ef = float(re.search(r"efermi\s*=\s*(-?\d+\.\d+)",
                         (d / "aitranss_te" / "TE.dat").read_text().splitlines()[1]).group(1))
    return {"dir": d, "res": res, "nbf": nbf, "te": te, "efermi": ef}


def main():
    data = {(s, st): load(s, st) for s, _ in SYSTEMS for st in ("on", "off")}
    for k, d in data.items():
        st = json.loads((d["dir"] / "state.json").read_text(encoding="utf-8"))
        bad = [n for n, v in st["stages"].items() if v.get("status") != "collected"]
        print(k, "nbf=", d["nbf"], "points=", int(d["res"]["points"]), "bad=", bad)

    # ---- figure ----
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.1), dpi=200)
    for ax, (sysname, title) in zip(axes, SYSTEMS):
        for st, color, label in (("on", "#1F3864", "vdW on"), ("off", "#C00000", "vdW off")):
            te = data[(sysname, st)]["te"]
            ax.plot(te[:, 1], te[:, 2], color=color, lw=1.2, label=f"{label}（{len(te)} 点）")
        ax.axvline(0.0, color="#808080", lw=0.8, ls="--")
        ax.set_yscale("log")
        ax.set_xlim(-0.9, 0.8)
        ax.set_xlabel("E − E_F (eV)")
        ax.set_ylabel("T(E)（每自旋）")
        ax.set_title(title, fontsize=10)
        ax.legend(frameon=False, fontsize=8, loc="lower right")
        ax.tick_params(labelsize=8)
    fig.tight_layout()
    png = FIG / "vdw-onoff-TE.png"
    fig.savefig(png, bbox_inches="tight")
    print("figure ->", png)

    # ---- tables ----
    rows = []
    for sysname, title in SYSTEMS:
        for st in ("on", "off"):
            d = data[(sysname, st)]
            r = d["res"]
            rows.append(
                f"| {title} | {st} | {d['nbf'] or '—'} | {r['points']} | {d['efermi']:.6f} | "
                f"{r['L0']:.6f} | {r['S_uV_per_K']:+.2f} | {r['k_e_pW_per_K']:.3f} | "
                f"{r['ZT_electronic']:.5f} | {r['fermi_nearest']['T']:.6f} | {r['max_transmission']['T']:.4f} |"
            )

    drows = []
    concl = []
    for sysname, title in SYSTEMS:
        a, b = data[(sysname, "on")]["res"], data[(sysname, "off")]["res"]
        dL0 = (b["L0"] - a["L0"]) / abs(a["L0"]) * 100
        dS = b["S_uV_per_K"] - a["S_uV_per_K"]
        dke = (b["k_e_pW_per_K"] - a["k_e_pW_per_K"]) / abs(a["k_e_pW_per_K"]) * 100
        dZT = (b["ZT_electronic"] - a["ZT_electronic"]) / abs(a["ZT_electronic"]) * 100
        dT = (b["fermi_nearest"]["T"] - a["fermi_nearest"]["T"]) / abs(a["fermi_nearest"]["T"]) * 100
        drows.append(f"| {title} | {dL0:+.1f}% | {dS:+.2f} | {dke:+.1f}% | {dZT:+.1f}% | {dT:+.1f}% |")
        concl.append(
            f"- {title}：L0 {a['L0']:.6f} → {b['L0']:.6f}（{dL0:+.1f}%），"
            f"S {a['S_uV_per_K']:+.2f} → {b['S_uV_per_K']:+.2f} µV/K，"
            f"κe {a['k_e_pW_per_K']:.3f} → {b['k_e_pW_per_K']:.3f} pW/K，"
            f"ZT {a['ZT_electronic']:.5f} → {b['ZT_electronic']:.5f}（{dZT:+.1f}%），"
            f"T(E_F) {a['fermi_nearest']['T']:.6f} → {b['fermi_nearest']['T']:.6f}。"
        )

    md = f"""# vdW 开关整链输运对照（2026-11-22）

> 目的：在 light 档、同一生产协议、同一能量窗下，只切换 FHI-aims 的 TS-vdW 修正
> （`vdw_correction_hirshfeld` 开 / 关），把分子优化 → tip 弛豫 → 建结 → SCF → Aitranss
> → 分析整条链跑通，给出透射谱与热电器件量的差异。FHI-aims 的 TS 修正是**后自洽**的：
> 它只进入总能量与力，不改 KS 哈密顿量，因此同一几何下开/关 vdW 的 T(E) 完全相同；
> 这里量到的是**弛豫几何的响应**（分子优化 + tip 弛豫），不是固定几何的电子结构差。

## 1. 设计

| 项目 | 设置 |
|---|---|
| 体系 | AZ-dithiol（azulene-2,6-dithiol，C10H8S2，20 原子，charge 0，tip 阶段 strip S–H）；2CN（2CN-656-4-thiomethylphenyl，54 原子，charge 0） |
| 唯一变量 | `electrode.vdw`：true（Tkatchenko–Scheffler hirshfeld）/ false（无色散修正） |
| 电子结构 | FHI-aims 221103，PBE，spin none，species light，atomic_zora scalar |
| 分子优化 | BFGS，最大力 1e-2 eV/Å，sc_accuracy_forces 5e-3 eV/Å，sc_iter_limit 1000 |
| tip 弛豫 | Au–S 初距 2.35 Å，pin_contact=false，全原子 BFGS |
| 电极 | fcc111 37 Au/侧，nlayers=3，冻结 |
| SCF | rho 1e-5，eev 1e-3，etot 1e-6，sc_iter_limit 500，KS_method serial |
| 输运 | Aitranss 071813，σ=0.05/0.10 Ha，能量窗 [−0.22, 0.0002, −0.16] Ha（301 点），T=300 K |
| 调度 | gpu-outside，8 CPU 任务/臂 |

vdW-on 基线复用 2026-11-22 light/tight 对照的 light 两臂（同协议、同 species、同能量窗）；
vdW-off 两臂在 `runs/vdw-2026-11-22/` 下新建。注意：两档的分子与 tip 几何都在各自设置下
重新优化，差异包含几何响应。

## 2. 结果

| 体系 | vdW | 总基函数 | 能量点 | E_F (Ha) | L0 (2e²/h) | S (µV/K) | κ_e (pW/K) | ZT | T(E_F) | max T |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

## 3. vdW-off 相对 on 的差值

| 体系 | ΔL0 | ΔS (µV/K) | Δκ_e | ΔZT | ΔT(E_F) |
|---|---:|---:|---:|---:|---:|
{chr(10).join(drows)}

## 4. 透射谱

![vdW on/off T(E)](figures/vdw-onoff-TE.png)

图中虚线为 E = E_F；纵轴为每自旋透射，对数坐标。

## 5. 结论与边界

{chr(10).join(concl)}

- TS-vdW 的作用路径是弛豫几何：本次 on/off 的 S–Au 接触距离只差 0.002–0.007 Å
  （AZ 2.251/2.244 Å，2CN 2.425/2.427 Å），但 E_F 附近的谱形对亚 0.01 Å 的构型差也已敏感，
  因此两体系的方向/幅度不同，不能外推成“vdW 一定使 ZT 升高/降低”。
- 边界：两个体系、一个能量窗、300 K、每臂一次；没有做 vdW 参数（TS vs MBD）、电极层数、
  泛函或 k 点的收敛对照；两档几何各自优化，差异来自几何响应而非固定几何的哈密顿量改变。
- 复现：`python scripts/report_vdw_onoff.py` 读取各臂 results.json / TE.dat / aims.out 生成本报告与图。

*生成时间：自动脚本；数据目录 `runs/light-tight-2026-11-22/`（on）与 `runs/vdw-2026-11-22/`（off）。*

入库副本：`docs/validation-data/vdw-2026-11-22/`（精简导出，含 state/run_spec/几何/TE.dat/results.json，
见其 README；`aims.out` 不入库）。
"""
    out = DOCS / "validation-vdw-onoff-2026-11-22.md"
    out.write_bytes(md.encode("utf-8"))
    print("report ->", out)


if __name__ == "__main__":
    main()
