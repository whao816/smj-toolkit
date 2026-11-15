"""烷烃二硫醇（HDT/ODT/DDT）构象散布：S–S 端距 → 等效碳数 → G 散布因子。

为什么做这个：判据族 `Au-alkane-dithiol-Au` 用"每碳指数衰减"的 β≈1.0/CH2 标定
（Xu & Tao, Science 2003），但那是针对 *全反式* 链长说的；柔性烷链的 gauche
缺陷会把 S–S 端距缩短 ~1 Å 量级，等效于少一个 CH2，G 就差 e^β ≈ 2.7 倍。
单个手画构象的 G 因此不能直接和族里的链长序列对齐。

口径
----
- 构象：ETKDGv3 + MMFF94s，重原子骨架 RMSD 0.25 Å 去重，按 MMFF 能量升序。
- d(S–S)：两个末端 S 的欧氏距离（去质子后 S 与 Au 成键，端距是结长的几何代理）。
- 等效碳数增量：Δn = (d_SS − d_SS,max) / 1.27 Å（全反式每 CH2 沿链轴增量 ~1.27 Å）。
- G 散布因子：exp(β·Δn)，β = 1.0 /CH2（Xu & Tao 2003 的量级）；β 只作量级估计，
  不是从本仓库算例拟合出来的。
- 只做几何与 MMFF，不含 DFT；DFT 级的 gap 散布要另跑 `smj ensemble`。

用法::

    python scripts/dithiol_conformer_spread.py --out runs/compare-2026-11-15
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
from pathlib import Path

from smj.conformers import select_conformers

#: 判据族 `Au-alkane-dithiol-Au` 的链长序列（criteria.py）。
SYSTEMS = {
    "HDT": dict(smiles="SCCCCCCS", carbons=6, note="1,6-己二硫醇"),
    "ODT": dict(smiles="SCCCCCCCCS", carbons=8, note="1,8-辛二硫醇"),
    "DDT": dict(smiles="SCCCCCCCCCCS", carbons=10, note="1,10-癸二硫醇"),
}
C_PER_CH2_A = 1.27   # 全反式烷链每 CH2 沿链轴的投影增量（Å）
BETA_PER_CH2 = 1.0   # Xu & Tao 2003，每 CH2 的隧穿衰减因子
BOLTZMANN_KCAL = 0.0019872043


def _s_indices(atoms) -> list[int]:
    return [i for i, a in enumerate(atoms) if a.element == "S"]


def _ss_distance(atoms) -> float | None:
    idx = _s_indices(atoms)
    if len(idx) != 2:
        return None
    a, b = atoms[idx[0]], atoms[idx[1]]
    return math.dist((a.x, a.y, a.z), (b.x, b.y, b.z))


def analyse(system: str, spec: dict, *, n_conf: int = 60, keep: int = 12) -> dict:
    conformers = select_conformers(spec["smiles"], n_conformers=n_conf, keep=keep)
    rows = []
    for c in conformers:
        d = _ss_distance(c.atoms)
        rows.append({
            "rank": c.rank,
            "energy_kcal": round(c.energy_kcal, 3),
            "rmsd_to_first_A": round(c.rmsd_to_first_A, 3),
            "ss_distance_A": None if d is None else round(d, 3),
        })

    distances = [r["ss_distance_A"] for r in rows if r["ss_distance_A"] is not None]
    if not distances:
        raise ValueError(f"{system}: 没有找到两个硫的构象")
    d_min, d_max = min(distances), max(distances)
    d_extended = d_max  # keep 集合里最长的那个即最接近全反式

    deltas = [(d - d_extended) / C_PER_CH2_A for d in distances]
    energies = [r["energy_kcal"] for r in rows]
    weights = []
    e_min = min(energies)
    for e in energies:
        weights.append(math.exp(-(e - e_min) / (BOLTZMANN_KCAL * 298.15)))
    z = sum(weights)
    weights = [w / z for w in weights]

    eff = []
    for delta, w in zip(deltas, weights):
        eff.append({
            "delta_carbons": round(delta, 3),
            "g_ratio": round(math.exp(BETA_PER_CH2 * delta), 3),
            "boltzmann_weight_298K": round(w, 3),
        })

    g_values = [e["g_ratio"] for e in eff]
    weighted_g = sum(e["g_ratio"] * e["boltzmann_weight_298K"] for e in eff)
    return {
        "system": system,
        "note": spec["note"],
        "smiles": spec["smiles"],
        "carbons": spec["carbons"],
        "n_conformers": len(rows),
        "ss_distance_A": {
            "min": round(d_min, 3),
            "max": round(d_max, 3),
            "mean": round(statistics.fmean(distances), 3),
            "spread": round(d_max - d_min, 3),
        },
        "delta_carbons_spread": round(max(deltas) - min(deltas), 3),
        "g_ratio": {
            "slackest": round(max(g_values), 3),   # 最短构象，G 最大
            "tightest": round(min(g_values), 3),
            "factor": round(max(g_values) / min(g_values), 3),
            "boltzmann_mean": round(weighted_g, 3),
        },
        "energy_kcal": {
            "min": round(min(energies), 3),
            "max": round(max(energies), 3),
            "spread": round(max(energies) - min(energies), 3),
        },
        "conformers": rows,
        "effective": eff,
    }


def markdown(results: list[dict]) -> str:
    lines = [
        "# HDT/ODT/DDT 构象散布（几何 + MMFF，非 DFT）",
        "",
        "口径：ETKDGv3 + MMFF94s，重原子骨架 RMSD 0.25 Å 去重，MMFF 能量升序取前 12；",
        "`d(S–S)` 是去质子后结长的几何代理；`Δn = (d−d_max)/1.27 Å`，`G 因子 = e^{βΔn}`，"
        "β = 1.0/CH2（Xu & Tao 2003 量级）。",
        "",
        "| 体系 | 碳数 | 构象 | d(S–S) 范围 (Å) | 极差 (Å) | Δn 极差 | G 因子范围 | G 极差倍数 | MMFF 能量极差 (kcal/mol) |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for r in results:
        lines.append(
            "| {s} | {c} | {n} | {a:.3f} – {b:.3f} | {sp:.3f} | {dn:.3f} | {g1:.2f} – {g2:.2f} | {gf:.2f}× | {e:.2f} |".format(
                s=r["system"], c=r["carbons"], n=r["n_conformers"],
                a=r["ss_distance_A"]["min"], b=r["ss_distance_A"]["max"],
                sp=r["ss_distance_A"]["spread"], dn=r["delta_carbons_spread"],
                g1=r["g_ratio"]["tightest"], g2=r["g_ratio"]["slackest"],
                gf=r["g_ratio"]["factor"], e=r["energy_kcal"]["spread"]))
    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="runs/compare-2026-11-15", help="输出目录")
    ap.add_argument("--n", type=int, default=60, help="每个体系的 ETKDG 候选数")
    ap.add_argument("--keep", type=int, default=12, help="去重后保留数")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    results = [analyse(s, spec, n_conf=args.n, keep=args.keep) for s, spec in SYSTEMS.items()]
    payload = {
        "generator": "ETKDGv3 + MMFF94s；RMSD 0.25 Å 重原子骨架去重",
        "caveat": "只做几何 + 力场；DFT 级 gap 散布需另跑 smj ensemble。β=1.0/CH2 为文献量级，未在本仓库拟合。",
        "systems": results,
    }
    (out / "dithiol_conformers.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    (out / "dithiol_conformers.md").write_text(markdown(results), encoding="utf-8")

    for r in results:
        print(f"{r['system']}: {r['n_conformers']} conf, d(S-S) {r['ss_distance_A']['min']:.3f}"
              f"-{r['ss_distance_A']['max']:.3f} A, G factor {r['g_ratio']['factor']:.2f}x,"
              f" E spread {r['energy_kcal']['spread']:.2f} kcal/mol")
    print(f"\n写出 {out/'dithiol_conformers.json'} 与 {out/'dithiol_conformers.md'}")


if __name__ == "__main__":
    main()
