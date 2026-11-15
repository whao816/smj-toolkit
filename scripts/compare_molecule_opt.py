"""汇总 `runs/compare-2026-11-15` 各臂 molecule_opt 结果成对照表。

每臂一个 run 目录（`<arm>/<run-name>/state.json` + `molecule_opt/opt.xyz`），
读 state 里的 results 与 opt.xyz 的几何量，输出 markdown + json。

用法::

    python scripts/compare_molecule_opt.py --root runs/compare-2026-11-15
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from smj import cube
from smj._vendor.junction.geometry import read_xyz


def _dist(a, b) -> float:
    return math.dist(a.coord(), b.coord())


def geometry_metrics(xyz_path: Path) -> dict:
    atoms = read_xyz(xyz_path)
    s_idx = [i for i, a in enumerate(atoms) if a.element == "S"]
    out: dict = {"atom_count": len(atoms), "s_atoms": len(s_idx)}
    if len(s_idx) == 2:
        out["ss_distance_A"] = round(_dist(atoms[s_idx[0]], atoms[s_idx[1]]), 3)
    cs = []
    sh = []
    for i in s_idx:
        s = atoms[i]
        for j, o in enumerate(atoms):
            if j == i:
                continue
            d = _dist(s, o)
            if o.element == "C" and d < 1.95:
                cs.append(d)
            elif o.element == "H" and d < 1.35:
                sh.append(d)
    if cs:
        out["c_s_mean_A"] = round(sum(cs) / len(cs), 3)
        out["c_s_count"] = len(cs)
    if sh:
        out["s_h_mean_A"] = round(sum(sh) / len(sh), 3)
    return out


#: 臂名 -> 体系标签（表里用，避免从 config 猜）。
MOLECULES = {
    "2cn-light": "2CN-656-硫醚 C30H20N2S2",
    "2cn-tight": "2CN-656-硫醚 C30H20N2S2",
    "az-thiol-light": "薁-2,6-二硫醇 C10H8S2",
    "az-thiol-tight": "薁-2,6-二硫醇 C10H8S2",
    "az-thiolate-light": "薁-2,6-二硫醇盐 C10H6S2^2-",
    "az-thiolate-tight": "薁-2,6-二硫醇盐 C10H6S2^2-",
    "az-aucap-light": "薁-2,6-二硫醇(Au 帽模型)",
}


def load_arm(arm_dir: Path) -> dict | None:
    runs = [p for p in arm_dir.iterdir() if p.is_dir() and (p / "state.json").is_file()]
    if not runs:
        return None
    run_dir = runs[0]
    state = json.loads((run_dir / "state.json").read_text(encoding="utf-8"))
    stage = state["stages"]["molecule_opt"]
    if stage.get("status") != "collected":
        return {"arm": arm_dir.name, "molecule": MOLECULES.get(arm_dir.name, ""),
                "status": stage.get("status")}
    res = stage["results"]["molecule_opt"]
    sd = run_dir / "molecule_opt"
    row = {
        "arm": arm_dir.name,
        "molecule": MOLECULES.get(arm_dir.name, ""),
        "status": "collected",
        "engine": res.get("engine"),
        "atom_count": res.get("atom_count"),
        "converged": res.get("converged"),
        "relax_steps": res.get("relax_steps"),
        "max_force_ev_per_A": res.get("max_force_ev_per_A"),
        "force_threshold_ev_per_A": res.get("force_threshold_ev_per_A"),
        "electron_count": res.get("electron_count"),
        "gap_ev": (res.get("frontier_orbitals") or {}).get("gap_ev"),
        "homo_ev": (res.get("frontier_orbitals") or {}).get("homo_ev"),
        "lumo_ev": (res.get("frontier_orbitals") or {}).get("lumo_ev"),
        "tier": (res.get("species") or {}).get("actual_tier"),
        "scf_forces": (res.get("scf") or {}).get("sc_accuracy_forces"),
    }
    xyz = sd / "opt.xyz"
    if xyz.is_file():
        row.update({k: v for k, v in geometry_metrics(xyz).items() if k not in row})
    weights = cube.frontier_weights(sd, anchor_elements=("S",))
    if weights:
        occ = row.get("electron_count")
        half = (occ if occ is not None else 0) // 2
        for w in weights:
            label = "HOMO" if w["eigenstate"] == half else (
                "LUMO" if w["eigenstate"] == half + 1 else None)
            row.setdefault("anchor_weights", []).append({
                "label": label, "eigenstate": w["eigenstate"],
                "anchor_fraction": w["anchor_fraction"]})
    return row


def markdown(rows: list[dict]) -> str:
    head = ("| 臂 | 体系 | 档位 | 原子 | charge 后电子数 | 弛豫步 | 末步最大力 (eV/Å) | "
            "gap (eV) | HOMO/LUMO (eV) | d(S–S) (Å) | C–S 均值 (Å) |")
    sep = "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"
    lines = [head, sep]
    for r in rows:
        if r["status"] != "collected":
            label = {"submitted": "未完成（见 REPORT §0）"}.get(r["status"], r["status"])
            lines.append(f"| {r['arm']} | {r.get('molecule','—') or '—'} | — | — | — | — | — | — | — | — | {label} |")
            continue
        hl = ""
        if r.get("homo_ev") is not None:
            hl = f"{r['homo_ev']:.2f} / {r['lumo_ev']:.2f}"
        lines.append("| {arm} | {mol} | {tier} | {n} | {e} | {steps} | {mf} | {gap} | {hl} | {ss} | {cs} |".format(
            arm=r["arm"], mol=r.get("molecule", ""), tier=r.get("tier"), n=r.get("atom_count"),
            e=r.get("electron_count"),
            steps=r.get("relax_steps") if r.get("relax_steps") is not None else "单点",
            mf=f"{r['max_force_ev_per_A']:.4g}" if r.get("max_force_ev_per_A") is not None else "—",
            gap=f"{r['gap_ev']:.3f}" if r.get("gap_ev") is not None else "—",
            hl=hl or "—",
            ss=f"{r['ss_distance_A']:.3f}" if r.get("ss_distance_A") else "—",
            cs=f"{r['c_s_mean_A']:.3f}" if r.get("c_s_mean_A") else "—"))
    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--root", default="runs/compare-2026-11-15")
    args = ap.parse_args()
    root = Path(args.root)
    rows = []
    for arm in sorted(p for p in root.iterdir() if p.is_dir() and p.name != "species-tight"):
        r = load_arm(arm)
        if r is not None:
            rows.append(r)
    # 单点臂（tight 单点落 light 优化几何）：没有 state.json，读 sp_result.json。
    for sp in sorted(root.glob("*/sp_result.json")):
        d = json.loads(sp.read_text(encoding="utf-8"))
        rows.append({
            "arm": d.get("arm", sp.parent.name), "molecule": d.get("molecule", ""),
            "status": "collected", "engine": "aims",
            "atom_count": d.get("atom_count"),
            "relax_steps": None, "max_force_ev_per_A": None, "tier": d.get("basis"),
            "electron_count": 246 if "2CN" in d.get("molecule", "") else 100,
            "gap_ev": d.get("gap_ev"), "homo_ev": d.get("homo_ev"), "lumo_ev": d.get("lumo_ev"),
            "note": "tight 单点 %s 几何" % ("light 优化" if "light" in d.get("geometry_from", "") else "?"),
        })
    payload = {"root": str(root), "arms": rows}
    (root / "molecule_opt_table.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    (root / "molecule_opt_table.md").write_text(markdown(rows), encoding="utf-8")
    print(markdown(rows))


if __name__ == "__main__":
    main()
