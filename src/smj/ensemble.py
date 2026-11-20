"""构象自动分支：把构象粗筛直接变成"每构象一个 run 目录"的批处理。

`smj conformers` 只给 top-k 结构；真正的科学问题是**构象对结果的影响量级**——
同一分子换一个构象，锚–锚距离可能差 0.3–1 Å，gap 与 T(E) 都会变。这里把
"生成 k 个构象 + 每个构象一条 `smj stage` 命令 + 跑完汇总散布"串起来：

* `write_ensemble(smiles, dir)` → `conf_XX.xyz`、`ensemble.json`、`run_ensemble.sh`；
  每条命令用 `--run-dir <root>/conf_XX` 隔离，互不覆盖 `run_spec.json` 与产物。
* `report_ensemble(dir)` → 读 `run/conf_*/molecule_opt/results.json`，给出 gap、
  最大力、电子数与构象能量的表，以及 gap/构象能量的极差（散布）。

口径：MMFF94s 只负责"哪些构象值得算"，不是能量判据；终判必须看 DFT 结果。
"""

from __future__ import annotations

import json
import shlex
from pathlib import Path
from typing import Any, Sequence

from . import conformers
from .junction.geometry import write_xyz

#: `ensemble.json` 里固定写出的口径说明。
CAVEAT = (
    "MMFF94s/ETKDGv3 只用于生成与去重候选构象，不是电子结构判据；"
    "每个构象各跑一次下游链，用 DFT 的 gap / 锚几何 / 接触距离比较。"
    "构象能量差在 kcal/mol 量级时应给出 Boltzmann 权重，而不是只取最低能。"
)


def _quote(path: Path) -> str:
    return shlex.quote(str(path))


def write_ensemble(smiles: str, outdir: str | Path, *, n_conformers: int = 20,
                   keep: int = 3, seed: int = conformers.DEFAULT_SEED,
                   rmsd_tol: float = conformers.DEFAULT_RMSD_TOL,
                   run_root: str = "run", config: str = "smj.yaml",
                   stage: str = "molecule_opt") -> dict[str, Any]:
    """生成 top-k 构象、批处理脚本与清单。返回清单字典。"""
    selected = conformers.select_conformers(
        smiles, n_conformers=n_conformers, seed=seed, rmsd_tol=rmsd_tol, keep=keep)
    if not selected:  # select_conformers 保证非空，这里只是防御
        raise ValueError("没有可用的构象（keep=%d）" % keep)

    root = Path(outdir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    run_dir_root = root / run_root
    config_arg = str(config)  # 原样写进脚本：默认 smj.yaml 相对于脚本执行目录

    entries: list[dict[str, Any]] = []
    commands: list[str] = []
    for item in selected:
        xyz = root / ("conf_%02d.xyz" % item.rank)
        write_xyz(xyz, item.atoms,
                  "conformer %d, MMFF94s %.3f kcal/mol, %s"
                  % (item.rank, item.energy_kcal, smiles))
        run_dir = run_dir_root / ("conf_%02d" % item.rank)
        entries.append({
            "rank": item.rank,
            "xyz": str(xyz),
            "run_dir": str(run_dir),
            "energy_kcal": item.energy_kcal,
            "rmsd_to_first_A": item.rmsd_to_first_A,
            "atom_count": item.atom_count,
        })
        # 脚本里的路径相对 ensemble 目录，便于整目录拷到集群再跑。
        commands.append(
            "smj stage %s --run --xyz %s --run-dir %s --config %s"
            % (stage, _quote(Path("conf_%02d.xyz" % item.rank)),
               _quote("%s/conf_%02d" % (run_root, item.rank)), _quote(config_arg)))

    script = root / "run_ensemble.sh"
    body = ["#!/usr/bin/env sh", "# 由 smj ensemble 生成：每个构象一个隔离的 run 目录。",
            "# 在脚本所在目录执行（路径都相对此处）；可用 SMJ_CONFIG 覆盖配置。",
            "set -eu", 'cd "$(dirname "$0")"']
    if not Path(config_arg).exists():
        body.append("# 注意：%s 不存在，先运行 `smj init`，或在集群上把 --config 指到真实配置"
                    % config_arg)
    body.extend(commands)
    body.append("")
    body.append("# 全部跑完后汇总：smj ensemble --dir . --report")
    script.write_text("\n".join(body) + "\n", encoding="utf-8", newline="\n")

    payload = {
        "generator": "smj.ensemble",
        "smiles": smiles,
        "seed": seed,
        "requested_conformers": n_conformers,
        "rmsd_tolerance_A": rmsd_tol,
        "kept": len(entries),
        "energy_unit": "kcal/mol (MMFF94s)",
        "stage": stage,
        "run_root": str(run_dir_root),
        "config": config_arg,
        "caveat": CAVEAT,
        "conformers": entries,
        "script": str(script),
        "commands": commands,
    }
    manifest = root / "ensemble.json"
    manifest.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
                        encoding="utf-8", newline="\n")
    return {"manifest": str(manifest), "script": str(script),
            "conformers": entries, "data": payload}


def _load_results(run_dir: Path) -> dict[str, Any] | None:
    path = run_dir / "molecule_opt" / "results.json"
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def report_ensemble(outdir: str | Path) -> dict[str, Any]:
    """汇总各构象 run 的 molecule_opt 结果，给出散布。"""
    root = Path(outdir)
    manifest_path = root / "ensemble.json"
    if not manifest_path.exists():
        raise FileNotFoundError("缺 %s：先跑 smj ensemble 生成清单" % manifest_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    rows: list[dict[str, Any]] = []
    missing: list[str] = []
    for entry in manifest.get("conformers", []):
        data = _load_results(Path(entry["run_dir"]))
        if data is None:
            missing.append(entry["xyz"])
            continue
        section = data.get("molecule_opt") or {}
        orbitals = section.get("frontier_orbitals") or {}
        rows.append({
            "rank": entry["rank"],
            "energy_kcal": entry.get("energy_kcal"),
            "atom_count": section.get("atom_count"),
            "converged": section.get("converged"),
            "max_force_ev_per_A": section.get("max_force_ev_per_A"),
            "gap_ev": orbitals.get("gap_ev"),
            "electron_count": section.get("electron_count"),
            "tier": (section.get("species") or {}).get("actual_tier"),
        })

    def spread(key: str) -> float | None:
        values = [row[key] for row in rows if isinstance(row.get(key), (int, float))]
        return (max(values) - min(values)) if len(values) >= 2 else None

    return {
        "dir": str(root),
        "smiles": manifest.get("smiles"),
        "kept": manifest.get("kept"),
        "rows": rows,
        "missing": missing,
        "spread": {"gap_ev": spread("gap_ev"),
                   "max_force_ev_per_A": spread("max_force_ev_per_A")},
        "caveat": manifest.get("caveat", CAVEAT),
    }


__all__ = ["CAVEAT", "write_ensemble", "report_ensemble"]
