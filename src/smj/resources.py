"""作业资源估算（墙钟 / 内存 / rank / 核时）。

数值与 AgenticDFT 引擎保持同一套（2026-11 的真实标定，见
``agenticdft/agenticdft/tools/resources.py`` 的注释）：

- SCF / Aitranss 是重阶段，按 132 原子 BDT 结的**端到端实测峰值**做幂律外推：
  wall ∝ N^3，mem ∝ N^2（SCF）/ N^1.5（TE）。
- Gaussian / tip-relax 是轻阶段，线性外推更稳。
- 墙钟/内存各留 1.25x / 1.20x 余量；``--time`` 再乘 1.9（SCF 收敛后写产物
  还要 ~25 min，132 原子实测曾因 `--time` 卡在产物写完之前被 Slurm 杀掉）。
- 上下限 72 h / 128 GiB。
"""

from __future__ import annotations

import math
import re

WALLTIME_CAP_S = 72 * 3600
MEM_CAP_GIB = 128.0

#: 重阶段：幂律锚点（n0 处的实测墙钟 / 峰值内存）
ANCHORS: dict[str, dict[str, float]] = {
    "aims_scf": {"n0": 132, "wall0": 4320.0, "peak0": 13.3, "wall_exp": 3.0, "mem_exp": 2.0},
    "aitranss_te": {"n0": 132, "wall0": 2552.0, "peak0": 8.0, "wall_exp": 3.0, "mem_exp": 1.5},
}

#: 轻阶段：实测锚点 + 每原子增量
LINEAR: dict[str, dict[str, float]] = {
    "gaussian_opt": {"wall0": 60.0, "per_atom_wall": 10.0, "mem0": 4.0, "per_atom_mem": 0.05},
    "aims_tip_relax": {"wall0": 180.0, "per_atom_wall": 30.0, "mem0": 8.0, "per_atom_mem": 0.10},
    # au_relax：背层冻结后只弛豫中间区，原子数约为 au60 全体系的 1/3 自由度，
    # 但每步力都要算全体系，墙钟按"整体系的一次 BFGS"估（比 SCF 略贵）。
    "au_relax": {"wall0": 300.0, "per_atom_wall": 40.0, "mem0": 10.0, "per_atom_mem": 0.10},
}

NTASKS: dict[str, int] = {
    "gaussian_opt": 1,
    "aims_tip_relax": 4,
    "au_relax": 4,
    "aims_scf": 4,
    "aitranss_te": 1,
}

#: Aitranss 在这个构建里不并行（2026-11-05 独占窗口 A/B：-c 16 反而慢 10%），
#: 只给 1 核；其它阶段用 MPI rank 表达并行度。
CPUS: dict[str, int] = {"aitranss_te": 1}

MARGIN_WALL = 1.25
MARGIN_MEM = 1.20
TIME_LIMIT_FACTOR = 1.9

FLOOR_WALL: dict[str, int] = {
    "gaussian_opt": 900,
    "aims_tip_relax": 900,
    "au_relax": 1800,
    "aims_scf": 1800,
    "aitranss_te": 1800,
}

SMJ_STAGES: tuple[str, ...] = (
    "gaussian_opt", "aims_tip_relax", "au_relax", "aims_scf", "aitranss_te",
)


def hhmmss(seconds: float) -> str:
    """秒 -> ``HH:MM:SS``，向上取整到分钟。"""
    secs = int(math.ceil(float(seconds) / 60.0) * 60)
    hh, rem = divmod(secs, 3600)
    mm, ss = divmod(rem, 60)
    return "%02d:%02d:%02d" % (hh, mm, ss)


def walltime_seconds(spec: str) -> int:
    """``HH:MM:SS`` / ``MM:SS`` -> 秒（非法输入返回 0）。"""
    try:
        parts = [int(x) for x in str(spec or "").split(":") if str(x).strip()]
    except ValueError:
        return 0
    if len(parts) == 3:
        return parts[0] * 3600 + parts[1] * 60 + parts[2]
    if len(parts) == 2:
        return parts[0] * 3600 + parts[1] * 60
    return 0


def estimate_smj_resources(stage: str, natoms: int = 1) -> dict:
    """smj 阶段资源建议。

    Returns
    -------
    dict
        ``{walltime, time_limit, mem, ntasks, cpus, cores, atoms, core_h, capped, basis}``；
        本地阶段返回 ``{'local': True, ...}``；未知阶段返回 ``{'unknown': True, ...}``。
    """
    if stage == "au60_build":
        return {"local": True, "basis": "au60_build 为纯本地阶段，无集群资源"}
    n = max(int(natoms or 1), 1)
    if stage in ANCHORS:
        a = ANCHORS[stage]
        ratio = n / float(a["n0"])
        wall = a["wall0"] * (ratio ** a["wall_exp"]) * MARGIN_WALL
        mem = a["peak0"] * (ratio ** a["mem_exp"]) * MARGIN_MEM
        basis = (
            "实测锚点 %d 原子 %.0f s / %.1f GiB；墙钟 x N^%.1f、内存 x N^%.1f 外推，余量 %.2fx/%.2fx"
            % (a["n0"], a["wall0"], a["peak0"], a["wall_exp"], a["mem_exp"], MARGIN_WALL, MARGIN_MEM)
        )
    elif stage in LINEAR:
        p = LINEAR[stage]
        wall = (p["wall0"] + p["per_atom_wall"] * n) * MARGIN_WALL
        mem = (p["mem0"] + p["per_atom_mem"] * n) * MARGIN_MEM
        basis = (
            "实测锚点线性外推：墙钟 %.0f+%.0f·N s，内存 %.1f+%.2f·N GiB，余量 %.2fx/%.2fx"
            % (p["wall0"], p["per_atom_wall"], p["mem0"], p["per_atom_mem"], MARGIN_WALL, MARGIN_MEM)
        )
    else:
        return {"unknown": True, "basis": "未知 smj 阶段: %s（可用: %s）" % (stage, list(SMJ_STAGES))}

    wall = max(wall, float(FLOOR_WALL.get(stage, 0)))
    capped: list[str] = []
    if wall > WALLTIME_CAP_S:
        wall = float(WALLTIME_CAP_S)
        capped.append("walltime->72h")
    wall_limit = min(wall * TIME_LIMIT_FACTOR, float(WALLTIME_CAP_S))
    if mem > MEM_CAP_GIB:
        mem = MEM_CAP_GIB
        capped.append("mem->128GiB")
    ntasks = NTASKS.get(stage, 1)
    cpus = CPUS.get(stage)
    cores = int(cpus or ntasks or 1)
    return {
        "walltime": hhmmss(wall),
        "time_limit": hhmmss(wall_limit),
        "mem": "%dG" % int(math.ceil(mem)),
        "ntasks": ntasks,
        "cpus": cpus,
        "cores": cores,
        "atoms": n,
        "core_h": round(cores * wall / 3600.0, 2),
        "capped": capped,
        "basis": "FHI-aims 内存随 rank 放大，rank 取 %d；%s；--time 取预期墙钟 x%.1f"
                 % (ntasks, basis, TIME_LIMIT_FACTOR),
    }


def check_smj_budget(est: dict, *, walltime_cap_s: int = WALLTIME_CAP_S,
                     mem_cap_gib: float = MEM_CAP_GIB) -> dict:
    """估算是否落在预算内：``{'ok': bool, 'reason': str}``。"""
    if not est:
        return {"ok": False, "reason": "无资源估算"}
    if est.get("local"):
        return {"ok": True, "reason": "本地阶段，不计集群预算"}
    if est.get("unknown"):
        return {"ok": False, "reason": est.get("basis", "未知阶段")}
    secs = walltime_seconds(est.get("walltime") or "")
    mem = int(re.sub(r"[^0-9]", "", str(est.get("mem") or "0")) or 0)
    over = []
    if secs > walltime_cap_s:
        over.append("walltime %s > %s" % (est.get("walltime"), hhmmss(walltime_cap_s)))
    if mem > mem_cap_gib:
        over.append("mem %s > %dG" % (est.get("mem"), int(mem_cap_gib)))
    if over:
        return {"ok": False, "reason": "; ".join(over)}
    return {"ok": True, "reason": "在预算内"}
