"""单分子结谱学解析:KS 本征值匹配 + Gaussian log 解析(优化结构/SCF 能量/频率)。

来源: 第三方工具 qmtransport/orbitals.py(47 行)+ qmtransport/gaussian.py 解析部分
(tmp/自动化/, 只读来源不删除不入库)按本项目协议适配迁入(AGD-2026-09-22-SMJ, P2)。
修改(§4.3):原子序数表从原工具 14 元素子集(缺 Li/Na/Mg/Al/Si 回退 X{编号})改为
**完整内置周期表映射**, 不依赖 pymatgen/ase optional 依赖。
纯本地工具: 只做解析, 不做真实计算。
"""
from __future__ import annotations

from pathlib import Path
import json
import re
from .geometry import Atom, write_xyz
from .landauer import nearest_energy, read_te


# ---- 完整周期表(1-92, 含常见 f 区; 超出报错而非 X{编号} 静默回退) ----
ATOMIC_SYMBOLS: dict[int, str] = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne",
    11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar",
    19: "K", 20: "Ca", 21: "Sc", 22: "Ti", 23: "V", 24: "Cr", 25: "Mn", 26: "Fe",
    27: "Co", 28: "Ni", 29: "Cu", 30: "Zn", 31: "Ga", 32: "Ge", 33: "As", 34: "Se",
    35: "Br", 36: "Kr", 37: "Rb", 38: "Sr", 39: "Y", 40: "Zr", 41: "Nb", 42: "Mo",
    43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn",
    51: "Sb", 52: "Te", 53: "I", 54: "Xe", 55: "Cs", 56: "Ba",
    57: "La", 58: "Ce", 59: "Pr", 60: "Nd", 61: "Pm", 62: "Sm", 63: "Eu", 64: "Gd",
    65: "Tb", 66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb", 71: "Lu",
    72: "Hf", 73: "Ta", 74: "W", 75: "Re", 76: "Os", 77: "Ir", 78: "Pt", 79: "Au",
    80: "Hg", 81: "Tl", 82: "Pb", 83: "Bi", 84: "Po", 85: "At", 86: "Rn",
    87: "Fr", 88: "Ra", 89: "Ac", 90: "Th", 91: "Pa", 92: "U",
}


def _symbol(atomic_number: int) -> str:
    sym = ATOMIC_SYMBOLS.get(atomic_number)
    if sym is None:
        raise ValueError(f"Unknown atomic number {atomic_number} (周期表 1-92 之外)")
    return sym


def parse_eigenvalues(out_path: str | Path) -> list[dict]:
    rows = []
    pattern = re.compile(r"^\s*(\d+)\s+([-+]?\d+(?:\.\d+)?)\s+([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)\s+([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)")
    for line in Path(out_path).read_text(encoding="utf-8", errors="ignore").splitlines():
        match = pattern.match(line)
        if match:
            state, occ, ha, ev = match.groups()
            rows.append({"state": int(state), "occupation": float(occ), "eigenvalue_Ha": float(ha), "eigenvalue_eV": float(ev)})
    if not rows:
        raise ValueError(f"No Kohn-Sham eigenvalue rows found in {out_path}")
    return rows


def match_orbital(
    te_path: str | Path,
    aims_out: str | Path,
    mode: str = "fermi",
    target_e_ev: float | None = None,
    output: str | Path | None = None,
) -> dict:
    te_rows = read_te(te_path)
    if mode == "fermi":
        target = nearest_energy(te_rows, 0.0)
    elif mode == "peak":
        r = max(te_rows, key=lambda row: row[2])
        target = {"E_Hartree": r[0], "E_minus_Ef_eV": r[1], "T": r[2]}
    elif mode == "energy":
        if target_e_ev is None:
            raise ValueError("--target-e-ev is required when --mode energy")
        target = nearest_energy(te_rows, target_e_ev)
    else:
        raise ValueError("mode must be fermi, peak, or energy")
    eigs = parse_eigenvalues(aims_out)
    best = min(eigs, key=lambda row: abs(row["eigenvalue_Ha"] - target["E_Hartree"]))
    result = {"target_transmission_point": target, "matched_state": best}
    if output:
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        Path(output).write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def parse_standard_orientation(text: str) -> list[Atom]:
    blocks: list[list[Atom]] = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        if "Standard orientation:" not in lines[i]:
            i += 1
            continue
        i += 1
        dash_count = 0
        while i < len(lines):
            if set(lines[i].strip()) == {"-"}:
                dash_count += 1
                i += 1
                if dash_count == 2:
                    break
                continue
            i += 1
        rows: list[Atom] = []
        while i < len(lines):
            line = lines[i]
            if set(line.strip()) == {"-"}:
                break
            fields = line.split()
            if len(fields) >= 6 and fields[0].isdigit():
                atomic_number = int(fields[1])
                element = _symbol(atomic_number)
                rows.append(Atom(element, float(fields[3]), float(fields[4]), float(fields[5])))
            i += 1
        if rows:
            blocks.append(rows)
        i += 1
    if not blocks:
        raise ValueError("No Standard orientation block found")
    return blocks[-1]


def parse_gaussian_log(log_path: str | Path, xyz_path: str | Path | None = None, summary_path: str | Path | None = None) -> dict:
    text = Path(log_path).read_text(encoding="utf-8", errors="ignore")
    scf = re.findall(r"SCF Done:\s+E\([^)]+\)\s+=\s+([-+]?\d+\.\d+)", text)
    freq_values = []
    for line in re.findall(r"Frequencies --\s+(.+)", text):
        freq_values.extend(float(x) for x in line.split())
    imaginary = [f for f in freq_values if f < 0]
    normal_termination = "Normal termination" in text
    atoms = parse_standard_orientation(text)
    # 2026-11-20 复核 P1：SCF 是否真的收敛、开壳层的自旋纯度。
    # Gaussian 在每次 SCF/稳定性分析后打印 "S**2 before annihilation X, after Y"，
    # 取**最后一段**（最终解）；没有该行说明是闭壳层单重态或没做自旋分析。
    s2_pairs = re.findall(
        r"S\*\*2\s+before annihilation\s+([-+]?\d+\.\d+)\s*,\s*after\s+([-+]?\d+\.\d+)",
        text)
    s2_before = float(s2_pairs[-1][0]) if s2_pairs else None
    s2_after = float(s2_pairs[-1][1]) if s2_pairs else None
    converged = bool(scf) and "Convergence failure" not in text
    summary = {
        "log": str(log_path),
        "normal_termination": normal_termination,
        "scf_converged": converged,
        "scf_done_blocks": len(scf),
        "final_scf_energy_hartree": float(scf[-1]) if scf else None,
        "imaginary_frequency_count": len(imaginary),
        "imaginary_frequencies_cm-1": imaginary,
        "atom_count": len(atoms),
        "s2_before_annihilation": s2_before,
        "s2_after_annihilation": s2_after,
    }
    if xyz_path:
        write_xyz(xyz_path, atoms, comment=f"extracted from {log_path}")
        summary["xyz"] = str(xyz_path)
    if summary_path:
        Path(summary_path).parent.mkdir(parents=True, exist_ok=True)
        Path(summary_path).write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary