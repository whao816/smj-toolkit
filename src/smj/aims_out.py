"""FHI-aims ``aims.out`` 的文本解析：几何收敛与前沿轨道。

纯正则扫文本，不读其它文件、不依赖 numpy；给 :mod:`smj.stages.molecule_opt`
的 collect 提供两类事实：

* :func:`relax_summary` —— 几何弛豫是否收敛、最后一步的最大力分量（eV/Å）、
  弛豫步数；
* :func:`frontier_orbitals` —— 最后一段 Kohn-Sham 特征值表的 HOMO/LUMO/gap（eV）。

口径（2026-11-15 实测 ``demo/2CN-656-4-thiomethylphenyl-aimsopt`` 的
``gaussian_opt/aims.out``）：

* 未收敛几何的每一步都会写 ``Present geometry is not yet converged.``，收敛时
  写 ``Present geometry is converged.``。判据用完整串 ``Present geometry is
  converged``——它不可能是 ``... not yet converged`` 的子串，不会误匹配。
* ``Maximum force component is  0.928904E-02 eV/A.`` 是**当步**的最大力，取
  日志里最后一个匹配即"最后一步"。
* ``relax_geometry bfgs 1.e-2`` 的 1.e-2 单位是 **eV/Å**（日志：
  ``Convergence accuracy for geometry relaxation: Maximum force < 0.100000E-01 eV/A.``），
  不是 Hartree/Bohr——旧注释与旧文档写错了，见 ``docs/smj-计算流程梳理.md``。
* 特征值表每步打印一次，表头 ``State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]``，
  故取**最后一张表**（最后一步几何的电子结构）。
"""

from __future__ import annotations

import re
from typing import Any

__all__ = ["RELAX_CONVERGED", "frontier_orbitals", "relax_summary", "scf_settings"]

RELAX_CONVERGED = "Present geometry is converged"
RELAX_NOT_CONVERGED = "Present geometry is not yet converged"
SCF_CONVERGED = "Self-consistency cycle converged"
#: 同一行/相邻行里带 not 的写法（不同版本措辞不完全一致，宁可宽一点）。
_SCF_FAIL_RE = re.compile(
    r"self-consistency cycle\b[^.\n]{0,40}\bnot\b[^.\n]{0,40}converged", re.IGNORECASE)

#   Maximum force component is  0.928904E-02 eV/A.
_FORCE_RE = re.compile(
    r"Maximum force component is\s+([-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)\s+eV/A")
_STEP_RE = re.compile(r"Relaxation step number\s+(\d+)")
_EIG_HEADER = "Eigenvalue [eV]"
_EIG_ROW_RE = re.compile(
    r"^\s*(\d+)\s+([-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)\s+"
    r"([-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)\s+"
    r"([-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)\s*$"
)

#: 占据数 ≥ 该值视为占据态；< 该值视为空态。PBE 分子计算里占据数是 0/1/2。
_OCCUPIED_MIN = 0.5
#: 0.05 < occ < 0.95 视为分数占据（金属性 / 费米面穿过能带）。
_FRACTIONAL_LOW = 0.05
_FRACTIONAL_HIGH = 0.95


def _to_float(token: str) -> float:
    """把 ``1.0E-02`` / ``1.0D-02`` 统一成 float（Fortran 的 D 指数）。"""
    return float(token.replace("D", "E").replace("d", "e"))


#: 直接从 ``control.in`` 读回、写进 results 的 SCF 精度键（顺序即输出顺序）。
SCF_KEYS = ("occupation_type", "mixer", "n_max_pulay", "charge_mix_param",
            "sc_accuracy_rho", "sc_accuracy_eev", "sc_accuracy_etot",
            "sc_accuracy_forces", "sc_iter_limit")


def scf_settings(control_text: str) -> dict[str, str]:
    """从 ``control.in`` 读回实际生效的 SCF 精度设置。

    由 control.in 本身解析而不是复制模板常量：旧 run（没有这些行的）会返回空
    dict，调用方据此判断"这步力没有显式精度保证"，而不是假装有。
    """
    found: dict[str, str] = {}
    for raw in control_text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        parts = line.split()
        if parts[0] in SCF_KEYS and len(parts) > 1:
            found[parts[0]] = " ".join(parts[1:])
    return found


def relax_summary(text: str) -> dict[str, Any]:
    """几何弛豫摘要。

    Returns:
        ``{"converged": bool, "max_force_ev_per_A": float | None,
        "relax_steps": int, "scf_converged_steps": int,
        "scf_not_converged": bool}``；``max_force_ev_per_A`` 取日志中最后一个
        ``Maximum force component`` 行，缺失时为 ``None``（不伪造 0.0）。

    SCF 每步写一行 ``Self-consistency cycle converged.``，所以
    ``scf_converged_steps`` 应 ≈ ``relax_steps``（demo 实测 106 / 105，首次 SCF
    在第一步之前就跑了）。若日志里出现 "not ... converged" 的 SCF 行，
    ``scf_not_converged`` 为真——那一步的力没有自洽保证。
    """
    forces = [_to_float(tok) for tok in _FORCE_RE.findall(text)]
    steps = [int(tok) for tok in _STEP_RE.findall(text)]
    return {
        "converged": RELAX_CONVERGED in text,
        "max_force_ev_per_A": forces[-1] if forces else None,
        "relax_steps": max(steps) if steps else 0,
        "scf_converged_steps": text.count(SCF_CONVERGED),
        "scf_not_converged": bool(_SCF_FAIL_RE.search(text)),
    }


def frontier_orbitals(text: str) -> dict[str, Any] | None:
    """最后一段 Kohn-Sham 特征值表给出的 HOMO/LUMO/gap（eV）。

    分子没有分数占据时 HOMO 是最高占据态、LUMO 是高于 HOMO 的最低空态。
    分数占据（金属性）时仍给这两个数，但置 ``fractional_occupation=True``，
    调用方不应把它当"整数占据的 HOMO-LUMO 带隙"来解读。

    Returns:
        ``None``（日志里没有特征值表，如任务在 SCF 前失败）或
        ``{"homo_ev", "lumo_ev", "gap_ev", "n_states", "fractional_occupation"}``。
    """
    tables = _eigen_tables(text)
    if not tables:
        return None
    rows = tables[-1]
    if not rows:
        return None
    occupied = [ev for _, occ, ev in rows if occ >= _OCCUPIED_MIN]
    if not occupied:
        return None
    homo_ev = max(occupied)
    above = [ev for _, occ, ev in rows
             if occ < _OCCUPIED_MIN and ev >= homo_ev]
    lumo_ev = min(above) if above else None
    return {
        "homo_ev": round(homo_ev, 6),
        "lumo_ev": None if lumo_ev is None else round(lumo_ev, 6),
        "gap_ev": None if lumo_ev is None else round(lumo_ev - homo_ev, 6),
        "n_states": len(rows),
        "fractional_occupation": any(
            _FRACTIONAL_LOW < occ < _FRACTIONAL_HIGH for _, occ, _ in rows),
    }


def _eigen_tables(text: str) -> list[list[tuple[int, float, float]]]:
    """切出所有 ``State/Occupation/Eigenvalue`` 表，返回 ``[(state, occ, eV), ...]``。"""
    tables: list[list[tuple[int, float, float]]] = []
    current: list[tuple[int, float, float]] | None = None
    for line in text.splitlines():
        if _EIG_HEADER in line and "State" in line:
            if current is not None:
                tables.append(current)
            current = []
            continue
        if current is None:
            continue
        match = _EIG_ROW_RE.match(line)
        if match:
            state, occ, _ha, ev = match.groups()
            current.append((int(state), _to_float(occ), _to_float(ev)))
        elif current and line.strip() and not line.startswith((" ", "\t")):
            # 表结束（下一段叙述文本的顶格行）；留空行不算结束。
            tables.append(current)
            current = None
    if current is not None:
        tables.append(current)
    return [table for table in tables if table]
