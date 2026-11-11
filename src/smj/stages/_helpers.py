"""阶段间共用的小工具——锚定原子解析、资源记录、dry-run 提示命令。

单独成模块而不是让六个阶段互相 import：阶段之间**只允许通过文件耦合**，
这一点对"任意阶段可以单独 prepare"是硬前提。
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

from .._vendor.junction.geometry import Atom, read_aims_geometry, read_xyz
from ..errors import StageBlocked

#: 写进 ``StageState.resources`` 的估算字段（与 resources.estimate_smj_resources 对齐）
EST_KEYS = ("walltime", "time_limit", "mem", "ntasks", "cpus", "cores",
            "atoms", "core_h", "capped", "basis", "local", "unknown")


def est_record(estimate: dict[str, Any]) -> dict[str, Any]:
    """只保留估算结果里可直接序列化/展示的字段（去掉内部说明字段以外的杂项）。"""
    return {key: estimate[key] for key in EST_KEYS if key in estimate}


def submit_command(script: str | None) -> list[str]:
    """dry-run 时给用户的手工提交命令（真正提交由 runner 负责，阶段自己不提交）。"""
    return ["sbatch %s" % script] if script else []


def vdw_mode(enabled: bool) -> str | None:
    """配置里的 ``electrode.vdw: true`` -> FHI-aims 的 ``vdw_correction_hirshfeld``。"""
    return "hirshfeld" if enabled else None


def _as_index(value: Any) -> int | None:
    text = str(value).strip()
    if text.isdigit():
        return int(text) or None
    return None


def _resolve_anchor(atoms: list[Atom], value: Any, label: str, *, stage: str) -> int:
    """把单个锚定配置（元素符号或 1-based 序号）解析成 1-based 原子序号。"""
    total = len(atoms)
    index = _as_index(value)
    if index is not None:
        if not 1 <= index <= total:
            raise StageBlocked(
                stage, "%s 锚定索引越界: %s（结构共 %d 原子）" % (label, value, total))
        return index
    symbol = str(value).strip()
    matches = [i + 1 for i, atom in enumerate(atoms) if atom.element == symbol]
    if not matches:
        available = ", ".join(sorted({atom.element for atom in atoms}))
        raise StageBlocked(
            stage,
            "%s 锚定元素 %r 在结构中找不到（现有元素: %s）" % (label, symbol, available))
    if len(matches) != 2:
        raise StageBlocked(
            stage,
            "%s 锚定元素 %r 出现 %d 次，需恰好 2 个（1-based 序号: %s）"
            % (label, symbol, len(matches), matches),
        )
    return matches[0] if label == "left" else matches[1]


def resolve_anchor_indices_atoms(atoms: Iterable[Atom], left_anchor: Any, right_anchor: Any,
                                 *, stage: str = "aims_tip_relax") -> tuple[int, int]:
    """在已解析的原子列表上解析左右锚定原子（1-based）。"""
    atoms = list(atoms)
    left = _resolve_anchor(atoms, left_anchor, "left", stage=stage)
    right = _resolve_anchor(atoms, right_anchor, "right", stage=stage)
    if left == right:
        raise StageBlocked(stage, "左右锚定原子解析成同一个序号 %d，无法建结" % left)
    return left, right


def resolve_anchor_indices(xyz_path: str | Path, left_anchor: Any, right_anchor: Any,
                           *, stage: str = "aims_tip_relax") -> tuple[int, int]:
    """在 Gaussian 的 ``opt.xyz`` 上解析左右锚定原子（1-based）。

    配置风格是 ``electrode.left_anchor: S``（两个 S 锚原子，按文件出现顺序分左右）；
    也接受 ``11`` 这类显式 1-based 序号，便于只拿到别人输入的阶段独立运行。
    """
    return resolve_anchor_indices_atoms(
        read_xyz(xyz_path), left_anchor, right_anchor, stage=stage)


def resolve_geometry_anchor_indices(geometry_path: str | Path, left_anchor: Any,
                                    right_anchor: Any, *,
                                    stage: str = "aims_tip_relax") -> tuple[int, int]:
    """在 FHI-aims ``geometry.in`` 上解析左右锚定原子（1-based，供 collect 使用）。"""
    return resolve_anchor_indices_atoms(
        read_aims_geometry(geometry_path), left_anchor, right_anchor, stage=stage)
