"""开壳层结（自由基/阳离子）的自旋口径与 control.in / geometry.in 后处理。

配置里只有 ``molecule.multiplicity`` 与 ``molecule.charge`` 两个物理量；FHI-aims
需要它们以 ``spin`` / ``fixed_spin_moment`` / ``charge`` 三个关键字出现。vendored
的 ``make_control`` 已支持 spin/charge（relax 模板还留了 ``spin_moment_line``
占位），但结级调用方此前没传——自由基从 tip 弛豫起就静默退化成闭壳层中性结。

本模块把"配置 → 关键字"集中在唯一一处，并补一条 FHI-aims 手册要求、但 vendored
模板没有的关键字：

* ``default_initial_moment``：``spin collinear`` 时 FHI-aims 需要初始自旋密度
  才起得来（手册原文：default_initial_moment / initial_moment 至少给一个）。

带电体系不能用 ``override_initial_charge_check`` 绕检查：项目使用的 FHI-aims
221103 不认这个关键字（2026-11-19 集群实测 parse_control 直接停算），而新版
手册要求 ``charge`` 等于 geometry.in 里 ``initial_charge`` 之和。因此这里反过来
写 per-atom ``initial_charge``：把总电荷均分到所有原子，让新旧版本的检查都满足。
它只是初始密度猜测，不影响 SCF 自洽结果。
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

#: 开壳层结的初始磁矩默认值（μB/原子）。净磁矩约束由 fixed_spin_moment 单独控制。
DEFAULT_INITIAL_MOMENT = 0.1


@dataclass(frozen=True)
class SpinState:
    """从配置翻译出的结级自旋状态（FHI-aims 口径）。"""

    multiplicity: int
    charge: int
    #: 是否约束 N↑−N↓=M−1。默认普通开壳层固定，BS 释放。
    fix_spin_moment: bool | None = None
    #: M=1 的破缺对称初猜使用共线自旋和反号局域磁矩。
    broken_symmetry: bool = False

    def __post_init__(self) -> None:
        if self.broken_symmetry and self.multiplicity > 1:
            raise ValueError(
                "broken_symmetry 只适用于 multiplicity=1 的开壳层单重态，"
                "当前 multiplicity=%d" % self.multiplicity)
        if self.fix_spin_moment is None:
            object.__setattr__(self, "fix_spin_moment", not self.broken_symmetry)

    @property
    def open_shell(self) -> bool:
        return self.multiplicity > 1 or self.broken_symmetry

    @property
    def spin(self) -> str:
        """FHI-aims ``spin`` 只接受 ``none`` / ``collinear``（不接受数字）。"""
        return "collinear" if self.open_shell else "none"

    @property
    def fixed_spin_moment(self) -> float:
        """固定时 N↑−N↓ 的目标值；BS 初猜对应 0，允许非零局域自旋密度。"""
        return float(self.multiplicity - 1)


def spin_state(config: Any) -> SpinState:
    """读 ``config.molecule`` 的 charge/multiplicity；缺 molecule 时按闭壳层中性。"""
    molecule = getattr(config, "molecule", None)
    multiplicity = int(getattr(molecule, "multiplicity", 1) or 1)
    charge = int(getattr(molecule, "charge", 0) or 0)
    fix = getattr(molecule, "fixed_spin_moment", None)
    broken = bool(getattr(molecule, "broken_symmetry", False))
    return SpinState(multiplicity=multiplicity, charge=charge,
                     fix_spin_moment=fix, broken_symmetry=broken)


def describe_spin_constraint(moment: float | None) -> str:
    """说明共线自旋的净磁矩约束；0 与释放有不同的物理含义。"""
    if moment is None:
        return "释放净磁矩约束，N↑−N↓ 由 SCF 自洽决定"
    return "fixed_spin_moment %g（固定 N↑−N↓，局域磁矩可自洽变化）" % moment


def spin_control_lines(state: SpinState, *, with_fixed_spin_moment: bool,
                       default_initial_moment: float | None = DEFAULT_INITIAL_MOMENT,
                       ) -> list[str]:
    """需要补进 control.in general 段的关键字行。

    ``with_fixed_spin_moment`` 只对 tip_relax 有意义的 ``control_relax.in`` 为真：
    relax 模板只有 ``spin_moment_line`` 占位（默认空），scf 模板则自带
    ``fixed_spin_moment {fixed_spin_moment}``，由 ``make_control`` 的参数写入。

    ``default_initial_moment=None`` 时不写全局初猜——用于 ``molecule.initial_moments``
    已经逐原子给出初猜的情形（FHI-aims 要求两者至少给一个）。
    """
    lines: list[str] = []
    if state.open_shell:
        if with_fixed_spin_moment and state.fix_spin_moment:
            lines.append("  fixed_spin_moment  %g" % state.fixed_spin_moment)
        if default_initial_moment is not None:
            lines.append("  default_initial_moment  %g" % float(default_initial_moment))
    return lines


def remove_control_lines(control: Path, keys: Sequence[str]) -> list[str]:
    """删除 control.in general 段里首个 token 命中 ``keys`` 的行；返回删除的行。

    scf 模板总是写 ``fixed_spin_moment {值}``；``molecule.fixed_spin_moment=false``
    时没有"不写"的模板参数，只能在写完后把这一行拿掉。
    """
    path = Path(control)
    wanted = set(keys)
    lines = path.read_text(encoding="utf-8").splitlines()
    removed = [line for line in lines
               if line.split() and line.split()[0] in wanted]
    if not removed:
        return []
    kept = [line for line in lines
            if not (line.split() and line.split()[0] in wanted)]
    path.write_text("\n".join(kept) + "\n", encoding="utf-8", newline="\n")
    return removed


def insert_control_lines(control: Path, lines: Sequence[str]) -> list[str]:
    """把缺少的关键字行插到第一个 ``species`` 行之前；返回实际插入的行。

    FHI-aims 的 general 段必须位于 species 块之前（之后会被当成 species 子 tag）。
    已经出现过的关键字跳过，所以对同一份 control.in 重复调用是幂等的。
    """
    path = Path(control)
    if not lines:
        return []
    existing = path.read_text(encoding="utf-8").splitlines()
    present = {line.split()[0] for line in existing if line.split()}
    inserted = [line for line in lines if line.split()[0] not in present]
    if not inserted:
        return []
    insert_at = len(existing)
    for position, line in enumerate(existing):
        if line.split()[:1] == ["species"]:
            insert_at = position
            break
    merged = existing[:insert_at] + inserted + existing[insert_at:]
    path.write_text("\n".join(merged) + "\n", encoding="utf-8", newline="\n")
    return inserted


def write_initial_charges(geometry: Path, charge: float) -> float | None:
    """在 geometry.in 每个 ``atom`` 行后插入 per-atom ``initial_charge``。

    FHI-aims 要求 ``charge`` 等于 geometry.in 里 ``initial_charge`` 之和（不写时
    q 必须为 0）；本包不知道分子内初始电荷该怎么分，就按原子数均分。返回每原子
    写入值；已存在 ``initial_charge`` 或没有 atom 行时返回 ``None``（不改文件）。
    """
    path = Path(geometry)
    if not charge:
        return None
    lines = path.read_text(encoding="utf-8").splitlines()
    if any(line.split()[:1] == ["initial_charge"] for line in lines):
        return None
    atoms = {i for i, line in enumerate(lines) if line.split()[:1] == ["atom"]}
    if not atoms:
        return None
    per_atom = float(charge) / len(atoms)
    merged: list[str] = []
    for i, line in enumerate(lines):
        merged.append(line)
        if i in atoms:
            merged.append("  initial_charge  %.10f" % per_atom)
    path.write_text("\n".join(merged) + "\n", encoding="utf-8", newline="\n")
    return per_atom


def map_molecule_moments(moments: dict[int, float], removed_h_indices: Sequence[int] = (),
                         *, molecule_atoms: int | None = None) -> dict[int, float]:
    """把分子坐标系（1-based）的初始磁矩映射到结几何的原子序号。

    tip_relax 只在 ``strip_thiol_h`` 触发时移除锚定 S–H 的氢，其余原子顺序不变；
    结几何（molecule + 2 tips + Au 团簇）里分子原子仍在最前面，所以映射就是
    "减掉它前面被删掉的氢"。``moments`` 指向被删掉的氢时报错（该原子已不存在）。
    """
    removed = sorted({int(index) for index in removed_h_indices or ()})
    mapped: dict[int, float] = {}
    for index, value in moments.items():
        index = int(index)
        if index in removed:
            raise ValueError(
                "molecule.initial_moments 指向原子 %d，但该原子在去质子时被移除" % index)
        new_index = index - sum(1 for removed_index in removed if removed_index < index)
        if molecule_atoms is not None and new_index > molecule_atoms:
            raise ValueError(
                "molecule.initial_moments 原子序号 %d 超出分子原子数 %d"
                % (index, molecule_atoms))
        mapped[new_index] = float(value)
    return dict(sorted(mapped.items()))


def write_initial_moments(geometry: Path, moments: dict[int, float]) -> dict[int, float]:
    """在 geometry.in 指定 ``atom`` 行后插入 per-atom ``initial_moment``。

    序号是**该 geometry 自身**的 1-based 原子序号（结几何的分子前缀）。
    已存在 ``initial_moment`` 行时原样返回（幂等），不做部分覆盖。
    """
    path = Path(geometry)
    if not moments:
        return {}
    lines = path.read_text(encoding="utf-8").splitlines()
    if any(line.split()[:1] == ["initial_moment"] for line in lines):
        return {}
    merged: list[str] = []
    written: dict[int, float] = {}
    atom_index = 0
    for line in lines:
        merged.append(line)
        if line.split()[:1] != ["atom"]:
            continue
        atom_index += 1
        if atom_index in moments:
            merged.append("  initial_moment  %g" % float(moments[atom_index]))
            written[atom_index] = float(moments[atom_index])
    missing = sorted(set(moments) - set(written))
    if missing:
        raise ValueError(
            "geometry.in 只有 %d 个 atom，initial_moment 序号越界: %s"
            % (atom_index, ", ".join(str(index) for index in missing)))
    path.write_text("\n".join(merged) + "\n", encoding="utf-8", newline="\n")
    return written


def _update_summary(summary_path: Path, values: dict[str, Any]) -> None:
    try:
        data = json.loads(Path(summary_path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return
    if not isinstance(data, dict):
        return
    data.update(values)
    Path(summary_path).write_text(
        json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")


def record_junction_spin(summary_path: Path, state: SpinState, *,
                         default_initial_moment: float | None = DEFAULT_INITIAL_MOMENT,
                         initial_moments: dict[int, float] | None = None) -> None:
    """把结级自旋设置补进 ``input_summary.json``（truth chain 可核对）。"""
    values: dict[str, Any] = {"spin": state.spin, "charge": float(state.charge)}
    if state.broken_symmetry:
        # BS 单重态在输入摘要里必须显式标出来：multiplicity=1 看起来像闭壳层，
        # 只看 multiplicity 的读者会把两套自旋通道误当成误配。
        values["broken_symmetry"] = True
    if state.open_shell:
        if state.fix_spin_moment:
            values["fixed_spin_moment"] = state.fixed_spin_moment
        else:
            values["fixed_spin_moment"] = None
            values["spin_constraint"] = "released"
        if default_initial_moment is not None:
            values["default_initial_moment"] = float(default_initial_moment)
        if initial_moments:
            values["initial_moments"] = {str(k): v for k, v in initial_moments.items()}
    _update_summary(summary_path, values)


def record_initial_charge(summary_path: Path, charge: float, per_atom: float) -> None:
    """记录写入 geometry.in 的初始电荷（均分值 + 总量）。"""
    _update_summary(summary_path, {
        "initial_charge_total": float(charge),
        "initial_charge_per_atom": per_atom,
    })
