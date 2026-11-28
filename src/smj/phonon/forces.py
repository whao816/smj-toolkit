"""从 ``aims.out`` 里读 FHI-aims 的逐原子力块。

为什么单独一个模块(方案 §5.3):

* FHI-aims 只在 control.in 里写了 ``compute_forces .true.`` 时才打印力块,
  且**每次 SCF / 每步弛豫各打印一块**。有限位移的每个目录都只有一块,但弛豫
  阶段会有几十块,必须取**最后**一块 —— 取第一块等于用未收敛的力算力常数。
* 数字是 Fortran 的 ``D`` 指数格式,``float()`` 直接解析会抛 ValueError。
* 力块只有 ``| 序号 fx fy fz``,没有元素名;元素顺序来自同一目录的 geometry.in,
  所以这里只返回数组,由调用方配对。

行格式(实测 ``runs/nap-sh/nap-sh/aims_tip_relax/aims.out``)::

      Total atomic forces (unitary forces cleaned) [eV/Ang]:
      |    1         -0.262905381321659E+00          0.954633163814336E+00   ...
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

__all__ = ["ForceBlockError", "count_force_blocks", "read_forces"]

#: 力块的标题行。FHI-aims 会加 "(unitary forces cleaned)" 之类的后缀,
#: 也可能换成 "(not cleaned)",所以只匹配前缀。
_TITLE = "Total atomic forces"
_ROW_RE = re.compile(
    r"^\s*\|\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$")


class ForceBlockError(ValueError):
    """aims.out 里没有可用的力块(通常是 control.in 少了 compute_forces)。"""


def _to_float(token: str, line_no: int) -> float:
    try:
        return float(token.replace("D", "E").replace("d", "e"))
    except ValueError:
        raise ForceBlockError(
            "第 %d 行的力分量无法解析: %r" % (line_no, token)) from None


def count_force_blocks(text: str) -> int:
    """力块个数(弛豫任务会有很多块,有限位移任务应该恰好 1 块)。"""
    return sum(1 for line in text.splitlines() if _TITLE in line)


def _parse_block(lines: list[str], start: int, line_offset: int) -> np.ndarray:
    rows: list[tuple[int, float, float, float]] = []
    index = start + 1
    while index < len(lines):
        line = lines[index]
        match = _ROW_RE.match(line)
        if match is None:
            # 力块后面可能是空行或下一段输出;遇到非空且不是数字行的就收工
            if line.strip():
                break
            index += 1
            continue
        rows.append((
            int(match.group(1)),
            _to_float(match.group(2), line_offset + index + 1),
            _to_float(match.group(3), line_offset + index + 1),
            _to_float(match.group(4), line_offset + index + 1),
        ))
        index += 1
    if not rows:
        raise ForceBlockError(
            "找到力块标题但没有解析到任何力行(标题在第 %d 行)"
            % (line_offset + start + 1))
    numbers = [row[0] for row in rows]
    expected = list(range(1, len(rows) + 1))
    if numbers != expected:
        raise ForceBlockError(
            "力块的原子序号不连续: 期望 1..%d,实际首尾 %r..%r"
            % (len(rows), numbers[0], numbers[-1]))
    return np.asarray([row[1:] for row in rows], dtype=float)


def read_forces(text_or_path: str | Path, *, n_atoms: int | None = None) -> np.ndarray:
    """返回**最后**一个力块,形状 ``(N, 3)``,单位 eV/Å。

    参数
    ----
    text_or_path:
        ``aims.out`` 路径,或已读入的文本(含换行即可)。
    n_atoms:
        期望原子数。给了就强制校验 —— 位矢和力揉错体系是这类流水线最贵的错误。

    异常
    ----
    ForceBlockError
        没有力块 / 序号不连续 / 数字不可解析 / 原子数不符。一律 fail-closed:
        **不要**返回全零数组当作"力为零"。
    """
    if isinstance(text_or_path, Path) or (
            isinstance(text_or_path, str) and "\n" not in text_or_path):
        path = Path(text_or_path)
        if not path.exists():
            raise FileNotFoundError("找不到 aims.out: %s" % path)
        text = path.read_text(encoding="utf-8", errors="replace")
    else:
        text = str(text_or_path)

    lines = text.splitlines()
    starts = [i for i, line in enumerate(lines) if _TITLE in line]
    if not starts:
        raise ForceBlockError(
            "aims.out 里没有 'Total atomic forces' 段:control.in 是否漏了 "
            "`compute_forces .true.`?(注意 `calculate_forces` 不是关键字)")
    forces = _parse_block(lines, starts[-1], line_offset=0)
    if n_atoms is not None and forces.shape[0] != int(n_atoms):
        raise ForceBlockError(
            "力块里 %d 个原子,期望 %d 个" % (forces.shape[0], int(n_atoms)))
    if not np.all(np.isfinite(forces)):
        raise ForceBlockError("力块含非有限数值(NaN/inf)")
    return forces
