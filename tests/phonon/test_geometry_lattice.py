"""``write_aims_geometry`` 的 ``lattice=`` 扩展 —— W0-1 逐字节兼容性守卫。

为什么这条测试值得单独一个文件(方案 §5.1 / devplan M2-T1):``geometry.py`` 是本仓
**唯一**需要改动的内核文件(``src/smj/junction/`` 由本仓维护,A 方案),而既有 25 条
smj 记录的锚点回归靠的就是 ``lattice=None`` 时的输出与 0.2.2 **逐字节相同**。
所以这里不但要测新功能,还要把旧输出的每一个空格钉死。
"""
from __future__ import annotations

import pathlib

import pytest

from smj.junction.geometry import Atom, write_aims_geometry

#: 0.2.2 的输出原样(含行内空格数量)。
_EXPECTED_NO_LATTICE = (
    "atom        0.00000000        0.00000000        0.00000000  Au\n"
    "constrain_relaxation .true.\n"
    "atom        0.00000000        0.00000000        2.35455000  S\n"
    "atom        1.40000000        0.00000000        3.90000000  C\n"
    "constrain_relaxation .true.\n"
)

_ATOMS = [
    Atom("Au", 0.0, 0.0, 0.0),
    Atom("S", 0.0, 0.0, 2.35455),
    Atom("C", 1.4, 0.0, 3.9),
]

_LATTICE = (
    (4.07823, 0.0, 0.0),
    (0.0, 4.07823, 0.0),
    (0.0, 0.0, 25.0),
)


def _write(workdir: pathlib.Path, **kwargs) -> str:
    target = workdir / "geometry.in"
    write_aims_geometry(target, _ATOMS, **kwargs)
    return target.read_text(encoding="utf-8")


def test_lattice_none_is_byte_identical(workdir):
    """``lattice=None`` 必须与 0.2.2 逐字节一致 —— 不是近似相等。"""
    assert _write(workdir, constrained=[1, 3]) == _EXPECTED_NO_LATTICE


def test_default_call_is_byte_identical_too(workdir):
    """连 ``constrained`` 都缺省时同样不能变。"""
    got = _write(workdir)
    assert got == (
        "atom        0.00000000        0.00000000        0.00000000  Au\n"
        "atom        0.00000000        0.00000000        2.35455000  S\n"
        "atom        1.40000000        0.00000000        3.90000000  C\n"
    )


def test_lattice_lines_come_first_with_ten_decimals(workdir):
    """三条 ``lattice_vector`` 在文件最前,精度 %.10f。"""
    got = _write(workdir, lattice=_LATTICE)
    lines = got.splitlines()
    assert lines[0] == "lattice_vector 4.0782300000 0.0000000000 0.0000000000"
    assert lines[1] == "lattice_vector 0.0000000000 4.0782300000 0.0000000000"
    assert lines[2] == "lattice_vector 0.0000000000 0.0000000000 25.0000000000"
    # 原子行紧接着三行晶格矢量
    assert lines[3].startswith("atom ")
    # 晶格矢量只出现在最前面
    assert sum(1 for line in lines if line.startswith("lattice_vector")) == 3


def test_constrained_order_unchanged_with_lattice(workdir):
    """加晶格后 ``constrain_relaxation`` 仍紧跟对应原子行(顺序不变)。"""
    got = _write(workdir, constrained=[2], lattice=_LATTICE)
    lines = got.splitlines()
    assert lines[3].startswith("atom ")          # Au
    assert lines[4].startswith("atom ")          # S(序号 2,冻结)
    assert lines[5] == "constrain_relaxation .true."
    assert lines[6].startswith("atom ")          # C
    assert len(lines) == 7


def test_lattice_rejects_bad_shapes(workdir):
    with pytest.raises(ValueError):
        _write(workdir, lattice=_LATTICE[:2])                       # 只有两条
    with pytest.raises(ValueError):
        _write(workdir, lattice=((1.0, 2.0), _LATTICE[1], _LATTICE[2]))  # 分量不足
    with pytest.raises(ValueError):
        _write(workdir, lattice=((float("nan"), 0.0, 0.0), _LATTICE[1], _LATTICE[2]))


def test_constrained_range_still_validated(workdir):
    with pytest.raises(ValueError):
        _write(workdir, constrained=[4])          # 只有 3 个原子
