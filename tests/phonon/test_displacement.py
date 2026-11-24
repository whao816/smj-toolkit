"""``apply_displacement.py`` 与 ``make_displacements`` 的符号约定。

这个文件的存在理由是一个真实的、代价很高的缺陷:早期 ``displace()`` 写成
``sign * delta_A``,而 ``make_displacements`` 写出的 ``delta_A`` 本身已带符号,
于是 (-1) * (-0.01) = +0.01 —— 正负两个位移都朝正向走,中心差分
``F+ - F-`` 恒为 0,力常数整矩阵变成 0。2026-11-24 在 gpu1 上用 54 个 SCF
的代价发现,而当时没有任何测试直接覆盖 ``apply_displacement``。

因此这里不只测"位移了多少",还测**符号配对之后能否恢复已知的力常数**:
后者才是真正会被这个缺陷打破的性质。
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from smj.phonon import apply_displacement as ap
from smj.phonon.fc import assemble_fc, make_displacements


GEOMETRY = """\
lattice_vector      5.0000000000       0.0000000000       0.0000000000
lattice_vector      2.5000000000       4.3301270189       0.0000000000
lattice_vector      0.0000000000       0.0000000000      20.0000000000
atom        0.00000000        0.00000000        0.00000000  Au
atom        2.88450000        0.00000000        0.00000000  S
atom        0.00000000        2.88450000        0.00000000  C
atom        1.44225000        1.44225000        2.35454990  H
"""


def _coords(text):
    out = []
    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 5 and parts[0] == "atom":
            out.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return np.array(out)


def test_plus_and_minus_entries_move_in_opposite_directions():
    """回归:``+`` 与 ``-`` 必须真的反向。"""
    disps = make_displacements(4, [0], 0.01)
    plus = [d for d in disps if d["sign"] > 0][0]
    minus = [d for d in disps if d["sign"] < 0][0]

    up, _ = ap.displace(GEOMETRY, plus)
    um, _ = ap.displace(GEOMETRY, minus)
    delta = _coords(up)[plus["atom"]] - _coords(um)[minus["atom"]]

    # 两个几何差恰好 2*amplitude,且方向与 axis 一致
    assert np.allclose(np.delete(delta, plus["axis"]), 0.0, atol=0.0)
    assert delta[plus["axis"]] == pytest.approx(0.02, abs=1e-12)


@pytest.mark.parametrize("slot", [0, 1, 2, 3, 4, 5])
def test_offset_matches_manifest_delta_exactly(slot):
    """位移量必须逐位等于 manifest 的 ``delta_A``(带符号)。"""
    disps = make_displacements(4, [0], 0.01)
    entry = disps[slot]
    moved, n_atoms = ap.displace(GEOMETRY, entry)
    assert n_atoms == 4

    before = _coords(GEOMETRY)
    after = _coords(moved)
    offset = after[entry["atom"]] - before[entry["atom"]]
    assert offset[entry["axis"]] == pytest.approx(entry["delta_A"], abs=1e-12)
    # 其余原子不动
    others = [i for i in range(4) if i != entry["atom"]]
    assert np.array_equal(after[others], before[others])


def test_lattice_and_constraint_lines_are_preserved():
    text = GEOMETRY + "constrain_relaxation .true.\n"
    moved, _ = ap.displace(text, make_displacements(4, [0], 0.01)[0])
    assert moved.count("lattice_vector") == 3
    assert "constrain_relaxation .true." in moved


def test_unknown_atom_index_raises():
    entry = {"slot": 0, "atom": 9, "axis": 0, "sign": 1, "delta_A": 0.01}
    with pytest.raises(ValueError, match="找不到"):
        ap.displace(GEOMETRY, entry)


def test_bad_axis_raises():
    entry = {"slot": 0, "atom": 0, "axis": 3, "sign": 1, "delta_A": 0.01}
    with pytest.raises(ValueError, match="axis"):
        ap.displace(GEOMETRY, entry)


def test_force_constants_are_recovered_from_signed_displacements():
    """端到端:已知 ``K`` → 用 manifest 造位移 → 谐波力 → 组装回 ``K``。

    这是会被"符号乘两次"直接打破的性质:两个位移都朝正向时
    ``F+ - F-`` 为 0,组装出的矩阵整片为 0,这里会以满量级偏差失败。
    """
    n_atoms = 4
    rng = np.random.default_rng(20261124)
    raw = rng.normal(size=(3 * n_atoms, 3 * n_atoms))
    k_true = 0.5 * (raw + raw.T)  # 对称,量级 ~1 eV/Ang^2

    device = list(range(n_atoms))
    disps = make_displacements(n_atoms, device, 0.01)

    forces = []
    for entry in disps:
        u = np.zeros(3 * n_atoms)
        u[3 * entry["atom"] + entry["axis"]] = entry["delta_A"]
        forces.append((-k_true @ u).reshape(n_atoms, 3))

    k_fit = assemble_fc(disps, forces, n_atoms)

    assert k_fit.shape == (3 * n_atoms, 3 * n_atoms)
    # 谐波模型下中心差分是精确的
    assert np.allclose(k_fit, k_true, atol=1e-10)


def test_cli_writes_displaced_geometry(workdir):
    """``main()`` 必须把结果落到 ``--out``,且 ``--index`` 真的被用上。"""
    root = Path(workdir)
    geometry = root / "geometry.in"
    geometry.write_text(GEOMETRY, encoding="utf-8")
    manifest = root / "displacements.json"
    manifest.write_text(
        json.dumps(make_displacements(4, [0], 0.01)), encoding="utf-8"
    )

    for index, expected_sign in ((0, +1), (1, -1)):
        out = root / ("out_%d.in" % index)
        rc = ap.main(
            [
                "--geometry", str(geometry),
                "--manifest", str(manifest),
                "--index", str(index),
                "--out", str(out),
            ]
        )
        assert rc == 0
        offset = _coords(out.read_text(encoding="utf-8"))[0][0]
        assert offset == pytest.approx(expected_sign * 0.01, abs=1e-12)


def test_cli_rejects_out_of_range_index(workdir):
    root = Path(workdir)
    geometry = root / "geometry.in"
    geometry.write_text(GEOMETRY, encoding="utf-8")
    manifest = root / "displacements.json"
    manifest.write_text(
        json.dumps(make_displacements(4, [0], 0.01)), encoding="utf-8"
    )
    rc = ap.main(
        [
            "--geometry", str(geometry),
            "--manifest", str(manifest),
            "--index", "99",
            "--out", str(root / "out.in"),
        ]
    )
    assert rc == 1


def test_module_stays_stdlib_only():
    """脚本要随作业上传到计算节点,不能依赖 numpy / smj 包。"""
    source = Path(ap.__file__).read_text(encoding="utf-8")
    assert "import numpy" not in source
    assert "from smj" not in source
    assert "import smj" not in source



