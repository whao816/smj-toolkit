"""``read_forces``:力块解析。

夹具 ``tests/data/aims_force_block_20atoms.out`` 是真实 run 里最后一个力块的
逐字节摘录(nap-sh 的 20 原子 tip_relax,第 27 块,原文件第 80992 行起)。
``runs/`` 被 gitignore,``docs/validation-data`` 里也没有任何 ``aims.out``,
所以只能把它裁剪进仓库(devplan T5)。
"""
from __future__ import annotations

import pathlib

import numpy as np
import pytest

from smj.phonon.forces import ForceBlockError, count_force_blocks, read_forces

DATA = pathlib.Path(__file__).resolve().parent.parent / "data"
FIXTURE = DATA / "aims_force_block_20atoms.out"


def test_reads_last_block_of_real_output():
    forces = read_forces(FIXTURE, n_atoms=20)
    assert forces.shape == (20, 3)
    assert np.all(np.isfinite(forces))
    # 首行逐字对应夹具里的 | 1 -0.768937654526694E-02 ...
    np.testing.assert_allclose(
        forces[0], [-7.68937654526694e-03, 2.49110346551588e-03,
                    1.92976622842031e-03], rtol=0, atol=0)
    # 末行 | 20 0.668789888210689E-02 ...
    np.testing.assert_allclose(
        forces[-1], [6.68789888210689e-03, 1.54041691885551e-03,
                     6.13196276381952e-03], rtol=0, atol=0)


def test_counts_one_block_in_fixture():
    assert count_force_blocks(FIXTURE.read_text(encoding="utf-8")) == 1


def test_takes_the_last_block_when_several():
    text = (
        "  Total atomic forces (unitary forces cleaned) [eV/Ang]:\n"
        "  |    1   1.0  2.0  3.0\n"
        "  |    2   4.0  5.0  6.0\n"
        "\n  Total atomic forces (unitary forces cleaned) [eV/Ang]:\n"
        "  |    1   -1.0  -2.0  -3.0\n"
        "  |    2   -4.0  -5.0  -6.0\n"
    )
    assert count_force_blocks(text) == 2
    np.testing.assert_allclose(read_forces(text),
                               [[-1, -2, -3], [-4, -5, -6]])


def test_fortran_d_exponent_is_parsed():
    text = (
        "  Total atomic forces (unitary forces cleaned) [eV/Ang]:\n"
        "  |    1   1.5D-02  -2.5D+00  0.0D+00\n"
    )
    np.testing.assert_allclose(read_forces(text), [[0.015, -2.5, 0.0]])


def test_missing_block_names_the_right_keyword():
    with pytest.raises(ForceBlockError) as excinfo:
        read_forces("  Have a nice day.\n")
    message = str(excinfo.value)
    assert "compute_forces" in message
    assert "Total atomic forces" in message


def test_atom_count_mismatch_is_rejected():
    with pytest.raises(ForceBlockError):
        read_forces(FIXTURE, n_atoms=21)


def test_non_contiguous_indices_are_rejected():
    text = (
        "  Total atomic forces (unitary forces cleaned) [eV/Ang]:\n"
        "  |    1   0.0  0.0  0.0\n"
        "  |    3   0.0  0.0  0.0\n"
    )
    with pytest.raises(ForceBlockError):
        read_forces(text)


def test_non_finite_force_is_rejected():
    text = (
        "  Total atomic forces (unitary forces cleaned) [eV/Ang]:\n"
        "  |    1   NaN  0.0  0.0\n"
    )
    with pytest.raises(ForceBlockError):
        read_forces(text)


def test_title_without_rows_is_rejected():
    text = "  Total atomic forces (unitary forces cleaned) [eV/Ang]:\n\n"
    with pytest.raises(ForceBlockError):
        read_forces(text)


def test_missing_file_raises_file_not_found():
    with pytest.raises(FileNotFoundError):
        read_forces(DATA / "no-such-file.out")


# --- 远程 stdlib 脚本必须与本模块逐位一致(M2-T6) -------------------------


def test_extract_forces_matches_read_forces():
    from smj.phonon import extract_forces

    text = FIXTURE.read_text(encoding="utf-8")
    payload = extract_forces.extract(text)
    assert payload["n_atoms"] == 20
    assert payload["n_blocks"] == 1
    assert payload["block_index"] == 0
    # 逐位一致:两条解析路径不能有 1e-16 的分叉
    np.testing.assert_array_equal(
        np.asarray(payload["forces"]), read_forces(text))


def test_extract_forces_picks_last_block():
    from smj.phonon import extract_forces

    text = (
        "  Total atomic forces (unitary forces cleaned) [eV/Ang]:\n"
        "  |    1   1.0  2.0  3.0\n"
        "\n  Total atomic forces (unitary forces cleaned) [eV/Ang]:\n"
        "  |    1   -1.0  -2.0  -3.0\n"
    )
    payload = extract_forces.extract(text)
    assert payload["n_blocks"] == 2
    assert payload["block_index"] == 1
    assert payload["forces"] == [[-1.0, -2.0, -3.0]]


def test_extract_forces_cli(workdir):
    from smj.phonon import extract_forces

    target = workdir / "forces.json"
    assert extract_forces.main([str(FIXTURE), str(target)]) == 0
    payload = __import__("json").loads(target.read_text(encoding="utf-8"))
    np.testing.assert_array_equal(np.asarray(payload["forces"]),
                                  read_forces(FIXTURE))
    assert extract_forces.main(["--help"]) == 0


def test_extract_forces_cli_reports_missing_block(workdir):
    from smj.phonon import extract_forces

    source = workdir / "aims.out"
    source.write_text("  Have a nice day.\n", encoding="utf-8")
    assert extract_forces.main([str(source), str(workdir / "f.json")]) == 1


def test_extract_forces_uses_stdlib_only():
    from smj.phonon import extract_forces

    source = pathlib.Path(extract_forces.__file__).read_text(encoding="utf-8")
    body = "\n".join(line for line in source.splitlines()
                     if not line.lstrip().startswith("#"))
    for forbidden in ("import numpy", "from smj", "import smj"):
        assert forbidden not in body
