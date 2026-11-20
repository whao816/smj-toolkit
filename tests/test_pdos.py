"""分组 PDOS（``src/smj/pdos.py``）的离线测试。

夹具按 FHI-aims 的真实形状造：``basis-indices.out`` 的 ``# fn. type at. n l m``
行、``mos.aims`` 的 ``<序号> <a|b> eigenvalue=… nsaos=…`` 头行 + 每行 4 个系数。
不依赖任何远端产物，也不碰 ``tmp_path``（沙箱下会 PermissionError）。
"""

from __future__ import annotations

import json
import math

import pytest

from smj.pdos import (
    DEFAULT_SIGMA_EV,
    HARTREE_TO_EV,
    PdosError,
    analyze_pdos,
    parse_basis_indices,
    parse_groups,
    parse_mos,
)

EFERMI_HA = -0.188452434085


def _basis(natoms: int, per_atom: int = 2) -> str:
    lines = ["$basis.indices.aims", "# amount of basis functions used:  %d" % (natoms * per_atom),
             "# fn.   type   at.   n   l   m"]
    fn = 0
    for atom in range(1, natoms + 1):
        for n in range(1, per_atom + 1):
            fn += 1
            lines.append("  %3d atomic      %d   %d   0   0" % (fn, atom, n))
    return "\n".join(lines) + "\n"


def _mos(state_energy_offsets_ev, natoms: int, per_atom: int = 2, spin: str = "a") -> str:
    """每个原子前一半系数给 alpha 组、后一半给 beta 组（方便手算权重）。"""
    nsaos = natoms * per_atom
    lines = ["$scfmo.aims   iteration:   12   format(4d20.14)"]
    for index, offset in enumerate(state_energy_offsets_ev, start=1):
        energy_ha = EFERMI_HA + offset / HARTREE_TO_EV
        lines.append("  %4d  %s      eigenvalue=%s   nsaos=  %d"
                     % (index, spin, ("%.14f" % energy_ha).replace("e", "D"), nsaos))
        coefficients = []
        for fn in range(1, nsaos + 1):
            # 只让每个原子的第一个基函数有系数：该原子在态里权重 = 1
            coefficients.append(math.sqrt(1.0 / natoms) if (fn - 1) % per_atom == 0 else 0.0)
        for start in range(0, nsaos, 4):
            chunk = coefficients[start:start + 4]
            lines.append("".join("%.14f" % value for value in chunk))
    return "\n".join(lines) + "\n"


def _te(grid_ev) -> str:
    lines = ["#bias = 0.000000 V efermi = %.12f H" % EFERMI_HA, "# energy[eV] T"]
    for energy in grid_ev:
        ha = EFERMI_HA + energy / HARTREE_TO_EV
        lines.append("%.10f %.10f %.10f" % (ha, energy, 0.5))
    return "\n".join(lines) + "\n"


def _run(workdir, *, natoms: int = 4, grid=(-0.5, -0.25, 0.0, 0.25, 0.5),
         offsets=(0.0, 0.25), mos: str | None = None):
    run = workdir / "run"
    scf = run / "aims_scf"
    te_dir = run / "aitranss_te"
    scf.mkdir(parents=True)
    te_dir.mkdir(parents=True)
    (scf / "basis-indices.out").write_text(_basis(natoms), encoding="utf-8")
    (scf / "mos.aims").write_text(mos if mos is not None else _mos(offsets, natoms), encoding="utf-8")
    (te_dir / "TE.dat").write_text(_te(grid), encoding="utf-8")
    return run


def test_parse_basis_indices_maps_functions_to_atoms(workdir):
    path = workdir / "basis-indices.out"
    path.write_text(_basis(3), encoding="utf-8")
    parsed = parse_basis_indices(path)
    assert parsed["nsaos"] == 6
    assert parsed["atom_of_fn"] == [1, 1, 2, 2, 3, 3]


def test_parse_basis_indices_rejects_broken_numbering(workdir):
    path = workdir / "basis-indices.out"
    path.write_text(_basis(2).replace("    2 atomic", "    5 atomic"), encoding="utf-8")
    with pytest.raises(PdosError, match="编号不连续"):
        parse_basis_indices(path)


def test_parse_mos_reads_energies_and_coefficients(workdir):
    path = workdir / "mos.aims"
    path.write_text(_mos((0.0, 0.5), natoms=2), encoding="utf-8")
    parsed = parse_mos(path)
    assert parsed["nsaos"] == 4
    assert [state["index"] for state in parsed["states"]] == [1, 2]
    first = parsed["states"][0]
    assert first["eigenvalue_Ha"] == pytest.approx(EFERMI_HA, abs=1e-12)
    assert len(first["coefficients"]) == 4
    assert first["coefficients"][0] == pytest.approx(math.sqrt(0.5), abs=1e-12)


def test_parse_mos_explains_truncated_file(workdir):
    """本地只回取了 head 的典型现场：头行在、系数行没有。"""
    path = workdir / "mos.aims"
    path.write_text(
        "$scfmo.aims   iteration:   90   format(4d20.14)\n"
        "     1  a      eigenvalue=-.18845243408500D+00   nsaos=  4364\n",
        encoding="utf-8")
    with pytest.raises(PdosError) as excinfo:
        parse_mos(path)
    message = str(excinfo.value)
    assert "0/4364" in message
    assert "head" in message


def test_parse_groups_accepts_ranges_and_rejects_overlap():
    groups = parse_groups("left=1-3,right:4-6", natoms=6)
    assert groups == [{"name": "left", "atoms": [1, 2, 3]},
                      {"name": "right", "atoms": [4, 5, 6]}]
    with pytest.raises(PdosError, match="同时出现在"):
        parse_groups("a=1-4,b=4-6", natoms=6)
    with pytest.raises(PdosError, match="超出"):
        parse_groups("a=1-9", natoms=6)


def test_analyze_pdos_groups_sum_to_one_state(workdir):
    """每个原子在态里权重 1/natoms，按原子分组后各组权重之和应为 1。"""
    grid = tuple(round(-0.9 + 0.05 * i, 6) for i in range(37))
    run = _run(workdir, natoms=4, offsets=(0.0,), grid=grid)
    result = analyze_pdos(run, groups="a=1-2,b=3-4", sigma_ev=0.3)
    channel = result["channels"]["mos"]
    weights = channel["integrated_weight"]
    assert weights["a"] == pytest.approx(0.5, abs=1e-9)
    assert weights["b"] == pytest.approx(0.5, abs=1e-9)
    # 展宽归一化：数值积分 ≈ 总权重（σ=0.3 eV、步长 0.05 eV）
    area = sum(channel["pdos"]["a"]) * (result["grid_eV"][1] - result["grid_eV"][0])
    assert area == pytest.approx(0.5, rel=0.01)


def test_analyze_pdos_aligns_state_at_fermi_to_zero(workdir):
    run = _run(workdir, natoms=2, offsets=(0.0,), grid=(-0.2, -0.1, 0.0, 0.1, 0.2))
    result = analyze_pdos(run, groups="mol=1-2", sigma_ev=0.05)
    peak_index = max(range(len(result["grid_eV"])),
                     key=lambda i: result["channels"]["mos"]["pdos"]["mol"][i])
    assert result["grid_eV"][peak_index] == pytest.approx(0.0, abs=1e-9)
    assert result["efermi_Ha"] == pytest.approx(EFERMI_HA, abs=1e-12)


def test_analyze_pdos_writes_json_and_sorts_grid(workdir):
    run = _run(workdir, natoms=2, grid=(0.3, -0.3, 0.0))
    target = workdir / "out" / "pdos.json"
    result = analyze_pdos(run, groups="mol=1-2", output=target)
    assert result["grid_eV"] == sorted(result["grid_eV"])
    payload = json.loads(target.read_text(encoding="utf-8"))
    assert payload["schema"] == "smj.pdos.v1"
    assert payload["sigma_eV"] == pytest.approx(DEFAULT_SIGMA_EV)
    assert payload["groups"][0]["functions"] == 4


def test_analyze_pdos_supports_spin_channels(workdir):
    run = _run(workdir, natoms=2, mos=_mos((0.0,), natoms=2))
    scf = run / "aims_scf"
    (scf / "alpha.aims").write_text(_mos((0.0,), natoms=2, spin="a"), encoding="utf-8")
    (scf / "beta.aims").write_text(_mos((0.1,), natoms=2, spin="b"), encoding="utf-8")
    result = analyze_pdos(run, groups="mol=1-2")
    assert set(result["channels"]) == {"alpha", "beta"}
    assert result["channels"]["beta"]["spin"] == "beta"
    assert result["channels"]["beta"]["state_energies_eV"][0] == pytest.approx(0.1, abs=1e-9)


def test_analyze_pdos_refuses_without_te(workdir):
    run = _run(workdir, natoms=2)
    (run / "aitranss_te" / "TE.dat").unlink()
    with pytest.raises(PdosError, match="找不到 TE.dat"):
        analyze_pdos(run, groups="mol=1-2")

def test_cli_pdos_prints_groups_and_channels(workdir, capsys):
    from smj.cli import main

    run = _run(workdir, natoms=2)
    rc = main(["pdos", str(run), "--groups", "mol=1-2", "--json"])
    assert rc == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["schema"] == "smj.pdos.v1"
    assert payload["groups"][0]["name"] == "mol"


def test_cli_pdos_requires_config_only_when_fetching(workdir, capsys):
    from smj.cli import main

    run = _run(workdir, natoms=2)
    rc = main(["pdos", str(run), "--groups", "mol=1-2", "--json"])
    assert rc == 0
    rc = main(["pdos", str(run), "--groups", "mol=1-2", "--fetch",
               "-c", str(workdir / "missing.yaml")])
    assert rc == 1
    assert "找不到配置文件" in capsys.readouterr().err
