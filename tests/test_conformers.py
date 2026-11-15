"""构象粗筛前端：ETKDGv3 + MMFF94s + 重原子 RMSD 去重。

这些用例只在装了 RDKit 的环境里跑（``pytest.importorskip``），因为
``smj`` 的核心链不依赖 RDKit——它只是 ``molecule_opt`` 之前的可选前端。
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

pytest.importorskip("rdkit")

from smj.conformers import (  # noqa: E402  (importorskip 之后才能 import)
    DEFAULT_RMSD_TOL,
    DEFAULT_SEED,
    select_conformers,
    write_conformers,
)

CHAIN = "CCCCCCCCCCCCCCCCS"   # 正十六硫醇，51 原子（C16H34S），柔性链
BDT = "Sc1ccc(S)cc1"          # 苯二硫醇，14 原子，刚性


def test_bdt_has_a_single_conformer():
    confs = select_conformers(BDT, n_conformers=10, keep=3)
    assert len(confs) == 1
    assert confs[0].rank == 1
    assert confs[0].atom_count == 14
    assert confs[0].rmsd_to_first_A == pytest.approx(0.0)


def test_flexible_chain_keeps_distinct_conformers_sorted_by_energy():
    confs = select_conformers(CHAIN, n_conformers=8, keep=5)
    assert len(confs) >= 2
    energies = [c.energy_kcal for c in confs]
    assert energies == sorted(energies)
    assert [c.rank for c in confs] == list(range(1, len(confs) + 1))
    assert all(c.atom_count == 51 for c in confs)
    # 去重后除第一个外，都应至少偏离"最低能构象"一个容差
    assert confs[0].rmsd_to_first_A == pytest.approx(0.0)
    assert all(c.rmsd_to_first_A >= DEFAULT_RMSD_TOL for c in confs[1:])


def test_selection_is_deterministic_for_a_fixed_seed():
    a = select_conformers(CHAIN, n_conformers=6, seed=DEFAULT_SEED)
    b = select_conformers(CHAIN, n_conformers=6, seed=DEFAULT_SEED)
    assert [(c.rank, round(c.energy_kcal, 6)) for c in a] == \
        [(c.rank, round(c.energy_kcal, 6)) for c in b]


def test_seed_zero_is_rejected():
    """RDKit randomSeed=0 会让所有构象用同一初猜（实测两两 RMS=0），直接拒绝。"""
    with pytest.raises(ValueError, match="seed"):
        select_conformers(BDT, seed=0)


def test_invalid_arguments_are_rejected():
    with pytest.raises(ValueError):
        select_conformers(BDT, n_conformers=0)
    with pytest.raises(ValueError):
        select_conformers(BDT, keep=0)
    with pytest.raises(ValueError):
        select_conformers(BDT, rmsd_tol=0.0)
    with pytest.raises(ValueError, match="SMILES"):
        select_conformers("这不是 SMILES")


def test_write_conformers_writes_xyz_manifest_and_primary(workdir):
    confs = select_conformers(BDT, n_conformers=6, keep=3)
    written = write_conformers(confs, workdir, smiles=BDT, n_conformers=6)

    out = workdir / "molecule.xyz"
    manifest = workdir / "conformers.json"
    assert out.exists() and manifest.exists()
    assert (workdir / "conformer_01.xyz").exists()
    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "14"

    data = json.loads(manifest.read_text(encoding="utf-8"))
    assert data["kept"] == len(confs)
    assert data["smiles"] == BDT
    assert data["seed"] == DEFAULT_SEED
    assert data["rmsd_tolerance_A"] == pytest.approx(DEFAULT_RMSD_TOL)
    assert data["energy_unit"].startswith("kcal/mol")
    assert "MMFF94s" in data["caveat"]
    assert len(data["conformers"]) == len(confs)
    assert Path(written["out"]) == out


def test_write_conformers_rejects_empty_selection(workdir):
    with pytest.raises(ValueError, match="没有可写的构象"):
        write_conformers([], workdir)

def test_cli_conformers_writes_rigid_molecule(workdir, capsys):
    from smj.cli import main

    assert main(["conformers", "--smiles", BDT, "--dir", str(workdir),
                 "--n", "6", "--keep", "3"]) == 0
    assert (workdir / "molecule.xyz").exists()
    assert (workdir / "conformers.json").exists()
    assert "去重保留 1 个" in capsys.readouterr().out


def test_cli_conformers_json_output_and_bad_smiles(workdir, capsys):
    from smj.cli import EXIT_FAILED, main

    assert main(["conformers", "--smiles", "not-a-smiles", "--dir", str(workdir)]) == EXIT_FAILED
    capsys.readouterr()

    assert main(["conformers", "--smiles", BDT, "--dir", str(workdir),
                 "--n", "6", "--json"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert data["kept"] == 1
    assert data["smiles"] == BDT