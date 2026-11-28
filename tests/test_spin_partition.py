"""S1：Hirshfeld 自旋矩解析与分区（2026-11-21 自由基支线改进计划）。

复核依据：``runs/ttf-cation`` 的 ``aims.out`` 最后一段 Hirshfeld 分析给出
88 个原子、总自旋矩 1.0000；按 1–14 分子 / 15–16 接触区 / 17–88 电极分区得
0.0056 / 0.0063 / 0.9881。旧实现只核验段存在（``aims_scf.py`` 的
``hirshfeld_mentions``），从不解析数值，因此"自旋留在分子上还是跑到电极上"
没有第一手证据。

这里锁住解析与分区的最小事実：
* 多段分析取最后一段（早期迭代表不能当终态）；
* 缺 ``Hirshfeld spin moment`` 行的原子容错为 ``None``，不被误当 0；
* 分区占比按组可加、和为 1；
* 接触区 = 距锚原子 ≤cutoff 的电极原子；
* 没有 Hirshfeld 段时 ``summarize`` 返回 ``None``（调用方据此跳过，而不是伪造 0）。
"""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from smj.spin_partition import (
    contact_atoms,
    default_groups,
    parse_hirshfeld,
    parse_hirshfeld_text,
    summarize,
)


def _section(
    rows: list[tuple[int, str, float | None]],
    *,
    header: str = "Performing Hirshfeld analysis of fragment charges and moments.",
) -> list[str]:
    """按真实 ``aims.out`` 版式造一段 Hirshfeld 分析。"""
    lines = ["  " + header, "  " + "-" * 70]
    for index, element, moment in rows:
        lines.append("  | Atom %5d: %s" % (index, element))
        lines.append("  |   Hirshfeld charge        :     -0.06377126")
        lines.append("  |   Free atom volume        :     37.22336203")
        lines.append("  |   Hirshfeld volume        :     31.59741609")
        if moment is not None:
            lines.append("  |   Hirshfeld spin moment   :      %.8f" % moment)
        lines.append("  |   Hirshfeld dipole moment :      0.05164383")
    return lines


def _write(workdir: Path, *sections: list[str]) -> Path:
    lines = ["  FHI-aims version      : 221103"]
    for block in sections:
        lines.extend(block)
        lines.append("  End of iteration")
    lines.append("  Have a nice day.")
    path = workdir / "aims.out"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_parse_takes_the_last_hirshfeld_section(workdir: Path):
    path = _write(workdir,
                  _section([(1, "S", 0.0), (2, "S", 0.0)]),
                  _section([(1, "S", 0.25), (2, "S", -0.25)]))
    table = parse_hirshfeld(path)
    assert table.analyses == 2
    assert [atom.index for atom in table.atoms] == [1, 2]
    assert [atom.element for atom in table.atoms] == ["S", "S"]
    assert [atom.spin_moment for atom in table.atoms] == pytest.approx([0.25, -0.25])
    assert table.total == pytest.approx(0.0)


def test_missing_moment_line_is_tolerated(workdir: Path):
    path = _write(workdir, _section([(1, "S", 0.5), (2, "S", None)]))
    table = parse_hirshfeld(path)
    assert table.atoms[0].spin_moment == pytest.approx(0.5)
    assert table.atoms[1].spin_moment is None
    assert table.missing == 1
    assert table.total == pytest.approx(0.5)


def test_no_section_returns_empty_table(workdir: Path):
    path = workdir / "aims.out"
    path.write_text("  FHI-aims\n  Have a nice day.\n", encoding="utf-8")
    table = parse_hirshfeld(path)
    assert table.atoms == []
    assert table.analyses == 0
    assert summarize(table, {}) is None


def test_summarize_partitions_by_group(workdir: Path):
    path = _write(workdir, _section([
        (1, "C", 0.60), (2, "C", 0.40), (3, "Au", 0.10), (4, "Au", -0.10),
    ]))
    table = parse_hirshfeld(path)
    summary = summarize(table, {"molecule": [1, 2], "electrode": [3, 4]})
    assert summary["total"] == pytest.approx(1.0)
    assert summary["atoms"] == 4
    assert summary["missing"] == 0
    molecule = summary["groups"]["molecule"]
    assert molecule["atoms"] == 2
    assert molecule["spin_moment"] == pytest.approx(1.0)
    assert molecule["fraction"] == pytest.approx(1.0)
    electrode = summary["groups"]["electrode"]
    assert electrode["spin_moment"] == pytest.approx(0.0)
    assert electrode["fraction"] == pytest.approx(0.0)
    assert electrode["sum_abs"] == pytest.approx(0.2)
    assert electrode["abs_fraction"] == pytest.approx(0.2 / 1.2)
    assert summary["n_positive"] == 3
    assert summary["n_negative"] == 1


def test_summarize_without_groups_reports_total_only(workdir: Path):
    table = parse_hirshfeld(_write(workdir, _section([(1, "S", 0.7), (2, "S", 0.3)])))
    summary = summarize(table, {})
    assert summary["total"] == pytest.approx(1.0)
    assert summary["groups"] == {}


def test_groups_must_cover_atoms_once(workdir: Path):
    table = parse_hirshfeld(_write(workdir, _section([(1, "S", 0.5), (2, "S", 0.5)])))
    with pytest.raises(ValueError, match="重复|越界|未覆盖"):
        summarize(table, {"molecule": [1, 2], "electrode": [2]})


def test_contact_atoms_are_electrode_atoms_near_anchor():
    geometry = [
        ("S", 0.0, 0.0, 0.0),    # 1 分子锚原子
        ("C", 3.0, 0.0, 0.0),    # 2 分子
        ("Au", 2.4, 0.0, 0.0),   # 3 电极，距锚 2.4 Å
        ("Au", 5.0, 0.0, 0.0),   # 4 电极，太远
        ("Au", 0.0, 2.9, 0.0),   # 5 电极，距锚 2.9 Å
    ]
    assert contact_atoms(geometry, molecule_atoms=2, anchors=(1,), cutoff=3.0) == [3, 5]
    # 锚原子自身不属于电极；越界的锚被忽略而不是抛错。
    assert contact_atoms(geometry, molecule_atoms=2, anchors=(1, 9), cutoff=3.0) == [3, 5]


def test_default_groups_split_contact_out_of_electrode():
    groups = default_groups(2, 5, contact=[3, 5])
    assert groups["molecule"] == [1, 2]
    assert groups["contact"] == [3, 5]
    assert groups["electrode"] == [4]
    plain = default_groups(2, 5)
    assert set(plain) == {"molecule", "electrode"}
    assert plain["electrode"] == [3, 4, 5]


def test_review_reports_spin_partition(workdir: Path):
    """复核里要能看到"自旋主要落在电极侧"，但只记 info、不阻断。"""
    from smj.review import ReviewResult, _spin_review

    run_dir = workdir / "run"
    (run_dir / "aims_scf").mkdir(parents=True)
    (run_dir / "aims_scf" / "spin_partition.json").write_text(json.dumps({
        "source": "aims.out",
        "total": 1.0,
        "atoms": 4,
        "atom_moments": [
            {"index": i, "spin_moment": m}
            for i, m in enumerate([0.0028, 0.0028, 0.4972, 0.4972], start=1)
        ],
        "groups": {
            "molecule": {"atoms": 2, "spin_moment": 0.0056, "fraction": 0.0056},
            "electrode": {"atoms": 2, "spin_moment": 0.9944, "fraction": 0.9944},
        },
    }), encoding="utf-8")
    (run_dir / "aitranss_te").mkdir()
    (run_dir / "aitranss_te" / "TE.dat").write_text(
        "#bias = 0.0 V efermi = -0.18 H\n"
        "-0.19 -0.2 0.20 0.19\n-0.18 0.0 0.25 0.24\n-0.17 0.2 0.20 0.19\n",
        encoding="utf-8")
    config = SimpleNamespace(molecule=SimpleNamespace(
        multiplicity=2, charge=0, fixed_spin_moment=True, broken_symmetry=False))
    result = ReviewResult(run_id="t")
    _spin_review(result, run_dir, config)
    assert "spin_localized_off_molecule" in result.codes()
    issue = next(i for i in result.issues if i["code"] == "spin_localized_off_molecule")
    assert issue["severity"] == "info"
    assert "0.56%" in issue["message"] and "0.0056" in issue["message"]


def test_summarize_flags_localized_antiparallel_moments(workdir: Path):
    """S4：破缺对称单重态的证据是"分子内正负自旋矩共存且净矩近零"。"""
    table = parse_hirshfeld(_write(workdir, _section([
        (1, "C", 0.25), (2, "C", -0.25), (3, "Au", 0.50), (4, "Au", 0.50),
    ])))
    summary = summarize(table, {"molecule": [1, 2], "electrode": [3, 4]})
    molecule = summary["groups"]["molecule"]
    assert molecule["n_positive"] == 1
    assert molecule["n_negative"] == 1
    assert molecule["net"] == pytest.approx(0.0)
    assert molecule["max_abs"] == pytest.approx(0.25)
    assert summary["localized_antiparallel"] is True
    assert summary["molecule_net"] == pytest.approx(0.0)
    assert summary["molecule_max_abs"] == pytest.approx(0.25)
    # 全表口径：净矩 1.0、最大绝对值 0.5，不能被分子组的"净零"覆盖。
    assert summary["net"] == pytest.approx(1.0)
    assert summary["max_abs"] == pytest.approx(0.5)


def test_summarize_does_not_flag_polarized_moments(workdir: Path):
    """分子内同号自旋矩是普通极化，不能报"局部反平行"。"""
    table = parse_hirshfeld(_write(workdir, _section([
        (1, "C", 0.30), (2, "C", 0.50), (3, "Au", 0.10), (4, "Au", 0.10),
    ])))
    summary = summarize(table, {"molecule": [1, 2], "electrode": [3, 4]})
    assert summary["localized_antiparallel"] is False
    assert summary["groups"]["molecule"]["n_negative"] == 0
    assert summary["molecule_net"] == pytest.approx(0.80)


def test_summarize_does_not_flag_numerical_noise(workdir: Path):
    """局部矩只有 1e-6 μB 时，净矩/幅度比再小也是噪声，不能报反平行。"""
    table = parse_hirshfeld(_write(workdir, _section([
        (1, "C", 1e-6), (2, "C", -1e-6), (3, "Au", 0.0), (4, "Au", 0.0),
    ])))
    summary = summarize(table, {"molecule": [1, 2], "electrode": [3, 4]})
    assert summary["groups"]["molecule"]["n_positive"] == 1
    assert summary["groups"]["molecule"]["n_negative"] == 1
    assert summary["localized_antiparallel"] is False


def test_render_shows_broken_symmetry_evidence():
    """BS run 的 §4.1 表下要有局部证据行；普通开壳层不出现该行。"""
    from smj.report import _render_spin_section

    partition = {
        "total": 1.0,
        "localized_antiparallel": True,
        "molecule_n_positive": 1,
        "molecule_n_negative": 1,
        "molecule_net": 0.0,
        "groups": {
            "molecule": {"atoms": 2, "spin_moment": 0.0, "fraction": 0.0,
                         "n_positive": 1, "n_negative": 1, "net": 0.0, "max_abs": 0.25},
            "electrode": {"atoms": 2, "spin_moment": 1.0, "fraction": 1.0},
        },
    }
    lines: list[str] = []
    _render_spin_section(lines, {"channels": 2, "broken_symmetry": True}, partition)
    text = "\n".join(lines)
    assert "破缺对称证据" in text
    assert "局部反平行" in text
    assert "⟨S²⟩" in text and "NOON" in text

    lines = []
    _render_spin_section(lines, {"channels": 2, "broken_symmetry": False}, partition)
    assert "破缺对称证据" not in "\n".join(lines)


def test_review_reports_broken_symmetry_local_moments(workdir: Path):
    """BS run 的复核要给局部自旋矩 info，且明确不是 ⟨S²⟩ 级别的证据。"""
    from smj.review import ReviewResult, _spin_review

    run_dir = workdir / "run"
    (run_dir / "aims_scf").mkdir(parents=True)
    (run_dir / "aims_scf" / "spin_partition.json").write_text(json.dumps({
        "source": "aims.out",
        "total": 1.0,
        "atoms": 4,
        "localized_antiparallel": True,
        "molecule_n_positive": 1,
        "molecule_n_negative": 1,
        "molecule_net": 0.0,
        "groups": {
            "molecule": {"atoms": 2, "spin_moment": 0.0, "fraction": 0.0,
                         "n_positive": 1, "n_negative": 1, "net": 0.0, "max_abs": 0.25},
            "electrode": {"atoms": 2, "spin_moment": 1.0, "fraction": 1.0},
        },
    }), encoding="utf-8")
    config = SimpleNamespace(molecule=SimpleNamespace(
        multiplicity=1, charge=0, fixed_spin_moment=False, broken_symmetry=True))
    result = ReviewResult(run_id="t")
    _spin_review(result, run_dir, config)
    assert "broken_symmetry_local_moments" in result.codes()
    assert "spin_localized_off_molecule" not in result.codes()
    issue = next(i for i in result.issues if i["code"] == "broken_symmetry_local_moments")
    assert issue["severity"] == "info"
    assert "局部反平行" in issue["message"]

    # 非 BS 的普通开壳层不产生该 code。
    config = SimpleNamespace(molecule=SimpleNamespace(
        multiplicity=2, charge=0, fixed_spin_moment=True, broken_symmetry=False))
    result = ReviewResult(run_id="t")
    _spin_review(result, run_dir, config)
    assert "broken_symmetry_local_moments" not in result.codes()
    assert "spin_localized_off_molecule" not in result.codes()


def test_bs_legacy_sidecar_uses_molecular_local_evidence(workdir: Path):
    from smj.report import _render_spin_section
    from smj.review import ReviewResult, _spin_review
    from smj.spin_report import build_spin_report, render_spin_report

    molecule = dict(multiplicity=1, broken_symmetry=True, charge=0)
    (workdir / "run_spec.json").write_text(json.dumps({"config": {"molecule": molecule}}))
    (workdir / "aims_scf").mkdir()
    # v1 sidecars saved by collect have raw atom moments and group counts.
    # A tiny residual junction moment makes the signed fraction misleading.
    partition = {
        "total": 0.001, "missing_moments": 0, "atoms": 3,
        "atom_moments": [
            {"index": i, "spin_moment": m}
            for i, m in enumerate([0.6, -0.6, 0.001], start=1)
        ],
        "groups": {
            "molecule": {"atoms": 2, "spin_moment": 0.0, "fraction": 0.0},
            "electrode": {"atoms": 1, "spin_moment": 0.001, "fraction": 1.0},
        },
    }
    (workdir / "aims_scf" / "spin_partition.json").write_text(json.dumps(partition))
    (workdir / "review.json").write_text(json.dumps({"issues": [
        {"code": "spin_localized_off_molecule", "message": "开壳层自旋主要在电极侧"}
    ]}))
    result = ReviewResult(run_id="bs")
    _spin_review(result, workdir, SimpleNamespace(molecule=SimpleNamespace(**molecule)))
    assert "spin_localized_off_molecule" not in result.codes()
    evidence = next(i for i in result.issues if i["code"] == "broken_symmetry_local_moments")
    lines = []
    _render_spin_section(lines, {"channels": 2, **molecule}, partition)
    cli_text = render_spin_report(build_spin_report(workdir))
    assert "spin_localized_off_molecule" not in cli_text
    for text in (evidence["message"], "\n".join(lines), cli_text):
        assert "局部反平行" in text
        assert "#1 +0.6000" in text and "#2 -0.6000" in text
        assert "1.2000" in text
        assert "100.00%" not in text

    molecule.update(multiplicity=2, broken_symmetry=False)
    (workdir / "run_spec.json").write_text(json.dumps({"config": {"molecule": molecule}}))
    cli_text = render_spin_report(build_spin_report(workdir))
    assert "spin_localized_off_molecule" not in cli_text
    assert "局部反平行判据" not in cli_text


def test_missing_local_moments_do_not_support_localization(workdir: Path):
    table = parse_hirshfeld(_write(workdir, _section([
        (1, "C", 0.25), (2, "C", -0.25), (3, "C", None), (4, "Au", 0.5),
    ])))
    partition = summarize(table, default_groups(3, 4))
    assert partition["localized_antiparallel"] is None
    assert partition["groups"]["molecule"]["abs_fraction"] is None


def test_bs_shows_mapped_seed_centers_when_other_sites_have_larger_moments(workdir: Path):
    from smj.report import _render_spin_section, _spin_partition_block
    from smj.review import ReviewResult, _spin_review
    from smj.spin_report import build_spin_report, render_spin_report

    molecule = {"multiplicity": 1, "broken_symmetry": True, "charge": 0,
                "initial_moments": {"3": 1.0, "5": -1.0}}
    (workdir / "run_spec.json").write_text(json.dumps({"config": {"molecule": molecule}}))
    stage = workdir / "aims_scf"
    stage.mkdir()
    (stage / "input_summary.json").write_text(json.dumps({"initial_moments": {"2": 1, "4": -1}}))
    table = parse_hirshfeld_text("\n".join(_section([
        (1, "C", 0.6), (2, "N", 0.12), (3, "C", -0.6), (4, "N", -0.12), (5, "Au", 0.001),
    ])))
    partition = {**table.as_dict(), "groups": {"molecule": {"atoms": 4}, "electrode": {"atoms": 1}}}
    (stage / "spin_partition.json").write_text(json.dumps(partition))
    result = ReviewResult(run_id="bs")
    _spin_review(result, workdir, SimpleNamespace(molecule=SimpleNamespace(**molecule)))
    message = next(i["message"] for i in result.issues if i["code"] == "broken_symmetry_local_moments")
    lines = []
    _render_spin_section(lines, molecule, _spin_partition_block(None, workdir))
    for text in (message, "\n".join(lines), render_spin_report(build_spin_report(workdir))):
        assert "N#2 +0.1200" in text and "N#4 -0.1200" in text
        assert "C#1 +0.6000" in text and "C#3 -0.6000" in text


def test_render_shows_spin_partition_rows():
    from smj.report import _render_spin_section

    lines: list[str] = []
    _render_spin_section(
        lines, {"channels": 2, "spin_polarization_at_ef": 0.0},
        {"total": 1.0, "groups": {
            "molecule": {"atoms": 14, "spin_moment": 0.0056, "fraction": 0.0056},
            "contact": {"atoms": 2, "spin_moment": 0.0063, "fraction": 0.0063},
            "electrode": {"atoms": 72, "spin_moment": 0.9881, "fraction": 0.9881},
        }},
    )
    text = "\n".join(lines)
    assert "自旋分布" in text
    assert "0.0056" in text and "0.9881" in text
    assert "接触区" in text
