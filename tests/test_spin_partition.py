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
