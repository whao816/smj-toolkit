"""2026-11-21 第三轮复审（master 6384ce0d）8 条断言的逐条入库。

源文件：``D:\\codex_sessions\\2026-11-19\\https-github-com-whao816-smj-toolkit\\outputs\\test_smj_review_6384ce0.py``
（141 行，8 个测试）。这里逐条移植为仓内回归，全部用 ``workdir`` 夹具
（``tmp_path`` 在沙箱下会 PermissionError），断言口径保留原附件：
粗/细网格常数谱、相对路径变体、基组轴、真实阶段结果、多轴逐切片、空指标、
闭壳层 PDOS 占据、同路径改坐标。

附件测试 → 问题 → 修复项 → 批次（详见 ``docs/fix-plan-2026-11-21-r3.md`` §12）：

===============================  ==========================  ======  ====
附件测试                          问题                         修复项   批次
===============================  ==========================  ======  ====
coordinate_change_at_same_path   结构变化被降级成 schema 演进   R5      3
relative_molecule_path_survives  相对 XYZ 失效                R2      2
documented_basis_axis_changes    基组轴被静默忽略             R2      2
compare_reads_normal_analysis    汇总恒 missing               R3      2
multiaxis_convergence_slices     多轴漏掉变化组合             R4      2
empty_metric_records             空指标判 converged            R4      2
closed_shell_pdos_occupations    闭壳层取旧占据表             R6      3
underresolved_constant_spectrum  粗网格虚假热电势             R1      1
===============================  ==========================  ======  ====

复跑附件原版（两个开关：``tmp_path`` 的 0700 目录在沙箱下不可写；pytest 的
cacheprovider 会在 sessionfinish 用 ``tempfile.mkdtemp`` 自旋）::

    PYTHONPATH=<tmpfix 目录> python -m pytest -q -p tmpfix -p no:cacheprovider \\
        D:\\codex_sessions\\2026-11-19\\https-github-com-whao816-smj-toolkit\\
        outputs\\test_smj_review_6384ce0.py        # 期望 8 passed
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

from smj.config import load_config
from smj.junction.criteria import evaluate_run_criteria, physics_verdict
from smj.junction.landauer import analyze_te
from smj.pdos import analyze_pdos
from smj.runner import StageRunner

from test_pdos import EFERMI_HA  # noqa: F401  （与附件同口径，便于对照）
from test_pdos import _aims_out, _run as _pdos_run
from test_review_provenance import _codes, _prepared_run

ROOT = Path(__file__).resolve().parents[1]
EF_HARTREE = -0.18
H_EV = 27.2114


# --------------------------------------------------------------------------- #
# 共用夹具：附件里的 write_constant_te / planned_variants / _base_config
# --------------------------------------------------------------------------- #
def write_constant_te(directory, energies):
    """常数透射谱 T=0.25：首列 Hartree（efermi=-0.18），第二列 E−E_F (eV)。"""
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "TE.dat"
    lines = ["#bias = 0.0000 V efermi = -0.18 H"]
    lines.extend("%.12f %.12f 0.25" % (EF_HARTREE + e / H_EV, e) for e in energies)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _scan():
    spec = importlib.util.spec_from_file_location(
        "scan_convergence", ROOT / "scripts" / "scan_convergence.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _base_config(workdir: Path) -> Path:
    target = workdir / "smj.yaml"
    target.write_text(
        "run:\n  name: base\n  mode: local\n"
        "molecule:\n  xyz: mol.xyz\n"
        "electrode:\n  nlayers: 3\n  au_per_electrode: 20\n",
        encoding="utf-8")
    return target


def _planned_variants(base: Path, axes) -> list[dict]:
    scan = _scan()
    plan = scan.build_plan(base, list(axes), base.parent / "scan")
    return scan.load_manifest(scan.write_plan(plan))["variants"]


# --------------------------------------------------------------------------- #
# 附件第 1 条（R5）：同一路径下改坐标 = 物理漂移
# --------------------------------------------------------------------------- #
def test_coordinate_change_at_same_path_is_protocol_drift(workdir):
    config, run_dir = _prepared_run(workdir)
    xyz = config.molecule_xyz()
    before = xyz.read_text(encoding="utf-8")
    after = before.replace("H 0 0 1\n", "H 0 0 1.5\n")
    assert before != after
    xyz.write_text(after, encoding="utf-8")

    codes = _codes(config, run_dir)
    assert "provenance_drift" in codes, codes
    assert "provenance_schema_evolved" not in codes, codes


# --------------------------------------------------------------------------- #
# 附件第 2 条（R2）：相对 molecule.xyz 在 scan 变体里仍可解析
# --------------------------------------------------------------------------- #
def test_relative_molecule_path_survives_scan_plan(workdir):
    base = _base_config(workdir)
    source_xyz = workdir / "mol.xyz"
    source_xyz.write_text("1\nexample\nH 0 0 0\n", encoding="utf-8")

    variant = _planned_variants(base, [("electrode.nlayers", [3])])[0]
    target = load_config(variant["config"]).molecule_xyz()
    assert target is not None and target.is_file(), target
    assert target.read_bytes() == source_xyz.read_bytes(), target


# --------------------------------------------------------------------------- #
# 附件第 3 条（R2）：文档里的基组轴必须真正改变生效基组
# --------------------------------------------------------------------------- #
def test_documented_basis_axis_changes_the_effective_basis(workdir):
    base = _base_config(workdir)
    tiers = [(workdir / tier).resolve() for tier in ("light", "tight")]
    for tier in tiers:
        tier.mkdir()
    with base.open("a", encoding="utf-8") as stream:
        stream.write("paths:\n  species_defaults: " + tiers[0].as_posix() + "\n")

    axis = ("paths.species_dir", [path.as_posix() for path in tiers])
    variants = _planned_variants(base, [axis])
    actual = [StageRunner(load_config(v["config"])).context.species_dir() for v in variants]
    assert actual == tiers, actual


# --------------------------------------------------------------------------- #
# 附件第 4 条（R3）：汇总读真实阶段写出的 results.json
# --------------------------------------------------------------------------- #
def test_compare_reads_the_normal_local_analysis_stage(workdir):
    scan = _scan()
    variant = _planned_variants(_base_config(workdir), [("electrode.nlayers", [3])])[0]
    config = load_config(variant["config"])
    write_constant_te(config.run_dir / "aitranss_te", (i / 100 for i in range(-50, 51)))

    runner = StageRunner(config, logger=lambda message: None)
    runner.prepare("junction_analyze")
    assert (config.run_dir / "junction_analyze" / "results.json").is_file()

    metrics = scan.variant_metrics(variant)
    assert Path(variant["run_dir"]) == config.run_dir
    assert metrics["status"] == "ok", metrics
    assert metrics["run_dir"] == str(config.run_dir), metrics
    assert metrics["values"]["S_uV_per_K"] is not None, metrics


# --------------------------------------------------------------------------- #
# 附件第 5 条（R4）：多轴时每个固定切片都要比较
# --------------------------------------------------------------------------- #
def test_multiaxis_convergence_checks_each_fixed_slice():
    scan = _scan()
    axis = {"key": "electrode.nlayers", "values": [3, 4]}
    other = "electrode.au_per_electrode"
    rows = []
    for nlayers in axis["values"]:
        for au in (20, 30):
            values = dict.fromkeys(scan.METRICS, 1.0)
            values["S_uV_per_K"] = 20.0 if (nlayers, au) == (4, 30) else 10.0
            rows.append({"name": "n%d-au%d" % (nlayers, au),
                         "values": {axis["key"]: nlayers, other: au},
                         "metrics": {"status": "ok", "values": values}})

    verdicts = [scan.compare_axis(a, rows, 0.05)["verdict"]
                for a in (axis, {"key": other, "values": [20, 30]})]
    assert verdicts == ["not-converged", "not-converged"], verdicts


# --------------------------------------------------------------------------- #
# 附件第 6 条（R4）：空指标不能证明收敛
# --------------------------------------------------------------------------- #
def test_empty_metric_records_do_not_certify_convergence():
    scan = _scan()
    axis = {"key": "electrode.nlayers", "values": [3, 4]}
    rows = [{"values": {axis["key"]: nlayers},
             "metrics": {"status": "ok", "values": {}}} for nlayers in (3, 4)]

    result = scan.compare_axis(axis, rows, 0.05)
    assert result["verdict"] != "converged", result
    assert result["verdict"] == "incomplete", result


# --------------------------------------------------------------------------- #
# 附件第 7 条（R6）：闭壳层 PDOS 不受上一轮迭代的占据表影响
# --------------------------------------------------------------------------- #
def test_closed_shell_pdos_is_independent_of_previous_occupations(workdir):
    run = _pdos_run(workdir, natoms=2, offsets=(-0.2, 0.2), grid=(-0.2, 0.0, 0.2))
    aims_out = run / "aims_scf" / "aims.out"
    aims_out.write_text(_aims_out([2.0, 0.0]), encoding="utf-8")
    expected = analyze_pdos(run, groups="mol=1-2", sigma_ev=0.01)["channels"]["mos"]["pdos"]["mol"]

    # 上一轮迭代：占据 [1, 1]、能量都在 0 附近；终态表在最后。
    aims_out.write_text(
        _aims_out([1.0, 1.0]) + "End of iteration\n" + _aims_out([2.0, 0.0]),
        encoding="utf-8")
    spectrum = analyze_pdos(run, groups="mol=1-2", sigma_ev=0.01)["channels"]["mos"]["pdos"]["mol"]
    assert spectrum == pytest.approx(expected)


# --------------------------------------------------------------------------- #
# 附件第 8 条（R1）：欠采样常数谱不能带着虚假 Seebeck 判 passed
# --------------------------------------------------------------------------- #
def test_underresolved_constant_spectrum_cannot_pass_with_spurious_seebeck(workdir):
    coarse = write_constant_te(workdir / "coarse", (-1.065 + i * 0.2 for i in range(12)))
    fine = write_constant_te(workdir / "fine", (-1.065 + i * 0.002 for i in range(1101)))

    reference = analyze_te(fine, temperature=300)
    assert reference["L0"] == pytest.approx(0.25, abs=1e-6)
    assert abs(reference["S_uV_per_K"]) < 1e-6

    result = analyze_te(coarse, temperature=300)
    checks = evaluate_run_criteria(coarse.parent, temperature=300, system="BDT")
    verdict = physics_verdict(checks)["physics_verdict"]

    correct = abs(result["L0"] - 0.25) < 1e-6 and abs(result["S_uV_per_K"]) < 1e-6
    assert correct or verdict != "passed", (
        verdict, result["thermal_window"]["coverage"], result["L0"], result["S_uV_per_K"])

    # 判据层面：求积不足必须是 error 级、missing（证据不足），而不是静默放过。
    by_name = {c["name"]: c for c in checks}
    quad = by_name["te_fermi_quadrature"]
    assert quad["severity"] == "error"
    assert quad["passed"] is False
    assert quad["missing"] is True
    assert "求积" in quad["detail"] and "kBT" in quad["detail"]

    summary = physics_verdict(checks)
    assert summary["physics_verdict"] == "inconclusive"
    assert "te_fermi_quadrature" in summary["missing_evidence_criteria"]
