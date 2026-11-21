"""2026-11-21 第三轮复审（master 6384ce0d）回归 —— 批 1（R1）+ 批 2（R2–R4）。

批 1 对应附件 ``test_smj_review_6384ce0.py`` 第 8 条：

    粗网格常数谱（12 点，步长 0.2 eV）在 300 K 下只有 1 个采样点落在 ±kBT 里，
    梯形求积把权重积分低估 42%、一阶矩变成 -1.96 kBT —— 于是常数谱也报出
    S=+169 uV/K。旧判据只看端点覆盖（f(E_min)-f(E_max)=1.0），照样判 passed。

R1 用"求积充分性"（权重积分/精确覆盖、权重一阶矩）来区分"网格够密"与
"网格太粗"，缺失时按证据不足（inconclusive）归约，不判物理否证。

批 2 对应附件第 2–6 条（扫描工具层）：

- R2 变体 yaml 把相对路径写进 `scan/vNNN/` 后失效；`paths.species_dir` 不是
  合法配置键，基组轴被静默忽略（两个变体都落到基线档位）。
- R3 汇总读 `<variant>/run/junction_analyze/analysis.json`（全仓无生产者），
  真实阶段写的是 `<variant>/run/<name>/junction_analyze/results.json`。
- R4 多轴只比较每个轴值的**第一行**，漏掉其余切片；空指标时 `worst` 保持 0.0
  直接判 converged。
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

from smj.config import ConfigError, load_config
from smj.junction.criteria import evaluate_run_criteria, physics_verdict
from smj.junction.landauer import analyze_te
from smj.runner import StageRunner

ROOT = Path(__file__).resolve().parents[1]
EF_HARTREE = -0.18
H_EV = 27.2114


def write_constant_te(directory, energies):
    """常数透射谱 T=0.25：首列 Hartree（efermi=-0.18），第二列 E-E_F (eV)。"""
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "TE.dat"
    lines = ["#bias = 0.0000 V efermi = -0.18 H"]
    lines.extend("%.12f %.12f 0.25" % (EF_HARTREE + e / H_EV, e) for e in energies)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


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


# --------------------------------------------------------------------------- #
# 批 2（R2–R4）：扫描工具层
# --------------------------------------------------------------------------- #
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


def test_relative_molecule_path_survives_scan_plan(workdir):
    base = _base_config(workdir)
    source_xyz = workdir / "mol.xyz"
    source_xyz.write_text("1\nexample\nH 0 0 0\n", encoding="utf-8")

    variant = _planned_variants(base, [("electrode.nlayers", [3])])[0]
    target = load_config(variant["config"]).molecule_xyz()
    assert target is not None and target.is_file(), target
    assert target.read_bytes() == source_xyz.read_bytes(), target


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


def test_empty_metric_records_do_not_certify_convergence():
    scan = _scan()
    axis = {"key": "electrode.nlayers", "values": [3, 4]}
    rows = [{"values": {axis["key"]: nlayers},
             "metrics": {"status": "ok", "values": {}}} for nlayers in (3, 4)]

    result = scan.compare_axis(axis, rows, 0.05)
    assert result["verdict"] != "converged", result
    assert result["verdict"] == "incomplete", result


def test_plan_rejects_unknown_axis_key(workdir):
    scan = _scan()
    with pytest.raises(scan.ScanError, match="electrode.typox"):
        scan.build_plan(_base_config(workdir), [("electrode.typox", [1])], workdir / "scan")


def test_paths_species_dir_alias_and_conflict(workdir):
    light = (workdir / "light").resolve()
    tight = (workdir / "tight").resolve()
    alias = workdir / "alias.yaml"
    alias.write_text(
        "run:\n  name: r\n  mode: local\nmolecule:\n  xyz: mol.xyz\n"
        "paths:\n  species_dir: " + light.as_posix() + "\n", encoding="utf-8")
    config = load_config(alias)
    assert config.paths.species_defaults is not None
    assert Path(config.paths.species_defaults) == light
    assert any("species_dir 已并入" in warning for warning in config.warnings), config.warnings
    assert not any("未知字段" in warning for warning in config.warnings), config.warnings

    clash = workdir / "clash.yaml"
    clash.write_text(
        "run:\n  name: r\n  mode: local\nmolecule:\n  xyz: mol.xyz\n"
        "paths:\n  species_dir: " + light.as_posix() + "\n"
        "  species_defaults: " + tight.as_posix() + "\n", encoding="utf-8")
    with pytest.raises(ConfigError, match="species_dir"):
        load_config(clash)
