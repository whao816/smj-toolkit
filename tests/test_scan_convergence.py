"""`scripts/scan_convergence.py` 的行为：矩阵展开、缺结果、收敛判定。

脚本不在包里（`scripts/`），所以用 importlib 从路径加载 —— 与
`tests/test_demo_regression.py` 的 `_demo_script()` 同一手法。
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "scan_convergence.py"


def _scan_script():
    spec = importlib.util.spec_from_file_location("scan_convergence", SCRIPT)
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


def _analysis(workdir: Path, variant: dict, s: float, g: float = -1.0) -> None:
    target = Path(variant["dir"]) / "run" / "junction_analyze"
    target.mkdir(parents=True, exist_ok=True)
    payload = {"S_uV_per_K": s, "log10_G_over_G0_like": g,
               "ZT_electronic": 0.03, "k_e_pW_per_K": 388.0}
    (target / "analysis.json").write_text(json.dumps(payload), encoding="utf-8")


def test_typed_and_axis_parsing():
    scan = _scan_script()
    assert scan.typed("3") == 3
    assert scan.typed("2.5") == 2.5
    assert scan.typed("true") is True
    assert scan.typed("None") is None
    assert scan.typed("/data/species/tight") == "/data/species/tight"
    assert scan.parse_axis("electrode.nlayers=3,4,5") == ("electrode.nlayers", [3, 4, 5])
    assert scan.parse_axis("a.b=x, x, y") == ("a.b", ["x", "y"])
    with pytest.raises(scan.ScanError):
        scan.parse_axis("electrode.nlayers")
    with pytest.raises(scan.ScanError):
        scan.parse_axis("nlayers=3")


def test_plan_expands_axes_and_writes_configs(workdir):
    scan = _scan_script()
    base = _base_config(workdir)
    plan = scan.build_plan(base, [("electrode.nlayers", [3, 4]),
                                  ("electrode.au_per_electrode", [20, 30])],
                           workdir / "scan")
    assert len(plan["variants"]) == 4
    manifest = scan.write_plan(plan)
    data = json.loads(manifest.read_text(encoding="utf-8"))
    assert data["schema"] == scan.SCHEMA
    assert [v["name"] for v in data["variants"]][0].startswith("v001__")
    first = json.loads(json.dumps(plan["variants"][0]["config"]))
    assert first["run"]["dir"] == "run"
    assert first["electrode"]["au_per_electrode"] == 20
    assert (workdir / "scan" / data["variants"][0]["name"] / "smj.yaml").is_file()


def test_plan_applies_constants_and_rejects_duplicate_axis(workdir):
    scan = _scan_script()
    base = _base_config(workdir)
    plan = scan.build_plan(base, [("electrode.nlayers", [3, 4])], workdir / "scan",
                           constants={"cluster.partition": "normal"})
    assert all(v["config"]["cluster"]["partition"] == "normal" for v in plan["variants"])
    with pytest.raises(scan.ScanError):
        scan.build_plan(base, [("electrode.nlayers", [3]), ("electrode.nlayers", [4])],
                        workdir / "scan")


def test_compare_marks_converged_only_when_successive_steps_are_small(workdir):
    scan = _scan_script()
    base = _base_config(workdir)
    plan = scan.build_plan(base, [("electrode.nlayers", [3, 4, 5])], workdir / "scan")
    manifest = scan.write_plan(plan)
    data = scan.load_manifest(manifest)
    for variant, s in zip(data["variants"], (10.0, 10.2, 10.3)):
        _analysis(workdir, variant, s)
    rows = [{"name": v["name"], "values": v["values"],
             "metrics": scan.variant_metrics(v)} for v in data["variants"]]
    result = scan.compare_axis(data["axes"][0], rows, 0.05)
    assert result["verdict"] == "converged"
    assert all(change["converged"] for change in result["changes"])


def test_compare_flags_a_real_step(workdir):
    scan = _scan_script()
    base = _base_config(workdir)
    plan = scan.build_plan(base, [("electrode.nlayers", [3, 4])], workdir / "scan")
    manifest = scan.write_plan(plan)
    data = scan.load_manifest(manifest)
    for variant, s in zip(data["variants"], (10.0, 13.0)):
        _analysis(workdir, variant, s)
    rows = [{"name": v["name"], "values": v["values"],
             "metrics": scan.variant_metrics(v)} for v in data["variants"]]
    result = scan.compare_axis(data["axes"][0], rows, 0.05)
    assert result["verdict"] == "not-converged"
    assert result["changes"][0]["worst_relative_change"] == pytest.approx(0.3)


def test_compare_reports_missing_results_instead_of_zero(workdir):
    scan = _scan_script()
    base = _base_config(workdir)
    plan = scan.build_plan(base, [("electrode.nlayers", [3, 4])], workdir / "scan")
    manifest = scan.write_plan(plan)
    data = scan.load_manifest(manifest)
    _analysis(workdir, data["variants"][0], 10.0)   # 第二个变体故意不写
    rows = [{"name": v["name"], "values": v["values"],
             "metrics": scan.variant_metrics(v)} for v in data["variants"]]
    assert rows[1]["metrics"]["status"] == "missing"
    assert rows[1]["metrics"]["values"] == {}
    assert scan.compare_axis(data["axes"][0], rows, 0.05)["verdict"] == "incomplete"


def test_run_without_yes_prints_commands_and_executes_nothing(workdir, capsys):
    scan = _scan_script()
    base = _base_config(workdir)
    manifest = scan.write_plan(scan.build_plan(
        base, [("electrode.nlayers", [3, 4])], workdir / "scan"))
    rc = scan.main(["run", str(manifest)])
    captured = capsys.readouterr().out
    assert rc == 0
    assert "尚未执行" in captured
    assert captured.count("smj.cli run") == 2
    # 没有 run/state.json 说明一个作业都没跑
    assert not list((workdir / "scan").glob("*/run/state.json"))
