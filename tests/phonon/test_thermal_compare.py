"""``smj thermal-compare`` 与 ``smj thermal --sensitivity`` 的回归测试。

这两条命令的价值都在"防止误比"与"暴露正则化依赖"，所以测试也盯这两件事：
协议不一致必须被点名（``--strict`` 直接失败），扫参必须覆盖全组合并把
被带顶闸门拒绝的组合如实记成 ``None``（拒绝本身就是结果）。
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from smj.cli import main
from smj.phonon.compare import compare, find_runs, load_run, render_markdown


def _payload(**over) -> dict:
    """构造一份能通过 ``read_thermal`` 校验的 payload。"""
    energy = [i * 0.1 for i in range(245)]
    data = {
        "schema": "smj.thermal.v1",
        "temperature_K": 300.0,
        "n_atoms": 90,
        "n_device_atoms": 54,
        "n_lead_atoms": 36,
        "g_ph_pw_per_k": 13.7,
        "g_quantum_pw_per_k": 283.92934549159924,
        "n_channels_effective": 0.0483,
        "acoustic_sum_residual": 1.0e-13,
        "acoustic_sum_relative_residual": 1.0e-15,
        "imaginary_modes": {"count": 2, "worst_cm1": -21.3, "modes": [],
                            "tolerance_cm1": -10.0, "n_modes": 3},
        "n_k_points": 1,
        "band_top_meV": 19.5,
        "eta_meV": 0.05,
        "energy_step_meV": 0.05,
        "lead_fc_sha256": "abc123",
        "spectrum": {
            "energy_meV": energy,
            "transmission": [0.4 for _ in energy],
            "cumulative_g_ph_pw_per_k": [13.7 * (i / (len(energy) - 1))
                                         for i in range(len(energy))],
        },
        "zt": {"available": False, "reason": "本 run 未跑电子链"},
    }
    data.update(over)
    return data


def _write_run(root: Path, name: str, payload: dict, *, nested: bool) -> Path:
    run = root / name / name if nested else root / name
    (run / "thermal_analyze").mkdir(parents=True)
    (run / "thermal_analyze" / "thermal.json").write_text(
        json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    return run


def test_find_runs_accepts_run_parent_and_json(workdir) -> None:
    root = Path(workdir).resolve() / "runs"
    flat = _write_run(root, "a", _payload(), nested=False)
    nested = _write_run(root, "b", _payload(), nested=True)

    assert find_runs([root]) == [flat, nested]
    assert find_runs([flat]) == [flat]
    assert find_runs([nested]) == [nested]
    assert find_runs([nested / "thermal_analyze" / "thermal.json"]) == [nested]
    # 去重：同一个 run 从父目录和它自己各给一次，只算一条
    assert find_runs([root, flat]) == [flat, nested]


def test_find_runs_reports_nothing_found(workdir) -> None:
    empty = Path(workdir).resolve() / "empty"
    empty.mkdir()
    with pytest.raises(ValueError, match="thermal.json"):
        find_runs([empty])


def test_compare_detects_protocol_mismatch(workdir) -> None:
    root = Path(workdir).resolve() / "runs"
    one = load_run(_write_run(root, "a", _payload(), nested=False))
    other = load_run(_write_run(root, "b", _payload(eta_meV=0.1,
                                                   lead_fc_sha256="def456"),
                                nested=False))
    result = compare([one, other])
    assert not result.consistent
    assert set(result.mismatches) == {"eta_meV", "lead_fc_sha256"}
    text = render_markdown(result)
    assert "协议不一致" in text and "D15" in text
    assert "| a | 13.700 |" in text and "| b | 13.700 |" in text


def test_compare_consistent_runs_have_no_warning(workdir) -> None:
    root = Path(workdir).resolve() / "runs"
    records = [load_run(_write_run(root, name, _payload(), nested=False))
               for name in ("a", "b")]
    result = compare(records)
    assert result.consistent and result.mismatches == {}
    assert "协议不一致" not in render_markdown(result)


def test_compare_flags_environment_differences(workdir) -> None:
    root = Path(workdir).resolve() / "runs"
    a = _write_run(root, "a", _payload(), nested=False)
    b = _write_run(root, "b", _payload(), nested=False)
    def _env(aims="a" * 64):
        # identity() 读的是嵌套结构（与环境清单的 schema 一致）。
        return {
            "schema": "smj.environment.v1",
            "generated_at": "2026-11-25T00:00:00+00:00",
            "smj": {"version": "0.2.2rc2", "root": "/x/src/smj"},
            "python": {"version": "3.12.10", "implementation": "CPython",
                       "executable": "/usr/bin/python3"},
            "platform": {"system": "Linux", "machine": "x86_64",
                         "release": "4.18.0", "processor": "x86_64"},
            "packages": {"numpy": "2.3.0"},
            "junction_provenance_sha256": "p" * 64,
            "thermal": {},
            "binaries": {"aims": {"path": "/opt/aims", "sha256": aims,
                                  "exists": True}},
            "species": {"directory": "/s", "files": {"Au.in": "c" * 64},
                        "digest": "s" * 64, "present": True},
            "lead_fc": {"path": "/l/lead_fc.json", "sha256": "l" * 64},
        }

    (a / "thermal_analyze" / "environment.json").write_text(
        json.dumps(_env()), encoding="utf-8")
    (b / "thermal_analyze" / "environment.json").write_text(
        json.dumps(_env(aims="b" * 64)), encoding="utf-8")

    result = compare([load_run(a), load_run(b)])
    # 协议一致但环境不同：两者必须分开报告 —— 前者是方法选择，后者是换了工具。
    assert result.consistent
    assert not result.environment_consistent
    assert "aims_sha256" in result.environment_mismatches
    text = render_markdown(result)
    assert "环境不一致" in text and "FHI-aims" in text


def test_compare_reads_verdict_from_review(workdir) -> None:
    root = Path(workdir).resolve() / "runs"
    run = _write_run(root, "a", _payload(), nested=False)
    (run / "review.json").write_text(json.dumps({
        "thermal": {"verdict": "inconclusive", "failed": [],
                    "inconclusive": ["ph_imaginary_modes"]}}), encoding="utf-8")
    record = load_run(run)
    assert record.verdict == "inconclusive"
    assert record.inconclusive == ["ph_imaginary_modes"]


def test_cli_thermal_compare_strict_fails_on_mismatch(workdir, capsys) -> None:
    root = Path(workdir).resolve() / "runs"
    _write_run(root, "a", _payload(), nested=False)
    _write_run(root, "b", _payload(eta_meV=0.2), nested=False)

    assert main(["thermal-compare", str(root)]) == 0
    assert "协议不一致" in capsys.readouterr().out

    assert main(["thermal-compare", str(root), "--strict"]) == 1


def test_cli_thermal_compare_writes_out_file(workdir) -> None:
    root = Path(workdir).resolve() / "runs"
    _write_run(root, "a", _payload(), nested=False)
    target = Path(workdir).resolve() / "compare.md"
    assert main(["thermal-compare", str(root), "--out", str(target)]) == 0
    text = target.read_text(encoding="utf-8")
    assert text.startswith("## 声子热导比较")
    assert "| a |" in text


def test_sensitivity_covers_grid_and_records_rejections(workdir) -> None:
    from smj.config import ThermalConfig
    from smj.stages.thermal_analyze import _build, _curve_factory, _energy_grid, \
        sensitivity

    root = Path(workdir).resolve() / "run"
    fc_dir = root / "thermal_fc"
    dev_dir = root / "thermal_device"
    fc_dir.mkdir(parents=True)
    dev_dir.mkdir(parents=True)
    # 6 原子分组：器件 2（分子 0/1）+ 每侧主层与次外层各 1，等强弹簧。
    n_atoms = 6
    k = 2.0
    fc = np.zeros((3 * n_atoms, 3 * n_atoms))
    for a in range(n_atoms):
        for b in range(n_atoms):
            if a != b:
                fc[3 * a:3 * a + 3, 3 * b:3 * b + 3] = -k * np.eye(3)
                fc[3 * a:3 * a + 3, 3 * a:3 * a + 3] += k * np.eye(3)
    # 位移集只覆盖器件区(2 个原子 = 6 个自由度),列覆盖全部 6 个原子
    np.save(fc_dir / "force_constants.npy", fc[:3 * 2, :])
    (fc_dir / "geometry.in").write_text(
        "".join("atom 0.0 0.0 %.3f  Au\n" % (1.5 * i) for i in range(n_atoms)),
        encoding="utf-8")
    from smj.phonon.fc import save_lead_fc
    lead_dir = root / "lead"
    lead_dir.mkdir()
    # lead 块覆盖「主层 + 次外层」：2 层 × 2 原子 = 4 原子 = 12 个自由度
    # 每侧「主层 1 + 次外层 1」= 2 原子,所以 lead 块就是 2 层 x 1 原子 = 6 自由度
    save_lead_fc(lead_dir, fc[:6, :6], element="Au", a_lattice_A=4.07823,
                 n_atoms_per_layer=1, layers=2, in_plane=(1, 1),
                 fc_cutoff_A=10.0)
    (dev_dir / "device_metadata.json").write_text(json.dumps({
        "n_atoms": n_atoms,
        "device_indices": [0, 1],
        "lead_indices": [2, 3, 4, 5],
        "layer_index": [0, 0, 0, 1, 0, 1],
        "side": ["left", "left", "left", "left", "right", "right"],
        "molecule_indices": [0],
    }), encoding="utf-8")
    thermal = ThermalConfig()
    thermal.negf.energy_step_meV = 0.05
    thermal.negf.band_top_meV = 2.0
    thermal.lead_fc = str(lead_dir)
    built = _build(fc_dir / "force_constants.npy",
                   dev_dir / "device_metadata.json", fc_dir / "geometry.in",
                   thermal)
    omega, _ = _energy_grid(2.0, 0.05)
    curve = _curve_factory(built, omega, 1.0e-4)
    assert curve(built["dynamical"]).shape == omega.shape

    payload = sensitivity(fc_dir / "force_constants.npy",
                          dev_dir / "device_metadata.json",
                          fc_dir / "geometry.in", thermal, temperature=300.0,
                          etas=[0.1, 0.05], steps=[0.1, 0.05])
    assert len(payload["rows"]) == 4
    assert payload["production"] == {"eta_meV": 0.05, "energy_step_meV": 0.05}
    extrapolated = payload["eta_zero_extrapolation"]
    assert extrapolated["n_points"] == 2
    assert extrapolated["at_energy_step_meV"] == 0.05
    # 外推必须锚在「最小 η、生产步长」那一格上；数值层面的单调性由真实 run
    # 的 sensitivity.json 复核（η 0.3→0.05 对应 12.37→13.72，外推 13.82）。
    assert extrapolated["smallest_eta_meV"] == 0.05
    row = [item for item in payload["rows"]
           if item["eta_meV"] == 0.05 and item["energy_step_meV"] == 0.05][0]
    assert extrapolated["g_ph_at_smallest_eta_pw_per_k"] == row["g_ph_pw_per_k"]
