"""``smj.phonon.environment`` 与 ``cutoff_sensitivity`` 的单元测试。

环境清单的意义是让"另一个团队能不能复现"变成可比对的字段，所以这里钉的是
**比对口径**（哪些键参与 identity、什么该被排除）而不只是字段存在性。
截断轴的测试则钉两件事：切到不超过矩阵实有跨度时结果不变、跨度要如实报出。
"""

from __future__ import annotations

import json

import numpy as np

from smj.phonon import environment as env


def _payload(**over):
    base = {
        "schema": env.SCHEMA,
        "generated_at": "2026-11-25T00:00:00+00:00",
        "smj": {"version": "0.2.2rc2", "root": "/x/src/smj"},
        "python": {"version": "3.12.10", "implementation": "CPython",
                   "executable": "/usr/bin/python3"},
        "platform": {"system": "Linux", "machine": "x86_64",
                     "release": "4.18.0", "processor": "x86_64"},
        "packages": {"numpy": "2.3.0", "scipy": "1.14.0",
                     "matplotlib": "3.9.0", "PyYAML": "6.0.1"},
        "junction_provenance_sha256": "a" * 64,
        "thermal": {"enabled": True, "temperature": None,
                    "electrode": {"in_plane": (3, 3), "layers": 4}},
        "binaries": {"aims": {"path": "/opt/aims", "sha256": "b" * 64,
                              "exists": True}},
        "species": {"directory": "/spec", "files": {"Au.in": "c" * 64},
                    "digest": "d" * 64, "present": True},
        "lead_fc": {"path": "/lead/lead_fc.json", "sha256": "e" * 64,
                    "in_plane": [3, 3], "layers": 2},
    }
    base.update(over)
    return base


def test_identity_keeps_only_the_six_keys_that_change_numbers():
    identity = env.identity(_payload())
    assert set(identity) == set(env.IDENTITY_KEYS)
    assert identity["aims_sha256"] == "b" * 64
    assert identity["species_digest"] == "d" * 64
    assert identity["lead_fc_sha256"] == "e" * 64


def test_identity_ignores_paths_so_two_installs_compare_equal():
    one = env.identity(_payload())
    other = env.identity(_payload(
        smj={"version": "0.2.2rc2", "root": "/somewhere/else/src/smj"},
        binaries={"aims": {"path": "/other/aims", "sha256": "b" * 64,
                           "exists": True}},
        species={"directory": "/other/spec", "files": {"Au.in": "c" * 64},
                 "digest": "d" * 64, "present": True}))
    # 装在不同目录但二进制与 species 内容一致 —— 这才是"同一套工具"。
    assert one == other


def test_digest_ignores_generated_at_but_catches_a_changed_binary():
    base = env.digest(_payload())
    assert base == env.digest(_payload(generated_at="2027-01-01T00:00:00+00:00"))
    changed = _payload(binaries={"aims": {"path": "/opt/aims",
                                          "sha256": "f" * 64, "exists": True}})
    assert env.digest(changed) != base


def test_canonical_is_stable_json():
    text = env.canonical(_payload())
    assert text == json.dumps(json.loads(text), sort_keys=True,
                              separators=(",", ":"))
    assert "generated_at" not in text


def test_describe_lines_reports_the_five_facts():
    lines = env.describe_lines(_payload())
    joined = "\n".join(lines)
    assert len(lines) == 5
    assert "CPython" in joined and "3.12.10" in joined
    assert "numpy" in joined
    assert "e" * 12 in joined            # lead FC 的 sha 前缀
    assert "b" * 12 in joined            # FHI-aims 的 sha 前缀


def test_describe_lines_flags_missing_aims_instead_of_pretending():
    lines = env.describe_lines(_payload(
        binaries={"aims": {"path": "/opt/aims", "sha256": None,
                           "exists": False}}))
    joined = "\n".join(lines)
    assert "未配置" in joined
    assert "—" in joined


# --- 截断轴 -------------------------------------------------------------

def _chain_built():
    """6 原子一维链：lead(0,1) — 中心(2,3) — lead(4,5)。

    最近邻弹簧给通带，另外在**两个 lead 原子之间**（0 与 5，相距 12.5 Å）
    加一条 1e-3 量级的长程耦合：截断轴的 ``effective_span_A`` 只统计 lead
    块内的非零元素，这一对正好把跨度顶到 12.5 Å。
    """
    from smj.phonon.units import ELEMENT_MASS_AMU, EV_PER_SQRT

    n = 6
    positions = np.zeros((n, 3))
    positions[:, 0] = np.arange(n) * 2.5
    mass = ELEMENT_MASS_AMU["Au"]

    k = 30.0
    K = np.zeros((3 * n, 3 * n))
    for i in range(n - 1):
        block = np.eye(3) * k
        K[3 * i:3 * i + 3, 3 * i + 3:3 * i + 6] -= block
        K[3 * i + 3:3 * i + 6, 3 * i:3 * i + 3] -= block
        K[3 * i:3 * i + 3, 3 * i:3 * i + 3] += block
        K[3 * i + 3:3 * i + 6, 3 * i + 3:3 * i + 6] += block
    matrix = K / mass * EV_PER_SQRT ** 2
    long_range = 1.0e-6
    matrix[0:3, 15:18] -= long_range * np.eye(3)
    matrix[15:18, 0:3] -= long_range * np.eye(3)
    return {
        "positions": positions,
        "dynamical": matrix,
        "lead_indices": [0, 1, 4, 5],
        "left_pl": [1], "left_outer": [0],
        "right_pl": [4], "right_outer": [5],
        "centre": [2, 3],
    }


def test_cutoff_above_the_real_span_changes_nothing():
    from smj.stages.thermal_analyze import _energy_grid, cutoff_sensitivity

    built = _chain_built()
    omega, _ = _energy_grid(40.0, 2.0)
    result = cutoff_sensitivity(built, omega, 0.5 * 1e-3, temperature=300.0,
                                band_top=40.0, cutoffs=(25.0, 40.0))
    values = [row["g_ph_pw_per_k"] for row in result["rows"]]
    # 跨度以内不切，所以两档必须逐位相同（含都被闸门拒的情形）。
    assert values[0] == values[1]
    assert result["relative_spread"] in (None, 0.0)


def test_effective_span_reports_the_farthest_nonzero_pair():
    from smj.stages.thermal_analyze import _energy_grid, cutoff_sensitivity

    built = _chain_built()
    omega, _ = _energy_grid(40.0, 2.0)
    result = cutoff_sensitivity(built, omega, 0.5 * 1e-3, temperature=300.0,
                                band_top=40.0)
    # 长程项只连 lead 原子 0 与 5，相距 12.5 Å。
    assert abs(result["effective_span_A"] - 12.5) < 1e-9


def test_cutting_below_the_span_drops_pairs():
    from smj.stages.thermal_analyze import _energy_grid, cutoff_sensitivity

    built = _chain_built()
    omega, _ = _energy_grid(40.0, 2.0)
    result = cutoff_sensitivity(built, omega, 0.5 * 1e-3, temperature=300.0,
                                band_top=40.0, cutoffs=(7.0, 20.0))
    tight, loose = result["rows"]
    # lead 共 4 个原子、12 个有序非对角对；7 Å 以内只有 (0,1) 与 (4,5) 两条
    # 各两个方向，所以 20 Å 档比 7 Å 档多 8 对 —— 掩码确实按距离在切。
    assert loose["pairs_kept"] == 12
    assert tight["pairs_kept"] == 4


def test_equivalent_cutoff_matrices_share_one_transmission(monkeypatch):
    from smj.stages import thermal_analyze as stage

    calls = []
    omega, _ = stage._energy_grid(40.0, 2.0)

    def curve(matrix):
        calls.append(matrix)
        return np.zeros_like(omega)

    monkeypatch.setattr(stage, "_curve_factory", lambda *args: curve)
    result = stage.cutoff_sensitivity(_chain_built(), omega, 0.0005, temperature=300,
                                      band_top=40, cutoffs=(7, 20, 25, 40))
    assert len(calls) == 2
    assert len(result["rows"]) == 4


def test_convergence_does_not_call_a_covering_cutoff_unconverged():
    """名义截断覆盖块内全部非零对时，不能判"未收敛"。

    ``effective_span_A`` 是组装矩阵里**实际存在**的最远非零 lead 对。它
    ≤ 名义 ``fc_cutoff_A`` 只说明"再按名义值切一刀不改变矩阵"，既不是
    收敛证据，也不是未收敛证据：更远的体相耦合本来就没进这份资产，
    后处理无法验证。以前把"没有 shortfall"当成 ``converged=False``，
    会在报告里给出一个没有理由的未收敛警告（2026-11-25 实测：span
    7.08 Å < cutoff 10 Å，G_ph 在 8/10/12 Å 三档逐位相同）。
    """
    from types import SimpleNamespace

    from smj.stages.thermal_analyze import convergence_axes

    built = {
        "spec": {"electrode_layers": 4, "device_layers": 2},
        "lead": SimpleNamespace(fc_cutoff_A=10.0),
    }
    block = {
        "values_A": [8.0, 10.0, 12.0],
        "rows": [
            {"cutoff_A": 8.0, "g_ph_pw_per_k": 15.92639, "rejected": None},
            {"cutoff_A": 10.0, "g_ph_pw_per_k": 15.92639, "rejected": None},
            {"cutoff_A": 12.0, "g_ph_pw_per_k": 15.92639, "rejected": None},
        ],
        "relative_spread": 0.0,
        "effective_span_A": 7.082,
    }
    axis = convergence_axes(built, block)["fc_cutoff"]
    assert axis["nominal_shortfall"] is None
    assert axis["nominal_covers_span"] is True
    assert axis["beyond_nominal_checked"] is False
    assert axis["converged"] is None


def test_cutoff_uses_the_uncut_result_for_a_nonmonotonic_scan():
    from types import SimpleNamespace
    from smj.stages.thermal_analyze import convergence_axes

    built = {"lead": SimpleNamespace(fc_cutoff_A=10.0)}
    block = {"effective_span_A": 11., "uncut_g_ph_pw_per_k": 10.,
             "rows": [{"cutoff_A": 8., "g_ph_pw_per_k": 20.},
                      {"cutoff_A": 10., "g_ph_pw_per_k": 10.},
                      {"cutoff_A": 12., "g_ph_pw_per_k": 10.}]}
    axis = convergence_axes(built, block)["fc_cutoff"]
    assert axis["nominal_shortfall"] == 0.
    assert axis["converged"] is True
    block["rows"][1]["g_ph_pw_per_k"] = 12.
    axis = convergence_axes(built, block)["fc_cutoff"]
    assert np.isclose(axis["nominal_shortfall"], 0.2)
    assert axis["converged"] is False


def test_curve_evaluates_frequency_batches_and_agrees_with_single_points(monkeypatch):
    from smj.stages import thermal_analyze as stage

    built = _chain_built()
    built["centre"] = stage._dof_indices([2, 3])
    omega = np.linspace(.001, .035, 16)
    matrix = built["dynamical"]
    reference = [stage._curve_factory(built, np.array([w]), .0005)(matrix)[0]
                 for w in omega]
    sizes = []
    original = stage.transmission

    def measured(w, *args, **kwargs):
        sizes.append(np.size(w))
        return original(w, *args, **kwargs)

    monkeypatch.setattr(stage, "transmission", measured)
    spectrum = stage._curve_factory(built, omega, .0005)(matrix)
    assert max(sizes) > 1
    np.testing.assert_allclose(spectrum, reference, atol=1e-9, rtol=1e-8)
