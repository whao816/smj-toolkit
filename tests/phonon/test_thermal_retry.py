"""数组作业的"失败可恢复"路径：收敛校验、缺槽位、单点重试。

为什么要单独一组测试：324 个位移点跑 9 小时，失败一两个点是常态，而失败有两种
形态且都**不会**被形状检查抓到：

1. 槽位没产出 ``forces.json``（作业被杀、写盘失败）；
2. 槽位产出了力，但 SCF 没收敛——数看着完全正常，只是不该用。

这两种都只能在远端从 ``aims.out`` 读出来，所以检查被放在数组作业的收尾步骤里，
本地 collect 只消费 ``fc_summary.json``。
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from smj.errors import StageError
from smj.phonon import concat_forces
from smj.stages.thermal_fc import ThermalFcStage

_AIMS_OK = """\
  Self-consistency cycle converged.
  Have a nice day.
"""

_AIMS_BAD = """\
  Self-consistency cycle NOT converged.
  Have a nice day.
"""

_AIMS_TWO_CYCLES = """\
  Self-consistency cycle 1 converged.
  Self-consistency cycle converged.
  Have a nice day.
"""


def test_scan_aims_distinguishes_converged_from_not():
    ok = concat_forces.scan_aims(_AIMS_OK)
    assert ok["converged"] is True
    assert ok["n_cycles"] == 1

    bad = concat_forces.scan_aims(_AIMS_BAD)
    assert bad["converged"] is False

    # "NOT converged" 里也含 "converged" 子串:两条标记必须分开判
    assert concat_forces.scan_aims(_AIMS_TWO_CYCLES)["n_cycles"] == 2


def _make_slot(root, slot, *, forces=None, aims=None):
    directory = root / "fc" / ("%04d" % slot)
    directory.mkdir(parents=True)
    if forces is not None:
        (directory / "forces.json").write_text(
            json.dumps({"n_atoms": len(forces), "forces": forces}),
            encoding="utf-8")
    if aims is not None:
        (directory / "aims.out").write_text(aims, encoding="utf-8")


def test_concat_with_missing_slots_records_them_instead_of_dying(workdir):
    """缺槽位时也要落盘汇总:否则一次失败就等于"没有任何可用产物"。"""
    slot_forces = [[[0.1, 0.0, 0.0], [0.0, 0.2, 0.0]]]
    _make_slot(workdir, 0, forces=slot_forces, aims=_AIMS_OK)
    _make_slot(workdir, 1, aims=_AIMS_OK)  # 有 aims.out,没 forces.json
    payload = concat_forces.concat(str(workdir / "fc"), 2,
                                   allow_missing=True, with_summary=True)
    assert payload["forces"][0] == slot_forces
    assert payload["forces"][1] is None
    summary = payload["summary"]
    assert summary["missing"] == [1]
    assert summary["not_converged"] == []


def test_concat_reports_not_converged_slots(workdir):
    slot_forces = [[[0.1, 0.0, 0.0], [0.0, 0.2, 0.0]]]
    _make_slot(workdir, 0, forces=slot_forces, aims=_AIMS_OK)
    _make_slot(workdir, 1, forces=slot_forces, aims=_AIMS_BAD)
    payload = concat_forces.concat(str(workdir / "fc"), 2,
                                   allow_missing=True, with_summary=True)
    assert payload["summary"]["not_converged"] == [1]
    assert payload["summary"]["missing"] == []


def test_concat_cli_exit_code_flags_incomplete_runs(workdir):
    slot_forces = [[[0.1, 0.0, 0.0], [0.0, 0.2, 0.0]]]
    _make_slot(workdir, 0, forces=slot_forces, aims=_AIMS_OK)
    out = workdir / "forces_all.json"
    summary = workdir / "fc_summary.json"
    code = concat_forces.main([
        str(workdir / "fc"), "2", str(out),
        "--summary", str(summary), "--allow-missing"])
    assert code == 1  # 不完整 → 非零,让 .dsh_done 记成失败
    assert json.loads(summary.read_text(encoding="utf-8"))["missing"] == [1]


@pytest.mark.parametrize("slots", [3])
def test_read_forces_rejects_unconverged_slots(workdir, slots):
    """未收敛的槽位即使力齐全、形状正确,也不能进力常数。"""
    n_atoms = 2
    good = [[0.1, 0.0, 0.0], [0.0, 0.2, 0.0]]
    for slot in range(slots):
        _make_slot(workdir, slot, forces=[good] * 1, aims=_AIMS_OK)
    (workdir / "fc_summary.json").write_text(json.dumps({
        "count": slots, "missing": [], "not_converged": [2],
        "slots": [{"slot": 2, "converged": False}],
    }), encoding="utf-8")
    with pytest.raises(StageError) as excinfo:
        ThermalFcStage._read_forces(workdir, slots, n_atoms)
    message = str(excinfo.value)
    assert "未收敛" in message
    assert "--retry 2" in message


def test_read_forces_points_at_retry_for_missing_slots(workdir):
    good = [[0.1, 0.0, 0.0], [0.0, 0.2, 0.0]]
    _make_slot(workdir, 0, forces=good, aims=_AIMS_OK)
    _make_slot(workdir, 1, forces=good, aims=_AIMS_OK)
    _make_slot(workdir, 2, aims=_AIMS_OK)
    with pytest.raises(StageError) as excinfo:
        ThermalFcStage._read_forces(workdir, 3, 2)
    message = str(excinfo.value)
    assert "缺 1 个位移点" in message
    assert "smj thermal --stage thermal_fc --retry 2" in message


def test_read_forces_prefers_individual_slots_over_stale_bundle(workdir):
    """重试后本地逐个文件比汇总文件新:必须以逐个为准。"""
    first = [[0.1, 0.0, 0.0], [0.0, 0.2, 0.0]]
    second = [[0.3, 0.0, 0.0], [0.0, 0.4, 0.0]]
    _make_slot(workdir, 0, forces=first, aims=_AIMS_OK)
    _make_slot(workdir, 1, forces=second, aims=_AIMS_OK)
    # 汇总文件里是旧的(全零)数据
    (workdir / "forces_all.json").write_text(json.dumps({
        "n_atoms": 2, "count": 2,
        "forces": [[[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]],
                   [[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]]]}),
        encoding="utf-8")
    stack = ThermalFcStage._read_forces(workdir, 2, 2)
    assert stack[0][0][0] == pytest.approx(0.1)
    assert stack[1][0][0] == pytest.approx(0.3)


def test_read_forces_falls_back_to_bundle_when_slots_absent(workdir):
    stack = np.arange(12, dtype=float).reshape(2, 2, 3)
    (workdir / "forces_all.json").write_text(json.dumps({
        "n_atoms": 2, "count": 2, "forces": stack.tolist()}), encoding="utf-8")
    read = ThermalFcStage._read_forces(workdir, 2, 2)
    assert np.array_equal(read, stack)


def test_retry_slots_accept_explicit_list_and_keywords(workdir):
    from smj.cli import _resolve_retry_slots

    (workdir / "fc_summary.json").write_text(json.dumps({
        "count": 10, "missing": [3, 7], "not_converged": [9], "slots": [],
    }), encoding="utf-8")
    assert _resolve_retry_slots("17, 203", workdir) == [17, 203]
    assert _resolve_retry_slots("missing", workdir) == [3, 7]
    assert _resolve_retry_slots("not-converged", workdir) == [9]


def test_retry_slots_reject_bad_input(workdir):
    from smj.cli import _resolve_retry_slots

    with pytest.raises(Exception) as excinfo:
        _resolve_retry_slots("17,x", workdir)
    assert "整数" in str(excinfo.value)

    # 关键词需要 fc_summary.json;没有它时给出"直接写下标"的出路
    with pytest.raises(Exception) as excinfo:
        _resolve_retry_slots("missing", workdir)
    assert "fc_summary.json" in str(excinfo.value)


def test_array_script_summarises_convergence_at_finalize():
    """收尾命令必须带 --summary/--allow-missing,否则失败时没有可读产物。"""
    from smj.stages._thermal_common import template_path

    text = template_path("sub_thermal_fc.sh").read_text(encoding="utf-8")
    assert "concat_forces.py ../ {count} ../../forces_all.json" in text
    assert "--summary ../../fc_summary.json" in text
    assert "--allow-missing" in text
