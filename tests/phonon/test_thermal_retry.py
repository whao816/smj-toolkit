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
        "slots": [{"slot": i, "converged": i != 2} for i in range(slots)],
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


def test_read_forces_prefers_retrieved_bundle_over_old_local_slots(workdir):
    first = [[0.1, 0.0, 0.0], [0.0, 0.2, 0.0]]
    second = [[0.3, 0.0, 0.0], [0.0, 0.4, 0.0]]
    _make_slot(workdir, 0, forces=first, aims=_AIMS_OK)
    _make_slot(workdir, 1, forces=second, aims=_AIMS_OK)
    # 远端重试只回取汇总文件；本地 fc/ 仍是旧数据。
    (workdir / "forces_all.json").write_text(json.dumps({
        "n_atoms": 2, "count": 2,
        "forces": np.zeros((2, 2, 3)).tolist()}),
        encoding="utf-8")
    stack = ThermalFcStage._read_forces(workdir, 2, 2)
    assert np.all(stack == 0.0)


def test_new_batch_without_scf_log_is_retryable(workdir):
    _make_slot(workdir, 0, forces=[[0, 0, 0]])
    out = workdir / "forces_all.json"
    summary = workdir / "fc_summary.json"
    assert concat_forces.main([str(workdir / "fc"), "1", str(out),
                               "--summary", str(summary), "--allow-missing"]) == 1
    assert json.loads(summary.read_text(encoding="utf-8"))["unknown"] == [0]
    with pytest.raises(StageError, match="--retry 0"):
        ThermalFcStage._read_forces(workdir, 1, 1)


@pytest.mark.parametrize("entries", [
    [],
    [{"slot": 0, "converged": None}],
    [{"slot": 0, "converged": "true"}],
    [{"slot": 0, "converged": True}, {"slot": 0, "converged": True}],
])
def test_incomplete_or_malformed_scf_evidence_is_rejected(workdir, entries):
    _make_slot(workdir, 0, forces=[[0, 0, 0]])
    (workdir / "fc_summary.json").write_text(json.dumps({
        "count": 1, "missing": [], "not_converged": [], "slots": entries,
    }), encoding="utf-8")
    with pytest.raises(StageError, match="fc_summary|--retry 0"):
        ThermalFcStage._read_forces(workdir, 1, 1)


@pytest.mark.parametrize("with_summary", [True, False])
def test_collect_reports_evidence_or_legacy_unknown(workdir, with_summary):
    from smj.phonon.fc import make_displacements

    disps = make_displacements(1, [0], 0.01)
    (workdir / "geometry.in").write_text("atom 0 0 0 Au\n", encoding="utf-8")
    (workdir / "displacements.json").write_text(json.dumps(disps), encoding="utf-8")
    for slot in range(len(disps)):
        _make_slot(workdir, slot, forces=[[0, 0, 0]], aims=_AIMS_OK)
    if with_summary:
        concat_forces.main([str(workdir / "fc"), str(len(disps)),
                            str(workdir / "forces_all.json"),
                            "--summary", str(workdir / "fc_summary.json")])
    result = ThermalFcStage().collect(None, workdir).results["thermal_fc"]
    assert result["scf_converged"] is (True if with_summary else None)


def test_multislot_retry_cleans_only_requested_results_before_submit(workdir):
    from types import SimpleNamespace
    from smj.remote import RemoteSession
    from smj.testing import FakeTransport
    from smj.transport import CommandResult

    class Cluster(FakeTransport):
        def on_run(self, command, cwd, timeout):
            if "&& sbatch" in command:
                return CommandResult(0, "Submitted batch job 77")

    cluster = Cluster()
    ctx = SimpleNamespace(remote_stage_dir=lambda stage: "/remote/fc")
    session = RemoteSession(ctx, transport=cluster)
    script = "sub_thermal_fc.sh"
    for name in (script, "forces_all.json", "fc_summary.json", "force_constants.npy"):
        (workdir / name).write_text("previous", encoding="utf-8")
    sub = session.submit("thermal_fc", workdir, script=script,
                         sbatch_args=("--array=0,2%2",))
    cleanup, submit = [cmd for kind, cmd in cluster.commands if kind == "run"]
    assert "--array=0,2%2" in submit
    for result in (".done", "forces.json", "aims.out"):
        assert f"fc/0000/{result}" in cleanup and f"fc/0002/{result}" in cleanup
        assert f"fc/0001/{result}" not in cleanup and f"fc/*/{result}" not in cleanup
    for name in ("forces_all.json", "fc_summary.json"):
        assert name in cleanup
    assert sub.uploaded == [script]


def test_read_forces_falls_back_to_bundle_when_slots_absent(workdir):
    stack = np.arange(12, dtype=float).reshape(2, 2, 3)
    (workdir / "forces_all.json").write_text(json.dumps({
        "n_atoms": 2, "count": 2, "forces": stack.tolist()}), encoding="utf-8")
    read = ThermalFcStage._read_forces(workdir, 2, 2)
    assert np.array_equal(read, stack)


def test_retry_slots_accept_explicit_list_and_keywords(workdir):
    from smj.cli import _resolve_retry_slots

    (workdir / "fc_summary.json").write_text(json.dumps({
        "count": 10, "missing": [3, 7], "not_converged": [9], "unknown": [6],
        "slots": [],
    }), encoding="utf-8")
    assert _resolve_retry_slots("17, 203", workdir) == [17, 203]
    assert _resolve_retry_slots("missing", workdir) == [3, 7]
    assert _resolve_retry_slots("not-converged", workdir) == [9]
    assert _resolve_retry_slots("unknown", workdir) == [6]


def test_retry_does_not_overwrite_an_active_array(workdir):
    from types import SimpleNamespace
    from smj.cli import _cmd_thermal_retry
    from smj.errors import SmjError

    (workdir / "displacements.json").write_text("[{}, {}]", encoding="utf-8")
    runner = SimpleNamespace(layout=SimpleNamespace(stage_dir=lambda name: workdir),
                             state=SimpleNamespace(status=lambda name: "waiting"))
    config = SimpleNamespace(mode="slurm-ssh", thermal=SimpleNamespace(fc=None))
    with pytest.raises(SmjError, match="--wait"):
        _cmd_thermal_retry(SimpleNamespace(retry="1", force=False),
                           config, runner, ("thermal_fc",))


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


def test_array_publishes_only_after_all_expected_retry_slots(workdir):
    import os
    from pathlib import Path
    import shlex
    import shutil
    import subprocess
    import sys

    from smj.phonon.fc import make_displacements
    from smj.stages._thermal_common import template_path

    bash = shutil.which("bash")
    if os.name == "nt" and shutil.which("git"):
        bundled = Path(shutil.which("git")).parent.parent / "bin" / "bash.exe"
        if bundled.is_file():
            bash = str(bundled)
    if not bash:
        pytest.skip("requires bash to execute the array script")
    stage = workdir.resolve()
    (stage / "control.in").write_text("", encoding="utf-8")
    (stage / "geometry.in").write_text("atom 0 0 0 Au\n", encoding="utf-8")
    (stage / "displacements.json").write_text(json.dumps(
        make_displacements(1, [0], 0.01)[:3]), encoding="utf-8")
    for name in ("concat_forces.py", "extract_forces.py", "apply_displacement.py"):
        shutil.copy(Path(concat_forces.__file__).with_name(name), stage / name)
    # 槽 1 有效，槽 8 是旧批次残留，不能用它顶替尚未结束的槽 2。
    for slot in (1, 8):
        _make_slot(stage, slot, forces=[[1, 0, 0]], aims=_AIMS_OK)
        (stage / "fc" / f"{slot:04d}" / ".done").write_text("0 old\n", encoding="utf-8")
    template = template_path("sub_thermal_fc.sh").read_text(encoding="utf-8")
    script = template.format(
        partition="test", ntasks=1, mem="1G", walltime="0:01:00", gres_line="",
        last=2, count=3, parallel=2, exe_define="",
        module_load="python() { " + shlex.quote(Path(sys.executable).as_posix()) + ' "$@"; }',
        aims_bin="printf 'Self-consistency cycle converged.\\nTotal atomic forces\\n"
                 "| 1 1 0 0\\n' > aims.out")
    (stage / "array.sh").write_text(script, encoding="utf-8", newline="\n")
    for slot in (0, 2):
        env = dict(os.environ, SLURM_SUBMIT_DIR=stage.as_posix(),
                   SLURM_ARRAY_TASK_ID=str(slot), SLURM_ARRAY_JOB_ID="771",
                   SLURM_ARRAY_TASK_COUNT="2")
        subprocess.run([bash, "array.sh"], cwd=stage, env=env, check=True)
        assert (stage / ".dsh_done").exists() == (slot == 2)
        assert (stage / "forces_all.json").exists() == (slot == 2)
    assert (stage / ".dsh_done").read_text().startswith("0 ")
    # Standalone sbatch --array=0,2 starts with the previous bundle and markers.
    for slot in (2, 0):
        env.update(SLURM_ARRAY_TASK_ID=str(slot), SLURM_ARRAY_JOB_ID="772")
        subprocess.run([bash, "array.sh"], cwd=stage, env=env, check=True)
        assert (stage / ".dsh_done").exists() == (slot == 0)
        assert (stage / "forces_all.json").exists() == (slot == 0)
    assert (stage / "fc" / "0001" / ".done").read_text() == "0 old\n"
