"""validate 的产物定位不能把弛豫日志当成静态 SCF 日志（2026-11-20 复核 R13b）。

`_ARTIFACT_STAGE_HOME['aims.out']` 原为 ``('aims_scf', 'aims_tip_relax')``：只要
``aims_scf/aims.out`` 不在，validate 会去读 ``aims_tip_relax/aims.out``。弛豫日志
每一步都有 ``Self-consistency cycle converged.``，因此"静态 SCF 没跑"的 run 也能
拿到一条 aims_scf=OK 的判据。修法是只认 ``aims_scf`` 阶段（扁平布局仍由
``run_dir/aims.out`` 命中）。
"""

from __future__ import annotations

from pathlib import Path

from smj.junction.validate import validate_junction_run

RELAX_LOG = (
    "Self-consistency cycle converged.\n"
    "Present geometry is converged.\n"
    "Have a nice day.\n"
)
SCF_LOG = "Self-consistency cycle converged.\nHave a nice day.\n"


def _aims_checks(verdict: dict) -> list[dict]:
    return [check for check in verdict['checks'] if check['name'] == 'aims_scf']


def test_tip_relax_log_does_not_satisfy_the_scf_check(workdir: Path):
    (workdir / 'aims_tip_relax').mkdir()
    (workdir / 'aims_tip_relax' / 'aims.out').write_text(RELAX_LOG, encoding='utf-8')

    verdict = validate_junction_run(workdir)
    assert _aims_checks(verdict) == [], (
        'aims_scf/aims.out 不存在时不应从弛豫日志拿到 SCF 判据: %s'
        % _aims_checks(verdict)
    )


def test_scf_stage_log_still_checked(workdir: Path):
    (workdir / 'aims_scf').mkdir()
    (workdir / 'aims_scf' / 'aims.out').write_text(SCF_LOG, encoding='utf-8')
    (workdir / 'aims_tip_relax').mkdir()
    (workdir / 'aims_tip_relax' / 'aims.out').write_text(RELAX_LOG, encoding='utf-8')

    verdict = validate_junction_run(workdir)
    checks = _aims_checks(verdict)
    assert len(checks) == 1 and checks[0]['passed'] is True
    # 证据必须来自 aims_scf 那一份
    assert '无失败标记' in checks[0]['detail']


def test_scf_stage_failure_markers_still_rejected(workdir: Path):
    (workdir / 'aims_scf').mkdir()
    (workdir / 'aims_scf' / 'aims.out').write_text(
        "SCF is not yet converged\nHave a nice day.\n", encoding='utf-8')
    verdict = validate_junction_run(workdir)
    checks = _aims_checks(verdict)
    assert len(checks) == 1 and checks[0]['passed'] is False
    assert checks[0]['severity'] == 'error'


def test_flat_layout_aims_out_still_found(workdir: Path):
    """扁平布局（产物直接放在 run 根）不能被这次收紧误伤。"""
    (workdir / 'aims.out').write_text(SCF_LOG, encoding='utf-8')
    checks = _aims_checks(validate_junction_run(workdir))
    assert len(checks) == 1 and checks[0]['passed'] is True
