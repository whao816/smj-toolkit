"""P2-5 (AGD-2026-09-25-SMJ): golden fixtures —— 把真跑冻结为可复算的回归基准。

动机: 本会话为验证改动反复重跑/手算; 有了冻结基准, 后续改管线只需 compare 一次即可
回答"有没有退化"。这是 P0-1/P0-3/P0-4/P1-x 各项改动的**共同安全网**。

数值来源(从**本地既有产物**实测, 非凭记忆):
- bdt-trans-smj-aae585 (opt=tight, 最终轮) 与 bdt-trans-smj-624059 (普通 opt) 均在
  work/smj_bdt/<run_id>/ 下, 5 阶段齐全, TE.dat 4001 点。
- 2cn-smj-g16t01(2026-10-14 更新): 该 run **已在本地** work/2cn-smj-g16t01/
  (5 阶段齐全 + TE.dat 4001 点), 现冻结为可本地复算的基准 —— 早前"仅在集群"的记录已作废。
  H₂O→Au60 体系按用户决定**不再作为验收体系**(其 harness 已删), 验收面由 2cn(主场)与
  bdt(基准对照)承担。

注: 会话叙述里 624059 曾记作 G=0.8611, 而本地产物现算 L0=0.861257(差 1.5e-4);
与已登记的"L0 第 7 位"同类, **以本地产物现算值为准**并留此说明。

注(2026-10-24 第二刀): G 理论带/Seebeck 判据登记了 applies_to(仅 Au-BDT-Au), 故 2cn
基准中这两条记 not_applicable; 2026-11-15 判据更新(t_ef_magnitude 取消、共振归因口径)
后 2cn 的 physics_verdict 为 passed。bdt 两条基准仍可复算(同族, 判据适用;
g_in_theory_band 现为 warn/派生包络, 故 inconclusive)。历史 run 的 report/DB 记录不重写。
"""
from __future__ import annotations

from pathlib import Path

FIXTURES = {
    'bdt-trans-smj-aae585': {
        'available_locally': True,
        'note': 'opt=tight 最终轮(本会话结论所依据的那一轮)',
        'te_sha256': 'aa961b1977e30aeb78f0c306cfa3d40f3a1f7b2502e4cd6163e269ab2c4249ab',
        'te_rows': 4001,
        't_ef': 0.8491, 'l0': 0.860733, 's_uV_per_K': 37.1507,
        'physics_verdict': 'inconclusive',
        'failed_criteria': [],
        'inconclusive_criteria': ['g_in_theory_band'],
        'stages': 5,
    },
    'bdt-trans-smj-624059': {
        'available_locally': True,
        'note': '普通 opt(因 Gaussian 虚频 -0.5188 cm^-1 而作废, 保留作对照)',
        'te_sha256': 'eac812ccf0f721dca53f8baa582304a92f1551aa43b85cef8d2131b3c35665c1',
        'te_rows': 4001,
        't_ef': 0.8496, 'l0': 0.861257, 's_uV_per_K': 37.1275,
        'physics_verdict': 'inconclusive',
        'failed_criteria': [],
        'inconclusive_criteria': ['g_in_theory_band'],
        'stages': 5,
    },
    '2cn-smj-g16t01': {
        'available_locally': True,
        'run_root': 'work/2cn-smj-g16t01',
        'system': '2CN',        # 判据适用域: 非 BDT 族 -> seebeck_sign 不作闸门(2026-10-24)
        'note': ('主场真实算例(2CN-656-4-thiomethylphenyl 54 原子 + au10 电极 -> 74 原子; '
                 '与 bdt 族分属不同 run 根, 故用 run_root 显式定位)。**阶段目录名为 au10_build**'
                 '(小电极变体)而非 au60_build, 故用 stage_names 显式声明, 不按 SMJ_STAGES 数数。'
                 '**2026-10-24 第二刀(判据适用域)**: seebeck_sign 与 G 理论带都按 Au-BDT-Au'
                 '标定(Reddy 2007 正号 / Thygesen 2005 派生带), 对 2cn 属**判据适用性**问题而非'
                 '物理错误 —— 由 criteria.applies_to 表达, 两条判据对 2CN 记 not_applicable、'
                 '不计入 failed/inconclusive。**2026-11-15(判据更新)**: t_ef_magnitude 取消'
                 '(无一手区间); t_ef_no_peak_near_ef 改为 t_ef_resonance_attribution —— 2CN 未'
                 '登记族, 按实测 S=-102.9 uV/K 自洽归入 LUMO 侧, 峰 +0.241 eV 落在 0.15-1.0 eV '
                 '-> pass; 于是 physics_verdict 由 inconclusive 变为 passed。历史 run 的 report/DB '
                 '记录不重写, 本 fixture 只冻结"当前口径下"的可复算基准'),
        'te_sha256': 'd9616f307f56cf51542cd53bc9af92d2941a44d0aa75c43c477a80e585d269e4',
        'te_rows': 4001,
        't_ef': 0.01863583574, 'l0': 0.0813315561957088, 's_uV_per_K': -102.91458413897126,
        'physics_verdict': 'passed',
        'failed_criteria': [],
        'inconclusive_criteria': [],
        'not_applicable_criteria': ['g_in_theory_band', 'seebeck_sign'],
        'stages': 5,
        'stage_names': ['gaussian_opt', 'aims_tip_relax', 'au10_build', 'aims_scf', 'aitranss_te'],
    },
}

#: 仓库根(fixtures.py -> junction -> tools -> agenticdft(内) -> agenticdft(包) -> 仓库根)
_REPO = Path(__file__).resolve().parents[4]

_TOL = {'t_ef': 5e-4, 'l0': 5e-6, 's_uV_per_K': 5e-3}


def verify_fixture(run_root, run_id: str, tol: dict | None = None) -> dict:
    """复算并比对冻结基准 -> {'ok', 'diffs', 'missing', 'skipped'}。"""
    tol = {**_TOL, **(tol or {})}
    spec = FIXTURES.get(run_id)
    if spec is None:
        return {'ok': False, 'diffs': [], 'missing': ['未登记的 fixture: %s' % run_id],
                'skipped': False}
    if not spec.get('available_locally', True):
        return {'ok': True, 'diffs': [], 'missing': [], 'skipped': True,
                'reason': '该 fixture 不可本地校验: %s' % spec.get('note')}
    root = Path(run_root)
    te = root / 'aitranss_te' / 'TE.dat'
    if not te.exists() and spec.get('run_root'):        # fixture 可自带 run 根(如 2cn 不在 bdt 族下)
        root = _REPO / spec['run_root']
        te = root / 'aitranss_te' / 'TE.dat'
    if not te.exists():
        return {'ok': False, 'diffs': [], 'missing': ['TE.dat'], 'skipped': False}
    from .landauer import analyze_te
    from .manifest import build_manifest, sha256_file
    diffs = []
    if sha256_file(te) != spec['te_sha256']:
        diffs.append('TE.dat sha256 不符(产物已变; 若为有意变更请更新 fixture)')
    rows = [l for l in te.read_text(encoding='utf-8', errors='ignore').splitlines()
            if l.strip() and not l.lstrip().startswith('#')]
    if len(rows) != spec['te_rows']:
        diffs.append('TE.dat 行数 %d != %d' % (len(rows), spec['te_rows']))
    a = analyze_te(te)
    pairs = (('t_ef', a['fermi_nearest']['T']), ('l0', a['L0']),
             ('s_uV_per_K', a['S_uV_per_K']))
    for k, got in pairs:
        if abs(got - spec[k]) > tol[k]:
            diffs.append('%s: 现算 %.6g != 冻结 %.6g (tol %g)' % (k, got, spec[k], tol[k]))
    man = build_manifest(root, stages=spec.get('stage_names'),    # 小电极变体用 au10_build 等异名阶段
                         system=spec.get('system'))               # 判据适用域(2CN != BDT 族)
    v = man.get('verdict') or {}
    if v.get('physics_verdict') != spec['physics_verdict']:
        diffs.append('physics_verdict: %s != %s' % (v.get('physics_verdict'), spec['physics_verdict']))
    if list(v.get('failed_criteria') or []) != spec['failed_criteria']:
        diffs.append('failed_criteria: %s != %s' % (v.get('failed_criteria'), spec['failed_criteria']))
    if 'not_applicable_criteria' in spec and \
            list(v.get('not_applicable_criteria') or []) != spec['not_applicable_criteria']:
        diffs.append('not_applicable_criteria: %s != %s'
                     % (v.get('not_applicable_criteria'), spec['not_applicable_criteria']))
    if man['totals']['stages'] != spec['stages']:
        diffs.append('阶段数 %s != %s' % (man['totals']['stages'], spec['stages']))
    return {'ok': not diffs, 'diffs': diffs, 'missing': [], 'skipped': False}
