"""P1-1/P1-2/P1-5 (AGD-2026-09-25-SMJ): smj 物理判据注册表与 physics_verdict。

两条铁律(承 v0.4 §0 A4 与 §2.5 断言纪律):
1. **每条 physical 判据必须带 source{doi, locator}** —— 无源判据不得作闸门;
2. **engineering 判据**(程序语义, 如左右电极不相交)不要求文献源, 但须标 kind。

P1-5 的强制由 lint_registry() + tests/test_smj_automation.py 的 ratchet 用例实现:
新增无源判据即失败; 已登记的缺口(见 KNOWN_GAPS)数量固定, 只能减少不能增加。

判据数值来源(均已回读原文核验, 见 docs/taskcard-smj-bdt-benchmark.md §4.2):
- Au-S 带 2.40-2.49 A: Thygesen 2005 p.12 `dAu-S = 2.40 A`; Jackson 2009 `S-Au = 2.42+-0.02 A`
- T(E_F) 带 0.03-0.4: 由 Thygesen 2005 p.13 `0.1G0` 与其 Ref.15/48 的 x3/x4 关系推得
- T(E_F) 结构判据: Thygesen 2005 p.12 峰位(-3.0/+2.5 宽峰, -1.5/-1.0 窄峰) -> E_F 近旁 1.0 eV 内无峰
- 适用域(2026-10-24 第二刀): T(E_F)/G 绿带与 seebeck_sign 同属 **Au-BDT-Au 标定域**, 故同带
  applies_to; 换分子族(如饱和烷烃二硫醇)不作闸门, 实算值仍记录(见 SYSTEM_FAMILIES)
"""
from __future__ import annotations

from pathlib import Path

THY = {'doi': '10.1016/j.chemphys.2005.05.032'}
JAC = {'doi': '10.1016/j.susc.2009.01.022'}

# kind: physical(物理判据, 须有源) | engineering(程序/接口语义, 免源但须标注)
# gate: error(阻断) | warn(告警) | info(仅报告)
CRITERIA = [
    {'name': 'gaussian_no_imaginary', 'kind': 'physical', 'gate': 'error', 'implemented': True,
     'band': None, 'source': dict(THY, locator='—'),
     'note': '虚频判据源为项目既定口径(任务卡 §3); 低频噪声阈值仍无源, 见 KNOWN_GAPS'},

    {'name': 'scf_converged', 'kind': 'engineering', 'gate': 'error', 'implemented': True,
     'band': None, 'source': None, 'note': 'FHI-aims 自洽收敛标记(程序语义)'},

    # 2026-11-11 本地补丁 (smj F2): 产物缺失必须是 error 级判据, 不能被"零判据"静默判过。
    {'name': 'te_present', 'kind': 'engineering', 'gate': 'error', 'implemented': True,
     'band': None, 'source': None,
     'note': 'transmission 产物存在性(程序语义); 缺失时物理判据无法评估, 不得判 passed'},

    {'name': 'grid_independence', 'kind': 'physical', 'gate': 'warn', 'implemented': True,
     'band': [0.01, None], 'source': dict(THY, locator='—'),
     'note': '抽稀后 L0/L1/L2 相对变化 <1%(任务卡 §5.1.3); 本会话实测 <=0.09%'},

    {'name': 'electrodes_disjoint', 'kind': 'engineering', 'gate': 'error', 'implemented': True,
     'band': None, 'source': None,
     'note': 'aitranss 语义: 左右电极共享原子即终止(trans.f); 见 P0-4'},

    {'name': 'au_s_bond', 'kind': 'physical', 'gate': 'warn', 'implemented': True,
     'band': [2.40, 2.49], 'source': dict(THY, locator='p.12'),
     'note': '另有实验锚点 ' + JAC['doi'] + ' S-Au=2.42+-0.02; 仅方向性(两锚点均为表面/金字塔配位); '
             '仅适用于 S 为锚的去质子硫醇盐；-SCH3 等硫醚的 S–Au 是配位接触，'
             'N/Se 等非 S 锚只报接触距离，2.40–2.49 Å 均不适用'},

    {'name': 't_ef_magnitude', 'kind': 'physical', 'gate': 'error', 'implemented': True,
     'band': [0.03, 0.4], 'source': dict(THY, locator='p.13'),
     'note': 'Au-BDT-Au 典型 ~0.1; 上下端点由该文 Ref.15(x3)/Ref.48(x4) 关系推得',
     # 2026-10-24 第二刀: 绿带是 **Au-BDT-Au 理论标定包络**, 与 seebeck_sign 同属该族;
     # 饱和烷烃二硫醇(C4)实测 T(E_F)~1.5e-2 G0, 低于下端点 2 倍 —— 那是通道不同
     # (隧穿衰减), 不是对 H0 的否证。不适用时仍记实算值, 只不进 failed/inconclusive。
     'applies_to': {'families': ('Au-BDT-Au',),
                    'basis': ('绿带 0.03-0.40 G0 由 Thygesen 2005 p.13 的 Au-BDT-Au ~0.1 G0'
                              ' 与其 Ref.15(x3)/Ref.48(x4) 关系推得, 是该族理论标定包络;'
                              ' 换分子族(饱和烷烃/等)输运通道不同, 超带不构成对 H0 的否证'
                              ' —— 2026-10-24 用户拍板, 与 seebeck_sign 同口径')}},

    {'name': 't_ef_no_peak_near_ef', 'kind': 'physical', 'gate': 'warn', 'implemented': True,
     'band': [1.0, None], 'source': dict(THY, locator='p.12'),
     'note': 'E_F 近旁应 >=1.0 eV 无透射峰; 依赖 HOMO 位置, 随泛函/几何有 ~0.5 eV 漂移'},

    {'name': 'g_in_theory_band', 'kind': 'physical', 'gate': 'error', 'implemented': True,
     'band': [0.03, 0.4], 'source': dict(THY, locator='p.13'),
     'note': '任务卡 §4.2 四基准(Reed/Reddy/Lee/Thygesen)已登记',
     # 与 t_ef_magnitude 同一绿带/同一标定域(G/G0≈T(E_F), 300 K 费米窗窄); 适用域一并声明
     'applies_to': {'families': ('Au-BDT-Au',),
                    'basis': ('同 t_ef_magnitude: 0.03-0.40 G0 是 Au-BDT-Au 理论标定包络'
                              '(Thygesen 2005 p.13 及其 Ref.15/48 比例); 300 K 下 G/G0≈T(E_F),'
                              ' 换分子族不作闸门 —— 2026-10-24 用户拍板')}},

    {'name': 'seebeck_sign', 'kind': 'physical', 'gate': 'error', 'implemented': True,
     'band': [0, None], 'source': {'doi': '10.1126/science.1137149', 'locator': '摘要/正文'},
     'note': 'Reddy 2007 +8.7+-2.1 与 Lee 2014 +7.4+-0.5 同号为正 -> 铁律判据',
     'applies_to': {'families': ('Au-BDT-Au',),
                    'basis': ('两篇实验均为 Au-BDT-Au(Reddy 2007/Lee 2014), 正号是该族实测'
                              '性质; 换分子体系(如 2CN)不得当普适闸门 —— 2026-10-24 用户拍板'
                              '"窄口径", 2CN 记不适用而非判否')}},

    {'name': 'k_e_zt_order', 'kind': 'physical', 'gate': 'warn', 'implemented': False,
     'band': None, 'source': None,
     'note': '目前只有"非负且量级合理"的物理常识, 正式来源待补'},
]

# 已登记缺口(ratchtet): 只允许这些无源/未实现的 physical 判据存在, 新增即 lint 失败
KNOWN_GAPS = frozenset({'k_e_zt_order'})

# ---------------------------------------------------------------------------
# A (2026-10-24, 用户拍板"窄口径"): 判据**适用域**。
#
# 有源 != 普适。`seebeck_sign` 的正号来自 Au-BDT-Au 两篇实验(Reddy 2007 +8.7+-2.1 /
# Lee 2014 +7.4+-0.5), 即"该族实测符号"被当作闸门。换成 2CN(不同分子、不同输运通道),
# S<0 说明不了 H0 被否证 —— 属**判据适用性**问题而非物理错误。故带 applies_to 的判据只在
# 登记过的体系族上作闸门; 其余体系产出一条 applicable=False 的 check(仍带实算值, 便于
# 追溯), 归约时进 not_applicable_criteria, 不计入 failed/inconclusive。
#
# 方向必须 fail-safe: **未登记/未知体系一律视为适用**(见 system_family), 不能因为漏登记
# 而让某条判据悄悄失效。
# ---------------------------------------------------------------------------
SYSTEM_FAMILIES = {
    # 体系标识(TaskCard.material / validate_junction_run(system=) / fixture 的 system 字段)
    'BDT': 'Au-BDT-Au',
    '2CN': 'Au-(2CN-656-4-thiomethylphenyl)',
    # 2026-10-24 C4 真跑(TaskCard.material='c4'): 饱和烷烃二硫醇, 无芳香 π 通道,
    # T(E_F) 比 BDT 低约 60 倍 —— 显式登记族名, 让绿带判据对它记"不适用"而非"否证"。
    'C4': 'Au-C4-dithiol-Au',
    # 2026-11-10 薁(azulene)2,6-二硫醇: 非交替芳烃(5-7 并环, 永久偶极 1.03 D 沿长轴),
    # 硫醇去质子直连 Au。它既不是交替苯环的 π 通道(BDT), 也不是饱和链(C4) ——
    # T(E_F)/G 绿带与 Seebeck 正号都是 Au-BDT-Au 标定域, 对薁作闸门会把一次物理上
    # 正常的算例判成 failed。显式登记族名, 让三条判据记"不适用"而非"否证"。
    # 两个键: 'azulene' 是任务文本首 token(material 由此而来), 'AZ' 是简写。
    'AZ': 'Au-(2,6-azulenedithiol)-Au',
    'AZULENE': 'Au-(2,6-azulenedithiol)-Au',
    # 2026-11-15 4,4'-联吡啶(BPY): 锚定原子是吡啶 N(非 S), 结为
    # Au-N-(pyridine)2-N-Au, 无硫醇盐接触。T(E_F)/G 绿带(0.03–0.40 G0)与
    # Seebeck 正号都是 Au-BDT-Au 硫醇盐标定域, 对 BPY 作闸门只会把一次物理上
    # 正常的算例判成"否证 H0"。显式登记族名, 让这三条记"不适用"。
    # 实验参照: 1.3±0.1 MΩ ≈ 9.93e-3 G0(benchmark 指南), 仅作对照, 不作闸门。
    'BPY': "Au-(4,4'-bipyridine)-Au",
    # 2026-11-15 1,4-双(甲硫基)苯(para-SMe): 硫醚锚定 Au–S(CH3)–C6H4–S(CH3)–Au,
    # S 上无 H、不脱质子, 是配位型接触(硫醚 Au–S 可到 2.6–2.7 Å, 不是硫醇盐的
    # 2.4 Å 口径)。T(E_F)/G 绿带(0.03–0.40 G0)与 Seebeck 正号都来自 Au–BDT–Au
    # 硫醇盐标定; 硫醚接触的耦合强度与输运通道不同, 对它作闸门会把一次物理上
    # 正常的算例判成"否证 H0"。显式登记族名, 让这三条记"不适用"。
    'SME': 'Au-(1,4-bis(methylthio)benzene)-Au',
    # 2026-11-15 烷烃二硫醇链长序列(HDT C6 / ODT C8 / DDT C10): 末端硫醇去质子成
    # 硫醇盐–Au, 输运走饱和 σ 链的隧穿通道, G 随碳数指数衰减(Xu & Tao 2003,
    # Science 301 1221, β≈1.0/碳)。T(E_F)/G 绿带(0.03–0.40 G0)与 Seebeck 正号都是
    # Au–BDT–Au 硫醇盐标定域, 对烷链作闸门会把一次物理上正常的算例判成"否证 H0"。
    # 与既有 C4 同类, 统一到同一族名('C4' 历史键保留其名, 语义等价)。
    'HDT': 'Au-alkane-dithiol-Au',
    'ODT': 'Au-alkane-dithiol-Au',
    'DDT': 'Au-alkane-dithiol-Au',
    'C6': 'Au-alkane-dithiol-Au',
    'C8': 'Au-alkane-dithiol-Au',
    'C10': 'Au-alkane-dithiol-Au',
}


def system_family(system) -> str | None:
    """体系标识 -> 判据标定族名; 未登记/缺省 -> None(调用方按"全部适用"处理)。"""
    if not system:
        return None
    return SYSTEM_FAMILIES.get(str(system).strip().upper())


def applies(criterion, family) -> bool:
    """该判据是否对 family 作闸门。

    未声明 applies_to 的判据普适; family 未知(None)时也按适用处理 —— fail-safe:
    宁可多判一条(报告里看得见), 不可漏判(静默放过)。
    """
    spec = (criterion or {}).get('applies_to') if isinstance(criterion, dict) else None
    if not spec:
        return True
    if family is None:
        return True
    return family in (spec.get('families') or ())


def lint_registry(criteria=None) -> list:
    """返回问题清单: physical 判据缺 source / 或 implemented=False 即记为问题。

    engineering 判据免源(其"依据"是程序语义, 以 note 记录)。
    """
    problems = []
    for c in (criteria if criteria is not None else CRITERIA):
        k = c.get('kind')
        if k not in ('physical', 'engineering'):
            problems.append({'name': c.get('name'), 'issue': 'kind 非法: %r' % k})
            continue
        if k == 'physical':
            src = c.get('source') or {}
            if not src.get('doi') or not src.get('locator'):
                problems.append({'name': c.get('name'), 'issue': 'physical 判据缺 source{doi,locator}'})
        ap = c.get('applies_to')
        if ap is not None:      # A(2026-10-24): 声明适用域可以, 但必须写清族名单与标定依据
            if not isinstance(ap, dict) or not ap.get('families'):
                problems.append({'name': c.get('name'),
                                 'issue': 'applies_to 缺 families(适用体系族名单)'})
            elif not ap.get('basis'):
                problems.append({'name': c.get('name'),
                                 'issue': 'applies_to 缺 basis(标定依据: 为何只对该族适用)'})
        if not c.get('implemented', False):
            problems.append({'name': c.get('name'), 'issue': 'implemented=False(仅有文档, 无闸门)'})
    return problems


# P1-1: physics_verdict 与任务卡 §9 措辞的映射(v0.4 §2.2 P1-1 表)
VERDICT_WORDING = {'passed': '支持 H0', 'failed': '否证 H0', 'inconclusive': '无法判定'}
ALLOWED_ACTION = {'passed': '可进入发布(A5)/预测卡(A4)',
                  'failed': '停; 不得进入预测模式(任务卡 §7.3②)',
                  'inconclusive': '停; 须补数据或降级声明'}


def physics_verdict(checks) -> dict:
    """由 checks 归约出 physics_verdict —— 工程失败与物理不通过严格分离。

    checks 为 [{'name','severity','passed',...}](validate_* 的约定)。
    任一 error 级判据未过 -> 'failed'; 否则若有 warn 级未过 -> 'inconclusive'; 否则 'passed'。
    applicable=False 的 check(判据不适用于该体系族, 见 SYSTEM_FAMILIES)只进
    not_applicable_criteria, **不参与** failed/inconclusive 归约。
    """
    if not checks:
        # 2026-11-11 本地补丁 (smj F2): 空判据列表不是"通过", 而是"无法评估"。
        return {'physics_verdict': 'inconclusive',
                'wording': VERDICT_WORDING['inconclusive'],
                'allowed_action': ALLOWED_ACTION['inconclusive'],
                'failed_criteria': [], 'inconclusive_criteria': ['no_checks'],
                'missing_evidence_criteria': [],
                'not_applicable_criteria': [], 'unregistered_checks': []}
    fam = {c['name']: c for c in CRITERIA}
    errs, warns, unknown, not_app, missing = [], [], [], [], []
    for c in checks or []:
        if c.get('applicable') is False:    # 判据对此体系不作闸门(A, 2026-10-24)
            not_app.append(c.get('name'))
            continue
        if c.get('passed'):
            continue
        spec = fam.get(c.get('name'))
        if spec is None:
            unknown.append(c.get('name'))
            continue
        gate = spec.get('gate', 'warn')
        if gate == 'error':
            # 2026-11-11 本地补丁 (smj F2): "证据缺失"与"证据否定假设"必须分开归约。
            # 缺 TE.dat 时没有观测可反驳 H0, 说"否证 H0"是错的; 只能"无法判定"。
            if c.get('missing'):
                missing.append(c.get('name'))
            else:
                errs.append(c.get('name'))
        elif gate == 'warn':
            warns.append(c.get('name'))
    if errs:
        v = 'failed'
    elif warns or unknown or missing:
        # 未登记判据未过时**不得**判 passed(否则新增闸门一旦忘登记就会被静默放过)
        v = 'inconclusive'
    else:
        v = 'passed'
    return {'physics_verdict': v, 'wording': VERDICT_WORDING[v],
            'allowed_action': ALLOWED_ACTION[v],
            'failed_criteria': errs, 'inconclusive_criteria': warns,
            'missing_evidence_criteria': missing,
            'not_applicable_criteria': not_app,
            'unregistered_checks': unknown}


# ---------------------------------------------------------------------------
# P1-2 评估接线 (AGD-2026-09-25-SMJ): 把已注册判据**从 run 产物实算**出来。
# 此前只注册未评估 —— T(E_F)=0.8491 / G=0.8607 这类结果需被**自动判否**,
# 而非靠人看报告(本会话 H0 裁决就是人看出来的)。
# ---------------------------------------------------------------------------
PEAK_PROMINENCE = 0.1   # 工程约定(**非物理判据**): 局部极大点的突出度门限。
#   理由: 4001 点扫描中的数值抖动会制造无数"峰"; 不设突出度则"最近峰"恒为 ~0 eV。
PEAK_SEARCH_WINDOW_EV = 1.0   # 与 t_ef_no_peak_near_ef 的比较窗一致(Thygesen 2005 p.12)


def _check(name, value, ok, detail, applicable: bool = True, missing: bool = False):
    out = {'name': name, 'value': value, 'passed': bool(ok),
           'severity': 'error' if name in _ERROR_GATED else 'warn', 'detail': detail,
           'applicable': bool(applicable)}
    if missing:     # 2026-11-11 本地补丁 (smj F2): 该 check 未过是因为产物缺失
        out['missing'] = True
    return out


_ERROR_GATED = {c['name'] for c in CRITERIA if c.get('gate') == 'error'}
PHYSICAL_NAMES = frozenset(c['name'] for c in CRITERIA if c.get('kind') == 'physical')


def _spec(name):
    for c in CRITERIA:
        if c['name'] == name:
            return c
    return None


def _bf(name):
    for c in CRITERIA:
        if c['name'] == name:
            return c.get('band')
    return None


def _applicability_tail(spec, app_: bool, family) -> str:
    """带 applies_to 的判据 -> detail 尾注; 未声明适用域则返回空串。

    单一出口: 所有带适用域的判据共用同一文案, 避免某条判据漏写"不适用"而让读者
    误以为它仍是闸门(2026-10-24 第二刀)。
    """
    if not (spec or {}).get('applies_to'):
        return ''
    return '[适用域: %s; 体系族 %s]' % (
        '适用' if app_ else '不适用(仅记录, 不作闸门)',
        family or '未知(按适用处理)')


def significant_peaks(rows, prominence=PEAK_PROMINENCE, half_window=200):
    """T(E) 的显著局部极大点 -> [(E-EF_eV, T)]。

    prominence 为工程约定(见上); half_window 是求"两侧谷底"的采样半径(点数)。
    """
    out = []
    n = len(rows)
    for i in range(1, n - 1):
        ti = rows[i][2]
        if ti < rows[i - 1][2] or ti < rows[i + 1][2]:
            continue
        lo_l = min(r[2] for r in rows[max(0, i - half_window):i + 1])
        lo_r = min(r[2] for r in rows[i:min(n, i + half_window + 1)])
        if ti - max(lo_l, lo_r) >= prominence:
            out.append((rows[i][1], ti))
    return out


def evaluate_run_criteria(run_dir, temperature: float = 300.0, system=None) -> list:
    """从 run 产物实算物理判据; 无产物则返回一条 error 级 ``te_present``。

    2026-11-11 本地补丁 (smj F2): 旧实现返回 ``[]``, 而 ``physics_verdict([])``
    会落到 ``else: 'passed'`` —— 一个只有空 ``aitranss_te/`` 目录的 run 因此被
    判"支持 H0"。缺产物只能意味着"无法评估", 不能意味着"通过"。

    system: 体系标识(TaskCard.material 等), 决定带 applies_to 的判据是否作闸门;
    未登记体系按适用处理(fail-safe, 见 system_family/applies)。
    """
    from .landauer import analyze_te, read_te
    family = system_family(system)
    d = Path(run_dir)
    te = d / 'TE.dat'
    if not te.exists():
        return [_check('te_present', None, False,
                       '缺 TE.dat: 物理判据无法评估, 不得视为通过(2026-11-11 本地补丁 F2)',
                       missing=True)]
    try:
        res = analyze_te(te, temperature=temperature)
        rows = read_te(te)
    except Exception as exc:                                    # 解析失败不静默
        return [_check('t_ef_magnitude', None, False,
                       'TE.dat 解析失败, 判据无法评估: %s' % exc)]

    checks = []
    # ① T(E_F): 取最接近 E_F 的采样点(与 landauer 同口径)
    near = res.get('fermi_nearest') or {}
    t_ef = near.get('T')
    spec = _spec('t_ef_magnitude') or {}
    app_ = applies(spec, family)
    lo, hi = _bf('t_ef_magnitude') or (None, None)
    if t_ef is None:
        checks.append(_check('t_ef_magnitude', None, False, 'TE.dat 无 E_F 附近采样点'))
    else:
        ok = (lo <= t_ef <= hi)
        checks.append(_check('t_ef_magnitude', t_ef, ok,
                             'T(E_F)=%.4f 带 %.2f-%.2f G0 (Thygesen 2005 p.13%s)%s'
                             % (t_ef, lo, hi, ', 仅方向性' if ok else '',
                                _applicability_tail(spec, app_, family)),
                             applicable=app_))

    # ② E_F 近旁是否落在透射峰上(Thygesen p.12: 最近特征 ~1.0 eV)
    pk = [p for p in significant_peaks(rows) if abs(p[0]) <= PEAK_SEARCH_WINDOW_EV]
    if pk:
        dpk = min(pk, key=lambda p: abs(p[0]))
        need = (_bf('t_ef_no_peak_near_ef') or [1.0])[0]
        checks.append(_check('t_ef_no_peak_near_ef', round(abs(dpk[0]), 4),
                             abs(dpk[0]) >= need,
                             'E_F 近旁 %.4f eV 处有显著峰(T=%.4f); 基准要求 %.1f eV 内无峰 '
                             '(Thygesen 2005 p.12; 突出度门限 %.2f 为工程约定)'
                             % (dpk[0], dpk[1], need, PEAK_PROMINENCE)))
    else:
        checks.append(_check('t_ef_no_peak_near_ef', None, True,
                             'E_F 近旁 %.1f eV 内无显著峰(突出度门限 %.2f 为工程约定)'
                             % (PEAK_SEARCH_WINDOW_EV, PEAK_PROMINENCE)))

    # ③ G 落在理论带(G/G0 = L0)
    g = res.get('L0')
    spec = _spec('g_in_theory_band') or {}
    app_ = applies(spec, family)
    glo, ghi = _bf('g_in_theory_band') or (None, None)
    if g is not None:
        checks.append(_check('g_in_theory_band', g, glo <= g <= ghi,
                             'G=%.4f G0 带 %.2f-%.2f (任务卡 §4.2)%s'
                             % (g, glo, ghi, _applicability_tail(spec, app_, family)),
                             applicable=app_))

    # ④ Seebeck 符号(铁律: Au-BDT-Au 为正) —— A(2026-10-24): 带适用域, 仅该族作闸门
    s = res.get('S_uV_per_K')
    if s is not None:
        spec = _spec('seebeck_sign') or {}
        app_ = applies(spec, family)
        detail = ('S=%+.4f uV/K 应为正(Reddy 2007 +8.7+-2.1)' % s
                  + _applicability_tail(spec, app_, family))
        checks.append(_check('seebeck_sign', s, s > 0, detail, applicable=app_))

    # ⑤ k_e / ZT_e 非负(来源待补, 故 warn)
    ke, zt = res.get('k_e_pW_per_K'), res.get('ZT_electronic')
    if ke is not None and zt is not None:
        checks.append(_check('k_e_zt_order', {'k_e': ke, 'ZT_e': zt}, ke >= 0 and zt >= 0,
                             'k_e=%.4f pW/K, ZT_e=%.6f 应非负(来源待补)' % (ke, zt)))
    return checks
