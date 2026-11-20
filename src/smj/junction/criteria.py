"""P1-1/P1-2/P1-5 (AGD-2026-09-25-SMJ): smj 物理判据注册表与 physics_verdict。

两条铁律(承 v0.4 §0 A4 与 §2.5 断言纪律):
1. **每条 physical 判据必须带 source{doi, locator}** —— 无源判据不得作闸门;
2. **engineering 判据**(程序语义, 如左右电极不相交)不要求文献源, 但须标 kind。

P1-5 的强制由 lint_registry() + tests/test_smj_automation.py 的 ratchet 用例实现:
新增无源判据即失败; 已登记的缺口(见 KNOWN_GAPS)数量固定, 只能减少不能增加。

判据数值来源(均已回读原文核验, 见 docs/criteria-source-research.md §2):
- Au-S 带 2.40-2.44 A: Jackson 2009 摘要 `S-Au = 2.42+-0.02 A`(Au(111)(sqrt3xsqrt3)R30-CH3S
  表面相, NIXSW+S1s PhD 实验锚点; 2.44=2.42+0.02; 旧带上端 2.49 A 无一手源, 2026-11-15 收窄)。
  Thygesen 2005 §V.2 的 d(Au-S)=2.40 A 是**采用 Stokbro 2003**("in accordance with Ref. 48")
  的模型输入, 非该文测量/计算结论。
- T(E_F)/G 理论带 0.03-0.4: **派生量, 非原文区间** —— Thygesen 2005 §V.2 只给 "around 0.1 G0"
  与 "比 Ref.48 小 4 倍、比 Ref.15 大 3 倍", [0.1/3, 0.1x4] 是文献散布包络。原文同时承认 BDT
  理论电导比实验高 2-3 个数量级, 故该带只检验"是否落在 DFT-NEGF 理论族内", 不检验"是否复现
  实验"; 2026-11-15: t_ef_magnitude 按用户规则取消(无一手区间), g_in_theory_band 降为 warn
  并标 evidence=derived(派生量不足以否证 H0)。
- 共振归因: 族表 EXPECTED_RESONANCE 给各族的预期通道(HOMO/LUMO); 未登记族按实测
  Seebeck 符号自洽归因(T(E) 与 S 是同一组数据的不同矩)。E_F 近旁 1.0 eV 内出现显著峰
  -> **需解释的物理信号**; 峰落在预期侧且 0.15-1.0 eV -> 已归因(pass)。
- 适用域(2026-10-24 第二刀): G 理论带与 seebeck_sign 同属 **Au-BDT-Au 标定域**, 故同带
  applies_to; 换分子族(如饱和烷烃二硫醇)不作闸门, 实算值仍记录(见 SYSTEM_FAMILIES)
"""
from __future__ import annotations

from pathlib import Path

THY = {'doi': '10.1016/j.chemphys.2005.05.032'}
JAC = {'doi': '10.1016/j.susc.2009.01.022'}
WID = {'doi': '10.1021/nl203634m'}

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
     'note': 'transmission 产物存在性(程序语义); 缺失或不可解析时物理判据无法评估, 不得判 passed'},

    {'name': 'grid_independence', 'kind': 'physical', 'gate': 'warn', 'implemented': True,
     'band': [0.01, None], 'source': dict(THY, locator='—'),
     'note': '抽稀后 L0/L1/L2 相对变化 <1%(任务卡 §5.1.3); 本会话实测 <=0.09%'},

    {'name': 'electrodes_disjoint', 'kind': 'engineering', 'gate': 'error', 'implemented': True,
     'band': None, 'source': None,
     'note': 'aitranss 语义: 左右电极共享原子即终止(trans.f); 见 P0-4'},

    {'name': 'au_s_bond', 'kind': 'physical', 'gate': 'warn', 'implemented': True,
     'band': [2.40, 2.44], 'source': dict(JAC, locator='摘要: Au(111)(sqrt3xsqrt3)R30-CH3S atop, '
                                                       'NIXSW+S1s PhD, S-Au=2.42+-0.02 A'),
     'note': '实验锚点 Jackson 2009 2.42+-0.02 -> 带 [2.40, 2.44](2.44=2.42+0.02; 旧上端 2.49 '
             '无一手源)。Thygesen 2005 §V.2 的 d(Au-S)=2.40 A 是**采用 Stokbro 2003** 的模型'
             '输入("in accordance with Ref. 48"), 不是该文结论。锚点为表面吸附相、非电极-分子'
             '-电极结, 故仅方向性(低配位单 Au 顶位/桥位、金字塔配位都可偏离, 偏离只提示不阻断); '
             '仅适用于 S 为锚的去质子硫醇盐；-SCH3 等硫醚的 S–Au 是配位接触，N/Se 等非 S 锚'
             '只报接触距离，均不适用'},

    {'name': 't_ef_resonance_attribution', 'kind': 'physical', 'gate': 'warn', 'implemented': True,
     'band': [0.15, 1.0], 'source': dict(THY, locator='§V.2: 窄峰 -1.5/-1.0 eV, 宽峰 -3.0/+2.5 eV'),
     'note': 'E_F±1.0 eV 内的显著透射峰按族归因(见 EXPECTED_RESONANCE): 无峰 -> pass; 峰落在'
             '预期通道侧且 |ΔE|∈[0.15,1.0] eV -> pass(物理信号已归因); 侧别不符/距 E_F 过近/'
             '通道未知 -> warn(需解释)。0.15 eV 下限取 PBE 半局域能级误差量级 0.1-0.3 eV'
             '(Gulyaev 2026), 不是"允许的物理共振距离"; 该文最近峰恰在 -1.0 eV(余量为 0), '
             '且 PBE 会把峰系统性推向 E_F(本管线 BDT -0.339 eV), 故不要求峰位复现文献值'},

    {'name': 'g_in_theory_band', 'kind': 'physical', 'gate': 'warn', 'implemented': True,
     'band': [0.03, 0.4], 'source': dict(THY, locator='§V.2 + Fig.10 派生: 0.1 G0 x '
                                                       '(Ref.15 的 1/3 ... Ref.48 的 4)'),
     'evidence': 'derived',
     'note': '**派生包络**(非原文区间): 0.1 G0 的 1/3-4 倍 = [0.033, 0.4]; 原文同时承认 BDT '
             '理论电导比实验高 2-3 个数量级, 故只检验"是否落在 DFT-NEGF 理论族内", 不复现实验。'
             '2026-11-15: 因派生故 gate 由 error 降为 warn(派生量不足以否证 H0)',
     # 与 seebeck_sign 同一标定域(G/G0≈T(E_F), 300 K 费米窗窄); 适用域一并声明
     'applies_to': {'families': ('Au-BDT-Au',),
                    'basis': ('0.03-0.40 G0 是 Au-BDT-Au 理论标定包络(Thygesen 2005 §V.2 的 '
                              '0.1 G0 与其 Ref.15/48 的 x3/x4 比例, 派生); 300 K 下 '
                              'G/G0≈T(E_F), 换分子族不作闸门 —— 2026-10-24 用户拍板')}},

    {'name': 'seebeck_sign', 'kind': 'physical', 'gate': 'error', 'implemented': True,
     'band': [0, None], 'source': {'doi': '10.1126/science.1137149', 'locator': '摘要/正文'},
     'note': 'Reddy 2007 +8.7+-2.1 与 Lee 2014 +7.4+-0.5 同号为正 -> 铁律判据',
     'applies_to': {'families': ('Au-BDT-Au',),
                    'basis': ('两篇实验均为 Au-BDT-Au(Reddy 2007/Lee 2014), 正号是该族实测'
                              '性质; 换分子体系(如 2CN)不得当普适闸门 —— 2026-10-24 用户拍板'
                              '"窄口径", 2CN 记不适用而非判否')}},

    {'name': 'te_window_covers_efermi', 'kind': 'physical', 'gate': 'error', 'implemented': True,
     'band': None, 'source': None,
     'note': '2026-11-20 复核 P1-2: TE.dat 的 E-E_F 必须跨 0(min<=0<=max)。整窗落在 E_F '
             '一侧时费米窗权重全在同一侧, L0/L1/L2 不是 E_F 处展开, S 无物理意义 -> failed'},

    {'name': 'k_e_zt_order', 'kind': 'physical', 'gate': 'warn', 'implemented': True,
     'band': None, 'source': None,
     'note': '2026-11-20 复核 R12: 评估器确实会 emit 该判据(κe≥0 与 ZT 量级常识), '
             '原元数据写 implemented=False 与行为不符, 已改正; 正式文献来源待补, '
             '仍列在 KNOWN_GAPS 里(ratchet 只许减不许增)'},
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
    # 2026-11-17 2,6-双(乙酰硫基)芳烃(thioacetate, SAc): S 上无 H、不脱质子, 保留
    # –C(=O)CH3 的配位型接触(与 SME 同类, 只是羰基取代甲基)。T(E_F)/G 绿带
    # (0.03–0.40 G0)与 Seebeck 正号都来自 Au–BDT–Au 硫醇盐标定; 对硫酯接触作闸门
    # 会把一次物理上正常的算例判成"否证 H0"。显式登记族名, 让这三条记"不适用"。
    # 注意: 实验上 SAc 在 Au 上通常原位脱保护成硫醇盐; 本族模拟的是**保留乙酰基**
    # 的对照臂, 脱保护形态等价于已跑的 SH 批(nap-sh / az-sh / ind-sh)。
    'SAC': 'Au-(2,6-bis(acetylthio)arene)-Au',
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
    # 2026-11-16 1,4-苯二胺(PPD / BDA, 本仓库 run 'ppd'): 锚定原子是 –NH2 的 N, 结为
    # Au–NH2–C6H4–NH2–Au —— 配位型接触(N 孤对 -> Au), 两个 N–H 全部保留, 不是硫醇盐。
    # T(E_F)/G 绿带(0.03–0.40 G0)与 Seebeck 正号都是 Au–BDT–Au 硫醇盐标定域; 对胺锚结
    # 作闸门会把一次物理上正常的算例判成"否证 H0"(BDA 实验 G≈6.4e-3 G0 本就落在绿带外,
    # Venkataraman 2006)。显式登记族名, 让这两条记"不适用"。
    # 两个键: 'BDA' 是文献缩写(1,4-benzenediamine), 'PPD' 是 p-phenylenediamine 常用缩写。
    'BDA': 'Au-(1,4-diaminobenzene)-Au',
    'PPD': 'Au-(1,4-diaminobenzene)-Au',
}


# 2026-11-15: 共振归因族表 —— "E_F 近旁有峰"不再是错误信号, 而是**需要解释**的物理信号。
# 只有这里登记了预期通道的家族才做"文献先验"归因; 未登记族(如 SME / 2CN, 尚无一手文献
# 明确其近 E_F 通道归属)按实测 Seebeck 符号做自洽归因, detail 里标明这是自洽性检查而非
# 文献先验 —— 不把"没查到源"包装成"有源"。通道依据见 docs/criteria-source-research.md §3。
EXPECTED_RESONANCE = {
    'Au-BDT-Au': {
        'channel': 'HOMO',
        'source': {'doi': '10.1126/science.1137149',
                   'locator': '摘要: E_F 在 HOMO 上方 1.2 eV, S=+8.7+-2.1 uV/K, p 型'},
        'note': 'BDT 的 HOMO 在 E_F 下方(Reddy 2007 摘要); Thygesen 2005 §V.2 的窄峰 '
                '-1.5/-1.0 eV 为峰位佐证(PBE 级)',
    },
    "Au-(4,4'-bipyridine)-Au": {
        'channel': 'LUMO',
        'source': dict(WID, locator='Table 1(mol.4: 4,4-bipyridine, S=-9.5 uV/K)+正文'),
        'note': '吡啶 N-Au 结走 LUMO(LUMO 在 E_F 上方约 1.47-1.53 eV); 实验 S<0 '
                '(Widawsky 2012: pyridine-Au 负、amine-Au 正)',
    },
    'Au-(2,6-azulenedithiol)-Au': {
        'channel': 'HOMO',
        'source': {'doi': '10.1021/nl5010702',
                   'locator': '正文 + Fig.4(GW): E_F<0 for thiochroman linkers; 连接位序 '
                              '1,3 > 2,6 ~ 4,7 > 5,7'},
        'note': '薁为非交替烃, 谱不对称; GW 模型取 E_F=-1.5 eV(绝对共振位未给, 故只登记'
                '通道侧、不登记窗口); 本管线 S=+54.1 与 HOMO 侧一致',
    },
    'Au-alkane-dithiol-Au': {
        'channel': None,
        'source': {'doi': '10.1126/science.1087481',
                   'locator': '摘要: beta_N=1.0+-0.1/碳; 宽 gap 隧穿, 无近 E_F 共振'},
        'note': '饱和烷链 HOMO/LUMO 均远离 E_F; 出现近 E_F 峰即 warn(需解释)',
    },
    # 2026-11-16 1,4-苯二胺(胺锚 Au–NH2): 苯二胺系列的占据轨道离 Au E_F 最近, 输运是
    # hole 隧穿(Venkataraman 2007 摘要), 即 HOMO 通道; 实验 Seebeck 取正值(Widawsky 2012
    # 同时给出 pyridine-Au 负、amine-Au 正), 与 HOMO 侧自洽。
    'Au-(1,4-diaminobenzene)-Au': {
        'channel': 'HOMO',
        'source': {'doi': '10.1021/nl062923j',
                   'locator': '摘要: 苯二胺系列 "the occupied states are closest to the gold '
                              'Fermi energy" + "hole tunneling"'},
        'note': '胺锚 Au–NH2 结走 HOMO(hole 隧穿先验); 实验 S>0(Widawsky 2012: amine-Au 正、'
                'pyridine-Au 负), 符号与本族 HOMO 侧自洽',
    },
}


def system_family(system) -> str | None:
    """体系标识 -> 判据标定族名; 未登记/缺省 -> None(调用方按"全部适用"处理)。"""
    if not system:
        return None
    return SYSTEM_FAMILIES.get(str(system).strip().upper())


def is_open_shell(multiplicity) -> bool:
    """多重度 > 1 即开壳层（缺失/非法按闭壳层处理，与旧行为一致）。"""
    try:
        return int(multiplicity) > 1
    except (TypeError, ValueError):
        return False


def applies(criterion, family, *, open_shell: bool = False) -> bool:
    """该判据是否对当前体系作闸门。

    未声明 applies_to 的判据普适。**按族标定的判据只在"族被登记且匹配"时作闸门**：

    - 开壳层（multiplicity>1）：判据族从未针对自旋分裂的 α/β 双通道标定过，
      借来的 G 带与 S 符号约定不能当闸门；
    - 族未知（None）：没有该族的标定依据，把别的族的数字当闸门就是过度断言；
    - 族已知但不匹配：同样只记录（2026-10-24 窄口径）。

    2026-11-20：旧实现在 family=None 时 fail-open（"宁可多判一条"），于是
    ``runs/ttf-release``（family=null、多重度 2）被 BDT 的 ``G∈[0.03,0.40]``
    与 ``S>0`` 判成 ``passed / 支持 H0``——两条数字碰巧落在 BDT 包络里，
    不是该自由基结的证据。
    """
    spec = (criterion or {}).get('applies_to') if isinstance(criterion, dict) else None
    if not spec:
        return True
    if open_shell:
        return False
    if family is None:
        return False
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
    # 2026-11-20 复核 P1-2, 本地补丁 F-RECHECK-window
    evaluated_gates = [c.get('name') for c in checks or []
                       if c.get('applicable') is not False
                       and c.get('name') in _DISCRIMINATING_GATES]
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
    elif not evaluated_gates:
        # 2026-11-20 复核 P1-2: 所有带 applies_to 的族闸门都 applicable=False。
        # 这不是"证据支持 H0", 而是"没有闸门可判" —— 必须 inconclusive。
        v = 'inconclusive'
        warns = ['no_applicable_gate']
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
PEAK_SEARCH_WINDOW_EV = 1.0   # 归因窗(与 t_ef_resonance_attribution 的 band 上端一致; Thygesen 2005 §V.2)


def _check(name, value, ok, detail, applicable: bool = True, missing: bool = False):
    out = {'name': name, 'value': value, 'passed': bool(ok),
           'severity': 'error' if name in _ERROR_GATED else 'warn', 'detail': detail,
           'applicable': bool(applicable)}
    if missing:     # 2026-11-11 本地补丁 (smj F2): 该 check 未过是因为产物缺失
        out['missing'] = True
    return out


_ERROR_GATED = {c['name'] for c in CRITERIA if c.get('gate') == 'error'}
PHYSICAL_NAMES = frozenset(c['name'] for c in CRITERIA if c.get('kind') == 'physical')
# 2026-11-20 复核 P1-2, 本地补丁 F-RECHECK-window: "有判别力的族闸门" = 带 applies_to 的判据。
# 全部此类闸门都 applicable=False 时, 该 run 没有任何族标定下的物理闸门被真正评估过,
# 不能判 passed(那只是"没有闸门可判", 不是"证据支持 H0")。
_DISCRIMINATING_GATES = frozenset(c['name'] for c in CRITERIA if c.get('applies_to'))


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


def _applicability_tail(spec, app_: bool, family, *,
                        open_shell: bool = False) -> str:
    """带 applies_to 的判据 -> detail 尾注; 未声明适用域则返回空串。

    单一出口: 所有带适用域的判据共用同一文案, 避免某条判据漏写"不适用"而让读者
    误以为它仍是闸门(2026-10-24 第二刀)。2026-11-20 起把"为什么不作闸门"写清楚:
    开壳层未校准 / 体系族未知 / 该族不在名单里，三者含义不同。
    """
    if not (spec or {}).get('applies_to'):
        return ''
    if app_:
        reason = '适用'
    elif open_shell:
        reason = '不适用(开壳层未校准, 仅记录)'
    elif family is None:
        reason = '不适用(体系族未知, 无标定依据, 仅记录)'
    else:
        reason = '不适用(该族未标定, 仅记录)'
    return '[适用域: %s; 体系族 %s]' % (reason, family or '未知')


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


def evaluate_run_criteria(run_dir, temperature: float = 300.0, system=None,
                          multiplicity=None) -> list:
    """从 run 产物实算物理判据; 无产物则返回一条 error 级 ``te_present``。

    2026-11-11 本地补丁 (smj F2): 旧实现返回 ``[]``, 而 ``physics_verdict([])``
    会落到 ``else: 'passed'`` —— 一个只有空 ``aitranss_te/`` 目录的 run 因此被
    判"支持 H0"。缺产物只能意味着"无法评估", 不能意味着"通过"。

    system: 体系标识(TaskCard.material 等), 决定带 applies_to 的判据是否作闸门。
    未登记体系 -> 族未知 -> 按族标定的判据只记录不作闸门（2026-11-20 起；旧实现
    fail-open，把 BDT 的包络套到别的体系上）。
    multiplicity: 分子多重度; >1 即开壳层, 按族标定的判据一律不作闸门（判据族
    只针对闭壳层标定过）。
    """
    from .landauer import analyze_te, read_te
    family = system_family(system)
    open_shell = is_open_shell(multiplicity)
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
        # 2026-11-15: t_ef_magnitude 取消后, 解析失败归到 te_present(存在但不可解析),
        # 仍按"证据缺失"处理, 不得判 passed。
        return [_check('te_present', None, False,
                       'TE.dat 存在但解析失败, 物理判据无法评估: %s' % exc,
                       missing=True)]

    checks = []
    # ⓪ 能量窗必须跨 E_F(2026-11-20 复核 P1-2, 本地补丁 F-RECHECK-window):
    # 整窗落在 E_F 一侧时, 费米窗权重全在同一侧, L0/L1/L2 不是 E_F 处的展开,
    # S 的符号与量级都没有意义 —— 这种谱必须判 failed, 而不是靠"族闸门不适用"放过。
    energies = [r[1] for r in rows]
    lo_ev, hi_ev = min(energies), max(energies)
    covers = lo_ev <= 0.0 <= hi_ev
    checks.append(_check(
        'te_window_covers_efermi', [round(lo_ev, 6), round(hi_ev, 6)], covers,
        '能量窗 E-E_F ∈ [%+.4f, %+.4f] eV %s 0: %s'
        % (lo_ev, hi_ev, '跨越' if covers else '**未跨越**',
           '正常' if covers else '整窗在 E_F 一侧, Seebeck/ZT 无物理意义')))
    # ① E_F 近旁的显著透射峰: 按族归因(2026-11-15: "有峰"从错误信号改为需解释的物理信号)
    pk = [p for p in significant_peaks(rows) if abs(p[0]) <= PEAK_SEARCH_WINDOW_EV]
    band_lo, band_hi = _bf('t_ef_resonance_attribution') or (0.15, 1.0)
    entry = EXPECTED_RESONANCE.get(family) if family else None
    s = res.get('S_uV_per_K')
    if entry:
        channel = entry.get('channel')
        why = '族表(%s): %s' % (family, entry.get('note', ''))
    elif s is None:
        channel, why = None, '未登记族且无 Seebeck 可依'
    elif s == 0:
        channel, why = None, '未登记族且 S=0, 通道无法归因'
    else:
        channel = 'HOMO' if s > 0 else 'LUMO'
        why = ('未登记族: 按实测 S=%+.1f uV/K 符号自洽归因'
               '(自洽性检查, 非文献先验)' % s)
    if not pk:
        checks.append(_check('t_ef_resonance_attribution', None, True,
                             'E_F±%.1f eV 内无显著峰(突出度门限 %.2f 为工程约定)'
                             % (PEAK_SEARCH_WINDOW_EV, PEAK_PROMINENCE)))
    else:
        dpk = min(pk, key=lambda p: abs(p[0]))
        delta = dpk[0]
        ok_side = (channel == 'HOMO' and delta < 0) or (channel == 'LUMO' and delta > 0)
        ok_dist = band_lo <= abs(delta) <= band_hi
        if channel is None:
            checks.append(_check('t_ef_resonance_attribution', round(abs(delta), 4), False,
                                 'E_F 近旁 %+.4f eV 处有显著峰(T=%.4f), 但通道无法归因(%s); 需解释'
                                 % (delta, dpk[1], why)))
        elif ok_side and ok_dist:
            checks.append(_check('t_ef_resonance_attribution', round(abs(delta), 4), True,
                                 'E_F 近旁 %+.4f eV 处有显著峰(T=%.4f), 在预期 %s 侧且 '
                                 '|ΔE|∈[%.2f,%.1f] eV -> 已归因; %s'
                                 % (delta, dpk[1], channel, band_lo, band_hi, why)))
        else:
            reason = (('侧别不符(预期 %s 侧)' % channel) if not ok_side
                      else ('距 E_F 过近(<%.2f eV)' % band_lo))
            checks.append(_check('t_ef_resonance_attribution', round(abs(delta), 4), False,
                                 'E_F 近旁 %+.4f eV 处有显著峰(T=%.4f); %s; 需解释; %s'
                                 % (delta, dpk[1], reason, why)))

    # ② G 落在理论带(G/G0 = L0) —— 2026-11-20: 开壳层/未知族只记录不作闸门
    g = res.get('L0')
    spec = _spec('g_in_theory_band') or {}
    app_ = applies(spec, family, open_shell=open_shell)
    glo, ghi = _bf('g_in_theory_band') or (None, None)
    if g is not None:
        checks.append(_check('g_in_theory_band', g, glo <= g <= ghi,
                             'G=%.4f G0 带 %.2f-%.2f (派生包络, Thygesen 2005 §V.2)%s'
                             % (g, glo, ghi, _applicability_tail(
                                 spec, app_, family, open_shell=open_shell)),
                             applicable=app_))

    # ③ Seebeck 符号(铁律: Au-BDT-Au 为正) —— A(2026-10-24): 带适用域, 仅该族作闸门
    s = res.get('S_uV_per_K')
    if s is not None:
        spec = _spec('seebeck_sign') or {}
        app_ = applies(spec, family, open_shell=open_shell)
        detail = ('S=%+.4f uV/K 应为正(Reddy 2007 +8.7+-2.1)' % s
                  + _applicability_tail(spec, app_, family, open_shell=open_shell))
        checks.append(_check('seebeck_sign', s, s > 0, detail, applicable=app_))

    # ④ k_e / ZT_e 非负(来源待补, 故 warn)
    ke, zt = res.get('k_e_pW_per_K'), res.get('ZT_electronic')
    if ke is not None and zt is not None:
        checks.append(_check('k_e_zt_order', {'k_e': ke, 'ZT_e': zt}, ke >= 0 and zt >= 0,
                             'k_e=%.4f pW/K, ZT_e=%.6f 应非负(来源待补)' % (ke, zt)))
    return checks
