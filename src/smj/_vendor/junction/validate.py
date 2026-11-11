"""单分子结程序族 run 校验(validate_junction_run, AGD-2026-09-22-SMJ P4)。

tools/validate.py 的 validate_run 按 family='smj' 分派到本模块。校验语义(24 章 §4/§5):
- gaussian_opt: mol.log 含 Normal termination(完成); 虚频 -> error(不得静默接受)
- aims_tip_relax/aims_scf: geometry.in.next_step / aims.out 存在(完成判定)
- aitranss_te: TE.dat 非空且有数据行; efermi 解析 TE.dat 表头优先, 缺失 fail-closed
- 7 点样本免责: 结果仅验证管线, _like 值不得当真实量报告
纯校验, 不做真实计算。
"""
from __future__ import annotations

import json
import math
from pathlib import Path


def _ok(name: str, value, detail: str) -> dict:
    return {'name': name, 'value': value, 'passed': True, 'severity': 'info', 'detail': detail}


def _err(name: str, value, detail: str) -> dict:
    return {'name': name, 'value': value, 'passed': False, 'severity': 'error', 'detail': detail}


def _warn(name: str, value, detail: str) -> dict:
    return {'name': name, 'value': value, 'passed': False, 'severity': 'warn', 'detail': detail}


# 产物归属: 引擎 run 根下产物按阶段分子目录存放(gaussian_opt/aims_scf/...),
# 而 driver / 手工 run 可能把同一批产物平铺在一个目录。以下映射给出每个产物
# 首选的阶段子目录顺序; 定位时先看 run_dir 本身, 再按序找子目录 —— 两种布局
# 都可校验, 且不靠阶段名反推产物(找不到就返回 None, 不臆断)。
_ARTIFACT_STAGE_HOME = {
    'mol.log': ('gaussian_opt',),
    'geometry.in.next_step': ('aims_tip_relax',),
    'aims.out': ('aims_scf', 'aims_tip_relax'),
    'TE.dat': ('aitranss_te',),
    'results.json': ('junction_analyze',),
}


def _find_artifact(run_dir: Path, name: str):
    """在 run_dir 本身或对应阶段子目录里定位产物; 找不到返回 None。"""
    p = Path(run_dir) / name
    if p.exists():
        return p
    for stage in _ARTIFACT_STAGE_HOME.get(name, ()):
        q = Path(run_dir) / stage / name
        if q.exists():
            return q
    return None


def validate_junction_run(run_dir, task_card=None, results=None, system=None) -> dict:
    """smj run 校验: 按 run_dir(扁平阶段目录或含阶段子目录的 run 根)内实际产物判定。

    system: 体系标识(缺省取 task_card.material), 供 criteria 的判据适用域使用;
    未登记体系按"全部判据适用"处理(A, 2026-10-24 窄口径)。
    """
    run_dir = Path(run_dir)
    if system is None and task_card is not None:
        system = getattr(task_card, 'material', None)
    checks = []
    mol_log = _find_artifact(run_dir, 'mol.log')
    if mol_log is not None:
        text = mol_log.read_text(encoding='utf-8', errors='ignore')
        if 'Normal termination' in text:
            checks.append(_ok('gaussian_termination', True, 'Gaussian Normal termination'))
            if 'imaginary frequencies' in text or 'Frequencies --' in text:
                import re
                freqs = []
                for line in re.findall(r'Frequencies --\s+(.+)', text):
                    freqs.extend(float(x) for x in line.split())
                imag = [f for f in freqs if f < 0]
                # 2026-10-24: emit 名必须与判据注册表(criteria.CRITERIA)一致 —— 旧名
                # 'gaussian_imaginary' 未注册, 归约时只进 unregistered_checks, 使虚频证据
                # 永远只能造成 inconclusive, 不进 failed_criteria(等于静默降级)。
                if imag:
                    checks.append(_err('gaussian_no_imaginary', imag,
                                       'Gaussian 虚频 %s cm^-1: 不得静默接受(24 章 §5)' % imag))
                else:
                    checks.append(_ok('gaussian_no_imaginary', 0, '无虚频(结构稳定)'))
        else:
            checks.append(_err('gaussian_termination', False, 'mol.log 无 Normal termination'))
    if _find_artifact(run_dir, 'geometry.in.next_step') is not None:
        checks.append(_ok('aims_tip_relax', True, 'geometry.in.next_step 存在(tip relax 完成)'))
    aims_out = _find_artifact(run_dir, 'aims.out')
    if aims_out is not None:
        aims_text = aims_out.read_text(encoding='utf-8', errors='ignore')
        # 24 章 §5: SCF 未收敛/异常终止不得静默接受(§9.2 失败演练)。
        # 只匹配特异失败标记(FHI-aims 常见 `** ERROR **`/`!ERROR!` 大写), 避免泛 error 误伤。
        if any(m in aims_text for m in ('SCF is not yet converged',
                                        'maximum number of scf iterations',
                                        'geometry optimization has not yet converged',
                                        '** ERROR **', '!ERROR!', '!!!ERROR')):
            checks.append(_err('aims_scf', True,
                                'FHI-aims 输出含失败标记(SCF 未收敛/异常): 不得静默降级, 见 caselib 24-§smj-2'))
        else:
            checks.append(_ok('aims_scf', True, 'aims.out 存在且无失败标记(SCF 完成)'))
    te = _find_artifact(run_dir, 'TE.dat')
    if te is not None:
        rows = [l for l in te.read_text(encoding='utf-8', errors='ignore').splitlines()
                if l.strip() and not l.lstrip().startswith('#')]
        if len(rows) >= 2:
            checks.append(_ok('aitranss_te', len(rows), 'TE.dat 含 %d 数据行' % len(rows)))
            # efermi 解析纪律: TE.dat 表头优先(24 章 §4/§5)
            header = te.read_text(encoding='utf-8', errors='ignore').splitlines()[:5]
            efermi = None
            for line in header:
                if 'efermi' in line.lower():
                    import re
                    m = re.search(r'efermi\s*=\s*([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)', line, re.IGNORECASE)
                    if m:
                        efermi = float(m.group(1))
                        break
            if efermi is not None:
                checks.append(_ok('efermi_source', efermi, 'efermi 取自 TE.dat 表头(aitranss 实际使用值)'))
            else:
                checks.append(_warn('efermi_source', None, 'TE.dat 表头无 efermi; 应回退 aims.out(仍缺失则 fail-closed)'))
        else:
            checks.append(_err('aitranss_te', len(rows), 'TE.dat 数据行不足(<2): 透射未完成'))
    # 7 点样本免责与 _like 语义(结果**若由示例数据产生**)。
    # P0-17(2026-10-24, 真跑缺口 #11 附带): 旧判据"results 里有 junction 就 warn"会把**真实**
    # 输运结果(1501 点 TE.dat 走管线)也当成 7 点示例 —— 实测 C4 真跑被挂上一条假告警,
    # 报告里像"用了示例数据"。改为按数据点数识别: 只有 <=7 点(原工具 examples/TE.dat 的口径)
    # 才是示例样本; 点数缺失(历史/外部)保守仍按示例处理(=旧行为, 不放过)。
    _jres = (results or {}).get('junction') if results else None
    if _jres:
        try:
            _npts = int((_jres or {}).get('points') or 0)
        except (TypeError, ValueError):
            _npts = 0
        if _npts <= 7:
            checks.append(_warn('anchors_disclaimer', True,
                                '示例 7 点样本: 仅验证解析→积分→命名管线, 非生产级科学结论(_like 值不得当真实量报告)'))
    from .landauer import analyze_te  # noqa: F401  (仅导入检查依赖)
    # P1-2: 已注册判据从产物实算并并入(此前只有注册表, 无评估)
    from .criteria import PHYSICAL_NAMES, evaluate_run_criteria, physics_verdict
    # 判据读取 TE.dat; 用实际定位到的目录(阶段子目录布局)而不是 run_dir 本身
    checks.extend(evaluate_run_criteria(te.parent if te is not None else run_dir,
                                        system=system))
    # A(2026-10-24): applicable=False(判据不适用于该体系族)的 check 不作闸门, 不计入
    # errs/warns —— 否则 2CN 的 seebeck_sign 会把工程判定与科学判定一起拉垮。
    _appl = lambda c: c.get('applicable') is not False      # noqa: E731
    errs = sum(1 for c in checks if c['severity'] == 'error' and not c['passed'] and _appl(c))
    warns = sum(1 for c in checks if c['severity'] == 'warn' and not c['passed'] and _appl(c))
    # A3 严格分离: passed 为**工程**判定(不含物理判据); physics_verdict 为科学判定。
    # 若让物理判据也拉垮 passed, 二者就被混成一体, 与本条要求相悖。
    _eng = [c for c in checks if c['name'] not in PHYSICAL_NAMES and _appl(c)]
    engineering_passed = not any(c['severity'] == 'error' and not c['passed'] for c in _eng)
    # P1-1: 另出 physics_verdict(科学判定), 与 passed(工程判定) 严格分离
    return {'passed': errs == 0, 'engineering_passed': engineering_passed,
            'errors': errs, 'warnings': warns, 'checks': checks,
            **physics_verdict(checks)}


# ---------------------------------------------------------------------------
# P0-4 (AGD-2026-09-25-SMJ): 建结几何自洽硬闸门(提交前机检)
#
# 动机(实测): 同一表面索引算法在三组几何上给出三组左表面索引
# ([55,56,60] / [58,52,57] / [50,51,56]); 其中一组直接使 aitranss 终止, 报
# "your LEFT and RIGHT electrodes share the same atoms"(v0.4 D11)。
# 依据: v0.4 §2.1 P0-4 四项; 元素对下限规则见 _bond_floor(替代 v0.1 过松的 "<1.0 A")。
# 注: "每侧电极沿输运轴的层数 >= N" 中的 N 在 v0.4 §6 标为**无源/待定**, 故此处
#     只报告延伸量(info), 不设硬阈值 —— 不臆造判据。
# ---------------------------------------------------------------------------
COVALENT_RADII = {'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66,
                  'S': 1.05, 'Au': 1.36, 'Ag': 1.45, 'Cu': 1.32, 'Ni': 1.24}
MIN_BOND_FACTOR = 0.7        # d < 0.7*(r_i+r_j) 判为重叠
AU_S_BAND = (2.40, 2.49)     # 有源带: Thygesen 2005 p.12 dAu-S=2.40; Jackson 2009 2.42+-0.02
PAIR_CUTOFF = 3.0            # 只对 3 A 内的原子对做重叠判定


def _bond_floor(e1: str, e2: str):
    """该元素对的键长下限; 半径未知返回 None(不判, 不臆断)。"""
    r1, r2 = COVALENT_RADII.get(e1), COVALENT_RADII.get(e2)
    if r1 is None or r2 is None:
        return None
    return MIN_BOND_FACTOR * (r1 + r2)


def read_geometry_in(path):
    """读 FHI-aims geometry.in -> {'atoms': [(x,y,z,el)], 'meta': {key: [tokens]}}。"""
    atoms, meta = [], {}
    for raw in Path(path).read_text(encoding='utf-8', errors='ignore').splitlines():
        parts = raw.split()
        if not parts:
            continue
        if parts[0] == 'atom' and len(parts) >= 5:
            atoms.append((float(parts[1]), float(parts[2]), float(parts[3]), parts[4]))
        elif len(parts) >= 2:
            meta[parts[0]] = parts[1:]
    return {'atoms': atoms, 'meta': meta}


def _dist(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def validate_junction_geometry(run_dir, metadata=None, geometry=None) -> dict:
    """建结几何硬闸门: ①表面索引合法/左右不相交 ②索引在输运轴正确一侧
    ③接触重叠(元素对+键长下限) ④锚定 S-Au 有源带; 另报告电极延伸量(info)。

    geometry/metadata: 显式给出产物路径(引擎提交 aims_scf 等下游阶段时, 结几何与
    electrode_metadata.json 都落在上游 au60_build 目录, 当前阶段目录里没有);
    缺省仍按 run_dir/geometry.in + run_dir/electrode_metadata.json 查找(旧行为不变)。"""
    run_dir = Path(run_dir)
    geo_p = Path(geometry) if geometry is not None else run_dir / 'geometry.in'
    if not geo_p.exists():
        return {'passed': False, 'errors': 1, 'warnings': 0,
                'checks': [_err('junction_geom_inputs', None,
                                'geometry.in 缺失: %s' % geo_p)]}
    if metadata is None:
        mp = run_dir / 'electrode_metadata.json'
        metadata = json.loads(mp.read_text(encoding='utf-8')) if mp.exists() else None
    if not metadata:
        return {'passed': False, 'errors': 1, 'warnings': 0,
                'checks': [_err('junction_geom_metadata', None,
                                'electrode_metadata.json 缺失: 无法校验表面索引')]}

    atoms = read_geometry_in(geo_p)['atoms']
    n = len(atoms)
    elems = [a[3] for a in atoms]
    checks = []

    # ① 表面索引: 非空 / 在范围 / 左右不相交 / 均为 Au
    L = list(metadata.get('left_surface_indices') or [])
    R = list(metadata.get('right_surface_indices') or [])
    bad = []
    if not L or not R:
        bad.append('表面索引为空')
    oob = [i for i in L + R if not (1 <= int(i) <= n)]
    if oob:
        bad.append('索引越界: %s' % oob)
    inter = sorted(set(L) & set(R))
    if inter:
        bad.append('左右表面索引相交: %s' % inter)
    nonau = [i for i in L + R if 1 <= int(i) <= n and elems[int(i) - 1] != 'Au']
    if nonau:
        bad.append('表面原子非 Au: %s' % nonau)
    if bad:
        checks.append(_err('junction_geom_surface_valid', {'left': L, 'right': R},
                           '; '.join(bad) + ' —— aitranss 会因此报 '
                           'LEFT and RIGHT electrodes share the same atoms'))
    else:
        checks.append(_ok('junction_geom_surface_valid', {'left': L, 'right': R},
                          '左右表面索引合法且不相交(左 %d / 右 %d 个, 共 %d 原子)'
                          % (len(L), len(R), n)))

    # ② 输运轴与左右侧一致性: 轴取两锚定原子连线的主方向(不预设 x/y/z)
    lk = metadata.get('left_linker_index')
    rk = metadata.get('right_linker_index')
    if lk and rk and 1 <= int(lk) <= n and 1 <= int(rk) <= n:
        sl = atoms[int(lk) - 1][:3]
        sr = atoms[int(rk) - 1][:3]
        d = [sr[i] - sl[i] for i in range(3)]
        ax = max(range(3), key=lambda i: abs(d[i]))
        u = [0.0, 0.0, 0.0]
        u[ax] = 1.0 if d[ax] >= 0 else -1.0      # 左锚 -> 右锚 的单位方向
        mid = [(sl[i] + sr[i]) / 2.0 for i in range(3)]

        def side(idx):                      # >0 表示偏向"右锚"一侧
            a = atoms[int(idx) - 1][:3]
            return (a[ax] - mid[ax]) * u[ax]

        wrong = [i for i in L if side(i) >= 0] + [i for i in R if side(i) <= 0]
        axis_name = 'xyz'[ax]
        if wrong:
            checks.append(_err('junction_geom_surface_side', {'axis': axis_name, 'wrong': wrong},
                               '表面索引落在输运轴(%s)错误一侧: %s —— 该索引不稳定正是 '
                               'aitranss "share the same atoms" 的成因' % (axis_name, wrong)))
        else:
            checks.append(_ok('junction_geom_surface_side', {'axis': axis_name},
                              '左右表面索引均位于输运轴(%s)正确一侧' % axis_name))
    else:
        checks.append(_warn('junction_geom_surface_side', None,
                            '缺 left/right_linker_index, 无法定向输运轴'))

    # ③ 接触重叠: 元素对 + 共价半径下限(替代过松的"任一原子对 <1.0 A")
    overlaps = []
    for i in range(n):
        for j in range(i + 1, n):
            floor = _bond_floor(elems[i], elems[j])
            if floor is None:
                continue
            dij = _dist(atoms[i][:3], atoms[j][:3])
            if dij < PAIR_CUTOFF and dij < floor:
                overlaps.append({'i': i + 1, 'j': j + 1, 'elements': elems[i] + '-' + elems[j],
                                 'd': round(dij, 3), 'floor': round(floor, 3)})
    if overlaps:
        checks.append(_err('junction_geom_contacts', overlaps[:8],
                           '%d 对原子短于元素对键长下限(0.7x共价半径和): %s'
                           % (len(overlaps), overlaps[:3])))
    else:
        checks.append(_ok('junction_geom_contacts', None,
                          '%d 原子无重叠(元素对下限 + 3 A 预筛; C-H 1.09 合法、Au-H 1.16 属重叠)' % n))

    # ④ 锚定几何: 锚定 S 与最近的 Au 距离 vs 有源带(仅方向性 -> warn)
    sulf = [i for i, e in enumerate(elems) if e == 'S']
    aus = [i for i, e in enumerate(elems) if e == 'Au']
    if sulf and aus:
        best = min(((_dist(atoms[s][:3], atoms[a][:3]), s + 1, a + 1)
                    for s in sulf for a in aus), key=lambda t: t[0])
        dmin, si, ai = best
        lo, hi = AU_S_BAND
        if lo <= dmin <= hi:
            checks.append(_ok('junction_geom_anchor_aus', round(dmin, 4),
                              'S%d-Au%d = %.4f A 落在有源带 %.2f-%.2f A (Thygesen 2005 p.12; '
                              'Jackson 2009)' % (si, ai, dmin, lo, hi)))
        else:
            checks.append(_warn('junction_geom_anchor_aus', round(dmin, 4),
                                'S%d-Au%d = %.4f A 偏离有源带 %.2f-%.2f A(仅方向性: 两锚点均为'
                                '表面/金字塔配位, 不覆盖单 Au 低配位)' % (si, ai, dmin, lo, hi)))
    else:
        checks.append(_warn('junction_geom_anchor_aus', None, '缺 S 或 Au, 无法判锚定几何'))

    # 附: 电极沿输运轴延伸量(info; 下限 N 待定, 见 v0.4 §6, 不设硬阈值)
    if lk and rk and 1 <= int(lk) <= n and 1 <= int(rk) <= n:
        tips = {metadata.get('left_tip_index'), metadata.get('right_tip_index')}
        ext = {}
        for name, sel in (('left', [i + 1 for i, e in enumerate(elems)
                                    if e == 'Au' and (i + 1) not in tips and side(i + 1) < 0]),
                          ('right', [i + 1 for i, e in enumerate(elems)
                                     if e == 'Au' and (i + 1) not in tips and side(i + 1) > 0])):
            ext[name] = round(max((abs(side(i)) for i in sel), default=0.0), 3)
        checks.append(_ok('junction_geom_electrode_extent', ext,
                          '电极沿输运轴(%s)延伸: 左 %.2f A / 右 %.2f A(仅报告; 下限 N 无源待定)'
                          % (axis_name, ext['left'], ext['right'])))

    # ⑤ P0-6: A2 人工确认与几何指纹绑定 —— 无记录 warn, 指纹不符 error(须重新确认)
    from .fingerprint import check_confirmation
    conf = check_confirmation(run_dir)
    if conf.get('valid'):
        checks.append(_ok('junction_a2_confirmation', conf.get('current'), conf['reason']))
    elif conf.get('present'):
        checks.append(_err('junction_a2_confirmation', conf.get('current'), conf['reason']))
    else:
        checks.append(_warn('junction_a2_confirmation', conf.get('current'), conf['reason']))

    errs = sum(1 for c in checks if c['severity'] == 'error' and not c['passed'])
    warns = sum(1 for c in checks if c['severity'] == 'warn' and not c['passed'])
    return {'passed': errs == 0, 'errors': errs, 'warnings': warns, 'checks': checks}