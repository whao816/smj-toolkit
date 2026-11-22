"""P2-2 (AGD-2026-09-25-SMJ): smj 运行前预检(preflight)。

动机(本会话五条一手证据):
1. species 需**手工**取齐 —— 我手工 scp 了 H/C/S/Au 的 light 文件; 若缺任一元素,
   FHI-aims 会在远端才失败(浪费一次提交 + 排队)。
2. `AIMS_MODULE` 非空时会**遮蔽** `aims_bin`(模板走 module 分支选二进制)
   -> 本会话据此误用了 GPU 版而非任务卡要求的 scalapack CPU 版。
   2026-10-24(缺口⑦): `AIMS_BIN`/`AIMS_LD_MODULE`/`MKL_LIB` 已进 `DEFAULT_ENV`;
   直接二进制仍须 `AIMS_LD_MODULE` 提供运行时库(该 scalapack 版链 libcublas/libmkl/libimf)。
3. 集群路径由 `inputs.py::DEFAULT_ENV` 在 **import 时**读 `os.environ`
   -> "先设 env 再 import" 是隐性前提, 值得显式检查。
4. `mos.aims`/`omat.aims` 体积大(0.96/0.48 GB), 本地只需表头(见 artifacts.py)。
5. 电极规模由 `au_per_electrode` 决定, 与 `self.energy.in` 的电极原子集必须一致
   (本会话曾报 `Expected 18 electrode Au, found 118`)。

纯检查, 不触网、不提交。
"""
from __future__ import annotations

import os
import re
from pathlib import Path

SPECIES_RE = re.compile(r'^(\d+)_([A-Z][a-z]?)_default$')


def species_available(species_dir, elements) -> dict:
    """检查 species 目录是否覆盖给定元素集合(命名约定 NN_El_default)。

    返回 {'ok', 'missing', 'found'}; 目录不存在 -> ok=False 且 missing=全部。
    """
    d = Path(species_dir) if species_dir else None
    have = {}
    if d and d.is_dir():
        for f in d.iterdir():
            m = SPECIES_RE.match(f.name)
            if m:
                have.setdefault(m.group(2), f.name)
    want = sorted(set(elements or []))
    missing = [e for e in want if e not in have]
    return {'ok': not missing, 'missing': missing,
            'found': {e: have[e] for e in want if e in have},
            'detail': ('species 齐备: %s' % ', '.join(want)) if not missing else
                      ('species 缺失: %s(需从集群 species_defaults 取齐; 本会话为手工 scp)'
                       % ', '.join(missing))}


def aims_executable(env=None) -> dict:
    """解析 FHI-aims 可执行: aims_bin(AIMS_BIN) 与 module(AIMS_MODULE) 的**遮蔽**关系。

    本会话实测: `AIMS_MODULE` 非空时模板优先走 module 分支, 使 `aims_bin` 被遮蔽
    -> 用了 GPU 版而非要求的 scalapack CPU 版。故此处显式报告二者, 并在同时存在时告警。
    2026-10-24(缺口⑦): 直接二进制(aims_bin 分支)仍须 module 提供运行时库, 由
    AIMS_LD_MODULE 单独声明 —— 它只设 LD_LIBRARY_PATH, 不选二进制, 不会遮蔽 aims_bin。
    """
    e = dict(os.environ if env is None else env)
    binp = (e.get('AIMS_BIN') or '').strip()
    mod = (e.get('AIMS_MODULE') or '').strip()
    ldmod = (e.get('AIMS_LD_MODULE') or '').strip()
    if binp and mod:
        return {'ok': False, 'aims_bin': binp, 'aims_module': mod, 'aims_ld_module': ldmod,
                'detail': ('同时设置了 AIMS_BIN 与 AIMS_MODULE —— module 分支会**遮蔽** '
                           'aims_bin(本会话实测); 若要用直接二进制请清空 AIMS_MODULE')}
    if not binp and not mod:
        return {'ok': False, 'aims_bin': '', 'aims_module': '', 'aims_ld_module': ldmod,
                'detail': 'AIMS_BIN 与 AIMS_MODULE 均未设置(至少给一个: AIMS_BIN=直接二进制路径)'}
    detail = 'FHI-aims 路径确定: %s' % (binp or ('module ' + mod))
    if binp and ldmod:
        detail += '(运行时库 module load %s)' % ldmod
    elif binp:
        detail += '(未设 AIMS_LD_MODULE: 若该二进制链 CUDA/MKL 运行时库, 作业内会缺库失败)'
    return {'ok': True, 'aims_bin': binp, 'aims_module': mod, 'aims_ld_module': ldmod,
            'detail': detail}


def cluster_paths(settings) -> dict:
    """集群路径/分区是否足以提交(remote_root + ssh_alias|ssh_host + partition)。"""
    c = getattr(settings, 'cluster', None)
    problems = []
    if c is None:
        return {'ok': False, 'detail': 'settings.cluster 缺失'}
    if not getattr(c, 'remote_root', None):
        problems.append('remote_root 缺失(远端工作根目录)')
    if not (getattr(c, 'ssh_alias', None) or getattr(c, 'ssh_host', None)):
        problems.append('ssh_alias/ssh_host 均缺失(无法连集群)')
    if not getattr(c, 'partition', None):
        problems.append('partition 缺失')
    return {'ok': not problems, 'remote_root': getattr(c, 'remote_root', None),
            'partition': getattr(c, 'partition', None),
            'detail': ('集群路径齐备: %s @ %s' % (getattr(c, 'remote_root', None),
                                                 getattr(c, 'partition', None)))
            if not problems else '; '.join(problems)}


# ---------------------------------------------------------------------------
# P1 前置 checklist(taskcard-smj-transport §5 第 1 条): 四项集群资源确认
#   Gaussian 许可 / FHI-aims species 路径 / aitranss 二进制 / sbatch。
# species+cluster+aims 见上方; gaussian/aitranss/sbatch 补充于下。
# 纯检查, 不触网、不提交 —— 目的: 在**提交前**发现资源未就绪, 而非排完队才失败。
# ---------------------------------------------------------------------------

#: inputs.py 在缺省时写入的占位值 —— 不得被当成「已配置」。
_PLACEHOLDER = 'UNCONFIGURED'


def _is_configured(v) -> bool:
    """值已配置(非空且非 UNCONFIGURED 占位)。"""
    s = str(v or '').strip()
    return bool(s) and s.upper() != _PLACEHOLDER


def gaussian_license(env=None) -> dict:
    """Gaussian 许可/安装是否足以提交(taskcard §5 第 1 项)。

    集群上 Gaussian 以两种方式可用(见 tools/junction/inputs.py sub_g16.sh 模板):
      ① module 方式: G16_MODULE 非空 -> module load;
      ② g16root 方式: 直接 source $g16root/bsd/g16.profile。
    许可文件位置由 GAUSS_LIC 或 G16_SCRATCH(g16 把 lic 放在 scratch 的部署惯例)指示。
    未配置一律 fail-closed; 两者并存**不是**错误(module 优先), 但必须显式告知。
    """
    e = dict(os.environ if env is None else env)
    root = (e.get('G16ROOT') or '').strip()
    mod = (e.get('G16_MODULE') or '').strip()
    lic = (e.get('GAUSS_LIC') or '').strip()
    scratch = (e.get('G16_SCRATCH') or '').strip()
    has_root, has_mod = _is_configured(root), _is_configured(mod)
    if not has_root and not has_mod:
        return {'ok': False, 'g16root': root, 'g16_module': mod,
                'license_from': '',
                'detail': 'G16ROOT 与 G16_MODULE 均未设置 —— Gaussian 无法执行(gjf/freq 阶段会失败)'}
    license_from = ('GAUSS_LIC' if _is_configured(lic)
                    else ('G16_SCRATCH' if _is_configured(scratch) else ''))
    note = ''
    if has_root and has_mod:
        note = ('G16_MODULE 非空时走 module 分支, g16root 路径方式被遮蔽'
                '(与 FHI-aims 的 AIMS_MODULE 遮蔽 AIMS_BIN 同类); 当前实际走: module ' + mod)
    return {'ok': True, 'g16root': root, 'g16_module': mod,
            'license_from': license_from, 'note': note,
            'detail': 'Gaussian 就绪: %s%s' % (('module ' + mod) if has_mod else ('g16root ' + root),
                                             ('' if license_from else '(未显式给许可路径)'))}


def aitranss_binary(env=None) -> dict:
    """aitranss 透射计算二进制是否配置(taskcard §5 第 1 项)。

    缺省值 UNCONFIGURED 是 inputs.py 的占位符, 会被 validate 拦截 —— 此处提前暴露。
    """
    e = dict(os.environ if env is None else env)
    b = (e.get('AITRANSS_BIN') or '').strip()
    if not _is_configured(b):
        return {'ok': False, 'aitranss_bin': b,
                'detail': ('AITRANSS_BIN 未设置或为占位值 %r —— aitranss_te 阶段无法执行'
                           % (b or ''))}
    return {'ok': True, 'aitranss_bin': b, 'detail': 'aitranss 就绪: ' + b}


def sbatch_available(settings=None, *, live: bool = False) -> dict:
    """sbatch 可提交(taskcard §5 第 1 项), **离线优先**。

    设计约束(本模块文档头): preflight 是「纯检查, 不触网、不提交」。而远程调度后端
    (`slurm-ssh`/`pbs-ssh`/`aws-batch`)的 probe 会 SSH 跑 `hostname && sbatch --version`,
    属真实网络往返 —— 故缺省**不调用**远程 probe, 按可确认性分级:

      * 本地调度后端(`slurm`/`pbs'): probe 只做本地 `shutil.which`, 直接调用, severity=error;
      * 非调度后端(`local`): 不经 sbatch, severity=info 且不阻断;
      * 远程后端: 缺省不触网 -> severity=warn(可确认缺口, 不阻断), 提示用 `live=True` 真探测;
        `live=True` 时才调用 probe 并据结果给 error/info。

    `live` 是显式开关: 调用方(如人工确认或诊断命令)自担一次网络往返。
    """
    from smj.compat import get_backend
    try:
        backend = get_backend(settings)
    except Exception as exc:  # noqa: BLE001 —— 后端不可用本身即未就绪
        return {'ok': False, 'severity': 'error', 'detail': '计算后端不可用: %s' % exc}
    name = getattr(backend, 'name', '')
    remote = bool(getattr(backend, 'remote', False))
    # 非调度后端: 本地直跑, 不涉及 sbatch
    if name not in ('slurm', 'slurm-ssh', 'pbs', 'pbs-ssh', 'aws-batch'):
        return {'ok': True, 'severity': 'info', 'backend': name, 'sbatch': '',
                'detail': '后端 %s 不经 sbatch(local 直跑), 本项不适用' % name}
    if remote and not live:
        return {'ok': False, 'severity': 'warn', 'backend': name, 'sbatch': '',
                'detail': ('远程后端 %s 的 sbatch 在集群侧; 离线预检不触网, 未探测 —— '
                           '如需真实确认请加 --live(会产生一次 SSH 往返)' % name)}
    probe = getattr(backend, 'probe', None)
    if not callable(probe):
        return {'ok': False, 'severity': 'error', 'backend': name,
                'detail': '后端 %s 无 probe 能力, 无法确认 sbatch' % name}
    try:
        r = probe() or {}
    except Exception as exc:  # noqa: BLE001
        return {'ok': False, 'severity': 'error', 'backend': name,
                'detail': 'sbatch 探测失败: %s' % exc}
    ok = bool(r.get('ok'))
    slurm = str(r.get('slurm') or r.get('sbatch') or '')
    return {'ok': ok, 'severity': 'error' if not ok else 'info', 'backend': name,
            'sbatch': slurm,
            'detail': ('sbatch 就绪: ' + (slurm or 'ok')) if ok
                      else ('sbatch 不可用(%s): %s' % (name,
                                                       r.get('error') or r.get('sbatch') or 'not found'))}


def _merge(name, r, default_severity='error') -> dict:
    """合并检查项: 检查函数可自带 severity(如 sbatch 的 warn/info), 未带则用默认。

    severity 语义: error=阻断(未就绪); warn/info=可确认缺口或说明, 不阻断。
    """
    d = dict(r or {})
    sev = d.pop('severity', None) or default_severity
    return {'name': name, 'severity': sev, **d}


def _checklist_checks(settings=None, *, elements=None, species_dir=None, env=None,
                      live=False) -> list:
    """taskcard §5 四项 + cluster/aims/species 的检查项汇总(顺序即报告顺序)。"""
    checks = [
        _merge('cluster_paths', cluster_paths(settings)),
        _merge('gaussian_license', gaussian_license(env)),
        _merge('aitranss_binary', aitranss_binary(env)),
        _merge('sbatch', sbatch_available(settings, live=live)),
        _merge('aims_executable', aims_executable(env)),
    ]
    if elements and species_dir:
        checks.append(_merge('species_complete', species_available(species_dir, elements)))
    else:
        # 未提供输入 -> warn(不阻断): 「没问」不等于「缺资源」, 不得据此判未就绪。
        # 提供了 species_dir 却缺元素才是 error(那才是真缺口)。
        checks.append({'name': 'species_complete', 'severity': 'warn', 'ok': False,
                       'detail': '未提供 elements/species_dir, 跳过 species 检查(不臆断)'})
    return checks


def resource_checklist(settings=None, *, elements=None, species_dir=None, env=None,
                       live=False) -> dict:
    """taskcard §5 四项集群资源确认的单一入口(供 CLI/driver/测试复用)。

    在提交前一次性判定: cluster 路径 / Gaussian 许可 / aitranss 二进制 / sbatch /
    FHI-aims 可执行 / species 齐备。任一项 error 未过即 ok=False(fail-closed)。
    """
    checks = _checklist_checks(settings, elements=elements, species_dir=species_dir, env=env,
                              live=live)
    errors = sum(1 for c in checks if c['severity'] == 'error' and not c['ok'])
    warns = sum(1 for c in checks if c['severity'] == 'warn' and not c['ok'])
    return {'ok': errors == 0, 'errors': errors, 'warnings': warns, 'checks': checks}


def smj_preflight(settings, elements=None, species_dir=None, env=None) -> dict:
    """smj 运行前预检: 汇总为 {'ok', 'errors', 'checks'}。

    elements/species_dir 由调用方按当前分子给出(未知则不检查, 不臆断)。

    P1 checklist 扩展(2026-10-13): 除原有 cluster/aims/species 外, 补齐 taskcard §5
    第 1 条的四项集群资源确认 —— Gaussian 许可 / FHI-aims species / aitranss 二进制 / sbatch。
    现为 resource_checklist 的别名(同一实现, 保留本名向后兼容既有调用方)。
    """
    return resource_checklist(settings, elements=elements, species_dir=species_dir, env=env)


# ---------------------------------------------------------------------------
# P2-3 (AGD-2026-09-25-SMJ): 最小输入 -> 完整规格的自动推导
#
# 目标形态(v0.4 §0): 用户只给 分子结构 + 电极 + 锚定点 [+ 温度], 其余自动推导。
# **纪律**: 每个默认值都注明来源; 凡无来源者标为"项目既往约定", 不得伪装成文献值。
# 特别地 `au_tip_distance` 默认 **2.5 A** —— 有源**判据带** 2.40-2.44 A 是**弛豫产物**
# 判据, **不是**输入放置距离的依据(这正是 v0.1 被否定的"用错对象")。
# ---------------------------------------------------------------------------
AU_TIP_DISTANCE_DEFAULT = 2.5
#: W1-4 (2026-11-05): 出账时的**实测利用率折扣**。理论槽位是"节点全空 + 队列不停"的上限,
#: 真实跑里有排队、审计等待、人工审批、被别人的作业挤占(实测 gpu1 上长期有 adft-*/l502
#: 常驻)。W0 四臂与 σ 扫描 12 作业的实际占用率约 0.7, 故日吞吐按 0.70 折算并显式标注。
ACCOUNTING_DERATE = 0.70
_DEFAULT_PROVENANCE = {
    'au_tip_distance': '2.5 A —— **项目既往约定, 非文献依据**(判据带 2.40-2.44 A 判的是弛豫产物, 不是输入放置距离)',
    'electrode_element': "Au —— 与改动前的 Electrode111 一致(默认路径逐字节不变)",
    'electrode_lattice_a': 'None —— 用 Electrode111 字面常量; 显式给 a 才按公式重算',
    'temperature_K': '300 —— 默认室温(判据 T(E_F)/S 均在室温口径)',
}


def derive_spec(molecule_geometry, electrode=None, anchors=None, temperature: float = 300.0,
                species_dir=None, settings=None, au_per_electrode: int | None = None,
                stages=None) -> dict:
    """从最小输入推导计算规格; 返回 {'spec', 'provenance', 'checks'}。

    molecule_geometry: 分子结构(xyz/geometry.in) —— 仅用于**读元素集合**与规模估计。
    electrode: {'element': 'Au', 'a_lattice': None|float, 'layers': 6, 'au_per_electrode': N}
    anchors: {'left': i, 'right': j}(可选; 校验在 P0-4 几何闸门做)
    """
    from .electrode import FCC_LATTICE_A
    from .preflight import aims_executable, cluster_paths, species_available
    from smj.compat import estimate_smj_resources, check_smj_budget

    if species_dir is None:
        species_dir = os.environ.get('AIMS_SPECIES_DIR') or None
    el = dict(electrode or {})
    element = el.get('element') or 'Au'
    a_lattice = el.get('a_lattice')
    layers = int(el.get('layers') or 6)
    aper = el.get('au_per_electrode', au_per_electrode)

    # 元素集合 + 规模: 直接读结构文件(不猜)
    elements, n_mol = [], 0
    gp = Path(molecule_geometry)
    if gp.exists():
        for line in gp.read_text(encoding='utf-8', errors='ignore').splitlines():
            f = line.split()
            if len(f) >= 5 and f[0] == 'atom':
                n_mol += 1
                if f[4] not in elements:
                    elements.append(f[4])
            elif len(f) == 4 and f[0] != 'lattice_vector':   # xyz: El x y z
                try:
                    float(f[1]); float(f[2]); float(f[3])
                except ValueError:
                    continue
                n_mol += 1
                if f[0] not in elements:
                    elements.append(f[0])

    checks = []
    if not gp.exists():
        checks.append({'name': 'molecule_geometry', 'severity': 'error', 'ok': False,
                       'detail': '分子结构文件不存在: %s' % gp})
    else:
        checks.append({'name': 'molecule_geometry', 'severity': 'info', 'ok': True,
                       'detail': '分子 %d 原子, 元素 %s' % (n_mol, ', '.join(elements))})
    cp = cluster_paths(settings) if settings is not None else {'ok': False, 'detail': '未给 settings'}
    checks.append({'name': 'cluster_paths', 'severity': 'error', **cp})
    ae = aims_executable()
    checks.append({'name': 'aims_executable', 'severity': 'error', **ae})
    sp = species_available(species_dir, elements) if (elements and species_dir) else {
        'ok': False, 'detail': '未提供 species_dir(或元素为空), 跳过(不臆断)'}
    checks.append({'name': 'species_complete', 'severity': 'warn' if not species_dir else 'error', **sp})

    if element not in FCC_LATTICE_A and a_lattice is None:
        checks.append({'name': 'electrode_element', 'severity': 'error', 'ok': False,
                       'detail': '未登记金属 %r 且未给 a_lattice(不猜晶格常数; 已登记: %s)'
                                 % (element, sorted(FCC_LATTICE_A))})

    # 资源: 逐阶段取**该阶段真实体系规模**。P1-5 教训: 旧实现把同一个 n_est 套给四个阶段,
    # 于是 12 原子分子的 SCF 按 132 原子定档(或反过来), 定档与实际规模无关。
    n_tip_relax = (n_mol + 2) if n_mol else 0              # +add_tip_au 的两颗 Au tip
    aper_eff = int(aper or 60)
    n_junction = (n_mol + 2 * aper_eff) if n_mol else 0     # +两侧电极(含 tip)
    per_stage_atoms = {'gaussian_opt': n_mol, 'aims_tip_relax': n_tip_relax,
                       'aims_scf': n_junction, 'aitranss_te': n_junction}
    resources = {}
    for st in (stages or ('gaussian_opt', 'aims_tip_relax', 'aims_scf', 'aitranss_te')):
        if per_stage_atoms.get(st):
            e = estimate_smj_resources(st, per_stage_atoms[st])
            resources[st] = {**e, 'budget': check_smj_budget(e)}

    # W1-4 (2026-11-05): 出账 —— 单结核时 / 建议并发 / 预计日吞吐。
    # 吞吐不是单作业水位的函数: 一个结要串行跑 SCF 与 TE, 两者资源画像差 4 倍,
    # 真瓶颈是 min(n_scf/t_scf, n_te/t_te)(见 concurrency.plan_slots 的推导)。
    accounting = {}
    if resources:
        from .concurrency import plan_slots
        from smj.compat import walltime_seconds
        specs = {}
        for st, e in resources.items():
            if not isinstance(e, dict) or e.get('local') or e.get('unknown'):
                continue
            # 稳态只对**批量循环里反复跑的两个重阶段**求并发: gaussian/tip_relax 是每分子一次性
            # 前置(0.25 h/1 槽), 若也按"常驻占内存"进水位, 会白扣 20 GiB 并把日吞吐从 ~45 压到 ~36。
            if st not in ('aims_scf', 'aitranss_te'):
                continue
            specs[st] = {'mem': e['mem'], 'cores': e.get('cores') or e.get('ntasks') or 1,
                         'wall_s': walltime_seconds(e['walltime'])}
        plan = plan_slots(specs, derate=ACCOUNTING_DERATE) if specs else {}
        accounting = {
            'per_junction_core_h': round(sum(float(e.get('core_h') or 0)
                                             for e in resources.values()
                                             if isinstance(e, dict)), 2),
            'stage_core_h': {st: e.get('core_h') for st, e in resources.items()
                             if isinstance(e, dict) and e.get('core_h') is not None},
            'slots': plan.get('slots', {}),
            'junctions_per_day': plan.get('junctions_per_day', 0.0),
            'cycle_hours': plan.get('cycle_hours'),
            'bound_by': plan.get('bound_by'),
            'usable_mem_gib': plan.get('usable_mem_gib'),
            'derate': ACCOUNTING_DERATE,
            'note': ('日吞吐 = 24*3600*derate / max(阶段墙钟/槽位); derate=%.2f 含排队/空转'
                     % ACCOUNTING_DERATE)}

    provenance = dict(_DEFAULT_PROVENANCE)
    provenance['species_dir'] = str(species_dir) if species_dir else '未提供'
    provenance['atom_count_estimate'] = ('结 %d = 分子 %d + 2x每侧电极 %d(含 tip); tip_relax %d '
                                        '= 分子 + 2 tip(**估算**, 实际以各阶段 geometry 为准)'
                                        % (n_junction, n_mol, aper_eff, n_tip_relax)
                                        if n_mol else '未知')

    errors = sum(1 for c in checks if c['severity'] == 'error' and not c['ok'])
    return {'ok': errors == 0, 'errors': errors, 'checks': checks,
            'spec': {'molecule': str(gp), 'elements': elements,
                     'anchors': dict(anchors or {}),
                     'electrode': {'element': element, 'a_lattice': a_lattice, 'layers': layers,
                                   'au_per_electrode': aper},
                     'au_tip_distance': AU_TIP_DISTANCE_DEFAULT,
                     'temperature_K': float(temperature)},
            'resources': resources, 'accounting': accounting, 'provenance': provenance}
