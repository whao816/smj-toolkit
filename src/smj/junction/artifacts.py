"""P1-7 (AGD-2026-09-25-SMJ): smj 产物回取与存储治理。

动机(实测): 本会话最终 run 的**集群侧**产物体积 >2.2 GB ——
  omat.aims 0.48 GB + mos.aims 0.96 GB + aims.restart 0.76 GB + Mulliken.out 0.089 GB。
其中 `mos.aims`/`omat.aims` 本地只为**读表头**(nsaos / 三角布局自检), 全量回取纯属浪费;
驱动一度卡在下载 1.4 GB。已验证的省法: 只取远端表头 `pull_head`(本会话实测可行)。

本模块只做**策略与判据**, 不做真实传输。
"""
from __future__ import annotations

from pathlib import Path

# 回取模式: full(全量) | head(只取表头 N 行) | skip(仅留远端)
# 依据(实测): mos.aims 仅用于读 nsaos; omat.aims 仅用于三角布局自检(行数=nsaos*(nsaos+1)/2);
# aims.restart / Mulliken.out 属"仅远端留存"(需要时可另取)。
SMJ_ARTIFACTS = {
    # 热导链(D12:独立于电子六阶段)。thermal_fc 只回取 forces_all.json:
    # 数组作业的收尾步骤已经把 N 份 fc/NNNN/forces.json 拼好,fc/ 整目录
    # (每个位移点一份 aims.out,300+ 份)留在远端。
    'thermal_device': {'must': ['geometry.in', 'device_metadata.json',
                                '.dsh_done'],
                       'head': {}, 'skip': []},
    'thermal_relax': {'must': ['geometry.in', 'geometry.in.next_step',
                               'control.in', 'aims.out', '.dsh_done'],
                      'head': {}, 'skip': []},
    'thermal_fc': {'must': ['displacements.json', 'control.in', 'geometry.in',
                            'forces_all.json', '.dsh_done'],
                   'head': {}, 'skip': ['fc']},
    'thermal_analyze': {'must': ['thermal.json', '.dsh_done'],
                        'head': {}, 'skip': []},
    'gaussian_opt': {'must': ['mol.gjf', 'mol.log', 'sub_g16.sh', '.dsh_done'],
                     'head': {}, 'skip': []},
    'aims_tip_relax': {'must': ['geometry.in', 'geometry.in.next_step', 'control.in',
                                'aims.out', 'sub_aims.sh', '.dsh_done'],
                       'head': {}, 'skip': ['aims.restart']},
    'au_relax': {'must': ['geometry.in', 'geometry.in.next_step', 'control.in',
                          'aims.out', 'sub_aims.sh', '.dsh_done'],
                 'head': {}, 'skip': ['aims.restart', 'Mulliken.out']},
    'aims_scf': {'must': ['geometry.in', 'control.in', 'aims.out', 'basis-indices.out',
                          'sub_aims.sh', '.dsh_done'],
                 'head': {}, 'skip': ['aims.restart', 'Mulliken.out']},
    'aitranss_te': {'must': ['tcontrol', 'TE.dat', 'te.out', 'self.energy.in', 'sub_te.sh',
                             '.dsh_done'],
                    'head': {}, 'skip': []},
}

# FHI-aims 的大文件: 本地只需要**表头**(nsaos 等), 全量回取无收益。
# 实测体积: mos.aims 0.96 GB / omat.aims 0.48 GB -> 取前 2 行即可满足本地解析需要。
AIMS_HEAD_ONLY = {'mos.aims': 2, 'omat.aims': 2}

# 只有真正产出这两个文件的阶段才应把它们列进 head；否则 gaussian_opt /
# aitranss_te 的回取计划会平白多出两个必然 missing 的文件（2026-11-11 本地补丁 F1）。
AIMS_HEAD_STAGES = ('aims_scf',)

# 单 run 体积告警门限(GiB)。**工程约定, 非物理判据** —— 目的是尽早暴露存储膨胀,
# 本会话单 run 集群侧实测已 >2.2 GB。
VOLUME_WARN_GIB = 2.0


def retrieval_plan(stage: str) -> dict:
    """该阶段的回取计划 {must, head, skip, unknown}。

    `mos.aims`/`omat.aims` 自动归入 head(只取表头), 不因阶段表未列出而漏掉。
    """
    spec = SMJ_ARTIFACTS.get(stage)
    if spec is None:
        return {'unknown': True, 'stage': stage,
                'reason': '未知 smj 阶段(可用: %s)' % sorted(SMJ_ARTIFACTS)}
    must = list(spec.get('must') or [])
    head = dict(spec.get('head') or {})
    if stage in AIMS_HEAD_STAGES:
        head.update(AIMS_HEAD_ONLY)
    skip = list(spec.get('skip') or []) + [k for k in AIMS_HEAD_ONLY if k not in head]
    unknown = [k for k in head if k in must]
    return {'must': must, 'head': head, 'skip': sorted(set(skip)), 'unknown': unknown,
            'stage': stage}


def retrieval_mode(stage: str, name: str) -> str:
    """单文件回取模式: 'full' | 'head' | 'skip' | 'unknown'。"""
    p = retrieval_plan(stage)
    if p.get('unknown'):
        return 'unknown'
    if name in p['head']:
        return 'head'
    if name in p['must']:
        return 'full'
    if name in p['skip']:
        return 'skip'
    return 'skip'          # 未登记文件默认仅留远端(避免误拉 GB 级中间件)


def estimate_transfer_bytes(sizes: dict, stage: str) -> dict:
    """按 sizes{name: bytes} 估算实际需传输字节, 并给出省下的量。

    只计 must(full) 与 head(按 head 行数比例粗估 —— **工程近似**, 非精确传输模型)。
    """
    plan = retrieval_plan(stage)
    if plan.get('unknown'):
        return {'unknown': True, 'reason': plan['reason']}
    total = sum(sizes.values())
    xfer = 0
    for name, n in sizes.items():
        mode = retrieval_mode(stage, name)
        if mode == 'full':
            xfer += n
        elif mode == 'head':
            # 表头占比极低; 按 min(1%, 100 KB) 粗估(工程近似)
            xfer += min(max(int(n * 0.01), 1024), 100 * 1024)
    return {'stage': stage, 'total_bytes': total, 'transfer_bytes': xfer,
            'saved_bytes': total - xfer,
            'saved_ratio': (1.0 - xfer / total) if total else 0.0}


def volume_warning(total_bytes: float, limit_gib: float = VOLUME_WARN_GIB) -> dict:
    """单 run 体积告警(**工程约定**, 非物理判据)。"""
    limit = limit_gib * (1024 ** 3)
    return {'ok': total_bytes <= limit, 'total_gib': round(total_bytes / (1024 ** 3), 3),
            'limit_gib': limit_gib,
            'detail': ('单 run 体积 %.3f GiB 超告警门限 %.1f GiB(工程约定)'
                       % (total_bytes / (1024 ** 3), limit_gib))
            if total_bytes > limit else
            ('单 run 体积 %.3f GiB 未超 %.1f GiB 门限'
             % (total_bytes / (1024 ** 3), limit_gib))}


# ---------------------------------------------------------------------------
# P0-5 余项 (AGD-2026-09-25-SMJ): smj 阶段产物契约 —— 供 resume 判定"可复用"
#
# 证据(本会话):
# - 阶段序: stages/驱动实际执行的 gaussian_opt -> aims_tip_relax -> au60_build
#   -> aims_scf -> aitranss_te;
# - au60_build 的输入就是上游的 geometry.in.next_step(stages.py:1447);
# - aims_scf 的 geometry.in 与 au60_build 产物 **sha256 完全相同**(本轮 manifest 实证
#   `1de1783a43bf`), 即 SCF 原样消费建结几何;
# - aitranss_te 消费 SCF 目录内的 omat/mos/basis + tcontrol + self.energy.in(驱动实测:
#   self.energy.in 必须在 upload 之前生成, 否则 aitranss 报左右电极共享原子)。
# ---------------------------------------------------------------------------
SMJ_STAGE_ORDER = ('gaussian_opt', 'aims_tip_relax', 'au60_build', 'aims_scf', 'aitranss_te')

SMJ_STAGE_CONTRACT = {
    'gaussian_opt': {'consumes': ['mol.gjf'],
                     'produces': ['mol.log']},
    'aims_tip_relax': {'consumes': ['geometry.in'],
                       'produces': ['geometry.in.next_step', 'aims.out']},
    # au60_build 是**本地**阶段(stages.py:1442 needs_submit=False) -> 不适用提交谓词
    'au60_build': {'consumes': ['geometry.in.next_step'],
                   'produces': ['geometry.in', 'electrode_metadata.json'],
                   'local': True},
    'aims_scf': {'consumes': ['geometry.in', 'control.in'],
                 'produces': ['aims.out', 'basis-indices.out']},
    'aitranss_te': {'consumes': ['tcontrol', 'self.energy.in'],
                    'produces': ['TE.dat', 'te.out']},
}


def reusable(stage: str, step_dir) -> dict:
    """该阶段产物是否可复用(resume 时跳过重算)。

    判据**复用 P0-3 的谓词**(junction.submit.smj_finished: 哨兵优先 + 成功标记 +
    扫到 $eend), **不另立一套**; 另要求契约里的 produces 全部存在。
    """
    d = Path(step_dir)
    spec = SMJ_STAGE_CONTRACT.get(stage)
    if spec is None:
        return {'ok': False, 'stage': stage, 'reason': '未知 smj 阶段(无契约)'}
    missing = [f for f in spec['produces'] if not (d / f).exists()]
    if missing:
        return {'ok': False, 'stage': stage, 'missing': missing,
                'reason': '产物缺失: %s' % ', '.join(missing)}
    if spec.get('local'):
        # 本地阶段无作业/无哨兵 -> 判据只能是"产物齐备"(此处已通过), 不得套用提交谓词
        return {'ok': True, 'stage': stage,
                'reason': '本地阶段(needs_submit=False): 产物齐备即视为可复用(不适用提交谓词)'}
    try:
        from .submit import smj_finished
        finished = bool(smj_finished(d))
    except Exception as exc:
        return {'ok': False, 'stage': stage, 'reason': '完成判定失败: %s' % exc}
    return {'ok': finished, 'stage': stage,
            'reason': '产物齐备且通过 P0-3 完成判定' if finished
                      else '产物齐备但**未通过**完成判定(哨兵/成功标记/终点检查)'}


def resume_plan(run_root, upto: str | None = None) -> dict:
    """按阶段序给出每个阶段的复用性; upto=从该阶段起续跑(之前的阶段须可复用)。"""
    root = Path(run_root)
    plan = {}
    for st in SMJ_STAGE_ORDER:
        d = root / st
        plan[st] = reusable(st, d) if d.is_dir() else {
            'ok': False, 'stage': st, 'reason': '阶段目录不存在'}
    if upto is None:
        return {'plan': plan}
    if upto not in SMJ_STAGE_ORDER:
        raise ValueError('upto=%r 不在阶段序中: %s' % (upto, list(SMJ_STAGE_ORDER)))
    idx = SMJ_STAGE_ORDER.index(upto)
    blocked = [st for st in SMJ_STAGE_ORDER[:idx] if not plan[st]['ok']]
    return {'plan': plan, 'resume_from': upto,
            'can_resume': not blocked,
            'reason': ('上游 %d 个阶段均可复用, 可从 %s 续跑' % (idx, upto)) if not blocked
                      else ('上游阶段不可复用, 不能直接从 %s 续跑: %s' % (upto, blocked))}
