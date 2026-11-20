"""P1-8 (AGD-2026-09-25-SMJ): smj 并发水位与内存申请纪律。

动机(两条实测/既有教训):
1. 任务卡 §7.1 **内存纪律**: "slurm `CR_CORE_MEMORY` 按**申请内存**记账(NEB 并发各 32G
   曾占满 96/122G) → 本族作业必须显式 `--mem` + 短 `--time`, 不做整节点默认申请"。
2. 本会话 gpu1 实测: 节点 **64 核 / 124703 MB**; smj 的 SCF 单作业用 `--mem=48G`,
   **内存而非核数**才是并发约束(16 ranks 连 64G 都不够, 降 rank 才过)。

故本模块按**内存**算并发水位, 并把"整节点默认申请"列为违规。
纯计算, 不做真实提交。
"""
from __future__ import annotations

# gpu1 实测(squeue/节点信息): 64 核 / 124703 MB
NODE_CORES = 64
NODE_MEM_MB = 124703
# 工程约定(**非物理判据**): 只用节点内存的 80%, 给调度与其它作业留余量。
MEM_HEADROOM = 0.8


def _gi(n_or_spec) -> float:
    """'48G' | '48G' | 48 | '48000M' | 48000 -> GiB。"""
    if isinstance(n_or_spec, (int, float)):
        return float(n_or_spec)
    s = str(n_or_spec).strip().upper()
    num = ''.join(ch for ch in s if ch.isdigit() or ch == '.')
    if not num:
        return 0.0
    val = float(num)
    if s.endswith('TB') or s.endswith('T'):
        return val * 1024.0
    if s.endswith('MB') or s.endswith('M'):   # 裸 'M' 亦按 MiB
        return val / 1024.0
    return val          # 'G'/'GB'/裸数字 -> GiB(本模块约定)


def safe_concurrency(mem_per_job, cores_per_job: int = 1, node_mem_mb: int = NODE_MEM_MB,
                     node_cores: int = NODE_CORES, headroom: float = MEM_HEADROOM) -> dict:
    """给定单作业资源, 返回**内存与核数双约束**下的并发水位(取小者)。

    本会话口径: aimu_scf 单作业 48G / 4 核 -> 内存侧紧, 核数侧宽 -> 水位由内存决定。
    """
    mem_gib = _gi(mem_per_job)
    if mem_gib <= 0:
        return {'safe': 0, 'bound_by': 'invalid', 'reason': '单作业内存未显式给出(违反 §7.1)'}
    usable = (node_mem_mb / 1024.0) * headroom
    by_mem = int(usable // mem_gib)
    by_cpu = int(node_cores // max(int(cores_per_job), 1))
    safe = max(0, min(by_mem, by_cpu))
    return {'safe': safe, 'by_mem': by_mem, 'by_cpu': by_cpu,
            'bound_by': 'mem' if by_mem <= by_cpu else 'cpu',
            'usable_mem_gib': round(usable, 2), 'mem_per_job_gib': mem_gib,
            'reason': '水位 %d(内存侧 %d / 核数侧 %d, 受限于 %s)'
                      % (safe, by_mem, by_cpu, 'mem' if by_mem <= by_cpu else 'cpu')}


def check_memory_discipline(mem_spec, node_mem_mb: int = NODE_MEM_MB,
                            headroom: float = MEM_HEADROOM) -> dict:
    """内存申请纪律: 必须显式 `--mem`, 且不得做整节点默认申请(任务卡 §7.1)。"""
    if mem_spec in (None, '', 0):
        return {'ok': False, 'reason': '缺显式 --mem(违反任务卡 §7.1: 不做整节点默认申请)'}
    mem_gib = _gi(mem_spec)
    limit = (node_mem_mb / 1024.0) * headroom
    if mem_gib > limit:
        return {'ok': False, 'mem_gib': mem_gib, 'limit_gib': round(limit, 2),
                'reason': ('单作业申请 %.1f GiB 超节点可用上限 %.1f GiB(124703MB x %.0f%%) '
                           '—— 属整节点默认申请' % (mem_gib, limit, headroom * 100))}
    return {'ok': True, 'mem_gib': mem_gib, 'limit_gib': round(limit, 2),
            'reason': '显式 --mem %.1f GiB 在纪律范围内' % mem_gib}


# ---------------------------------------------------------------------------
# W1-3 (2026-11-05, AGD-2026-11-03-SMJ-HTP): 多阶段二维水位 —— 从"单作业水位"升级到
# "一个结要占多少槽位"。原先 safe_concurrency 只回答"一个作业能开几个"(全包无调用者),
# 但一批结的吞吐不是单作业水位的函数: 一个结要**串行**跑 SCF 与 TE, 两者资源画像差 4 倍,
# 于是"SCF 开满"与"TE 开满"会互相抢内存/核, 真瓶颈是 min(n_scf/t_scf, n_te/t_te)。
# 例(gpu1, 132 原子): SCF 12G/4 核/1.5 h, TE 11G/16 核/0.9 h —— 全给 SCF 得 8 槽,
# 但 TE 抢不到槽, 吞吐被 TE 卡死; 最优是 5 SCF + 2 TE ≈ 80 结/日(实测口径)。
# ---------------------------------------------------------------------------
def plan_slots(specs: dict, node_cores: int = NODE_CORES, node_mem_mb: int = NODE_MEM_MB,
               headroom: float = MEM_HEADROOM, derate: float = 1.0) -> dict:
    """给定一个结的**各阶段资源画像**, 求节点上的最优并发分配与日吞吐。

    specs: {阶段: {'mem': '12G', 'cores': 4, 'wall_s': 5400}}(wall_s 为秒; mem/cores 必给)
    derate: 有效利用率折扣(留调度空转/排队/人工介入), 0–1。

    返回 {'slots': {阶段: n}, 'junctions_per_day': x, 'bound_by': 'mem'|'cpu'|'mix',
          'usable_mem_gib':, 'cores_used':, 'mem_used_gib':, 'reason':}
    口径: 一个结 = 每个阶段各跑 1 次(串行) -> 日吞吐 = 24*3600*derate / max(墙钟/槽位)。
    不可行(某阶段连 1 槽都放不下)时 slots 该阶段为 0、junctions_per_day=0, 原因写明。
    """
    if not specs:
        return {'slots': {}, 'junctions_per_day': 0.0, 'bound_by': 'invalid',
                'reason': '未给阶段画像(不猜)'}
    usable = (node_mem_mb / 1024.0) * headroom
    stages, caps = [], []
    for name, spec in specs.items():
        mem_gib = _gi(spec.get('mem'))
        cores = max(int(spec.get('cores') or 1), 1)
        wall_s = float(spec.get('wall_s') or 0)
        if mem_gib <= 0 or wall_s <= 0:
            return {'slots': {name: 0}, 'junctions_per_day': 0.0, 'bound_by': 'invalid',
                    'reason': '%s 缺 mem 或 wall_s(不猜)' % name}
        stages.append((name, mem_gib, cores, wall_s))
        caps.append(max(0, min(int(usable // mem_gib), node_cores // cores)))

    # 每个阶段各自最多能开几槽; 再暴力枚举组合(阶段数 ≤4、单阶段槽位 ≤64, 组合数可控)
    best = None
    def _rec(i, chosen):
        nonlocal best
        if i == len(stages):
            cores_used = sum(c * chosen[k] for k, (_, _m, c, _w) in enumerate(stages))
            mem_used = sum(chosen[k] * stages[k][1] for k in range(len(stages)))
            if cores_used > node_cores or mem_used > usable + 1e-9:
                return
            if any(n <= 0 for n in chosen):
                return
            cycle_s = max(w / n for (_nm, _m, _c, w), n in zip(stages, chosen))
            jpd = 24 * 3600 * derate / cycle_s
            key = (round(jpd, 4), -mem_used, -cores_used)
            if best is None or key > best[0]:
                best = (key, list(chosen), cores_used, mem_used, cycle_s)
            return
        for n in range(caps[i] + 1):
            _rec(i + 1, chosen + [n])

    _rec(0, [])
    if best is None:
        return {'slots': {nm: 0 for nm, _m, _c, _w in stages}, 'junctions_per_day': 0.0,
                'bound_by': 'infeasible', 'usable_mem_gib': round(usable, 2),
                'reason': '任一阶段都放不下 1 槽: 单结需求 %s'
                          % {nm: '%s/%d 核' % (specs[nm].get('mem'), specs[nm]['cores'])
                             for nm, _m, _c, _w in stages}}
    _key, chosen, cores_used, mem_used, cycle_s = best
    slots = {nm: n for (nm, _m, _c, _w), n in zip(stages, chosen)}
    core_room = node_cores - cores_used
    mem_room = usable - mem_used
    if mem_room < min(m for _nm, m, _c, _w in stages) and core_room < min(
            c for _nm, _m, c, _w in stages):
        bound_by = 'mix'
    elif mem_room < min(m for _nm, m, _c, _w in stages):
        bound_by = 'mem'
    else:
        bound_by = 'cpu'
    return {'slots': slots, 'junctions_per_day': round(24 * 3600 * derate / cycle_s, 1),
            'cycle_hours': round(cycle_s / 3600.0, 2),
            'bound_by': bound_by, 'usable_mem_gib': round(usable, 2),
            'cores_used': cores_used, 'mem_used_gib': round(mem_used, 2),
            'reason': ('最优分配 %s: 单结周期 %.2f h(瓶颈阶段取最慢), 余量 %.1f GiB / %d 核; '
                       '受限 %s' % (slots, cycle_s / 3600.0, mem_room, core_room, bound_by))}
