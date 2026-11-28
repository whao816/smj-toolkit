"""哨兵式批量等待: 一条 ssh 覆盖 N 个作业的"完成 + 队列"判定。

## 为什么需要它

`backends/slurm_ssh.py` 的 `wait_many` 已对 **vasp** 家族做了批量探针
(`tools/submit.py::_batch_probe`, 读 OUTCAR 完成标记) —— 但它的批量分支硬性要求
`family == 'vasp'`, 其余家族整批回退逐句柄 `_wait_remote`, 也就是**每作业每 tick 一次 ssh**。

SMJ 线的四个真实阶段(Gaussian / aims_tip_relax / aims_scf / aitranss_te)都没有 OUTCAR,
它们的完成标记是阶段脚本自己写的 `.dsh_done`(内容 `<退出码> <ISO8601 时间>`)。于是每个
单结脚本都得自己写一遍轮询, 而且写错一次的代价很实在 —— 2026-11-03 第一轮 σ 扫描的 ssh
瞬断让 `squeue` 返回空 stdout, 三个仍在跑的作业被同时判成"已完成"。

本模块把这件事收敛成一份实现: **一个轮询 tick = 一次 rc.run**, 内里既读哨兵又查队列,
并给三种"不知道"留出安全出口(见下"语义")。campaign runner(W4)与单结实验脚本共用它。

## 语义(与 `_batch_probe` 对齐, 绝不因一次查询失败误杀)

- **finished**: 该目录存在 `.dsh_done`。内容首列是退出码; 非 0 -> `failed=True` 且
  `note` 带上退出码(作业跑完了但失败了, 与"没跑起来"不是一回事)。
- **queue unknown**: `squeue` 退出码非 0(含 ssh 瞬断、`Invalid job id` 的 stderr 场景)
  -> 本轮**不累计缺席**, 只继续等。把"查不到"当"已离队"正是上面那次误判的根因。
- **absent 计数**: 队列查询成功、该 job_id 不在队列、且无哨兵 -> 记一次; **连续两次**
  才判 `failed`(作业消失且没有任何完成标记)。中间任何一次查到它在队列里 -> 计数清零。
- **timeout**: 到点返回 `finished=False` 且**不置 failed** —— 作业可能还在跑, 由调用方
  决定续等还是取消(单结脚本的经验: 驱动自设的 8 h 上限曾在 9.6 h 的 SCF 上误杀等待)。

只读远端、不改任何远端状态; 不 sleep 在 import 时, 时钟与 sleep 可注入以便测试。
"""

from __future__ import annotations

import shlex
import time
from typing import Callable, Iterable, Mapping

__all__ = ['SENTINEL', 'build_probe_command', 'parse_probe', 'wait_many']

#: 阶段脚本写在 stage 目录里的完成哨兵(见 `templates/sub_aims.sh` / `sub_te.sh` / 驱动脚本)。
SENTINEL = '.dsh_done'

_MARK_SENT = '@@SENT'
_MARK_QUEUE = '@@Q'
_MARK_QUEUE_RC = '@@QRC'
_MARK_ID = '@@ID'
_MARK_END = '@@END'


def build_probe_command(dirs: Iterable[str], job_ids: Iterable[str]) -> str:
    """拼一条 shell 命令: 逐个目录读哨兵, 最后一次 squeue 查全部 job_id。

    队列查询的退出码单独回传(`@@QRC`), 解析侧据此区分"没有作业在排队"与"这次没查成"。
    """
    dirs = [str(d) for d in dirs]
    ids = [str(j) for j in job_ids]
    parts = ['set +e']
    for i, d in enumerate(dirs):
        q = shlex.quote(d)
        # `$(cat ...)` 缺失时展开为空串 -> 该行尾部无内容, 解析成"无哨兵"。
        parts.append('printf "%s %d %%s\\n" "$(cat %s/%s 2>/dev/null | tr -d "\\r\\n")"'
                     % (_MARK_SENT, i, q, SENTINEL))
    if ids:
        parts.append('_adft_q=$(squeue -h -j %s -o "%%i" 2>/dev/null); _adft_rc=$?'
                     % shlex.quote(','.join(ids)))
        parts.append('printf "%s %%s\\n" "$_adft_rc"' % _MARK_QUEUE_RC)
        # 必须 `printf "%s\n" "$_adft_q"` —— 写成 `printf "@@Q\n" "$_adft_q"` 时格式串没有
        # 占位符, printf 会**丢掉**参数, 于是队列信息整段消失、解析侧把每个作业都当成
        # "不在队列"(首版就是这么误判的: 7 个正在跑的作业在 1.5 分钟内全被判失败)。
        parts.append('printf "%%s\\n" "$_adft_q" | sed "s/^/%s /"' % _MARK_ID)
    parts.append('printf "%s\\n"' % _MARK_END)
    return '; '.join(parts)


def parse_probe(out: str) -> dict:
    """解析 `build_probe_command` 的输出。

    返回 `{'sentinels': {idx: '内容'}, 'queue_known': bool, 'queued': set[str], 'complete': bool}`。
    `queue_known=False` 表示这次队列查询不可用(非零退出码或该行缺失), 调用方**不得**据此累计缺席。
    判定的唯一依据是 `@@QRC` 行本身: 它是本模块的命令才会打印的, 不必再要求 `@@Q` 段头出现
    (空队列时 `@@Q` 后面一行都没有, 早先按段头判会把"队列明确的空"误当成"没查成")。
    """
    sentinels: dict[int, str] = {}
    queued: set[str] = set()
    queue_known = False
    complete = False
    for line in (out or '').splitlines():
        line = line.rstrip('\r')
        if line.startswith(_MARK_SENT + ' '):
            rest = line[len(_MARK_SENT) + 1:]
            idx_s, _, payload = rest.partition(' ')
            try:
                sentinels[int(idx_s)] = payload.strip()
            except ValueError:
                continue
        elif line.startswith(_MARK_QUEUE_RC + ' '):
            queue_known = line[len(_MARK_QUEUE_RC) + 1:].strip() == '0'
        elif line.startswith(_MARK_ID + ' '):
            jid = line[len(_MARK_ID) + 1:].strip()
            if jid:
                queued.add(jid)
                # squeue %i lists array elements as 712_3 or 712_[0-11%2].
                # sbatch returns their parent 712, which stays active for both.
                parent, separator, _ = jid.partition('_')
                if separator and parent.isdigit():
                    queued.add(parent)
        elif line.startswith(_MARK_END):
            complete = True
    return {'sentinels': sentinels, 'queue_known': queue_known,
            'queued': queued, 'complete': complete}


def _sentinel_exit_code(payload: str) -> int:
    """`.dsh_done` 内容首列是退出码; 解析不出时按 -1(未知失败)处理, 不当作成功。"""
    head = (payload or '').split()[0] if (payload or '').split() else ''
    try:
        return int(head)
    except (TypeError, ValueError):
        return -1


def wait_many(rc, jobs: Iterable[Mapping], *, poll_s: float = 60.0,
              timeout_s: float = 6 * 3600.0,
              sleep: Callable[[float], None] = time.sleep,
              clock: Callable[[], float] = time.time,
              on_tick: Callable[[int, int], None] | None = None) -> list[dict]:
    """批量等待哨兵式作业(见模块 docstring 的四条语义)。

    `jobs`: 每项需有 `job_id` 与 `remote_dir`(其它键原样忽略)。返回与输入等长的列表,
    元素形如 `{'job_id', 'finished', 'failed', 'exit_code', 'note', 'elapsed_s', 'polls'}`。
    """
    jobs = [dict(j) for j in jobs]
    if not jobs:
        return []
    dirs = []
    for j in jobs:
        if not str(j.get('job_id') or '').strip():
            raise ValueError('wait_many: 作业缺少 job_id: %r' % (j.get('remote_dir'),))
        if not str(j.get('remote_dir') or '').strip():
            raise ValueError('wait_many: 作业缺少 remote_dir: %r' % (j.get('job_id'),))
        dirs.append(str(j['remote_dir']))
    ids = [str(j['job_id']) for j in jobs]

    poll_s = float(poll_s)
    timeout_s = float(timeout_s)
    start = clock()
    results: list[dict | None] = [None] * len(jobs)
    absent = [0] * len(jobs)
    ticks = 0
    while True:
        pending = [i for i, r in enumerate(results) if r is None]
        if not pending:
            break
        cmd = build_probe_command([dirs[i] for i in pending], [ids[i] for i in pending])
        try:
            code, out, _err = rc.run(cmd, timeout=120)
        except Exception as exc:                      # ssh 抛异常 = 未知, 继续等
            code, out = 1, '%s: %s' % (type(exc).__name__, exc)
        probe = parse_probe(out) if code == 0 else {'sentinels': {}, 'queue_known': False,
                                                    'queued': set(), 'complete': False}
        elapsed = round(clock() - start, 1)
        ticks += 1
        for pos, i in enumerate(pending):
            payload = probe['sentinels'].get(pos)
            if payload:
                code_i = _sentinel_exit_code(payload)
                results[i] = {'job_id': ids[i], 'finished': True,
                              'failed': code_i != 0, 'exit_code': code_i,
                              'note': '' if code_i == 0 else '哨兵退出码 %d' % code_i,
                              'elapsed_s': elapsed, 'polls': ticks}
                continue
            if probe['queue_known'] and ids[i] not in probe['queued']:
                absent[i] += 1
                if absent[i] >= 2:
                    results[i] = {'job_id': ids[i], 'finished': False, 'failed': True,
                                  'exit_code': None,
                                  'note': ('作业 %s 连续 2 轮不在队列且无 %s'
                                           '(疑似未跑起来/被取消; 查 slurm-*.out 与阶段日志)'
                                           % (ids[i], SENTINEL)),
                                  'elapsed_s': elapsed, 'polls': ticks}
            elif probe['queue_known']:
                absent[i] = 0
            # queue_known=False: 本轮不知道, 缺席计数原地不动
        if on_tick is not None:
            on_tick(ticks, len([r for r in results if r is None]))
        if all(r is not None for r in results):
            break
        if clock() - start >= timeout_s:
            break
        sleep(poll_s)

    elapsed = round(clock() - start, 1)
    for i, r in enumerate(results):
        if r is None:
            results[i] = {'job_id': ids[i], 'finished': False, 'failed': False,
                          'exit_code': None, 'note': '超时(%ds)未完成' % int(timeout_s),
                          'elapsed_s': elapsed, 'polls': ticks}
    return results  # type: ignore[return-value]
