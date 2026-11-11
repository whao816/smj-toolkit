"""结构人工确认检查点(2026-10-24, 用户 m02353): 两处 xyz 输出 + 确认后才继续算。

动机(用户指令): 「增加人工检查, 一个是分子左右各加一个 Au 原子后, 一个是贴装电极后,
要输出 xyz 文件, 供我检查结构, 结构确认没问题后再继续计算」。

为什么值得停: smj 链只有两处几何定稿点 —— 一旦相对位置/取向/键长错, 后面 aims_scf
(实测单次 3886 s) + aitranss 整轮机时全废, 而且错在几何上时输运判据只会给出"物理不达标"
的误导结论。此前两处都只是流水线里的中间文件(aims_tip_relax/geometry.in、
<auN>_build/geometry.in), 人要看得自己翻 work 目录, 也没有"确认"这一步 —— 引擎一口气跑到
提交。

做法(沿用 fingerprint.py 的教训, **不做**"确认快照"):
- 每个检查点由引擎写一份**人读 xyz**(review_tips.xyz / review_junction.xyz; 标准 XYZ 格式,
  常见可视化软件可直接打开);
- 人工确认时把**写出时这份 xyz 的字节 sha256** 落进 structure_review.json;
- 续跑/校验时**现算**当前 xyz 指纹再比对: 不符 -> 确认失效(几何变了/阶段重跑了), 必须
  重新看。绑定的是"人看到的这份文件", 不是外部状态。

纯文件操作, 不做真实计算(检查点本身不提交任何作业)。
"""
from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

from .geometry import read_aims_geometry

REVIEW_FILE = 'structure_review.json'
GEOMETRY_FILE = 'geometry.in'
# 引擎在哪些阶段停下等人看结构 -> 检查点 key(顺序即链上顺序)。
# 电极阶段名随电极规模变(au60_build / au10_build ...), 故除字面名外还按 `^au\d+_build$` 归一,
# 与 manifest.discover_smj_stages 的阶段识别同口径 —— 否则小电极变体(2CN 的 au10_build)会漏掉检查点。
REVIEW_STAGES = {'aims_tip_relax': 'tips', 'au60_build': 'electrodes'}
_ELECTRODE_STAGE_RE = re.compile(r'^au\d+_build$')
REVIEW_XYZ = {'tips': 'review_tips.xyz', 'electrodes': 'review_junction.xyz'}
REVIEW_TITLES = {
    'tips': '分子 + 左右各 1 个 Au tip(FHI-aims 弛豫前)',
    'electrodes': '贴装左右电极后的结区(进入 SCF/输运前)',
}
# 引擎写进 rec.results['error'] 的前缀 —— 确认通过后据此清掉这条"待确认"叙述,
# 免得成功续跑的 run 报告里还挂着旧提示。
ERROR_PREFIX = '结构检查点未确认: '


def review_key(stage) -> str | None:
    r"""阶段名 -> 检查点 key(tips/electrodes); 不是检查点阶段 -> None。

    `au\d+_build`(au60_build / au10_build / ...)一律归一到 electrodes。
    """
    s = str(stage or '')
    if s in REVIEW_STAGES:
        return REVIEW_STAGES[s]
    return 'electrodes' if _ELECTRODE_STAGE_RE.match(s) else None


def xyz_path(step_dir, key) -> Path:
    return Path(step_dir) / REVIEW_XYZ[key]


def confirmation_path(step_dir) -> Path:
    return Path(step_dir) / REVIEW_FILE


def discover_review_stages(work_root) -> list:
    """按链上顺序列出 work_root 下**实际存在**的检查点阶段目录(tips 先, 电极后)。"""
    root = Path(work_root)
    out = []
    if (root / 'aims_tip_relax').is_dir():
        out.append('aims_tip_relax')
    if root.is_dir():
        out.extend(sorted(p.name for p in root.iterdir()
                          if p.is_dir() and _ELECTRODE_STAGE_RE.match(p.name)))
    return out


def stage_for_key(key, work_root=None) -> str:
    """检查点 key -> 实际阶段目录名(电极阶段随规模可能是 au60_build / au10_build)。"""
    if key not in REVIEW_XYZ:
        raise KeyError(key)
    if work_root is not None:
        for s in discover_review_stages(work_root):
            if review_key(s) == key:
                return s
    return next((s for s, k in REVIEW_STAGES.items() if k == key), key)


def _sha256(text: str) -> str:
    return hashlib.sha256(text.encode('utf-8')).hexdigest()


def render_xyz(step_dir, key, stage=None, system=None) -> dict:
    """由当前 geometry.in 渲染检查点 xyz 文本(不落盘)。

    几何文件读不到/为空 -> ValueError(调用方按失败处理, 不能"没几何也放行")。
    """
    step_dir = Path(step_dir)
    gp = step_dir / GEOMETRY_FILE
    try:
        atoms = read_aims_geometry(gp)              # 空/无 atom 行 -> ValueError
        geom_bytes = gp.read_bytes()
    except OSError as exc:                          # 文件缺失/不可读 -> 同一契约(ValueError)
        raise ValueError('无法读取检查点几何 %s: %s' % (gp, exc))
    geom_fp = hashlib.sha256(geom_bytes).hexdigest()
    stage = stage or next((s for s, k in REVIEW_STAGES.items() if k == key), key)
    comment = ('smj 检查点 %s | %s | 体系 %s | 原子 %d | 来源 %s | geometry sha256 %s'
               % (key, REVIEW_TITLES.get(key, stage), system or '未知', len(atoms),
                  GEOMETRY_FILE, geom_fp[:12]))
    lines = [str(len(atoms)), comment]
    lines.extend('%2s %16.8f %16.8f %16.8f' % (a.element, a.x, a.y, a.z) for a in atoms)
    text = '\n'.join(lines) + '\n'
    return {'text': text, 'atoms': len(atoms), 'comment': comment,
            'fingerprint': _sha256(text), 'geometry_fingerprint': geom_fp,
            'geometry_path': str(gp), 'stage': stage, 'key': key}


def write_review_xyz(step_dir, key, stage=None, system=None) -> dict:
    """渲染并写出人读 xyz; 返回 {'path','atoms','fingerprint','comment',...}。

    fingerprint = **写出文件的字节 sha256**(即人看到的这份内容), 与 confirmation 绑定。
    """
    step_dir = Path(step_dir)
    info = render_xyz(step_dir, key, stage=stage, system=system)
    p = xyz_path(step_dir, key)
    p.write_text(info['text'], encoding='utf-8', newline='\n')
    out = {k: info[k] for k in ('atoms', 'comment', 'fingerprint', 'geometry_fingerprint',
                                'geometry_path', 'stage', 'key')}
    out['path'] = str(p)
    return out


def check_confirmation(step_dir, key, system=None) -> dict:
    """检查点确认有效性: {'present','valid','reason',...}。

    无记录 -> present=False, valid=False(调用方按"需确认"处理);
    xyz 指纹不符 -> present=True, valid=False(**必须重新确认**, 几何已变)。
    与 fingerprint.check_confirmation 同原则: 以**现算**为准, 不用记录里的快照。

    **system 必须与写入时一致**: 人读 xyz 的注释行带体系名, 指纹是整份文件字节的 sha256
    —— 传 None 会渲染成"未知", 与 record_confirmation/evaluate 写入时的体系名不符, 于是
    同一几何也会被判"已失效"(2026-10-24 真实 run 的 tips 检查点因此永远无法确认)。
    """
    step_dir = Path(step_dir)
    cp = confirmation_path(step_dir)
    if not cp.exists():
        return {'present': False, 'valid': False, 'current': None,
                'reason': '尚无结构确认记录(%s 不存在)' % REVIEW_FILE}
    try:
        data = json.loads(cp.read_text(encoding='utf-8'))
    except (OSError, ValueError) as exc:
        return {'present': True, 'valid': False, 'current': None,
                'reason': '结构确认记录不可解析: %s' % exc}
    ent = ((data.get('reviews') if isinstance(data, dict) else None) or {}).get(key)
    if not isinstance(ent, dict):
        return {'present': False, 'valid': False, 'current': None,
                'reason': '确认记录里没有检查点 %s 的条目' % key}
    try:
        cur = render_xyz(step_dir, key, system=system)['fingerprint']
    except (ValueError, OSError) as exc:
        return {'present': True, 'valid': False, 'current': None,
                'reason': '无法据当前 %s 重算 xyz 指纹: %s' % (GEOMETRY_FILE, exc)}
    rec_fp = str(ent.get('fingerprint') or '')
    if rec_fp != cur:
        return {'present': True, 'valid': False, 'current': cur, 'recorded': rec_fp,
                'atoms': ent.get('atoms'), 'who': ent.get('who'), 'when': ent.get('when'),
                'reason': ('结构确认已失效: 记录指纹 %s != 当前结构指纹 %s —— 几何已变'
                           '(或该阶段已重跑), 旧确认不得沿用, 须重新检查 xyz 并确认'
                           % (rec_fp[:12] or '(空)', cur[:12]))}
    return {'present': True, 'valid': True, 'current': cur, 'recorded': rec_fp,
            'atoms': ent.get('atoms'), 'who': ent.get('who'), 'when': ent.get('when'),
            'comment': ent.get('comment') or '',
            'reason': ('结构确认有效(指纹 %s, %s @ %s)'
                       % (cur[:12], ent.get('who'), ent.get('when')))}


def record_confirmation(step_dir, key, who: str, when: str | None = None,
                        comment: str = '', stage=None, system=None) -> dict:
    """登记人工结构确认(绑定**当前**写出 xyz 的指纹)。

    先重新渲染/写出 xyz 再取指纹: 保证记录的就是确认时人看到的那份内容。
    """
    step_dir = Path(step_dir)
    info = write_review_xyz(step_dir, key, stage=stage, system=system)
    cp = confirmation_path(step_dir)
    data = {}
    if cp.exists():
        try:
            data = json.loads(cp.read_text(encoding='utf-8'))
        except ValueError:
            data = {}                      # 坏文件不足以拒绝确认: 覆盖重写(下面整体落盘)
        if not isinstance(data, dict):
            data = {}
    reviews = data.get('reviews')
    if not isinstance(reviews, dict):
        reviews = {}
    when = when or datetime.now(timezone.utc).isoformat(timespec='seconds')
    reviews[key] = {'stage': info['stage'], 'xyz': REVIEW_XYZ[key],
                    'fingerprint': info['fingerprint'], 'atoms': info['atoms'],
                    'geometry_fingerprint': info['geometry_fingerprint'],
                    'who': who, 'when': when, 'comment': comment}
    data['reviews'] = reviews
    data['updated_at'] = when
    cp.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding='utf-8')
    return {'key': key, 'stage': info['stage'], 'xyz': info['path'],
            'fingerprint': info['fingerprint'], 'atoms': info['atoms'],
            'who': who, 'when': when, 'comment': comment}


def evaluate(step_dir, stage, system=None) -> dict | None:
    """检查点当前状态: None(非检查点阶段) 或 {'key','stage','xyz','atoms','fingerprint',...}。

    总是重写一次 xyz(幂等): 保证磁盘上的人读文件与当前几何一致, 且与现实算指纹同源。
    """
    key = review_key(stage)
    if key is None:
        return None
    step_dir = Path(step_dir)
    info = write_review_xyz(step_dir, key, stage=stage, system=system)
    # 写入与校验必须同源(system 参与注释行 -> 参与指纹): 只给写不给查会把同一几何判成
    # "已失效", 真实 run 的检查点永远等不到确认(2026-10-24 缺口⑥)。
    chk = check_confirmation(step_dir, key, system=system)
    return {'key': key, 'stage': stage, 'xyz': info['path'], 'atoms': info['atoms'],
            'fingerprint': info['fingerprint'], 'comment': info['comment'],
            'approved': bool(chk['valid']), 'who': chk.get('who'), 'when': chk.get('when'),
            'reason': chk.get('reason') or ''}


def pending_message(run_id: str, info: dict) -> str:
    """停下来的可操作提示(写进 rec.results['error'], CLI 原样打印)。

    2026-11-11 本地补丁 F4: 原提示指向主仓的 `agenticdft smj-review <run> --approve`,
    独立包 smj-toolkit 里没有这个命令。且独立包**没有**接线检查点 xyz
    (structure_review.json)的 CLI —— 唯一接线的类似门是 auN_build 的 A2 表面
    三元组确认, 由 `smj review --approve --who <名字>` 登记(fingerprint 机制,
    审计/出图不通过时拒绝签署, 退出码 2)。这里如实说明, 不再给不存在的命令。
    """
    return (ERROR_PREFIX
            + '阶段 %s 已写出待检结构 %s(原子 %d, 指纹 %s)。请先人工检查该 xyz; '
              '独立包不提供检查点 xyz 的签署命令, 与它对应的真实门是 A2 表面确认: '
              '`smj review --approve --who <你的名字>`(审计/出图不通过会拒绝, 退出码 2)。'
              '若结构有问题, 调整旋钮(左右锚定原子/Au tip '
              '距离/每侧电极原子数等)后按原任务重跑, 不要带着错几何进集群。'
              % (info['stage'], info['xyz'], info['atoms'], info['fingerprint'][:12]));