"""P0-6 (AGD-2026-09-25-SMJ): 建结几何指纹与 A2 人工确认绑定。

问题(实测): 本会话同一表面索引算法在三组几何上给出三组左表面索引
([55,56,60] / [58,52,57] / [50,51,56]), 而 A2(表面原子人工确认)只在**占位几何**上
做过一次 —— 几何一变, 旧确认**静默失效**却仍在被引用(v0.4 D12; 三组索引见
docs/smj-bdt-benchmark-results.md §9.4)。

做法: au60_build 记录 surface_fingerprint(geometry.in 的字节 sha256);
A2 确认落成 a2_confirmation.json{fingerprint, indices, who, when};
校验时**指纹不符即确认失效**, 必须按新几何重新确认。

纯文件操作, 不做真实计算。
"""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

FINGERPRINT_KEY = 'surface_fingerprint'
CONFIRM_FILE = 'a2_confirmation.json'
SURFACE_KEYS = ('left_surface_indices', 'right_surface_indices')


def surface_fingerprint(geometry_path) -> str:
    """geometry.in 的字节级 sha256 —— 人工确认的就是这份文件本身。"""
    return hashlib.sha256(Path(geometry_path).read_bytes()).hexdigest()


def attach_fingerprint(step_dir, metadata_path=None) -> dict:
    """把 surface_fingerprint 与汇总的 surface_indices 写进 electrode_metadata.json。"""
    step_dir = Path(step_dir)
    mp = Path(metadata_path) if metadata_path else step_dir / 'electrode_metadata.json'
    meta = json.loads(mp.read_text(encoding='utf-8'))
    meta[FINGERPRINT_KEY] = surface_fingerprint(step_dir / 'geometry.in')
    meta['surface_indices'] = {k: meta.get(k) for k in SURFACE_KEYS}
    mp.write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding='utf-8')
    return meta


def record_confirmation(step_dir, who: str, when: str | None = None,
                        metadata_path=None) -> dict:
    """登记 A2 人工确认(绑定当前几何指纹)。"""
    step_dir = Path(step_dir)
    mp = Path(metadata_path) if metadata_path else step_dir / 'electrode_metadata.json'
    meta = json.loads(mp.read_text(encoding='utf-8'))
    fresh = surface_fingerprint(step_dir / 'geometry.in')   # 现算(不用建结时快照)
    meta[FINGERPRINT_KEY] = fresh                            # 同步快照, 保持 metadata 真实
    meta['surface_indices'] = {k: meta.get(k) for k in SURFACE_KEYS}
    mp.write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding='utf-8')
    rec = {'fingerprint': fresh,
           'indices': {k: meta.get(k) for k in SURFACE_KEYS},
           'who': who,
           'when': when or datetime.now(timezone.utc).isoformat(timespec='seconds')}
    (step_dir / CONFIRM_FILE).write_text(json.dumps(rec, ensure_ascii=False, indent=2),
                                         encoding='utf-8')
    return rec


def check_confirmation(step_dir) -> dict:
    """A2 确认有效性: {'present','valid','reason',...}。

    无记录 -> present=False, valid=False(调用方按 warn 处理);
    指纹不符 -> present=True, valid=False(**须重新确认**, 调用方按 error 处理)。
    """
    step_dir = Path(step_dir)
    cp = step_dir / CONFIRM_FILE
    meta = {}
    mp = step_dir / 'electrode_metadata.json'
    if mp.exists():
        try:
            meta = json.loads(mp.read_text(encoding='utf-8'))
        except ValueError:
            meta = {}
    # 关键: 以**当前 geometry.in 现算**的指纹为准(metadata 里的只是建结时快照)。
    # 若优先用快照, 几何改变后仍会判"有效" —— 恰与本机制目的相反(测试已抓到该缺陷)。
    gp = step_dir / 'geometry.in'
    if gp.exists():
        cur = surface_fingerprint(gp)
    else:
        cur = meta.get(FINGERPRINT_KEY)
        if not cur:
            return {'present': False, 'valid': False, 'current': None,
                    'reason': 'geometry.in 缺失, 无法计算几何指纹'}
    if not cp.exists():
        return {'present': False, 'valid': False, 'current': cur,
                'reason': '尚无 A2 确认记录(几何指纹 %s); 需人工确认表面原子' % cur[:12]}
    try:
        rec = json.loads(cp.read_text(encoding='utf-8'))
    except (OSError, ValueError) as exc:
        return {'present': True, 'valid': False, 'current': cur,
                'reason': 'A2 确认记录不可解析: %s' % exc}
    if rec.get('fingerprint') != cur:
        return {'present': True, 'valid': False, 'current': cur,
                'recorded': rec.get('fingerprint'),
                'reason': ('A2 确认已失效: 记录指纹 %s != 当前几何指纹 %s —— 几何已变, '
                           '旧确认不得沿用, 须按新几何重新确认'
                           % (str(rec.get('fingerprint'))[:12], cur[:12]))}
    return {'present': True, 'valid': True, 'current': cur,
            'recorded': rec.get('fingerprint'), 'who': rec.get('who'), 'when': rec.get('when'),
            'reason': 'A2 确认有效(指纹 %s, %s @ %s)'
                      % (cur[:12], rec.get('who'), rec.get('when'))}
