"""P2-4 (AGD-2026-09-25-SMJ): run manifest —— 可追溯的计算过程单一出口。

动机: 本会话那份「§8 交付物」是我**手工拼**的(收敛表/结果表/对比表/图), 溯源信息则散落在
驱动日志、results.db、任务卡与报告里。本模块把一条 run 的可追溯信息收成单一 manifest.json。

设计上复用既有判据(不另起炉灶):
- 是否对某文件算 sha256, 由 P1-7 的 `artifacts.retrieval_mode` 决定 —— 标记为 head/skip 的
  大文件(如 mos.aims 0.96 GB / omat.aims 0.48 GB)只记 size, **不做哈希**(避免为 2.2 GB
  产物算校验而白耗时); must/full 的文件才哈希。
- 判定直接采纳 P1-1/P1-2 的 `validate_junction_run`(含 physics_verdict 与逐条判据)。

纯本地聚合, 不触网。
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
from pathlib import Path

# 已知 smj 阶段(与 STEPS 一致; 用于识别 run 目录下的阶段子目录)
SMJ_STAGES = ('gaussian_opt', 'aims_tip_relax', 'au60_build', 'aims_scf', 'aitranss_te')
_HASH_MAX_BYTES = 64 * 1024 * 1024      # 超过 64 MiB 且未被策略标为 full 的不哈希(工程约定)

# 建结阶段的目录名由**电极规模**决定, 不是固定字面量: 全量电极是 au60_build,
# 小电极变体(2cn 主场真实 run)是 au10_build。若照 SMJ_STAGES 字面量数数, au10_build
# 会在 manifest/report/入库里整个消失(产物看不见、阶段数 5 变 4)。故按模式发现。
_BUILD_STAGE_RE = re.compile(r'^au\d+_build$')
_STAGE_PREFIX = ('gaussian_opt', 'aims_tip_relax')
_STAGE_SUFFIX = ('aims_scf', 'aitranss_te', 'junction_analyze')


def discover_smj_stages(run_root) -> list:
    """按 run 目录里**实际存在**的阶段子目录, 按执行顺序返回阶段名。

    顺序: gaussian_opt → aims_tip_relax → <auN_build> → aims_scf → aitranss_te
    → junction_analyze。不存在的阶段不出现在返回值里(与 build_manifest 逐目录
    跳过的旧行为一致), 因此调用方拿到的是"这条 run 真跑过哪些阶段"。
    """
    root = Path(run_root)
    out = [s for s in _STAGE_PREFIX if (root / s).is_dir()]
    if root.is_dir():
        out.extend(sorted(p.name for p in root.iterdir()
                          if p.is_dir() and _BUILD_STAGE_RE.match(p.name)))
    out.extend(s for s in _STAGE_SUFFIX if (root / s).is_dir())
    return out


def sha256_file(path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def _file_entry(stage: str, f: Path) -> dict:
    e = {'name': f.name, 'bytes': f.stat().st_size}
    mode = 'full'
    try:
        from .artifacts import retrieval_mode
        mode = retrieval_mode(stage, f.name)
    except Exception:
        pass
    e['retrieval'] = mode
    if mode == 'full' or e['bytes'] <= _HASH_MAX_BYTES:
        e['sha256'] = sha256_file(f)
    else:
        e['sha256'] = None
        e['note'] = '按 P1-7 策略免哈希(retrieval=%s); 需要时单独取' % mode
    return e


def build_manifest(run_root, run_id: str | None = None, spec: dict | None = None,
                   stages=None, system=None) -> dict:
    """汇总一条 run 的可追溯信息 -> manifest dict(可 json.dump)。

    run_root: 形如 work/<material>/<run_id>/ 的目录(内含各阶段子目录)。
    system: 体系标识('BDT'/'2CN'), 传给判据适用域判定; 缺省未知 -> 全部判据适用。
    """
    root = Path(run_root)
    rid = run_id or root.name
    man = {'run_id': rid, 'root': str(root), 'spec': spec or {},
           'stages': {}, 'verdict': None, 'totals': {}}
    total = 0
    names = stages or discover_smj_stages(root)
    for st in names:
        d = root / st
        if not d.is_dir():
            continue
        files = [_file_entry(st, f) for f in sorted(d.iterdir()) if f.is_file()]
        sub = sum(f['bytes'] for f in files)
        total += sub
        man['stages'][st] = {'dir': str(d), 'n_files': len(files),
                             'bytes': sub, 'files': files}

    # 判定: 有 aitranss_te 阶段就按 P1-1/P1-2 的口径出 verdict
    te_dir = root / 'aitranss_te'
    if te_dir.is_dir():
        try:
            from .validate import validate_junction_run
            # 传 run 根: 校验层按阶段子目录解析产物(gaussian/aims/te 全链),
            # 只传 te_dir 会看不到 mol.log/aims.out, 工程判定不完整
            v = validate_junction_run(root, system=system)
            man['verdict'] = {k: v.get(k) for k in
                              ('passed', 'engineering_passed', 'errors', 'warnings',
                               'physics_verdict', 'wording', 'allowed_action',
                               'failed_criteria', 'inconclusive_criteria',
                               'not_applicable_criteria')}
            man['checks'] = v.get('checks', [])
        except Exception as exc:
            man['verdict'] = {'error': 'validate_junction_run 失败: %s' % exc}
    man['totals'] = {'bytes': total, 'gib': round(total / (1024 ** 3), 3),
                     'stages': len(man['stages'])}
    return man


def write_manifest(run_root, out_path=None, **kw) -> dict:
    man = build_manifest(run_root, **kw)
    out = Path(out_path) if out_path else (Path(run_root) / 'manifest.json')
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(man, ensure_ascii=False, indent=2)
    # 2026-11-22 复审 P2：与 state.save 同样的原子写；旧实现直接 write_text，
    # report/外部工具可能读到半截 JSON。
    fd, tmp_name = tempfile.mkstemp(
        prefix=out.name + '.', suffix='.tmp', dir=str(out.parent))
    try:
        with os.fdopen(fd, 'w', encoding='utf-8', newline='\n') as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_name, out)
    except BaseException:
        Path(tmp_name).unlink(missing_ok=True)
        raise
    man['manifest_path'] = str(out)
    return man


# ---------------------------------------------------------------------------
# P2-4 余项: 从 manifest 机械生成报告(任务卡 §8 的可自动化部分)
#
# **边界必须说清**(不假装全自动): 本函数只产出 manifest 里**确有的**内容 ——
# 逐阶段产物/校验/体积、判据表、verdict。以下内容 manifest 里没有, 故**不由本函数生成**,
# 报告里显式列出为"待人工/其它工具补充": 图(T(E)/收敛/对比)、收敛测试结论、
# 独立复核(手算 Landauer)、与文献基准的对照叙述。
# ---------------------------------------------------------------------------
_REPORT_UNMECHANIZED = (
    '图: T(E) 曲线 / 积分网格收敛 / 计算 vs 基准 (需绘图输入, 不在 manifest 内)',
    '收敛测试结论: 基组 / SCF / 网格独立性 (须由各自测试产物给出)',
    '独立复核: 脱离管线锚点手算 Landauer L0/L1/L2 (须独立计算)',
    '文献基准对照叙述与 H0 裁决论证 (须人工撰写)',
)


def render_report(man: dict, title: str | None = None) -> str:
    """manifest -> Markdown 报告(机械可生成部分)。"""
    rid = man.get('run_id', '?')
    L = []
    L.append('# %s' % (title or ('smj run 报告: %s' % rid)))
    L.append('')
    L.append('> 由 `tools/junction/manifest.py::render_report` 从 `manifest.json` 机械生成；')
    L.append('> **不含**图、收敛测试与独立复核（见文末"未机械化部分"）。')
    L.append('')
    tot = man.get('totals') or {}
    L.append('## 0. 概览')
    L.append('')
    L.append('| 项 | 值 |')
    L.append('|---|---|')
    L.append('| run_id | `%s` |' % rid)
    L.append('| 阶段数 / 文件总量 | %s / %.3f GiB |' % (tot.get('stages'), tot.get('gib', 0)))
    v = man.get('verdict') or {}
    L.append('| 工程判定 | `engineering_passed=%s` |' % v.get('engineering_passed'))
    L.append('| 科学判定 | **`%s`** (%s) |' % (v.get('physics_verdict'), v.get('wording')))
    L.append('| 允许动作 | %s |' % v.get('allowed_action'))
    if v.get('failed_criteria'):
        L.append('| 未过判据 | %s |' % ', '.join('`%s`' % x for x in v['failed_criteria']))
    if v.get('inconclusive_criteria'):
        L.append('| 存疑判据 | %s |' % ', '.join('`%s`' % x for x in v['inconclusive_criteria']))
    if v.get('not_applicable_criteria'):
        L.append('| 未适用判据 | %s |' %
                 ', '.join('`%s`' % x for x in v['not_applicable_criteria']))
    L.append('')

    L.append('## 1. 逐阶段产物（可追溯）')
    L.append('')
    L.append('| 阶段 | 文件数 | 体积 (MiB) | 关键产物 sha256(前12) |')
    L.append('|---|---|---|---|')
    for st, info in (man.get('stages') or {}).items():
        key = [f for f in info['files'] if f['name'] in
               ('TE.dat', 'geometry.in', 'mol.log', 'aims.out')]
        shas = ', '.join('%s:%s' % (f['name'], (f['sha256'] or '未哈希')[:12]) for f in key) or '—'
        L.append('| `%s` | %d | %.2f | %s |' % (st, info['n_files'],
                                               info['bytes'] / (1024 * 1024), shas))
    L.append('')

    checks = man.get('checks') or []
    if checks:
        L.append('## 2. 判据逐条（%d 条）' % len(checks))
        L.append('')
        L.append('| 判据 | 级别 | 通过 | 说明 |')
        L.append('|---|---|---|---|')
        for c in checks:
            ok_cell = ('➖' if c.get('applicable') is False
                       else ('✅' if c.get('passed') else '❌'))
            L.append('| `%s` | %s | %s | %s |'
                     % (c.get('name'), c.get('severity'), ok_cell,
                        str(c.get('detail', '')).replace('|', '/')))
        L.append('')

    L.append('## 3. 未机械化部分（manifest 内没有，须另行补）')
    L.append('')
    for x in _REPORT_UNMECHANIZED:
        L.append('- %s' % x)
    L.append('')
    return chr(10).join(L)


def write_report(run_root, out_path=None, **kw) -> dict:
    """建 manifest 并落 report.md; 返回 {'manifest_path', 'report_path', 'manifest'}。"""
    man = write_manifest(run_root, **kw)
    out = Path(out_path) if out_path else (Path(run_root) / 'report.md')
    out.write_text(render_report(man), encoding='utf-8')
    return {'manifest_path': man.get('manifest_path'), 'report_path': str(out), 'manifest': man}
