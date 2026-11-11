"""单分子结程序族完成/收敛谓词(谓词注入接缝, AGD-2026-09-22-SMJ P4)。

submit.py 的 VASP 谓词(_finished/_oszar_converged/_remote_finished/_remote_converged)
按 family='smj' 分派到本模块。

P0-3 完成判定 sentinel 化(AGD-2026-09-25-SMJ): 旧谓词只看"文件非空", 实测被误判两次
- job 24500: 失败的 FHI-aims 把动态库报错写进 aims.out(16 行, 非空) -> 被判"完成";
- job 24513: aitranss 边算边写 TE.dat, 32 行(扫到 E-EF~-4.2 eV) 时被判"完成"。
现按三级判定:
  1. 哨兵 `.dsh_done`(阶段脚本末尾写 "退出码 时间戳", 见 templates/sub_*.sh) -> 优先;
  2. 成功标记: aims.out 含 'Have a nice day'(实测: 成功 tip_relax/SCF 各 1 次,
     失败作业 0 次); mol.log 含 'Normal termination';
  3. 结构完整: TE.dat 的末行能量须扫到 tcontrol 声明的 $eend(±0.05 eV)。
无哨兵的旧 run 走 2/3 兜底(向后兼容)。

P0-13(2026-10-24, 真跑缺口 #8) 追加一条更强的约束: **工作目录里声明了阶段脚本、
且脚本正文写 DONE_SENTINEL 时, 完成/收敛只认哨兵**。原因是 smj 的 aitranss_te 与
aims_scf 共用远端目录(见 smj_remote_dir 的 R3), 目录级成功标记属于**上一阶段**, 
会让刚提交 22 s 的 aitranss 被判"完成", 回取到 0 字节 TE.dat.part(实测 run
c4-trans-smj-692b6c)。哨兵是本轮作业自己写的, 缺失即"本轮未结束"。

纯判定, 不做真实计算。
"""
from __future__ import annotations

import re
from pathlib import Path

DONE_SENTINEL = '.dsh_done'
# 阶段脚本名(tools/submit.py 的 SMJ_STAGE_SCRIPTS 值域; 这里只做"目录里声明了哪个阶段")。
SMJ_SCRIPTS = ('sub_g16.sh', 'sub_aims.sh', 'sub_te.sh')
AIMS_SUCCESS = 'Have a nice day'
SCF_CONVERGED = 'Self-consistency cycle converged'
RELAX_CONVERGED = 'Present geometry is converged'
HARTREE_EV = 27.211386245988
EEND_TOL_EV = 0.05


def read_sentinel(cwd: Path):
    """返回 (found, rc): 哨兵首字段为退出码。损坏/不可解析按 (True, None)。"""
    p = Path(cwd) / DONE_SENTINEL
    if not p.exists():
        return (False, None)
    try:
        first = p.read_text(encoding='utf-8', errors='ignore').split()[0]
        return (True, int(first))
    except (IndexError, ValueError):
        return (True, None)


def _read(path: Path) -> str:
    try:
        return Path(path).read_text(encoding='utf-8', errors='ignore')
    except OSError:
        return ''


def _te_data_rows(text: str) -> list:
    out = []
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        parts = s.split()
        if len(parts) >= 3:
            try:
                out.append((float(parts[0]), float(parts[1]), float(parts[2])))
            except ValueError:
                continue
    return out


def _te_scan_target_ev(cwd: Path):
    """由 tcontrol 的 $eend 与 TE.dat 表头的 efermi 求扫描终点(相对 E_F, eV)。"""
    tc = _read(Path(cwd) / 'tcontrol')
    m = re.search(r'^\s*\$eend\s+([-\d.eEdD+]+)', tc, re.M)
    if not m:
        return None
    try:
        eend = float(m.group(1).replace('D', 'E').replace('d', 'e'))
    except ValueError:
        return None
    hdr = _read(Path(cwd) / 'TE.dat')
    m2 = re.search(r'efermi\s*=\s*([-\d.eEdD+]+)\s*H', hdr)
    if not m2:
        return None
    try:
        efermi = float(m2.group(1).replace('D', 'E').replace('d', 'e'))
    except ValueError:
        return None
    return (eend - efermi) * HARTREE_EV


def te_reached_end(cwd: Path, rows=None) -> bool:
    """TE.dat 是否扫到 tcontrol 声明的 $eend(±EEND_TOL_EV)。

    无 tcontrol/表头时无法判终点 -> 退化为"至少 2 个数据行"(弱判据, 显式记录)。
    """
    cwd = Path(cwd)
    if rows is None:
        rows = _te_data_rows(_read(cwd / 'TE.dat'))
    if len(rows) < 2:
        return False
    target = _te_scan_target_ev(cwd)
    if target is None:
        return True
    return max(r[1] for r in rows) >= target - EEND_TOL_EV


def smj_finished(cwd: Path) -> bool:
    """smj 阶段完成判定: 哨兵优先, 否则成功标记 + 结构完整。"""
    cwd = Path(cwd)
    found, rc = read_sentinel(cwd)
    if found:
        return rc == 0

    mol_log = cwd / 'mol.log'
    if mol_log.exists() and 'Normal termination' in _read(mol_log):
        return True
    aims_out = cwd / 'aims.out'
    if aims_out.exists() and AIMS_SUCCESS in _read(aims_out):
        return True
    if (cwd / 'TE.dat').exists() and te_reached_end(cwd):
        return True
    return False


def smj_converged(cwd: Path) -> bool:
    """smj 阶段收敛判定: Gaussian 无虚频 / SCF 自洽收敛或弛豫收敛 / TE 扫到终点。"""
    cwd = Path(cwd)
    found, rc = read_sentinel(cwd)
    if found and rc != 0:
        return False

    mol_log = cwd / 'mol.log'
    if mol_log.exists():
        text = _read(mol_log)
        if 'Normal termination' in text:
            freqs = []
            for line in re.findall(r'Frequencies --\s+(.+)', text):
                freqs.extend(float(x) for x in line.split())
            # 虚频 = Frequencies 行含负值(24 章 §5: 不静默接受)
            return not any(f < 0 for f in freqs)

    aims_out = cwd / 'aims.out'
    if aims_out.exists():
        text = _read(aims_out)
        if AIMS_SUCCESS in text and (SCF_CONVERGED in text or RELAX_CONVERGED in text):
            return True

    if (cwd / 'TE.dat').exists():
        return te_reached_end(cwd)
    return False


# --------------------------------------------------------------------------
# 远端变体: 同一判定的 shell 实现(哨兵优先, 再成功标记/结构完整)
# --------------------------------------------------------------------------
def stage_script(run_dir) -> str | None:
    """run_dir 里声明的 smj 阶段脚本名(sub_g16/sub_aims/sub_te;无/多命中 -> None)。"""
    if not run_dir:
        return None
    try:
        hits = [n for n in SMJ_SCRIPTS if (Path(run_dir) / n).is_file()]
    except OSError:
        return None
    return hits[0] if len(hits) == 1 else None


def sentinel_strict_script(run_dir) -> str | None:
    """阶段脚本自己写哨兵时返回该脚本名 —— 该阶段的完成判定**只认哨兵**。

    P0-13(2026-10-24, 真跑缺口 #8): 六个 smj 阶段里有四个共用远端目录
    (aitranss_te 与 aims_scf 同目录, 见 smj_remote_dir 的 R3 说明)。成功标记回退
    ("Have a nice day"/"Normal termination"/TE.dat 终点) 读的是**目录**而不是**本轮作业**,
    于是 aitranss 提交 22 s 就被上一阶段 SCF 的 aims.out 判成"完成": 回取到 0 字节的
    TE.dat.part 与半截 te.out, 下游 junction_analyze 硬停「缺上游 TE.dat」。
    判据必须来自本轮作业自己写的那份东西 —— 阶段脚本末尾写哨兵时, 哨兵缺失只可能意味
    "还没跑完"(或作业被杀), 不得再用任何目录级成功标记回退。
    """
    sc = stage_script(run_dir)
    if not sc:
        return None
    try:
        txt = (Path(run_dir) / sc).read_text(encoding='utf-8', errors='ignore')
    except OSError:
        return None
    return sc if DONE_SENTINEL in txt else None


def _sh_finished(rc, d: str, script: str | None = None) -> bool:
    sent = repr(d + '/' + DONE_SENTINEL)
    # 哨兵存在 -> 只认退出码(首字段为 0)
    if rc.run_ok('test -s ' + sent + ' && awk \'NR==1{exit ($1==0)?0:1}\' ' + sent):
        return True
    if rc.run_ok('test -s ' + sent + ' 2>/dev/null'):
        return False
    if script is not None:
        # 本阶段脚本写哨兵 -> 哨兵没出现 = 本轮还没结束(缺口 #8): 不回退到目录级成功标记。
        return False
    if rc.run_ok('grep -q "Normal termination" ' + repr(d + '/mol.log') + ' 2>/dev/null'):
        return True
    if rc.run_ok('grep -q "Have a nice day" ' + repr(d + '/aims.out') + ' 2>/dev/null'):
        return True
    return _sh_te_end(rc, d)


def _sh_te_end(rc, d: str) -> bool:
    """TE.dat 是否扫到 tcontrol 声明的 $eend(与本地 te_reached_end 同判据)。

    取参一律用 awk(实测 `awk '/^\\$eend/{print $2}'` 稳健; grep -E 的 `$` 转义
    在多层引号下会静默失配, 曾导致 $eend 取空后被当 0 —— 若某 run 的 $eend≠0
    即会静默判错, 故此处不依赖 grep)。
    """
    cmd = (
        'e=$(awk \'/^\\$eend/{print $2; exit}\' ' + repr(d + '/tcontrol') + ' 2>/dev/null); '
        'ef=$(awk \'/efermi *=/{sub(/.*efermi *= */,""); sub(/ *H.*/,""); print; exit}\' '
        + repr(d + '/TE.dat') + ' 2>/dev/null); '
        '[ -n "$e" ] && [ -n "$ef" ] || exit 1; '
        'awk -v e="$e" -v ef="$ef" \'BEGIN{t=(e-ef)*27.211386-0.05; ok=0} '
        '/^[^#]/ && NF>=3 { if ($2+0 >= t) ok=1 } END{exit ok?0:1}\' ' + repr(d + '/TE.dat')
    )
    return rc.run_ok(cmd)


def smj_remote_finished(rc, rdir: str, run_dir=None) -> bool:
    # run_dir 给定时按阶段脚本决定判据(哨兵严格); 旧 run 无脚本 -> 旧回退链(向后兼容)。
    return _sh_finished(rc, rdir, sentinel_strict_script(run_dir))


def smj_remote_converged(rc, rdir: str, run_dir=None) -> bool:
    d = rdir
    sc = sentinel_strict_script(run_dir)
    sent = repr(d + '/' + DONE_SENTINEL)
    if rc.run_ok('test -s ' + sent + ' 2>/dev/null'):
        if not rc.run_ok('awk \'NR==1{exit ($1==0)?0:1}\' ' + sent):
            return False
    elif sc is not None:
        return False      # 哨兵没出现: 本轮还没结束, 谈不上收敛(缺口 #8)
    # 阶段脚本决定收敛判据: 共用目录里别的阶段的 aims.out/mol.log 不能冒充本轮结果。
    if sc == 'sub_te.sh':
        return _sh_te_end(rc, d)
    if sc == 'sub_aims.sh':
        return (rc.run_ok('grep -q "Have a nice day" ' + repr(d + '/aims.out') + ' 2>/dev/null')
                and (rc.run_ok('grep -q "Self-consistency cycle converged" '
                               + repr(d + '/aims.out') + ' 2>/dev/null')
                     or rc.run_ok('grep -q "Present geometry is converged" '
                                  + repr(d + '/aims.out') + ' 2>/dev/null')))
    if sc == 'sub_g16.sh':
        if not rc.run_ok('grep -q "Normal termination" ' + repr(d + '/mol.log') + ' 2>/dev/null'):
            return False
        return not rc.run_ok('grep -E "Frequencies --" ' + repr(d + '/mol.log')
                             + ' 2>/dev/null | grep -q -- "-"')
    if rc.run_ok('grep -q "Normal termination" ' + repr(d + '/mol.log') + ' 2>/dev/null'):
        # 虚频 = Frequencies 行含负值(24 章 §5); 无 Frequencies 行视为无虚频
        return not rc.run_ok('grep -E "Frequencies --" ' + repr(d + '/mol.log')
                             + ' 2>/dev/null | grep -q -- "-"')
    if rc.run_ok('grep -q "Have a nice day" ' + repr(d + '/aims.out') + ' 2>/dev/null'):
        return (rc.run_ok('grep -q "Self-consistency cycle converged" ' + repr(d + '/aims.out') + ' 2>/dev/null')
                or rc.run_ok('grep -q "Present geometry is converged" ' + repr(d + '/aims.out') + ' 2>/dev/null'))
    return _sh_te_end(rc, d)
