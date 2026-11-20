"""P2-5 (AGD-2026-09-25-SMJ): mock 集群 —— 无集群时回放远端产物。

缺口(实测): P0-3 改的是**远端**谓词(`junction/submit.py::smj_remote_finished/converged`),
但本会话只**在集群上手工验过**; 本地用例覆盖的是 local 分支。没有替身, 远端分支就永远
进不了 CI。

设计取舍:
- 只实现 smj 远端谓词实际发出的那几种命令形态(见下面 SENTINEL/GREP/TEST 分支);
  **不认识的一律抛 UnsupportedCommand** —— 宁可响, 不可静默给错答案。
- `_sh_te_end` 的复合命令交由**本地** `te_reached_end` 求值: 即"远端终点判据"由本地
  实现代理。这样 mock 下的**其余远端流程**是真测的, 而终点判据属"以本地为参考" ——
  **此边界必须记住**, 不要据此宣称远端终点判据已被端到端覆盖。
"""
from __future__ import annotations

import re
from pathlib import Path


class UnsupportedCommand(RuntimeError):
    """mock 不认识的命令: 显式报错(不静默给错答案)。"""


def _unq(s: str) -> str:
    s = s.strip()
    if len(s) >= 2 and s[0] == s[-1] and s[0] in ("'", '"'):
        return s[1:-1]
    return s


def _split_paths(s: str) -> list:
    """按 shell 的引号切分路径参数; **不做反斜杠转义**(Windows 路径会被 shlex.split 吃掉)。"""
    out = []
    for m in re.finditer(r"""'([^']*)'|"([^"]*)"|(\S+)""", s):
        out.append(next(g for g in m.groups() if g is not None))
    return out


class MockTransport:
    """把本地目录当远端的最小替身。"""

    def __init__(self, root):
        self.root = Path(root)
        self.calls = []

    # ---- 命令执行 -------------------------------------------------------
    def run_ok(self, cmd: str) -> bool:
        return self.run(cmd)[0] == 0

    def run(self, cmd: str, timeout=None):
        self.calls.append(cmd)
        c = cmd.replace(' 2>/dev/null', '').replace(' 2>&1', '').strip()
        # P0-3 哨兵: awk 'NR==1{exit ($1==0)?0:1}' <sentinel>
        m = re.match(r"awk 'NR==1\{exit \(\$1==0\)\?0:1\}' (.+)$", c)
        if m:
            p = Path(_unq(m.group(1)))
            if not p.exists():
                return (1, '', '')
            try:
                return (0 if int(p.read_text().split()[0]) == 0 else 1, '', '')
            except (IndexError, ValueError):
                return (1, '', '')
        # 复合: test -s <p> && awk 'NR==1{exit ($1==0)?0:1}' <p>   (哨兵优先判定)
        m = re.match(r"test -s (.+?) && awk 'NR==1\{exit \(\$1==0\)\?0:1\}' (.+)$", c)
        if m:
            p1, p2 = Path(_unq(m.group(1))), Path(_unq(m.group(2)))
            if not (p1.exists() and p1.stat().st_size > 0):
                return (1, '', '')
            try:
                return (0 if int(p2.read_text().split()[0]) == 0 else 1, '', '')
            except (IndexError, ValueError):
                return (1, '', '')
        # test -s <path>
        if c.startswith('test -s '):
            rest = c[len('test -s '):]
            if ' && ' in rest or '|' in rest:
                # 未实现的复合形态必须**响亮失败**, 不得静默误解析(本会话已踩过一次)
                raise UnsupportedCommand('mock 未实现的复合 test 命令: %s' % cmd[:120])
            p = Path(_unq(rest))
            return (0 if (p.exists() and p.stat().st_size > 0) else 1, '', '')
        # grep -q "<pat>" <path>
        m = re.match(r'grep -q "([^"]*)" (.+)$', c)
        if m:
            pat, p = m.group(1), Path(_unq(m.group(2)))
            if not p.exists():
                return (1, '', '')
            return (0 if pat in p.read_text(encoding='utf-8', errors='ignore') else 1, '', '')
        # grep -E "Frequencies --" <path> | grep -q -- "-"
        m = re.match(r'grep -E "Frequencies --" (.+?) \| grep -q -- "-"$', c)
        if m:
            p = Path(_unq(m.group(1)))
            if not p.exists():
                return (1, '', '')
            hits = [l for l in p.read_text(encoding='utf-8', errors='ignore').splitlines()
                    if 'Frequencies --' in l]
            return (0 if any('-' in l for l in hits) else 1, '', '')
        # _sh_te_end 复合命令: 交给本地 te_reached_end(见模块 docstring 的边界说明)
        if '/^$eend/' in c or '\\$eend' in c:
            from .submit import te_reached_end
            # 命令里有多个单引号段(awk 脚本本身也被引起来): 只取**不含引号**、且以
            # /TE.dat 结尾的那一段, 否则会把 awk 脚本的开引号当成路径起点(静默取到错目录)。
            cands = re.findall(r"'([^']*?/TE\.dat)'", c)
            d = Path(cands[-1]).parent if cands else None
            try:
                return (0 if (d and te_reached_end(d)) else 1, '', '')
            except Exception:
                return (1, '', '')
        # R4(2026-10-23 smj 提交/回取接线): head -N <path> -> 前 N 行(与真实 head 同形态)。
        # smj 回取用它对 mos.aims/omat.aims 只取头部(GB 级文件不整取)。
        m = re.match(r'head -(\d+) (.+)$', c)
        if m:
            n, p = int(m.group(1)), Path(_unq(m.group(2)))
            if not p.exists():
                return (1, '', '')
            try:
                lines = p.read_text(encoding='utf-8', errors='ignore').splitlines(True)
            except OSError:
                return (1, '', '')
            return (0, ''.join(lines[:n]), '')
        # mkdir -p <path> [<path>...]
        m = re.match(r'mkdir -p (.+)$', c)
        if m:
            ok = True
            for p in _split_paths(m.group(1)):
                try:
                    Path(p).mkdir(parents=True, exist_ok=True)
                except OSError:
                    ok = False
            return (0 if ok else 1, '', '')
        # sbatch 提交形态: cd <dir> && sbatch <script> [args]
        # 只认本仓提交层实际发出的这一种(带 cd 且脚本**必须已落在该目录**), 借此把
        # "aitranss 的 cd 指向错误远端目录"变成可见失败, 而不是静默返回 job id。
        m = re.match(r"cd (.+?) && sbatch (\S+)(?: .*)?$", c)
        if m:
            d = Path(_unq(m.group(1)))
            script = _unq(m.group(2))
            if not (d.is_dir() and (d / script).is_file()):
                return (1, '', 'mock: 远端目录缺阶段脚本 %s' % (d / script))
            return (0, 'Submitted batch job 100001', '')
        # cd <dir> && rm -f <names...>  (smj 阶段提交前清远端残留, 见 submit._smj_stale_names)
        m = re.match(r"cd (.+?) && rm -f (.+)$", c)
        if m:
            d = Path(_unq(m.group(1)))
            if not d.is_dir():
                return (1, '', 'mock: 远端目录不存在 %s' % d)
            for tok in m.group(2).split():
                try:
                    victims = sorted(d.glob(tok)) if any(ch in tok for ch in '*?[') else [d / tok]
                except (OSError, ValueError):
                    return (1, '', 'mock: 无法解析清理目标 %r' % tok)
                for p in victims:
                    try:
                        if p.is_file() or p.is_symlink():
                            p.unlink()
                    except OSError:
                        return (1, '', 'mock: 清理失败 %s' % p)
            return (0, '', '')
        raise UnsupportedCommand('mock 未实现的命令形态: %s' % cmd[:120])

    # ---- 文件面(供 stage 复用) -------------------------------------------
    def exists(self, path) -> bool:
        return Path(path).exists()

    def read_bytes(self, path, timeout=None):
        """读远端字节; 不存在/不可读返回 None(与 Transport.read_bytes 契约一致)。"""
        try:
            return Path(path).read_bytes()
        except OSError:
            return None

    def mkdir_p(self, path) -> None:
        Path(path).mkdir(parents=True, exist_ok=True)

    def write_text(self, text: str, path) -> bool:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding='utf-8')
        return True

    def write_bytes(self, data: bytes, path) -> bool:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(data)
        return True
