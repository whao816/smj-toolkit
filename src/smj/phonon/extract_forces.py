"""``aims.out`` → ``forces.json``(远程阵列作业里用,**只依赖 stdlib**)。

集群上跑阵列作业时不保证装了 ``smj``(而且不该为了一个正则去装整条链),
所以这个脚本自成一体:只用 ``re``/``json``/``sys``。解析逻辑与
:func:`smj.phonon.forces.read_forces` **刻意重复**一份,由
``tests/phonon/test_forces.py`` 用同一夹具逐位交叉验证 —— 重复可以,分叉不行。

用法::

    python extract_forces.py aims.out forces.json
"""
from __future__ import annotations

import json
import re
import sys

_TITLE = "Total atomic forces"
_ROW_RE = re.compile(r"^\s*\|\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$")

USAGE = "用法: python extract_forces.py <aims.out> <forces.json>"


def _to_float(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


def extract(text: str) -> dict:
    """返回 ``{"n_atoms", "forces", "block_index", "n_blocks"}``。"""
    lines = text.splitlines()
    starts = [i for i, line in enumerate(lines) if _TITLE in line]
    if not starts:
        raise ValueError(
            "aims.out 里没有 'Total atomic forces' 段:control.in 是否漏了 "
            "`compute_forces .true.`?(注意 `calculate_forces` 不是关键字)")
    start = starts[-1]
    rows: list[tuple[int, float, float, float]] = []
    index = start + 1
    while index < len(lines):
        match = _ROW_RE.match(lines[index])
        if match is None:
            if lines[index].strip():
                break
            index += 1
            continue
        rows.append((int(match.group(1)), _to_float(match.group(2)),
                     _to_float(match.group(3)), _to_float(match.group(4))))
        index += 1
    if not rows:
        raise ValueError("力块标题在第 %d 行,但没有任何力行" % (start + 1))
    expected = list(range(1, len(rows) + 1))
    if [row[0] for row in rows] != expected:
        raise ValueError(
            "力块原子序号不连续: 首尾 %d..%d,共 %d 行"
            % (rows[0][0], rows[-1][0], len(rows)))
    return {
        "n_atoms": len(rows),
        "forces": [[row[1], row[2], row[3]] for row in rows],
        "block_index": len(starts) - 1,
        "n_blocks": len(starts),
    }


def main(argv=None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in ("-h", "--help"):
        print(USAGE)
        return 0 if args else 2
    if len(args) != 2:
        print(USAGE, file=sys.stderr)
        return 2
    source, target = args
    try:
        with open(source, "r", encoding="utf-8", errors="replace") as handle:
            payload = extract(handle.read())
    except (OSError, ValueError) as error:
        print("extract_forces: %s" % error, file=sys.stderr)
        return 1
    with open(target, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=1)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
