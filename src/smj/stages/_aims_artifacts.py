"""校验 FHI-aims SCF 产物的"真实性"（非空 + 头部标识 + 可解析规模）。

2026-11-20 复核 P1-4：``aims_scf.collect`` 只用 ``is_file()`` 判存在，零字节的
``basis-indices.out`` / ``omat.aims`` / ``mos.aims`` 也能通过，而 ``aitranss_te``
的 ``tcontrol`` 会把它们当输入直接读。这里按文件头做最小真实性校验：

===================  ===========================================================
文件                 期望头部（FHI-aims 自己写的标识行）
===================  ===========================================================
basis-indices.out    ``$basis.indices.aims``（后随 ``# amount of basis functions used: N``）
omat.aims            ``$omat.aims   format(4d20.14)   nsaos=  4364``
mos.aims             ``$scfmo.aims``（闭壳层；自旋极化时改为 alpha/beta.aims）
alpha.aims / beta.aims 同上
===================  ===========================================================

只读每个文件的前 ``_HEAD_BYTES`` 字节，不把大文件读进内存。
"""

from __future__ import annotations

import re
from pathlib import Path

_HEAD_BYTES = 4096

#: 文件名 → 必须出现在头部的标识串。
HEADERS = {
    "basis-indices.out": "$basis.indices.aims",
    "omat.aims": "$omat.aims",
    "mos.aims": "$scfmo.aims",
    "alpha.aims": "$scfmo.aims",
    "beta.aims": "$scfmo.aims",
}

_BASIS_COUNT = re.compile(r"#\s*amount of basis functions used:\s*(\d+)")
_NSAOS = re.compile(r"nsaos\s*=\s*(\d+)")


def probe(path: str | Path) -> dict:
    """读头部并返回 ``{name, bytes, header, header_ok, size_value}``。

    文件不存在时抛 ``FileNotFoundError``；存在但读失败时抛 ``OSError``。
    """
    path = Path(path)
    size = path.stat().st_size
    with path.open("rb") as fh:
        head = fh.read(_HEAD_BYTES)
    text = head.decode("utf-8", "replace")
    expected = HEADERS.get(path.name)
    header_ok = True
    if expected is not None:
        header_ok = expected in text
    size_value = None
    if path.name == "basis-indices.out":
        match = _BASIS_COUNT.search(text)
        if match:
            size_value = int(match.group(1))
    elif path.name == "omat.aims":
        match = _NSAOS.search(text)
        if match:
            size_value = int(match.group(1))
    return {
        "name": path.name,
        "path": str(path),
        "bytes": size,
        "header": expected,
        "header_ok": header_ok,
        "size_value": size_value,
    }


def problem(info: dict) -> str | None:
    """返回不合格原因；合格时返回 ``None``。"""
    if info["bytes"] <= 0:
        return "%s 是零字节文件" % info["name"]
    if not info["header_ok"]:
        return ("%s 头部缺少 %r 标识（文件可能被截断或不是 FHI-aims 产物）"
                % (info["name"], info["header"]))
    if info["size_value"] is not None and info["size_value"] <= 0:
        return "%s 里记录的规模为 %d（不可能）" % (info["name"], info["size_value"])
    return None
