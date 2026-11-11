"""测试公共夹具。

本机 %TEMP% 处在 DSH 沙箱的写限制之外（在 ``%TEMP%\\dsh-*`` 下 mkdir 报
``WinError 5``），pytest 自带的 ``tmp_path`` 会因此在 ``tmp_path_factory.mktemp``
处失败，且它建的 basetemp 目录事后连 ``Get-Acl`` 都无权访问、删不掉。
所以所有测试统一用仓库内 ``.tmp/`` 下自建自删的 ``workdir``（已在 .gitignore）。
"""
from __future__ import annotations

import shutil
import uuid
from pathlib import Path
from typing import Iterator

import pytest


@pytest.fixture()
def workdir() -> Iterator[Path]:
    base = Path(".tmp") / ("t-" + uuid.uuid4().hex[:8])
    base.mkdir(parents=True)
    try:
        yield base
    finally:
        shutil.rmtree(base, ignore_errors=True)
