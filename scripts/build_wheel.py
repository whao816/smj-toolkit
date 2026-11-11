"""用标准库构建 ``dist/smj_toolkit-<version>-py3-none-any.whl``。

为什么不用 ``python -m build`` / ``pip wheel``：本机系统 Python 没有 setuptools，
``ensurepip`` 在沙箱里也装不进去，而本包是纯 Python、无编译步骤 —— wheel 就是
一个按 PEP 427 布局的 zip。本脚本只依赖标准库，输出**可复现**（固定时间戳、
按路径排序、deflate），RECORD 里的 sha256 与文件一一对应，`pip install` 可校验。

用法::

    python scripts/build_wheel.py            # 构建到 dist/
    python scripts/build_wheel.py --check    # 构建后对照 RECORD 自检
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
DIST = ROOT / "dist"

#: 打包进 wheel 的文件模式（相对 src/）。模板与 VENDOR.json 都挂在
#: smj/_vendor/junction/ 下（`smj verify` 的模板检查也读那里），不是 smj/templates/。
INCLUDE_PATTERNS = (
    "smj/**/*.py",
    "smj/_vendor/junction/templates/*",
    "smj/_vendor/junction/VENDOR.json",
)

AUTHOR = "AgenticDFT contributors"
SUMMARY = "Single-molecule junction transport workflow (standalone)"
FIXED_DATE = (1980, 1, 1, 0, 0, 0)


def _version() -> str:
    import tomllib

    data = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    return str(data["project"]["version"])


def _package_version() -> str | None:
    import re

    text = (SRC / "smj" / "__init__.py").read_text(encoding="utf-8")
    match = re.search(r'^__version__\s*=\s*"([^"]+)"', text, re.MULTILINE)
    return match.group(1) if match else None


def _package_files() -> list[tuple[str, Path]]:
    """返回 [(wheel 内路径, 磁盘路径)]，顺序确定。"""
    out: dict[str, Path] = {}
    for pattern in INCLUDE_PATTERNS:
        for path in SRC.glob(pattern):
            if not path.is_file() or "__pycache__" in path.parts:
                continue
            rel = path.relative_to(SRC).as_posix()
            out[rel] = path
    return sorted(out.items())


def _metadata(version: str) -> str:
    return "\n".join([
        "Metadata-Version: 2.4",
        "Name: smj-toolkit",
        "Version: %s" % version,
        "Summary: %s" % SUMMARY,
        "Author: %s" % AUTHOR,
        "Requires-Python: >=3.10",
        "Requires-Dist: PyYAML<7,>=6",
        "Requires-Dist: numpy>=1.24",
        "Provides-Extra: plot",
        "Requires-Dist: matplotlib>=3.7; extra == \"plot\"",
        "Provides-Extra: dev",
        "Requires-Dist: pytest>=8; extra == \"dev\"",
        "Provides-Extra: all",
        "Requires-Dist: matplotlib>=3.7; extra == \"all\"",
        "Requires-Dist: pytest>=8; extra == \"all\"",
        "",
    ])


WHEEL_META = "\n".join([
    "Wheel-Version: 1.0",
    "Generator: smj-build-wheel (stdlib, 0.1.1 local)",
    "Root-Is-Purelib: true",
    "Tag: py3-none-any",
    "",
])

ENTRY_POINTS = "[console_scripts]\nsmj = smj.cli:main\n"


def _record_line(name: str, data: bytes) -> str:
    digest = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=")
    return "%s,sha256=%s,%d" % (name, digest.decode("ascii"), len(data))


def build() -> Path:
    version = _version()
    package_version = _package_version()
    if package_version != version:
        raise SystemExit("版本不一致：pyproject.toml=%s，smj.__version__=%s"
                         % (version, package_version))
    dist_info = "smj_toolkit-%s.dist-info" % version
    entries: list[tuple[str, bytes]] = []

    package_files = _package_files()
    if not any("/templates/" in name for name, _ in package_files):
        raise SystemExit("打包文件里没有 templates/*：检查 "
                         "src/smj/_vendor/junction/templates/ 是否存在")
    if not any(name.endswith("VENDOR.json") for name, _ in package_files):
        raise SystemExit("打包文件里没有 VENDOR.json：先运行 python scripts/gen_vendor_json.py")

    for rel, path in package_files:
        entries.append((rel, path.read_bytes()))

    entries.extend([
        (dist_info + "/METADATA", _metadata(version).encode("utf-8")),
        (dist_info + "/WHEEL", WHEEL_META.encode("utf-8")),
        (dist_info + "/entry_points.txt", ENTRY_POINTS.encode("utf-8")),
        (dist_info + "/top_level.txt", b"smj\n"),
    ])

    record = [_record_line(name, data) for name, data in entries]
    record.append(dist_info + "/RECORD,,")
    entries.append((dist_info + "/RECORD", ("\n".join(record) + "\n").encode("utf-8")))

    DIST.mkdir(parents=True, exist_ok=True)
    wheel = DIST / ("smj_toolkit-%s-py3-none-any.whl" % version)
    with zipfile.ZipFile(wheel, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name, data in entries:
            info = zipfile.ZipInfo(name, date_time=FIXED_DATE)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            zf.writestr(info, data)
    return wheel


def check(wheel: Path) -> int:
    """按 RECORD 逐条校验 sha256 与大小。"""
    problems: list[str] = []
    with zipfile.ZipFile(wheel) as zf:
        names = set(zf.namelist())
        record_name = next(n for n in names if n.endswith(".dist-info/RECORD"))
        for line in zf.read(record_name).decode("utf-8").splitlines():
            if not line.strip():
                continue
            name, digest, size = line.rsplit(",", 2)
            if name not in names:
                problems.append("RECORD 指向不存在的条目: %s" % name)
                continue
            data = zf.read(name)
            if digest:
                actual = "sha256=" + base64.urlsafe_b64encode(
                    hashlib.sha256(data).digest()).rstrip(b"=").decode("ascii")
                if actual != digest:
                    problems.append("sha256 不符: %s" % name)
            if size and len(data) != int(size):
                problems.append("大小不符: %s" % name)
        for name in sorted(names):
            if name.endswith("/RECORD"):
                continue
            if name.endswith(".dist-info/RECORD"):
                continue
    if problems:
        for item in problems:
            print("FAIL %s" % item)
        return 1
    print("RECORD 自检通过：%d 个条目" % len(names))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="构建后按 RECORD 自检")
    args = parser.parse_args(argv)

    vendor_json = SRC / "smj" / "_vendor" / "junction" / "VENDOR.json"
    if not vendor_json.is_file():
        print("缺少 %s：先运行 python scripts/gen_vendor_json.py" % vendor_json)
        return 1
    wheel = build()
    print("构建 %s（%d 字节）" % (wheel, wheel.stat().st_size))
    return check(wheel) if args.check else 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
