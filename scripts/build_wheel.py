"""用标准库构建 ``dist/smj_toolkit-<version>-py3-none-any.whl``（可选 sdist）。

为什么不用 ``python -m build`` / ``pip wheel``：本机系统 Python 没有 setuptools，
``ensurepip`` 在沙箱里也装不进去，而本包是纯 Python、无编译步骤 —— wheel 就是
一个按 PEP 427 布局的 zip，sdist 就是一个按 PEP 517 布局的 tar.gz。本脚本只依赖
标准库，输出**可复现**（固定时间戳、按路径排序、deflate/gzip），RECORD 里的 sha256
与文件一一对应，`pip install` 可校验。

用法::

    python scripts/build_wheel.py             # 只构建 wheel 到 dist/
    python scripts/build_wheel.py --check     # 构建后对照 RECORD 自检
    python scripts/build_wheel.py --sdist     # 同时构建 sdist（tar.gz）
    python scripts/build_wheel.py --check --sdist
"""
from __future__ import annotations

import argparse
import base64
import gzip
import hashlib
import io
import sys
import tarfile
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
FIXED_EPOCH = 315532800  # 1980-01-01T00:00:00Z，sdist 的 tar/gzip 时间戳

#: 兜底 classifiers：pyproject 读不到时使用；正常路径以 pyproject.toml 为准，避免元数据漂移。
CLASSIFIERS = (
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Operating System :: POSIX :: Linux",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Programming Language :: Python :: 3.13",
    "Programming Language :: Python :: 3.14",
    "Topic :: Scientific/Engineering :: Chemistry",
    "Topic :: Scientific/Engineering :: Physics",
)


def _classifiers() -> tuple[str, ...]:
    """优先从 pyproject.toml 读 classifiers（Python 3.10 无 tomllib 时回退常量）。"""
    try:
        import tomllib
    except ModuleNotFoundError:  # pragma: no cover - Python 3.10
        return CLASSIFIERS
    data = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    items = data.get("project", {}).get("classifiers")
    if isinstance(items, list) and items:
        return tuple(str(item) for item in items)
    return CLASSIFIERS


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
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    return "\n".join([
        "Metadata-Version: 2.4",
        "Name: smj-toolkit",
        "Version: %s" % version,
        "Summary: %s" % SUMMARY,
        "Author: %s" % AUTHOR,
        "License-Expression: MIT",
        "License-File: LICENSE",
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
        "Keywords: single-molecule junction, transport, FHI-aims, aitranss, DFT",
    ] + ["Classifier: %s" % c for c in _classifiers()] + [
        "Description-Content-Type: text/markdown",
        "",
        readme,
        "",
    ])


WHEEL_META = "\n".join([
    "Wheel-Version: 1.0",
    "Generator: smj-build-wheel (stdlib)",
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
    license_file = ROOT / "LICENSE"
    if not license_file.is_file():
        raise SystemExit("缺少 LICENSE：发布包必须带许可证")

    for rel, path in package_files:
        entries.append((rel, path.read_bytes()))

    entries.extend([
        (dist_info + "/licenses/LICENSE", license_file.read_bytes()),
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


def _sdist_files(version: str) -> list[tuple[str, bytes]]:
    files: list[tuple[str, bytes]] = []
    for name in ("pyproject.toml", "README.md", "LICENSE", "CHANGELOG.md"):
        files.append((name, (ROOT / name).read_bytes()))
    for rel, path in _package_files():
        files.append(("src/" + rel, path.read_bytes()))
    for path in sorted((ROOT / "tests").glob("*.py")):
        files.append(("tests/" + path.name, path.read_bytes()))
    for path in sorted((ROOT / "scripts").glob("*.py")):
        files.append(("scripts/" + path.name, path.read_bytes()))
    files.append(("PKG-INFO", _metadata(version).encode("utf-8")))
    return files


def build_sdist(version: str) -> Path:
    """按 PEP 517 布局写 tar.gz（固定 mtime/uid/gid，可复现）。"""
    root = "smj_toolkit-%s" % version
    DIST.mkdir(parents=True, exist_ok=True)
    archive = DIST / ("%s.tar.gz" % root)
    with gzip.GzipFile(str(archive), "wb", mtime=FIXED_EPOCH) as gz:
        with tarfile.open(fileobj=gz, mode="w", format=tarfile.PAX_FORMAT) as tf:
            for name, data in sorted(_sdist_files(version)):
                info = tarfile.TarInfo("%s/%s" % (root, name))
                info.size = len(data)
                info.mtime = FIXED_EPOCH
                info.mode = 0o644
                info.uid = info.gid = 0
                info.uname = info.gname = ""
                tf.addfile(info, io.BytesIO(data))
    return archive


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
    if problems:
        for item in problems:
            print("FAIL %s" % item)
        return 1
    print("RECORD 自检通过：%d 个条目" % len(names))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="构建后按 RECORD 自检")
    parser.add_argument("--sdist", action="store_true", help="同时构建 sdist（tar.gz）")
    args = parser.parse_args(argv)

    vendor_json = SRC / "smj" / "_vendor" / "junction" / "VENDOR.json"
    if not vendor_json.is_file():
        print("缺少 %s：先运行 python scripts/gen_vendor_json.py" % vendor_json)
        return 1
    wheel = build()
    print("构建 %s（%d 字节）" % (wheel, wheel.stat().st_size))
    code = check(wheel) if args.check else 0
    if code:
        return code
    if args.sdist:
        archive = build_sdist(_version())
        print("构建 %s（%d 字节）" % (archive, archive.stat().st_size))
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
