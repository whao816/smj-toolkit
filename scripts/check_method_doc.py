# -*- coding: utf-8 -*-
"""Rebuild the method document in a temp dir and compare with the repo copy.

    python scripts/check_method_doc.py

Checks both files:

* ``docs/AutoSMJ_DFT方法说明_v1.2.docx`` — rebuilt bytes (the builder fixes
  core-property dates and zip timestamps, so the same toolchain is
  byte-identical; if another zlib emits different bytes the member contents
  are compared instead);
* ``docs/AutoSMJ_DFT方法说明_v1.2.md`` — derived from the rebuilt DOCX, must be
  byte-identical.

Exit code 0 when everything matches, 1 otherwise.
"""
from __future__ import annotations

import hashlib
import shutil
import subprocess
import sys
import uuid
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRATCH = ROOT / ".tmp"
DOCX = ROOT / "docs" / "AutoSMJ_DFT方法说明_v1.2.docx"
MD = ROOT / "docs" / "AutoSMJ_DFT方法说明_v1.2.md"
BUILD_DOCX = ROOT / "scripts" / "build_method_doc_docx.py"
BUILD_MD = ROOT / "scripts" / "build_method_doc_md.py"


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _docx_entries(path: Path) -> dict[str, bytes]:
    with zipfile.ZipFile(path) as zf:
        return {name: zf.read(name) for name in zf.namelist()}


def main() -> int:
    for path in (DOCX, MD):
        if not path.is_file():
            print(f"[fail] missing {path.relative_to(ROOT)}")
            return 1

    SCRATCH.mkdir(parents=True, exist_ok=True)
    # NOTE: tempfile.mkdtemp sets mode 0o700, whose explicit DACL some Windows
    # sandboxes deny to child processes.  A plainly created scratch dir inherits
    # the workspace ACL and writes fine.
    tmp = SCRATCH / f"smj-method-doc-{uuid.uuid4().hex}"
    tmp.mkdir()
    try:
        rebuilt_docx = tmp / "rebuilt.docx"
        rebuilt_md = tmp / "rebuilt.md"
        subprocess.run(
            [sys.executable, str(BUILD_DOCX), str(rebuilt_docx)], check=True
        )
        subprocess.run(
            [sys.executable, str(BUILD_MD), str(rebuilt_docx), str(rebuilt_md)],
            check=True,
        )

        ok = True

        committed_docx = DOCX.read_bytes()
        fresh_docx = rebuilt_docx.read_bytes()
        if committed_docx == fresh_docx:
            print(
                "[ok] DOCX byte-identical "
                f"({len(committed_docx)} bytes, sha256 {_sha256(committed_docx)[:16]}…)"
            )
        elif _docx_entries(DOCX) == _docx_entries(rebuilt_docx):
            print(
                "[warn] DOCX zip bytes differ but all member contents match "
                "(different zlib/toolchain); rebuild on this machine to refresh"
            )
        else:
            ok = False
            old, new = _docx_entries(DOCX), _docx_entries(rebuilt_docx)
            changed = sorted(
                name for name in set(old) | set(new) if old.get(name) != new.get(name)
            )
            print("[fail] DOCX member contents differ:", ", ".join(changed))

        committed_md = MD.read_bytes()
        fresh_md = rebuilt_md.read_bytes()
        if committed_md == fresh_md:
            print(
                "[ok] MD byte-identical "
                f"({len(committed_md)} bytes, sha256 {_sha256(committed_md)[:16]}…)"
            )
        else:
            ok = False
            print(
                "[fail] MD differs: committed "
                f"{len(committed_md)} bytes vs rebuilt {len(fresh_md)} bytes"
            )
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print("[pass] method document is reproducible" if ok else "[fail] method document drifted")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
