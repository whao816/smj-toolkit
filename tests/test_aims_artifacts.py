"""``_aims_artifacts`` 头部校验：自旋极化写出的 alpha/beta.aims 必须被接受。

2026-11-20 修：FHI-aims 221103 在 spin collinear 下写 ``alpha.aims`` /
``beta.aims``，头部标识是 ``$uhfmo_alpha.aims`` / ``$uhfmo_beta.aims``；
旧表只认 ``$scfmo.aims``，真实产物会在 collect 时被判成"文件被截断"。
"""
from __future__ import annotations

from pathlib import Path




from smj.stages._aims_artifacts import HEADERS, probe, problem


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_real_fhi_aims_uhfmo_headers_accepted(workdir) -> None:
    alpha = _write(workdir / "alpha.aims",
                   "$uhfmo_alpha.aims   iteration:   90   format(4d20.14)\n"
                   "     1  a      eigenvalue=-.3256D+04   nsaos=  4934\n")
    beta = _write(workdir / "beta.aims",
                  "$uhfmo_beta.aims   iteration:   90   format(4d20.14)\n"
                  "     1  b      eigenvalue=-.3256D+04   nsaos=  4934\n")
    for path in (alpha, beta):
        info = probe(path)
        assert info["header_ok"] is True
        assert problem(info) is None


def test_legacy_scfmo_alpha_still_accepted(workdir) -> None:
    info = probe(_write(workdir / "alpha.aims",
                        "$scfmo.aims   alpha spin channel\n  1.0\n"))
    assert info["header_ok"] is True
    assert problem(info) is None


def test_truncated_alpha_rejected(workdir) -> None:
    info = probe(_write(workdir / "alpha.aims", "garbage\n"))
    assert info["header_ok"] is False
    reason = problem(info)
    assert reason is not None and "头部缺少" in reason


def test_zero_byte_rejected(workdir) -> None:
    empty = workdir / "omat.aims"
    empty.write_bytes(b"")
    assert "零字节" in (problem(probe(empty)) or "")


def test_header_table_covers_every_tested_name() -> None:
    assert set(HEADERS) == {
        "basis-indices.out", "omat.aims", "mos.aims", "alpha.aims", "beta.aims"}
    for markers in HEADERS.values():
        assert isinstance(markers, tuple) and markers


