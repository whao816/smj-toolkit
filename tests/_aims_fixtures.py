"""FHI-aims 输入/输出片段的最小可信样本（供各测试复用）。

2026-11-20 复核 P1-3/P1-4 之后，``aims_scf.collect`` 与弛豫阶段开始校验
**内容**而不只是存在性：

* ``basis-indices.out`` / ``omat.aims`` / ``mos.aims`` 必须带 FHI-aims 的
  头部标识（见 ``src/smj/stages/_aims_artifacts.py``）；
* ``aims.out`` 必须含 ``Self-consistency cycle converged``（SCF）或
  ``Present geometry is converged.``（弛豫），并且都以 ``Have a nice day.`` 结束。

以前各测试统一写 ``"placeholder\\n"``，新口径下会被正确拒绝。这里集中提供
"真实形状"的最小片段，避免每个测试各写一份。
"""
from __future__ import annotations

from pathlib import Path

BASIS_INDICES = (
    "$basis.indices.aims\n"
    "# amount of basis functions used:   4364\n"
    "# fn.   type   at.   n   l   m\n"
    "   1     1     1    0   0   0\n"
)

OMAT = (
    "$omat.aims   format(4d20.14)   nsaos=  4364\n"
    "#column     1\n"
    "  1.00000000000000E+00\n"
)

MOS = "$scfmo.aims   iteration    1\n  1.00000000000000E+00\n"
ALPHA = "$scfmo.aims   alpha spin channel\n  1.00000000000000E+00\n"
BETA = "$scfmo.aims   beta spin channel\n  1.00000000000000E+00\n"

SCF_CONVERGED = "Self-consistency cycle converged.\nHave a nice day.\n"
RELAX_CONVERGED = "Present geometry is converged.\nHave a nice day.\n"


def write_scf_artifacts(stage_dir: Path, *, spin: bool = False) -> None:
    """写一套可被 ``aims_scf.collect`` 接受的 basis/omat/波函数文件。"""
    stage_dir.mkdir(parents=True, exist_ok=True)
    (stage_dir / "basis-indices.out").write_text(BASIS_INDICES, encoding="utf-8")
    (stage_dir / "omat.aims").write_text(OMAT, encoding="utf-8")
    if spin:
        (stage_dir / "alpha.aims").write_text(ALPHA, encoding="utf-8")
        (stage_dir / "beta.aims").write_text(BETA, encoding="utf-8")
    else:
        (stage_dir / "mos.aims").write_text(MOS, encoding="utf-8")


def write_scf_sidecars(stage_dir: Path) -> None:
    """只写 ``basis-indices.out`` + ``omat.aims``（自旋通道由调用方自行落盘）。"""
    stage_dir.mkdir(parents=True, exist_ok=True)
    (stage_dir / "basis-indices.out").write_text(BASIS_INDICES, encoding="utf-8")
    (stage_dir / "omat.aims").write_text(OMAT, encoding="utf-8")


def write_scf_log(stage_dir: Path, *, converged: bool = True,
                  normal_end: bool = True) -> None:
    lines = ["  FHI-aims version      : 221103\n",
             "  aims_uuid : TEST-UUID-1234\n"]
    if converged:
        lines.append(SCF_CONVERGED)
    elif normal_end:
        lines.append("Have a nice day.\n")
    (stage_dir / "aims.out").write_text("".join(lines), encoding="utf-8")


def write_relax_log(path: Path, *, converged: bool = True,
                    normal_end: bool = True) -> None:
    lines = ["  FHI-aims version      : 221103\n",
             "  No maximum number of relaxation steps, defaulting to     1000\n"]
    if converged:
        lines.append(RELAX_CONVERGED)
    else:
        lines.append("Present geometry is not yet converged.\n")
        if normal_end:
            lines.append("Have a nice day.\n")
    path.write_text("".join(lines), encoding="utf-8")
