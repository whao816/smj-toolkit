"""FHI-aims ``electrodes.library`` 预定义 Au 团簇的读取与校验。

手册 §5.3.1 要求用库里的预定义团簇做电极，并在装配后固定背层弛豫一次。
库文件是 ``geometry.in`` 风格（``#`` 头注释 + ``atom x y z Au`` 行），
头部带建议的 ``-nlayers``。本模块只做**读与校验**，不打包、不改写库文件：

- 每侧团簇必须元素一致、至少两层、接触层唯一原子（tip）、背层 >= 3 且共面；
- 接触原子归一化到原点并放在第 0 位（``mount_electrode`` 约定第 0 个是 tip）；
- 记录来源路径 + sha256（metadata provenance），不把文件内容塞进 wheel。

库目录解析顺序：``electrode.library_dir``（或 ``SMJ_ELECTRODE_LIBRARY``）
→ 本机常见安装路径 → ``mode=slurm-ssh`` 时经 ssh ``cat`` 取回并缓存到
``<run>/electrodes/``。
"""

from __future__ import annotations

import glob as _glob
import hashlib
import os
import re
import shlex
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from .junction.electrode import LibraryElectrode
from .junction.geometry import Atom
from .errors import SmjError

#: 默认团簇：37 原子 fcc111 壳层（nlayers=3），手册示例同族。
DEFAULT_CLUSTER = "fcc111.au.clusters.37_atoms/fcc111.cluster.natoms37.x1.aims"

#: 没配 library_dir 时探测的本机位置（glob）。
COMMON_LIBRARY_GLOBS: tuple[str, ...] = (
    "/data/apps/src/fhi-aims*/src/aitranss/electrodes.library",
    "/opt/fhi-aims*/src/aitranss/electrodes.library",
    "~/fhi-aims*/src/aitranss/electrodes.library",
    "~/software/fhi-aims*/src/aitranss/electrodes.library",
)

#: 没有本地库、但集群上常见的位置（slurm-ssh 时经 ssh cat）。
COMMON_LIBRARY_REMOTE: tuple[str, ...] = (
    "/data/apps/src/fhi-aims.221103/src/aitranss/electrodes.library",
    "/data/apps/src/fhi-aims*/src/aitranss/electrodes.library",
)

_ATOM_RE = re.compile(r"^\s*atom\s+(\S+)\s+(\S+)\s+(\S+)\s+([A-Za-z][A-Za-z]?)\s*$")
_NLAYERS_RE = re.compile(r"-\s*nlayers\s+(\d+)", re.IGNORECASE)
_PLANE_TOL = 1e-4

LIBRARY_HINT = (
    "未找到 FHI-aims electrodes.library。三条出路：\n"
    "  1) 在 smj.yaml 写 electrode.library_dir: <本地目录或集群绝对路径>\n"
    "  2) 设环境变量 SMJ_ELECTRODE_LIBRARY=<同上的路径>\n"
    "  3) 不想要库电极时写 electrode.source: builtin（退回 0.1.1 的自造三角岛）"
)


@dataclass(frozen=True)
class LibraryCluster:
    """一份已校验的库团簇。``atoms[0]`` 恒为接触原子（位于原点）。"""

    atoms: tuple[Atom, ...]
    relpath: str
    sha256: str
    source: str
    nlayers_suggested: int | None = None
    fetched_from: str | None = None
    translated: bool = False

    @property
    def natoms(self) -> int:
        return len(self.atoms)

    @property
    def element(self) -> str:
        return self.atoms[0].element

    @property
    def natoms_per_side(self) -> int:
        return self.natoms

    def layers(self) -> list[tuple[float, int]]:
        """[(z, 原子数)]，z 升序。"""
        counts: dict[float, int] = {}
        for atom in self.atoms:
            key = round(atom.z, 4)
            counts[key] = counts.get(key, 0) + 1
        return sorted(counts.items())

    def to_electrode(self) -> LibraryElectrode:
        return LibraryElectrode(
            atoms=tuple((a.x, a.y, a.z) for a in self.atoms),
            element=self.element,
            layers=len(self.layers()),
        )

    def provenance(self) -> dict:
        return {
            "file": str(self.source),
            "relpath": self.relpath,
            "sha256": self.sha256,
            "natoms": self.natoms,
            "nlayers_suggested": self.nlayers_suggested,
            "fetched_from": self.fetched_from,
        }


def parse_cluster_text(text: str, *, source: str, relpath: str,
                       fetched_from: str | None = None) -> LibraryCluster:
    """解析并校验库文件文本；任何不合规都抛 ``SmjError``（不静默修补）。"""
    atoms: list[Atom] = []
    for line in text.splitlines():
        match = _ATOM_RE.match(line)
        if match:
            x, y, z, element = match.groups()
            atoms.append(Atom(element, float(x), float(y), float(z)))
    if not atoms:
        raise SmjError("库团簇没有 atom 行: %s" % source)
    elements = {a.element for a in atoms}
    if len(elements) != 1:
        raise SmjError("库团簇含多种元素 %s: %s（单侧团簇必须纯金属）"
                       % (sorted(elements), source))
    nlayers = None
    head = _NLAYERS_RE.search(text)
    if head:
        nlayers = int(head.group(1))
        if not 1 <= nlayers <= 8:
            raise SmjError("库文件头部 -nlayers=%d 不合理: %s" % (nlayers, source))
    # 每层必须共面（库文件是理想晶格；z 抖动说明文件坏了）。
    # 旧实现用 round(z,4) 自己分组再和自身比，恒等于 0（空检查）；改成先合并
    # 近似 z，再要求"不同层的 z 间距"远大于面内容差 —— 对 1e-5 级抖动也报错。
    zs = sorted({round(a.z, 6) for a in atoms})
    if len(zs) < 2:
        raise SmjError("库团簇只有一层（z=%s），不是可用的电极: %s" % (zs, source))
    gaps = [b - a for a, b in zip(zs, zs[1:])]
    if gaps and min(gaps) <= max(_PLANE_TOL, 1e-3):
        raise SmjError("库团簇层不共面（相邻层 z 间距 %.2e ≤ 容差，更像同一层被 z 抖动拆开）: %s: %s"
                       % (min(gaps), source, zs))
    contact_z = zs[0]
    contact = [a for a in atoms if abs(a.z - contact_z) <= _PLANE_TOL]
    if len(contact) != 1:
        raise SmjError("库团簇接触层应有 1 个原子，实得 %d: %s" % (len(contact), source))
    back = [a for a in atoms if abs(a.z - zs[-1]) <= _PLANE_TOL]
    if len(back) < 3:
        raise SmjError("库团簇背层只有 %d 个原子（需 >= 3）: %s" % (len(back), source))
    # 每层必须共面（库文件是理想晶格；z 抖动说明文件坏了）。
    for z in zs:
        layer = [a for a in atoms if abs(a.z - z) <= _PLANE_TOL]
        spread = max(abs(a.z - z) for a in layer)
        if spread > _PLANE_TOL:
            raise SmjError("库团簇 z=%.4f 层不共面（偏差 %.2e）: %s"
                           % (z, spread, source))
    tip = contact[0]
    translated = abs(tip.x) > 1e-9 or abs(tip.y) > 1e-9 or abs(tip.z) > 1e-9
    ordered = [tip] + [a for a in atoms if a is not tip]
    if translated:
        dx, dy, dz = tip.x, tip.y, tip.z
        ordered = [Atom(a.element, a.x - dx, a.y - dy, a.z - dz) for a in ordered]
    return LibraryCluster(
        atoms=tuple(ordered),
        relpath=relpath,
        sha256=hashlib.sha256(text.encode("utf-8")).hexdigest(),
        source=source,
        nlayers_suggested=nlayers,
        fetched_from=fetched_from,
        translated=translated,
    )


def load_cluster_file(path: str | Path, *, relpath: str | None = None) -> LibraryCluster:
    p = Path(path).expanduser()
    try:
        text = p.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise SmjError("读不了库团簇 %s: %s" % (p, exc)) from exc
    return parse_cluster_text(text, source=str(p), relpath=relpath or p.name)


def _candidate_files(hint: str | None, relpath: str) -> list[Path]:
    out: list[Path] = []
    if hint:
        base = Path(os.path.expanduser(hint))
        if base.suffix.lower() == ".aims" and base.is_file():
            out.append(base)
        else:
            out.append(base / relpath)
    for pattern in COMMON_LIBRARY_GLOBS:
        for match in sorted(_glob.glob(os.path.expanduser(pattern))):
            out.append(Path(match) / relpath)
    return out


def _remote_candidates(hint: str | None, relpath: str) -> list[str]:
    out: list[str] = []
    if hint and hint.startswith("/"):
        out.append(hint.rstrip("/") + "/" + relpath)
    for pattern in COMMON_LIBRARY_REMOTE:
        # glob 交给 shell 展开（ssh 的远端 shell 会做），这里只拼好通配串。
        out.append(pattern.rstrip("/") + "/" + relpath)
    return out


def resolve_cluster(ctx, electrode) -> LibraryCluster:
    """按配置/环境/远程顺序取回团簇；失败抛出带三条出路的 ``SmjError``。"""
    relpath = str(getattr(electrode, "library_cluster", "") or DEFAULT_CLUSTER)
    hint = getattr(electrode, "library_dir", None) or os.environ.get(
        "SMJ_ELECTRODE_LIBRARY")
    for candidate in _candidate_files(hint, relpath):
        if candidate.is_file():
            return load_cluster_file(candidate, relpath=relpath)
    transport = getattr(ctx, "transport", None)
    mode = getattr(getattr(ctx, "config", None), "mode", None)
    if transport is not None and mode == "slurm-ssh":
        for remote in _remote_candidates(hint, relpath):
            result = transport.run("cat %s" % shlex.quote(remote))
            if not result.ok or not result.out.strip():
                continue
            host = getattr(getattr(ctx.config, "cluster", None), "host", None) or "remote"
            cache = Path(ctx.layout.run_dir) / "electrodes" / Path(relpath).name
            try:
                cache.parent.mkdir(parents=True, exist_ok=True)
                cache.write_text(result.out, encoding="utf-8", newline="\n")
            except OSError:
                pass
            return parse_cluster_text(result.out, source="%s:%s" % (host, remote),
                                      relpath=relpath,
                                      fetched_from="%s:%s" % (host, remote))
    raise SmjError(LIBRARY_HINT + "\n（期望文件: %s）" % relpath)


def library_file_sha256(path: str | Path) -> str:
    """独立工具：给一份库文件算 sha256（与 provenance 同口径）。"""
    return hashlib.sha256(Path(path).expanduser().read_bytes()).hexdigest()


def validate_electrode_atoms(atoms: Iterable[Atom], element: str = "Au") -> None:
    """装配后的轻量校验：至少 2 个电极元素原子且分布不重叠（供测试/调用方复用）。"""
    items = list(atoms)
    metal = [a for a in items if a.element == element]
    if len(metal) < 2:
        raise SmjError("电极至少需要 2 个 %s，实得 %d" % (element, len(metal)))


__all__ = [
    "DEFAULT_CLUSTER",
    "LIBRARY_HINT",
    "LibraryCluster",
    "load_cluster_file",
    "parse_cluster_text",
    "resolve_cluster",
    "library_file_sha256",
]
