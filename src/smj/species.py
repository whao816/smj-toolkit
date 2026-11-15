"""FHI-aims species 目录审计：实际用了哪套基组，而不是配置里的标签。

`molecule.opt_basis` 只是写在配置里的标签（`light` / `tight` …），真正生效的是
`cluster.species_dir`（或 `paths.species_defaults`）里逐元素的 species 文件。两者
可以不一致，且 species 文件本身决定了径向基函数、角动量与积分网格——所以每次
`molecule_opt` 都应把"实际用了什么"记进 `species.json` 与 `results.json`。

本模块只读 species 文件，不启动任何计算；它给的是**口径可核对性**，不是敏感性结果。
真正的基组敏感性要拿同一分子分别用 light / tight 目录各跑一次再比（见 docs）。
"""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

from ._vendor.junction.tcontrol import find_species_file

#: 文件头里的基组档位：`Suggested "light" defaults for S atom`。
_TIER_RE = re.compile(r'Suggested\s+"([A-Za-z][A-Za-z0-9_-]*)"\s+defaults', re.IGNORECASE)
#: 生效的基组行：`  hydro 2 p 1.8` / `  ionic 3 d auto`（注释行 `#     hydro …` 不算）。
_BASIS_RE = re.compile(r"^\s*(hydro|ionic)\s+(\d+)\s+([spdfgh])\b", re.IGNORECASE)
#: 目录名里可能出现的档位词（FHI-aims 发行版通常按 light/tight 分目录）。
_DIR_TIERS = ("light", "tight", "intermediate", "really_tight", "very_tight")

#: 角动量序号 → 字母，用于报告"最高角动量"。
_ANGULAR = ("s", "p", "d", "f", "g", "h")


@dataclass(frozen=True)
class SpeciesInfo:
    """单个元素的 species 文件审计结果。"""

    element: str
    path: str
    sha256: str
    tier: str            # 文件头声明；没有则 "unknown"
    basis_functions: int  # 生效的基组行数
    max_angular: str     # 最高角动量（s/p/d/f/g/h）；没有基组行则 ""

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_species(text: str) -> dict[str, Any]:
    """从 species 文件文本里取档位与基组统计。"""
    lines = text.splitlines()
    tier = "unknown"
    for line in lines:
        match = _TIER_RE.search(line)
        if match:
            tier = match.group(1).lower()
            break
    functions = 0
    highest = -1
    for line in lines:
        stripped = line.lstrip()
        if not stripped or stripped.startswith("#"):
            continue
        match = _BASIS_RE.match(line)
        if not match:
            continue
        functions += 1
        highest = max(highest, _ANGULAR.index(match.group(3).lower()))
    return {
        "tier": tier,
        "basis_functions": functions,
        "max_angular": _ANGULAR[highest] if highest >= 0 else "",
    }


def species_info(species_dir: str | Path, element: str) -> SpeciesInfo:
    """审计一个元素：文件位置、inode 级内容指纹、档位、基组规模。"""
    path = find_species_file(species_dir, element)
    parsed = parse_species(path.read_text(encoding="utf-8", errors="ignore"))
    return SpeciesInfo(element=element, path=str(path), sha256=_sha256(path), **parsed)


def species_dir_tier(species_dir: str | Path) -> str:
    """从目录路径猜档位（不依赖文件头；找不到返回空串）。

    只在文件头缺失时作为补充线索：`.../species_defaults/light` → `light`。
    """
    parts = [part.lower() for part in Path(species_dir).parts]
    for part in reversed(parts):
        if part in _DIR_TIERS:
            return part
    return ""


def audit_species_dir(species_dir: str | Path,
                      elements: Iterable[str] | None = None) -> list[SpeciesInfo]:
    """审计一组元素（缺省：目录里所有 `*_default` 文件）。"""
    root = Path(species_dir)
    if elements is None:
        names = [path.name for path in root.glob("*_default")]
        elements = sorted({name.rsplit("_", 2)[-2] for name in names})
    return [species_info(root, str(element)) for element in elements]


def elements_of(atoms: Sequence[Any]) -> list[str]:
    """按首次出现顺序去重元素符号（以几何里的顺序为准）。"""
    seen: list[str] = []
    for atom in atoms:
        symbol = str(getattr(atom, "element", atom))
        if symbol and symbol not in seen:
            seen.append(symbol)
    return seen


def summary(species_dir: str | Path, atoms: Sequence[Any] | None = None,
            elements: Iterable[str] | None = None, *,
            label: str | None = None) -> dict[str, Any]:
    """生成写进 `species.json` 的摘要，并核对配置标签与实际档位是否一致。"""
    if elements is not None:
        wanted = [str(item) for item in elements]
    elif atoms:
        wanted = elements_of(atoms)
    else:
        wanted = [info.element for info in audit_species_dir(species_dir)]
    infos = [species_info(species_dir, element) for element in wanted]
    tiers = sorted({info.tier for info in infos if info.tier != "unknown"})
    dir_tier = species_dir_tier(species_dir)
    actual = tiers[0] if len(tiers) == 1 else ("mixed" if tiers else dir_tier)
    warnings: list[str] = []
    if len(tiers) > 1:
        warnings.append("species 目录里混用了多个档位 %s：基组口径不统一"
                        % ", ".join(tiers))
    if label:
        expected = str(label).strip().lower()
        if actual and expected and actual != expected and actual not in ("mixed", ""):
            warnings.append(
                "molecule.opt_basis=%s 只是配置标签，实际 species 档位是 %s"
                "（%s）；以 species 文件为准"
                % (expected, actual, Path(species_dir).name))
    return {
        "species_dir": str(species_dir),
        "declared_label": label or "",
        "actual_tier": actual or "unknown",
        "dir_tier": dir_tier,
        "elements": [info.as_dict() for info in infos],
        "warnings": warnings,
    }


def write_species_json(payload: dict[str, Any], path: str | Path) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
                      encoding="utf-8", newline="\n")
    return target


__all__ = [
    "SpeciesInfo",
    "parse_species",
    "species_info",
    "species_dir_tier",
    "audit_species_dir",
    "elements_of",
    "summary",
    "write_species_json",
]
