"""把多个 run 的 ``thermal.json`` 摆到一张表里（``smj thermal-compare``）。

这个模块存在的唯一理由，是**防止拿不同方法算出来的数互比**。``G_ph`` 对本仓
自己的正则化参数（D15 冻结的 η 与能量步长）、温度、以及电极力常数极为敏感：

* η 从 0.5 meV 降到 0.05 meV，nap-sh 的 ``G_ph`` 从 10.9 涨到 13.7 pW/K；
* 换一份 lead 力常数会改变带顶附近的透射权重。

所以本模块把「协议」当成和一阶结果同等级的产物：先比对 ``temperature_K`` /
``eta_meV`` / ``energy_step_meV`` / ``band_top_meV`` / ``lead_fc_sha256``，
不一致就显式列出来，由调用方决定是警告还是拒绝。
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

from ..analysis import read_thermal
from .environment import IDENTITY_KEYS, identity

__all__ = [
    "PROTOCOL_KEYS",
    "IDENTITY_KEYS",
    "RunRecord",
    "Comparison",
    "find_runs",
    "load_run",
    "compare",
    "render_markdown",
]

#: 必须一致才能横向比较的字段；``temperature_K`` 允许通过 normalize 对齐。
PROTOCOL_KEYS = ("eta_meV", "energy_step_meV", "band_top_meV", "lead_fc_sha256")


@dataclass
class RunRecord:
    """一个 run 的热导结果与其协议。"""

    name: str
    path: Path
    payload: dict
    verdict: str | None = None
    failed: list[str] = field(default_factory=list)
    inconclusive: list[str] = field(default_factory=list)
    environment: dict | None = None

    @property
    def g_ph(self) -> float:
        return float(self.payload["g_ph_pw_per_k"])

    @property
    def protocol(self) -> dict:
        return {key: self.payload.get(key) for key in PROTOCOL_KEYS}

    @property
    def environment_identity(self) -> dict:
        return identity(self.environment)


@dataclass
class Comparison:
    """多 run 的比较结果。"""

    records: list[RunRecord]
    mismatches: dict[str, list[str]]
    consistent: bool
    environment_mismatches: dict[str, list[str]] = field(default_factory=dict)

    @property
    def environment_consistent(self) -> bool:
        return not self.environment_mismatches

    def rows(self) -> list[dict]:
        return [
            {
                "name": record.name,
                "g_ph_pw_per_k": record.g_ph,
                "n_channels_effective": record.payload.get("n_channels_effective"),
                "g_quantum_pw_per_k": record.payload.get("g_quantum_pw_per_k"),
                "temperature_K": record.payload.get("temperature_K"),
                "eta_meV": record.payload.get("eta_meV"),
                "energy_step_meV": record.payload.get("energy_step_meV"),
                "n_imaginary": (record.payload.get("imaginary_modes") or {}).get(
                    "count"),
                "verdict": record.verdict,
            }
            for record in self.records
        ]


def find_runs(paths) -> list[Path]:
    """把用户给的路径展开成「含 ``thermal_analyze/thermal.json`` 的目录」。

    接受三种写法：run 目录本身、run 目录的父目录（扫一层）、直接指向
    ``thermal.json``。父目录只扫一层，因为本仓的 run 目录布局是
    ``runs/<name>/<name>/thermal_analyze``（``smj init`` 会多套一层同名目录）。
    """
    found: list[Path] = []
    for raw in paths:
        path = Path(raw)
        if path.is_file() and path.name == "thermal.json":
            found.append(path.parent.parent)
            continue
        if (path / "thermal_analyze" / "thermal.json").is_file():
            found.append(path)
            continue
        for child in sorted(path.iterdir()) if path.is_dir() else []:
            if (child / "thermal_analyze" / "thermal.json").is_file():
                found.append(child)
            elif (child / child.name / "thermal_analyze" / "thermal.json").is_file():
                found.append(child / child.name)
    seen: set[Path] = set()
    unique: list[Path] = []
    for item in found:
        resolved = item.resolve()
        if resolved not in seen:
            seen.add(resolved)
            unique.append(item)
    if not unique:
        raise ValueError(
            "没有找到任何 thermal_analyze/thermal.json；请给 run 目录或其父目录")
    return unique


def load_run(path: Path, *, name: str | None = None) -> RunRecord:
    """读一个 run 的 ``thermal.json`` 与（若存在）``review.json`` 的热导判定。"""
    path = Path(path)
    thermal_path = path / "thermal_analyze" / "thermal.json"
    payload, _warnings = read_thermal(thermal_path)
    record = RunRecord(name=name or path.name, path=path, payload=payload)
    environment_path = path / "thermal_analyze" / "environment.json"
    if environment_path.is_file():
        try:
            record.environment = json.loads(
                environment_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            record.environment = None
    review_path = path / "review.json"
    if review_path.is_file():
        try:
            review = json.loads(review_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            review = {}
        thermal = review.get("thermal") or {}
        record.verdict = thermal.get("verdict")
        record.failed = list(thermal.get("failed") or [])
        record.inconclusive = list(thermal.get("inconclusive") or [])
    return record


def compare(records: list[RunRecord]) -> Comparison:
    """找出协议不一致的字段（每个字段列出「run 名: 取值」的差异清单）。"""
    mismatches: dict[str, list[str]] = {}
    for key in ("temperature_K",) + PROTOCOL_KEYS:
        values = {}
        for record in records:
            value = record.payload.get(key)
            values.setdefault(repr(value), []).append(record.name)
        if len(values) > 1:
            mismatches[key] = ["%s: %s" % (", ".join(names), raw)
                               for raw, names in values.items()]

    # 环境单独比：它不属于"协议"，因为环境差异往往不是选择而是意外
    # （换了 numpy、换了 aims 二进制）。但差异同样能让 G_ph 变几个百分点，
    # 所以也要摆出来。单 run 时不比 —— 没有第二个东西可对照。
    environment_mismatches: dict[str, list[str]] = {}
    if len(records) > 1:
        for key in IDENTITY_KEYS:
            values: dict[str, list[str]] = {}
            for record in records:
                values.setdefault(repr(record.environment_identity.get(key)),
                                  []).append(record.name)
            if len(values) > 1:
                environment_mismatches[key] = [
                    "%s: %s" % (", ".join(names), raw)
                    for raw, names in values.items()]

    return Comparison(records=records, mismatches=mismatches,
                      consistent=not mismatches,
                      environment_mismatches=environment_mismatches)


def render_markdown(result: Comparison) -> str:
    """把比较结果渲染成 markdown：协议警告在前，结果表在后。"""
    lines = ["## 声子热导比较", ""]
    if not result.consistent:
        lines.append("> **协议不一致**：下面这些字段在不同 run 之间取值不同。"
                     "跨协议的 `G_ph` 不能直接横向比较（见方案 D15），"
                     "请先对齐 η、能量步长、温度与电极力常数。")
        lines.append("")
        for key, detail in result.mismatches.items():
            lines.append("- `%s`：%s" % (key, "；".join(detail)))
        lines.append("")
    if result.environment_mismatches:
        lines.append("> **环境不一致**：决定数值的输入或工具在 run 之间不同。"
                     "这不是协议选择，通常是换了机器、换了 FHI-aims 二进制、"
                     "species 目录或 numpy 版本 —— 差异可能来自这里而不是体系。")
        lines.append("")
        for key, detail in result.environment_mismatches.items():
            lines.append("- `%s`：%s" % (key, "；".join(detail)))
        lines.append("")
    lines.append("| run | G_ph (pW/K) | N_eff | G_Q (pW/K) | T (K) | "
                 "η (meV) | 步长 (meV) | 虚频 | 判定 |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for row in result.rows():
        lines.append(
            "| %s | %.3f | %s | %s | %s | %s | %s | %s | %s |" % (
                row["name"],
                row["g_ph_pw_per_k"],
                "—" if row["n_channels_effective"] is None
                else "%.4f" % row["n_channels_effective"],
                "—" if row["g_quantum_pw_per_k"] is None
                else "%.2f" % row["g_quantum_pw_per_k"],
                "—" if row["temperature_K"] is None else "%.0f" % row["temperature_K"],
                "—" if row["eta_meV"] is None else "%.3f" % row["eta_meV"],
                "—" if row["energy_step_meV"] is None
                else "%.3f" % row["energy_step_meV"],
                "—" if row["n_imaginary"] is None else str(row["n_imaginary"]),
                row["verdict"] or "—",
            ))
    lines.append("")
    return "\n".join(lines)
