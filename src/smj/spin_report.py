"""``smj spin``：一份 run 的开壳层自旋汇总（2026-11-21 批 D / S8）。

把散在四处的自旋信息压成一份可读表格或 JSON：

- ``run_spec.json`` 的 ``config.molecule``：多重度、电荷、是否钉住总自旋、BS；
- ``aitranss_te/TE.dat``：α/β 双通道 → P(E_F)、P(L0)、Sα/Sβ/S_total、P_S
  （旧 run 的 ``results.json`` 里的 spin 块可能没有分通道 S，直接从谱重算）；
- ``aims_scf/spin_partition.json``：Hirshfeld 分区（分子/接触区/电极）与
  ``localized_antiparallel``；
- ``review.json`` 的自旋相关 code（不重跑复核，只汇总已有结论）。

只读、不写任何文件；缺哪块就明说缺哪块，不补零、不猜测。
"""
from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from .landauer_spin import analyze_te_spin
from .spin_summary import (
    FIXED_CONSTRAINT_CAVEAT,
    SINGLE_CHANNEL_CAVEAT,
    SPIN_CAVEAT,
    findings as spin_findings,
    normalize,
    percent,
)

_SPIN_CODE_PREFIXES = ("open_shell", "spin_", "broken_symmetry")

_GROUP_ORDER = ("molecule", "contact", "electrode")


def _load_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _number(value: Any, fmt: str = "%.6f") -> str:
    try:
        return fmt % float(value)
    except (TypeError, ValueError):
        return "—"


def _te_path(run: Path) -> Path | None:
    for candidate in (run / "aitranss_te" / "TE.dat", run / "TE.dat"):
        if candidate.is_file():
            return candidate
    return None


def _molecule_block(spec: Mapping[str, Any]) -> dict[str, Any]:
    config = spec.get("config") or {}
    molecule = config.get("molecule") or {}
    return {
        "multiplicity": molecule.get("multiplicity"),
        "charge": molecule.get("charge"),
        "fixed_spin_moment": molecule.get("fixed_spin_moment", True),
        "broken_symmetry": bool(molecule.get("broken_symmetry", False)),
        "family": molecule.get("family"),
    }


def _temperature(spec: Mapping[str, Any], spin: Mapping[str, Any] | None,
                 override: float | None) -> float:
    if override is not None:
        return float(override)
    config = spec.get("config") or {}
    transport = config.get("transport") or {}
    for value in (transport.get("temperature"), (spin or {}).get("temperature_K")):
        try:
            if value is not None:
                return float(value)
        except (TypeError, ValueError):
            continue
    return 300.0


def _partition(run: Path) -> dict[str, Any] | None:
    payload = _load_json(run / "aims_scf" / "spin_partition.json")
    if not isinstance(payload, Mapping) or not isinstance(payload.get("groups"), Mapping):
        return None
    return dict(payload)


def _review_codes(run: Path) -> list[dict[str, Any]]:
    review = _load_json(run / "review.json")
    if not isinstance(review, Mapping):
        return []
    codes: list[dict[str, Any]] = []
    for issue in review.get("issues") or ():
        if not isinstance(issue, Mapping):
            continue
        code = str(issue.get("code") or "")
        if code.startswith(_SPIN_CODE_PREFIXES):
            codes.append({"severity": issue.get("severity"),
                          "code": code, "message": issue.get("message")})
    return codes


def build_spin_report(run_dir: str | Path, *, temperature: float | None = None) -> dict:
    """汇总一份 run 的自旋信息；run 目录不存在时抛 FileNotFoundError。"""
    run = Path(run_dir)
    if not run.is_dir() or not ((run / "run_spec.json").is_file()
                                or (run / "state.json").is_file()
                                or _te_path(run) is not None):
        raise FileNotFoundError(
            "找不到 run 目录（需要 run_spec.json / state.json / TE.dat 至少一项）: %s"
            % run)

    spec = _load_json(run / "run_spec.json")
    spec = spec if isinstance(spec, Mapping) else {}
    molecule = _molecule_block(spec)
    open_shell = bool(molecule["broken_symmetry"]) or bool(
        molecule["multiplicity"] and int(molecule["multiplicity"]) > 1)

    te = _te_path(run)
    spin_source = None
    raw_spin = None
    if te is not None:
        try:
            raw_spin = analyze_te_spin(
                te, temperature=_temperature(spec, None, temperature))
            spin_source = "TE.dat"
        except Exception:  # 谱坏了就退回 results.json 的 spin 块
            raw_spin = None
    if raw_spin is None:
        results = _load_json(run / "junction_analyze" / "results.json")
        if isinstance(results, Mapping) and isinstance(results.get("spin"), Mapping):
            raw_spin = dict(results["spin"])
            spin_source = "junction_analyze/results.json"

    temp = _temperature(spec, raw_spin, temperature)
    summary = normalize(
        raw_spin,
        multiplicity=molecule["multiplicity"],
        charge=molecule["charge"],
        fixed_spin_moment=molecule["fixed_spin_moment"],
        broken_symmetry=molecule["broken_symmetry"],
    )
    if summary is not None and raw_spin is not None:
        summary.setdefault("temperature_K", raw_spin.get("temperature_K"))

    partition = _partition(run)
    findings = spin_findings(
        summary,
        multiplicity=molecule["multiplicity"],
        open_shell=open_shell,
        te_present=te is not None,
        broken_symmetry=molecule["broken_symmetry"],
    )
    caveats: list[str] = []
    if open_shell:
        caveats.append(SPIN_CAVEAT)
        if molecule["fixed_spin_moment"]:
            caveats.append(FIXED_CONSTRAINT_CAVEAT)
        if summary is not None and summary.get("channels") == 1:
            caveats.append(SINGLE_CHANNEL_CAVEAT)

    run_id = spec.get("run_id") or run.name
    return {
        "run": str(run),
        "run_id": run_id,
        "molecule": molecule,
        "open_shell": open_shell,
        "temperature_K": temp,
        "spin_source": spin_source,
        "spin": summary,
        "partition": partition,
        "findings": findings,
        "review_codes": _review_codes(run),
        "caveats": caveats,
    }


def _spin_lines(payload: Mapping[str, Any]) -> list[str]:
    spin = payload.get("spin")
    if not payload.get("open_shell"):
        return ["自旋: 闭壳层（multiplicity=1）；无 α/β 通道。"]
    if not isinstance(spin, Mapping):
        return ["自旋: 开壳层，但没有可用的 TE.dat / spin 块。"]
    lines = ["自旋: %s 通道 / %s 点 / %s K（来源 %s）" % (
        spin.get("channels"), spin.get("points"),
        _number(payload.get("temperature_K"), "%.1f"), payload.get("spin_source"))]
    lines.append("Tα(E_F) / Tβ(E_F): %s / %s" % (
        _number(spin.get("T_alpha_at_ef")), _number(spin.get("T_beta_at_ef"))))
    lines.append("P(E_F) / P(L0): %s / %s" % (
        percent(spin.get("spin_polarization_at_ef")),
        percent(spin.get("spin_polarization_L0"))))
    lines.append("Gα / Gβ / Gtotal (G0): %s / %s / %s" % (
        _number(spin.get("g_over_g0_alpha")), _number(spin.get("g_over_g0_beta")),
        _number(spin.get("g_over_g0_total"))))
    lines.append("Sα / Sβ / S_total [μV/K]: %s / %s / %s" % (
        _number(spin.get("S_alpha_uV_per_K"), "%+.3f"),
        _number(spin.get("S_beta_uV_per_K"), "%+.3f"),
        _number(spin.get("S_total_uV_per_K"), "%+.3f")))
    lines.append("P_S = (Sα−Sβ)/(Sα+Sβ): %s" % percent(spin.get("spin_polarization_S")))
    if spin.get("note"):
        lines.append("备注: %s" % spin["note"])
    return lines


def _partition_lines(partition: Mapping[str, Any] | None) -> list[str]:
    if not isinstance(partition, Mapping):
        return ["自旋分布: 未采集（aims_scf/spin_partition.json 不存在；"
                "开壳层重跑 aims_scf 会自动生成）。"]
    total = partition.get("total")
    head = "自旋分布 (Hirshfeld, μB): 总 %s" % _number(total, "%.4f")
    groups = partition.get("groups") or {}
    parts = []
    for name in _GROUP_ORDER:
        group = groups.get(name)
        if not isinstance(group, Mapping):
            continue
        parts.append("%s %s (%s)" % (
            {"molecule": "分子", "contact": "接触区", "electrode": "电极"}[name],
            _number(group.get("spin_moment"), "%+.4f"),
            percent(group.get("fraction"))))
    if parts:
        head += "；" + "；".join(parts)
    lines = [head]
    if partition.get("molecule_n_positive") is not None:
        verdict = partition.get("localized_antiparallel")
        tag = ("局部反平行（正负共存且净矩很小）" if verdict is True
               else "未满足局部反平行判据" if verdict is False
               else "分区数据不完整")
        lines.append("分子内 +/− 自旋矩: %s / %s；净矩 %s μB（%s）" % (
            partition.get("molecule_n_positive"),
            partition.get("molecule_n_negative"),
            _number(partition.get("molecule_net"), "%+.4f"), tag))
    return lines


def render_spin_report(payload: Mapping[str, Any]) -> str:
    """把 :func:`build_spin_report` 的结果渲染成人类可读文本。"""
    molecule = payload.get("molecule") or {}
    lines = ["自旋汇总: %s" % payload.get("run_id"),
             "目录: %s" % payload.get("run"),
             "分子: multiplicity=%s charge=%s fixed_spin_moment=%s "
             "broken_symmetry=%s" % (
                 molecule.get("multiplicity"), molecule.get("charge"),
                 molecule.get("fixed_spin_moment"),
                 molecule.get("broken_symmetry"))]
    lines.extend(_spin_lines(payload))
    lines.extend(_partition_lines(payload.get("partition")))

    review_codes = payload.get("review_codes") or []
    findings = payload.get("findings") or []
    if review_codes:
        lines.append("review 已有自旋提示:")
        for item in review_codes:
            lines.append("  - [%s] %s: %s" % (
                item.get("severity"), item.get("code"), item.get("message")))
    if findings:
        lines.append("本地复核提示（未写回 review.json）:")
        for item in findings:
            lines.append("  - [%s] %s: %s" % (
                item.get("severity"), item.get("code"), item.get("message")))
    caveats = payload.get("caveats") or []
    if caveats:
        lines.append("口径提醒:")
        lines.extend("  - %s" % text for text in caveats)
    return "\n".join(lines)
