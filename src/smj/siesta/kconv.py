"""Transverse k-point convergence for TranSIESTA junction runs.

The comparison is deliberately narrow: it reads the *collected* transport
result of two runs (``junction_analyze/results.json``) and asks whether the
linear-response conductance changed by less than a caller-supplied tolerance.
It does not pretend that agreement between two k meshes is a proof of
convergence; the report keeps that scope statement.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

SCHEMA = "smj.siesta-kconv.v1"
DEFAULT_TOLERANCE = 0.01


def _results_path(path: str | Path) -> Path:
    candidate = Path(path)
    if candidate.is_file():
        if candidate.name != "results.json":
            raise ValueError("expected junction_analyze/results.json, got %s" % candidate)
        return candidate
    for child in (candidate / "junction_analyze" / "results.json",
                  candidate / "results.json"):
        if child.is_file():
            return child
    raise ValueError("no junction_analyze/results.json under %s" % candidate)


def load_run(path: str | Path) -> dict[str, Any]:
    """Load the transport scalars needed for a k-mesh comparison."""
    source = _results_path(path)
    data = json.loads(source.read_text(encoding="utf-8"))
    missing = [key for key in ("G_over_G0", "G_S", "fermi_eV") if key not in data]
    if missing:
        raise ValueError("%s is missing %s" % (source, ", ".join(missing)))
    return {
        "path": str(source),
        "run": str(source.parent.parent),
        "name": source.parent.parent.name,
        "G_over_G0": float(data["G_over_G0"]),
        "G_S": float(data["G_S"]),
        "fermi_eV": float(data["fermi_eV"]),
        "physical_verdict": data.get("physical_verdict"),
    }


def _relative_change(base: float, refined: float) -> float | None:
    if base == 0.0:
        return None if refined == 0.0 else float("inf")
    return (refined - base) / base


def compare(base: str | Path, refined: str | Path,
            tolerance: float = DEFAULT_TOLERANCE) -> dict[str, Any]:
    """Compare two collected runs; ``passed`` means |ΔG/G| <= tolerance."""
    if tolerance <= 0:
        raise ValueError("tolerance must be positive")
    left = load_run(base)
    right = load_run(refined)
    relative = _relative_change(left["G_over_G0"], right["G_over_G0"])
    passed = relative is not None and abs(relative) <= tolerance
    return {
        "schema": SCHEMA,
        "base": left,
        "refined": right,
        "tolerance": tolerance,
        "metrics": {
            "G_over_G0": {
                "base": left["G_over_G0"],
                "refined": right["G_over_G0"],
                "relative_change": relative,
                "passed": passed,
            },
            "fermi_eV": {
                "base": left["fermi_eV"],
                "refined": right["fermi_eV"],
                "delta_eV": right["fermi_eV"] - left["fermi_eV"],
            },
        },
        "passed": passed,
        "scope": (
            "two transverse k meshes agree within tolerance; this is a "
            "numerical convergence check, not evidence of experimental validity"
            if passed else
            "two transverse k meshes differ beyond tolerance; this is a "
            "numerical convergence check, not evidence of experimental validity"
        ),
    }


def render_markdown(payload: Mapping[str, Any]) -> str:
    """Render a compact Markdown table for the validation record."""
    left = payload["base"]
    right = payload["refined"]
    metric = payload["metrics"]["G_over_G0"]
    rel = metric["relative_change"]
    rel_text = "n/a" if rel is None else "%+.3f%%" % (100.0 * rel)
    delta_ef = payload["metrics"]["fermi_eV"]["delta_eV"]
    verdict = "passed" if payload["passed"] else "failed"
    return "\n".join([
        "| k 网格对比 | G/G0 | G/S | E_F/eV |",
        "| --- | ---: | ---: | ---: |",
        "| %s | %.6e | %.6e | %.9f |" % (
            left.get("name", left["run"]), left["G_over_G0"], left["G_S"], left["fermi_eV"]),
        "| %s | %.6e | %.6e | %.9f |" % (
            right.get("name", right["run"]), right["G_over_G0"], right["G_S"], right["fermi_eV"]),
        "",
        "相对变化：%s；容差：%.2f%%；判定：**%s**。" % (
            rel_text, 100.0 * payload["tolerance"], verdict),
        "E_F 偏移：%+.6f eV。" % delta_ef,
        "",
        "> %s" % payload["scope"],
    ])
