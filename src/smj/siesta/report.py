"""Reports for the six-stage SIESTA route, using collected evidence only."""
from __future__ import annotations

from pathlib import Path

from ..errors import StageError
from ..state import RunLayout, now_iso
from . import get_stage
from .io import read_json
from .workflow import check_outputs, stage_protocol


def _kconv(run_dir):
    path = Path(run_dir) / "kconv.json"
    if not path.is_file():
        return None
    data = read_json(path)
    if data.get("schema") != "smj.siesta-kconv.v1":
        raise StageError("kconv.json 的 schema 不受支持")
    if not isinstance(data.get("passed"), bool):
        raise StageError("kconv.json 缺少 passed 判定")
    return data


def _physics_verdict(scf, kconv):
    screening = (scf or {}).get("screening")
    if not isinstance(screening, dict):
        screening = {}
    if screening.get("status") == "failed" or (kconv and kconv.get("passed") is False):
        return "failed"
    if screening.get("status") == "passed" and kconv and kconv.get("passed") is True:
        return "passed"
    return "inconclusive"


def _verdict_reason(scf, kconv):
    screening = (scf or {}).get("screening")
    status = screening.get("status") if isinstance(screening, dict) else None
    if status == "failed":
        return "屏蔽判据未通过"
    if kconv and kconv.get("passed") is False:
        return "横向 k 点相对变化超过容差"
    if status != "passed":
        return "屏蔽证据缺失"
    return "横向 k 点收敛证据尚未生成"


def evidence(run_dir, state, config=None):
    layout = RunLayout(run_dir)
    issues = {}
    detail = {}
    contacts = scf = result = None
    for name in state.order:
        stage, record = get_stage(name), state.stage(name)
        detail[name] = {"status": record.status, "job_id": record.job_id, "error": record.error,
                        "outputs_sha256": record.outputs_sha256, "protocol_sha256": record.extra.get("siesta_protocol")}
        try:
            if record.status != "collected":
                raise StageError(f"状态为 {record.status}" + (f"，{record.error}" if record.error else ""))
            recorded_parents = record.extra.get("siesta_parents") or {}
            failed_parents = [p for p in recorded_parents if p in issues]
            if failed_parents:
                raise StageError("上游待核实 " + ", ".join(failed_parents))
            if recorded_parents != {p: state.stage(p).outputs_sha256 for p in recorded_parents}:
                raise StageError("上游产物摘要已变化")
            protocol = record.extra.get("siesta_protocol")
            if not protocol or (config and protocol != stage_protocol(config, name)):
                raise StageError("运行协议与本阶段记录不符")
            paths = {Path(file).name: layout.artifact(name, file) for file in record.outputs_sha256}
            if not set(stage.outputs).issubset(paths):
                raise StageError("标准产物缺少采集摘要")
            check_outputs(layout, name, record)
            if name == "siesta_junction_relax":
                contacts = read_json(paths["contacts.json"])
                if contacts.get("passed") is not True:
                    raise StageError("最终接触或电极匹配未通过")
            elif name == "transiesta_scf":
                scf = read_json(paths["scf.json"])
                if scf.get("open_boundary_converged") is not True:
                    raise StageError("缺少开放边界自洽收敛记录")
            elif name == "junction_analyze":
                saved = read_json(paths["results.json"])
                if saved.get("backend") != "transiesta" or saved.get("contacts") != contacts:
                    raise StageError("输运结果与本次结构或路线不符")
                result = saved
        except (StageError, OSError, ValueError) as exc:
            issues[name] = str(exc)
    kconv = None
    try:
        kconv = _kconv(run_dir)
    except (StageError, OSError, ValueError) as exc:
        issues["kconv"] = str(exc)
    return {"stages": detail, "issues": issues, "contacts": contacts, "scf": scf,
            "transport": result, "kconv": kconv,
            "physics_verdict": _physics_verdict(scf, kconv)}


def manifest(run_dir, state, config, rid, spec, run_spec, provenance, include_checks, warnings):
    audit = evidence(run_dir, state, config)
    issues = audit["issues"]
    verdict = {"engineering_passed": not issues, "physics_verdict": audit["physics_verdict"],
               "errors": len(issues), "warnings": len((audit["transport"] or {}).get("warnings", []))}
    return {"schema": "smj.report.v1", "backend": "transiesta", "run_id": rid,
            "root": str(run_dir), "generated_at": now_iso(), "run_spec": run_spec, "spec": spec,
            "provenance": provenance, "state": {"stages": state.statuses(), "is_done": not issues,
                                                "first_incomplete": next(iter(issues), None)},
            "stages": audit["stages"], "stage_issues": issues, "contacts": audit["contacts"],
            "scf": audit["scf"], "transport": audit["transport"], "kconv": audit["kconv"],
            "warnings": warnings, "verdict": verdict if include_checks else None}


def _screening_text(screening):
    if not isinstance(screening, dict):
        return str(screening or "unverified")
    parts = [f"状态 {screening.get('status', 'unverified')}"]
    for side, values in (screening.get("sides") or {}).items():
        outer = values.get("outermost_delta_e")
        second = values.get("second_outermost_delta_e")
        if outer is not None and second is not None:
            parts.append(f"{side} 最外两层 Δq={outer:+.2e}/{second:+.2e} e")
    spread = screening.get("electrode_layer_spread_e")
    if spread is not None:
        parts.append(f"电极层间 spread={spread:.2e} e")
    if screening.get("reason"):
        parts.append(f"原因：{screening['reason']}")
    if screening.get("criterion"):
        parts.append(f"判据：{screening['criterion']}")
    return "；".join(parts)


def render(data):
    lines = [f"# {data['run_id']} 零偏压电子输运报告", "",
             "计算路线为 SIESTA → TranSIESTA → TBtrans。", "",
             "## 阶段与采集证据", "", "| 阶段 | 状态 | 作业号 | 复核 |", "|---|---|---|---|"]
    for name, entry in data["stages"].items():
        issue = data["stage_issues"].get(name, "已核对")
        lines.append(f"| {name} | {entry['status']} | {entry['job_id'] or '—'} | {issue} |")
    contacts = data["contacts"] or {}
    if contacts.get("contacts"):
        lines.extend(["", "## 最终接触结构", "",
                      "| 电极 | 形貌 | 初始位点 | 最终位点 | 锚点编号 | 配位数 |",
                      "|---|---|---|---|---|---|"])
        for side, contact in contacts["contacts"].items():
            values = [contact.get(k, "—") for k in ("motif", "initial_site", "final_site",
                                                     "anchor_id", "coordination")]
            lines.append("| " + " | ".join(map(str, [side, *values])) + " |")
    scf = data.get("scf")
    if scf:
        lines.extend(["", f"开放边界 SCF 电荷偏差为 {scf.get('charge_deviation_e', '—')} e。"
                      f"电极边界屏蔽证据为 {_screening_text(scf.get('screening'))}。"])
    lines.extend(["", "## 线性电导", ""])
    result = data["transport"]
    if result:
        lines.extend(["| 量 | 值 |", "|---|---|",
                      f"| G/G₀ | {result['G_over_G0']:.10g} |",
                      f"| G (S) | {result['G_S']:.10g} |",
                      f"| 分析温度 (K) | {result['temperature_K']} |",
                      f"| 占据展宽 (eV) | {result['smearing_eV']} |",
                      f"| 共同费米能级 (eV) | {result['fermi_eV']} |",
                      f"| 自旋约定 | {result['spin_convention']} |",
                      ""])
        entry = data["stages"]["junction_analyze"]
        image = next((Path("junction_analyze") / p for p in entry["outputs_sha256"]
                      if Path(p).name == "T_E.png"), None)
        if image:
            lines.append(f"![相对费米能级的透射谱](<{image.as_posix()}>)")
        lines.extend(f"\n积分提示：{warning}" for warning in result.get("warnings", []))
    else:
        lines.append("本次计算缺少完整、可核实的阶段记录。")
    kconv = data.get("kconv")
    if kconv:
        left, right = kconv["base"], kconv["refined"]
        metric = kconv["metrics"]["G_over_G0"]
        relative = metric.get("relative_change")
        relative_text = "n/a" if relative is None else f"{100.0 * relative:+.3f}%"
        lines.extend(["", "## 横向 k 点收敛", "",
                      "| 网格 | G/G₀ | E_F (eV) |", "|---|---:|---:|",
                      f"| {left.get('name', left['run'])} | {left['G_over_G0']:.10g} | {left['fermi_eV']:.9f} |",
                      f"| {right.get('name', right['run'])} | {right['G_over_G0']:.10g} | {right['fermi_eV']:.9f} |",
                      "",
                      f"相对变化 {relative_text}，容差 {100.0 * kconv['tolerance']:.2f}%；"
                      f"判定 **{'passed' if kconv['passed'] else 'failed'}**。"])
    if data.get("verdict"):
        verdict = data["verdict"]["physics_verdict"]
        if verdict == "passed":
            lines.extend(["", "数值判据（横向 k 点与屏蔽）通过；这不构成实验有效性证据。"])
        elif verdict == "failed":
            reason = _verdict_reason(data.get("scf"), data.get("kconv"))
            lines.extend(["", f"物理判定为 **failed**：{reason}。"])
        else:
            reason = _verdict_reason(data.get("scf"), data.get("kconv"))
            lines.extend(["", f"物理判定为 **inconclusive**：{reason}。"])
    lines.extend(["", "## 产物与来源", ""])
    for stage, entry in data["stages"].items():
        for file in entry["outputs_sha256"]:
            path = (Path(stage) / file).as_posix()
            lines.append(f"- [{path}](<{path}>)")
    return "\n".join(lines) + "\n"
