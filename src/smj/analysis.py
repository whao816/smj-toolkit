"""TE.dat 输运分析：Landauer L0/L1/L2 与 G / S / k_e / ZT 摘要。

本模块是 :mod:`smj.junction.landauer`（以及可选的
:mod:`smj.junction.independent_check`）的**薄封装**，只做三件事：

1. 调 vendored ``analyze_te`` 得到固定字段的统计摘要；
2. 可选跑独立复核（独立解析器 + Simpson 求积 + CODATA 常数），抓实现级单位/
   口径事故（见 vendored ``independent_check`` 的独立性边界声明）；
3. 落 ``analysis.json``，可选出图（``matplotlib`` 属 extra ``plot``，未装时
   优雅跳过并说明），可选与参考值 compare（默认容差：log10(G/G0) 绝对 0.1、
   S 相对 10%）。

口径提醒
--------
aitranss 的 ``T(E)`` 以 **G0 = 2e^2/h 为单位（per spin）**，非极化时单列、
极化时双列（vendored 解析器按通道取均值）。因此 ``L0`` 的数值直接就是
``G/G0``，而 ``log10_G_over_G0_like`` 里的 ``_like`` 表示命名是项目近似，
不得把示例数据的值当生产级科学结论（7 点样本免责见 vendored landauer 头注）。

所有 vendored 模块都在函数内 **lazy import**：``import smj.analysis`` 不会
触发 vendored 依赖，也就不会因缺可选依赖而崩。
"""

from __future__ import annotations

import json
import math
import os
import re
from pathlib import Path
from typing import Any, Mapping

#: analysis.json 的 schema 标识。
ANALYSIS_SCHEMA = "smj.analysis.v1"

#: vendored ``analyze_te`` 必须给出的字段（缺失即视为口径事故）。
REQUIRED_KEYS: tuple[str, ...] = (
    "source",
    "temperature_K",
    "points",
    "L0",
    "L1",
    "L2",
    "log10_G_over_G0_like",
    "S_uV_per_K",
    "k_e_pW_per_K",
    "ZT_electronic",
    "fermi_nearest",
    "max_transmission",
)

#: 参考值默认容差。``abs`` = 绝对容差；``rel`` = 相对容差（对参考值的模）。
DEFAULT_TOLERANCES: dict[str, dict[str, float]] = {
    "log10_G_over_G0_like": {"abs": 0.1},
    "S_uV_per_K": {"rel": 0.10},
    "ZT_electronic": {"rel": 0.10},
    "T_E_F": {"rel": 0.05},
    "k_e_pW_per_K": {"rel": 0.10},
}

#: 参考值常用别名 -> 结果字段（或结果内的取值路径）。
REFERENCE_ALIASES: dict[str, str] = {
    "log10_G": "log10_G_over_G0_like",
    "log10_G_over_G0": "log10_G_over_G0_like",
    "G_over_G0": "L0",
    "S": "S_uV_per_K",
    "S_uV_per_K": "S_uV_per_K",
    "ZT": "ZT_electronic",
    "ZT_e": "ZT_electronic",
    "T_E_F": "T(E_F)",
    "T_EF": "T(E_F)",
    "T(E_F)": "T(E_F)",
    "T_EF_G0": "T(E_F)",
    "k_e": "k_e_pW_per_K",
}


# --------------------------------------------------------------------------- #
# 小工具
# --------------------------------------------------------------------------- #
def _write_json(payload: Mapping[str, Any], path: str | Path) -> Path:
    """原子写 JSON（``.tmp`` -> ``os.replace``），父目录自动创建。

    ``allow_nan=False``：NaN/Infinity 不是 JSON 数字，写出去会让严格解析器
    （以及下游判据）把损坏证据读成正常值。非有限量应在构造 payload 时转 None。
    """
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(dict(payload), indent=2, ensure_ascii=False,
                      allow_nan=False) + "\n"
    tmp = target.with_name(target.name + ".tmp")
    tmp.write_text(text, encoding="utf-8", newline="\n")
    os.replace(tmp, target)
    return target


def resolve_te(path: str | Path) -> Path:
    """定位 TE.dat：接受文件本身、run 根目录或阶段目录。"""
    p = Path(path)
    if p.is_file():
        return p
    candidates = (
        p / "TE.dat",
        p / "aitranss_te" / "TE.dat",
        p / "aims_scf" / "TE.dat",
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        "找不到 TE.dat：%s（试过 %s）" % (p, ", ".join(str(c) for c in candidates)))


def _lookup(result: Mapping[str, Any], key: str) -> Any:
    """按参考值键取结果值；支持别名与 ``T(E_F)`` 这类派生键。"""
    target = REFERENCE_ALIASES.get(key, key)
    if target == "T(E_F)":
        nearest = result.get("fermi_nearest") or {}
        return nearest.get("T") if isinstance(nearest, Mapping) else None
    return result.get(target)


def _tolerance_for(key: str, tolerances: Mapping[str, Any] | None) -> dict[str, Any]:
    merged = dict(DEFAULT_TOLERANCES)
    if tolerances:
        for k, v in tolerances.items():
            if isinstance(v, Mapping):
                merged[k] = dict(v)
            elif isinstance(v, (int, float)):
                merged[k] = {"rel": float(v)}
    spec = dict(merged.get(key) or {"rel": 0.05})
    if "abs" not in spec and "rel" not in spec:
        spec = {"rel": 0.05}
    return spec


# --------------------------------------------------------------------------- #
# 独立复核 / 参考值对比
# --------------------------------------------------------------------------- #
def independent_check(result: Mapping[str, Any], te_path: str | Path, *,
                      temperature: float | None = None) -> dict:
    """独立的 L0/L1/L2 复核（不同解析器 + Simpson 求积 + CODATA 常数）。

    返回 ``{ok, L0, L1, L2, G_over_G0, S_uV_per_K, ZT_electronic, comparison}``；
    独立复核失败不影响主管线结果，调用方在 ``warnings`` 里记录原因。
    """
    from smj.junction import independent_check as ind  # lazy: vendored

    temp = float(temperature if temperature is not None else result.get("temperature_K") or 300.0)
    rows = ind.parse_te(te_path)
    independent = ind.l0_l1_l2(rows, temp)
    derived = ind.derived(independent, temp)
    pipeline = {k: result.get(k) for k in ("L0", "L1", "L2")}
    comparison = ind.compare(pipeline, independent, temperature=temp)
    return {
        "ok": bool(comparison.get("ok")),
        "L0": independent["L0"],
        "L1": independent["L1"],
        "L2": independent["L2"],
        **derived,
        "comparison": comparison,
        "temperature_K": temp,
        "points": independent.get("points"),
    }


def compare_reference(result: Mapping[str, Any], reference: Mapping[str, Any], *,
                      tolerance: Mapping[str, Any] | None = None) -> dict:
    """结果 vs 参考值逐量对比。

    ``reference`` 支持 ``log10_G_over_G0_like`` / ``S_uV_per_K`` / ``ZT_electronic``
    / ``T_E_F`` 及其常见别名。默认容差见 :data:`DEFAULT_TOLERANCES`：
    log10(G/G0) 绝对 0.1、S 相对 10%。返回 ``{ok, tolerance, items}``，
    不抛异常（缺字段即 ``ok=False`` 并给原因）。
    """
    items: list[dict[str, Any]] = []
    for key, expected in dict(reference or {}).items():
        actual = _lookup(result, key)
        spec = _tolerance_for(key, tolerance)
        item: dict[str, Any] = {
            "key": key,
            "field": REFERENCE_ALIASES.get(key, key),
            "reference": expected,
            "actual": actual,
            "tolerance": spec,
            "ok": False,
        }
        if actual is None:
            item["reason"] = "结果中没有该字段"
            items.append(item)
            continue
        try:
            a, e = float(actual), float(expected)
        except (TypeError, ValueError):
            item["reason"] = "参考值/实际值非数值"
            items.append(item)
            continue
        delta = a - e
        if "abs" in spec:
            limit = float(spec["abs"])
            ok = abs(delta) <= limit
            item["limit"] = limit
        else:
            limit = abs(float(spec.get("rel", 0.05)) * e)
            ok = abs(delta) <= limit
            item["limit"] = limit
            item["rel"] = (abs(delta) / abs(e)) if e else None
        item.update({"delta": delta, "ok": bool(ok)})
        if not ok:
            item["reason"] = "超出容差"
        items.append(item)
    return {"ok": all(i["ok"] for i in items), "tolerance": dict(tolerance or {}),
            "items": items}


# --------------------------------------------------------------------------- #
# 出图（可选依赖）
# --------------------------------------------------------------------------- #
#: T(E) 图纵轴一律对数；T<=0 的点压到这个下限，避免 log(0) 丢整条曲线。
T_LOG_FLOOR = 1e-12


def _log_floor(values: list[float]) -> float:
    """对数纵轴下界：最小正值的一半（不低于 :data:`T_LOG_FLOOR`）。"""
    positive = [float(v) for v in values if float(v) > 0.0]
    if not positive:
        return 1.0
    return max(min(positive) * 0.5, T_LOG_FLOOR)


def plot_transmission(te_path: str | Path, result: Mapping[str, Any],
                      out_png: str | Path, *, title: str | None = None) -> dict:
    """T(E) 曲线（E-E_F vs T，**纵轴对数**）。兼容入口：版式与导出统一走
    :func:`smj.plotting.plot_te`（需要 matplotlib，extra ``plot``），缺失即抛错。
    """
    from .plotting import plot_te

    meta = plot_te(te_path, out_png, scale="log", title=title, result=result,
                   formats=("png",), annotate=True)
    return {"path": meta["path"], "points": meta["points"],
            "log_scale": True, "floor": meta["floor"]}


# --------------------------------------------------------------------------- #
# 主入口
# --------------------------------------------------------------------------- #
def analyze_te(te_path: str | Path, *, temperature: float = 300.0,
               output: str | Path | None = None, plot: str | Path | bool | None = None,
               reference: Mapping[str, Any] | None = None,
               tolerance: Mapping[str, Any] | None = None,
               independent: bool = True) -> dict:
    """对一个 TE.dat 做输运分析，返回固定字段摘要（可选写 JSON / 出图 / 对比）。

    Parameters
    ----------
    te_path:
        ``TE.dat`` 路径（或包含它的 run/阶段目录，见 :func:`resolve_te`）。
    temperature:
        费米窗温度（K）。``analyze_te`` 的 ``T(E)`` 单位是 G0，温度只影响权重。
    output:
        若给出则写 ``analysis.json``（父目录自动创建）。
    plot:
        ``None``/``False`` 不出图；``True`` 写到 TE.dat 同目录 ``T_E_plot.png``；
        字符串/Path 则写到该路径。matplotlib 缺失时**跳过**并在 ``plot.skipped``
        说明，不抛错。
    reference / tolerance:
        参考值 dict 与容差覆盖；结果写入 ``comparison``（不通过不抛错）。
    independent:
        是否跑独立复核（默认 True）；失败只记 ``warnings``，不影响主结果。
    """
    te = resolve_te(te_path)
    from smj.junction.landauer import analyze_te as _analyze  # lazy: vendored

    core = _analyze(te, temperature=float(temperature))
    missing = [k for k in REQUIRED_KEYS if k not in core]
    if missing:  # vendored 口径变了 —— 显式失败，不静默降级
        raise RuntimeError("vendored analyze_te 缺字段: %s" % ", ".join(missing))

    result: dict[str, Any] = dict(core)
    result["schema"] = ANALYSIS_SCHEMA
    result["source"] = str(te)
    # 保留 vendored 的口径告警（如 TE.dat 行序非升序被排序），不要在这里清空：
    # 旧实现直接覆盖成 []，会让"倒序谱被静默排序"这条唯一提示消失。
    result["warnings"] = list(core.get("warnings") or [])

    if independent:
        try:
            result["independent"] = independent_check(core, te, temperature=temperature)
        except Exception as exc:  # noqa: BLE001 —— 复核失败不得拖垮主结果
            result["independent"] = {"ok": False, "error": str(exc)}
            result["warnings"].append("独立复核失败: %s" % exc)

    if reference:
        result["comparison"] = compare_reference(core, reference, tolerance=tolerance)

    if plot:
        target = (te.parent / "T_E_plot.png") if plot is True else Path(str(plot))
        try:
            result["plot"] = plot_transmission(te, core, target)
        except Exception as exc:  # noqa: BLE001 —— matplotlib 是可选依赖
            result["plot"] = {
                "skipped": True,
                "reason": ("matplotlib 不可用（extra \"plot\" 未安装）：%s；"
                           "analysis.json 仍可用" % exc),
            }

    if output:
        result["analysis_path"] = str(Path(output))
        result["analysis_path"] = str(_write_json(result, output))
    return result


def analyze_run(run_dir: str | Path | Any, *, temperature: float | None = None,
                output: str | Path | None = None,
                plot: str | Path | bool | None = None,
                reference: Mapping[str, Any] | None = None,
                tolerance: Mapping[str, Any] | None = None,
                independent: bool = True, force: bool = False) -> dict:
    """分析一个 run 目录（默认 ``<run>/aitranss_te/TE.dat`` 与
    ``<run>/junction_analyze/analysis.json``）。

    ``run_dir`` 也接受 :class:`smj.config.RunConfig`（CLI 直接传配置对象），
    此时温度取 ``config.transport.temperature``。
    """
    if hasattr(run_dir, "run_dir") and not isinstance(run_dir, (str, Path)):
        cfg = run_dir
        run = Path(cfg.run_dir)
        if temperature is None:
            temperature = getattr(cfg.transport, "temperature", None)
    else:
        run = Path(run_dir)
    te = resolve_te(run)
    temp = float(temperature if temperature is not None else 300.0)

    target = Path(output) if output else (run / "junction_analyze" / "analysis.json")
    if target.exists() and not force:
        return {
            "skipped": True,
            "reason": "analysis.json 已存在（用 force=True 覆盖）: %s" % target,
            "analysis_path": str(target),
            "source": str(te),
        }
    result = analyze_te(te, temperature=temp, output=target, plot=plot,
                        reference=reference, tolerance=tolerance,
                        independent=independent)
    result["run_id"] = run.name
    return result


def write_analysis(result: Mapping[str, Any], path: str | Path) -> Path:
    """把分析结果写到指定 JSON（幂等；供 CLI/测试直接调用）。"""
    return _write_json(result, path)


def orbital_alignment(te_path: str | Path, aims_out: str | Path, *,
                      mode: str = "fermi", target_e_ev: float | None = None,
                      output: str | Path | None = None) -> dict:
    """可选：把 TE.dat 的费米窗与 FHI-aims 轨道本征值对齐（vendored spectra）。"""
    from smj.junction.spectra import match_orbital  # lazy: vendored
    return match_orbital(te_path, aims_out, mode=mode, target_e_ev=target_e_ev,
                         output=output)


def read_analysis(path: str | Path) -> dict:
    """读回 analysis.json（报告模块复用）。"""
    return json.loads(Path(path).read_text(encoding="utf-8"))


def te_header(te_path: str | Path) -> dict:
    """从 TE.dat 表头/数据行提取 ``bias`` / ``efermi`` / ``points``（不臆断）。"""
    te = resolve_te(te_path)
    text = Path(te).read_text(encoding="utf-8", errors="ignore")
    header = text.splitlines()[:8]
    out: dict[str, Any] = {"source": str(te)}
    for line in header:
        m = re.search(r"efermi\s*=\s*([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)", line,
                      re.IGNORECASE)
        if m:
            out["efermi_Ha"] = float(m.group(1))
        m = re.search(r"bias\s*=\s*([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)", line,
                      re.IGNORECASE)
        if m:
            out["bias_V"] = float(m.group(1))
    points = 0
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) >= 3:
            try:
                float(parts[0].replace("D", "E"))
                float(parts[1].replace("D", "E"))
                float(parts[2].replace("D", "E"))
            except ValueError:
                continue
            points += 1
    out["points"] = points
    return out


#: CLI/兼容别名（首选名是 analyze_run）。
analyze = analyze_run

# --------------------------------------------------------------------------- #
# 热导（声子）产物的 schema 校验（2026-11-23, M5-T4）
#
# `thermal.json` 由 `thermal_analyze` 阶段写、由 `review`/`report` 读。这里的
# 作用与 `analysis.json` 一致：**只做契约校验，不做物理判据** —— 判据在
# `smj.phonon.checks` 与 `smj.junction.criteria`，物理结论不在这里下。
# --------------------------------------------------------------------------- #

#: thermal.json 的 schema 标识（与 `thermal_analyze.SCHEMA` 同值）。
THERMAL_SCHEMA = "smj.thermal.v1"

#: 必需键。缺任何一个都 raise —— 静默补默认值会让"没算"看起来像"算过了"。
THERMAL_REQUIRED_KEYS = (
    "schema", "temperature_K", "n_atoms", "n_device_atoms", "n_lead_atoms",
    "g_ph_pw_per_k", "g_quantum_pw_per_k", "n_channels_effective",
    "acoustic_sum_residual", "acoustic_sum_relative_residual",
    "imaginary_modes", "n_k_points", "band_top_meV", "eta_meV",
    "energy_step_meV", "spectrum", "zt",
)

#: spectrum 段的三个等长数组。
THERMAL_SPECTRUM_KEYS = ("energy_meV", "transmission", "cumulative_g_ph_pw_per_k")

#: 允许出现在 thermal.json 里但不在必需键中的附加键（多出来的保留并告警）。
THERMAL_OPTIONAL_KEYS = ("relax_max_force_ev_per_A", "lead_fc_sha256", "fc_asymmetry", "notes",
                         # D13：ASR 投影记录。旧文件没有这三个键，故列为可选。
                         "asr_projected", "asr_projection_norm",
                         "asr_residual_after_ev_per_A2", "asr_region",
                         # 软模归属（modes.describe_soft_modes）。解释性字段，
                         # 缺了不影响任何判据，故可选。
                         "soft_modes",
                         # 免费的温度曲线（同一份透射谱，只换玻色因子）
                         "g_ph_by_temperature",
                         # 把负本征值投影掉后 G_ph 的相对变化（None = 无虚频）
                         "soft_mode_impact_relative_g_ph",
                         "g_ph_symmetric_pw_per_k", "g_ph_soft_removed_pw_per_k",
                         # 三条收敛轴的状态（截断已测，另两条登记为未测）。
                         "convergence",
                         # 环境清单指纹（正文在 environment.json）。
                         "environment")


def read_thermal(path: str | Path) -> tuple[dict, list[str]]:
    """读并校验 ``thermal.json``，返回 ``(payload, warnings)``。

    校验项：schema 值、必需键齐备、spectrum 三数组等长且长度为 >=2、
    数值有限（``allow_nan=False`` 语义）。任何一项不过都 raise ``ValueError``
    —— 与 :func:`read_analysis` 的口径一致：**契约坏掉时宁可炸，不要猜**。
    """
    p = Path(path)
    if not p.is_file():
        raise ValueError("thermal.json 不存在: %s" % p)
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError("thermal.json 不是合法 JSON: %s" % exc) from exc
    if not isinstance(data, dict):
        raise ValueError("thermal.json 顶层必须是对象")
    if data.get("schema") != THERMAL_SCHEMA:
        raise ValueError("thermal.json 的 schema 是 %r，期望 %r"
                         % (data.get("schema"), THERMAL_SCHEMA))
    missing = [k for k in THERMAL_REQUIRED_KEYS if k not in data]
    if missing:
        raise ValueError("thermal.json 缺必需键: %s" % ", ".join(missing))

    warnings: list[str] = []
    extra = sorted(set(data) - set(THERMAL_REQUIRED_KEYS) - set(THERMAL_OPTIONAL_KEYS))
    if extra:
        warnings.append("thermal.json 含未知键，已保留: %s" % ", ".join(extra))

    spectrum = data["spectrum"]
    if not isinstance(spectrum, dict):
        raise ValueError("thermal.json 的 spectrum 必须是对象")
    absent = [k for k in THERMAL_SPECTRUM_KEYS if k not in spectrum]
    if absent:
        raise ValueError("thermal.json 的 spectrum 缺: %s" % ", ".join(absent))
    lens = {k: len(spectrum[k]) for k in THERMAL_SPECTRUM_KEYS}
    if len(set(lens.values())) != 1:
        raise ValueError("thermal.json 的 spectrum 三数组长度不一致: %r" % lens)
    if lens["energy_meV"] < 2:
        raise ValueError("thermal.json 的 spectrum 至少要有 2 个能量点")
    for key in THERMAL_SPECTRUM_KEYS:
        for value in spectrum[key]:
            if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
                raise ValueError("thermal.json 的 spectrum.%s 含非有限数值" % key)
    for key in ("g_ph_pw_per_k", "g_quantum_pw_per_k", "n_channels_effective"):
        try:
            number = float(data[key])
        except (TypeError, ValueError) as exc:
            raise ValueError("thermal.json 的 %s 不是数值" % key) from exc
        if not math.isfinite(number):
            raise ValueError("thermal.json 的 %s 非有限" % key)
    if float(data["g_ph_pw_per_k"]) < 0.0:
        raise ValueError("thermal.json 的 g_ph_pw_per_k 为负：%r"
                         % data["g_ph_pw_per_k"])
    return data, warnings
