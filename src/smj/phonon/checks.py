"""热导（声子）判据的机械评估（M5-T2 的数据侧）。

分工：
- ``smj.junction.criteria`` 登记**判据的定义**（名称/门限/源），是唯一的口径来源；
- 本模块把 ``thermal.json`` 的数值喂进这些定义，产出一份 ``checks`` 清单；
- ``smj.review`` 只负责把清单翻译成 review.json 的 issues 与结论。

三条设计约束（都来自方案 §11）：
1. **缺证据不判 pass**。拿不到的量（例如没留 ``thermal_relax/aims.out`` 时的
   ``max|F|``）记成 ``missing``，整体归 ``inconclusive``，不当作通过。
2. **error 门 → failed，warn 门 → inconclusive**，与 ``physics_verdict`` 同构。
3. **不在这里改数据**。评估是只读的；任何"调参让它过"的冲动都应该写进报告。
"""
from __future__ import annotations

import math
from typing import Any, Mapping

import numpy as np

from ..junction import criteria as _criteria
from .landauer import (ACOUSTIC_ENERGY_EV, ACOUSTIC_TOL, DEFAULT_MAX_CHANNELS,
                       MAX_ABOVE_BAND_T, NEGATIVE_TRANSMISSION_TOL, above_band,
                       conductance_bound, thermal_conductance)

__all__ = ["evaluate_thermal", "summarize"]

#: 力收敛门限（eV/Å）。与 ``thermal_relax`` 的 ``relax_force_threshold`` 默认一致。
MAX_RELAX_FORCE = 1.0e-2

#: 网格独立性门限：粗化一半后 G_ph 的相对变化上限。
#:
#: (D14, 2026-11-25) 1% 不可达 —— 本方法的被积函数在峰宽 ~eta 的窄结构上求和,
#: 梯形积分的求积噪声实测 2.0%-2.5%(nap-sh, eta=0.05/0.1/0.2/0.3 meV 各测一遍
#: 都是这个量级), 加密网格反而在 1.0%-2.0% 之间摆动, 压不到 1% 以下。
#: 真正的收敛控制量是 eta(见 `thermal.negf.eta_meV`), 不是步长;
#: 因此门限取实测底噪的两倍, 只用来抓"网格明显不够" (>5%) 这一档。
MAX_GRID_DRIFT = 5.0e-2

#: 声学支在 ω→0 的透射上限（k 平均口径，三条声学支）。
MAX_ACOUSTIC_T = 3.0

#: ASR 的**相对**残差门限（``max|Σ_j K_ij| / max|K|``，无量纲）。
#:
#: 绝对门限（1e-6 eV/Å²）不可用：有限差分 + 基组噪声本身就在 1e-3–1e-1 eV/Å²
#: 量级，拿绝对值判会让每个 run 都 inconclusive。ASR 违反本来就该与 ``max|K|``
#: 比，这也是 phonopy 社区的口径。超过本门限时 ``thermal_analyze`` 会做一次
#: ASR 投影，并把投影前后的残差与投影范数一并写进 ``thermal.json``（D13）。
MAX_ASR_RELATIVE = 1.0e-3


def _check(name: str, *, passed: bool | None, detail: str) -> dict:
    spec = _criteria._spec(name)
    return {
        "name": name,
        "kind": spec.get("kind"),
        "gate": spec.get("gate"),
        "passed": passed,
        "missing": passed is None,
        "detail": detail,
    }


def evaluate_thermal(
    thermal: Mapping[str, Any],
    *,
    relax_max_force: float | None = None,
    relax_sc_accuracy_forces: str | None = None,
    lead_cell_match: bool | None = None,
    fc_asymmetry: float | None = None,
    max_channels: float = DEFAULT_MAX_CHANNELS,
) -> dict:
    """按 ``CRITERIA`` 的十条热导判据评估一份 ``thermal.json``。

    返回 ``{"checks": [...], "verdict": "passed"|"inconclusive"|"failed",
    "failed": [names], "inconclusive": [names]}``。
    """
    spectrum = thermal["spectrum"]
    energy = np.asarray(spectrum["energy_meV"], dtype=float) * 1.0e-3   # eV
    t_ph = np.asarray(spectrum["transmission"], dtype=float)
    g_ph = float(thermal["g_ph_pw_per_k"])
    band_top_mev = float(thermal["band_top_meV"])
    temperature = float(thermal["temperature_K"])

    checks: list[dict] = []

    checks.append(_check(
        "ph_transmission_present", passed=bool(np.all(t_ph >= -NEGATIVE_TRANSMISSION_TOL)),
        detail="spectrum 有 %d 个能量点，最小透射 %.3g" % (energy.size, np.min(t_ph))))

    if relax_max_force is None:
        detail = "没找到 thermal_relax/aims.out 的力块，无法判定收敛"
        passed: bool | None = None
    elif relax_sc_accuracy_forces is None:
        detail = "max|F|=%.3e eV/A，但 control.in 里没有 sc_accuracy_forces" % relax_max_force
        passed = False
    else:
        passed = relax_max_force <= MAX_RELAX_FORCE
        detail = "max|F|=%.3e eV/A（门限 %.1e），sc_accuracy_forces=%s" % (
            relax_max_force, MAX_RELAX_FORCE, relax_sc_accuracy_forces)
    checks.append(_check("ph_relax_converged", passed=passed, detail=detail))

    above = above_band(spectrum["energy_meV"], band_top_mev)
    if np.any(above):
        t_above = float(np.max(np.abs(t_ph[above])))
        passed = t_above <= MAX_ABOVE_BAND_T
        detail = "带顶 %.1f meV 之上 %d 点，max|T|=%.1e（门限 %.0e）" % (
            band_top_mev, int(above.sum()), t_above, MAX_ABOVE_BAND_T)
    else:
        passed = None
        detail = "带顶 %.1f meV 以上缺少采样" % band_top_mev
    checks.append(_check("ph_band_top_zero", passed=passed, detail=detail))

    bound = conductance_bound(temperature, max_channels)
    checks.append(_check(
        "ph_quantum_bound", passed=g_ph <= bound,
        detail="G_ph=%.3f pW/K，%g 通道上限为 %.3f pW/K（%.1f K）"
               % (g_ph, max_channels, bound, temperature)))

    low = energy <= ACOUSTIC_ENERGY_EV
    t_low = float(t_ph[0]) if t_ph.size else float("nan")
    checks.append(_check(
        "ph_acoustic_limit",
        passed=(t_low <= MAX_ACOUSTIC_T + ACOUSTIC_TOL) if low.any() else None,
        detail="T_ph(%.3f meV)=%.4f（k 平均口径，上限 %.0f）"
               % (energy[0] * 1e3, t_low, MAX_ACOUSTIC_T)))

    asr = float(thermal["acoustic_sum_residual"])
    asr_rel = thermal.get("acoustic_sum_relative_residual")
    projected = bool(thermal.get("asr_projected", False))
    if isinstance(asr_rel, (int, float)) and not math.isnan(float(asr_rel)):
        passed = float(asr_rel) <= MAX_ASR_RELATIVE
        detail = ("相对残差 %.2e（门限 %.0e；绝对 %.3e eV/A^2）"
                  % (float(asr_rel), MAX_ASR_RELATIVE, asr))
    else:                          # 旧文件没有相对残差,退回绝对口径
        passed = asr <= 1.0e-6
        detail = "max|sum_j K_ij|=%.3e eV/A^2（无相对残差记录）" % asr
    if projected:
        norm = thermal.get("asr_projection_norm")
        detail += ("；已做 ASR 投影（相对投影量 %s）"
                   % ("%.2e" % float(norm) if isinstance(norm, (int, float)) else "未记录"))
    transport = thermal.get("transport_asr") or {}
    retained = transport.get("relative_residual")
    if isinstance(retained, (int, float)) and math.isfinite(retained):
        passed = passed and retained <= MAX_ASR_RELATIVE
        detail += "；输运分区投影前相对残差 %.2e" % retained
        if thermal.get("asr_projection_reason") == "incompatible_fixed_coupling":
            detail += "；固定界面耦合与平移约束不相容"
    else:
        if passed:
            passed = None
        detail += "；输运分区 ASR 未记录"
    checks.append(_check("ph_asr_residual", passed=passed, detail=detail))

    if fc_asymmetry is None:      # 未显式给出时退回 thermal.json 里记录的值
        recorded = thermal.get("fc_asymmetry")
        if isinstance(recorded, (int, float)):
            fc_asymmetry = float(recorded)
    if fc_asymmetry is None:
        checks.append(_check(
            "ph_fc_symmetry", passed=None, detail="thermal.json 没记录 fc_asymmetry"))
    else:
        checks.append(_check(
            "ph_fc_symmetry", passed=fc_asymmetry < 5.0e-2,
            detail="||K-K^T||/||K||=%.2f%%（门限 5%%）" % (fc_asymmetry * 100.0)))

    modes = thermal["imaginary_modes"]
    count = int(modes.get("count", 0))
    worst = modes.get("worst_cm1")
    checks.append(_check(
        "ph_imaginary_modes", passed=count == 0,
        detail=("%d 个虚频，最差 %s cm^-1（容限 %.0f）"
                % (count, "%.2f" % worst if worst is not None else "—",
                   _criteria.PH_IMAGINARY_TOL_CM1)) if count else "无虚频"))

    if lead_cell_match is None:
        checks.append(_check(
            "ph_lead_cell_match", passed=None, detail="没读到 lead_fc.json 的胞参数"))
    else:
        checks.append(_check(
            "ph_lead_cell_match", passed=lead_cell_match,
            detail="冻结 lead FC 的面内胞与本次器件配置一致" if lead_cell_match
                   else "冻结 lead FC 的面内胞与本次器件配置**不一致**"))

    passed = None
    detail = "能量网格步长 %.3f meV，粗化一半后重算" % float(thermal["energy_step_meV"])
    try:
        coarse = slice(None, None, 2)
        if energy[coarse].size >= 2 and g_ph > 0.0:
            g_coarse = thermal_conductance(
                energy[coarse], t_ph[coarse], temperature,
                band_top_mev=band_top_mev, max_channels=max_channels)
            drift = abs(g_coarse - g_ph) / g_ph
            passed = drift < MAX_GRID_DRIFT
            detail += "：ΔG/G=%.3e（门限 %.0e）" % (drift, MAX_GRID_DRIFT)
        else:
            detail += "：点数不足或 G_ph=0，无法比较"
    except (ValueError, RuntimeError) as exc:      # 粗网格触发闸门 -> 证据不足
        detail += "：粗化重算失败（%s）" % exc
    checks.append(_check("ph_grid_independence", passed=passed, detail=detail))

    return summarize(checks)


def summarize(checks: list[dict]) -> dict:
    """把 checks 归约成 verdict —— error 级未过即 failed，其余未过即 inconclusive。"""
    failed = [c["name"] for c in checks
              if c["passed"] is False and c.get("gate") == "error"]
    inconclusive = [c["name"] for c in checks
                    if c["passed"] is not True and c.get("gate") != "error"]
    # warn 门未过也归 inconclusive；error 门 missing 同样是证据不足而非通过
    inconclusive += [c["name"] for c in checks
                     if c["passed"] is None and c.get("gate") == "error"]
    if failed:
        verdict = "failed"
    elif inconclusive:
        verdict = "inconclusive"
    else:
        verdict = "passed"
    return {"checks": checks, "verdict": verdict,
            "failed": failed, "inconclusive": inconclusive}
