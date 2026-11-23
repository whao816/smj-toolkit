"""thermal_analyze（声子输运）：力常数 → 动力学矩阵 → T_ph(ω) → G_ph(T)（本地阶段）。

这一步是全链唯一产出物理量的地方，也是**唯一**做单位换算的地方：

1. ``force_constants`` 组装成完整的 ``3N×3N``（器件行来自 DFT，铅行来自冻结铅 FC）；
2. ``D = K/√(m_i m_j) × (0.0646542)²`` —— 因子把 ``eV/(Å²·amu)`` 换成 ``eV²``
   （方案 §3.2 与附录 A-1：报告里的 0.2543 是把正确系数又开了一次方，偏大 3.93 倍）；
3. 铅侧取**一个主层**（D11），做 López-Sancho decimation 得表面格林函数；
4. Caroli 公式 ``T_ph = Tr[G^r Γ_L G^a Γ_R]``，``N_k=1``（D10：3×3 超胞的 Γ）；
5. Landauer 积分，三条 fail-closed（带顶非零 / 声学支 > 3 / G > 3G_Q）。

``ZT_total = ZT_e / (1 + G_ph/k_e)`` 只在本阶段**附加**给出，**不参与门控**：
它是 engineering 参考量，``thermal.json`` 里与 ``G_ph`` 分字段存放。
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ..errors import StageBlocked, StageError
from ..phonon.dynamical import (acoustic_sum_relative_residual,
                               acoustic_sum_residual, dynamical_matrix,
                               imaginary_modes, masses_from_elements)
from ..phonon.fc import FullFCMeta, assemble_full_fc, load_lead_fc
from ..phonon.landauer import cumulative_thermal_conductance, thermal_conductance
from ..phonon.negf import transmission
from ..phonon.units import gq_pw_per_k
from .base import Collected, Prepared, Stage, register

SCHEMA = "smj.thermal.v1"


def _lead_fc_path(value: str | Path) -> Path:
    """``thermal.lead_fc`` 可以写成目录(用其中的 ``lead_fc.json``)或直接写文件。

    ``save_lead_fc`` 落盘的是 ``lead_fc.npy`` + ``lead_fc.json`` 两个文件,
    用户按目录配置更自然,这里两种都接受。
    """
    path = Path(value)
    if path.is_dir():
        candidate = path / "lead_fc.json"
        if candidate.is_file():
            return candidate
    return path


@register
class ThermalAnalyzeStage(Stage):
    name = "thermal_analyze"
    depends = ("thermal_fc",)
    needs_submit = False
    input_hint = "消费 thermal_fc/force_constants.npy 与 thermal.lead_fc"
    outputs = ("thermal.json",)

    def prepare(self, ctx, stage_dir: Path) -> Prepared:
        stage_dir = Path(stage_dir)
        thermal = ctx.config.thermal
        fc_path = ctx.artifact("thermal_fc", "force_constants.npy")
        if not Path(fc_path).is_file():
            raise StageBlocked(
                self.name, "缺 thermal_fc/force_constants.npy：请先跑 finite_fc 并 collect",
                [str(fc_path)])
        if not thermal.lead_fc:
            raise StageBlocked(
                self.name,
                "配置缺 thermal.lead_fc：声子 NEGF 需要冻结铅力常数（方案 D4）。"
                "先用一份 bulk/厚板算一次并 save_lead_fc 落盘。",
                [])
        path = _lead_fc_path(thermal.lead_fc)
        if not path.is_file():
            raise StageBlocked(self.name, "thermal.lead_fc 指向的文件不存在: %s" % path,
                               [str(path)])
        return self.prepared(
            ctx, stage_dir,
            inputs=[Path(fc_path), path],
            notes=[
                "本阶段为纯本地后处理：读力常数、算透射、做 Landauer 积分，"
                "不提交任何作业",
                "温度 %.3f K；N_k=1（3×3 超胞的 Γ，Γ 已把 9 个原胞 q 点折叠进来）"
                % thermal.resolved_temperature(ctx.config.transport.temperature),
            ],
            commands=[],
            parent_stage="thermal_fc",
        )

    def collect(self, ctx, stage_dir: Path) -> Collected:
        stage_dir = Path(stage_dir)
        thermal = ctx.config.thermal
        temperature = thermal.resolved_temperature(ctx.config.transport.temperature)

        k_dev = np.load(ctx.artifact("thermal_fc", "force_constants.npy"))
        device_meta = json.loads(
            ctx.artifact("thermal_device", "device_metadata.json").read_text(
                encoding="utf-8"))
        lead = load_lead_fc(_lead_fc_path(thermal.lead_fc))

        n_atoms = int(device_meta["n_atoms"])
        elements = self._elements(ctx.artifact("thermal_fc", "geometry.in"))
        if len(elements) != n_atoms:
            raise StageError(
                "thermal_analyze 的元素数(%d)与 device_metadata 的 n_atoms(%d) 不一致"
                % (len(elements), n_atoms))
        masses = masses_from_elements(elements)

        device_indices = [int(i) for i in device_meta["device_indices"]]
        lead_indices = [int(i) for i in device_meta["lead_indices"]]
        meta = FullFCMeta(n_atoms=n_atoms, device_indices=device_indices,
                          left_indices=[int(i) for i in device_meta["left_lead_indices"]],
                          right_indices=[int(i) for i in device_meta["right_lead_indices"]])
        full = assemble_full_fc(k_dev, lead, meta)
        residual = acoustic_sum_residual(full)
        relative = acoustic_sum_relative_residual(full)
        dynamical = dynamical_matrix(full, masses)
        imag = imaginary_modes(dynamical)

        layers = device_meta["layer_index"]
        sides = device_meta["side"]
        left_pl = [i for i in lead_indices if layers[i] == 1 and sides[i] == "left"]
        right_pl = [i for i in lead_indices if layers[i] == 1 and sides[i] == "right"]
        left_outer = [i for i in lead_indices if layers[i] == 0 and sides[i] == "left"]
        right_outer = [i for i in lead_indices if layers[i] == 0 and sides[i] == "right"]
        if not (left_pl and right_pl and left_outer and right_outer):
            raise StageError(
                "thermal_analyze 找不到两侧的主层或最外层：layers=%s 时每侧至少需要 2 层"
                % sorted(set(layers)))

        centre = np.array(device_indices)
        d_c = dynamical[np.ix_(centre, centre)]
        d_l = dynamical[np.ix_(left_pl, left_pl)]
        d_r = dynamical[np.ix_(right_pl, right_pl)]
        # negf.transmission 的 v 约定是「中心区 × 铅」,即 (n_c, n_lead)
        v_l = dynamical[np.ix_(centre, left_pl)]
        v_r = dynamical[np.ix_(centre, right_pl)]
        v_lead_l = dynamical[np.ix_(left_pl, left_outer)]
        v_lead_r = dynamical[np.ix_(right_pl, right_outer)]

        negf = thermal.negf
        band_top = float(negf.band_top_meV)
        step = float(negf.energy_step_meV)
        n_points = int(round(band_top / step)) + 1
        energy_mev = np.linspace(0.0, band_top * 1.25, max(n_points, 2))
        omega = energy_mev * 1.0e-3
        eta = float(negf.eta_meV) * 1.0e-3
        t_ph = np.array([
            transmission(float(w), d_c, d_l, d_r, v_l, v_r, eta ** 2,
                        v_lead_l=v_lead_l, v_lead_r=v_lead_r)
            for w in omega])
        g_ph = thermal_conductance(omega, t_ph, temperature, band_top_mev=band_top)
        cumulative = cumulative_thermal_conductance(omega, t_ph, temperature)

        payload = {
            "schema": SCHEMA,
            "temperature_K": temperature,
            "n_atoms": n_atoms,
            "n_device_atoms": len(device_indices),
            "n_lead_atoms": len(lead_indices),
            "g_ph_pw_per_k": g_ph,
            "g_quantum_pw_per_k": gq_pw_per_k(temperature),
            "n_channels_effective": g_ph / gq_pw_per_k(temperature),
            "acoustic_sum_residual": residual,
            "acoustic_sum_relative_residual": relative,
            "imaginary_modes": imag,
            "n_k_points": 1,
            "band_top_meV": band_top,
            "eta_meV": float(negf.eta_meV),
            "energy_step_meV": step,
            "spectrum": {
                "energy_meV": energy_mev.tolist(),
                "transmission": t_ph.tolist(),
                "cumulative_g_ph_pw_per_k": cumulative.tolist(),
            },
        }
        payload["zt"] = self._zt_payload(ctx, g_ph)
        (stage_dir / "thermal.json").write_text(
            json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8")
        results = {"thermal": {
            "g_ph_pw_per_k": g_ph,
            "temperature_K": temperature,
            "acoustic_sum_relative_residual": relative,
            "n_imaginary": imag["count"],
            "worst_imaginary_cm1": imag["worst_cm1"],
            "zt_total": payload["zt"].get("zt_total"),
        }}
        notes = [
            "G_ph(%.0f K) = %.3f pW/K（有效通道数 %.3f，占单通道上限的 %.1f%%）"
            % (temperature, g_ph, g_ph / gq_pw_per_k(temperature),
               100.0 * g_ph / gq_pw_per_k(temperature)),
            "ASR 残差 %.3e eV/Å²（相对 %.3e）；虚频 %d 个（最差 %.1f cm⁻¹）"
            % (residual, relative, imag["count"], imag["worst_cm1"]),
        ]
        if payload["zt"].get("zt_total") is not None:
            notes.append(
                "ZT_total = %.5f（电子 ZT %.5f，G_ph/k_e = %.3f）；仅作参考，不参与门控"
                % (payload["zt"]["zt_total"], payload["zt"]["zt_electronic"],
                   payload["zt"]["g_ph_over_k_e"]))
        return self.collected(
            ctx, stage_dir, artifacts=["thermal.json"], results=results,
            notes=notes)

    # -- 内部 -------------------------------------------------------------- #
    @staticmethod
    def _zt_payload(ctx, g_ph: float) -> dict:
        """从 ``junction_analyze/results.json`` 读电子量,给出总 ZT。

        缺电子结果时返回 ``{"available": False, ...}`` —— 声子链不依赖电子链,
        单独跑 thermal 全链是允许的。
        """
        path = ctx.artifact("junction_analyze", "results.json")
        try:
            data = json.loads(Path(path).read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return {"available": False,
                    "reason": "缺 junction_analyze/results.json（电子链未跑）"}
        block = data.get("junction", data)
        try:
            zt_e = float(block["ZT_electronic"])
            k_e = float(block["k_e_pW_per_K"])
        except (KeyError, TypeError, ValueError):
            return {"available": False, "reason": "results.json 缺 ZT_electronic/k_e"}
        if k_e <= 0:
            return {"available": False, "reason": "k_e 非正，无法定义 ZT_total"}
        return {
            "available": True,
            "zt_electronic": zt_e,
            "k_e_pW_per_K": k_e,
            "g_ph_over_k_e": g_ph / k_e,
            "zt_total": zt_e / (1.0 + g_ph / k_e),
        }

    @staticmethod
    def _elements(geometry: Path) -> list[str]:
        elements = []
        for line in Path(geometry).read_text(encoding="utf-8").splitlines():
            parts = line.split()
            if len(parts) >= 5 and parts[0] == "atom":
                elements.append(parts[4])
        return elements
