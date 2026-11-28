"""thermal_analyze（声子输运）：力常数 → 动力学矩阵 → T_ph(ω) → G_ph(T)（本地阶段）。

这一步是全链唯一产出物理量的地方，也是**唯一**做单位换算的地方：

1. ``force_constants`` 组装成完整的 ``3N×3N``（器件行来自 DFT，lead 行来自冻结 lead FC）；
2. ``D = K/√(m_i m_j) × (0.0646542)²`` —— 因子把 ``eV/(Å²·amu)`` 换成 ``eV²``
   （方案 §3.2 与附录 A-1：报告里的 0.2543 是把正确系数又开了一次方，偏大 3.93 倍）；
3. lead 侧取**一个主层**（D11），做 López-Sancho decimation 得表面格林函数；
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
from ..phonon.assets import (lead_boundary, read_device_structure,
                             resolve_lead_fc_path as _lead_fc_path)
from ..phonon.dynamical import (acoustic_sum_residuals, dynamical_matrix, imaginary_modes,
                                masses_from_elements)
from ..phonon.checks import MAX_ASR_RELATIVE
from ..phonon.fc import (FullFCMeta, asr_diagnostics, assemble_full_fc,
                         load_lead_fc)
from ..phonon.geometry import atom_distances as _atom_distances
from ..phonon.landauer import cumulative_thermal_conductance, thermal_conductance
from ..phonon.modes import describe_soft_modes
from ..phonon.negf import transmission
from ..phonon.units import gq_pw_per_k
from .base import Collected, Prepared, Stage, register

SCHEMA = "smj.thermal.v1"


def _fc_asymmetry(force_constants) -> float:
    """`||K - K^T|| / ||K||`（Frobenius）。供判据 `ph_fc_symmetry`。"""
    matrix = np.asarray(force_constants, dtype=float)
    scale = float(np.linalg.norm(matrix))
    if scale <= 0.0:
        return 0.0
    return float(np.linalg.norm(matrix - matrix.T) / scale)


def _dof_indices(atoms) -> np.ndarray:
    """原子下标 → 自由度下标。

    ``dynamical_matrix`` 返回的是 3N×3N 矩阵,第 ``3*a+k`` 行对应原子 ``a``
    的第 ``k`` 个方向。整条热导链上凡是"取某个原子子集的块",都必须先过这道
    换算;直接拿原子下标去 ``np.ix_`` 会静默切出错位的块(见 collect 里的注释)。
    """
    return np.array([3 * int(a) + k for a in atoms for k in range(3)], dtype=int)


def _energy_grid(band_top: float, step: float) -> tuple[np.ndarray, np.ndarray]:
    """能量网格：到 ``1.25 × 带顶``，给闸门留出"带顶以上必须为零"的取样点。"""
    if not np.isfinite([band_top, step]).all() or min(band_top, step) <= 0:
        raise ValueError("带顶和能量步长必须是正的有限值")
    n_intervals = int(np.ceil(1.25 * band_top / step))
    energy_mev = np.arange(n_intervals + 1) * step
    return energy_mev * 1.0e-3, energy_mev


def _transport_matrix(full, masses, centre, omitted_left, omitted_right):
    """Use reciprocal FCs and retain only device contacts entering the lead self-energies."""
    reciprocal = (full + full.T) * 0.5
    discarded = {}
    for side, atoms in (("left", omitted_left), ("right", omitted_right)):
        columns = _dof_indices(atoms)
        contact = np.ix_(centre, columns)
        discarded[side] = float(np.max(np.abs(full[contact]), initial=0.0))
        reciprocal[contact] = 0.0
        reciprocal[np.ix_(columns, centre)] = 0.0
    diagnostics = asr_diagnostics(dynamical_matrix(reciprocal, masses, symmetrize=False),
                                  masses, threshold=MAX_ASR_RELATIVE, rows=centre)
    return diagnostics, discarded


def _build(fc_path, device_meta_path, geometry_path, thermal) -> dict:
    """把"力常数 → 完整矩阵 → 动力学矩阵 → 切块下标"这段装配抽出来。

    抽成函数不只是为了让 ``collect`` 短一点：``--sensitivity`` 要用**同一份**
    矩阵、**同一套**切块把 η/步长扫一遍。如果两边各写一份，就会重演 2026-11-24
    那个 bug（原子下标当自由度下标切 D，形状恰好撞上，T_ph 恒为 0 而检查全过）。
    """
    k_dev = np.load(fc_path)
    device_meta = json.loads(Path(device_meta_path).read_text(encoding="utf-8"))
    lead = load_lead_fc(_lead_fc_path(thermal.lead_fc))

    n_atoms = int(device_meta["n_atoms"])
    elements, positions, lattice = read_device_structure(Path(geometry_path))
    try:
        (left_outer, left_pl), (right_outer, right_pl) = lead_boundary(
            lead, device_meta, elements, lattice, thermal.fc.cutoff)
    except ValueError as exc:
        raise StageError(f"thermal_analyze: {exc}") from exc
    masses = masses_from_elements(elements)
    device_indices = [int(i) for i in device_meta["device_indices"]]
    lead_indices = [int(i) for i in device_meta["lead_indices"]]
    left_block = left_outer + left_pl
    right_block = right_outer + right_pl
    meta = FullFCMeta(n_atoms=n_atoms, device_indices=device_indices,
                      left_indices=left_block, right_indices=right_block)
    full = assemble_full_fc(k_dev, lead, meta)
    # D13:ASR 只在**器件区行**上判 —— lead 行来自被截断成两个主层的 lead
    # 力常数,它们与第三层及更远的耦合不在矩阵里,行和按构造不可能闭合
    # (实测 lead 块自身相对残差 0.6)。器件区的邻居集合是完整的,才可判。
    centre = _dof_indices(device_indices)
    residual, relative = acoustic_sum_residuals(full, rows=centre)
    sides = device_meta["side"]
    left_omitted = [i for i in lead_indices if sides[i] == "left" and i not in left_pl]
    right_omitted = [i for i in lead_indices if sides[i] == "right" and i not in right_pl]
    asr, discarded = _transport_matrix(full, masses, centre, left_omitted, right_omitted)
    dynamical = asr["dynamical"]
    imag = imaginary_modes(dynamical)
    # 虚频只报波数没有意义：把本征矢拆成「谁在动、动的是哪类键」。
    # 区域划分是判断物理软模 vs lead 截断假象的第一道筛子。
    molecule_indices = [int(i) for i in device_meta.get("molecule_indices") or []]
    molecule_set = set(molecule_indices)
    regions = {
        "molecule": molecule_indices,
        "device_au": [i for i in device_indices if i not in molecule_set],
        "lead": lead_indices,
    }
    soft_modes = describe_soft_modes(dynamical, masses=masses, elements=elements,
                                     positions=positions, regions=regions)

    # D 是 3N×3N,按**自由度**编号。用原子下标去切会得到一个形状"看起来
    # 差不多"的子块(例如 54 个原子的 54×54 恰好与 18 个原子的 3N=54
    # 撞上),transmission 的形状检查照样通过,但结果是垃圾(2026-11-24
    # 本地干跑时 T_ph 恒为 0 才发现)。
    arrived = dynamical[np.ix_(centre, centre)]
    expected = 3 * len(device_indices)
    if arrived.shape != (expected, expected):
        raise StageError(
            "中心区动力学矩阵形状 %r,期望 (%d,%d);请检查自由度下标换算"
            % (arrived.shape, expected, expected))

    return {
        "n_atoms": n_atoms,
        "elements": elements,
        "positions": positions,
        "lattice": device_meta.get("lattice"),
        "masses": masses,
        "device_indices": device_indices,
        "lead_indices": lead_indices,
        "left_pl": left_pl,
        "right_pl": right_pl,
        "left_outer": left_outer,
        "right_outer": right_outer,
        "centre": centre,
        "full": full,
        "dynamical": dynamical,
        "asr": asr,
        "residual": residual,
        "relative": relative,
        "discarded": discarded,
        "imag": imag,
        "soft_modes": soft_modes,
        "lead": lead,
        # 器件规格（层数、面内、截断）——收敛轴的"当前值"从这里读，
        # 免得在报告里再写一遍可能与 geometry 不符的常数。
        "spec": device_meta.get("spec") or {},
    }


def cutoff_sensitivity(built: dict, omega: np.ndarray, eta: float, *,
                       temperature: float, band_top: float,
                       cutoffs: tuple[float, ...] | None = None,
                       uncut_g_ph: float | None = None) -> dict:
    """电极力常数截断半径的敏感性（零机时：只切现有矩阵）。

    切**只切 lead 的 lead×lead 块**（``d_L``/``d_R`` 与 ``V_L``/``V_R``）。
    为什么不整体切：器件区的力来自周期 DFT 胞，长程耦合已经折进胞内，按距离
    再切等于凭空造一个不同的哈密顿量（实测整体切 5 Å 会把分子切断，``G_ph``
    掉到 1e-19，那是伪结果而不是敏感性）。真正带截断参数的是电极力常数
    —— ``gen_au_lead_fc.py`` 就是按 ``fc_cutoff_A`` 切的，所以这条轴只该切它。

    距离必须用**面内最小像**：lead 块本身虽只有两层，但它来自周期体相，器件坐标
    里 11.9 Å 的一对在体相里可能只隔 3.1 Å（nap-sh 的 (18,35) 就是这样）。用裸距
    会在名义 10 Å 处切掉真实近邻，造出 4.3% 的伪差异。
    """
    positions = np.asarray(built["positions"], dtype=float)
    # 切的是**动力学矩阵**（eV²），不是力常数：``_curve_factory`` 要的是前者，
    # 拿力常数当输入会让每个频率都跑到带顶以上，透射恒为 0（第一版就踩了这个）。
    # 零化的元素位置两者一一对应，D_ij = K_ij/√(m_i m_j)，所以切法等价。
    matrix = np.asarray(built["dynamical"], dtype=float)
    curve = _curve_factory(built, omega, eta)

    lead_dofs = _dof_indices(built["lead_indices"])
    lead_atoms = list(built["lead_indices"])
    lead_mask = _atom_distances(positions[lead_atoms], built.get("lattice"))

    if cutoffs is None:
        cutoffs = (4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0)
        nominal = getattr(built.get("lead"), "fc_cutoff_A", None)
        if nominal is not None:
            cutoffs = tuple(sorted(set((*cutoffs, float(nominal)))))

    # 组装矩阵里 lead 块最远的**非零**元素（面内最小像）有多远。名义截断
    # 10 Å 通常大于两层块的面内跨度，所以这个数只用来说明"这个块里到底有哪些
    # 距离在起作用"，不能拿它证明 10 Å 对**体相**力常数已经足够。
    lead_block_indices = np.ix_(lead_dofs, lead_dofs)
    block = matrix[lead_block_indices]
    n_lead = len(lead_atoms)
    coupled = block.reshape(n_lead, 3, n_lead, 3).any(axis=(1, 3))
    span = float(lead_mask[coupled].max(initial=0.0))

    def conductance(candidate):
        try:
            value = thermal_conductance(omega, curve(candidate), temperature,
                                        band_top_mev=band_top)
            return float(value), None
        except ValueError as exc:
            return None, str(exc)

    uncut = conductance(matrix) if uncut_g_ph is None else (float(uncut_g_ph), None)
    evaluated = {np.ones(np.count_nonzero(coupled), dtype=bool).tobytes(): uncut}
    rows = []
    for cutoff in cutoffs:
        keep = lead_mask <= float(cutoff)
        # 掩码只在原本为零的耦合上不同时，动力学矩阵和整条透射谱都相同。
        key = keep[coupled].tobytes()
        if key not in evaluated:
            truncated = matrix.copy()
            keep_dof = np.repeat(np.repeat(keep, 3, axis=0), 3, axis=1)
            truncated[lead_block_indices] = np.where(keep_dof, block, 0.0)
            evaluated[key] = conductance(truncated)
        g_ph, rejected = evaluated[key]
        rows.append({"cutoff_A": float(cutoff), "g_ph_pw_per_k": g_ph,
                     "rejected": rejected, "pairs_kept": int(keep.sum() - n_lead)})

    accepted = [(row["cutoff_A"], row["g_ph_pw_per_k"]) for row in rows
                if row["g_ph_pw_per_k"] is not None]
    spread = None
    if len(accepted) >= 2:
        values = [value for _, value in accepted]
        largest = max(values)
        spread = float((largest - min(values)) / largest) if largest else None
    return {
        "values_A": [float(c) for c in cutoffs],
        "rows": rows,
        "relative_spread": spread,
        "effective_span_A": span,
        "uncut_g_ph_pw_per_k": uncut[0],
        "uncut_rejected": uncut[1],
        "method": "只对 lead×lead 块按面内最小像距离置零后重算 NEGF（不碰 DFT）",
    }


def convergence_axes(built: dict, cutoff_block: dict) -> dict:
    """三条收敛轴的状态：一条已测（后处理可得），两条未测（需重跑 DFT）。

    为什么要把"未测"写进产物：读 ``thermal.json`` 的人只看到 13.7，不会知道
    lead 层数与器件区大小从没扫过。把它们连同所需作业写在同一份文件里，
    比写在文档里更难被忽略。
    """
    spec = built.get("spec") or {}
    layers = int(spec.get("electrode_layers") or 0)
    device_layers = int(spec.get("device_layers") or 0)
    # 生产截断半径来自 lead 力常数本身（gen_au_lead_fc.py 按它切的），
    # 不在 device spec 里。
    cutoff = float(getattr(built.get("lead"), "fc_cutoff_A", 0.0) or 0.0)
    span = float(cutoff_block.get("effective_span_A") or 0.0)
    # ``effective_span_A`` 用的是**面内最小像距离**（lead 块来自周期体相，
    # 裸坐标距离会把跨边界的近邻算成 11.9 Å，实测最小像只有 3.1 Å）。
    # 它 ≤ 名义截断，只说明"再按名义值切一刀不改变矩阵"，既不是收敛证据、
    # 也不是未收敛证据：超出名义截断的体相耦合本来就没进这份资产，后处理
    # 无法验证。所以三种状态要分开：
    #   span > cutoff 且按名义值切会改变结果 → converged=False（有确凿短差）
    #   span > cutoff 但影响 <1%            → converged=True
    #   span ≤ cutoff                       → converged=None（beyond unverified）
    covers = bool(cutoff and span <= cutoff)
    gap = None
    converged = None
    if cutoff and not covers:
        nominal = next((row for row in cutoff_block.get("rows", [])
                        if abs(row["cutoff_A"] - cutoff) < 1e-9), None)
        full = cutoff_block.get("uncut_g_ph_pw_per_k")
        value = nominal.get("g_ph_pw_per_k") if nominal else None
        if full is not None and value is not None:
            if full > 0:
                gap = abs(value - full) / full
                converged = gap <= 0.01
            else:
                converged = value == 0
                gap = 0.0 if converged else None
    return {
        "fc_cutoff": dict(cutoff_block, production_A=cutoff,
                          nominal_shortfall=gap,
                          nominal_covers_span=covers,
                          beyond_nominal_checked=False,
                          converged=converged),
        "lead_layers": {
            "status": "unmeasured",
            "current": layers,
            "candidates": [layers, layers + 2] if layers else [],
            "requires": "新建 %d 层版热导器件并重跑位移集（约 1 天/档）"
                        % (layers + 2 if layers else 0),
            "why_not_post_hoc": "层数改变参与 DFT 的原子，力常数矩阵必须重算",
        },
        "device_layers": {
            "status": "unmeasured",
            "current": device_layers,
            "candidates": [1, device_layers] if device_layers else [],
            "requires": "device_layers=1 时器件区约 36 原子、216 个位移"
                        "（比生产档少 1/3 机时），可复用同一 relax 几何",
            "why_not_post_hoc": "器件区决定哪些行来自 DFT、哪些来自冻结 lead 块",
        },
    }


def _curve_factory(built: dict, omega: np.ndarray, eta: float):
    """给一份动力学矩阵算整条 ``T_ph(ω)``（η 在这里定，便于扫参）。"""
    centre = built["centre"]
    left, right = (_dof_indices(built[key]) for key in ("left_pl", "right_pl"))
    left_outer, right_outer = (_dof_indices(built[key])
                               for key in ("left_outer", "right_outer"))
    pairs = [(centre, centre), (left, left), (right, right), (centre, left),
             (centre, right), (left, left_outer), (right, right_outer)]
    blocks = [np.ix_(rows, cols) for rows, cols in pairs]
    # Bound the working set of complex matrices used in the batched inversions.
    width = max(len(centre), len(left), len(right))
    batch_size = max(1, (32 * 1024 ** 2) // (8 * 16 * width ** 2))

    def curve(matrix) -> np.ndarray:
        dc, dl, dr, vl, vr, vl_lead, vr_lead = (matrix[index] for index in blocks)
        values = np.empty(len(omega))
        for start in range(0, len(omega), batch_size):
            batch = slice(start, start + batch_size)
            values[batch] = transmission(omega[batch], dc, dl, dr, vl, vr, eta ** 2,
                                         v_lead_l=vl_lead, v_lead_r=vr_lead)
        return values

    return curve


def sensitivity(fc_path, device_meta_path, geometry_path, thermal, *, temperature,
                etas: list[float] | None = None,
                steps: list[float] | None = None) -> dict:
    """用同一份力常数扫 η 与能量步长（不碰 DFT、不花机时）。

    为什么需要它：这条谱在 η=0.05 meV 下有法布里–珀罗针状结构，η 从 0.5 降到
    0.05 会让 ``G_ph`` 从 10.9 升到 13.7（+25%），而梯形积分的点数从 0.1 加到
    0.05 meV 只改变 1–2%。也就是说**η 是方法参数、步长只是求积精度**；不把这条
    曲线摆出来，横向比较就变成了"不同 η 的两次计算谁高谁低"。
    """
    built = _build(fc_path, device_meta_path, geometry_path, thermal)
    band_top = float(thermal.negf.band_top_meV)
    etas = [0.3, 0.2, 0.1, 0.05] if etas is None else list(etas)
    steps = [0.1, 0.05, 0.025] if steps is None else list(steps)

    rows = []
    for eta_mev in etas:
        eta = float(eta_mev) * 1.0e-3
        for step in steps:
            omega, _ = _energy_grid(band_top, float(step))
            curve = _curve_factory(built, omega, eta)
            entry = {"eta_meV": float(eta_mev), "energy_step_meV": float(step)}
            try:
                t_ph = curve(built["dynamical"])
                entry["g_ph_pw_per_k"] = float(
                    thermal_conductance(omega, t_ph, temperature,
                                        band_top_mev=band_top))
                entry["rejected"] = None
            except ValueError as exc:
                # 带顶闸门拒绝本身就是结果：η 太大时泄漏被采到，说明该 η 下
                # T_ph 已经不可信，不是"算错了"。
                entry["g_ph_pw_per_k"] = None
                entry["rejected"] = str(exc)
            rows.append(entry)

    by_eta: dict[float, dict[float, float | None]] = {}
    for row in rows:
        by_eta.setdefault(row["eta_meV"], {})[row["energy_step_meV"]] = \
            row["g_ph_pw_per_k"]

    # η → 0 外推。**不能用全部点做一次全局线性拟合**：G(η) 是饱和的
    # （实测 0.3/0.2/0.1/0.05 → 12.37/13.12/13.61/13.72），全局拟合被低端
    # 的陡坡拖出 13.67，反而**低于**最小 η 处的 13.72 —— 不合物理。
    # 渐近斜率只由**最小的两个 η**给出,且固定在**生产步长**那一行:不同步长
    # 之间本来就有 2% 的非单调摆动(13.449 vs 13.717),拿最细步长外推会把
    # 求积底噪当成 η 趋势。
    production_step = float(thermal.negf.energy_step_meV)
    if not by_eta or production_step not in next(iter(by_eta.values())):
        production_step = min(steps)
    extrapolated = None
    usable = sorted(((eta, values[production_step]) for eta, values in by_eta.items()
                     if values.get(production_step) is not None),
                    key=lambda item: item[0])
    if len(usable) >= 2:
        eta_small, g_small = usable[0]
        eta_next, g_next = usable[1]
        slope = (g_small - g_next) / (eta_small - eta_next)
        extrapolated = {
            "g_ph_pw_per_k": float(g_small + slope * (0.0 - eta_small)),
            "slope_pw_per_k_per_meV": float(slope),
            "n_points": 2,
            "at_energy_step_meV": production_step,
            "smallest_eta_meV": float(eta_small),
            "g_ph_at_smallest_eta_pw_per_k": float(g_small),
            "method": "生产步长下、两个最小 η 的渐近线性外推",
        }

    step_spread = {}
    for eta, values in by_eta.items():
        finite = [v for v in values.values() if v is not None]
        if len(finite) >= 2 and max(finite) > 0:
            step_spread[str(eta)] = float((max(finite) - min(finite)) / max(finite))

    return {
        "temperature_K": temperature,
        "band_top_meV": band_top,
        "etas_meV": [float(e) for e in etas],
        "energy_steps_meV": [float(s) for s in steps],
        "rows": rows,
        "eta_zero_extrapolation": extrapolated,
        "step_relative_spread": step_spread,
        "production": {"eta_meV": float(thermal.negf.eta_meV),
                       "energy_step_meV": float(thermal.negf.energy_step_meV)},
    }


@register
class ThermalAnalyzeStage(Stage):
    name = "thermal_analyze"
    depends = ("thermal_fc",)
    needs_submit = False
    input_hint = "消费 thermal_fc/force_constants.npy 与 thermal.lead_fc"
    outputs = ("thermal.json", "environment.json")

    def prepare(self, ctx, stage_dir: Path) -> Prepared:
        stage_dir = Path(stage_dir)
        thermal = ctx.config.thermal
        fc_path = ctx.artifact("thermal_fc", "force_constants.npy")
        if not Path(fc_path).is_file():
            raise StageBlocked(
                self.name, "缺 thermal_fc/force_constants.npy：请先跑 finite_fc 并 collect",
                [str(fc_path)])
        path = _lead_fc_path(thermal.lead_fc)
        if not path.is_file():
            raise StageBlocked(
                self.name,
                "thermal.lead_fc 指向的文件不存在: %s（若用 auto，"
                "请用 `smj lead-fc --show` 确认随包资产是否完整）" % path,
                [str(path)])
        try:
            lead = load_lead_fc(path)
            device_meta = json.loads(Path(ctx.artifact(
                "thermal_device", "device_metadata.json")).read_text(encoding="utf-8"))
            elements, _, lattice = read_device_structure(Path(ctx.artifact("thermal_fc", "geometry.in")))
            lead_boundary(lead, device_meta, elements, lattice, thermal.fc.cutoff)
        except (KeyError, IndexError, OSError, TypeError, ValueError) as exc:
            raise StageBlocked(self.name, "与 thermal.lead_fc 自检失败: %s" % exc,
                               [str(path)]) from exc
        return self.prepared(
            ctx, stage_dir,
            inputs=[Path(fc_path), path],
            notes=[
                "本阶段为纯本地后处理：读力常数、算透射、做 Landauer 积分，"
                "不提交任何作业",
                "温度 %.3f K；N_k=1（面内超胞的 Γ 点）"
                % thermal.resolved_temperature(ctx.config.transport.temperature),
                "电极力常数：%s（面内 %s，%d 层，sha256 %s）"
                % (path, "×".join(str(v) for v in lead.in_plane), lead.layers,
                   lead.sha256[:12]),
            ],
            commands=[],
            parent_stage="thermal_fc",
        )

    def collect(self, ctx, stage_dir: Path) -> Collected:
        stage_dir = Path(stage_dir)
        thermal = ctx.config.thermal
        temperature = thermal.resolved_temperature(ctx.config.transport.temperature)

        built = _build(
            ctx.artifact("thermal_fc", "force_constants.npy"),
            ctx.artifact("thermal_device", "device_metadata.json"),
            ctx.artifact("thermal_fc", "geometry.in"),
            thermal)
        n_atoms = built['n_atoms']
        device_indices = built['device_indices']
        lead_indices = built['lead_indices']
        full = built['full']
        dynamical = built['dynamical']
        asr = built['asr']
        residual = built['residual']
        relative = built['relative']
        imag = built['imag']
        soft_modes = built['soft_modes']

        negf = thermal.negf
        band_top = float(negf.band_top_meV)
        step = float(negf.energy_step_meV)
        omega, energy_mev = _energy_grid(band_top, step)
        eta = float(negf.eta_meV) * 1.0e-3
        transmission_curve = _curve_factory(built, omega, eta)

        t_ph = transmission_curve(dynamical)
        g_ph = thermal_conductance(omega, t_ph, temperature, band_top_mev=band_top)
        cumulative = cumulative_thermal_conductance(omega, t_ph, temperature)

        # 虚频对热导有没有影响：把全部负本征值投影掉再算一遍 G_ph。这是一次
        # 同一套切块的对照（多花一遍 NEGF 时间，不碰 DFT），比"虚频不可传播
        # 所以无所谓"这种断言可靠。
        #
        # 生产矩阵已经互易,直接以生产热导为基准。
        soft_impact = None
        g_symmetric = None
        g_removed = None
        eigenvalues, eigenvectors = np.linalg.eigh(dynamical)
        if (eigenvalues < 0.0).any():
            negative = eigenvalues < 0.0
            basis = eigenvectors[:, negative]
            cleaned = dynamical - basis @ np.diag(eigenvalues[negative]) @ basis.T
            g_symmetric = g_ph
            g_removed = thermal_conductance(
                omega, transmission_curve(cleaned), temperature,
                band_top_mev=band_top)
            if g_symmetric:
                soft_impact = (g_removed - g_symmetric) / g_symmetric
        # 谐波 Landauer 的 T_ph(ω) 与温度无关，换 T 只是换玻色因子 —— 这条
        # G_ph(T) 曲线**不额外花机时**。它回答两件事：300 K 的值离高温平台
        # 还有多远（差多少算"准"），以及低温端是否回到 G∝T 的量子极限。
        temperature_grid = np.array(
            [10.0, 20.0, 30.0, 50.0, 75.0, 100.0, 150.0,
             200.0, 250.0, 300.0, 400.0, 500.0, 600.0, 800.0])
        g_ph_by_temperature = [
            thermal_conductance(omega, t_ph, float(t), band_top_mev=band_top)
            for t in temperature_grid]

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
            # 判据只对器件区行成立：lead 行是截断块，行和按构造不闭合。
            "asr_region": "device_rows",
            "transport_matrix_protocol": "reciprocal_principal_layer_v1",
            "transport_asr": {key: asr[key] for key in (
                "residual_ev_per_A2", "relative_residual", "residual_after_ev_per_A2")},
            "discarded_device_lead_max_ev_per_A2": built["discarded"],
            # 投影应用于输运分区；判据保留完整器件行与输运分区的投影前残差。
            "asr_projected": bool(asr["projected"]),
            "asr_projection_reason": asr["projection_reason"],
            "asr_projection_norm": asr["projection_norm"],
            "asr_residual_after_ev_per_A2": asr["residual_after_ev_per_A2"],
            "imaginary_modes": imag,
            # 软模归属（modes.describe_soft_modes）：原子/区域/键三类权重。
            # 判据只看 imaginary_modes 的 count 与 worst_cm1，这里是解释材料。
            "soft_modes": soft_modes,
            # 把负本征值全部投影掉后 G_ph 的相对变化（None = 无虚频）。这是
            # 「虚频到底要不要紧」的定量答复，比定性论证可复核。两个基准值
            # 一并落盘，免得报告只能给一个百分数而看不出绝对量级。
            "soft_mode_impact_relative_g_ph": soft_impact,
            "g_ph_symmetric_pw_per_k": g_symmetric,
            "g_ph_soft_removed_pw_per_k": g_removed,
            # 免费的温度曲线（同一份 T_ph，只换玻色因子）。300 K 是否已到平台、
            # 低温是否 G∝T，都从这两列直接读出。
            "g_ph_by_temperature": {
                "temperature_K": [float(t) for t in temperature_grid],
                "g_ph_pw_per_k": [float(v) for v in g_ph_by_temperature],
            },
            "n_k_points": 1,
            "band_top_meV": band_top,
            "eta_meV": float(negf.eta_meV),
            "energy_step_meV": step,
            # 力常数的对称性偏差：判据 ph_fc_symmetry 用它（M5-T2）。
            # 中心差分的力应当对称；明显不对称说明力噪声大或正负没配对。
            "fc_asymmetry": _fc_asymmetry(full),
            "lead_fc_sha256": built["lead"].sha256,
            "right_lead_transform": "reflect_z",
            "spectrum": {
                "energy_meV": energy_mev.tolist(),
                "transmission": t_ph.tolist(),
                "cumulative_g_ph_pw_per_k": cumulative.tolist(),
            },
        }
        payload["zt"] = self._zt_payload(ctx, g_ph, temperature)
        # 三条收敛轴：截断半径可后处理（零机时），另两条只登记为未测并写清
        # 代价。它回答的是"13.7 这个数对哪几个开关最不敏感"。
        cutoffs = cutoff_sensitivity(built, omega, eta, temperature=temperature,
                                     band_top=band_top, uncut_g_ph=g_ph)
        payload["convergence"] = convergence_axes(built, cutoffs)
        # 环境清单：决定数值的输入与工具。单独落盘是因为它比 thermal.json
        # 长得快（species 目录有几十份文件），digest 才进主文件。
        from ..phonon.environment import collect as collect_environment
        from ..phonon.environment import digest as environment_digest
        environment = collect_environment(ctx.config)
        (stage_dir / "environment.json").write_text(
            json.dumps(environment, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8")
        payload["environment"] = {
            "schema": environment["schema"],
            "digest": environment_digest(environment),
            "file": "environment.json",
        }
        (stage_dir / "thermal.json").write_text(
            json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8")
        # 出图是演示品，不是判定依据：没有 matplotlib 时只记 note（同 junction_analyze）。
        plot_note = None
        try:
            from ..plotting import plot_thermal, plot_thermal_vs_temperature
            plotted = plot_thermal(stage_dir / "thermal.json",
                                   stage_dir / "thermal_spectrum.png",
                                   formats=("png", "pdf"))
            if not plotted.get("ok"):
                plot_note = "热导谱未出图：%s" % plotted.get("reason")
            curve = plot_thermal_vs_temperature(
                stage_dir / "thermal.json",
                stage_dir / "thermal_conductance_vs_temperature.png",
                formats=("png", "pdf"))
            if not curve.get("ok"):
                plot_note = "G_ph(T) 曲线未出图：%s" % curve.get("reason")
        except Exception as exc:  # noqa: BLE001 —— 出图失败不得影响数值结论
            plot_note = "热导出图失败：%s" % exc
        results = {"thermal": {
            "g_ph_pw_per_k": g_ph,
            "temperature_K": temperature,
            "acoustic_sum_relative_residual": relative,
            "n_imaginary": imag["count"],
            "transport_asr_relative_residual": asr["relative_residual"],
            "worst_imaginary_cm1": imag["worst_cm1"],
            "zt_total": payload["zt"].get("zt_total"),
        }}
        notes = [
            "G_ph(%.0f K) = %.3f pW/K（有效通道数 %.3f，占单通道上限的 %.1f%%）"
            % (temperature, g_ph, g_ph / gq_pw_per_k(temperature),
               100.0 * g_ph / gq_pw_per_k(temperature)),
            "完整器件行 ASR 残差 %.3e eV/Å²；输运分区残差 %.3e eV/Å²"
            "（投影前相对 %.3e，门限 %.0e）；虚频 %d 个（最差 %.1f cm⁻¹）"
            % (residual, asr["residual_ev_per_A2"], asr["relative_residual"],
               MAX_ASR_RELATIVE, imag["count"], imag["worst_cm1"]),
            "G_ph(%.0f K) = %.3f pW/K" % (temperature_grid[-1], g_ph_by_temperature[-1]),
        ]
        spread = cutoffs.get("relative_spread")
        span = cutoffs.get("effective_span_A")
        nominal = payload["convergence"]["fc_cutoff"].get("production_A")
        shortfall = payload["convergence"]["fc_cutoff"].get("nominal_shortfall")
        if spread is None:
            notes.append("截断扫描缺少可比较的非零热导，见 convergence.fc_cutoff.rows")
        else:
            notes.append(
                "收敛轴：力常数截断扫描的 G_ph 相对散布 %.1f%%；"
                "lead 层数与器件区大小仍未测（需重跑位移集）"
                % (100.0 * spread))
        if shortfall is not None and shortfall > 0.01:
            notes.append(
                "收敛轴警告：名义截断 %.1f Å，但组装矩阵里最远非零元素在 %.2f Å，"
                "按名义值截断后的热导与当前完整矩阵相差 %.1f%%"
                % (nominal or 0.0, span or 0.0, 100.0 * shortfall))
        if plot_note:
            notes.append(plot_note)
        if payload["zt"].get("zt_total") is not None:
            notes.append(
                "ZT_total = %.5f（电子 ZT %.5f，G_ph/k_e = %.3f）；仅作参考，不参与门控"
                % (payload["zt"]["zt_total"], payload["zt"]["zt_electronic"],
                   payload["zt"]["g_ph_over_k_e"]))
        return self.collected(
            ctx, stage_dir, artifacts=["thermal.json", "environment.json"],
            results=results, notes=notes)

    # -- 内部 -------------------------------------------------------------- #
    @staticmethod
    def _zt_payload(ctx, g_ph: float, temperature: float) -> dict:
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
        block = data.get("junction", data) if isinstance(data, dict) else {}
        try:
            electronic_temperature = float(block["temperature_K"])
        except (KeyError, TypeError, ValueError):
            return {"available": False, "reason": "results.json 缺有效电子温度 temperature_K"}
        if not np.isfinite(electronic_temperature) or electronic_temperature <= 0:
            return {"available": False, "reason": "results.json 电子温度无效"}
        if not np.isclose(electronic_temperature, temperature, atol=0, rtol=1e-10):
            return {"available": False, "reason": "电子与声子温度不一致",
                    "electronic_temperature_K": electronic_temperature,
                    "phonon_temperature_K": temperature}
        try:
            zt_e = float(block["ZT_electronic"])
            k_e = float(block["k_e_pW_per_K"])
        except (KeyError, TypeError, ValueError):
            return {"available": False, "reason": "results.json 缺 ZT_electronic/k_e"}
        if not np.isfinite([zt_e, k_e]).all() or k_e <= 0:
            return {"available": False, "reason": "电子 ZT 或 k_e 无效，无法定义 ZT_total"}
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
