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
from ..phonon.dynamical import (acoustic_sum_relative_residual,
                               acoustic_sum_residual, dynamical_matrix,
                               imaginary_modes, masses_from_elements)
from ..phonon.checks import MAX_ASR_RELATIVE
from ..phonon.fc import (FullFCMeta, asr_diagnostics, assemble_full_fc,
                         load_lead_fc)
from ..phonon.landauer import cumulative_thermal_conductance, thermal_conductance
from ..phonon.modes import describe_soft_modes
from ..phonon.negf import transmission
from ..phonon.units import gq_pw_per_k
from .base import Collected, Prepared, Stage, register

SCHEMA = "smj.thermal.v1"


def _fc_asymmetry(force_constants) -> float:
    """`||K - K^T|| / ||K||`（Frobenius）。供判据 `ph_fc_symmetry`。"""
    import numpy as np

    matrix = np.asarray(force_constants, dtype=float)
    scale = float(np.linalg.norm(matrix))
    if scale <= 0.0:
        return 0.0
    return float(np.linalg.norm(matrix - matrix.T) / scale)



def _lead_fc_path(value: str | Path | None) -> Path:
    """``thermal.lead_fc`` 可以写成目录(用其中的 ``lead_fc.json``)或直接写文件。

    ``save_lead_fc`` 落盘的是 ``lead_fc.npy`` + ``lead_fc.json`` 两个文件,
    用户按目录配置更自然,这里两种都接受。

    ``None`` / ``"auto"`` 解析到**随包发布的 Au(111) 资产**（``smj.phonon.assets``）：
    电极是金、与分子无关，那份力常数只需算一次，不该让每个使用者先搭一套
    体相计算才能跑分子。
    """
    if value is None or str(value) == "auto":
        from ..phonon.assets import lead_fc_dir
        return lead_fc_dir() / "lead_fc.json"
    path = Path(value)
    if path.is_dir():
        candidate = path / "lead_fc.json"
        if candidate.is_file():
            return candidate
    return path


def _structure(geometry: Path) -> tuple[list[str], np.ndarray]:
    """读 FHI-aims ``geometry.in`` 的元素序列与坐标（Å，只取 ``atom`` 行）。"""
    elements: list[str] = []
    coords: list[list[float]] = []
    for line in Path(geometry).read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) >= 5 and parts[0] == "atom":
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
            elements.append(parts[4])
    return elements, np.asarray(coords, dtype=float)


def _dof_indices(atoms) -> np.ndarray:
    """原子下标 → 自由度下标。

    ``dynamical_matrix`` 返回的是 3N×3N 矩阵,第 ``3*a+k`` 行对应原子 ``a``
    的第 ``k`` 个方向。整条热导链上凡是"取某个原子子集的块",都必须先过这道
    换算;直接拿原子下标去 ``np.ix_`` 会静默切出错位的块(见 collect 里的注释)。
    """
    return np.array([3 * int(a) + k for a in atoms for k in range(3)], dtype=int)


def _energy_grid(band_top: float, step: float) -> tuple[np.ndarray, np.ndarray]:
    """能量网格：到 ``1.25 × 带顶``，给闸门留出"带顶以上必须为零"的取样点。"""
    n_points = int(round(band_top / step)) + 1
    energy_mev = np.linspace(0.0, band_top * 1.25, max(n_points, 2))
    return energy_mev * 1.0e-3, energy_mev


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
    elements, positions = _structure(Path(geometry_path))
    if len(elements) != n_atoms:
        raise StageError(
            "thermal_analyze 的元素数(%d)与 device_metadata 的 n_atoms(%d) 不一致"
            % (len(elements), n_atoms))
    masses = masses_from_elements(elements)

    device_indices = [int(i) for i in device_meta["device_indices"]]
    lead_indices = [int(i) for i in device_meta["lead_indices"]]
    layers = device_meta["layer_index"]
    sides = device_meta["side"]
    # 主层(principal layer)= 紧挨器件区的那一层;它外面还有一层用来给
    # decimation 提供层间耦合 V_lead。**这两层**才是 lead 力常数要覆盖的
    # 范围 —— lead_fc 存的是「两个主层」的块(18 原子 × 2 层 = 54 维),
    # 而不是器件侧全部 4 层冻结 lead(36 原子)。传全部 36 个会让
    # _place_lead_block 直接报「lead 区原子数 36 与 lead 力常数原子数 18 不符」。
    left_pl = [i for i in lead_indices if layers[i] == 1 and sides[i] == "left"]
    right_pl = [i for i in lead_indices if layers[i] == 1 and sides[i] == "right"]
    left_outer = [i for i in lead_indices if layers[i] == 0 and sides[i] == "left"]
    right_outer = [i for i in lead_indices if layers[i] == 0 and sides[i] == "right"]
    if not (left_pl and right_pl and left_outer and right_outer):
        raise StageError(
            "thermal_analyze 找不到两侧的主层或最外层：layers=%s 时每侧至少需要 2 层"
            % sorted(set(layers)))
    left_block = left_pl + left_outer
    right_block = right_pl + right_outer
    if len(left_block) != lead.n_atoms or len(right_block) != lead.n_atoms:
        raise StageError(
            "lead 力常数覆盖 %d 个原子，但器件侧的主层+次外层是 %d 个（左）/%d 个（右）："
            "两者必须一致，否则拼不出完整矩阵"
            % (lead.n_atoms, len(left_block), len(right_block)))
    meta = FullFCMeta(n_atoms=n_atoms, device_indices=device_indices,
                      left_indices=left_block, right_indices=right_block)
    full = assemble_full_fc(k_dev, lead, meta)
    # D13:ASR 只在**器件区行**上判 —— lead 行来自被截断成两个主层的 lead
    # 力常数,它们与第三层及更远的耦合不在矩阵里,行和按构造不可能闭合
    # (实测 lead 块自身相对残差 0.6)。器件区的邻居集合是完整的,才可判。
    centre = _dof_indices(device_indices)
    residual = acoustic_sum_residual(full, rows=centre)
    relative = acoustic_sum_relative_residual(full, rows=centre)
    # symmetrize=False:器件↔lead 交叉块已由 assemble_full_fc 补了转置,
    # 再走一次对称化只会把 lead 侧那两块按构造不闭合的行平均进器件行,
    # 凭空造出 0.1% 的"合规则违反"(实测 2.6e-13 → 7.5e-2 eV/Å²)。
    asr = asr_diagnostics(dynamical_matrix(full, masses, symmetrize=False),
                          masses, threshold=MAX_ASR_RELATIVE, rows=centre)
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
        "imag": imag,
        "soft_modes": soft_modes,
        "lead": lead,
    }


def _curve_factory(built: dict, omega: np.ndarray, eta: float):
    """给一份动力学矩阵算整条 ``T_ph(ω)``（η 在这里定，便于扫参）。"""
    centre = built["centre"]
    dofs = _dof_indices
    left_pl, right_pl = built["left_pl"], built["right_pl"]
    left_outer, right_outer = built["left_outer"], built["right_outer"]

    def curve(matrix) -> np.ndarray:
        return np.array([
            transmission(float(w),
                         matrix[np.ix_(centre, centre)],
                         matrix[np.ix_(dofs(left_pl), dofs(left_pl))],
                         matrix[np.ix_(dofs(right_pl), dofs(right_pl))],
                         matrix[np.ix_(centre, dofs(left_pl))],
                         matrix[np.ix_(centre, dofs(right_pl))],
                         eta ** 2,
                         v_lead_l=matrix[np.ix_(dofs(left_pl), dofs(left_outer))],
                         v_lead_r=matrix[np.ix_(dofs(right_pl), dofs(right_outer))])
            for w in omega])

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
    outputs = ("thermal.json",)

    def prepare(self, ctx, stage_dir: Path) -> Prepared:
        stage_dir = Path(stage_dir)
        thermal = ctx.config.thermal
        fc_path = ctx.artifact("thermal_fc", "force_constants.npy")
        if not Path(fc_path).is_file():
            raise StageBlocked(
                self.name, "缺 thermal_fc/force_constants.npy：请先跑 finite_fc 并 collect",
                [str(fc_path)])
        if thermal.lead_fc is None:
            raise StageBlocked(
                self.name,
                "配置缺 thermal.lead_fc：声子 NEGF 需要冻结 lead 力常数（方案 D4）。"
                "把 thermal.lead_fc 写成 auto（用随包资产），或用 "
                "`smj lead-fc --export DIR` 导出后指向该目录。",
                [])
        path = _lead_fc_path(thermal.lead_fc)
        if not path.is_file():
            raise StageBlocked(
                self.name,
                "thermal.lead_fc 指向的文件不存在: %s（若用 auto，"
                "请用 `smj lead-fc --show` 确认随包资产是否完整）" % path,
                [str(path)])
        from ..phonon.assets import lead_fc_info
        try:
            lead_info = lead_fc_info(path.parent)
        except (OSError, ValueError) as exc:
            raise StageBlocked(self.name, "与 thermal.lead_fc 自检失败: %s" % exc,
                               [str(path)]) from exc
        return self.prepared(
            ctx, stage_dir,
            inputs=[Path(fc_path), path],
            notes=[
                "本阶段为纯本地后处理：读力常数、算透射、做 Landauer 积分，"
                "不提交任何作业",
                "温度 %.3f K；N_k=1（3×3 超胞的 Γ，Γ 已把 9 个原胞 q 点折叠进来）"
                % thermal.resolved_temperature(ctx.config.transport.temperature),
                "电极力常数：%s（面内 %s，%d 层，sha256 %s）"
                % (lead_info["path"], "×".join(str(v) for v in
                                               lead_info.get("in_plane", [])),
                   int(lead_info.get("layers", 0)),
                   str(lead_info.get("sha256", ""))[:12]),
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
        left_pl, right_pl = built['left_pl'], built['right_pl']
        left_outer, right_outer = built['left_outer'], built['right_outer']
        centre = built['centre']
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
        # 基准取**对称化**后的矩阵：本征分解要求厄米，而有限差分给的 K 不是
        # 严格对称的（互易性残差 ~0.5%）。若拿非对称 D 当基准，这个残差会混
        # 进 "+1.1%"，看起来像虚频的影响（实测对称化后是 +0.58%）。
        soft_impact = None
        g_symmetric = None
        g_removed = None
        symmetric = 0.5 * (dynamical + dynamical.T)
        eigenvalues = np.linalg.eigvalsh(symmetric)
        if (eigenvalues < 0.0).any():
            _, eigenvectors = np.linalg.eigh(symmetric)
            negative = eigenvalues < 0.0
            basis = eigenvectors[:, negative]
            cleaned = symmetric - basis @ np.diag(eigenvalues[negative]) @ basis.T
            g_symmetric = thermal_conductance(
                omega, transmission_curve(symmetric), temperature,
                band_top_mev=band_top)
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
            # D13：投影是"超阈值才做"，所以这三个字段要能表达"没投影"。
            # 判据 ph_asr_residual 判的是上面那两个**投影前**的值。
            "asr_projected": bool(asr["projected"]),
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
        # 跨 run 比较时要确认用的是同一份电极力常数(smj thermal-compare 会查)。
        "lead_fc_sha256": built["lead"].sha256,
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
            "worst_imaginary_cm1": imag["worst_cm1"],
            "zt_total": payload["zt"].get("zt_total"),
        }}
        notes = [
            "G_ph(%.0f K) = %.3f pW/K（有效通道数 %.3f，占单通道上限的 %.1f%%）"
            % (temperature, g_ph, g_ph / gq_pw_per_k(temperature),
               100.0 * g_ph / gq_pw_per_k(temperature)),
            "ASR 残差 %.3e eV/Å²（相对 %.3e）；虚频 %d 个（最差 %.1f cm⁻¹）"
            % (residual, relative, imag["count"], imag["worst_cm1"]),
            "G_ph(T)：%.0f K 时是高温平台（%.3f pW/K）的 %.1f%%"
            % (temperature, g_ph_by_temperature[-1],
               100.0 * g_ph / g_ph_by_temperature[-1]),
        ]
        if plot_note:
            notes.append(plot_note)
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
