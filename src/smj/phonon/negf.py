"""声子 NEGF:表面格林函数(decimation)与 Caroli 透射。

**公式口径**(方案 §3.4/§3.5),与本仓电子链的 Caroli 完全同构,只是把
``E - H`` 换成 ``ω² - D``、把电子自能换成声子自能:

    G^r(ω) = [ (ω² + iη) I - D_c - Σ_L - Σ_R ]^{-1}
    Σ_α    = V_{cα} g_α(ω) V_{αc}
    Γ_α    = i (Σ_α - Σ_α†)
    T_ph   = Tr[ G^r Γ_L G^a Γ_R ]

**单位**:``d_layer``/``v_coupling``/``eta_ev2`` 都是 **eV²**(质量加权动力学矩阵
的本征值平方根才是 eV),``omega_ev`` 是 eV。``eta_ev2`` 加在 ``ω²`` 上,必须是
eV² 量级 —— 配置里的 ``eta_mev`` 按 ``(eta_mev*1e-3)**2`` 换算(1 meV 对应
1e-6 eV²),不是把 meV 直接当 eV² 用。

**为什么 η 加在 ω² 上而不是 ω 上**:本仓的物理展宽来自电极自能(金属在 ω→0 有
声学支,Im Σ 不为零),η 只负责矩阵条件数。若写 ``(ω+iη)²``,ω→0 时虚部消失,
``T_ph(ω→0)`` 会变成 0/0;写成 ``ω²+iη`` 则在整个频段都稳定。
"""
from __future__ import annotations

import numpy as np

__all__ = ["surface_gf", "transmission"]


def surface_gf(
    omega_ev,
    d_layer,
    v_coupling,
    eta_ev2: float = 1.0e-9,
    *,
    max_iter: int = 1000,
    tol: float = 1.0e-12,
):
    """半无限电极的表面格林函数 ``g_α(ω)``,decimation 迭代(López-Sancho)。

    Parameters
    ----------
    omega_ev : float or array_like (Nω,)
        声子能量 ħω,单位 eV。**支持批量** —— 生产一次要 251 个频点,逐点循环做
        2000 次 27×27 求逆会慢一个数量级。
    d_layer : (n, n) array
        一个主层的动力学矩阵,eV²。n = 3×主层原子数(每个原子 3 个自由度)。
    v_coupling : (n, n) array
        相邻主层之间的耦合矩阵 V,eV²(非厄米本身,``V_{i,i+1}``)。
    eta_ev2 : float
        加在 ω² 上的虚部,eV²(见模块 docstring)。
    max_iter, tol : int, float
        迭代上限与收敛判据 ``‖V_new‖ < tol``。

    Returns
    -------
    (n, n) 或 (Nω, n, n) complex
        表面格林函数,``Im g < 0``(推迟传播子)。输入标量时返回 (n,n)。

    Notes
    -----
    迭代式::

        ε_s ← ε_s + α g β
        ε   ← ε   + α g β + β g α
        α   ← α g α ,  β ← β g β ,  g ← (z - ε)^{-1} ,  z = ω² + iη

    收敛后 ``g_s = (z - ε_s)^{-1}``。1D 单通道链在通带中部收敛是代数式的
    (误差 ~ 1/N²),所以 ``max_iter`` 默认给到 1000;真实三维电极因为耦合矩阵
    非对角元分散,通常几十步就到 1e-12。收敛不足由调用方的残差判据
    (``ph_lead_cell_match`` 一族)兜底,不在这里静默放宽。
    """
    omega = np.atleast_1d(np.asarray(omega_ev, dtype=float))
    # atleast_2d:标量(单通道链)按 1×1 矩阵处理,避免调用方为了一个数套两层括号
    d = np.atleast_2d(np.asarray(d_layer, dtype=complex))
    v = np.atleast_2d(np.asarray(v_coupling, dtype=complex))
    if d.ndim != 2 or d.shape[0] != d.shape[1]:
        raise ValueError("d_layer 必须是方阵,得到 %r" % (d.shape,))
    if v.shape != d.shape:
        raise ValueError(
            "v_coupling 形状 %r 与 d_layer %r 不一致" % (v.shape, d.shape))
    if eta_ev2 < 0:
        raise ValueError("eta_ev2 必须非负(推迟格林函数),得到 %r" % (eta_ev2,))

    n = d.shape[0]
    nw = omega.size
    eye = np.eye(n, dtype=complex)
    # z = ω² + iη,形状 (Nω,1,1) 便于与 (Nω,n,n) 广播
    z = (omega ** 2)[:, None, None] + 1j * eta_ev2 * np.ones((nw, 1, 1))

    eps_s = np.broadcast_to(d, (nw, n, n)).copy()
    eps = eps_s.copy()
    alpha = np.broadcast_to(v, (nw, n, n)).copy()
    beta = np.broadcast_to(v.conj().T, (nw, n, n)).copy()

    residual = np.inf
    for _ in range(int(max_iter)):
        # 注意必须写 z*eye - eps:z 是 (Nω,1,1),若直接 z-eps 会把标量广播成
        # **满矩阵**(而不是单位矩阵的平移),得到 [[z-d00, z],[z, z-d11]] 这种
        # 非对角垃圾。1×1 情形下两者恰好相同,所以标量链测试抓不到这个 bug。
        g = np.linalg.inv(z * eye - eps)
        agb = alpha @ g @ beta
        eps_s = eps_s + agb
        eps = eps + agb + beta @ g @ alpha
        alpha_next = alpha @ g @ alpha
        residual = float(np.max(np.abs(alpha_next)))
        if residual < tol:
            alpha = alpha_next
            break
        beta = beta @ g @ beta
        alpha = alpha_next

    # fail-closed:没收敛就报错,不静默返回一个「看起来像解」的数(本仓 F2 口径)。
    # 真实的 Au 电极通常 20–50 步到 1e-12,走到 max_iter 说明输入有问题(死带、
    # 主层切错、耦合矩阵非厄米等),必须让调用方看见。
    # 残差是整个 ω 批量上的最大值,无法归因到单个频点,故只报范围。
    if residual >= tol:
        raise RuntimeError(
            "decimation 未收敛:‖V‖ 残差 %.3e > tol %.1e(%d 步,%d 个频点,ω∈[%.6g, %.6g] eV);"
            "检查主层划分与耦合矩阵,或提高 max_iter"
            % (residual, tol, int(max_iter), nw, float(omega.min()), float(omega.max())))

    g_surf = np.linalg.inv(z * eye - eps_s)
    # 标量输入返回 (n,n);数组输入返回 (Nω,n,n) —— 与 np 的常规约定一致
    return g_surf[0] if np.ndim(omega_ev) == 0 else g_surf


def _center_batched(d_c, nk: int, nw: int, has_k_axis: bool):
    """中心区矩阵 → (nk, nω, n_c, n_c) 视图。

    两种口径不能混:N_k=1(本仓生产口径,D10)时第三维是**频率**,矩阵随 ω 变
    (``ω²I - D_c`` 里的 ω² 是标量平移,但 D_c 本身也可以逐频点给);显式给
    ``k_weights`` 时第三维是 **k**,不叠加频率轴。
    """
    arr = np.asarray(d_c, dtype=complex)
    if arr.ndim == 0:                      # 单自由度中心区(1D 链基准)
        arr = arr.reshape(1, 1)
    if has_k_axis:
        if arr.ndim == 2:
            arr = arr[None, ...]
        if arr.ndim != 3 or arr.shape[0] != nk:
            raise ValueError(
                "d_c 形状 %r 与 k_weights 长度 %d 不匹配(应为 (%d,n,n))"
                % (arr.shape, nk, nk))
        return np.broadcast_to(arr[:, None, :, :], (nk, nw) + arr.shape[1:])
    if arr.ndim == 2:
        arr = arr[None, ...]
    if arr.ndim != 3:
        raise ValueError(
            "d_c 形状 %r:单 k 时应为 (n,n) 或 (Nω,n,n)" % (arr.shape,))
    if arr.shape[0] not in (1, nw):
        raise ValueError(
            "d_c 的 ω 批量 %d 与 ω 网格 %d 不一致" % (arr.shape[0], nw))
    return np.broadcast_to(arr[None, ...], (1, nw) + arr.shape[1:])


def _lead_batched(d_lead, nk: int, name: str):
    """电极一层矩阵 → (nk, n, n)。电极矩阵与 ω 无关,不支持频率批次。"""
    arr = np.asarray(d_lead, dtype=complex)
    if arr.ndim == 0:
        arr = arr.reshape(1, 1)
    if arr.ndim == 2:
        return np.broadcast_to(arr[None, ...], (nk,) + arr.shape)
    if arr.ndim != 3 or arr.shape[0] != nk:
        raise ValueError(
            "%s 形状 %r:应为 (n,n),或给 k 轴时 (%d,n,n);电极矩阵不带频率批次"
            % (name, arr.shape, nk))
    return arr


def transmission(
    omega_ev,
    d_c,
    d_l,
    d_r,
    v_l,
    v_r,
    eta_ev2: float = 1.0e-9,
    *,
    v_lead_l=None,
    v_lead_r=None,
    k_weights=None,
    **surface_kw,
):
    """声子透射 ``T_ph(ω) = Tr[G^r Γ_L G^a Γ_R]``(Caroli),支持 k 权重归一化。

    Parameters
    ----------
    omega_ev : float or (Nω,) array
        声子能量 ħω,eV。
    d_c : (n_c,n_c)、(Nω,n_c,n_c) 或 (Nk,n_c,n_c)
        中心区(分子 + 两侧若干 Au 层)的动力学矩阵,**eV²**。单 k 时第三维是
        频率(逐频点给矩阵),显式给 ``k_weights`` 时第三维是 k。
    d_l, d_r : (n_l,n_l)
        一个**主层**的电极动力学矩阵,eV²(两侧可以不同 → 允许不对称的物理,
        但不产生整流,见 Notes)。
    v_l, v_r : (n_c,n_l) / (n_c,n_r)
        中心区 ↔ 电极表面层的耦合矩阵,eV²。
    eta_ev2 : float
        加在 ω² 上的虚部,eV²。
    v_lead_l, v_lead_r : (n_l,n_l) / (n_r,n_r), optional
        电极**内部**相邻主层的耦合。缺省等于 ``v_l``/``v_r`` —— 当中心区最外层
        与电极层在化学上一致(本方案的 Au(111) 板正是如此)时这是正确的;
        若两者不同(例如中心区用了弛豫后的层),必须显式给出。
    k_weights : (Nk,) array, optional
        横向 k 权重,必须归一化到 1(``Σ w_k = 1``)。缺省 = 单 k(``N_k=1``),
        对应 3×3 超胞的 Γ 点(方案 D10)。**求和不除 N_k 是这里最容易犯的错**:
        会把 ``T_ph(0)`` 放大成 ``3N_k``、把 ``G_ph`` 放大 ``N_k`` 倍。

    Returns
    -------
    float 或 (Nω,) ndarray
        ``T_ph``。矢量输入时逐频点。

    Notes
    -----
    **线性响应下 ``T_ph`` 对 L↔R 严格对称**(方案 A-9):交换两侧只是把
    ``Tr[G^r Γ_L G^a Γ_R]`` 里两个 Γ 的位置对调,由迹的循环性结果不变。
    所以本函数**给不出热整流**;要方向性必须上非谐 NEGF 或带温差的 NEMD。
    两侧不对称只改变**绝对值**(``Σ_L ≠ Σ_R``),不产生方向依赖。
    """
    omega = np.atleast_1d(np.asarray(omega_ev, dtype=float))
    nw = omega.size

    has_k_axis = k_weights is not None
    if not has_k_axis:
        nk = 1
        weights = np.ones(1)
    else:
        weights = np.asarray(k_weights, dtype=float).reshape(-1)
        nk = weights.size
        if abs(weights.sum() - 1.0) > 1.0e-9:
            raise ValueError(
                "k_weights 必须归一化到 1(Σw=%.6f);求和而不除 N_k 会放大 T_ph"
                % (weights.sum(),))

    dc = _center_batched(d_c, nk, nw, has_k_axis)
    dl = _lead_batched(d_l, nk, "d_l")
    dr = _lead_batched(d_r, nk, "d_r")
    nc = dc.shape[-1]
    if dc.shape[-2] != nc:
        raise ValueError("d_c 必须是方阵,得到 %r" % (dc.shape,))
    for name, arr in (("d_l", dl), ("d_r", dr)):
        if arr.shape[1] != arr.shape[2]:
            raise ValueError("%s 必须是方阵" % (name,))

    vl = np.atleast_2d(np.asarray(v_l, dtype=complex))
    vr = np.atleast_2d(np.asarray(v_r, dtype=complex))
    if vl.shape != (nc, dl.shape[1]):
        raise ValueError("v_l 形状 %r,期望 (%d,%d)" % (vl.shape, nc, dl.shape[1]))
    if vr.shape != (nc, dr.shape[1]):
        raise ValueError("v_r 形状 %r,期望 (%d,%d)" % (vr.shape, nc, dr.shape[1]))

    vl_lead = np.atleast_2d(np.asarray(vl if v_lead_l is None else v_lead_l, dtype=complex))
    vr_lead = np.atleast_2d(np.asarray(vr if v_lead_r is None else v_lead_r, dtype=complex))
    if vl_lead.shape != (dl.shape[1], dl.shape[1]):
        raise ValueError("v_lead_l 形状 %r,期望 (%d,%d)" % (vl_lead.shape, dl.shape[1], dl.shape[1]))
    if vr_lead.shape != (dr.shape[1], dr.shape[1]):
        raise ValueError("v_lead_r 形状 %r,期望 (%d,%d)" % (vr_lead.shape, dr.shape[1], dr.shape[1]))

    eye_c = np.eye(nc, dtype=complex)
    z = (omega ** 2)[:, None, None] + 1j * eta_ev2 * np.ones((nw, 1, 1))

    out = np.zeros((nw, nk), dtype=float)
    for ik in range(nk):
        g_l = surface_gf(omega, dl[ik], vl_lead, eta_ev2, **surface_kw)
        g_r = surface_gf(omega, dr[ik], vr_lead, eta_ev2, **surface_kw)
        sigma_l = vl @ g_l @ vl.conj().T
        sigma_r = vr @ g_r @ vr.conj().T
        gamma_l = 1j * (sigma_l - sigma_l.conj().transpose(0, 2, 1))
        gamma_r = 1j * (sigma_r - sigma_r.conj().transpose(0, 2, 1))
        g_c = np.linalg.inv(z * eye_c - dc[ik] - sigma_l - sigma_r)
        g_c_dag = g_c.conj().transpose(0, 2, 1)
        t_k = np.real(np.trace(g_c @ gamma_l @ g_c_dag @ gamma_r, axis1=-2, axis2=-1))
        out[:, ik] = t_k

    t_avg = out @ weights
    return float(t_avg[0]) if np.ndim(omega_ev) == 0 else t_avg
