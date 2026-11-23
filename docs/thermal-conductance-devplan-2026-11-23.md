# 单分子结热导（`G_ph`）开发计划（2026-11-23）

> **配套文档**
> - 方案（设计/物理/决策）：[`thermal-conductance-plan-2026-11-23.md`](thermal-conductance-plan-2026-11-23.md)
> - 文献调研（三支技术栈）：[`thermal-conductance-literature-2026-11-23.md`](thermal-conductance-literature-2026-11-23.md)
>
> 本文只回答「怎么做、按什么顺序、每步怎么算做完」。**不重复方案的物理推导**，引用处直接给 `§x` 或 `T1`。

**基线（2026-11-23 实测）**：`git 608368e`；`python -m pytest -q` → **620 passed, 2 deselected**（22.7 s）；Python 3.12.10；集群 FHI-aims **221103**（无 DFPT，`atomic_pert_response` 0 命中）、Slurm `MaxArraySize=1001`、单节点 64 核 / 124.7 GB / 4×RTX 4090。

---

## 0. 一页速览

| 项 | 值 |
|---|---|
| 任务数 | **30 个**（M1 6 / M2 6 / M3 7 / M4 6 / M5 5） |
| 新增源码 | 约 **1,450 行**（`src/smj/phonon/` 7 个模块 + 1 个远程脚本 + 4 个阶段 + 配置 + 判据 + 报告） |
| 新增测试 | 约 **900 行**（`tests/phonon/` 10 个文件 + 1 个夹具 + 1 个集成用例） |
| 内核改动 | **1 处、≤15 行**（`write_aims_geometry(lattice=)`），需刷新 provenance + demo 对账 |
| 人日 | **18.75 人日**（不含等待机时；关键路径 14.5） |
| 机时 | **约 52–93 h 墙钟**（8 路并发，含开发期收敛矩阵；见 §6.2） |
| 关键路径 | M1 → M2 → M3 → M4 → M5（M1 与 M2 可并行） |
| 电子链影响 | **零**。`thermal.enabled` 默认 `false`，`STAGE_ORDER` 不变，`smj run` 逐字节不变 |

**顺序原则**：先把**不碰 DFT 的数值内核**（M1）做完并用解析解钉死，再做几何与 IO（M2）——这样后面所有物理数字都有一个独立的对照物。反过来先跑 DFT 会得到一个「看起来对」的数，但没有办法判断它对不对。

---

## 1. 交付物与文件清单

### 1.1 新增

| 文件 | 行数估 | 内容 | 首次出现 |
|---|---|---|---|
| `src/smj/phonon/__init__.py` | 15 | 导出面 + 版本常量 | M1-T1 |
| `src/smj/phonon/units.py` | 70 | 唯一换算常量、`gq_pw_per_k` | M1-T1 |
| `src/smj/phonon/negf.py` | 210 | `surface_gf`（decimation）、`transmission` | M1-T2/T3 |
| `src/smj/phonon/landauer.py` | 130 | `thermal_conductance`（唯一公式源） | M1-T4 |
| `src/smj/phonon/dynamical.py` | 100 | `dynamical_matrix`、`imaginary_modes` | M1-T5 |
| `src/smj/phonon/forces.py` | 95 | `read_forces` | M2-T3 |
| `src/smj/phonon/fc.py` | 260 | 位移生成、FC 组装、ASR、lead FC 装载与拼接 | M2-T4 |
| `src/smj/phonon/device.py` | 380 | 周期板构造、角色/层号、S–H 去除 | M2-T5 |
| `src/smj/phonon/extract_forces.py` | 60 | **远程**用（仅 stdlib），阵列作业里把 `aims.out` → `forces.json` | M2-T6 |
| `src/smj/junction/templates/control_thermal_relax.in` | 40 | 弛豫模板 | M2-T2 |
| `src/smj/junction/templates/control_thermal_force.in` | 30 | 单点力模板 | M2-T2 |
| `src/smj/stages/thermal_device.py` | 150 | 阶段 | M3-T2 |
| `src/smj/stages/thermal_relax.py` | 140 | 阶段 | M3-T3 |
| `src/smj/stages/thermal_fc.py` | 260 | 阶段（数组作业） | M3-T4 |
| `src/smj/stages/thermal_negf.py` | 180 | 阶段（纯后处理） | M3-T5 |
| `src/smj/junction/data/au111_3x3_l4.json` | 数据 | 冻结 lead 力常数（约 200 KB） | M4-T2 |
| `scripts/gen_au_lead_fc.py` | 200 | 一次性生成 lead FC | M4-T2 |
| `tests/phonon/`（10 个文件） | 约 900 | 见 §2 | M1…M5 |
| `tests/data/aims_force_block_20atoms.out` | 约 40 行 | 力解析夹具 | M2-T3 |
| `tests/test_thermal_integration.py` | 120 | `@pytest.mark.integration` | M5-T5 |

### 1.2 修改

| 文件 | 改动 | 风险 |
|---|---|---|
| `src/smj/junction/geometry.py:68` | `write_aims_geometry(path, atoms, constrained=(), lattice=None)` | **唯一内核改动**；`lattice=None` 必须逐字节不变 |
| `src/smj/junction/tcontrol.py:102` | `make_control(mode=...)` 增加 `thermal_relax`/`thermal_force` 两个 mode | 低，纯分支 |
| `src/smj/config.py:469/485/718` | `_KNOWN_SECTIONS` + `_THERMAL_KEYS` + `ThermalConfig` + 解析 | 低 |
| `src/smj/state.py:41-48` | 新增 `THERMAL_STAGES`，**`STAGE_ORDER` 不动** | 低 |
| `src/smj/stages/__init__.py` | import 4 个新阶段完成注册 | 低 |
| `src/smj/junction/criteria.py` | 10 条判据 + 4 条 `physical` 的 `source` | 中：`lint_registry()` 会校验 |
| `src/smj/report.py:1010` 之后 | 新增 `## 5. 热导` 节 + `_thermal_summary` | 低 |
| `src/smj/cli.py` | 新增 `thermal` 子命令（照 `cmd_stage`:637） | 中：argparse 冲突 |
| `pyproject.toml` | 无需改（numpy 已在依赖里） | — |

---

## 2. 任务分解

每个任务五段：**目标 / 文件 / 实现要点 / 测试 / DoD**。估时为「专注人日」。

### 2.1 M1：数值内核（纯离线，不碰 DFT）

> 门禁：**T1–T4 全绿**。这一批做完之前不要碰集群。

#### M1-T1 `units.py`：唯一换算常量

- **文件**：`src/smj/phonon/units.py`、`tests/phonon/test_units.py`
- **实现要点**
  ```python
  EV_PER_SQRT: Final[float] = 0.0646542      # ħω[eV] = EV_PER_SQRT·√λ, λ in eV/(Å²·amu)
  THZ_PER_SQRT: Final[float] = 15.6333       # ν[THz] = THZ_PER_SQRT·√λ
  EV_TO_CM1: Final[float] = 8065.544         # 1 eV = 8065.544 cm⁻¹
  KB_EV_PER_K: Final[float] = 8.617333262e-5
  KB_J_PER_K: Final[float] = 1.380649e-23
  H_J_S: Final[float] = 6.62607015e-34
  H_EV_S: Final[float] = 4.135667696e-15
  E_CHARGE: Final[float] = 1.602176634e-19

  def lambda_to_ev(lam) -> np.ndarray          # sqrt(lam)*EV_PER_SQRT
  def ev_to_cm1(ev) -> np.ndarray             # ev*EV_TO_CM1
  def gq_pw_per_k(temperature_K: float) -> float   # π²k_B²T/(3h)·1e12
  ```
- **易错点**：`EV_TO_CM1` 必须是 **8065.544**（= 1/1.239841984e-4），不是 8065.73（旧 CODATA）。`gq_pw_per_k` 必须用 **SI 的 `k_B`/`h`**，再用 `1e12` 换成 pW——用 eV 版的 `k_B` 会差 `E_CHARGE` 倍，这是本方案第三类单位陷阱（前两类见附录 A-1/A-2）。
- **测试**（`test_units.py`，6 条）
  1. `EV_PER_SQRT**2 * 1e20` 与 CODATA 独立推导一致（相对 1e-12）；
  2. `THZ_PER_SQRT` 与 `EV_PER_SQRT / H_EV_S / 1e12` 自洽，且与 phonopy `VaspToTHz=15.633302` 相对误差 < 1e-6；
  3. `lambda_to_ev(0.1) = 20.445 meV ± 1e-3`（≈ Au 带顶 19.5 meV）；
  4. C–C 双原子地标：`k=30/50 eV/Å²`、`m=12.011` → **824 / 1064 cm⁻¹**（±2）；
  5. `gq_pw_per_k(300) = 283.93 ± 0.01`、`(100)=94.64`、`(4)=3.786`；
  6. `ev_to_cm1(0.1)=806.554 ± 1e-2`。
- **DoD**：T1 绿；模块无 numpy 以外依赖；所有常量带 `Final` 与推导注释。

#### M1-T2 `negf.surface_gf`：decimation

- **文件**：`src/smj/phonon/negf.py`、`tests/phonon/test_negf_analytic.py`（与 T3/T4 共用）
- **签名**
  ```python
  def surface_gf(omega_ev, d_layer, v_coupling, eta_ev,
                 *, max_iter: int = 500, tol: float = 1e-12) -> np.ndarray:
      """返回半无限 lead 的表面格林函数 g_surf(ω)，形状 (n,n)，n = d_layer 的维数。
      d_layer: 主层动力学矩阵，单位 eV²；v_coupling: 相邻主层耦合，单位 eV²。"""
  ```
- **实现要点**：标准 López-Sancho 迭代 `g ← (ω²I + iη − d − v·g·v†)⁻¹`，收敛判据 `‖g_new − g_old‖_F < tol`；**必须支持批量 ω**（输入 `(Nω, n, n)`，避免 251 次 Python 循环）。`η` 用 `iη` 加在 `ω²` 上。
- **测试**
  1. 一维单原子链（`d=2k`、`v=−k`）解析：`g_surf = (−1/√(4k²−(ω²+iη)²))` 或等价形式，相对误差 < 1e-10；
  2. **T4**：200 层有限堆叠直接求逆 vs decimation，`g_surf` 相对误差 < 1e-8（取 η=1e-3 抑制病态）；
  3. 批量与逐点结果逐元素相等（< 1e-14）。
- **DoD**：T4 绿；`omega_ev` 标量/数组两种输入都工作；不收敛时 raise（不静默返回）。

#### M1-T3 `negf.transmission`：Caroli 公式 + **`1/N_k`**

- **文件**：`src/smj/phonon/negf.py`
- **签名**
  ```python
  def transmission(omega_ev, d_c, d_l, d_r, v_l, v_r, eta_ev,
                   *, k_weights: np.ndarray | None = None) -> np.ndarray:
      """T_ph(ω) = (1/N_k) Σ_k w_k · Tr[G^r Γ_L G^a Γ_R]。
      k_weights 为 None 时等价于单个 k 点（w=1）。Σ_k w_k 必须为 1，否则 raise。"""
  ```
- **实现要点**
  - `G^r = inv((ω²+iη)I − d_c − Σ_L − Σ_R)`；`Σ_α = v_α g_α v_α†`；`Γ_α = i(Σ_α − Σ_α†)`；
  - **`(ω²)` 不是 `ω`**：`d_c` 必须已是 **eV²** 量纲（由 `dynamical_matrix` 负责换算），函数 docstring 明写这一前置条件；
  - `np.real(trace(G @ Γ_L @ G.conj().T @ Γ_R))`；结果裁掉 `|T| < 1e-12` 的负噪声；
  - **`k_weights` 校验**：`abs(sum - 1) > 1e-9` → raise；`len(w) != n_k` → raise。这是 §3.4 警示框的代码化。
- **测试**
  1. 单原子链通带内 `T_ph=1`、带外 `0`（**T2**），与解析色散逐点吻合 < 1e-10；
  2. `k_weights=[1.0]` 与 `None` 结果相同；`k_weights=[1/3,1/3,1/3]` 与单点结果相同（同构重复）；`k_weights=[1,1,1]` → **raise**；
  3. 对称性：交换 L/R 后 `T_ph` 不变（< 1e-12）——这条同时验证了「线性响应无整流」（附录 A-9）；
  4. `T_ph` 对 η 的连续性：η=0.5/1.0 meV 下同一点差 < 1e-3。
- **DoD**：T2 绿；`k_weights` 归一化有 raise 用例。

#### M1-T4 `landauer.thermal_conductance`：唯一公式源

- **文件**：`src/smj/phonon/landauer.py`、`tests/phonon/test_landauer.py`
- **签名**
  ```python
  def thermal_conductance(omega_ev, t_ph, temperature_K: float,
                          *, band_top_meV: float = 19.5) -> dict:
      """返回 {G_ph_pW_per_K, GQ_pW_per_K, N_eff, g_cum_pW_per_K,
               t_ph_at_zero, t_ph_above_band_max, x_max}。"""
  ```
- **实现要点**
  - `x = ħω/(k_B T)`（`k_B` 用 **eV/K**）；核 `x²eˣ/(eˣ−1)²`；`x→0` 时用极限 1.0 避免 0/0（`np.where(x < 1e-8, 1.0, ...)`）；
  - `G = (k_B_eV²·T/h_eV_s) · E_CHARGE · ∫...` → **W/K**，再 `×1e12` 得 pW/K；
  - `g_cum` 用累积梯形积分；
  - fail-closed：`t_ph[ω > band_top] > 1e-6` → raise；`t_ph[0] > 3` → raise；`G > 3·G_Q` → raise（三条对应 §3.6 前三行）。
- **测试**
  1. **T3**：`T_ph≡3` → `G_ph = 3·G_Q(T)`，300/100/4 K 三档相对误差 < 1e-10；
  2. 恒等式 `∫₀^∞ x²eˣ/(eˣ−1)² dx = π²/3`（数值 3.28987，误差 < 1e-8）；
  3. 低温标度：`T_ph≡1` 时 `G(10 K)/G(20 K) = 0.5 ± 1e-3`；
  4. 三条 raise 各一个用例；
  5. `g_cum[-1] == G_ph`（< 1e-12）。
- **DoD**：T3 绿；`∫x²…=π²/3` 进测试；返回值键与 `thermal.json` 契约（§6.4）一致。

#### M1-T5 `dynamical.py`：D 矩阵与虚频

- **文件**：`src/smj/phonon/dynamical.py`、`tests/phonon/test_dynamical.py`
- **签名**
  ```python
  def dynamical_matrix(K, masses, *, scale_to_ev2: bool = True) -> np.ndarray
  def imaginary_modes(D_ev2, *, cm1_tol: float = -10.0) -> list[float]
  def acoustic_sum_residual(K) -> float          # eV/Å²，max|Σ_j K_ij|
  ```
- **实现要点**
  - `D = K/√(M_iM_j) × EV_PER_SQRT²`（`scale_to_ev2=True`），**返回 eV²**，docstring 明写；
  - `imaginary_modes`：`eigvalsh(D)`，取 `<0` 的，`cm⁻¹ = sqrt(-λ)/EV_PER_SQRT×EV_TO_CM1`… 更直接：`cm1 = sqrt(-λ_ev2) * EV_TO_CM1`（因为 `sqrt(λ_ev2)` 就是 ħω 的 eV 值）；
  - `acoustic_sum_residual` 在 **K** 上算（不在 D 上），单位 eV/Å²——这是 §5.1「保留自作用对角块」的理由所在。
- **测试**
  1. 双原子链：解析声学/光学支频率，误差 < 1e-8；
  2. 刚体平移：`K` 满足 ASR 时 `residual < 1e-12`，`D` 有三个零模；
  3. 缩放：`scale_to_ev2=False` 时 `sqrt(eig)` 是 λ 不是 ħω（相对 `EV_PER_SQRT` 倍数）；
  4. 虚频：人为加一个负模 → `imaginary_modes` 返回对应 cm⁻¹，符号正确。
- **DoD**：D 的 eV² 约定在 docstring 与测试里双钉；虚频返回 cm⁻¹。

#### M1-T6 解析基准收口

- **文件**：`tests/phonon/test_negf_analytic.py`
- **内容**：把 T2/T3/T4 组装成一个「端到端人造体系」用例：1D 链（`d`/`v` 解析给出）→ `transmission` → `thermal_conductance`，与闭式解 `G = (k_B²T/h)∫T_ph(x)x²…` 逐点比。
- **DoD**：`pytest tests/phonon -q` 全绿；无网络、无 DFT、< 5 s。

---

### 2.2 M2：几何、力与模板

> 门禁：**T5 / T5b 绿** + provenance 刷新 + demo 不漂。

#### M2-T1 `write_aims_geometry(lattice=)`

- **文件**：`src/smj/junction/geometry.py:68`
- **实现要点**
  ```python
  def write_aims_geometry(path, atoms, constrained=(), lattice=None):
      """lattice=None 时输出与 0.2.2 逐字节相同（W0-1 兼容约束）。"""
  ```
  `lattice` 为 `((a1x,a1y,a1z),(a2x,...),(a3x,...))` 或 `None`；非 None 时在文件**最前**写三行 `lattice_vector`，然后是原子与 `constrain_relaxation`。
- **测试**（`tests/phonon/test_geometry_lattice.py`）
  1. **逐字节回归**：用 `runs/nap-sh/.../aims_scf/geometry.in` 的原子列表，`lattice=None` 输出与现有实现 `==`（不是近似）；
  2. `lattice` 三行格式与 FHI-aims 一致（`%.10f`，与原子行同精度）；
  3. 带 `constrained` 时顺序不变。
- **DoD**：既有 25 条锚点回归不漂；`gen_demo_reference.py --check` 通过；**这是唯一需要刷新 `PROVENANCE.json` 的提交**。

#### M2-T2 两个 control 模板

- **文件**：`src/smj/junction/templates/control_thermal_{relax,force}.in`、`tcontrol.py:102`
- **实现要点**：以 `control_scf.in` 为底，**加** `compute_forces .true.` / `final_forces_cleaned .true.` / `sc_accuracy_forces {sc_accuracy_forces}` / `k_grid {kx} {ky} {kz}`；**去掉** `output aitranss`、`output mulliken`、`restart aims.restart`（§5.4）。relax 版再加 `relax_geometry bfgs {force_threshold}`。
  `make_control(mode=...)` 增加两个 mode，其余签名不变。
- **测试**：渲染出的 `control.in` 里 `compute_forces` 恰出现一次；**不含** `calculate_forces`（拼错会静默无效，附录 A-3）；`k_grid` 行存在；`output aitranss` 不存在。
- **DoD**：两种 mode 的渲染快照测试；既有 mode 输出不变。

#### M2-T3 `forces.read_forces`

- **文件**：`src/smj/phonon/forces.py`、`tests/phonon/test_forces.py`、夹具 `tests/data/aims_force_block_20atoms.out`
- **实现要点**：锚 `Total atomic forces (unitary forces cleaned) [eV/Ang]:`，**取最后一个块**；逐行 `^\s*\|\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)$`；原子数不符 / 非有限值 / 无块 → raise。
- **夹具**：从 `runs/nap-sh/nap-sh/aims_tip_relax/aims.out` 裁出表头 + **末 2 个力块**（20 原子，约 40 行）。**必须进 git**——`runs/` 被 gitignore、`docs/validation-data/` 无任何 `aims.out`。
- **测试**（T5）：20 原子、取第 27 块；逐分量与原文一致（< 1e-15）；截断文件 → raise；把某分量改成 `nan` → raise；原子数改成 19 → raise。
- **DoD**：T5 绿；夹具 < 5 KB。

#### M2-T4 `fc.py`：位移集与 FC 组装

- **文件**：`src/smj/phonon/fc.py`、`tests/phonon/test_fc.py`
- **签名**
  ```python
  def make_displacements(n_atoms, device_indices, amplitude_A) -> list[dict]
  def assemble_fc(disps, forces_by_disp, n_atoms) -> np.ndarray   # (3*N_dev, 3*N)
  def apply_asr(D_ev2, masses) -> tuple[np.ndarray, float]
  def load_lead_fc(path) -> LeadFC                              # 含 sha256 校验
  def assemble_full_fc(k_dev, lead, meta) -> np.ndarray         # (3N,3N) eV/Å²
  ```
- **实现要点**
  - 位移集：`device_indices` × 3 轴 × ±，共 `6·N_dev`；每条记录 `{index, atom, axis, sign, delta_A}`；顺序确定（先原子、再轴、再符号），便于数组下标映射。
  - `K[3i+a, 3j+b] = −(F[+,i,a][3j+b] − F[−,i,a][3j+b]) / (2δ)`；只有 `i ∈ device` 的行有值。
  - `apply_asr`：先算 `acoustic_sum_residual(K)`（**在 K 上**），再对 D 做正交投影 `D ← (I−P)D(I−P)`，`P = T(TᵀT)⁻¹Tᵀ`，T 为三条质量加权平移模。返回 `(D_proj, residual)`。
  - `load_lead_fc`：读 JSON，校验 `in_plane`/`layers`/`fc_cutoff_A`/`sha256`；`sha256` 不符或面内胞与配置不符 → raise（`ph_lead_cell_match` 的数据源）。
  - `assemble_full_fc`：把 `k_dev` 填进器件行，lead 行由 `LeadFC` 的层间块填充；重叠区（器件区与 lead 的交界）用**器件侧实测值**优先，并记录冲突数。
- **测试**：人造 3 原子体系的解析 K；中心差分对 O(δ²) 误差的收敛（δ→δ/2 时误差 /4）；ASR 投影后 `residual < 1e-12`；`sha256` 篡改 → raise；`k_dev` 与 `lead` 重叠冲突 > 0 时 warning 而非 raise。
- **DoD**：`fc.py` 全部纯函数、无 IO 副作用（除 `load_lead_fc`）。

#### M2-T5 `device.build_device`

- **文件**：`src/smj/phonon/device.py`、`tests/phonon/test_device.py`
- **签名**
  ```python
  def build_device(molecule_atoms, cfg: ThermalConfig) -> tuple[Atoms, DeviceMeta]
  ```
- **实现要点**
  1. **先 `strip_thiol_h`**（`geometry.py:160`）——否则直接踩报告陷阱⑥；
  2. 生成 Au(111) 3×3 面内胞（`a=4.07823`，`a/√2=2.885` 最近邻），每侧 `layers` 层，`delta_z=2.35455`；
  3. 分子置于两板之间，`S···S` 沿 z，真空 `vacuum_A`；
  4. 输出 `layer_index`（0=最外，向分子递增）、`principal_layer_id`、`role ∈ {molecule, device_au, lead_au}`、`device_indices`、`constrained_indices`；
  5. `DeviceMeta` 写 `device_metadata.json`（**不改**电子链的 `electrode_metadata.json`）。
- **测试**（**T5b**）：nap-sh 分子 18 原子；S 的邻居恰 1 个 Au、无 H；总原子 90（4 层）/ 126（6 层）；层号单调；`device_indices` 长度 54；`constrained` 与 `device` 互补。
- **DoD**：T5b 绿；`device_metadata.json` 字段与 `thermal_device` 阶段契约一致。

#### M2-T6 `extract_forces.py`（远程 stdlib 脚本）

- **文件**：`src/smj/phonon/extract_forces.py`
- **实现要点**：**只用 stdlib**（`re`/`json`/`sys`），因为集群上不一定装了 `smj`。`python extract_forces.py aims.out forces.json`；输出 `{"n_atoms": N, "forces": [[fx,fy,fz],...], "block_index": k}`。
- **测试**：与 `forces.read_forces` 在夹具上**逐位一致**（同一正则、同一块选择逻辑，用同一夹具交叉验证）。
- **DoD**：两处解析结果一致；脚本 `--help` 可用；无第三方 import。

---

### 2.3 M3：阶段与配置

> 门禁：`smj thermal --prepare` 在 nap-sh 上 dry-run 通过；620 既有用例全绿。

#### M3-T1 `ThermalConfig`

- **文件**：`src/smj/config.py`（三处：`:469` `_KNOWN_SECTIONS`、`:485` 附近 `_THERMAL_KEYS`、`:718` 附近解析）
- **实现要点**：照 `TransportConfig` 的写法，含嵌套 `electrode`/`fc`/`negf` 三段 + `lead_fc` + `temperature_K`（`None` = 继承）。所有未知键走 `_warn_unknown`（不报错，只警告，与既有风格一致）。
- **测试**：完整 `smj.yaml` 解析；缺省值；未知键产生 warning；`temperature_K: null` 时取 `transport.temperature`。
- **DoD**：`RunConfig.thermal` 可用；既有配置解析零变化。

#### M3-T2 `thermal_device` 阶段

- **文件**：`src/smj/stages/thermal_device.py`
- **实现要点**：`depends=("molecule_opt",)`、`needs_submit=False`、`local_only=True`；prepare 读 `molecule_opt` 的 `molecule.xyz` → `build_device` → 写 `thermal/geometry.in`（带 `lattice`）+ `device_metadata.json`。
- **测试**：`smj thermal --prepare` 在 nap-sh 上产出文件；`geometry.in` 含 3 行 `lattice_vector`；原子数 90。
- **DoD**：dry-run 通过；产物可被 `aims` 直接读（用集群上 `aims --dry-run` 或 `preflight` 校验）。

#### M3-T3 `thermal_relax` 阶段

- **文件**：`src/smj/stages/thermal_relax.py`
- **实现要点**：`depends=("thermal_device",)`、`needs_submit=True`、1 个作业；prepare 写 `geometry.in`（含 `constrain_relaxation`）+ `control_thermal_relax.in` + `sub_relax.sh`；collect 用 `aims_out.parse` 取 `max_force_ev_per_A`，缺则 raise（不补 0），并**保留 `aims.restart`** 供 `thermal_fc` 复用。
- **测试**：脚本含 `--mem=12G` 且**不含** `--gres`；collect 对缺 `aims.out` 的 fail-closed；`max|F| > 1e-2` 时阶段报 blocked。
- **DoD**：阶段可 prepare/collect；`aims.restart` 被记录在 artifacts 里。

#### M3-T4 `thermal_fc` 阶段（数组作业）

- **文件**：`src/smj/stages/thermal_fc.py`
- **实现要点**
  - prepare：读 `relaxed_geometry.in` + `device_metadata.json` → `make_displacements` → 写 `disps.json`、`disp_XXXX/{geometry.in,control.in}` × 324、`extract_forces.py`、`sub_fc.sh`（`--array=0-323%8`、`--mem=12G`、**无 gres**、`--time=06:00:00`）；每个 `disp_*` 里**软链或复制 `aims.restart`**（从 `thermal_relax` 拷来，省 60% 的 SCF）。
  - collect：要求 324 个 `forces.json` **齐全**，缺一个 → raise；`assemble_fc` → `apply_asr` → 写 `fc_device.npz` + `fc_report.json`（`asr_residual`、`fc_asymmetry`、`n_disps`、`missing`）。
- **测试**：`make_displacements` 数量 = 6×54；数组脚本字符串含 `--array=0-323%8`；缺 1 个 forces → raise；`aims.restart` 复制逻辑用假文件验证。
- **DoD**：`--prepare` 产出 324 个目录（本地 dry-run）；`jobwait` 零改动可用。

#### M3-T5 `thermal_negf` 阶段

- **文件**：`src/smj/stages/thermal_negf.py`
- **实现要点**：`needs_submit=False`；读 `fc_device.npz` + 冻结 lead FC → `assemble_full_fc` → 取器件块 `d_c`、lead 块 `d_l/d_r`、耦合 `v_l/v_r` → 251 个 ω → `transmission`（`k_weights=None`，**N_k=1**，见 §8 D10）→ `thermal_conductance` → 写 `thermal.json`、`transmission_ph.csv`、`conductance_cum.csv`；`plotting` 出 png/pdf（失败只记 note，与 `junction_analyze` 同风格）；只读 `junction_analyze/results.json` 取 `k_e`/`ZT_e`，取不到写 `null` + warning。
- **测试**：对合成 FC（人造）跑通；`thermal.json` 通过 schema 校验；缺 `results.json` 时 `ZT_total=null` 且不失败。
- **DoD**：纯后处理 < 10 s；`--analyze-only` 可单独重跑。

#### M3-T6 `smj thermal` 子命令

- **文件**：`src/smj/cli.py`
- **实现要点**：照 `cmd_stage`（`:637`）的形状；`--run/--prepare/--submit/--wait/--collect/--analyze-only`；内部用 `THERMAL_STAGES` 按序调用 `StageRunner`。
- **测试**：`smj thermal --help`；`--prepare` 只跑 `thermal_device`；未知 stage 报错清晰。
- **DoD**：argparse 无冲突；既有子命令零变化。

#### M3-T7 注册隔离与回归

- **文件**：`src/smj/stages/__init__.py`、`src/smj/state.py`
- **实现要点**：4 个阶段注册；`STAGE_ORDER` **不动**；`smj status`/`smj plan`/`smj run` 的默认输出**不得出现** thermal 阶段。
- **测试**：`get_stage("thermal_fc")` 可用；`STAGE_ORDER` 长度不变；`smj plan` 输出快照不含 `thermal_`；**620 全绿**。
- **DoD**：T11 的电子链零回归。

---

### 2.4 M4：端到端与标定（要机时）

> 门禁：`thermal.json` 产出 + T6–T9 有记录。**这一批开始产生不可回退的机时成本，先看 §6 的降本方案。**

#### M4-T1 全链路 dry-run

- **内容**：本地 `smj thermal --prepare`（4 阶段）→ 检查 324 个目录、脚本、元数据；用 `preflight` 校验 `geometry.in`。
- **DoD**：无 `UNCONFIGURED` 占位；`--time`/`--mem` 符合 §5.2 三条实测约束。

#### M4-T2 冻结 lead 力常数

- **文件**：`scripts/gen_au_lead_fc.py`、`src/smj/junction/data/au111_3x3_l4.json`
- **内容**：体相 Au、3×3 面内胞、沿 z 周期（≥8 层）、位移中央主层 9 个原子 × 3 轴 × ± → **54 个 SCF**；提取层内与层间 FC 到 10 Å；写 `in_plane`/`layers`/`method`/`fc_cutoff_A`/`sha256`。
- **DoD**：`ph_lead_cell_match` 可通过；文件 < 500 KB；生成脚本可复跑（带 `--check` 比对 sha256）。

#### M4-T3 单位移标定

- **内容**：只跑 `disp_0000`，看 `aims.out` 的 **self-consistency cycles** 与 `| Number of k-points`。
- **判据**：`cycles ≤ 30` → restart 生效，按 §5.5「有 restart」档估（20–36 h）；`cycles > 40` → restart 没生效，改用无 restart 档（85 h）并排查；k 点数回填 §5.5 的 `N_kpts`。
- **DoD**：成本估算从「估算」变成「实测外推」。

#### M4-T4 全量 324 位移

- **内容**：提交数组作业，`wait`，`collect`。
- **DoD**：324 个 `forces.json` 齐全；`asr_residual < 1e-6 eV/Å²`；`fc_asymmetry < 5%`。

#### M4-T5 开发期收敛矩阵（**降本版**，见 §6.2）

- **T6 位移幅度**：只做**子集**（分子 18 原子 + 界面 6 原子 = 24 原子 → 144 位移）在 `0.005/0.02 Å` 下重算，比 `‖K(δ)−K(0.01)‖/‖K(0.01)‖`；若 < 5% 判过，**不跑全量**。若 > 5% 才升级为全量（+1 个完整 FC 集）。
- **T7 η 与网格**：纯后处理，零机时。`η ∈ {0.25,0.5,1.0}`、步长 0.1→0.05，`ΔG/G < 3%` / `< 1%`。
- **T8 电极层数**：**几何级检查**——建 6 层器件、弛豫，比器件区坐标与 4 层版的差；`< 0.01 Å` 判过（不必重算 FC）。若超差再升级。
- **T9 基组**：`tight` 在**同一子集**（144 位移）上重算，报 `ΔG/G`；写进报告，**不作为通过/失败判据**。
- **DoD**：四项各有数字与结论；T6/T8 若升级则单独记录。

#### M4-T6（可选）文献对账 T14

- **前置**：用户明确同意（AGENTS.md：对照臂须先说明理由与代价）。
- **内容**：`octanedithiol` 或 `BDT` 跑一次完整链。
- **判据**：`G_ph(300 K) ∈ [15, 25] pW/K`（**理论对理论**，见方案 §9 的警告）。
- **DoD**：结果单列，不进 nap/az/ind 的生产结论。

---

### 2.5 M5：判据、报告与收口

> 门禁：**T10–T13 全绿**。

#### M5-T1 十条判据进 `CRITERIA`

- **文件**：`src/smj/junction/criteria.py`
- **实现要点**：10 条（方案 §6.5 表）；4 条 `physical` 必须带 `source{doi, locator}`（`Lynn PRB 8, 3493`、`Schwab Nature 404, 974`、报告 §6 陷阱④）；`ph_acoustic_limit` 的 note 里写明「k 平均口径」；`KNOWN_GAPS` 只减不增。
- **测试**：`lint_registry()` 通过；每条判据有正/负用例；`ph_transmission_present` 的 **name 里不含 `te`**（回归用例）。
- **DoD**：`smj review` 对热导 run 能出 `pass/fail/inconclusive`。

#### M5-T2 review 集成

- **文件**：`src/smj/review.py`
- **实现要点**：热导 run 走独立分支；`ZT_total=null` → `inconclusive`（不冒充 `pass`）；有软模但 > −10 cm⁻¹ → `inconclusive` + warning。
- **测试**：三种结论各一个用例。
- **DoD**：方案 §11「`inconclusive` 不得包装成 `passed`」有代码与用例。

#### M5-T3 报告节

- **文件**：`src/smj/report.py`（`:1010` 之后）、`_thermal_summary`
- **实现要点**：`## 5. 热导` 表（`G_ph`/`G_Q`/`N_eff`/`ZT_e`/`ZT_total`/`T_ph(ω>带顶)`/ASR/虚频数/lead FC sha）+ `T_ph(ω)` 与 `G_cum(ω)` 双轴图；**必须显式写明**「热导器件几何（周期板 ~90 原子 3×3）与电子链几何（金字塔团簇 92 原子）不是同一体系，`G_ph` 与 `G/G0` 不可并列」。
- **测试**：报告快照含该警告句；缺 thermal.json 时报告不崩。
- **DoD**：`smj report` 在 nap-sh 的既有 run 上零变化（无 thermal.json）。

#### M5-T4 `thermal.json` schema 校验

- **文件**：`src/smj/analysis.py`（照 `ANALYSIS_SCHEMA` 的写法）
- **实现要点**：`THERMAL_SCHEMA = "smj.thermal.v1"`；必需键校验；`allow_nan=False`。
- **测试**：缺键 → raise；多键 → 保留并 warning。
- **DoD**：契约与方案 §6.4 一致（含新增的 `n_kpoints`）。

#### M5-T5 集成回归

- **文件**：`tests/test_thermal_integration.py`（`@pytest.mark.integration`，默认不跑）
- **内容**：对 `docs/validation-data/` 风格的冻结 `fc_device.npz`（新造一个小的）跑 `thermal_negf`，比对 `G_ph` 到 1e-9。
- **DoD**：默认 `pytest -q` 仍 620 + 新增用例；`-m integration` 时跑通。

---

## 3. 关键路径与并行度

```
M1 (4.00 人日) ──┐
                 ├─→ M3 (5.00) ─→ M4 (2.50 + 机时) ─→ M5 (3.00)
M2 (4.25 人日) ──┘
```

- **M1 与 M2 可完全并行**：M1 只碰新包与纯数值，M2 碰几何/IO/模板。两者文件不相交（M2 的 `fc.py` 依赖 M1 的 `dynamical_matrix` 签名，但可先按签名写）。
- **M3 必须等 M2 的 `build_device` 与模板**，但配置（M3-T1）与 CLI 骨架（M3-T6）可提前。
- **M4 是唯一要排队等机时的**；M5 的判据（T1）与报告（T3）可在 M4 跑数期间并行做。
- **关键路径长度**：M1 → M3 → M4 → M5 = 4.00 + 5.00 + 2.50 + 3.00 = **14.5 人日** + 机时；M2 被 M1 吸收。
- **总人日** 18.75（= 14.5 关键路径 + 4.25 的 M2）。

**建议排法**：第 1 周 M1 + M2 并行（两名开发者；一人串行则 8.25 人日）；第 2 周 M3；第 3 周 M4（机时在后台跑，人同时做 M5）；第 4 周 M5 收口。

---

## 4. 里程碑门禁

| 门 | 条件（全部满足才进下一批） | 证据 |
|---|---|---|
| **M1** | T1–T4 全绿；`pytest tests/phonon -q` < 5 s；无网络/无 DFT | CI 日志 |
| **M2** | T5/T5b 绿；`lattice=None` 逐字节回归绿；`PROVENANCE.json` 已刷新；`gen_demo_reference.py --check` 不漂 | 两个脚本输出 |
| **M3** | `smj thermal --prepare` 在 nap-sh 上 dry-run 通过；`STAGE_ORDER` 不变；620 全绿 | CLI 输出 + pytest |
| **M4** | `thermal.json` 产出；10 条判据全绿；T6–T9 有数字 | 真实 run 目录 |
| **M5** | T10–T13 全绿；`smj review` 出结论；报告含几何差异警告 | 全量回归 |

**任一门不过就不进下一批**——尤其是 M1：解析基准没立住就上 DFT，后面所有的「看起来合理」都无从判断。

---

## 5. 提交与回滚

| 阶段 | 提交粒度 | 回滚方式 |
|---|---|---|
| M1 | 每任务一提交（`units`/`negf`/`landauer`/`dynamical`/`tests`） | 删文件即可，无外部依赖 |
| M2-T1 | **单独一个提交**（唯一内核改动 + provenance 刷新 + demo 对账） | `git revert` 一个提交，基线立刻回来 |
| M2 其余 | 每任务一提交 | 删文件 |
| M3 | 配置/阶段/CLI 各一提交 | 删注册行即完全回退 |
| M4 | 不提交数据（`runs/` 已 gitignore）；lead FC 单独一个提交 | 删 JSON 文件 |
| M5 | 判据/报告/schema 各一提交 | 删判据条目（`KNOWN_GAPS` 只减不增，注意不要反向） |

**整体回滚点**：`thermal.enabled` 默认 `false` + 阶段不进 `STAGE_ORDER` + 新包独立 → 把 `src/smj/phonon/` 与 4 个 stage 文件删掉、`criteria.py` 的 10 条去掉，电子链回到 `608368e` 的行为。

---

## 6. 工时与机时

### 6.1 人日

| 批次 | 任务 | 人日 |
|---|---|---|
| M1 | T1 0.5 / T2 1.0 / T3 0.75 / T4 0.75 / T5 0.5 / T6 0.5 | **4.0** |
| M2 | T1 0.5 / T2 0.5 / T3 0.5 / T4 1.25 / T5 1.25 / T6 0.25 | **4.25** |
| M3 | T1 0.5 / T2 0.75 / T3 0.75 / T4 1.25 / T5 1.0 / T6 0.5 / T7 0.25 | **5.0** |
| M4 | T1 0.25 / T2 0.5 / T3 0.25 / T4 0.5 / T5 0.5 / T6 0.5 | **2.5** |
| M5 | T1 0.75 / T2 0.5 / T3 0.75 / T4 0.5 / T5 0.5 | **3.0** |
| | | **18.75** |

### 6.2 机时（**开发期收敛矩阵是主要成本，不是生产 run**）

方案 §9 的 T6/T8/T9 若各按「完整 FC 集」做，需要 **4 个额外全量集 ≈ 120 h 墙钟**。本计划改用**子集 + 几何级检查**，把开发期机时压到 **1/4**：

| 项 | 内容 | 位移数 | 墙钟（8 路） |
|---|---|---|---|
| lead FC | 体相 54 SCF | 54 | 2–4 h |
| `thermal_relax` | 1 个弛豫作业 | — | 1–2 h |
| 单位移标定 | `disp_0000` | 1 | 0.5–1 h |
| **生产 FC** | 全量 | **324** | **20–36 h** |
| T6 子集 | 24 原子 × 6，δ=0.005/0.02 | 288 | 18–32 h |
| T8 几何级 | 6 层器件弛豫（不跑 FC） | — | 1–2 h |
| T9 子集 | 同上子集，`tight` 基组 | 144 | 9–16 h |
| T7 | 纯后处理 | 0 | 0 |
| **合计** | | | **约 52–93 h（2–4 天）** |

> 若 T6/T9 的子集判据不通过（>5%），才升级为全量：T6 +20–36 h、T9 +20–36 h。**升级必须单独记录理由**，不能默认跑。

### 6.3 成本的三条实测锚点（写进报告）

- 单 SCF 锚点：92 原子非周期 **1563.9 s / 68 cycles ≈ 23 s/cycle**；
- 有 restart 时 15–30 cycles → **0.5–0.9 h/SCF**（§5.5）；
- 内存峰值 **871.9 MB/task** → `--mem=12G` 有 13× 余量。

---

## 7. 风险处置映射

| 风险 | 触发点 | 处置 | 落到哪个任务 |
|---|---|---|---|
| R1 SCF 比外推慢 | M4-T3 标定 | 退到 §5.5 最小档（120 位移） | M4-T3 |
| R2 力噪声大 | M4-T4 `fc_asymmetry > 5%` | 收紧 `sc_accuracy_forces` 到 1e-4；或 δ→0.02 Å | M4-T4 |
| R3 器件区虚频 | M5-T1 判据 | 报告列出；`<−10 cm⁻¹` 才 fail | M5-T1/T2 |
| R4 `G_ph` 差 3× | M4-T4 后 | 按方案 §3.5 逐条查（带顶归零→ASR→层数→接触→基组），**不调参** | M4-T5 |
| R5 队列配额 | M4-T4 | `%8` → `%4` | M3-T4 |
| R6 整条线失败 | 任意 | 删新包 + 注册行 | §5 |
| R7 无先例可对账 | 全程 | 解析基准 + 内部一致性 + （可选）T14 + 报告标注「首次实现」 | M1-T6 / M4-T6 / M5-T3 |

---

## 8. 待决项（**开工前需要拍板**）

| # | 问题 | 建议 | 影响 |
|---|---|---|---|
| **D9** | 要不要跑 T14 文献对账（C8 或 BDT）？ | **建议跑一次**，但须用户明确同意 | +1 个完整 FC 集（20–36 h）；结果单列 |
| **D10** | NEGF 的横向 k：`N_k=1`（3×3 超胞的 Γ）还是真 k 分辨？ | **`N_k=1`**。3×3 超胞的 Γ 已经把 9 个原胞 q 点折叠进来，是标准做法；真 k 分辨需要跨胞 FC，而有限位移只给胞内力常数。API 保留 `k_weights`，`N_k=1` 时归一化自动为 1 | 若坚持 k 分辨，M3-T5 增加 ~2 人日且需要 3×3×3 超胞的 DFT |
| **D11** | lead 主层厚度：1 层（27 DOF，只保留 V₁）还是 4 层（108 DOF，耦合只到相邻主层）？ | **先 1 层**，M4-T5 加一条 T8b 量 `V₂` 的影响；若 `ΔG/G > 5%` 改 4 层（NEGF 代码不变，只是矩阵变大） | 影响 `G_ph` 的 5–10% 量级 |
| **D12** | 热导阶段要不要出现在 `smj status`/`smj plan` 的默认输出里？ | **不出现**，只在 `smj thermal` 下可见 | M3-T7 的回归用例 |

> D10 与 D11 是**本计划对方案的两处技术收紧**，都写在任务里；若用户不同意，按备选走，人日与机时相应增加。

---

## 9. 验收映射（T → 任务）

| T | 由哪个任务交付 | 批次 |
|---|---|---|
| T1 | M1-T1 | M1 |
| T2 | M1-T3 | M1 |
| T3 | M1-T4 | M1 |
| T4 | M1-T2 | M1 |
| T5 | M2-T3 | M2 |
| T5b | M2-T5 | M2 |
| T6 | M4-T5（子集） | M4 |
| T7 | M4-T5（纯后处理） | M4 |
| T8 | M4-T5（几何级） | M4 |
| T9 | M4-T5（子集） | M4 |
| T10 | M4-T4 + M5-T1/T2 | M4/M5 |
| T11 | M3-T7 | M3 |
| T12 | M2-T1 | M2 |
| T13 | M2-T1 | M2 |
| T14 | M4-T6（可选，须同意） | M4 |

---

## 10. 与项目规则的关系

- **不默认加跑**：M4 的 T6/T8/T9 是**开发期一次性验收**，不是每次 run 的控制臂；T14 是**对照臂**，须用户明确同意、结果单列。`thermal.enabled` 默认 `false`，`smj run` 完全不受影响。
- **生产默认口径不变**：S 质子化（`strip_thiol_h`）、吸附位、`pin_contact`/`tip_direction`/`anchor_au_s_target` 全部沿用电子链，热导支线不引入新的结构自由度（除几何形状本身，见方案 §4）。
- **判定是机械的**：`ph_imaginary_modes` 有软模、`ZT_total=null` → `inconclusive`，不得包装成 `passed`。
- **协议不一致必须标注**：热导器件（周期板 ~90 原子 3×3）与电子链器件（金字塔团簇 92 原子）**不是同一个体系**，M5-T3 的 report 必须显式写这句话。
- **内核改动**：只有 `geometry.py:68` 一处，M2-T1 单独提交，刷新 `PROVENANCE.json` + `gen_demo_reference.py --check`；`scripts/crosscheck_agenticdft.py` 不进 CI。
