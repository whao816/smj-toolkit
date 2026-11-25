# 单分子结热导（G_ph）功能开发方案（2026-11-23）

> 基线：`master @ 608368e`，`python -m pytest -q` → **620 passed, 2 deselected**（22.7 s，Python 3.12.10）。
> 上游依据：`D:\WorkBuddy\单分子结计算\docs\单分子结热导率计算方法调研.md`（369 行；2026-11-23 全文核对，**勘误 13 条见附录 A**）。
> 本方案是**开发方案**，不含任何热导计算结果。文中所有「实测」数字都来自本仓真实 run 或一手文献；凡未核实的一律标 **UNVERIFIED**。
> 软件与机器现状（gpu1 逐项实测，2026-11-23）见 **§12**。
> 已发表工作的技术栈/软件栈调研见 `docs/thermal-conductance-literature-2026-11-23.md`（三条主线 + 四个反直觉核实结果；本方案 §4/§6 的选型依据）。

---

## 0. 一句话结论

| 现状 | 缺口 | 本方案 |
|---|---|---|
| 六阶段链只算电子：`TE.dat` → `results.json`（`log10_G_over_G0_like`/`S`/`k_e`/`ZT_electronic`）。全仓 `phonon`/`hessian`/`force_constant` 零命中；`aims_out.py` 只解析一个「Maximum force component」标量，**没有力矢量解析器**；SCF 不开力关键字（实测 `aims_scf/aims.out` 无 `Total atomic forces` 段）。 | ① 没有力常数（FC）生产链路；② 没有声子 NEGF / Landauer 积分；③ **现有金字塔团簇几何在物理上做不了声子 NEGF**（无横向周期性，无法构造主层与半无限 lead）；④ 配置/几何/判据/报告/CLI 全链路都没有热导位置。 | 新增**独立的「热导支线」**：`thermal_device`（周期板器件几何）→ `thermal_relax` → `thermal_fc`（6N 位移 + 力，数组作业）→ `thermal_negf`（自写声子 NEGF + Landauer）。**不进 `STAGE_ORDER`，默认关闭**，`smj run` 行为逐字节不变。产出 `thermal.json`（`G_ph`、`T_ph(ω)`、`N_eff`、`ZT_total`）。 |

**最大的科学增值不是「补一个 κ_ph」，而是把 `ZT` 从电子口径改成总口径。** 本仓已有 `k_e_pW_per_K` 与 `ZT_electronic`，`ZT_total = ZT_e / (1 + G_ph/k_e)` 是恒等式。用实测 `k_e` 代入可见：强耦合结（nap-sh 116 pW/K、az 224 pW/K）热导只是 10–20% 的修正；**弱耦合结（2CN 1.24 pW/K、hdt 0.91 pW/K）声子项压倒电子项，`ZT` 掉一个数量级**。这把「萘/薁/茚 对比」从电子输运问题变成了**热-电联合优化**问题（§2.3）。

---

## 1. 现状与缺口

### 1.1 已经有的（可直接复用）

| 资产 | 位置 | 热导支线怎么用 |
|---|---|---|
| 阶段框架（prepare/collect 分离、fail-closed、`UNCONFIGURED` 拒提） | `src/smj/stages/base.py`、`src/smj/runner.py` | 新阶段照 `junction_analyze.py` 骨架写 |
| 一个阶段 = 一个 sbatch 作业 + 哨兵等待 | `src/smj/junction/jobwait.py:49 build_probe_command` | 数组作业的 `squeue -j <id>` 天然列出全部 element，`wait_many` **不用改** |
| FHI-aims 输入渲染与模板 | `src/smj/junction/tcontrol.py:102 make_control`、`templates/control_scf.in` | 新增 `control_thermal.in`（加 `compute_forces`/`final_forces_cleaned`/`k_grid`） |
| 几何读写、`constrain_relaxation` 冻结 | `src/smj/junction/geometry.py:68 write_aims_geometry`、`73-84` | 加可选 `lattice=` 参数（见 §6.1） |
| 电极元数据（tip/linker/surface/backplane 索引） | `electrode_build/electrode_metadata.json` | 热导支线**不复用**金字塔团簇，但复用「索引化元数据」这个模式 |
| 唯一公式源 + 非有限值 fail-closed 的写法 | `src/smj/junction/landauer.py:174 transport_moments` | `G_ph` 的公式也集中在一个函数里 |
| 判据注册表（physical 必须带 `source`，engineering 免源） | `src/smj/junction/criteria.py` | 新增 10 条热导判据（§6.5） |
| 报告 / review / results 契约 | `report.py:1010-1043`、`review.py`、`analysis.py:33 ANALYSIS_SCHEMA` | 加一节 + 加 schema |
| 真实 run 回归（`@pytest.mark.integration`） | `tests/test_demo_regression.py` | 热导侧加一个同构的冻结参考 |

### 1.2 缺口

| # | 缺口 | 证据 |
|---|---|---|
| G1 | 没有力矢量解析 | `src/smj/aims_out.py:41` 只有 `_FORCE_RE`（`Maximum force component`）；`calculate_forces`/`compute_forces` 全仓零命中 |
| G2 | SCF 不产出力 | 实测 `runs/nap-sh/nap-sh/aims_scf/aims.out` 第 23484 行直接 `Have a nice day.`，**无 `Total atomic forces` 段** |
| G3 | 没有 FC 生产链路 | 全仓无位移生成、无力常数矩阵、无 ASR 投影 |
| G4 | 没有声子 NEGF | 无动力学矩阵、无表面格林函数/decimation、无声子 Landauer 积分 |
| G5 | **几何不可用** | 全部 60+ 真实 run 的 `geometry.in` 中 `lattice_vector`/`k_grid` 计数为 **0**；`write_aims_geometry` 不支持晶格矢量；`config.py` 无 `k_grid` |
| G6 | 电极模型是唯象的 | `src/smj/junction/self_energy.py` 给每个电极原子一个**常数虚部** `DEFAULT_SIGMA_NEAR=0.05`/`DEFAULT_SIGMA_FAR=0.10` Ha，docstring 明写「是**模型常数**, 不是收敛参数」——这不是半无限 lead，套到声子上会让 `G_ph` 变成 σ 的函数 |
| G7 | 电极元数据无层号 | `electrode_metadata.json` 无逐原子 `layer_index`；电极随分子倾斜，按 z 聚类得不到层（实测 z 连续分布） |
| G8 | 配置/报告/CLI 无位置 | `config.py:469 _KNOWN_SECTIONS`、`report.py:1010` 输运表、`cli.py:1319` 起 parser |

### 1.3 本次核查的实测证据

- **体系规模与电子基线**（`runs/<n>/<n>/aims_scf/geometry.in`、`junction_analyze/results.json`）：nap-sh 92 原子 / Au 74 / `G/G0=0.1903` / `S=42.08 µV/K`；az 92 / 0.3806 / 54.07；ind-sh 91 / 0.3170 / 39.83；bdt 86 / 0.2789 / 34.49；2cn 128 / 0.00132 / −33.83；hdt 94 / 0.0016 / 13.61。
- **电子热导（关键）**：nap-sh `k_e=116.2 pW/K`、`ZT_e=0.0674`；az `224.0`/`0.1155`；bdt `165.6`/`0.0466`；**2cn `1.24`/`0.0284`；hdt `0.91`/`0.0075`**。与 Wiedemann–Franz 自洽：`L0·T·G = 2.44e-8 × 300 × (0.1903 × 7.748e-5) = 108 pW/K`（代码 116，差在 `L2−L1²/L0` 项）。
- **SCF 成本锚点**（`docs/validation-0.2.0.md` §4，8 核 + 1 GPU）：92 原子 **1563.9 s**；128 原子 **2896.4 s**；TE 阶段与体系大小几乎无关（218–251 s）。
- **力块格式**（`runs/nap-sh/nap-sh/aims_tip_relax/aims.out`，20 原子，共 **27 个力块**，首块第 8305 行、末块第 80992 行）：
  ```
    Total atomic forces (unitary forces cleaned) [eV/Ang]:
    |    1         -0.262905381321659E+00          0.954633163814336E+00        0.560751120457043E+00
  ```
  正则 `^\s*\|\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)$`；relax 每步一块，**解析必须取最后一块**。
- **FHI-aims 版本**：集群是 **221103**（`tcontrol.py:58` 注释、README「221103 不认 override_initial_charge_check」）。DFPT 自 **240507** 起文档化（[arXiv:2501.16091](https://arxiv.org/abs/2501.16091)），**当前装不上**（§8 D3）。

### 1.4 为什么这不是「再加一个后处理」

电子链路的 `junction_analyze` 是纯后处理：读 `TE.dat`、积分、写 JSON，零新作业。热导不是——它需要**新的物理对象（力常数）**，而力常数只能由**新的 DFT 作业**产出，且这些作业必须跑在一个**新的几何**上。工作量分布：FC 生产占 70–80%（与调研报告 §8 的判断一致），NEGF 数值只占小头。方案必须按这个比例分配。

---

## 2. 目标与边界

### 2.1 目标（本方案完成时）

1. `smj thermal --run` 能在**一个已验证体系**（nap-sh）上从分子几何一路跑到 `thermal.json`，产出 `G_ph(T)`、`T_ph(ω)`、`G_cum(ω)`、`N_eff`、`ZT_total`。
2. NEGF 数值有**解析基准**：一维单原子链的 `T_ph` 与解析色散逐点吻合（§9 T2）；`T_ph≡3` 时 `G_ph` 与 `3G_Q` 吻合到 1e-10。
3. 单位换算只有**一处常量**，且有跨软件回归（phonopy `VaspToTHz=15.633302`）。
4. 十条热导判据进 `CRITERIA`，`smj review` 能对热导 run 出 `pass`/`fail`/`inconclusive`。
5. `smj report` 多一节「热导」，含 `T_ph(ω)` 图与 `ZT_total`。
6. **电子链路零变化**：620 个既有用例全绿，`gen_demo_reference.py --check` 数值不漂，`smj run` 默认链不变。

### 2.2 不做（避免越界）

- **不做非谐**（NEGF-PT / SCLMD / RNEMD）。调研报告 §2 把非谐列为 M3/M5/M6，结论是简谐 NEGF 在 T>100 K 高估 ~15%（Li et al., PRB 104, 245413 (2021)，[arXiv:2202.06490](https://arxiv.org/abs/2202.06490)，ReaxFF/烷链体系）。本方案在报告里**显式声明 +15% 系统偏差**，不实现。
- **不做力场 / MLP 快速档**（报告 §7 的 A 档、C 档）。理由：与本仓「FHI-aims 单一引擎 + 冻结基线」的定位冲突，且引入第二套力常数来源会让 `G_ph` 失去可比性。**文献侧证据**：C 支（LAMMPS + 经验势 + NEMD）是与 A 支并列的独立路线，但其 Au–S 参数是文献传承的 3 参数 Morse（`D_e=0.38 eV`、`α=1.47 Å⁻¹`、`r₀=2.65 Å`，引自 Ong, *JPCC* **118**, 7288 (2014)），**没有针对具体分子结的 DFT 拟合**；同一路线在 C60 二聚体上与 DFT-NEGF 差 **6 倍**（41.5 vs 7.0 pW/K，Li/Wang/Segal, [arXiv:2402.06172](https://arxiv.org/abs/2402.06172)）。要做力场档，前提是先有本方案的 DFT 基准作对账目标——顺序不能反。
- **不做电子-声子耦合、不做热电子的 Peltier 项**。
- **不做热整流（方向性）**。简谐 Landauer 框架下 `T_ph(ω)` 对 L↔R 严格对称，线性响应下两端热导必然相等（Onsager 互易）——**薁的 A→B / B→A 方向性差异这个方法在原理上给不出来**。见附录 A-9。若将来要做整流，需要非谐 NEGF 或带温差偏置的 NEMD，是另一个立项。
- **不接第三方声子输运包**（GOLLUM / PHtrans / QuantumATK）。**文献核实后理由比原先更硬**：GOLLUM 2014 原文**没有声子模块**（结论段写「including the incorporation of phonon transport」为 future work）；PHtrans 的实际输入是 **GULP** 而非 SIESTA（Papior 本人说明）；QuantumATK 模块存在但**检索范围内没有一篇论文点名用它算单分子结 `κ_ph`**。详见 §4.2、§12.2 与附录 A-10～A-12。
- **不默认加跑收敛扫描**（电极层数、位移幅度、η、基组）。这些是**开发期一次性的验收项**（§9 T6–T9），不是每次 run 的控制臂；符合 AGENTS.md「不要默认加跑对照」。
- **不自动改用户几何**：S 的质子化状态、吸附位一律沿用电子链已确认的口径（AGENTS.md 生产默认 `pin_contact: false` / `tip_direction: axis` / `anchor_au_s_target: null`）。

### 2.3 科学增值：`ZT_total`

本仓 `ZT_electronic = 1/((L0·L2)/(L1·L1) − 1)`（`landauer.py:194`）只含电子热导。加进 `G_ph` 后：

```
ZT_total = ZT_e / (1 + G_ph / k_e)          # k_e、G_ph 同为 W/K（或同为 pW/K）
```

这是恒等式，不是近似。以 `G_ph = 20 pW/K`（附录 B 的实验量级）做的敏感度估算——**仅为量级演示，不是计算结果**：

| run | `k_e` (pW/K) | `ZT_e` | `ZT_total` | 声子修正 |
|---|---|---|---|---|
| az | 224.0 | 0.1155 | 0.1060 | −8% |
| bdt | 165.6 | 0.0466 | 0.0415 | −11% |
| nap-sh | 116.2 | 0.0674 | 0.0575 | −15% |
| **2cn** | **1.24** | 0.0284 | **0.0017** | **−94%** |
| **hdt** | **0.91** | 0.0075 | **0.0003** | **−96%** |

**结论性判断（可证伪）**：`G_ph` 由质量-弹簧网络与 Au–S 接触刚度决定，与 HOMO–LUMO 能隙无关；因此 2cn 的 `G/G0` 比 az 小 300 倍，`G_ph` 却应与之**同量级**。于是「声子底噪 ~20 pW/K」在弱耦合结上成为 `ZT` 的天花板。这直接接上调研报告 §7 的三个预期科学产出（π 共轭 vs 非共轭、锚定基团质量效应、以及非对称结的接触刚度差），并把它们统一到「`ZT_total` 随电子耦合强度的标度」这一个可检验命题上。

> ⚠️ 报告 §7 的第 ② 项「薁非对称性带来 10–30% **方向性**差异」**这个方法给不出来**：简谐 Landauer 的 `T_ph(ω)` 在 L↔R 交换下严格不变，线性响应热导必然对称（Onsager）。非对称只影响**绝对值**（两端 Au–S 接触刚度不同 → `Σ_L ≠ Σ_R`），不产生整流。见附录 A-9。

---

## 3. 物理与数值规格

### 3.1 量：`G_ph`（W/K），不是 κ

单分子结没有确定截面积，**禁止**输出 `W/(m·K)`。输出量固定为 `G_ph_pW_per_K`（端到端热导）。同时输出 `N_eff = G_ph / G_Q`（有效通道数），因为这是跨体系可比的量。

### 3.2 单位与唯一换算常量

内部单位：**eV、Å、amu**；能量轴用 `ħω`（单位 eV），角频率轴用 `ω`（单位 rad/s）仅在积分核里出现。

**换算常量（本方案唯一一处，`src/smj/phonon/units.py`）**：

```python
# λ 单位 eV/(Å²·amu) → ħω 单位 eV
EV_PER_SQRT = 0.064654     # = ħ·sqrt(1 eV/(Å²·amu))，CODATA 2018
```

推导（可复核）：`1 eV/(Å²·amu) = 1.602176634e-19 / (1e-20 × 1.66053907e-27) = 9.6485e27 s⁻²`，`√λ = 9.8227e13 rad/s`，`ħω = 1.054571817e-34 × 9.8227e13 = 1.03587e-20 J = 0.064654 eV`。交叉核对：`ν = ω/2π = 15.6333·√λ THz`，正是 phonopy 的 `VaspToTHz = 15.633302`（[官方 `units.py`](https://github.com/phonopy/phonopy/blob/develop/phonopy/units.py)）；`hν = 64.654·√λ meV`。

> ⚠️ **调研报告 §3.4 与 §9 的 `0.2543` 是错的**，偏大 3.933 倍（`0.2543² = 0.064668`，即把本该开一次方的系数又开了一次方）。照抄会把 Au 声子带上移到 57.6 meV、C–C 伸缩移到 3200–4200 cm⁻¹（X–H 区）。详见附录 A-1。

### 3.3 Landauer 积分

用与 `ħω` 直接配对的等价形式（与报告 §1 的 `ħ²/(2πk_BT²)` 形式逐项等价，推导见附录 A-2）：

```
G_ph(T) = (k_B²T/h) ∫₀^∞ T_ph(x) · x² eˣ/(eˣ−1)² dx,     x = ħω/(k_BT)
```

- `T_ph ≡ 3`（三条声学支全透射）时给出 `3·π²k_B²T/(3h) = 3G_Q`，即 `∫₀^∞x²eˣ/(eˣ−1)²dx = π²/3`。**这条恒等式进单元测试。**
- 积分上限 `ħω_max` 取 `25 meV`（Au 声子带顶 19.5 meV + 余量，见附录 B-1）。超过带顶 `T_ph` 必须恒为 0；把它设成分子最高频率（100–400 meV）是报告 §6 的第 ⑨ 号陷阱。
- 低温行为：核函数只在 `x ~ O(1)` 采样，故 `G ∝ T`；`G_Q(300 K)=283.9`、`G_Q(100 K)=94.64`、`G_Q(4 K)=3.786 pW/K`（[Schwab et al., Nature 404, 974 (2000)](https://doi.org/10.1038/35010065)）。

### 3.4 动力学矩阵与力常数

```
K_ij       = −∂F_j/∂u_i                # eV/Å²，3×3 块，实对称
D_ij       = K_ij / √(M_i M_j)         # eV/(Å²·amu)，对称
D          → ASR 投影 → 声子谱
Σ_L/R      = V_c,α g_α(ω) V_α,c        # g_α 用 López-Sancho/decimation
G^r(ω)     = [(ħω)²I − D_c − Σ_L − Σ_R]⁻¹      # 注意：D 已换算到 eV² 量纲
Γ_α        = i(Σ_α^r − Σ_α^a)
T_ph(ω)    = (1/N_k) Σ_k w_k · Tr[G^r_k Γ_L,k G^a_k Γ_R,k]
```

> ⚠️ **`1/N_k` 这个归一化必须写死在函数里**（`negf.transmission` 的 docstring 与单元测试都要钉）。每个横向 k 点代表 BZ 的 `w_k` 份额，`Σ_k w_k = 1`。若改成「对所有 k 求和」而不除 `N_k`，`T_ph(ω→0)` 会变成 `3·N_k` 而不是 `3`，`G_ph` 直接放大 `N_k` 倍（本配置 `3×3×1` 是 3–9 倍）——**这正是 §3.6 那两条 raise 想抓的错，却抓不到自己**。判据 `ph_acoustic_limit`、`ph_quantum_bound` 都只在「k 平均」口径下成立。

实现要点（全部写进 `src/smj/phonon/negf.py` 的 docstring）：

- **`(ħω)²` 而不是 `ω²`**：一旦 `ω` 用 eV 表示，`D` 必须乘 `EV_PER_SQRT²`，否则量纲错。这是最容易错的一行，单元测试直接钉死（§9 T1）。
- **ASR**：对 `D` 做投影 `D ← (I − P)D(I − P)`，`P` 由三个平动零模张成；同时记录残差 `‖D·1‖`。
- **对称化**：`K ← (K + Kᵀ)/2` 前先记录 `‖K − Kᵀ‖/‖K‖`（中心差分下应为 O(δ²) 小量，见 §5.1）。
- **η**：`η = 0.5 meV` 默认；周期性 lead 使 `T_ph(ω)` 连续，η 只需抹掉数值锯齿。**文献锚点**：A 支用 `η = 1e-5 a.u. ≈ 0.27 meV`（Klöckner et al., [PRB **95**, 245404 (2017)](https://arxiv.org/abs/1704.00591) 是唯一给全 NEGF 参数的一篇），本方案 0.5 meV 与其同量级。开发期验证 `η ∈ {0.25, 0.5, 1.0} meV` 对 `G_ph` 影响 < 3%（§9 T7）。**注意（D14，2026-11-25 实测）**：η 不是小量修正而是**主敏感量**——同一份力常数下 η 由 0.5 降到 0.05 meV，`G_ph` 从 10.9 升到 13.7 pW/K（+25%），单调外推到 η→0 约 13.8。文献锚点 0.27 meV 因此不只是一个数值选择，它决定报出来的数落在哪一档。
- **网格**：`ħω ∈ [0, 25] meV`，步长 `0.1 meV` → 251 点；梯形积分。收敛判据：步长减半 `ΔG/G < 5%`（**D14**，见 §8）。

### 3.5 电极与自能

**必须是真半无限 lead**，不能沿用 `self_energy.py` 的常数虚部模型（G6）。理由不是洁癖：常数 σ 会把 `G_ph` 变成 `σ_near/σ_far` 的函数，而这两个数在声子问题里没有任何标定依据（电子的 0.05/0.10 Ha 来自冻结的 BDT/C4/2CN 基准，与声子无关）。

主层（principal layer）取 **一层 Au(111)（3×3 = 9 原子）**，`D_L`/`V_L` 为 27×27，decimation 每次迭代 `O(27³)`，251 个频点 < 1 s。lead 的力常数用**冻结数据文件**（§5.4），不每次现算。

**这套结构与文献 A 支完全同构**：Pauly 组的电极力常数也是**另算一次体相**（原文「a big cluster of several hundred atoms」）再做 decimation，器件区与 lead 的化学环境不一致是他们**显式接受的近似**（[Cui et al., Nature **572**, 628 (2019)](https://doi.org/10.1038/s41586-019-1420-z)；[Klöckner et al., PRB **95**, 245404 (2017)](https://arxiv.org/abs/1704.00591)）。我们唯一的差别是：lead 的 FC 来自**二维周期板**而不是团簇，因此带横向 k 采样；团簇电极做不到这一点，只能给分立谱。

### 3.6 数值纪律（fail-closed，与本仓现有风格一致）

| 检查 | 阈值 | 动作 |
|---|---|---|
| `T_ph(ħω > 19.5 meV)` | > 1e-6 | **raise**（lead 里没有模式，非零即错） |
| `T_ph(ω→0)` | > 3 | **raise**（超过声学支数） |
| `G_ph` | > `3G_Q` | **raise**（物理不可能） |
| 器件区虚频 | < −10 cm⁻¹ | **raise**；`−10…0` 记 warning（金属界面软模噪声） |
| `‖D·1‖`（ASR 残差） | 相对残差 `max\|Σ_j K_ij\| / max\|K\|` > 1e-3 | **warn**，超过则做 ASR 投影并把投影前后残差都记下来（D13；绝对门限 1e-6 eV/Å² 被实测噪声否掉） |
| `‖K − Kᵀ‖/‖K‖` | > 5% | **raise**（中心差分失效或力噪声过大） |
| 力收敛 | `max|F|` 缺失或 `sc_accuracy_forces` 未达 | **raise**（不猜、不补 0） |

---

## 4. 器件几何：本方案的核心决策

### 4.1 为什么现有金字塔团簇不能直接做声子 NEGF

三个独立理由，任何一个都足以否决：

1. **主层不存在**。`Electrode111` 的层原子数是三角数 `1, 3, 6, 10, 15, …`（`electrode.py`），**层与层不等价**。decimation 要求 lead 是同一主层的无限重复；金字塔既无横向周期性，也无重复单元，`Σ_L/R` 无法定义。
2. **谱是离散的，而且正好缺掉最要紧的那一段**。有限团簇的振动谱是一堆分立峰，没有 0–19.5 meV 的连续声子带。而全部热流都走这一段。用 η 把分立峰抹成连续谱，`G_ph` 就变成 η 的函数——**这是把数值伪影当物理**。
3. **唯一的替代方案是唯象的**。若沿用 `self_energy.py` 的常数虚部吸收边界，`G_ph` 由 `σ_near=0.05`/`σ_far=0.10` 决定，而这两个数是为电子链的冻结基线标定的，与声子无任何关系。得到的数无法写进论文。

> 这不是「本仓做得不够」，而是**电子链路和声子链路对电极的要求本来就不一样**：aitranss 用常数 σ 是有意的模型选择（并在 docstring 里声明为模型常数），声子 NEGF 没有对应的合法简化。

### 4.2 三个候选

| | A. 周期板 + 自写 NEGF（**推荐**） | B. 现成包（GOLLUM / PHtrans / QuantumATK） | C. 团簇 + 唯象 σ（否决） |
|---|---|---|---|
| 几何 | Au(111) 3×3、每侧 4 层、真空 ≥15 Å、`k_grid 3 3 1` | 同 A | 复用现有 92 原子团簇 |
| lead 模型 | 真半无限（decimation） | 真半无限 | 常数虚部 |
| 力常数来源 | 本仓 FHI-aims 有限位移 | **必须换电子结构引擎**：GOLLUM 的声子模式读 SIESTA `.FC`（[Gollum2 manual](http://www.physics.lancs.ac.uk/gollum/images/manual/Gollum2_manual.pdf)，**且 gpu1 上没装**）；PHtrans 不吃原始 `.FC`，要 TSHS 类容器（[TBtrans/PHtrans 手册](https://siesta-project.gitlab.io/-/documentation/siesta-docs/-/jobs/14622928758/artifacts/public/reference/tbtrans.html)，SIESTA 5.4.2 与 phtrans 二进制**已装在 gpu1**，见 §12）；QuantumATK 商业 | 本仓 FHI-aims |
| 与现有链的关系 | 加新几何，**电子链零改动** | 需要引入第二套 DFT 引擎与格式转换 | 复用，但结果不可发表 |
| **文献先例** | ✅ **本方案就是 A 支本身**：Pauly 组（TURBOMOLE + DFPT + 自研 NEGF + decimation）2016–2019 连续四篇；Sadeghi–Lambert 组（SIESTA + 0.01 Å 中心差分 + 自研 NEGF）2015–2021。**两组的输运都是自写码** | ❌ **检索范围内没有一篇用它算过单分子结 `κ_ph`**：GOLLUM 2014 原文无声子模块；PHtrans 实际输入是 GULP；QuantumATK 模块零文献（附录 A-10～A-12） | ❌ 无 |
| 成熟度 | 自写，靠解析基准兜底；**无 FHI-aims 先例**（§10 R7） | PHtrans 官方原话「**not as tested as TBtrans**」；GOLLUM 的 Au–分子–Au 声子示例未找到（**UNVERIFIED**） | — |
| 结论 | ✅ | ❌（除非已有 QuantumATK 许可） | ❌ |

**候选 B 的致命点**：三条路都要求**放弃 FHI-aims 产力常数**。SIESTA 需要重做整个电子结构链路（基组、赝势、`MD.TypeOfRun fc`），等于把已验证的 BDT/AZ/2CN 基线推倒重来；QuantumATK 是另一套闭源引擎。本仓刚在 2026-11-20 把 junction 内核收回自维护（A 方案），再引入第二套引擎与方向相反。

**文献核实把 B 的每一条退路都堵掉了**（2026-11-23，详见 [thermal-conductance-literature-2026-11-23.md](thermal-conductance-literature-2026-11-23.md)）：

1. **GOLLUM 的声子能力被高估**。[Ferrer et al., NJP **16**, 093029 (2014)](https://arxiv.org/abs/1502.04966) 结论段原文：「Further developments are in the pipeline, **including the incorporation of phonon transport**」——标题里的 "thermal" 指的是**电子热导**（Wiedemann–Franz）。声子模式是后续版本才加的（手册 mode 5），但**没有一篇分子结 `κ_ph` 论文点名用它**，包括 Mosso 2019（他们只写「following the method described in [2,3]」）。报告 §8 把它列为「直接出热导」的那条证据链是断的。
2. **PHtrans 的输入不是 SIESTA**。Papior 在 siesta-l（2019-03-31）原话：用 SIESTA 做声子透射「**extremely costly** since it requires a huge system to ensure the phonons are 'relaxed' in the electrode region. **This is why we (in the paper) did phonon transmissions using GULP**」；sisl 的官方例子就是 GULP→PHtrans。他 2019 年还在说从 DFT 提取动力学矩阵「**still work to do**」。
3. **QuantumATK 有模块、零文献**。`PhononTransmissionSpectrum` + `thermalConductance()` 确实存在且文档完整，但官方唯一教程是 CNT 同位素热电；多轮检索（web / arXiv / OpenAlex / Crossref / DBLP）**未出现任何一篇点名该模块算单分子结 `κ_ph`**。反例：Markussen（后来去了 QuantumWise）做分子结声子时用的是 **GPAW + 自写码**。
4. **VASP + phonopy 在这条线上同样缺席**。能确证的「DFT 出力常数 + 自写输运」实例，DFT 侧是 TURBOMOLE / SIESTA / GPAW / QE / Gaussian03，**没有 VASP**。原因是 phonopy 的超胞 FC 只解决中心区，分子结还需要**半无限电极的表面格林函数与力常数拼接**——这条链在 LCAO / 量子化学码上现成，在平面波框架里要另写一套。

**结论**：候选 B 不是「少写点代码」，而是「换一套生态，并且换完之后仍然没有先例可对账」。**A 支是文献里唯一被反复走通的路**，而且不需要任何现成包。

**候选 A 的新增成本**：几何层要支持 `lattice_vector` + `k_grid`（G5）。这是本方案唯一需要触碰 `src/smj/junction/` 的地方，改动面 ≤ 15 行且向后兼容（§6.1）。

### 4.3 推荐配置

```
面内：3×3 Au(111) = 9 原子/层（报告要求 ≥3×3）
层数：每侧 4 层（报告 4–6；4 层是成本/收敛折中，§9 T8 验证）
真空：≥15 Å（FHI-aims 对 ≥50 Å 真空几乎不计成本，可直接给 20 Å）
器件区（参与位移）：分子 18 + 每侧 2 层 36 = 54 原子
lead（半无限）：第 3 层起，力常数取冻结文件
k 点：3×3×1
总原子数：≈ 90（72 Au + 18 原子分子；nap-sh 实测分子为 92−74=18 原子）
```

分子必须**在板内重新弛豫**（`thermal_relax`），不能直接搬电子链的团簇几何——S···S 距离与 Au–S 接触角都不同。

---

## 5. 力常数生产

### 5.1 位移集

- **中心差分**：每个器件区原子 `±0.01 Å` × 3 方向 → `6·N_dev` 个位移。`N_dev = 54`（nap-sh）→ **324 个 SCF**。
- 为什么用中心差分而不是前向（3N）：本仓没有 DFPT 交叉验证，`K_ij = K_ji` 是**唯一可用的内部一致性检验**；前向差分的对称性误差是 O(δ)，中心差分是 O(δ²)。代价 2×，用数组并发摊掉。
- 位移幅度 `0.01 Å` 默认（报告 §3.3）；开发期验 `0.005/0.02 Å` 对 `∫T_ph dω` 的影响 < 3%（§9 T6）。
- 只位移器件区，其余原子冻结（`constrain_relaxation .true.`）——这就是「extended molecule + frozen leads」的标准做法。力要记录**全部原子**，这样 `V_{c,lead}`（器件区 ↔ lead 第一层）自然落在同一个矩阵里。
- **文献先例**：`0.01 Å` 中心差分正是 B 支的标准做法——Mosso 2019 的 SI 给出 `K_ij = [F_i(δ) − F_i(−δ)]/2δ`，Sadeghi 2015 把整个流程写全。A 支用 DFPT 解析 Hessian，**根本没有「位移幅度」这个参数**；我们与 B 支同构，与 A 支的差别**纯粹是引擎能力**（FHI-aims 221103 无 DFPT，§12.1），不是方法选择。
- **一个刻意不采用的技巧**：Mosso 把对角块用动量守恒补齐（`K_ii = −Σ_{j≠i}K_ij`）。本方案**保留自作用力测得的对角块**，因为 §3.6 要拿 `‖D·1‖` 当**独立的质量指标**；若对角块由 ASR 构造性给出，这个残差恒为 0，检查就退化成空检查。ASR 投影放在**测完残差之后**（§3.4），两者不冲突。

### 5.2 作业结构：数组作业

一个阶段 = 一个 sbatch 作业的假设在这里失效（324 个 SCF）。采用 **Slurm 数组**（`MaxArraySize = 1001`，够用）：

```bash
#SBATCH -p main
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=12G               # 必须降到 12G：见下
#SBATCH --time=06:00:00
#SBATCH --array=0-323%8         # 8 路并发 = 64 核 / 8 核
# 第 $SLURM_ARRAY_TASK_ID 个位移：写 disp_XXXX/{geometry.in,control.in} → aims → forces.json
```

**三条实测约束**（2026-11-23 在 gpu1 上核实，见 §12）：

1. **不要 `--gres=gpu:1`**。生产 run 的 `sub_aims.sh` 里带这个请求，但 `input_summary.json` 的 `aims_gpus` 是 `"cpu"`——实际是 `mpirun -np 8` 纯 CPU 跑，GPU 申请是模板残留。带着它，4 张卡只能开 4 路，白白砍掉一半并发。
2. **`--mem` 必须从 24G 降到 12G**。`aims.out` 实测峰值 tracked memory **871.9 MB/task**（92 原子）；24G × 8 路 = 192G > 节点 124.7G，Slurm 会只放行 5 个。12G × 8 = 96G 安全。
3. **只有一台机器**：`gpu-outside` 与 `gpu-inside` 是同一节点 `gpu1`（64 核 / 124.7G / 4×RTX 4090），不存在第二台可分摊。

好处：`jobwait.build_probe_command` 的 `squeue -h -j <array_id> -o "%i"` 天然列出全部 element，**等待逻辑零改动**；`collect` 只要求 324 个 `forces.json` 齐全，缺一个就 fail-closed。

### 5.3 力解析

`src/smj/phonon/forces.py` 新增 `read_forces(aims_out_path) -> np.ndarray (N,3)`：

- 锚 `Total atomic forces (unitary forces cleaned) [eV/Ang]:`，**取最后一个块**；
- 逐行 `^\s*\|\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)$`；
- 原子数必须等于 geometry 原子数，否则 raise（不截断、不补零）；
- 力不是有限值 → raise。

### 5.4 收敛与自检

`control_thermal.in` 在 `control_scf.in` 基础上**只加**：

```
compute_forces .true.
final_forces_cleaned .true.
sc_accuracy_forces 1.0e-3
k_grid 3 3 1
```

（**不是 `calculate_forces`**——手册关键词是 `compute_forces`，且明确要求「such as a finite-difference calculation of vibrational frequencies」时配合 `final_forces_cleaned`；见 [FHI-aims 手册](https://fhi-aims.org/uploads/documents/manuals/FHI-aims.071914.pdf) 与 [ASE 接口](https://docs.ase-lib.org/_modules/ase/calculators/aims.html)。报告此处写错，附录 A-3。）

去掉 `output aitranss`/`output mulliken`/`restart aims.restart`（热导不需要，省时）。

**lead 力常数**：`src/smj/junction/data/au111_3x3_l4.json` 冻结文件，由 `scripts/gen_au_lead_fc.py` 一次性生成（体相 Au 同一面内胞、周期 z、位移中央主层的 9 个原子 = 54 个 SCF）。文件带 `in_plane`/`layers`/`method`/`sha256` 元数据；配置里的面内胞与之不匹配 → **fail closed**。理由：lead 的力常数是**材料常数**，不是体系属性，与冻结 BDT/C4/2CN 基线同一设计哲学。

**FC 截断半径取 10 Å**（在文件里记 `fc_cutoff_A`）。依据：[Müller et al., PRB **106**, 195401 (2022)](https://arxiv.org/abs/2208.09171) 把 `r_c` 从 2 扫到 52 Å，最终取 **10 Å**；本器件 4 层沿 z 的跨度是 `3 × 2.3545 = 7.06 Å`，面内 Au(111) 最近邻 `a/√2 = 2.885 Å`，10 Å 覆盖约 3.5 个近邻壳层，足够。截断小于 7 Å 会切断层间耦合，必须 raise。

### 5.5 成本

单 SCF 成本按原子数线性外推（锚点：92 原子非周期 1563.9 s、68 个 SCF cycle，实测）：

```
t_SCF ≈ 1563.9 s × (N_atoms/92) × (N_kpts/1) × (cycles/68)
90 原子 × 5 个不可约 k 点，从零起收敛 ≈ 2.1 h/SCF
```

> ⚠️ **`N_kpts=5` 是未标定的上界**。`3×3×1` 的 9 个网格点在 C3v + 时间反演下大概率只剩 **3** 个不可约点（spglib 一跑就知道，B1 阶段顺手核）。若为 3，单 SCF 从零起 ≈ **1.28 h**，有 restart 时 ≈ 0.3–0.55 h，生产墙钟 **12–22 h** 而不是 20–36 h。**本节的成本按 5 个点给，是刻意取的保守侧**；标定时用 `disp_0000` 实测 cycle 数与 `aims.out` 里的 `| Number of k-points` 回填。

**但位移 SCF 不该从零起收敛**：`control_scf.in` 已带 `restart aims.restart`，集群上 **60 个 run 的 `aims_scf/` 都产出了 `aims.restart`**（291 MB，实测）。0.01 Å 位移对密度的扰动很小，从 restart 起步预计 15–30 cycle 而不是 68 → **0.5–0.9 h/SCF**。

| 场景 | 单 SCF | 324 个 / 8 路 | 说明 |
|---|---|---|---|
| 有 restart 复用 | 0.5–0.9 h | **20–36 h（≈1–1.5 天）** | 预期情形，需实测标定 |
| 无 restart（从零起） | 2.1 h | **85 h（≈3.5 天）** | 保守上界 |

**这是估算，不是实测**（k 点与周期的实际标度可能偏离线性，见 §10 R1）。**并发上限是 8，不是 16**——节点只有 64 核，8 核/作业（§5.2）。可用的降本档（按有 restart 的 0.5–0.9 h 估）：

| 档 | `device_layers` | `N_dev` | 位移数 | 估算墙钟（8 路） |
|---|---|---|---|---|
| 生产 | 2 | 54 | 324 | 20–36 h |
| 快速 | 1 | 36 | 216 | 14–24 h |
| 最小 | 1（且分子只位移 C/H，Au 只位移 3 个接触原子） | ~20 | 120 | 8–14 h |

> **第一个位移就是标定实验**：跑 `disp_0000` 后看 `aims.out` 的「Number of self-consistency cycles」，若 > 40 说明 restart 没生效（检查文件是否随位移目录复制过去），此时退回上表的保守上界。

> **不确定度的诚实边界**：`device_layers=1` 时器件区与 lead 的耦合只到一层，界面 FC 的截断误差未标定。默认 `2`；降到 `1` 必须在报告里显式标注。

**与文献的成本对照（唯一的公开 core-hour 数据）**：[Müller, Cuevas, Pauly & Nielaba, PRB **106**, 195401 (2022)](https://arxiv.org/abs/2208.09171) 逐字写「the dynamical matrix ... **is typically the most time-demanding computational step**」，并给了三档实测：NEMD 单次 288 core-h × 20 构型 × 14 温度 ≈ **80,640 core-h**；温度依赖 FC ≈ **13,440 core-h**；固定温度 FC ≈ **960 core-h**。

本方案：`324 SCF × 8 核 × 0.5–0.9 h ≈ 1,300–2,300 core-h`，与他们的**固定温度 FC 档同量级**（960 core-h），慢约 1.4–2.4×——考虑到我们是 PBE + 金属 + k 点、他们是同样的重元素体系，这个偏差是合理的。**这条对照的意义**：如果实测 core-hour 落到 10⁴ 量级，说明 restart 没生效或 k 点/层数设错了，应当停下来查，而不是硬等。

> 检索范围内**没有找到 Au–硫醇结声子 FC 的位移数或机时的公开数字**（UNVERIFIED）——A 支用 DFPT 一次拿全 Hessian，成本结构不可比；B 支论文都不报机时。所以上表是我们自己的估算，文献只能给量级锚点。

---

## 6. 代码落点

### 6.1 新增包 `src/smj/phonon/`（不进 `smj.junction` 内核）

| 文件 | 内容 | 关键签名 |
|---|---|---|
| `units.py` | 唯一换算常量 + `lambda_to_ev`、`ev_to_cm1`、`gq_pw_per_k` | `EV_PER_SQRT = 0.064654`；`def gq_pw_per_k(temperature_K: float) -> float` |
| `forces.py` | 力矢量解析 | `def read_forces(aims_out: Path) -> np.ndarray` |
| `fc.py` | 位移生成、FC 组装、ASR 投影、对称化 | `def make_displacements(n_atoms, device_indices, amplitude) -> list[dict]`；`def assemble_fc(disps, forces) -> np.ndarray`；`def apply_asr(K) -> tuple[np.ndarray, float]` |
| `dynamical.py` | `D = K/√(MMᵀ)`、特征值、虚频检测 | `def dynamical_matrix(K, masses) -> np.ndarray`；`def imaginary_modes(D, cm1_tol=-10.0) -> list[float]` |
| `negf.py` | 表面格林函数（decimation）、`T_ph(ω)` | `def surface_gf(omega_ev, d_layer, v_coupling, eta_ev) -> np.ndarray`；`def transmission(omega_ev, d_c, d_l, d_r, v_l, v_r, eta_ev) -> float` |
| `landauer.py` | `G_ph` 积分（唯一公式源） | `def thermal_conductance(omega_ev, t_ph, temperature_K) -> dict` |
| `device.py` | 周期板几何构造、原子角色划分、**S–H 去除** | `def build_device(molecule, cfg) -> tuple[Atoms, DeviceMeta]` |

> **`build_device` 内部必须先调 `geometry.strip_thiol_h`**（`src/smj/junction/geometry.py:160`），再挂板。`molecule_opt` 产出的是**中性二硫醇**（S 上带 H，nap-sh 的 18 个分子原子里就含 2 个 H）；直接把它贴到 Au(111) 上就是报告 §6 陷阱⑥「用 S–H 建模」的原样复现——Au–S 距离会多出 ~1.3 Å，`K` 与 `T_ph` 全错。电子链的 `electrode_build` 已经这么做（`strip_thiol_h: true`），热导支线必须复用同一口径，这也是 §2.2「不引入新的结构自由度」的具体含义。**验收挂在 T5 之外单列一条**：`build_device` 输出里 S 的邻居数必须恰为 1（Au），不含 H。

**与内核的边界**：只 `from smj.junction.geometry import read_aims_geometry, write_aims_geometry`。`smj/junction/*.py` **唯一改动**是 `write_aims_geometry` 加一个可选参数（现签名 `write_aims_geometry(path, atoms, constrained=())`，见 `geometry.py:68`）：

```python
def write_aims_geometry(path, atoms, constrained=(), lattice=None):
    """lattice=None 时输出与 0.2.2 逐字节相同。

    既有 docstring 已把这一点写成 W0-1 兼容约束（"既有 25 条 smj 记录的锚点回归靠它"），
    新增参数必须沿用该约束，并加一条 `lattice=None` 时逐字节相等的用例。
    """
```

`lattice` 非 None 时在文件头写 `lattice_vector`。此改动需重跑 `python scripts/gen_junction_provenance.py` 刷新 `PROVENANCE.json` 的 sha256，并跑 `gen_demo_reference.py --check` 确认数值基线不漂（AGENTS.md 硬要求）。

### 6.2 新增阶段（**不进 `STAGE_ORDER`**）

| 阶段 | 位置 | 依赖 | 提交 | 产出 |
|---|---|---|---|---|
| `thermal_device` | local | `molecule_opt` | 否 | `thermal/geometry.in`、`device_metadata.json`（逐原子角色 + `layer_index` + `principal_layer_id`） |
| `thermal_relax` | remote | `thermal_device` | 是（1 作业） | `relaxed_geometry.in`、`max_force_ev_per_A` |
| `thermal_fc` | remote | `thermal_relax` | 是（**数组作业**） | `disps.json`、`disp_*/forces.json` × 6N |
| `thermal_negf` | local | `thermal_fc` | 否 | `thermal.json`、`transmission_ph.csv`、`transmission_ph.png/pdf` |

`thermal_device` 补 G7：`electrode_metadata.json` 缺层号的问题在这里以新元数据文件解决（不改电子链的元数据）。

### 6.3 配置（新增 `thermal` 段）

`config.py` 三处：`_KNOWN_SECTIONS`（:469）加 `"thermal"`、新增 `_THERMAL_KEYS`（:485 附近）、新增 `ThermalConfig` 与解析（:718 附近，照 `transport` 段的 `_warn_unknown` 写法）。

```yaml
thermal:
  enabled: false            # 默认关闭；smj run 不受影响
  temperature_K: null       # null = 继承 transport.temperature
  electrode:
    in_plane: [3, 3]
    layers: 4
    device_layers: 2
    vacuum_A: 20.0
    k_grid: [3, 3, 1]
  fc:
    displacement_A: 0.01
    central_difference: true
    sc_accuracy_forces: 1.0e-3
    array_max_parallel: 8
  negf:
    eta_meV: 0.5
    omega_max_meV: 25.0
    omega_step_meV: 0.1
    asr: true
  lead_fc: "au111_3x3_l4.json"
```

### 6.4 契约：`thermal.json`（schema `smj.thermal.v1`）

```json
{
  "schema": "smj.thermal.v1",
  "source": "thermal_negf",
  "temperature_K": 300.0,
  "G_ph_pW_per_K": 0.0,
  "GQ_pW_per_K": 283.9,
  "N_eff": 0.0,
  "ZT_total": null,
  "k_e_pW_per_K": null,
  "omega_max_meV": 25.0,
  "points": 251,
  "eta_meV": 0.5,
  "asr_residual": 0.0,
  "fc_asymmetry": 0.0,
  "imaginary_modes_cm1": [],
  "band_top_meV": 19.5,
  "t_ph_above_band_max": 0.0,
  "lead_fc_file": "au111_3x3_l4.json",
  "lead_fc_sha256": "…",
  "warnings": []
}
```

配套 `transmission_ph.csv`（`omega_meV,t_ph`）与 `conductance_cum.csv`（`omega_meV,g_cum_pW_per_K`）。写入沿用 `analysis.py:_write_json` 的 `.tmp` + `os.replace` + `allow_nan=False`。

`ZT_total` 与 `k_e_pW_per_K` 从 `junction_analyze/results.json` **只读**取得；取不到就写 `null` + warning，**不阻塞**（热导可以独立于电子结果成立）。

### 6.5 判据（`criteria.py` 新增 10 条）

| name | kind | band | source |
|---|---|---|---|
| `ph_transmission_present` | engineering | — | — |
| `ph_relax_converged` | engineering | `max\|F\| ≤ 1e-2 eV/Å` 且 `sc_accuracy_forces` 已设 | — |
| `ph_band_top_zero` | physical | `t_ph_above_band_max < 1e-6` | Lynn et al., PRB **8**, 3493 (1973), DOI [10.1103/PhysRevB.8.3493](https://doi.org/10.1103/PhysRevB.8.3493)（Au 带顶 4.71 THz = 19.5 meV） |
| `ph_quantum_bound` | physical | `G_ph < 3·G_Q(T)` | [Schwab et al., Nature 404, 974 (2000)](https://doi.org/10.1038/35010065) |
| `ph_acoustic_limit` | physical | `T_ph(ω→0) ≤ 3`（**k 平均口径**，见 §3.4 警示框） | 同上（三条声学支） |
| `ph_asr_residual` | engineering | 相对残差 `< 1e-3`（判投影前的值，D13） | — |
| `ph_fc_symmetry` | engineering | `‖K−Kᵀ‖/‖K‖ < 5%` | — |
| `ph_imaginary_modes` | physical | 无 `< −10 cm⁻¹` | 报告 §6 陷阱④（−10 cm⁻¹ 噪声容限） |
| `ph_lead_cell_match` | engineering | 冻结 lead FC 的面内胞 == 器件配置 | — |
| `ph_grid_independence` | engineering | 步长减半 `ΔG/G < 5%`（D14，原 1% 被实测求积噪声否掉） | — |

`physical` 四条都带 `source`，满足 `lint_registry()` 的铁律；`KNOWN_GAPS` 只减不增。

> **命名说明**：不要叫 `ph_te_present`——`te` 在本仓专指电子的 `TE.dat`（`aitranss_te` 阶段产物），热导侧没有这个文件，`T_ph(ω)` 写在 `transmission_ph.csv` 里。`ph_relax_converged` 对应 §3.6 那条「力收敛 `max|F|` 缺失或未达 → raise」——这是**力常数精度的上游闸门**，`thermal_relax` 的几何没收敛，后面 324 个位移全部白算，必须有独立判据，不能靠 `scf_converged` 顶替（后者只管 SCF 残差）。

### 6.6 报告与 CLI

- `report.py`：新增 `## 5. 热导`（现有 `## 4. 输运结果` 在 1010-1043 之后）。表行：`G_ph`、`G_Q`、`N_eff`、`ZT_e`、`ZT_total`、`T_ph(ω>带顶)`、ASR 残差、虚频数、lead FC 文件 sha。图：`T_ph(ω)` + 累积 `G_cum(ω)` 双轴。
- `report.py` 的 `_transport_summary`（:391）不动；热导单独一个 `_thermal_summary`。
- `cli.py`：新增子命令 `smj thermal`（`--run/--prepare/--submit/--wait/--collect`，照 `cmd_stage`:637 的形状），以及 `smj thermal --analyze-only`（只重跑 NEGF，不重跑 DFT）。
- `state.py`：`STAGE_ORDER` **不变**；新增 `THERMAL_STAGES = ("thermal_device", "thermal_relax", "thermal_fc", "thermal_negf")` 供 `smj thermal` 按序调用。

---

## 7. 批次与里程碑

| 批次 | 内容 | 涉及 | 完成判据 |
|---|---|---|---|
| **B1 数值内核** | `units.py`/`negf.py`/`landauer.py`/`dynamical.py` + 解析基准测试 | 新包；无 DFT | T1–T4 全绿；**不碰 FHI-aims** |
| **B2 几何与力** | `device.py`、`write_aims_geometry(lattice=)`、`forces.py`、`control_thermal.in` | 内核 1 处 + 新包 | T5；`gen_demo_reference --check` 不漂；provenance 刷新 |
| **B3 阶段与配置** | 4 个 stage、`ThermalConfig`、数组作业脚本 | `stages/`、`config.py`、`state.py` | `smj thermal --prepare` 在 nap-sh 上 dry-run 通过 |
| **B4 端到端** | nap-sh 真实跑通；η/位移/层数/基组四项开发期验收；（可选）D9 文献对账 | 集群 | T6–T9；`thermal.json` 产出；T14 视用户是否同意而定 |
| **B5 报告与判据** | `thermal.json` schema、10 条判据、report 节、CLI 打磨 | `criteria.py`、`report.py`、`analysis.py` | T10；620 既有用例全绿 |

B1 可以完全离线做，是**风险最低、信息量最大**的一批：解析基准一旦立住，后面所有的物理数字都有了对照。建议先做 B1，用它产出的 `T_ph(ω)` 在**人造力常数**（如 1D 链、二聚体）上验证，再进入 B2。

---

## 8. 决策点

| # | 问题 | 建议 | 理由 |
|---|---|---|---|
| **D1** | 热导用新几何（周期板）还是复用金字塔团簇？ | **新几何** | §4.1；复用会让 `G_ph` 变成唯象 σ 的函数 |
| **D2** | 自写 NEGF 还是接 GOLLUM/PHtrans/QuantumATK？ | **自写** | ① 三个包都要求换电子结构引擎（§4.2）；② **文献侧：这个领域没有一支主线用现成声子输运包**——Pauly 组与 Sadeghi 组的输运都是自写码，GOLLUM 2014 原文无声子模块，PHtrans 实际输入是 GULP，QuantumATK 模块零文献（§4.2、附录 A-10～A-12）；③ 自写 NEGF 的核心 < 200 行，且有解析基准可钉（T2–T4） |
| **D3** | 是否申请把集群 FHI-aims 从 221103 升到 ≥240507 用 DFPT？ | **不升**（结论已从「并行提申请」改为否定，理由见 §12.1） | ① 官方文档：周期体系的 DFPT「only limited support (**LDA functional**, and systems that do not require large memory)」——我们的 slab 用 PBE，**升了也用不上**；② 官方对「其他泛函 + 大体系」推荐的恰恰是 FHI-vibes 有限位移，即本方案路线；③ 升级风险不对称：`aitranss 071813` 与全部冻结基线绑在 221103 上。若将来要做，只能**旁装**（不替换 221103）并在小团簇上单独验证 |
| **D4** | lead 力常数冻结入库还是每次现算？ | **冻结入库** | 是材料常数；与 BDT/C4/2CN 冻结基线同一哲学；省 54 个 SCF/次 |
| **D5** | 中心差分（6N）还是前向（3N）？ | **中心差分** | 唯一可用的内部一致性检验是 `K_ij=K_ji`（无 DFPT 交叉验证）；2× 成本由数组并发吸收 |
| **D6** | 基组用 `light` 还是 `tight`？ | **生产 `light`（与电子链一致），开发期做一次 light vs tight 对比** | 换基组会让热导与电子结果不可比；但 FC 对基组比能量敏感，必须量一次偏差并写进报告 |
| **D7** | `G_ph` 与 `ZT_total` 进不进 `smj review` 的 gate？ | **进，但 `ZT_total` 为 `engineering` 非门控** | `G_ph` 的物理判据可以门控；`ZT_total` 依赖电子结果，取不到时为 `null`，不能门控 |
| **D8** | 温度：`G_ph` 用 `transport.temperature` 还是允许列表？ | **默认继承，允许 `thermal.temperature_K` 覆盖** | 报告 §6 陷阱⑪：实验是电极温差，理论用平均温度；单一温度避免混淆 |
| **D9** | 要不要跑一个「文献对账」控制体系（octanedithiol 或 BDT）？ | **建议做一次，但必须用户明确同意，且结果不混进生产结论** | 我们是**唯一用 FHI-aims 做分子结 `κ_ph` 的组**，没有第三方数值可对账（§10 R7）。文献基准很窄但一致：Au–烷二硫醇–Au 在 300 K 应落在 **15–25 pW/K**（Cui 实验 ~20；Mosso 理论 19、实验 20±6；Klöckner 烷链 15–45）。跑一个 C8 或 BDT 就能一次性检验整条链（力常数→NEGF→积分）。**这是对照臂，按 AGENTS.md 必须先说明理由和代价、经用户同意才跑**；且报告里单列，不得与 nap/az/ind 的结果并列 |
| **D13** | ASR 残差怎么算、阈值多少？ | **判相对残差 `max\|Σ_j K_ij\|/max\|K\| < 1e-3`，先测再投影（C 方案）** | 绝对门限 1e-6 eV/Å² 被实测噪声否掉：有限差分 + 基组噪声本身就在 1e-3–1e-1 eV/Å²。且投影必须作用在**器件行子集**上——lead 两层块按构造就不闭合（完整 6 层胞 1.05e-3、2 层块 0.60），拿它判会永远失败 |
| **D14** | 网格独立性门限 1% 达不到，怎么办？ | **重标为 5%，并承认 η 才是主敏感量** | 2026-11-25 实测（nap-sh，同一份 `force_constants.npy`）：粗化一半的 `ΔG/G` 在 η=0.05/0.1/0.2/0.3 meV 下分别是 1.99%/1.58%/2.46%/2.21%，加密网格不降反在 1–2% 摆动 → 梯形积分在峰宽 ~η 的窄共振上有 2.0%–2.5% 的求积底噪。真正决定数值的是 η：η 0.5→0.05 meV 使 `G_ph` 由 10.9 升到 13.7 pW/K（+25%），单调外推 η→0 ≈ 13.8。故 `G_ph` 一律按「η→0 外推值 + η 敏感性作为误差带」报告，门限只保留来抓「网格明显不够」（>5%） |

---

## 9. 验收

| # | 项 | 命令/证据 | 期望 |
|---|---|---|---|
| T1 | 单位常量 | `pytest tests/phonon/test_units.py` | `EV_PER_SQRT=0.064654`；`0.064654/4.135667696 = 15.6333 THz`，与 phonopy `VaspToTHz=15.633302` 相对误差 < 1e-6；`λ=0.1 → 20.4 meV`（≈Au 带顶 19.5 meV） |
| T2 | 一维链解析基准 | `pytest tests/phonon/test_negf_analytic.py` | 单原子链 `T_ph(ω)` 在通带内 = 1、带外 = 0，与解析色散逐点吻合（< 1e-10） |
| T3 | 量子极限恒等式 | 同上 | `T_ph≡3` 时 `G_ph(T) = 3G_Q(T)`，300/100/4 K 三档 < 1e-10 相对 |
| T4 | decimation vs 直接求逆 | 同上 | 200 层有限堆叠直接求逆 vs decimation，`g_surf` 相对误差 < 1e-8 |
| T5 | 力解析 | `pytest tests/phonon/test_forces.py`，夹具 `tests/data/aims_force_block_20atoms.out`（**新增**：从 `aims_tip_relax/aims.out` 裁出表头 + 末 2 个力块，约 40 行） | 20 原子、取第 27 块（末块）；逐分量与文件一致；缺块/非有限值/原子数不符 → raise |
| T5b | S–H 去除 | `pytest tests/phonon/test_device.py` | `build_device` 后 S 的邻居恰为 1 个 Au、无 H；分子原子数 18（nap-sh） |
| T6 | 位移幅度 | 开发期一次性（nap-sh） | `0.005/0.01/0.02 Å` 下 `∫T_ph dω` 变化 < 3% |
| T7 | η 与网格 | 开发期一次性 | `η ∈ {0.05,0.1,0.2,0.5}` 与步长减半；**实测结论（2026-11-25）**：η 单调收敛、是主敏感量（10.9→13.7 pW/K），网格恒有 2% 底噪 → 门限改 5%（D14） |
| T8 | 电极层数 | 开发期一次性 | 4 vs 6 层，`G_ph` 变化 < 5%（报告 §3.5 要求） |
| T9 | 基组 | 开发期一次性 | light vs tight 的 `G_ph` 偏差记录在案，写进报告 |
| T10 | 端到端 | `smj thermal --run`（nap-sh） | `thermal.json` 产出，10 条判据全绿，`ZT_total` 与 `ZT_e/(1+G_ph/k_e)` 自洽 |
| T11 | 电子链零回归 | `python -m pytest -q` | 620 + 新增用例全绿；`2 deselected` 不变 |
| T12 | 数值基线 | `python scripts/gen_demo_reference.py --check` | 不漂 |
| T13 | provenance | `python scripts/gen_junction_provenance.py` | `PROVENANCE.json` 的 sha256 与 `git status` 一致 |
| T14 | **文献对账**（可选，**须用户明确同意**，见 D9） | `smj thermal --run`（octanedithiol 或 BDT） | `G_ph(300 K)` 落在 **15–25 pW/K**。超出 3× 按 R4 逐条查，**不许调参去凑** |

> ⚠️ **T14 的对账带只对「DFT 谐波理论值」成立，不要拿实验值去比**。15–25 pW/K 是从 B-4 里那批**同为简谐 Landauer** 的计算结果取的（Klöckner 烷链 15–45、苯二胺 25.24、Sadeghi 烷烃 25.4/33.4/30.3、Mosso 理论 22/19）。Mosso 的**实验**值是 ODT **37±5**、OPE3 **23±5**——比他们自己的 DFT 理论值高 1.7×。我们的方法在结构上就少了非谐（Li et al. 2021 给 ~15%），预期**系统性偏低**，所以拿 37 当靶子会必然「差 2 倍」并误判为 bug。若要与实验比，必须在报告里写明「简谐 NEGF 与实验的已知偏差来源」，而不是当成验收线。

---

## 10. 风险与回退

| # | 风险 | 概率 | 回退 |
|---|---|---|---|
| R1 | 周期板 SCF 比线性外推慢很多（金属 + k 点 + 电荷晃动） | 中 | 先跑 1 个位移做**实测标定**，再定并发数与档位；最坏退到 §5.5「最小档」 |
| R2 | 力噪声压不到 1e-3 eV/Å，`K_ij=K_ji` 超 5% | 中 | 收紧 `sc_accuracy_forces`/`sc_accuracy_rho`；或改中心差分 + 提高 δ 到 0.02 Å |
| R3 | 器件区出现大量虚频（金属界面软模） | 中 | 报告里显式列出；`<−10 cm⁻¹` 才 fail；同时检查 `thermal_relax` 是否收敛到 |F|<0.01 |
| R4 | `G_ph` 与实验 20 pW/K 差 3× 以上 | 中 | **不要调参去凑**。按报告 §3.5 逐条查：带顶归零、ASR、电极层数、接触构型、基组。差得离谱就是有 bug，不是「理论就该差这么多」 |
| R5 | 数组作业打爆队列配额 | 低 | `array_max_parallel` 默认已定 8（§6.3，= 64 核 / 8 核每作业）；队列仍紧张时降到 4，或改单作业内循环 + 逐位移哨兵 |
| R6 | 整个热导支线失败 | 低 | 全部为新增文件；`thermal.enabled` 默认 false；删 stage 注册即可完全回退，电子链零影响 |
| R7 | **没有 FHI-aims 基的分子结 `κ_ph` 先例，无第三方数值可对账** | **高（确定存在）** | 这是文献核实后的**新风险**，无法消除，只能分层兜底：① 解析基准（T2–T4）钉死数值内核；② 内部一致性（ASR、带顶归零、电极层数、弹道不变性）逐条过；③ D9 的文献对账（需用户同意）；④ 报告里**显式标注「无同类先例，本结果为首次实现」**，不冒充「与文献吻合」。**A 支与 B 支的公式相同、引擎不同**——这既是风险来源，也是可信度依据 |

---

## 11. 与项目规则的关系

- **「只跑用户明说的」**：热导支线默认关闭，`smj run` 的六阶段链逐字节不变；只有用户显式 `smj thermal --run` 才跑。T6–T9 四项是**开发期一次性验收**，不是每次 run 的控制臂，且只在 nap-sh 上跑一次。
- **生产默认口径**：S 的质子化、吸附位、`pin_contact`/`tip_direction`/`anchor_au_s_target` 一律沿用电子链配置，热导支线**不引入新的结构自由度**。
- **`inconclusive` 不得包装成 `passed`**：`ph_imaginary_modes` 有软模但 < −10 cm⁻¹ 之外、或 `ZT_total` 因缺电子结果而为 `null` 时，review 出 `inconclusive`，不冒充 `pass`。
- **协议不一致必须标注**：热导器件几何（周期板、约 90 原子、3×3）与电子链几何（金字塔团簇、nap-sh 92 原子）**不是同一个体系**。`report.md` 必须在热导节显式写明这一点，禁止把两者的 `G_ph` 与 `G/G0` 当作同一几何的两个性质并列。
- **内核改动**：本方案只改 `src/smj/junction/geometry.py` 一处（可选 `lattice=`），必须刷新 `PROVENANCE.json` 并跑 demo 对账；`scripts/crosscheck_agenticdft.py` 不进 CI。
- **文献对账是控制臂，不是默认结果**：D9/T14 的 C8 或 BDT 属于 AGENTS.md 明令「不要默认加跑」的对照计算，必须**先说明理由与代价、经用户明确同意**才跑；结果在 `report.md` 里**单列一节**，不得与 nap/az/ind 并列，更不得混进「默认结果」。
- **目标分子无文献值，不得声称「与文献吻合」**：萘/薁/茚二硫醇的 `κ_ph` 在检索范围内**没有任何已发表数值**（§B-4 末注）。它们的可信度只能建立在内部一致性上；报告必须写明「无同类先例，本结果为首次实现」。

---

## 12. 软件与运行环境（2026-11-23 在 gpu1 上实测）

**机器就是一台**：`gpu-outside` 与 `gpu-inside` 都解析到 hostname `gpu1`（`4.18.0-553`，CentOS 8.10），64 核 / 124.7 GB / 4×RTX 4090，Slurm 分区 `main`（唯一节点，`MaxArraySize=1001`，7 天上限）。

| 软件 | 状态 | 位置 / 版本 |
|---|---|---|
| **FHI-aims** | ✅ **221103**（CUDA 12.8 + oneAPI 2025.3） | `/data/apps/fhi-aims.221103-cuda12.8-oneapi2025.3/bin/aims.221103.scalapack.mpi.x` |
| ↳ `compute_forces` / `final_forces_cleaned` | ✅ 支持 | 二进制内确认（`read_control_module_mp_compute_forces_`、`runtime_choices_mp_final_forces_cleaned_`） |
| ↳ `lattice_vector` / `k_grid` / `use_dipole_correction` | ✅ 支持 | 二进制内 69 / 74 / 2 处命中 → 周期板可行 |
| ↳ **DFPT**（`atomic_pert_response`） | ❌ **无** | 二进制零命中，印证 221103 < 240507（附录 A-4） |
| aitranss | ✅ 071813 | 同目录 `aitranss.071813.openmp.ifx2025.3.x` |
| **Python** | ✅ 3.12.12 | miniforge 25.11.0-0 |
| numpy / scipy | ✅ 2.5.3 / 1.17.1 | `~/.local` |
| **phonopy** | ✅ 4.5.0，**含 `interface/aims.py`** | 原生 FHI-aims 接口可用（`phonopy -d --aims`） |
| spglib / ase | ✅ 2.7.0 / 3.29.0 | — |
| fhi-vibes | ❌ 未装 | **不需要**：phonopy + ase 已覆盖 |
| **SIESTA** | ✅ 5.4.2 | `/data/apps/siesta-5.4.2-oneapi2025.3-impi2021.17-mkl-dftd3` |
| ↳ `phtrans` / `tbtrans` | ✅ 二进制存在 | 直接跑会缺 MKL → 需先 `module load` oneAPI（`env.mkl_lib` 已配 `/opt/intel/oneapi/2025.3/lib`） |
| ↳ `MD.TypeOfRun FC` / `MD.FCDispl` | ✅ 支持 | 二进制内确认 |
| **GOLLUM** | ❌ **未装** | 报告 §8 的「零预算路线」需先自行编译 |
| QuantumATK | ❌ 无（商业） | — |
| LAMMPS / GPUMD | ✅ 2025.12（含 kokkos-cuda）/ 5.1 | 将来做 RNEMD/HNEMD 校验用 |
| Gaussian 16 / VASP 6.5/6.6 | ✅ | 与本方案无关 |

**对本方案的三条直接结论**：

1. **B1（数值内核）零依赖**：numpy + scipy 已就位，可在登录节点直接跑完 T1–T4，不占作业、不等队列。这也是把 B1 排在最前面的又一理由。
2. **B2–B4 需要的 FHI-aims 能力全部具备**：`compute_forces` + `final_forces_cleaned` + `lattice_vector` + `k_grid` 都在 221103 里。**唯一的缺口是 DFPT**，而方案本来就走有限位移（D3），所以**不阻塞**。
3. **候选 B 的软件障碍比预期低**：SIESTA 5.4.2 + phtrans 竟然已经装在这台机器上（§4.2 曾假设要自己搭）。但改用 B 仍要重建第二套电子结构链路（基组、赝势、`MD.TypeOfRun FC` 产物），并把已验证的 BDT/AZ/2CN 基线推倒重来；加上 PHtrans 官方自认「not as tested as TBtrans」，**结论不变：仍选 A**。这条只是让「万一 A 的自写 NEGF 走不通」时有一个已装好的退路。

### 12.1 要不要把 FHI-aims 升到 ≥240507 用 DFPT？

**结论：不升。** 四条理由，前两条是「升了也没用」，后两条是「升了有风险」。

**① 周期体系的 DFPT 只支持 LDA。** 官方 DFPT 论文（arXiv:2501.16091，FHI-aims Roadmap 章节）原话：

> "For periodic phonons, only limited support (**LDA functional**, and systems that do not require large memory) is provided in the current implementation."

我们的器件是**二维周期 slab**，用的是 **PBE + vdW**（与电子链一致，D6）。换 LDA 会让 `G_ph` 与电子结构口径脱节，而 Au–S 键刚度对泛函是有响应的。**升到 240507 也过不了这一关。**

**② 官方对「其他泛函 + 大体系」推荐的正是本方案路线。** 同一段紧接着写：

> "A more comprehensive and very efficient handling of the calculation of phonons with different functionals and larger system sizes is provided through the FHI-aims connection to the **FHI-vibes** infrastructure."

FHI-vibes 就是 phonopy + 有限位移。也就是说：对 PBE、对 90 原子 slab，**官方给的答案就是有限位移**，不是 DFPT。

**③ 金属的扩展没点名声子。** 论文提到「Recently, we have extended the support to allow 'electric response' and 'electronic friction' calculations of fractionally occupied systems (including metallic systems)」——只列了电响应与电子摩擦，**没列 atomic response（声子）**。Au 电极是金属，这条路是否通还得单独验证。

**④ 收益本来就不重要，风险却不对称。** 就算 DFPT 可用，它把 `6N` 次 SCF 换成 `3N` 次线性响应，是 2× 的差别；而 §5.5 显示**重启密度后 324 个位移本身只要 20–36 h**，2× 不是决定性收益。反过来：`aitranss 071813` 与全部冻结基线（BDT/AZ/2CN 的 `G/G0`、`docs/validation-data/`、`PROVENANCE.json`）都绑在 221103 上；`/data/apps` 下的安装属 `wanghc`，升级要动别人的东西，且要重新验证 aitranss 的输入格式兼容性。**这是拿确定的基线风险换不确定的 2× 提速。**

> 若将来确实要用 DFPT：只能**旁装**（`fhi-aims.24xxxx` 与 221103 并存，用 `env.aims_bin` 切换），并在**小团簇**（非周期、PBE）上先验证 `atomic_pert_response DFPT` 能否给出正确的 Hessian；不要替换默认安装。

### 12.2 这些软件在文献里实际被怎么用

文献调研（[thermal-conductance-literature-2026-11-23.md](thermal-conductance-literature-2026-11-23.md)）与 §12 的机器清单逐项对照后，结论比预期整齐：**gpu1 的配置覆盖了本方案需要的全部软件，也覆盖了 B 支（SIESTA）与 C 支（LAMMPS/GPUMD）的全部软件**；缺的只有 A 支的 TURBOMOLE（其能力由 FHI-aims 等价提供）和两个在文献里最不成立的选项（GOLLUM、QuantumATK）。

| 软件 | gpu1 状态 | 文献里的**实际**用法 | 对我们的含义 |
|---|---|---|---|
| **FHI-aims** | ✅ 221103 | **检索范围内没有一篇分子结 `κ_ph` 用它** | 我们是首次实现；这是 §10 R7 的来源 |
| TURBOMOLE | ❌ 未装 | A 支主力（Pauly 组），**DFPT 解析 Hessian** | 引擎不同、**方法相同**：团簇/周期 + 自写 NEGF + decimation |
| **SIESTA** | ✅ 5.4.2 | B 支主力（Sadeghi–Lambert 组），`0.01 Å` 中心差分；**声子输运不点名软件** | 公式与流程和本方案 §5 **逐字相同** → 最好的方法对照 |
| ↳ PHtrans | ✅ 二进制在 | **实际输入是 GULP**；Papior 自认 SIESTA 路线「extremely costly」 | 不是「A 的退路」，是「另一套生态」（§4.2） |
| **GOLLUM** | ❌ 未装 | 电子热导（Wiedemann–Franz）；**2014 原文无声子模块** | 报告 §8 的「零预算首选」不成立（附录 A-10） |
| QuantumATK | ❌ 商业 | 模块存在，**零篇分子结 `κ_ph` 论文**；官方教程是 CNT | 「买来就能用」是空集（附录 A-11） |
| **phonopy** | ✅ 4.5.0，含 `interface/aims.py` | 体相 / 超胞 FC；**VASP+phonopy 在分子结这条线上缺席** | 只当**交叉校验工具**（T1 的 `VaspToTHz`），不作主链 |
| LAMMPS / GPUMD | ✅ 2025.12 / 5.1 | C 支：LAMMPS + NEMD 家族，**没有一篇 HNEMD**；力场是文献 Morse，无 DFT 拟合 | 将来做非谐/校验档的入口，不是本版（§2.2） |
| Gaussian 16 | ✅ C02 | 历史上有 Gaussian03 + 自写 NEGF（Lü & Wang, PRB **78**, 235436 (2008)） | 同样受「无 PBC」限制，不构成更好的选择 |

**三条结论**：

1. **技术选型不变**。文献核实没有推翻任何一个决策，只是把每条否决理由从「推断」升级为「引文」：GOLLUM 无声子、PHtrans 吃 GULP、QuantumATK 零文献、VASP+phonopy 无先例。
2. **唯一实质变化是风险显性化**。我们成了「唯一用 FHI-aims 做分子结 `κ_ph` 的组」，没有第三方数值可对账 → 新增 §10 R7 与决策点 D9（文献对账控制体系）。
3. **不需要装任何新软件**。B1 用 numpy/scipy；B2–B4 用 FHI-aims 221103 + Slurm；唯一可选的是 `module load` oneAPI（只在走 PHtrans 退路时才需要）。

---

## 附录 A：调研报告勘误（13 条）

| # | 报告位置 | 报告写法 | 核实结论 | 依据 |
|---|---|---|---|---|
| **A-1** | §3.4、§9 | `ħω[eV] = 0.2543·√λ`，校验点「λ=0.010→25.4 meV≈Au 带顶」 | **错**。正确 `0.064654·√λ`；`0.2543² = 0.064668`，即多开了一次方，偏大 **3.933×**。正确的校验点是 `λ=0.1 → 20.4 meV ≈ Au 带顶 19.5 meV`。照抄会把 C–C 伸缩算成 3200–4200 cm⁻¹（X–H 区） | CODATA 2018 自推 + phonopy `VaspToTHz=15.633302` 交叉核对 + 两个物理地标 |
| **A-2** | §1 | `G=(ħ²/2πk_BT²)∫ω²T_ph eˣ/(eˣ−1)²dω` | **对**（`ω` 用 rad/s 时）。等价形式 `G=(k_B²T/h)∫T_ph(x)x²eˣ/(eˣ−1)²dx`。但 §9 若写成「`(k_B/h)∫T_ph x²…`」则量纲错（少一个 `k_BT`，得 J/K 而非 W/K） | 自推：`ħω·∂n/∂T = k_Bx²eˣ/(eˣ−1)²`，`dω=(k_BT/ħ)dx` |
| **A-3** | §3.3、§8 | FHI-aims 用 `calculate_forces` | **关键词是 `compute_forces`**；做有限差分力常数时还须 `final_forces_cleaned .true.`。手册原文：「This is necessary for use with external tools that require forces, such as a finite-difference calculation of vibrational frequencies」 | [FHI-aims 手册](https://fhi-aims.org/uploads/documents/manuals/FHI-aims.071914.pdf)、[ASE Aims 接口](https://docs.ase-lib.org/_modules/ase/calculators/aims.html)、官方 tutorial 的 `control.in` |
| **A-4** | §3.3 | 金属不支持 DFPT，只能用有限位移 | **过时，但结论仍成立**。FHI-aims 自 **240507** 起有 DFPT（`atomic_pert_response DFPT`），然而「For periodic phonons, only limited support (LDA functional, and systems that do not require large memory)」——**周期体系只支持 LDA**，且官方对 PBE/大体系推荐的就是 FHI-vibes 有限位移。本仓集群是 221103（二进制内 `atomic_pert_response` 零命中），且**不该升**（§12.1） | [arXiv:2501.16091](https://arxiv.org/abs/2501.16091)、手册 §3.39/§3.40 |
| **A-5** | §5 | Mosso 2019 的「拉伸区 16–21 pW/K、松弛区 22–33 pW/K」 | **归因错**。这两个区间是 **Cui 2019 自己的计算结果**，对应两个**都拉伸**的构型（AB 更接近断裂、G 更低；AG 拉伸较少），不是「拉伸 vs 松弛」，也不是实验值。Mosso 的实测值：ODT (37±5)、OPE3 (23±5) pW/K（来自作者学位论文，正文付费墙） | Cui et al., [Nature **572**, 628 (2019)](https://doi.org/10.1038/s41586-019-1420-z)；Mosso et al., [Nano Lett. **19**, 7614 (2019)](https://doi.org/10.1021/acs.nanolett.9b02089)；Mosso 学位论文 DOI 10.5283/epub.41047 |
| **A-6** | §5、§6⑧ | 「电子热导通常 <5 pW/K」 | **只对弱耦合烷链成立**。本仓实测：nap-sh **116.2**、az **224.0**、bdt **165.6 pW/K**（G/G0 0.19–0.38）；只有 2cn **1.24**、hdt **0.91 pW/K** 才 <5。Cui 的 5.7 pW/K 对应 C2 的 G/G0≈0.01。**这条直接决定 `ZT` 的框架**（§2.3） | 本仓 `results.json` 实测 + WF 自洽核对 |
| **A-7** | §0、§6⑨ | Au 声子带顶 ~25 meV（6 THz），同时又说 θ_D≈170 K→14.6 meV | **两个都偏**。实测 fcc-Au 最高支在 L 点 **4.71 THz = 19.5 meV**；θ_D 是低温德拜平均截止，**不是带顶**，两者并列本身自相矛盾。代码必须保证 `T_ph(ħω>19.5 meV) ≡ 0` | Lynn, Smith & Nicklow, [PRB **8**, 3493 (1973)](https://doi.org/10.1103/PhysRevB.8.3493) |
| **A-8** | §8 | `blaschke/phonon_transport`（Markussen 宽带电极模型） | **仓库不存在**：`github.com/blaschke` 是无关用户，`/blaschke/phonon_transport` 返回 404；未找到该模型的公开实现。该条应删或标 **UNVERIFIED** | 直接访问核实（2026-11-23） |
| **A-9** | §7 产出② | 「薁非对称性可能带来 10–30% **方向性**差异」 | **原理上不成立**。Caroli 公式 `T_ph=Tr[G^rΓ_LG^aΓ_R]` 在 L↔R 交换下不变（`Γ_α` 厄米、迹的循环性），且线性响应热导由 Onsager 互易必然对称。非对称结只能改变**绝对值**，不能整流。若确实要方向性数据，须上非谐 NEGF 或带温差的 NEMD | 自证 + 线性响应互易；报告 §0 的「两端接触刚度差」表述是**对的**，§7 把它升级成「方向性」就错了 |
| **A-10** | §8 | GOLLUM 能直接出热导（含声子） | **2014 原始文献里没有声子模块**。结论段逐字：「Further developments are in the pipeline, **including the incorporation of phonon transport**」——标题中的 "thermal" 指**电子热导**（Wiedemann–Franz）。声子模式是后续版本才加的（手册 mode 5，读 SIESTA `.FC`），但**没有一篇分子结 `κ_ph` 论文点名用它**（含 Mosso 2019） | [Ferrer et al., NJP **16**, 093029 (2014)](https://arxiv.org/abs/1502.04966) 结论段；多轮检索 |
| **A-11** | §8 | QuantumATK 用 `PhononTransmissionSpectrum` 可直接算 | **模块确实存在**（`thermalConductance(phonon_temperature=...)` 返回 nW/K），但**检索范围内零篇论文**用它算单分子结 `κ_ph`；官方唯一教程是 CNT 同位素热电。反例：Markussen 做分子结声子时用的是 **GPAW + 自写码** | QuantumATK 官方文档 + 论坛；web / arXiv / OpenAlex / Crossref / DBLP 多轮检索 |
| **A-12** | §8 | 「已有 VASP → VASP+phonopy 转 SIESTA `.FC` → GOLLUM」 | **这条链在文献里没有先例**。能确证的「DFT 出力常数 + 自写输运」实例，DFT 侧是 TURBOMOLE / SIESTA / GPAW / QE / Gaussian03，**没有 VASP**。原因是 phonopy 超胞 FC 只解决中心区，分子结还需要**半无限电极的表面格林函数与力常数拼接**。报告已提到 `FORCE_CONSTANTS` 与 `.FC` 格式不同，但**未把转换脚本 + 对称性/ASR 自检计入工作量** | 检索（arXiv / OpenAlex / Crossref）+ phonopy/SIESTA 文档 |
| **A-13** | §8 | 「SIESTA+PHtrans」可行 | **包存在，但实际输入不是 SIESTA**。Papior 在 siesta-l（2019-03-31）原话：用 SIESTA 做声子透射「**extremely costly** since it requires a huge system to ensure the phonons are 'relaxed' in the electrode region. **This is why we (in the paper) did phonon transmissions using GULP**」；sisl 官方例子即 GULP→PHtrans。他 2019 年还在说从 DFT 提取动力学矩阵「**still work to do**」。官方另有「**NOTE: PHtrans is not as tested as TBtrans.**」 | [siesta-l 2019-03-31](https://lists.ugr.es/pipermail/siesta-l/2019-March/)；sisl `examples/ex_03.py`；[PHtrans 手册](https://siesta-project.gitlab.io/-/documentation/siesta-docs/-/jobs/14622928758/artifacts/public/reference/tbtrans.html) |

补充：**A-10～A-13 全部集中在报告 §8「软件」一节**——这不是巧合，而是该节依赖厂商文档与二手转述、没有逐条回原文核实；A-1～A-9 则集中在 §3（数值）与 §5（实验归因）。完整证据链、每篇的 DFT 设置与未核实清单见 [thermal-conductance-literature-2026-11-23.md](thermal-conductance-literature-2026-11-23.md)。

---

## 附录 B：可引用的物理地标与基准数据

**B-1 能量标尺**

| 量 | 值 | 出处 |
|---|---|---|
| Au 声子带顶（L 点 LA） | **4.71 THz = 19.5 meV** | Lynn et al., PRB 8, 3493 (1973) |
| Au 德拜温度 / `k_Bθ_D` | 165 K / 14.2 meV（平均截止，**非**带顶） | 同上文献体系 |
| 300 K 的 `k_BT` | 25.85 meV | — |
| 1 THz | 4.135667696 meV | CODATA |
| `VaspToTHz` | 15.633302 THz·√(eV/Å²/amu) | phonopy `units.py` |
| C–C 伸缩（校验用） | 824–1064 cm⁻¹（`k=30–50 eV/Å²`, `m=12.011`） | 由 `0.064654·√λ` 反推 |

**B-2 热导基准**

| 量 | 值 | 出处 |
|---|---|---|
| `G_Q(300 K)` | **283.9 pW/K**（`π²k_B²/(3h)·T`，每通道） | Schwab et al., Nature 404, 974 (2000) |
| `G_Q(100 K)` / `G_Q(4 K)` | 94.64 / 3.786 pW/K | 同上 |
| Au–烷二硫醇–Au（C2–C10，室温） | **~20 pW/K，几乎不随长度变化** | Cui et al., Nature 572, 628 (2019) |
| 同上，Cui 的计算值 | 16–21（AB，更拉伸）/ 22–33（AG，拉伸较少）pW/K | 同上 |
| ODT / OPE3（实测） | (37±5) / (23±5) pW/K | Mosso et al., Nano Lett. 19, 7614 (2019) + 学位论文 |
| 同上（DFT 饱和值） | 22 / 19 pW/K | 同上 |
| 电子项：C2 / C4 | 5.7 / 1.1 pW/K（更长可忽略） | Cui 2019 |
| 电子项：OPE3 / ODT | 0.07 / 0.01–0.02 pW/K | Mosso 2019 |
| Au 单原子结 | ~500 pW/K（近量子化，WF 成立） | Cui et al., [Science **355**, 1192 (2017)](https://doi.org/10.1126/science.aam6622) |
| 非谐修正（简谐 NEGF 高估） | **~15%**（T>100 K，烷链，ReaxFF） | Li et al., PRB 104, 245413 (2021)，[arXiv:2202.06490](https://arxiv.org/abs/2202.06490) |

**B-3 室温热输运的主导模式（权威表述）**

- 「热输运由分子与金属电极之间的**声子失配**支配」——Mosso et al. 摘要。
- 「非谐会散射**携带大部分热量的低频声子**」——Li et al. 2021。
- 分子内 100–400 meV 的振动**弹性透射严格为零**（远高于 19.5 meV 的 lead 带顶）；玻色因子压制（300 K、100 meV 时 `e^{-E/k_BT}≈0.02`）只是次要原因。**报告 §0 的结论对，但理由不完整**——主因是能带失配，不是玻色因子。

**B-4 同族体系的更多数值（文献调研补充，供 D9 对账带用）**

| 体系 | `G_ph` (pW/K, 300 K) | 方法 | 出处 |
|---|---|---|---|
| 烷二硫醇 n=2–30 CH₂ | **15–45** | TURBOMOLE DFPT + 自写 NEGF | Klöckner et al., [PRB **94**, 205425 (2016)](https://arxiv.org/abs/1609.03827) |
| 苯二胺 / 2,6-二溴 | **25.24 / 10.46** | 同上（**唯一明写检查过无虚频**） | Klöckner et al., [PRB **96**, 245419 (2017)](https://arxiv.org/abs/1708.06021) |
| C60 单体 / 二聚体 | **25.5–46.3 / 7.3–7.4** | 同上（**唯一给全 NEGF 参数**：声子 `η=1e-5 a.u.`、32×32 横向 k） | Klöckner et al., [PRB **95**, 245404 (2017)](https://arxiv.org/abs/1704.00591) |
| 烷烃 N=1,2,4,8 | **25.4 / 33.4 / 30.3 / 5.6** | SIESTA 有限位移 + 自写 NEGF | Sadeghi et al., [arXiv:1510.05945](https://arxiv.org/abs/1510.05945) |
| 同上（力场对照） | **40–45**（强接触）/ **20–22**（弱接触） | LAMMPS EAM + Morse，RNEMD | Wang et al., [arXiv:2209.06266](https://arxiv.org/abs/2209.06266) |
| C60 二聚体（力场对照） | **41.5**（DFT-NEGF 为 7.0，**差 6×**） | LAMMPS AIREBO，无 Au–S 键 | Li/Wang/Segal, [arXiv:2402.06172](https://arxiv.org/abs/2402.06172) |

**对账带的取法**：取**同一化学族**（烷二硫醇 / ODT / OPE3）的 **DFT 谐波计算值**，交集是 **15–25 pW/K**（Cui 计算 16–21 / 22–33、Mosso 理论 19、Klöckner 15–45 的下半段）。**不把 C60（25–46）与力场（20–45）混进来**——前者是 π 体系且电极耦合更强，后者差 6 倍说明唯象势不可用于定量。这就是 D9 用 **C8 或 BDT** 而不是用 nap/az/ind 对账的原因。

> **口径警告（与 T14 那条一致）**：Cui 的**实验**值 ~20 pW/K 与 Mosso 的**实验**值 ODT 37±5 / OPE3 23±5 **不在**这个带里——它们比各自的 DFT 理论值高 1.7×（Mosso）到持平（Cui）。本仓方法在结构上缺非谐（预期偏低 ~15%），因此 15–25 这条线是**「理论对理论」的验收带**；若某天要报「与实验吻合」，必须先扣掉非谐与接触构型的系统偏差，否则是拿两个不同物理量在比。

> 本仓三个目标分子（萘/薁/茚二硫醇）**没有任何文献 `κ_ph` 值**，无法外部对账；它们的可比性只能靠**内部一致性**（同族、同接触、同协议下的相对比较）。这一点必须写进报告，不能包装成「与文献吻合」。
