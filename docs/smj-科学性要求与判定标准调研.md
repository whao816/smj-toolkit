# SMJ 计算流程的科学性要求与"结果正确可用"判定标准调研

> 范围：AutoSMJ（`smj-toolkit`）单分子结输运计算链。  
> 证据分级：**【项目】** = 本仓库文档/代码中的既有约定；**【文献】** = 已检索并核验的公开文献（含 DOI/arXiv）；**【建议】** = 本次调研提出的、尚未在 SMJ 中实现的检查或阈值。  
> 文献检索源：Zotero 本地库、arXiv、OpenAlex、Crossref、Europe PMC；付费全文未读的条目标注为"仅摘要/题录"。

---

## 0. 结论速览

1. **"结果正确"不是一个布尔值，而是四个层次**：实现正确 → 数值收敛 → 物理自洽 → 外部有效。SMJ 现有 10 条判据主要覆盖前两层和部分第三层；第四层（实验/高阶理论基准）目前只有一条派生性质的理论带，没有实验对照，也没有不确定度合成。
2. **通过全部闸门 ≠ 结果正确**。现有闸门是"必要不充分"条件：它们能排除已知的失败模式（未收敛、几何断裂、电极重叠、明显错误的量级/符号），但不能证明模型本身适用。报告与判定语义（`failed` / `inconclusive` / `passed`）应显式写明这一点。
3. **三类系统性偏差最危险**，且都不在现有闸门覆盖内：
   - **泛函能级对齐**：半局域泛函把分子共振放在离 E_F 过近的位置，系统性高估线性响应电导（PBE 的已知偏差方向）；
   - **电极模型**：有限团簇尺寸、缺少镜像电荷/表面极化、E_F 定义与对齐；
   - **接触几何**：吸附位点、锚定 denticity、tip 形貌、结内弛豫——同一分子不同接触可差数倍。
4. **操作性定义**：一个 SMJ 结果可称为"**在声明适用域内正确可用**"，当且仅当同时满足：
   (a) 全链证据完整、可复现、协议未漂移；
   (b) 所有数值收敛检查落在阈值内（基组/电极尺寸/能量网格/窗口/SCF/几何）；
   (c) 物理自洽检查无冲突（符号、量级、趋势、对称性、极限）；
   (d) 与外部基准的偏差落在声明的不确定度内，且适用域（分子家族、锚定化学、电极、温度、偏压线性区）被明确写出。  
   四条任一不满足，只能判 `failed` 或 `inconclusive`，不得称"正确"。
5. **建议优先补的检查**（详见 §6）：基组收敛、电极尺寸收敛、能量窗口收敛、能级对齐报告、不确定度块、ZT 语义修正（电子项≠总 ZT）、构象统计、线性响应适用域声明。

---

## 1. "正确"的四个层次与判定空间

| 层次 | 问题 | 典型证据 | SMJ 现状 |
|---|---|---|---|
| L1 实现正确 | 代码是否做了它声称的事？ | 单元测试、解析极限（G0、T≤通道数）、独立解析器交叉、输入输出哈希 | 部分：解析器交叉校验、哈希链、几何闸门 |
| L2 数值收敛 | 离散近似是否收敛到方程的解？ | 基组/k 点/团簇尺寸/能量网格/窗口/SCF/几何力 | 部分：网格独立性、SCF、力阈值；缺基组、电极尺寸、窗口 |
| L3 物理自洽 | 模型内的结果是否自洽、是否落在已知物理边界内？ | 符号规则、量级带、同系物趋势、左右对称/互易、极限行为 | 有：符号、量级带、峰判据；缺趋势与对称性检查 |
| L4 外部有效 | 模型对该问题是否足够好？ | 实验基准、GW/高阶理论、方法间散布、误差棒、适用域 | 弱：仅 BDT 理论带，无实验对照、无不确定度 |

**判定语义（沿用 SMJ 现约定，正确且应保留）**：证据缺失 → `inconclusive`，不等于"否证"；`gate=error` 未过 → `failed`；`gate=warn` 未过 → `inconclusive`；空 `checks` 不算通过。  
**必须写进报告的一句话**：`passed` 的含义是"在声明的适用域内、证据齐全、无已知冲突"，**不是**"与实验一致"或"物理上唯一正确"。

---

## 2. 必要条件（按计算阶段）

### 2.1 分子与结的化学定义

| 必要条件 | 为什么必要 | 证据 | SMJ 现状 |
|---|---|---|---|
| 明确质子化/电荷态 | 硫醇（-SH）与硫醇盐（-S⁻）的锚定化学与键长不同；去质子化与否改变能级对齐 | 【项目】仅对去质子硫醇盐标定 2.40–2.49 Å；锚定化学不可混 | 已有（几何闸门按锚定类型分流） |
| 明确自旋态/开壳层 | 开壳层分子需自旋极化；闭壳层假设会给出错误的近 E_F 谱 | 【建议】通用量子化学实践；开壳层/双自由基体系的输运对自旋态敏感 | 未见显式检查【建议】 |
| 构象与异构体 | 结电导对构象/扭转角敏感；单构象结果不能代表分子 | 【文献】[Li 2007](https://doi.org/10.1021/ja0762386)（烷基二硫醇的配位几何与构象自由度）；[Yasini 2026](https://doi.org/10.1002/anie.202522129)（室温下结在构象空间涨落） | `smj ensemble` 存在，但不属于阶段链；生产 SCF 用单一几何【建议：报告构象选择依据或分布】 |
| 几何优化收敛 | 未收敛几何的力/键长会污染接触与透射 | 【项目】BFGS 力阈值 1e-2 eV/Å | 已有 |
| 无虚频（孤立分子） | 虚频 = 鞍点，不是极小 | 【项目】`gaussian_no_imaginary` 仅 Gaussian 引擎 | **缺口：FHI-aims 路径无振动检查**【建议】 |
| 结内弛豫（不是气相几何直接移植） | 分子在电极间的键长/取向与气相不同 | 【文献】FHI-aims 手册 AITRANSS 流程：固定两个背层原子弛豫、调分子–Au 距离再弛豫 | 已有（`aims_tip_relax`） |
| 环境（气相/溶剂/温度） | 溶剂通过介电屏蔽移动能级；0 K 真空计算与室温实验有系统差 | 【文献】[Yasini 2026](https://doi.org/10.1002/anie.202522129)：多数计算在真空 0 K 基态几何上做，实验在室温、有溶剂/气氛，差异可达数量级 | 未建模【建议：至少在适用域中声明"真空、0 K 电子结构、线性响应"】 |

### 2.2 结几何与接触

| 必要条件 | 为什么必要 | 证据 | SMJ 现状 |
|---|---|---|---|
| 两侧电极表面/表面指数明确且不重叠 | 左右表面索引错误会让"结"变成同一块电极 | 【项目】`surface_valid`、`electrodes_disjoint` | 已有（error） |
| 左右表面分居分子两侧 | 防止几何镜像/坐标错误 | 【项目】`surface_side` | 已有（error） |
| 锚定原子–电极距离在化学键范围 | 断裂的接触 = 隧穿结而非化学键结 | 【项目】`contacts`（0.7×共价半径和下限）、`anchor_aus`（2.40–2.49 Å） | 已有（error / warn） |
| 吸附位点与 denticity 明确 | 同一分子 top/bridge/hollow、单齿/双齿接触电导可差数倍，且决定可重复性 | 【文献】[Kaneko 2016](https://doi.org/10.1021/jacs.5b11559)（吸附位点选择决定电导与可重复性）；[Li 2007](https://doi.org/10.1021/ja0762386) | **缺口：未显式记录/收敛吸附位点**【建议】 |
| tip 形貌与电极层数 | 尖端原子配位、层数改变电极 DOS 与耦合；SMJ 的团簇模型更接近 STM 针尖–表面构型，不等同于机械断裂结的电极几何，跨实验比较需说明 | 【文献】[Strange 2008](https://doi.org/10.1063/1.2839275)（超胞尺寸与 k∥ 收敛）；【项目】README §9：Au60 与 20 Au/侧不可直接比较 | **缺口：无电极尺寸/层数收敛检查**；模型–实验构型差异未声明【建议】 |
| 结几何的拉伸/压缩自由度 | 实验是动态拉伸过程，理论取单一几何；需说明取的是哪种构型 | 【文献】[Yasini 2026](https://doi.org/10.1002/anie.202522129) | 未声明【建议】 |

### 2.3 电极模型与能级对齐

| 必要条件 | 为什么必要 | 证据 | SMJ 现状 |
|---|---|---|---|
| E_F 的定义与两侧对齐 | 透射、S、ZT 都在 E_F 处取值；E_F 偏移 0.1 eV 可显著改变近共振结果 | 【项目】`surface_fingerprint` 校验、TE 数据 | 部分（有指纹校验，无 E_F 敏感性检查） |
| 镜像电荷/表面极化 | 金属表面对分子能级的极化重整化可达 eV 量级；PBE 不含 | 【文献】[Neaton 2006](https://doi.org/10.1103/physrevlett.97.216405)（苯/石墨 gap 10.5→7.2 eV）；[García-Lastra 2009](https://arxiv.org/abs/0910.5304)；[Perrin 2013](https://arxiv.org/abs/1403.4266)（实验可调镜像电荷）；[Tamblyn 2011](https://doi.org/10.1103/physrevb.84.201402)（GW 与光电子能谱对齐） | **缺口：无极化修正，也无标志位**【建议】 |
| 电极 DOS 与 tip 原子角色 | 透射峰高与宽度受电极 DOS、耦合原子数控制 | 【文献】[Gulyaev 2026](https://doi.org/10.1021/acs.nanolett.6c01462)（Breit-Wigner + 轨道能 + 金属 DOS 修正） | 未报告【建议】 |
| 团簇 vs 周期模型的选择与边界 | 有限团簇无 k∥，但边缘态/尺寸效应需收敛；周期模型需 k∥ 收敛 | 【文献】[Strange 2008](https://doi.org/10.1063/1.2839275) | 团簇模型；尺寸收敛未做【建议】 |
| 相对论（Au） | scalar ZORA 影响 Au 5d/6s 位置 | 【项目】`relativistic atomic_zora scalar` | 已有 |

### 2.4 电子结构计算

| 必要条件 | 为什么必要 | 证据 | SMJ 现状 |
|---|---|---|---|
| 泛函的已知偏差被承认并记录 | GGA 过度离域、SIE、共振位置偏差；不同泛函可差数倍 | 【文献】[Toher & Sanvito 2008](https://arxiv.org/abs/0712.1747)（LDA/GGA 过度离域、人为增大展宽；SIC 使共振远离 E_F）；[Koentopp 2008](https://arxiv.org/cond-mat/0703591)（弱耦合时 KS 能级错位）；[Gulyaev 2026](https://doi.org/10.1021/acs.nanolett.6c01462)（PBE 把共振放得离 E_F 太近→高估电导） | 固定 PBE；**无泛函敏感性检查**【建议】 |
| 基组收敛 | 基组不足会系统性移动透射；定量需 DZP 级别 | 【文献】[Strange 2008](https://doi.org/10.1063/1.2839275)：DZP 才足以定量，SIESTA 透射系统性下移；【项目】默认 `light` | **缺口：light/tight 收敛未做**【建议】 |
| SCF 收敛且残差一致 | "收敛"标志不等于结果稳定；近简并/近共振时残差敏感 | 【项目】`sc_accuracy_rho/eev/etot = 1E-5/1E-3/1E-6`、`sc_iter_limit 500` | 已有（`scf_converged`，engineering/error） |
| 占据/展宽设置 | `gaussian 0.01` 的展宽影响 E_F 与 E_F 附近 T(E) | 【项目】`occupation_type gaussian 0.01` | 已有；**无展宽敏感性检查**【建议】 |
| 真空/偶极修正（若用 slab） | 周期性 slab 的偶极相互作用移动能级 | 通用 DFT 实践 | 团簇模型下不适用；若未来切周期需补【建议】 |
| 定量精度需超越 DFT | GW 对烷烃 benchmark 可定量；DFT 不行 | 【文献】[Strange & Thygesen 2011](https://arxiv.org/abs/1108.3687)：GW 对 Au/烷烃二胺 benchmark 与实验定量一致；[FHI-aims Roadmap 2025](https://arxiv.org/abs/2505.00125) 承认半局域泛函 KS 谱偏差，计划 scissor/hybrid + image-charge 修正 | 不适用/不计划【建议：在适用域中声明"趋势与量级，不是定量实验预测"】 |

### 2.5 输运计算（NEGF / AITRANSS）

| 必要条件 | 为什么必要 | 证据 | SMJ 现状 |
|---|---|---|---|
| 自能模型参数固定且被视为"模型常数" | σ 决定耦合展宽；换参数结果不可比 | 【项目】σ_near 0.05 Ha / σ_far 0.10 Ha，明确写"不是收敛参数"，生产必须原样使用 | 已有（且文档已声明） |
| 能量网格收敛 | T(E)、dT/dE、S 都对网格敏感 | 【项目】0.2→0.05 meV 测试：L0/S/k_e/ZT_e 相对差 <1e-8；T(E_F) 差异仅来自最近格点 | 已有（`grid_independence`，warn） |
| 能量窗口收敛（积分截断） | ±2 eV 窗口若截断积分尾部，L0/L1/L2 不稳 | 【建议】重算 ±1.5/±2.5 eV 窗口，要求 <1% | **缺口**【建议】 |
| 温度与费米展宽 | 300 K 展宽影响 E_F 邻域；S 是导数型量，更敏感 | 【项目】T=300 K | 固定值；**无温度敏感性检查**【建议】 |
| 线性响应适用域 | AITRANSS 是零偏压后处理；G、S 是 V→0 极限；有限偏压非线性不覆盖 | 【文献】[Arnold 2007](https://doi.org/10.1063/1.2716664)（形式框架）；[Brandbyge 2002](https://doi.org/10.1103/physrevb.65.165401)；[Papior 2017](https://arxiv.org/abs/1607.04464) | 未声明【建议：写入 run_spec/报告】 |
| 左右互易/对称性 | 对称结应有 T_LR(E)=T_RL(E)；不对称需解释 | 【建议】交换左右电极重算并比较 | **缺口**【建议】 |
| 电流/粒子数守恒（有限偏压时） | 非自洽 NEGF 在有限偏压下可能不守恒 | 通用 NEGF 实践；本项目零偏压不涉及 | 不适用（零偏压）【建议：若扩展有限偏压需补】 |

### 2.6 后处理与物理量

| 必要条件 | 为什么必要 | 证据 | SMJ 现状 |
|---|---|---|---|
| Landauer 积分与常数正确 | 单位/常数错会让所有量整体偏移 | 【项目】k_B=8.6173324e-5 eV/K、h=4.135667662e-15 eV·s | 已有 |
| S 的符号约定与机制一致 | HOMO 传导 S>0、LUMO 传导 S<0；符号错=机制判断错 | 【文献】[Widawsky 2012](https://arxiv.org/abs/1201.1837)：amine-Au（HOMO）S>0，pyridine-Au（LUMO）S<0；同时测 G 与 S | 已有（BDT 专属 `seebeck_sign`，error） |
| ZT 的语义 | `ZT_electronic` 只是电子项；总 ZT 需 κ_ph，电子项是上限 | 【项目】README §9 已声明；【文献】[Klöckner & Pauly 2019](https://arxiv.org/abs/1910.02443)（热导的接触/取向/构型变异性） | 已声明；**`k_e_zt_order` 未实现** |
| 峰判据的物理定义与敏感性 | "E_F±1 eV 内有峰"依赖展宽/网格/定义；近共振可能是真实物理而非错误 | 【项目】`t_ef_resonance_attribution`（warn，按族归因 0.15–1.0 eV） | 有；**峰宽/来源与敏感性分析仍缺**【建议】 |
| 量级带与趋势 | 单点落带不足以证明正确；同系物趋势更稳健 | 【建议】至少 2–3 个同系物/锚定基团的一致趋势 | **缺口**【建议】 |
| 解析极限/边界 | T(E)≥0、T≤通道数、G<G0（单分子隧穿典型远小于 G0） | 通用；【项目】BDT 带 0.03–0.40 G0 | 部分（带检查隐含）【建议显式】 |

### 2.7 溯源、可复现与报告

| 必要条件 | 为什么必要 | 证据 | SMJ 现状 |
|---|---|---|---|
| 输入/输出哈希与父阶段指纹 | 防止上游被静默替换；可定位偏差来源 | 【项目】`run_spec.json`、每阶段 `inputs_sha256/outputs_sha256/parent_fingerprint` | 已有 |
| 协议哈希与漂移拒绝 | 相当于"计算预注册"：判据与参数在跑之前固定，防止事后调参 | 【项目】`protocol_sha256` 只覆盖会改物理结果的子集；漂移默认拒绝，`--force` 留档 | 已有（设计良好） |
| 软件/版本/二进制哈希 | 跨机器复现的前提 | 【项目】SMJ/Python 版本、git commit/dirty、二进制 sha256 | 已有 |
| 失败与超时的语义 | 超时不算失败，但也不能算通过 | 【项目】README §9 | 已有 |
| 不确定度报告 | 单点值无法表达数值/模型误差 | 【建议】见 §3.3 | **缺口** |
| 适用域声明 | 防止把"趋势正确"当成"定量预测" | 【文献】[Evers 2020](https://arxiv.org/abs/1906.10449)：理论预测趋势好、定量一致性差，需定义良好的 benchmark | 部分（判据族 fail-safe）【建议：显式字段】 |

---

## 3. 认定"正确可用"的量化标准

### 3.1 L2 数值收敛（建议阈值）

| 检查 | 建议阈值 | 依据 |
|---|---|---|
| 几何力 | < 1e-2 eV/Å（项目既有） | 【项目】 |
| SCF | 电子密度/总能/本征值残差达到设定值且不贴限 | 【项目】1E-5/1E-3/1E-6 |
| 能量网格 | L0/L1/L2 相对变化 < 1%（项目既有）；另单独检查 S | 【项目】0.2→0.05 meV 实测 <1e-8；【建议】S 是导数，需单独看 |
| 能量窗口 | ±1.5/±2.0/±2.5 eV 下 L0/L1/L2 变化 < 1% | 【建议】 |
| 基组 | `light`→`tight`（或更高）G 变化 < 5%，否则报告差异并降级结论 | 【文献】[Strange 2008](https://doi.org/10.1063/1.2839275)（DZP 才定量）【建议】 |
| 电极团簇尺寸 | +1 层或更大团簇（如 37→61 atoms）G 变化 < 10% | 【文献】[Strange 2008](https://doi.org/10.1063/1.2839275)（超胞收敛）【建议】 |
| k 点（若周期） | ΔG < 1–2% | 【文献】[Strange 2008](https://doi.org/10.1063/1.2839275)【建议】 |
| 展宽/占据 | 0.005/0.01/0.02 Ha 下 G、S 变化报告 | 【建议】 |
| 温度 | 100/200/300 K 下 S 变化报告（实验低温下 S 近温度无关） | 【文献】[Widawsky 2012](https://arxiv.org/abs/1201.1837)【建议】 |
| E_F 敏感性 | E_F±0.05 eV 内移动时 G、S 的变化 | 【建议】 |

### 3.2 L3 物理自洽（建议检查）

| 检查 | 判据 | 依据 |
|---|---|---|
| 符号规则 | HOMO 传导 S>0；LUMO 传导 S<0；与锚定基团/传输机制一致 | 【文献】[Widawsky 2012](https://arxiv.org/abs/1201.1837) |
| 量级带 | 按分子家族给出理论带（现有 BDT 0.03–0.40 G0），并注明"理论带，非实验带" | 【项目】；【文献】[Evers 2020](https://arxiv.org/abs/1906.10449)（定量一致性差） |
| 同系物趋势 | 至少 2–3 个同系物/锚定基团的 G、S 单调性或已知趋势一致 | 【建议】 |
| 对称/互易 | 对称结 T_LR=T_RL；不对称结差异需可由几何解释 | 【建议】 |
| 极限行为 | T(E)≥0；T≤通道数；G<G0（单分子隧穿）；远离共振时指数衰减 | 通用 |
| 近共振判据 | E_F±1 eV 内有峰 → 只作 warn，不作否证；按族预期通道归因，需报告峰宽/来源（分子轨道 vs 电极态） | 【项目】`t_ef_resonance_attribution`；【文献】[Gulyaev 2026](https://doi.org/10.1021/acs.nanolett.6c01462)（PBE 近共振是已知偏差来源） |
| 独立实现交叉 | 与另一代码/另一基组/另一泛函或 Breit-Wigner 解析拟合对照 | 【文献】[Gulyaev 2026](https://doi.org/10.1021/acs.nanolett.6c01462)（BW 拟合 + 轨道能 + DOS 修正）【建议】 |

### 3.3 L4 外部有效（建议报告格式）

- **基准表**：对每个 benchmark 体系（烷烃二硫醇/二胺、BDT 等）给出：本文值、文献理论值范围、实验值范围、偏差倍数、方法差异来源。
- **不确定度块**：把数值误差（网格/基组/电极尺寸/窗口）与模型误差（泛函、σ 参数、极化缺失）分开列；给出"合成区间"而非单点。
- **适用域声明**（建议字段）：分子家族、锚定化学、电极材料/表面、tip 形貌、温度、偏压（线性响应）、环境（真空/溶剂）、构象处理。
- **不能宣称的结论**：定量实验预测、有限偏压非线性、总 ZT、溶剂/温度效应、开壳层自旋效应（除非已建模）。

---

## 4. SMJ 现有判据的覆盖度评估（2026-11-15 修订：10 条）

> 本节在调研报告交付后随判据口径修订同步更新：删除 `t_ef_magnitude`，
> `t_ef_no_peak_near_ef` 改为 `t_ef_resonance_attribution`，`g_in_theory_band` 降为
> `warn`；证据核验见 [`docs/criteria-source-research.md`](criteria-source-research.md)。

| 判据 | 层次 | 评价 |
|---|---|---|
| `gaussian_no_imaginary` | L2 | 仅 Gaussian 路径；FHI-aims 缺振动检查 |
| `scf_converged` | L2 | 必要；但"收敛标志"不等于残差稳定 |
| `te_present` | L1 | 数据完整性 |
| `grid_independence` | L2 | 已做且实测极稳；建议补 S 与窗口 |
| `electrodes_disjoint` | L1 | 必要 |
| `au_s_bond` | L3 | 带收窄为 2.40–2.44 Å（Jackson 2009 摘要 S–Au=2.42±0.02 Å）；仅 thiolate；标定来源已写入 note |
| `t_ef_resonance_attribution` | L3 | 由"无峰"改为"按族归因"：BDT=HOMO、BPY=LUMO、AZ=HOMO 侧、烷烃=无；2CN/SME 按 Seebeck 自洽回退；0.15 eV 下限是 PBE 能级误差量级的工程容差 |
| `g_in_theory_band` | L3 | 0.03–0.40 G0 是**派生包络**（Thygesen §V.2 的 0.1 G0 与 ×1/3–×4），gate 降为 warn、evidence=derived；BDT 专属 |
| `seebeck_sign` | L3 | 有文献支撑（Widawsky）；BDT 专属 |
| `k_e_zt_order` | L3/L4 | 未实现；ZT_e 语义需修正 |

**结构性缺口**：基组收敛、电极尺寸收敛、窗口收敛、能级对齐、极化/镜像电荷、泛函敏感性、不确定度、构象统计、实验基准、互易/对称性、线性响应适用域。

**原设计问题与已落地修订（2026-11-15）**：`docs/validation-0.2.0.md` 中三体系物理判定原先均为 `inconclusive`，触发项都是 `t_ef_no_peak_near_ef`（E_F±1 eV 内有透射峰）；该判据对全链适用，而 2CN（硫醚锚定、腈基取代，LUMO 传导）与 AZ 的近 E_F 共振可能是**真实机制**。已按本报告建议改为 `t_ef_resonance_attribution`：无峰 `pass`；峰落在族预期通道侧且 0.15≤|ΔE|≤1.0 eV → 「物理信号已归因」；侧别不符/过近/通道未知 → `warn`（需解释）；每条族规则均带 `source{doi,locator}`。0.2.0 三体系重评均 `passed`。仍未覆盖：展宽/网格敏感性、峰宽与来源（分子轨道 vs 电极态）的报告字段。

---

## 5. 建议的判定流程（决策树）

```
输入完整性/哈希链 → 失败：failed
        ↓ 通过
数值收敛（SCF/力/网格/窗口/基组/电极尺寸）→ 任一未过：inconclusive
        ↓ 通过
物理自洽（符号/量级/趋势/对称/极限）→ 冲突：inconclusive（或 failed，若属硬错误）
        ↓ 通过
外部基准与适用域（benchmark 表 + 不确定度 + 域声明）
        ↓
passed（含义：域内证据齐全、无已知冲突；不等于实验正确）
```

**输出三态之外，建议增加一个字段** `confidence`（或 `evidence_level`）：`numeric_only` / `physically_checked` / `benchmarked`，让读者一眼看出"通过"到了哪一层。

---

## 6. 落地建议（按优先级）

**P0（文档与语义，零代码）**
1. 在 `README` 与 `docs/smj-计算流程梳理.md` §6 增加一段"判定语义与适用域"：明确 `passed` ≠ 实验正确；明确零偏压线性响应、真空/0 K 电子结构、单构象、团簇电极、σ 为模型常数、PBE/light。
2. `run_spec.json` 增加 `applicability` 块（分子家族/锚定/电极/温度/偏压/环境/构象），随报告输出。
3. 修正 ZT 表述：`ZT_electronic` 为电子项上限；`k_e_zt_order` 要么实现 κ_ph 估算，要么降级为 `info` 并写明"仅电子项"。

**P1（新增检查，可先 warn）**
4. `basis_set_convergence`（physical, warn）：light→tight（或更高）重算 G；变化 >5% 记 `inconclusive` 并附差异。来源：[Strange 2008](https://doi.org/10.1063/1.2839275)。
5. `electrode_size_convergence`（physical, warn）：37→61 atoms 或 +1 层；变化 >10% 记 `inconclusive`。来源同上 + 项目 README §9。
6. `window_sensitivity`（physical, warn）：±1.5/±2.0/±2.5 eV 窗口下 L0/L1/L2 变化 <1%。
7. `level_alignment`（physical, info/warn）：报告 ε_HOMO/ε_LUMO（或透射峰位）相对 E_F 的位置；若 |E_res−E_F| < 0.3 eV 且结论依赖线性响应，标 warn。来源：[Gulyaev 2026](https://doi.org/10.1021/acs.nanolett.6c01462)、[Neaton 2006](https://doi.org/10.1103/physrevlett.97.216405)、[Tamblyn 2011](https://doi.org/10.1103/physrevb.84.201402)。
8. `uncertainty_block`（info）：数值/模型误差分列，输出合成区间。

**P2（需要更多工作）**
9. 构象集成：`ensemble` 结果进入报告，给 Boltzmann 权重下的 G/S 分布（中位数 + IQR），而非单构象点值。
10. 互易/对称性检查：交换左右电极重算，比较 T(E)。
11. 实验 benchmark 表：烷烃二硫醇/二胺 + BDT；给出理论/实验比值分布。来源：[Evers 2020](https://arxiv.org/abs/1906.10449)、[Strange & Thygesen 2011](https://arxiv.org/abs/1108.3687)。
12. tip 形貌/吸附位点枚举：至少报告候选位点及其 G 散布。来源：[Kaneko 2016](https://doi.org/10.1021/jacs.5b11559)。
13. 正/负对照：同一流程跑一对已知对照（烷烃类隧穿低电导 vs BDT 近共振），确认流程能复现"低–高"趋势，而不仅是单点落带。

---

## 7. 关键文献证据摘录（与上述论断对应）

- **[Evers 2020]** 单分子电子输运的理论预测在趋势上可靠，但定量一致性仍是核心挑战；需要实验定义良好的 benchmark 体系；现有模型可能缺失动态镜像电荷、电子关联、极化子等机制。
- **[Koentopp 2008]** DFT-NEGF 在弱耦合下 KS 能级错位；标准方法漏掉 XC 对电流的修正——给出"何时可信、何时应失效"的边界。
- **[Toher & Sanvito 2008]** LDA/GGA 过度离域分子能级、人为增大展宽；SIC 使共振远离 E_F、电流下降。
- **[Neaton 2006]** 金属表面镜像电荷使分子 gap 显著重整化（苯/石墨 10.5→7.2 eV）。
- **[García-Lastra 2009]** 金属/半导体表面极化系统性重整化分子能级。
- **[Perrin 2013]** 实验观测到可调镜像电荷效应，量级与理论一致。
- **[Tamblyn 2011]** GW 能级对齐与光电子能谱一致；常规 DFT 对齐有系统偏差。
- **[Strange 2008]** 五个结的 benchmark：平面波+MLWF vs SIESTA；已对超胞尺寸与 k∥ 收敛；DZP 基组才足以定量；SIESTA 透射系统性下移。
- **[Strange & Thygesen 2011]** 对 Au/烷烃二胺 benchmark，GW 与实验定量一致——定量精度需要超越 DFT。
- **[Gulyaev 2026]** PBE 把共振放得离 E_F 太近→线性响应电导高估；DFT+Σ 同时修正 SIE 与表面极化；BW 拟合 + 轨道能 + 金属 DOS 修正可逼近高阶泛函与实验；修正"largely systematic"。
- **[Widawsky 2012]** 同时测 G 与 S：amine-Au（HOMO）S>0，pyridine-Au（LUMO）S<0；自能修正 DFT 与实验定量一致。
- **[Yasini 2026]** 综述：计算多在真空 0 K 基态几何，实验在室温、有溶剂/气氛、构象涨落；理论-实验可差数量级；溶剂通过介电屏蔽移动能级；吸附位点/denticity/取向强烈影响电导。
- **[Kaneko 2016]** 吸附位点选择决定电导与可重复性。
- **[Li 2007]** Au|烷基二硫醇|Au 的配位几何与构象自由度导致多组电导峰。
- **[Klöckner & Pauly 2019]** Au–烷烃–Au 热导的接触/取向/构型变异性——ZT 需 κ_ph，电子项不是全部。
- **[FHI-aims Roadmap 2025]** 承认半局域泛函 KS 谱偏差，计划 scissor/hybrid + image-charge 修正。
- **[Kim 2008]** 题名即针对 NEGF 器件输运的数值精度（烷烃单分子结）；仅题录，未读全文。
- **[Brandbyge 2002] / [Papior 2017] / [Arnold 2007]** DFT-NEGF 与 AITRANSS 的方法学基础。

---

## 8. 参考文献

[1] EVERS F, KORYTÁR R, TEWARI S, 等. Advances and challenges in single-molecule electron transport[J/OL]. arXiv:1906.10449, 2019. https://arxiv.org/abs/1906.10449. （亦见 Rev. Mod. Phys. 92, 035001 (2020)）

[2] KOENTOPP M, CHANG C, BURKE K, 等. Density functional calculations of nanoscale conductance[J/OL]. arXiv:cond-mat/0703591, 2007. https://arxiv.org/abs/cond-mat/0703591. （亦见 J. Phys.: Condens. Matter 20, 083203 (2008)）

[3] STRANGE M, KRISTENSEN I S, THYGESEN K S, 等. Benchmark density functional theory calculations for nanoscale conductance[J]. The Journal of Chemical Physics, 2008, 128(11): 114714. DOI:10.1063/1.2839275.

[4] TOHER C, SANVITO S. Effects of self-interaction corrections on the transport properties of phenyl-based molecular junctions[J/OL]. arXiv:0712.1747, 2007. https://arxiv.org/abs/0712.1747. （亦见 Phys. Rev. B 77, 155402 (2008)）

[5] NEATON J B, HYBERTSEN M S, LOUIE S G. Renormalization of molecular electronic levels at metal-molecule interfaces[J]. Physical Review Letters, 2006, 97(21): 216405. DOI:10.1103/PhysRevLett.97.216405.

[6] GARCIA-LASTRA J M, ROSTGAARD C, RUBIO A, 等. Polarization-induced renormalization of molecular levels at metallic and semiconducting surfaces[J/OL]. arXiv:0910.5304, 2009. https://arxiv.org/abs/0910.5304. （亦见 Phys. Rev. B 80, 245427 (2009)）

[7] PERRIN M L, VERZIJL C J O, MARTIN C A, 等. Large tunable image-charge effects in single-molecule junctions[J/OL]. arXiv:1403.4266, 2014. https://arxiv.org/abs/1403.4266. （亦见 Nat. Nanotechnol. 8, 282–287 (2013)）

[8] TAMBLYN I, DARANCET P, QUEK S Y, 等. Electronic energy level alignment at metal-molecule interfaces with a GW approach[J]. Physical Review B, 2011, 84(20): 201402. DOI:10.1103/PhysRevB.84.201402.

[9] GULYAEV A, HAZARIKA J, LIU Z-F, 等. A computationally efficient and accurate method for predicting conductance of single-molecule junctions[J]. Nano Letters, 2026, 26(22): 7429–7434. DOI:10.1021/acs.nanolett.6c01462.

[10] STRANGE M, THYGESEN K S. Towards quantitative accuracy in first-principles transport calculations: the GW method applied to alkane/gold junctions[J/OL]. arXiv:1108.3687, 2011. https://arxiv.org/abs/1108.3687.

[11] WIDAWSKY J R, DARANCET P, NEATON J B, 等. Simultaneous determination of conductance and thermopower of single molecule junctions[J/OL]. arXiv:1201.1837, 2012. https://arxiv.org/abs/1201.1837. （亦见 Nano Lett. 2012, DOI:10.1021/nl203634m）

[12] YASINI P, BATZINGER K, SMEU M, 等. Fundamentals, measurement and regulation of the conductance of single molecule junctions[J]. Angewandte Chemie International Edition, 2026, 65(35): e202522129. DOI:10.1002/anie.202522129.

[13] KANEKO S, MURAI D, MARQUÉS-GONZÁLEZ S, 等. Site-selection in single-molecule junction for highly reproducible molecular electronics[J]. Journal of the American Chemical Society, 2016, 138(4): 1294–1300. DOI:10.1021/jacs.5b11559.

[14] LI C, POBELOV I, WANDLOWSKI T, 等. Charge transport in single Au|alkanedithiol|Au junctions: coordination geometries and conformational degrees of freedom[J]. Journal of the American Chemical Society, 2007, 130(1): 318–326. DOI:10.1021/ja0762386.

[15] KLÖCKNER J C, PAULY F. Variability of the thermal conductance of gold-alkane-gold single-molecule junctions studied using ab-initio and molecular dynamics approaches[J/OL]. arXiv:1910.02443, 2019. https://arxiv.org/abs/1910.02443. DOI:10.48550/arXiv.1910.02443.

[16] KIM Y-H. Toward numerically accurate first-principles calculations of nano-device charge transport characteristics: the case of alkane single-molecule junctions[J]. Journal of the Korean Physical Society, 2008, 52: 1181–1186. DOI:10.3938/jkps.52.1181.

[17] BRANDBYGE M, MOZOS J-L, ORDEJÓN P, 等. Density-functional method for nonequilibrium electron transport[J]. Physical Review B, 2002, 65(16): 165401. DOI:10.1103/PhysRevB.65.165401.

[18] PAPIOR N, LORENTE N, FREDERIKSEN T, 等. Improvements on non-equilibrium and transport Green function techniques: the next-generation transiesta[J/OL]. arXiv:1607.04464, 2016. https://arxiv.org/abs/1607.04464. （亦见 Comput. Phys. Commun. 212, 8–24 (2017)）

[19] ARNOLD A, WEIGEND F, EVERS F. Quantum chemistry calculations for molecules coupled to reservoirs: formalism, implementation, and application to benzenedithiol[J]. The Journal of Chemical Physics, 2007, 126(17): 174101. DOI:10.1063/1.2716664.

[20] ABBOTT J W, ACOSTA C M, AKKOUSH A, 等. Roadmap on advancements of the FHI-aims software package[J/OL]. arXiv:2505.00125, 2025. https://arxiv.org/abs/2505.00125. （亦见 Electronic Structure, 2026, DOI:10.1088/2516-1075/ae8067）

**补充方法学/软件文档（非 DOI 条目）**

[21] FHI-aims 官方手册（AITRANSS 章节）：extended molecule 构建、固定背层弛豫、分子–电极距离调整。https://fhi-aims.org

---

## 9. 结语

SMJ 的判定体系在"工程可核验 + fail-closed + 协议哈希"上是扎实的，这部分已经超出了多数同类脚本的水平；真正需要补的是**数值收敛的完整证据链**（基组、电极尺寸、窗口）、**模型偏差的显式承认**（泛函近共振、镜像电荷缺失）、以及**不确定度与适用域的报告**。把这三件事补上，"结果正确可用"才能从一句判断变成一个可审计的结论。
