# 单分子结热导计算：已发表工作的技术栈与软件栈调研（2026-11-23）

配套文档：`docs/thermal-conductance-plan-2026-11-23.md`（功能开发方案）。
本文只回答一个问题：**别人算单分子结热导，实际用什么软件、什么流程**。物理结论不在此展开。

标注约定：**【明】** = 论文/SI/手册原文逐字读到；**【推】** = 从作者谱系、致谢或相邻文献推断；**【缺】** = 检索范围内未找到（不等于不存在）。

---

## §0 一句话结论

这个领域**没有"装好就能算"的软件**。能算单分子结 κ_ph 的，实际是 4–5 个理论组的自研码；DFT 软件只负责产力常数，输运一律自己解。三条主线互不相通：

| 主线 | DFT/力常数 | 输运 | 代表 |
|---|---|---|---|
| **A 量子化学团簇 + 解析 Hessian** | TURBOMOLE（DFPT） | 自研 NEGF + decimation | Pauly 组（Cui 2019 Nature） |
| **B 局域基 DFT + 有限位移** | SIESTA（±0.01 Å 中心差分） | 自研 NEGF（电子走 GOLLUM） | Sadeghi–Lambert 组（Mosso 2019） |
| **C 经典 MD** | LAMMPS + EAM/Morse/ReaxFF/QM-FF | NEMD 家族 | Segal 组、Cui–Prampolini（Yelishala 2025） |

**没有一篇已发表工作点名 QuantumATK `PhononTransmissionSpectrum` 算单分子结 κ_ph；没有一篇用 VASP+phonopy 做分子结弹道声子透射；GOLLUM 2014 原文里根本没有声子模块。**

---

## §1 主线 A：量子化学团簇 + DFPT + 自研 NEGF（Pauly 组）

### 1.1 技术栈（三篇 PRB + 一篇 Nature 逐字相同）

| 环节 | 具体 |
|---|---|
| DFT 引擎 | **TURBOMOLE**（2016 篇 6.5，2019 Nature 篇 **7.1**）【明】 |
| 泛函/基组 | **PBE / def2-SV(P)** + 对应 Coulomb fitting basis（RI-DFT）；全电子、无赝势、**无 k 点**（体系是团簇）【明】 |
| 力常数 | **DFPT 解析二阶导数**（Deglmann/Furche/Ahlrichs 2002/2004 实现）——**不是有限位移，没有位移幅度这个参数**【明】 |
| 收敛判据 | 总能量 1e-9 a.u.；几何优化至笛卡尔梯度最大范数 < 1e-5 a.u.【明】 |
| 电极 | **半无限完美晶体**；表面格林函数由 **decimation** 得到【明】 |
| 电极力常数 | **另算一次**：先用"a big cluster of several hundred atoms"算 bulk FC，Fourier 变换 + 对称化后 decimation【明】 |
| 输运码 | **自研 NEGF**（Pauly NJP 10, 125019 (2008) 的 cluster-based 框架；Bürkle PRB 91, 165419 (2015)）【明】 |
| 公式 | `G_ph = (1/h)∫dE·E·τ_ph(E)·∂n/∂T`，`τ = Tr[G^r Γ_L G^a Γ_R]`【明】 |

### 1.2 具体论文

| 论文 | 体系 | 关键数字 | 备注 |
|---|---|---|---|
| Klöckner, Bürkle, Cuevas, Pauly, [PRB 94, 205425 (2016)](https://arxiv.org/abs/1609.03827) | Au–烷二硫醇–Au，n=2–30 个 CH₂（~4.5 nm） | 15–45 pW/K @300 K | 长度依赖主线 |
| Klöckner, Siebler, Cuevas, Pauly, [PRB 95, 245404 (2017)](https://arxiv.org/abs/1704.00591) | C60 单体/二聚体 | κ_ph 25.5/46.3/20.0（单体，三种尖端组合） | **唯一给出完整 NEGF 数值参数的一篇**：电子 η=1e-2 a.u.、声子 **η=1e-5 a.u.**、均 **32×32 横向 k 点**；额外加 Grimme 色散；还算了光子通道（fluctuational electrodynamics + proximity approximation） |
| Klöckner, Cuevas, Pauly, [PRB 96, 245419 (2017)](https://arxiv.org/abs/1708.06021) | Au–苯二胺–Au 及卤素取代，**氨基（–NH₂）锚定** | **Table 1** 的 `κ_pn` 列：1,4-二氨基苯(para) **25.24**；2-F 24.40 / 2-Cl 22.17 / 2-Br 16.94 / 2-I 15.24 / 2,5-二Br 17.85 / 2,6-二Br **10.46** pW/K；同行 `κ_el`=2.61 pW/K、`G`=4.7×10⁻³ G₀（理论）/6.4×10⁻³ G₀（实验 [70]） | **唯一写了 HPC 资源**（bwUniCluster / JUSTUS）；**唯一明写"检查过无虚频"**；干涉分析用 η=1e-3 meV 的谱表示（展示用，非 κ 计算值）；文中「Au 的德拜能约 20 meV」可作带顶佐证。**锚定是 –NH₂，不能直接对标本仓的 –S⁻ 结** |
| Cui et al., [Nature 572, 628 (2019)](https://doi.org/10.1038/s41586-019-1420-z) | Au–烷二硫醇–Au，C2–C10（**硫醇锚定**） | **实验（Fig. 3c，C-SThM，295 K）**：「nearly length-independent … approximately **20 pW K⁻¹**」，每分子 108–155 条痕（C2 155 / C4 133 / C8 110 / C10 108）；**同文理论（Landauer–Büttiker）**：AG **22–33**、AB **16–21** pW/K，且「nearly independent of molecular length」 | 理论由 **J. C. Klöckner 在 F. Pauly 指导下**完成；Code availability 原文："The corresponding **custom-developed code for the description of phonon transport** implements the procedures outlined in ref. 31 and is available from F.P. on reasonable request."。**实验与理论是两个量，不要混用** |
| Blaschke & Pauly, [ACS Nano 19, 32093 (2025)](https://arxiv.org/abs/2505.19158) | 遗传算法筛分子 | 0.07–33 pW/K（600 K） | **已不是 DFT**：弛豫与 Hessian 用 **xTB/GFN1**；电极退化成 **Debye 解析表面 GF**；每侧只留 1 个 Au 原子，γ(Au–S)=−7.0 eV/Å²；代码在 Zenodo |

### 1.3 这一支的结构特征

- **没有 6N 问题**。DFPT 一次拿到全 Hessian，所以"位移幅度多大""超胞多大"这两个参数在 A 支根本不存在。取而代之的收敛旋钮是 SCF 能量、几何梯度、以及"有没有虚频"。
- **电极与器件区的化学环境不一致是显式接受的近似**：器件区 FC 来自真实团簇，电极 FC 来自另一次 bulk 计算。
- 复现不了的细节九成落在 **η 与横向 k 点**（只有 PRB 95 给了）。
- 2025 年开始自我修正：用 xTB 换 DFT（为了 10²⁷ 的化学空间），说明"简谐 NEGF-DFT"正在退化成**筛查器**。

---

## §2 主线 B：SIESTA + 有限位移 + 自研 NEGF（Sadeghi–Lambert 组）

### 2.1 技术栈（Mosso 2019 SI 逐字）

| 环节 | 具体 |
|---|---|
| DFT 引擎 | **SIESTA**【明】 |
| 泛函/基组 | **GGA/PBE + DZP**，实空间网格等效 **250 Ry**【明】（同组 2021 篇用 LDA/350 Ry，见 §2.2） |
| 几何收敛 | 力容差 **10 meV/Å**【明】 |
| 力常数 | **有限差分，δ = 0.01 Å**，沿 ±x/±y/±z 六个方向；`K_ij^qq' = [F_i^q(δq'_j) − F_i^q(−δq'_j)] / 2δq'_j`【明】 |
| 对角项 | **用动量守恒补齐**：`K_ii = −Σ_{j≠i} K_ij`【明】——省掉 3N 次自项位移，很实用 |
| 电极 | 半无限【明】 |
| 电子输运 | **GOLLUM**【明："the underlying mean-field Hamiltonian H was combined with our quantum transport code, Gollum"】 |
| 声子输运 | **未点名**，只写 "Following the method described in [2,3]"（= Sadeghi/Sangtarash/Lambert Nano Lett 15, 7467 (2015)；Sadeghi, Nanotechnology 29, 373001 (2018)）→【推】同族自研 NEGF |
| 公式 | `T_ph = Tr[Γ_L G Γ_R G†]`，`G = (ω²I − D − Σ_L − Σ_R)^(−1)`【明】 |

### 2.2 具体论文

| 论文 | DFT 设置 | 关键数字 |
|---|---|---|
| Mosso et al., [Nano Lett. 19, 7614 (2019)](https://doi.org/10.1021/acs.nanolett.9b02089) | SIESTA / PBE / DZP / 250 Ry / 力 10 meV/Å | 理论 κ_ph 在 T>200 K 饱和于 **19 pW/K**；实验 **20 ± 6 pW/K**；电子热导 OPE3 0.07、ODT 0.01–0.02 pW/K（可忽略） |
| Noori, Sangtarash, Sadeghi, [Appl. Sci. 11, 1066 (2021)](https://doi.org/10.3390/app11031066) | SIESTA / LDA(Ceperley–Alder) / DZP / 350 Ry / 力 5 meV/Å | OPE3 计算 19 pW/K vs 实验 20±6 |
| Zheng et al., [Nano Lett. (2025)](https://doi.org/10.1021/acs.nanolett.5c00225) | 同族 SIESTA + NEGF | 连接方式依赖的声子干涉 |
| Sadeghi, Sangtarash, Lambert, [arXiv:1510.05945](https://arxiv.org/abs/1510.05945) (2015) | SIESTA / PBE / DZP / 250 Ry / 力 20 meV/Å | **唯一把有限位移流程写全的一篇**：δ=0.01 Å 中心差分 + 动量守恒补对角项；烷烃 25.4/33.4/30.3/5.6 pW/K（N=1,2,4,8）；寡炔 15.6/9.2/7.7 |

### 2.3 这一支的结构特征

- **这就是本仓计划要走的流程的"隔壁版本"**：0.01 Å 中心差分、动量守恒补对角项、`G=(ω²I−D−Σ_L−Σ_R)^(−1)`、`T_ph=Tr[Γ_L G Γ_R G†]`，与本方案 §5 完全一致。
- 差别在电极：B 支用**有限团簇**（金字塔/团簇）而非周期板，且没有明确做 decimation（正文未给电极自能的构造细节）。
- 前因子写法 `(1/2π)∫dω ħω τ ∂n/∂T` 与 A 支的 `(1/h)∫dE E τ ∂n/∂T` **代数等价**（E=ħω），不是两套物理。

---

## §3 主线 C：经典 MD + 经验势 + NEMD

### 3.1 技术栈

| 环节 | 具体 |
|---|---|
| MD 引擎 | **LAMMPS**（8 篇里 7 篇明确；**没有 GROMACS、没有 GPUMD、没有自研全栈**）；**版本号一律不披露** |
| 势函数分层 | Au 用 **EAM**（Grochola 2005 常见）；分子内谐波键角 + 余弦二面角 或 ReaxFF/REBO/AIREBO；界面用 **Morse 或 LJ**；**H 普遍隐式（联合原子）**；**几乎不写库仑项** ⇒ 无极化、无电荷转移 |
| 方法 | **NEMD 家族**（RNEMD 反向热流 / 直接 NEMD / Langevin-NEMD / AEMD Newton 冷却）。**没有一篇用 HNEMD**；GK 只在不含结的 EMD 姊妹篇出现 |
| 规模 | 双侧电极各 ~2×10³ 原子，总体系 10³–10⁴ 原子，2–50 ns |
| 后处理 | OVITO；G = (恒温器累计交换能量对时间求斜率) / ΔT |

### 3.2 具体论文

| 论文 | 力场 | 关键数字 |
|---|---|---|
| Wang, Gong, McGaughey, Segal, [arXiv:2209.06266](https://arxiv.org/abs/2209.06266) / JCP (2022) | **Morse Au–S：D_e=0.38 eV, α=1.47 Å⁻¹, r₀=2.65 Å**；EAM Au；无电荷项 | G_Au-Au = **40–45 pW/K @300 K**（RNEMD 与 AEMD 一致）；弱接触 20–22 |
| Li, Wang, Segal, [arXiv:2402.06172](https://arxiv.org/abs/2402.06172) / JPCM (2024) | AIREBO(C–C) + EAM(Au) + LJ(Au–C)，**物理吸附无 Au–S 键** | C60 单体 ≈50、二聚体 ≈41.5、三聚体 10.0 pW/K；**DFT-NEGF 给二聚体 7.0 pW/K → 差 6 倍**，作者归因于唯象势 |
| Wang, Gerry, Segal, [arXiv:2311.07830](https://arxiv.org/abs/2311.07830) / JCP (2024) | 同上 | **软件实现审计**：多原子 Langevin 恒温器算出的热交换分布违反涨落对称性，换单原子恒温器才恢复 ⇒ "平均值容易对，涨落可能系统性错" |
| Li, Hu, Yang, Lü, [PRB 104, 245413 (2021)](https://arxiv.org/abs/2202.06490) | **ReaxFF（Järvi 2011 + Bae & Aikens 2013）**，含 QEq 电荷转移；石墨烯用 REBO | 自写 MD 主程序 + LAMMPS 提供力；SCLMD（给恒温器喂有色噪声强制 Bose–Einstein）；**非谐使 κ 比同势的谐波 NEGF 低 ~15%（T>100 K）** |
| Luo & Lloyd, [arXiv:0905.2154](https://arxiv.org/abs/0905.2154) / IJHMT (2010) | Hautman–Klein 联合原子 + Morse(Au–Au, Au–S) | 速度标定 NEMD；界面热导 340 MW/m²K；**做了 Au–S Morse 阱深 ±10/20/30% 的敏感性扫描** ⇒ 等于承认 Au–S 键强是可调参数 |
| Yelishala et al., [Nat. Mater. 24, 1258 (2025)](https://doi.org/10.1038/s41563-025-02195-w) | **GAUSSIAN16（B3LYP/cc-pVDZ + D3(BJ)）参数化 QM 力场 JOYCE**；Sheng EAM(Au–Au)；Au–S Morse（Mahaffy 1997 / Ventura-Macias 2024 vdW-DFT 拟合）；INTERFACE FF | LAMMPS NEMD，每块 5600 原子/40 层，Langevin 2 ps⁻¹，0.25 fs；**明确主张"DFT 谐波忽略了非谐"**，是当前与 A 支唯一明显分歧点 |
| Wang & Segal, [arXiv:2512.12478](https://arxiv.org/abs/2512.12478) (2026) | **SIESTA AIMD → DeePMD-kit → LAMMPS** | 承认经验势不可扩展；**代价是把金属电极整个删掉**，只算分子本征热导 |

### 3.3 Au–S 力场参数到底哪来的（三条谱系）

- **谱系 A：文献传承的 3 参数 Morse，没有 DFT 拟合**。Segal 组用 D_e=0.38 eV / α=1.47 Å⁻¹ / r₀=2.65 Å，引自 Ong et al. JPCC 118, 7288 (2014)。**注意一处不一致**：同组 Sierra Suarez et al., J. Appl. Phys. 119, 145306 (2016) 表 I 写同一组 Morse 但 **α=1.67 Å⁻¹**，与 1.47 差 ~14%，引用前必须核对。Luo & Lloyd 则引 Sung & Kim (2005) 与 Mahaffy (1997)。
- **谱系 B：ReaxFF（真正 DFT 参数化）**，Järvi et al. JPCA 115, 10315 (2011) + Bae & Aikens JPCA 117, 10438 (2013)，含 QEq 电荷转移。只有 Li/Hu/Yang/Lü 一线用。
- **谱系 C：MLIP**，把 Au–S 交给神经网络势学出来，但目前丢掉了金属电极。

---

## §4 零星但重要的其他栈

| 工作 | DFT/FC | 输运 | 说明 |
|---|---|---|---|
| Markussen, [JCP 139, 244101 (2013)](https://arxiv.org/abs/1308.2543) | **GPAW**（PAW/DZP/PBE）；电极 FC 用 GULP/Brenner 或 Mingo 解析表面 GF | **自写 Green 函数码**（非 ATK） | 作者当时在 DTU，之后才去 QuantumWise |
| Zhang, Polanco, Ghosh, [arXiv:1612.04807](https://arxiv.org/abs/1612.04807) | **Quantum ESPRESSO**（LDA，35/300 Ry，8×8×8 k；DFPT→实空间 FC，Au 到第 2 近邻、聚乙烯到第 7 近邻） | 自写 NEGF | 少见地把"FC 截断到第几近邻"写清楚 |
| Li, Duchemin, Xiong, Solomon, Donadio, [JPCC 119, 24636 (2015)](https://arxiv.org/abs/1511.06058) | **QE**（PBE，ultrasoft，±0.02 Å 有限差分） | **elastic scattering kernel**（非 NEGF，Duchemin & Donadio PRB 84, 115423） | 另一条数学等价路线 |
| Nikolić, Saha, Markussen, Thygesen, [J. Comput. Electron. 11, 78 (2012)](https://doi.org/10.1007/s10825-012-0386-y) | **GPAW**；对照 GULP+Brenner | NEGF | |
| Lü & Wang, [PRB 78, 235436 (2008)](https://arxiv.org/abs/0807.1193) | **Gaussian03** | 自写 NEGF + mode matching | **最早的"量子化学团簇 + 自写 NEGF"实例**；因 Gaussian 不能做 PBC，他们优化 7 个 SWCNT 周期、只取中心周期的 FC |

---

## §5 四个"反直觉"的核实结果

### 5.1 GOLLUM 2014 原文没有声子模块

[Ferrer et al., NJP 16, 093029 (2014)](https://arxiv.org/abs/1502.04966) 结论段逐字：**"Further developments are in the pipeline, including the incorporation of phonon transport."**

标题里的 "thermal" 指的是**电子热导**（Wiedemann–Franz，由透射函数给），不是声子。声子模式是后续版本才加的（手册 mode 5，读 SIESTA `.FC`），但**检索范围内没有一篇分子结 κ_ph 论文明确写"we used GOLLUM to compute the phonon transmission"**——包括 Mosso 2019（只写 "following the method described in [2,3]"）。

⇒ 本仓调研报告 §8 把 GOLLUM 列为"直接出热导"，**证据链是断的**；应改成"GOLLUM 后续版本有声子模式（手册可查），但已发表文献未见点名使用"。

### 5.2 SIESTA + PHtrans 的真实输入是 GULP，不是 SIESTA

- PHtrans 确实存在，官方原话 "**NOTE: PHtrans is not as tested as TBtrans.**"
- Papior 在 siesta-l（2019-03-31）原话：用 SIESTA 做声子透射 "**extremely costly** since it requires a huge system to ensure the phonons are 'relaxed' in the electrode region. **This is why we (in the paper) did phonon transmissions using GULP**"。
- sisl 的 `examples/ex_03.py` 就是 GULP → PHtrans（ZZ 石墨烯）。
- Papior 在 Matter Modeling SE："**We are working on extracting the DFT related dynamical matrices, but still work to do.**"

⇒ 本仓计划 D2 选"自写 NEGF"而非 PHtrans，这条证据是**加分**：PHtrans 连 SIESTA 自己的动力学矩阵都还没接上。

### 5.3 QuantumATK 有模块、零文献

`PhononTransmissionSpectrum(dynamical_matrix, configuration, ...)` + `thermalConductance(phonon_temperature=300*Kelvin)`（返回 nW/K）确实存在、文档完整。但：

- 官方唯一相关教程是 **CNT 同位素热电**（ForceField 声子 + SE 电子），不是分子结；
- QuantumATK 论坛 2016 年员工引的自家声子输运文献清单**以 NEMD + Tersoff 为主**；
- 多轮检索（web / arXiv / OpenAlex / Crossref / DBLP）**未出现任何一篇点名该模块算单分子结 κ_ph 的论文**。

⇒ "用商业软件算分子结热导"在已发表文献里**基本是空集**。商业软件更多出现在"分子结电子输运/热电系数"。

### 5.4 VASP + phonopy 在这条线上缺席

能确证的"DFT 出力常数 + 自写输运码"实例，DFT 侧全是 **TURBOMOLE / SIESTA / GPAW / QE / Gaussian**，**没有 VASP**。

原因：phonopy 产出的超胞 FC 只解决"中心区"；分子结弹道声子透射还需要**半无限电极的表面格林函数 + 力常数拼接**（Mosso 的有限差分、Klöckner 的 K_LL 单独体相计算、Markussen 的解析表面 GF 都在处理这件事）。这条链在 LCAO/量子化学码上现成，在 VASP 的平面波框架上要另写一套。

⇒ 调研报告 §8 的"已有 VASP → VASP+phonopy 转 SIESTA .FC → GOLLUM"组合，**在文献里找不到先例**，风险被低估了。

---

## §6 力常数成本：文献怎么说

- 所有谱系一致：**动力学矩阵/力常数是最贵的一步**。Müller, Cuevas, Pauly, Nielaba, [PRB 106, 195401 (2022)](https://arxiv.org/abs/2208.09171) 逐字："the dynamical matrix ... **is typically the most time-demanding computational step**"。
- 这是**唯一给出 core-hour 的论文**（LAMMPS + EAM，MD 关联矩阵反演 FC）：NEMD 单次 288 core-h × 20 次 × 14 温度 ≈ **80,640 core-h**；温度依赖 FC 路线 ≈ **13,440 core-h**（快 ~6×）；固定温度 FC ≈ **960 core-h**（快 ~84×）。
- FC 收敛检验做得最扎实的也是它：截断半径 r_c 从 2 扫到 52 Å（monomer 6–14 Å 平台、8 原子链 6–28 Å 平台 ⇒ 取 10 Å），并做 3/6/12/18/36 层有限尺寸检验（变化仅几 pW/K）。
- **DFPT vs 有限位移的成本差是数量级的**：A 支一次 DFPT 拿全 Hessian；B 支 6N 次力评估。这个差别来自**DFT 引擎能力**，不是方法要求。

---

## §7 对本仓方案（`docs/thermal-conductance-plan-2026-11-23.md`）的含义

1. **本方案在结构上就是 A 支**：量子化学/团簇基 DFT 产力常数 + 自研 NEGF + 半无限电极自能 + decimation。只是把 TURBOMOLE 换成 FHI-aims。**这条路有 2016–2025 的连续先例**，不是自创。
2. **唯一的结构性差异**：FHI-aims 221103 无 DFPT（且 ≥240507 的 DFPT 对周期体系只支持 LDA），所以必须 6N 有限位移。A 支一次 DFPT 就拿到全 Hessian，我们要 324 次 SCF。**这是引擎限制，不是方法限制**——方案 §12.1 的结论（不升级）由此得到文献侧的支持：连官方 DFPT 论文都建议大体系/其他泛函走有限位移（FHI-vibes）。
3. **B 支有两个细节可以直接抄**：δ=0.01 Å 中心差分；**对角项用动量守恒补齐**（省 3N 次位移）。
4. **没有 FHI-aims 基的分子结 κ_ph 先例（【缺】）**——这既是新颖性，也是风险：**没有第三方数值可对账**。因此方案 §9 的 T1–T13 内部一致性检验（声学求和规则、带顶截断、电极层数收敛、弹道不变性）必须做足，不能靠"对文献"验收。
5. **对账基准的选择**：最可比的是 Mosso 2019（SIESTA/PBE/DZP，19 pW/K）与 Cui 2019（TURBOMOLE/PBE/def2-SV(P)，~20 pW/K）——两者泛函相同、基组相当、数值一致。本仓算 Au–烷二硫醇时应落在 15–25 pW/K 区间；若偏离一个量级，先查单位换算（附录 A-1）与带顶截断（A-7）。

---

## §8 检索范围与未核实项

**已读全文（开放获取）**：Cui 2019（Konstanz OA 副本）、Mosso 2019（Zenodo submitted version + SI）、Yelishala 2025（CSIC OA PDF）、Klöckner PRB 94/95/96、Blaschke & Pauly 2025、Sadeghi 2015、Müller 2022、Lü & Wang 2008、Ferrer 2014、Papior CPC 2017、Noori 2021、Zheng 2025、Markussen 2013、Zhang 2018、Li 2015、Nikolić 2012、Wang 2022/2024、Li 2024、Wang & Segal 2026、Li/Hu/Yang/Lü 2021、Luo & Lloyd 2010a/b。

**未取得全文**：Diao et al. IJHMT 2020（付费墙，只有摘要；且其体系是 Au–CNT 而非 Au–分子–Au）、Gotsmann/Gemma/Segal APL 120, 160503 (2022) 综述、Sadeghi JPCC 123, 12556 (2019)、Bro-Jørgensen et al. ACS Phys. Chem. Au 2025。

**负面结论的边界**：§5.3/§5.4 的"未找到"是 web + arXiv + OpenAlex + Crossref + DBLP + Europe PMC + 机构库 + Zenodo/Unpaywall 检索范围内的负面，**不等价于穷尽证明**。
