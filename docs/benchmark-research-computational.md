# 单分子结输运计算 Benchmark 表

> 文献计算值、实验锚点与 smj 0.2.0 的对齐口径  
> 工作稿 v0.1 · 2026-11-15 · 维护：bench-comp

## 0. 用途与一句话结论

本表为 `smj 0.2.0` 的 PBE/light + 37 Au/侧团簇 + 3 层 + T=300 K + σ_near/σ_far=0.05/0.10 Ha 结果提供**可追溯的文献坐标**，不是“唯一正确值”排行榜。所有文献数值均从一手论文的正文、表格或图注抄录，并给出 DOI 与 locator；拿不到全文的条目标“摘要级/未核验”，不进入数值列。领域背景与实验/理论差异的综述见 Evers et al. 2020（DOI 10.1103/RevModPhys.92.035001）。

一句话结论：

- **BDT 是唯一可以“同锚定 + 同泛函 + 同结构类”对齐的体系**：smj `bdt` L0=0.2789，Strange et al. 2008 的 PBE/DZP 周期板值为 0.24 G0（+16%）；量级一致，但 smj 的透射峰在 −0.339/−0.687 eV，Strange 的共振更靠近 E_F，G 相同不证明机制相同。
- **AZ 没有可用的绝对理论锚点**：Xia et al. 2014 的 2,6-azulene 实验值 3.2×10⁻⁴ G0 用的是 thiochroman 连接子，不是直接硫醇盐，不能拿来判 smj `az` 的 0.3806 对错。
- **2CN 没有同分子锚点**：Gulyaev et al. 2026 的 SMePh（硫醚锚）为 1.3×10⁻² G0，但分子只有苯环长度；smj `2cn` 是 C30 共轭分子，0.001325 低约 10× 可由长度/耦合解释，不能直接判误差。
- **PBE/light 对胺、腈基、吡啶锚定量不可靠**：文献中 PBE 对这些体系高估 3–226×；GW/DFT+Σ 通常可把误差压到 2× 内（BDCN 的 DFT+Σ 仍高 5.6×）。对裸硫醇盐 BDT，GW 反而把 G 推到 0.83 G0，说明“实验测到的 BDT”很可能不是裸硫醇盐，而是 H 饱和/Au 加合物结构。

## 1. 口径、单位与缩写

| 项 | 定义 |
|---|---|
| G0 | 2e²/h ≈ 77.48 μS（CODATA；Li et al. 用 77.5 μS） |
| G/G0 | 零偏压 Landauer 电导；在无特别说明时等于 T(E_F) |
| T(E_F) | 费米能级处的透射系数；理论值通常是 T=0 K、单一几何 |
| L0 | smj 的 Landauer 积分 ∫T(E)(−∂f/∂E)dE，T=300 K；`G/G0 = L0` |
| 计算/实验 | 本表比值均为“计算值 ÷ 实验值”，标“≈”；实验本身有分布，不是精确靶 |
| [全文] | 已读到原文表格/图/正文，数值可定位 |
| [摘要级] | 只有摘要，未见数值表 |
| [未核验] | 未拿到全文，不用于结论 |
| 拟合值 | 原文的系列拟合参数（如 Gc、β），不是单一结的 G |

**缩写**：PBE/GGA、LDA 为交换关联泛函；GW 为多体自能（不是泛函）；DFT+Σ 为“DFT 透射 + 轨道能量/展宽修正”（常含镜像电荷）；ASIC/SIC 为原子自相互作用修正。PAO/NAO 为数值原子轨道；DZP 为双 ζ 加极化；SZ 为单 ζ；TZVP、def-SVP 为量子化学基组。NEGF 为非平衡格林函数；Landauer 为散射态/透射公式（不一定是 NEGF）。周期板指二维周期 Au(111) 电极；团簇指有限 Au 团簇电极。

**结构类型计数**：本表主表覆盖 **7 类**——①芳香二硫醇（BDT、BDT+H）、②胺锚（BDA）、③腈基（BDCN）、④吡啶锚（BP）、⑤饱和烷烃（HDA）、⑥富勒烯（C60）、⑦硫醚（SMePh）。非交替芳烃（azulene）与 OPE 另见第 3 节说明：前者没有可用的绝对计算值，后者全文未取到，故不强行补数值行。

## 2. 文献计算 benchmark（Table 1）

**读表说明**：每行给出一个“主方法”；同一行用分号列出的多个 G 值均来自同一篇论文的同一张表/同一组计算。G 列保留原文有效数字；实验锚点列给出典型实验值与计算/实验比。电极模型列只写论文明确给出的层数、超胞和接触位点；原文未写的，不猜。

| # | 体系（分子 + 电极） | 结构类型 | 计算方法：泛函 / 基组 / 代码 / NEGF | 电极模型 | 计算 G/G0（原文值） | 实验锚点与计算/实验比 | 一手文献（DOI + locator） | 备注（自能/镜像电荷、基组/尺寸收敛） |
|---|---|---|---|---|---|---|---|---|
| 1 | Au–BDT–Au(111)（1,4-苯二硫醇，硫醇盐） | 芳香二硫醇 | PBE；SIESTA 模守恒 PAO，DZP（S/C/H）；平面波 + MLWF 交叉验证；NEGF 双探针 | 2D 周期 Au(111) 板；3–4 原子层收敛；S 置 fcc hollow；Au 晶格 4.18 Å，面间距 9.68 Å；d(Au–S)=2.45 Å | **0.24**（DZP，4×4 k∥）；Γ-only **0.37**；SZ 明显偏低 | 实验典型 0.011（范围 10⁻⁴–10⁻¹）；0.24/0.011 ≈ 22 | Strange et al. 2008, *J. Chem. Phys.* **128**, 114714；DOI 10.1063/1.2839275；Sec. VI、Fig. 6(b) | 4×4 k∥ 与 Γ-only 差约 50%；DZP 才定量；3–4 层收敛；纯 PBE，无自能/镜像修正 |
| 2 | Au–BDT–Au(111)（硫醇盐，fcc hollow 略倾 bridge） | 芳香二硫醇 | 自洽 GW；GPAW PAW/LCAO，分子 DZ、Au DZP；NEGF + GW 自能 | 7 层 4×4 Au(111) 周期板；4×4 k；表面格林函数用 4×4 超胞、k⊥=(0.1,0.2)；面间距 10.20 Å | **0.83**（GW）；同文 DFT-PBE **0.28** | 实验 0.011；0.83/0.011 ≈ 75（GW 反而更差） | Strange et al. 2011, *Phys. Rev. B* **83**, 115108；DOI 10.1103/PhysRevB.83.115108；Table 2、Sec. III.3 | GW 含动态屏蔽；硫醇盐共振钉在 E_F 附近，GW 不降反升；说明裸硫醇盐不是实验结构 |
| 3 | Au–BDT–Au(111)（hollow / bridge / top / adatom） | 芳香二硫醇 | LDA/GGA vs ASIC（原子自相互作用修正）；SIESTA PAO；SMEAGOL NEGF | Au(111) 周期电极（SMEAGOL 双探针）；hollow 的 S–表面距离 d_SS=1.9 Å（对应 Au–S≈2.53 Å）；另测 bridge/top/adatom | **0.23**（LDA hollow）；**0.06**（ASIC hollow）；bridge LDA **0.10** / ASIC **0.06** | 实验 0.011（Xiao）；0.06/0.011 ≈ 5.5；0.23/0.011 ≈ 21 | Toher & Sanvito 2008, *Phys. Rev. B* **77**, 155402；DOI 10.1103/PhysRevB.77.155402；Sec. III.1–III.3、Fig. 5/9 | ASIC 打开 HOMO–LUMO gap、压低 G；对 BDMT 影响小；电极层数原文未在正文强调 |
| 4 | Au–BDT+H–Au(111)（氢饱和硫醇，经 Au 尖端） | 芳香二硫醇（thiol，H 饱和） | 自洽 GW；GPAW；NEGF | 7 层 4×4 Au(111)；两侧 Au 尖端；Au–Au 9.80 Å，Au–S 2.47 Å | **0.010**（GW）；同文 DFT-PBE **0.054** | 实验 0.011；0.010/0.011 ≈ 0.9 | Strange et al. 2011, *Phys. Rev. B* **83**, 115108；Table 2、Sec. III.7 | H 饱和 S 模拟 Au 加合物/化学钝化；GW 与实验吻合；锚定化学比泛函更关键 |
| 5 | Au–BDA–Au(111)（1,4-苯二胺） | 胺锚 | PBE/GGA；SIESTA；Au d 单 ζ、其余 DZP；NEGF/相干散射态 | 6 层 Au(111) 板 ×2；16 Au/层；中心区两侧各 4 Au 层；8×8 k∥（部分 16×16）；胺 atop 欠配位 Au | **0.046**（15 个几何平均；范围 0.032–0.062） | 实验 0.0064±0.0005；0.046/0.0064 ≈ 7（原文明确“7×”） | Quek et al. 2007, *Nano Lett.* **7**, 3477–3482；DOI 10.1021/nl072058i；Fig. 4(a)、摘要 | HOMO 共振 −1.27…−0.94 eV、宽 0.34–0.56 eV；GGA 自能误差使共振过浅；镜像电荷估计需更深 |
| 6 | Au–BDA–Au（经 Au 尖端） | 胺锚 | 自洽 GW；GPAW；NEGF | 7 层 4×4 Au(111)；两侧 Au 尖端；Au–Au 9.90 Å，Au–N 2.55 Å | **0.0042**（GW）；同文 DFT **0.021** | 实验 0.007（0.0064–0.007）；0.0042/0.0064 ≈ 0.66 | Strange et al. 2011, *Phys. Rev. B* **83**, 115108；Table 2 | GW 比 DFT 低 5×，与实验同量级；尖端几何 |
| 7 | Au–BDA–Au | 胺锚 | PBE + DFT+Σ（Σ_GP + Σ_CORR 镜像电荷）；GPAW LCAO DZP；NEGF | 8 层 Au(111)，16 原子/层；(4,4,1) k；tip/adatom/trimer 三种接触 | **5.7×10⁻³**（+1 Å 像面）；2.9×10⁻³（−1 Å）；DFT-PBE **24×10⁻³** | 实验 6.4×10⁻³；5.7/6.4 ≈ 0.89；24/6.4 ≈ 3.8 | Markussen et al. 2013, *Phys. Status Solidi B* **250**, 2394–2402；DOI 10.1002/pssb.201349217；Table 2 | 像面 ±1 Å 为不确定度；DFT+Σ 把误差压到 2× 内；BDA 的 DFT 24×10⁻³ 与 Quek 46×10⁻³ 的差异来自接触 motif/平均 |
| 8 | Au–BDCN–Au（1,4-苯二甲腈，腈基锚） | 腈基 | PBE / GW / DFT+Σ；GPAW；NEGF | 8 层 Au(111)，16 原子/层；(4,4,1) k；tip/adatom/trimer | **1.9×10⁻²**（DFT-PBE）；**6.3×10⁻⁵**（GW）；**4.7×10⁻⁴**（DFT+Σ +1 Å）/ 2.2×10⁻⁴（−1 Å） | 实验 8.4×10⁻⁵；1.9×10⁻²/8.4×10⁻⁵ ≈ 226；6.3×10⁻⁵/8.4×10⁻⁵ ≈ 0.75；4.7×10⁻⁴/8.4×10⁻⁵ ≈ 5.6 | Markussen et al. 2013, *Phys. Status Solidi B* **250**, 2394–2402；Table 2、Sec. 3.2 | DFT 因 LUMO 严重低估而高估 >2 个数量级；GW 低 25%；DFT+Σ 仍高 5.6× |
| 9 | Au–BP–Au（4,4′-联吡啶） | 吡啶锚 | PBE / DFT+Σ；SIESTA DZP；NEGF | 7 层 Au(111) ×16 Au；最后 4 层固定；3 种接触基元；BP 取 HighG | **2.2×10⁻³–9.4×10⁻³**（DFT-PBE）；**1.5×10⁻⁴–2.8×10⁻⁴**（DFT+Σ） | 实验 6×10⁻⁴；DFT/实验 ≈ 3.7–15.7；DFT+Σ/实验 ≈ 0.25–0.47 | Darancet et al. 2012, *Nano Lett.* **12**, 6250–6254；DOI 10.1021/nl3033137；Table 1、Fig. 1 | DFT 对 DAS/APA/BP 高估 >10×；DFT+Σ 把误差压到 2× 内；共振移动 ≥1.5 eV |
| 10 | Au–HDA–Au（1,6-己二胺） | 饱和烷烃（二胺） | PBE / DFT+Σ；SIESTA DZP；NEGF | 7 层 Au(111) ×16 Au；3 种接触基元 | **5.1×10⁻⁴–5.4×10⁻⁴**（DFT-PBE）；**0.7×10⁻⁴**（DFT+Σ） | 实验 1.2×10⁻⁴；DFT/实验 ≈ 4.3–4.5；DFT+Σ/实验 ≈ 0.58 | Darancet et al. 2012, *Nano Lett.* **12**, 6250–6254；Table 1 | 烷烃无低能共振，DFT 误差仅 4×；DFT+Σ 略低；与 BDT/BDA 的误差方向不同 |
| 11 | Au–C60–Au | 富勒烯 | BP86；TURBOMOLE 6.1 def-SVP；有限团簇 NEGF（Pauly et al. 2008） | 有限 Au 团簇：先 20/19 Au 优化，再扩到约 63 Au/侧；top vs hollow C60 结合；无 vdW | **1.85**（hollow）；**0.55**（top）；拉伸平台 **0.07–0.2**；Au–BDC60–Au **9.0×10⁻³**（hollow）/ 2.5×10⁻⁶（top） | 实验约 0.1（10 K）；1.85/0.1 ≈ 18.5；0.55/0.1 ≈ 5.5；平台与 0.1 相符；BDC60 实验 3×10⁻⁴ | Bilan et al. 2012, *Phys. Rev. B* **85**, 205403；DOI 10.1103/PhysRevB.85.205403；Fig. 1(b–d)、Fig. 2、Sec. III | LUMO 主导；有限团簇、无 vdW；室温 STM 另有约 5×10⁻⁴（大展宽）；团簇 20→63 Au 影响能级对齐 |
| 12 | Au–SMePh–Au（双(硫甲基)苯，两个 SMe 锚；硫醚锚） | 硫醚 | PBE 几何 FHI-aims；AITRANSS 透射；PBE→PBE0/B3LYP 轨道移动 + Au DOS 缩放（DFT+Σ 近似） | Au21 金字塔电极/侧；FHI-aims NAO；具体基组 tier 原文未强调 | **1.3×10⁻²**（PBE 透射 + PBE0 修正）；直接 PBE0 **1.6×10⁻²** | 实验 1.2×10⁻²；1.3/1.2 ≈ 1.08；1.6/1.2 ≈ 1.33 | Gulyaev et al. 2026, *Nano Lett.* **26**, 7429–7434；DOI 10.1021/acs.nanolett.6c01462；Table 1、Fig. 3c | 与 smj 同代码链（FHI-aims + AITRANSS）；Breit–Wigner 拟合 + 高 rung 泛函位移 + Au DOS 缩放；无显式镜像电荷；修正为经验型 |

### Table 1b 同一体系的方法学横向对照（单位 G0）

把同一体系、不同方法的数放在一起，才能看出“方法误差”而不是“体系差异”。下表全部来自 Table 1 的同一批论文。

| 体系 | 纯 DFT（PBE/LDA） | GW | DFT+Σ | ASIC | 实验 |
|---|---|---|---|---|---|
| BDT 硫醇盐 | 0.23–0.28 | **0.83** | — | 0.06 | 0.011 |
| BDT+H 硫醇 | 0.054 | **0.010** | — | — | 0.011 |
| BDA 胺 | 0.021–0.046 | **0.0042** | 0.0057 | — | 0.0064 |

三个要点：

1. **同一体系不同方法可差 10× 以上**：BDT 硫醇盐从 ASIC 0.06 到 GW 0.83，跨越 14×；BDA 从 PBE 0.046 到 GW 0.0042，跨越 11×。
2. **“GW 更好”不是普遍规律**：BDA/BDT+H 的 GW 与实验吻合；BDT 硫醇盐的 GW 反而把 G 推高到 0.83，因为硫醇盐的共振被钉在 E_F 附近，GW 只移动了能级、没有去掉共振。
3. **锚定化学比泛函更关键**：BDT 硫醇盐 PBE 0.24、BDT+H PBE 0.054，同一泛函下差 4.4×；GW 下 0.83 vs 0.010，差 83×。实验 0.011 更像 BDT+H/Au 加合物，而不是裸硫醇盐。

### 读表注（逐行补充）

- **第 1 行 BDT/PBE**：Strange 的 0.24 是 DZP、4×4 k∥ 的值；Γ-only 0.37，说明 k∥ 收敛不能省。SZ 基组明显偏低，DZP 才定量。这是 smj `bdt` 最直接的同锚定、同泛函参照。
- **第 2 行 BDT/GW**：同文 DFT 0.28、GW 0.83；GW 把硫醇盐的透射峰钉得更靠近 E_F。该文自己结论是“裸硫醇盐不能解释实验 0.01”，需要 H 饱和或 Au 加合物。
- **第 3 行 BDT/SIC**：LDA hollow 0.23、ASIC hollow 0.06；bridge 0.10/0.06。ASIC 打开 gap，对 BDT 影响大、对 BDMT 影响小。d_SS=1.9 Å 对应 Au–S≈2.53 Å，比 smj 的 2.253 Å 长，耦合更弱。
- **第 4 行 BDT+H/GW**：0.010，与实验 0.011 几乎一致；PBE 0.054。这个条目是“锚定化学决定 G”的最强证据。
- **第 5 行 BDA/Quek**：0.046 是 15 个几何的平均（0.032–0.062），实验 0.0064，原文明确 7× 高估。作者用镜像电荷模型说明共振应更深，GGA 自能误差是主因。
- **第 6 行 BDA/GW**：0.0042，实验 0.0064；GW 比 DFT 低 5×，方向正确。
- **第 7 行 BDA/DFT+Σ**：+1 Å 像面 5.7×10⁻³，−1 Å 2.9×10⁻³；实验 6.4×10⁻³。像面位置是不确定度，但都在 2× 内。BDA 的 DFT 24×10⁻³ 与 Quek 46×10⁻³ 的差异来自接触 motif 和几何平均，不是方法矛盾。
- **第 8 行 BDCN**：PBE 1.9×10⁻²，实验 8.4×10⁻⁵，高估 226×；GW 6.3×10⁻⁵（低 25%）；DFT+Σ 4.7×10⁻⁴（高 5.6×）。腈基锚的 LUMO 问题比胺更严重。
- **第 9 行 BP**：PBE 2.2×10⁻³–9.4×10⁻³（接触 motif 造成 4× 展宽），实验 6×10⁻⁴；DFT+Σ 1.5–2.8×10⁻⁴。DFT 对 BP 高估 >10×，DFT+Σ 压到 2× 内。
- **第 10 行 HDA**：PBE 5.1–5.4×10⁻⁴，实验 1.2×10⁻⁴（仅 4× 高估）；DFT+Σ 0.7×10⁻⁴。烷烃没有低能共振，误差比胺/腈小。
- **第 11 行 C60**：hollow 1.85、top 0.55；拉伸平台 0.07–0.2 与 10 K 实验约 0.1 相符。BDC60 的 hollow 9.0×10⁻³、top 2.5×10⁻⁶，实验优先值 3×10⁻⁴；说明 C60 的接触几何决定一切。
- **第 12 行 SMePh**：PBE 透射 + PBE0 位移 + Au DOS 缩放给 1.3×10⁻²，实验 1.2×10⁻²；直接 PBE0 1.6×10⁻²。与 smj 同代码链，但 smj 没有这步修正，不能直接拿 1.3×10⁻² 当“PBE/light 应该等于”的靶。

## 3. 补充表（Table 2）：系列拟合、长度序列与实验锚点

以下数值**不进入主表**：要么是系列拟合参数，要么是实验值，要么是同一方法的长度序列。列在这里是为了给“趋势/量级”参照，避免把它们误当成单结计算 G。

| 体系/系列 | 类型 | 方法/电极 | 结果（原文值） | 文献 + locator | 备注 |
|---|---|---|---|---|---|
| Au–(CH2)n–Au，n=4–10（烷烃二硫醇） | 饱和烷烃（二硫醇） | BP86/TZVP，TURBOMOLE；Au 金字塔 40–55 原子；Landauer（非 NEGF） | 接触电导 **Gc=0.24 G0=18.5 μS**；**βN=0.83/CH2**；HOMO–LUMO gap 7.5 eV；ΦB=2.14 eV；实验 β=0.94–0.96 | Li et al. 2007, *J. Am. Chem. Soc.* **130**, 318–326；DOI 10.1021/ja0762386；arXiv:0802.2407 p.10、Fig. 8A | Gc/β 是系列拟合，不是单一结 G；S atop/bridge/hollow 的 d(Au–S)=2.28/2.36/2.72 Å；all-trans 与 gauche 差 4×；B3LYP 把 β 推到 1.0 |
| Au–BDA/DBDA/TBDA/TetBDA–Au（胺系列） | 胺锚，长度序列 | PBE 几何 FHI-aims；Au21 金字塔；AITRANSS；PBE0 修正 | 实验：7.0×10⁻³ / 1.2×10⁻³ / 3.1×10⁻⁴ / 6.0×10⁻⁵；PBE 派生修正：5.4×10⁻³ / 1.2×10⁻³ / 2.7×10⁻⁴ / 4.6×10⁻⁵；直接 PBE0：5.5×10⁻³ / 1.1×10⁻³ / 1.8×10⁻⁴ / 2.8×10⁻⁵ | Gulyaev et al. 2026, *Nano Lett.* **26**, 7429–7434；Table 1 | 同一方法跨 4 个长度：PBE 派生修正误差 0.8–1.0×；直接 PBE0 为 0.5–0.9×；说明“方法 + 长度趋势”可迁移，但绝对 G 仍依赖修正 |
| Au–azulene–Au（1,3-/1,3-Py-/2,6-/4,7-/5,7-） | 非交替芳烃 | 实验 STM-BJ；GW 只算相对趋势、省略 thiochroman linker | **1,3Az 3.2×10⁻⁴；1,3PyAz 9×10⁻⁵；2,6Az 3.2×10⁻⁴；4,7Az 8×10⁻⁵；5,7Az 2×10⁻⁵**（G0） | Xia et al. 2014, *Nano Lett.* **14**, 2941–2945；DOI 10.1021/nl5010702；Table 1 | 实验值，不是计算值；连接子是 thiochroman（含饱和 S 环），不是直接硫醇盐；不能拿来判 smj `az` 的绝对 G |

**未形成数值行的类型（说明而非凑数）**：

- **非交替芳烃/薁（azulene）**：Xia et al. 2014 有实验 Table 1，但 GW 只做相对趋势、模型省略了 thiochroman linker；Dutta et al. 2008（arXiv:0807.3510）只有 I–V/透射图，正文没有绝对 G 数值；El-Nahas et al. 2016/2017 的摘要没有 G 数值且正文付费。**没有可用的“同分子 + 绝对 G”计算值**，因此不设计算行。smj `az` 只能作为“预测”，不能作为“已验证”。
- **OPE**：Manrique et al. 2015（*Nat. Commun.*，DOI 10.1038/ncomms7389）全文读取失败（arXiv:1509.00990 的 HTML/PDF 均不可读），无法抄录原始 G；因主表已达 7 类，不强行补位，避免编数。若后续拿到全文，可补一行 pyridyl-OPE/Au。
- **Thygesen & Jacobsen 2005**（*Phys. Rev. B* **72**, 033401）是 k 点采样/超胞干涉的方法论文，不含分子结绝对 G；本表在 5.2 引用它说明 Γ-only 的 van Hove 伪峰。**Brandbyge 2002** 只算 C 线/Au 线/CNT 缺陷，不含分子结 G，故只作为 transiesta/NEGF 方法源头引用，不进入数值行。

## 4. 本项目 smj 0.2.0 实测值（Table 3）

> **标注：本工作/未发表。** 以下三行是 `smj 0.2.0` 的发布验证样本（`docs/validation-0.2.0.md`），不是生产结论；按 2026-11-15 口径，三套物理判定均为 `passed`（近 E_F 峰已按族归因：BDT/AZ 在 HOMO 侧、2CN 在 LUMO 侧）。方法为 FHI-aims 221103 + AITRANSS 071813，PBE/light，电极库 `fcc111.au.clusters.37_atoms`（37 Au/侧，3 层），T=300 K，σ_near/σ_far=0.05/0.10 Ha，能量窗 `[-0.22, 0.0002, -0.16]` Ha（301 点）。L0 为 300 K Landauer 积分；T(E_F) 为最近格点值。

| 体系（本工作/未发表） | 结构类型 | 方法 | 电极模型 | L0 = G/G0 | T(E_F)（最近点） | 文献对照 | 备注 |
|---|---|---|---|---|---|---|---|
| Au–BDT–Au（1,4-苯二硫醇→硫醇盐；输入 C6H6S2 14 原子，去质子后 C6H4S2 12 原子；总结 86 原子） | 芳香二硫醇 | FHI-aims 221103 + AITRANSS 071813；PBE/light；T=300 K；σ 0.05/0.10 Ha | 37 Au/侧 fcc111 团簇，3 层；S 三元组接触（[40,41,44]/[76,77,80]）；Au–S≈2.253 Å | **0.278853973** | **0.272587** @ −1.59 meV；efermi −0.188142 Ha | Strange 2008 PBE 0.24：+16%；实验 0.011：25× | 1.0 eV 内有两个显著峰：−0.339 eV（T=1.010）与 −0.687 eV（T=1.055），均在 BDT 预期 HOMO 侧；G 由更靠近 E_F 的 −0.339 eV 峰尾贡献 |
| Au–azulene-2,6-dithiol–Au（输入 C10H8S2 20 原子，去质子后 C10H6S2 18 原子；总结 92 原子） | 非交替芳烃（薁） | 同上 | 同上；Au–S 2.2511/2.2549 Å | **0.380578051** | **0.349526** @ +2.36 meV；efermi −0.188086748 Ha | 无同分子绝对理论值；Xia 2,6Az 实验 3.2×10⁻⁴（thiochroman 连接，不可直接比） | max T=1.1795 @ −0.406 eV；E_F 附近由更近的 −0.177 eV 峰尾贡献；与 Xia 的比值仅作数量级提示，不作判据 |
| Au–2CN-656-4-thiomethylphenyl–Au（C30H20N2S2，54 原子；S–CH3 硫醚锚；总结 128 原子） | 硫醚 + 腈基取代共轭分子 | 同上 | 同上；Au–S 2.4248 Å（配位接触，弱耦合） | **0.001324925** | **0.001253** @ +1.76 meV；efermi −0.179464717 Ha | Gulyaev SMePh 1.3×10⁻²：约低 10×；实验 1.2×10⁻²：约低 9×；BDCN DFT+Σ 4.7×10⁻⁴：高 2.8×（锚定不同） | 显著峰在 +0.241 eV（T=0.871），落在 LUMO 侧；2CN 无一手通道文献，按实测 S<0 自洽归入 LUMO 侧；弱耦合是硫醚配位接触的预期行为 |

**Table 3 的读法**：

- smj 的 `L0` 与 `T(E_F)` 不是同一个数：BDT 差 2.3%、AZ 差 8.9%、2CN 差 5.7%。与文献的 T=0 `T(E_F)` 比，应优先用 smj 的 `T(E_F)` 最近点；与 300 K 实验 G 比，才用 `L0`。
- `bdt` 的 0.2789 与 Strange 2008 的 0.24 只差 16%，但两者 T(E) 峰位不同（smj −0.339/−0.687 eV；Strange 更靠近 E_F）。**G 接近不等于机制一致**；如果后续要拿 smj 做“实验预测”，必须同时对齐峰位/峰宽，而不只是 G。
- `az` 的 0.3806 是“直接硫醇盐”的强耦合结果；Xia 的 3.2×10⁻⁴ 是 thiochroman 连接子的弱耦合结果，两者不是同一分子/同一锚定，**1190× 的差不能当成误差**。
- `2cn` 的 0.001325 与 SMePh 的 1.3×10⁻² 同为硫醚锚，但分子长度差很多；低 10× 与“更长、更弱耦合”的预期一致。它是**量级自洽**，不是**定量验证**。

## 5. 使用注意：方法差异与对齐口径

### 5.1 泛函/自能：误差方向取决于“E_F 落在峰上还是峰尾”

- **PBE/GGA 低估 HOMO–LUMO gap**，把分子共振拉向 E_F。若 E_F 落在共振峰附近，G 被高估；若 E_F 落在峰尾，误差方向可能反转。BDA 的 PBE 0.021–0.046 vs 实验 0.0064（高 3–7×）；BDCN 的 PBE 1.9×10⁻² vs 实验 8.4×10⁻⁵（高 226×）；BP 的 PBE 高 3.7–15.7×；HDA 只高 4.3–4.5×——同样是 PBE，误差从 4× 到 226×，不能用一个“PBE 修正因子”统一修。
- **GW/DFT+Σ 不是万能**：BDA 的 GW 0.0042、DFT+Σ 0.0057，都把误差压到 2× 内；BDCN 的 GW 0.75×、DFT+Σ 5.6×；但 BDT 硫醇盐的 GW 0.83 反而比 PBE 0.24 更远离实验 0.011。原因是硫醇盐的共振被钉在 E_F，GW 只重正化能级，没有去掉共振；ASIC 打开 gap 才能把 G 压到 0.06。
- **SIC/ASIC 是另一条修正路线**：BDT hollow 0.23（LDA）→0.06（ASIC），factor 3.8；BDMT 0.006→0.004，影响小。SIC 对“共振是否在 E_F”极敏感，不能与 GW/DFT+Σ 的结果混为一谈。
- **结论**：表中每一行的“方法”必须与 G 一起读；不要拿 PBE 的 0.24 与 GW 的 0.83 直接说“文献不一致”。它们是不同物理近似下的不同模型。

### 5.2 基组与 k 点：DZP 与 Γ-only 的差距不能忽略

- Strange 2008：SZ 明显偏低，DZP 才定量；4×4 k∥ 给 0.24，Γ-only 给 0.37（差 54%）。如果只做 Γ 点，BDT 会被高估约 50%。
- **Thygesen & Jacobsen 2005**（*Phys. Rev. B* **72**, 033401）用紧束缚模型 + Pt 接触证明：Γ-only 超胞会在透射中引入源于 lead 的 van Hove 伪峰，必须做 k∥ 采样；这正是 Strange 2008 的 Γ-only 0.37 与 4×4 的 0.24 差异来源。smj 没有 k∥，但团簇尺寸/取向有等价的采样问题。
- Quek 2007：Au d 用单 ζ、其余 DZP；8×8 k∥，部分体系需 16×16。Markussen 2013：DZP，(4,4,1) k。Darancet 2012：DZP，7 层 ×16 Au。
- smj 用 `light` 基组，收敛程度低于文献的 DZP/TZVP；因此 smj 与文献的 10–20% 一致可能是“误差抵消”，不能当成“light 足够”。若要把 smj 用于定量预测，至少要做一次 tight 基组对照（项目已有 light/tight 的 gap 差 0.01–0.02 eV 数据，但 G 的 tight 对照应单独记录）。
- **规则**：先确认基组/k 点收敛，再比较 G；未收敛时，把 <2× 的差异视为“同量级”，不要解读为物理差异。

### 5.3 电极模型：周期板、团簇、尖端、接触位点

- **周期板 vs 团簇**：Strange 用 3–4 层、Quek 6 层、Darancet 7 层、Markussen 8 层的二维周期 Au(111) 板；smj 用 37 Au/3 层团簇；Bilan 用约 63 Au/侧团簇；Gulyaev 用 Au21 金字塔。团簇的离散能级 + σ 自能展宽会改变 E_F 附近的 DOS 和能级对齐；Bilan 明确说 20→63 Au 的扩展是为了正确描述电荷转移和能级对齐。smj 的 37 Au/3 层比 Bilan 的 63 Au 更小，**可能欠屏蔽**。
- **尖端 vs 平面**：BDA 的 Quek 用欠配位 atop Au；Strange 2011 的 BDA 用 Au 尖端（Au–Au 9.90 Å）；两者都叫“胺锚”，但耦合不同。BDT 的 Strange 2008 用平面 fcc hollow；Toher 2008 还测了 bridge/top/adatom，G 从 0.23 到 0.06 不等。
- **接触位点**：BDT hollow 0.23（LDA）、bridge 0.10、ASIC 后都 0.06；C60 hollow 1.85、top 0.55（差 3.4×）；BDC60 hollow 9.0×10⁻³、top 2.5×10⁻⁶（差 3600×）。**接触几何是第一阶变量**，比泛函更敏感。
- **k∥**：4×4、8×8、16×16 不等；Γ-only 可差 50%。smj 的团簇没有 k∥ 概念，但团簇尺寸/取向等价地影响采样。
- **规则**：比较前先问“电极是周期还是团簇、几层、什么接触位点、k∥ 多少”；不同电极模型的差异通常 10–100%，不要用单一修正因子。

### 5.4 几何与锚定：硫醇盐、硫醇、硫醚不是一回事

- **BDT 的三种锚定**：裸硫醇盐（S–Au，PBE 0.24；GW 0.83；ASIC 0.06）、H 饱和硫醇（BDT+H，PBE 0.054；GW 0.010）、Au 加合物（SR–Au(I)–SR）。实验 0.011 更像后两者；Strange 2011 明确说“裸硫醇盐不能解释实验”。smj `bdt` 是去质子硫醇盐（`strip_thiol_h: true`），所以它与 Strange 2008 的 0.24 同锚定，但与“实验 0.011”不同结构。
- **Au–S 距离**：smj bdt 2.253 Å、az 2.2511/2.2549 Å、2cn 2.4248 Å；Strange 2008 BDT 2.45 Å；Toher 2008 d_SS=1.9 Å（Au–S≈2.53 Å）；Strange 2011 BDT+H 2.47 Å。smj 的 Au–S 更短，耦合更强；bdt 的 +16% 可能部分来自这个几何差异，而不是 PBE 本身。
- **2CN 是硫醚（S–CH3）**，不是硫醇盐；2.4248 Å 是配位接触，弱耦合。Gulyaev 的 SMePh 也是硫醚，但分子是苯环；Markussen 的 BDCN 是腈基锚（Au–NC），不是硫醚。**2CN 不能拿 BDCN 当同锚定对照**。
- **分子长度/构象**：Li 2007 的 all-trans vs gauche 差 4×；2CN（C30）vs SMePh（C6）长度差很多；azulene 的 2,6 vs 1,3 连接方式不同。长度/构象差异可以轻易盖过泛函差异。
- **规则**：先对齐“锚定化学 + 接触位点 + Au–S 距离”，再谈泛函；不同锚定的 G 不可直接比。

### 5.5 温度、展宽与 G 的定义

- **理论 T=0 vs 实验 300 K**：文献的 T(E_F) 多是 T=0、单几何；实验是 300 K、多几何平均，且有 10–50% 的半宽（BDA 半宽 47±8%）。smj 的 L0 是 300 K 积分，T(E_F) 是最近格点；两者相差 2.3–8.9%（见 Table 3）。与文献 T=0 的 T(E_F) 比，用 smj 的 T(E_F)；与实验 300 K 的 G 比，用 L0。
- **smj 的 σ_near/σ_far=0.05/0.10 Ha**：这是电极原子自能虚部的**模型常数**（最近 3 个电极原子 0.05 Ha，其余 6 个 0.10 Ha），不是收敛参数；0.05 Ha≈1.36 eV、0.10 Ha≈2.72 eV，远大于文献中常用的 meV 量级 η。它会把电极 DOS/透射峰展宽，对“峰尾主导”的 G 影响尤其大。文献的 T(E_F) 通常用很小的 η，**绝对 G 不能无条件对比**；要么用同样的展宽重算，要么只比量级。
- **G 的定义**：文献的“G”有的是零偏压微分电导，有的是 T(E_F)，有的是有限偏压 I–V 拟合；本表统一按“零偏压 G/G0 = T(E_F)”理解，例外在备注中说明（如 Li 2007 的 Gc/β 是系列拟合）。
- **规则**：报 G 时必须同时报温度、展宽/η、电极模型、G 的定义；缺一项，比较的可信度就降一档。

### 5.6 与 smj 0.2.0 对齐的具体规则

1. **同泛函**：smj 是 PBE/light，优先与 PBE/GGA 行比（Strange 2008 BDT 0.24；Quek 2007 BDA 0.046；Darancet 2012 的 PBE；Toher 2008 的 LDA 0.23）。不要拿 GW 0.83、DFT+Σ 5.7×10⁻³、ASIC 0.06 当成“PBE 应该等于”的靶；它们只用来判断误差方向。
2. **同锚定**：`bdt`/`az` 是硫醇盐（去质子），与 Strange 2008/2011 的 BDT 同锚定；`2cn` 是硫醚，与 Gulyaev 2026 的 SMePh 同锚定类，但与 Markussen 的 BDCN（腈基）不同。**不要跨锚定比绝对 G**。
3. **同电极模型**：smj 是 37 Au/3 层团簇，文献多是 3–8 层周期板或 63 Au 团簇。可用 Gulyaev（FHI-aims + AITRANSS、Au21 金字塔）作为“同代码链”的桥梁，但它仍不是同电极；团簇/层数差异通常 10–30%，没有通用修正因子。
4. **同能量点**：与文献 T=0 的 T(E_F) 比，用 smj 的 T(E_F) 最近点（bdt 0.2726、az 0.3495、2cn 0.00125）；与 300 K 实验 G 比，用 L0（0.2789、0.3806、0.001325）。不要拿 L0 去比文献 T(E_F) 而不说明 2–9% 的差异。
5. **同几何**：smj 的 Au–S 2.253 Å（bdt）比 Strange 的 2.45 Å 短，耦合更强；bdt 的 +16% 可能来自几何而非泛函。比较时记录 Au–S 距离、接触位点、是否去质子。
6. **同展宽**：smj 的 σ=0.05/0.10 Ha 是模型常数；文献 η 多为 meV。三套的 G 都由 E_F 附近共振峰的尾巴主导（BDT/AZ 峰在 HOMO 侧，2CN 在 LUMO 侧，见 `t_ef_resonance_attribution`），展宽会显著改变 G；要么用相同 σ 重算文献几何，要么只比量级。
7. **同分子**：`2cn` 与 SMePh 不是同分子（C30 vs C6）；`az` 与 Xia 的 2,6Az 不是同锚定（硫醇盐 vs thiochroman）。这两行只能做“量级自洽”判断，不能做“定量验证”。
8. **结论**：`bdt` 可以做“同锚定 + 同泛函 + 同结构类”的定量对标（0.2789 vs 0.24，+16%，但峰位不同）；`az`/`2cn` 目前只能做“预测 + 量级检查”，要变成“验证”，必须补同分子/同锚定的文献值或自己做 periodic/GW 对照。

### 5.7 逐体系结论（给 smj 用户）

- **bdt（0.2789）**：与 Strange 2008 的 PBE 0.24 同锚定、同泛函、同结构类，差 +16%，在团簇/周期、3 层/3–4 层、light/DZP、Au–S 2.253/2.45 Å 的合理差异范围内。**这是三套里唯一可以写“与文献一致（量级/16%）”的**。但 smj 的 T(E) 在 −0.339/−0.687 eV 有两个峰、E_F 尾由 −0.339 eV 峰主导，Strange 的共振更靠近 E_F；如果审稿人问“为什么 G 对但峰位不对”，要如实说明。与实验 0.011 差 25×，是 PBE 裸硫醇盐的典型高估；GW 0.83 更差；BDT+H 0.010 才接近实验——所以 smj 的 bdt 是“干净硫醇盐模型”，不是“实验结构复现”。
- **az（0.3806）**：无同分子绝对理论值；Xia 的 3.2×10⁻⁴ 是 thiochroman 连接的实验值，不能直接比。smj 的 0.3806 是“直接硫醇盐 + 强耦合”的预测，T(E_F) 是更近的 −0.177 eV 峰的尾巴（max T 在 −0.406 eV）。**在找到同分子计算值之前，只能标“预测/未验证”**；不要用 1190× 的比值证明对或错。
- **2cn（0.001325）**：硫醚锚、弱耦合、C30 长度；与 Gulyaev 的 SMePh 1.3×10⁻² 同锚定类但不同长度，低 10× 与“更长/更弱”的预期一致；与 BDCN 4.7×10⁻⁴ 的 2.8× 差异主要来自锚定（硫醚 vs 腈基）。**量级自洽，但没有同分子定量锚点**；若要验证，需找 2CN-656-4-thiomethylphenyl 的同分子计算/实验，或做长度序列（SMePh→2CN）对照。

### 5.8 什么时候必须重算，而不是引用文献

- **要绝对 G 误差 <2×**：胺/腈基/吡啶锚必须上 GW/DFT+Σ；纯 PBE 的 3–226× 误差不可接受。
- **要解释 BDT 实验 0.011**：裸硫醇盐 PBE/GW 都不够（0.24/0.83）；必须显式建 H 饱和或 Au 加合物结构（BDT+H GW 0.010）。
- **要把 smj 与文献逐字对齐**：用 smj 的电极/基组/展宽重算文献几何，或用文献的周期板/DZP 重算 smj 几何；否则只能比量级。
- **要给 az/2cn 找锚点**：目前没有同分子绝对计算值；必须自己算同分子 periodic/GW 对照，或找到同分子实验（注意锚定/连接子一致）。

## 6. 数据核验与限制

- **Table 1 全部为 [全文]**：Strange 2008、Strange 2011、Quek 2007、Toher 2008、Bilan 2012、Li 2007、Xia 2014 读的是 arXiv/PubMed Central 全文；Darancet 2012、Markussen 2013、Gulyaev 2026 读的是论文全文/PMC；数值来自正文、Table 1/2 或图注。
- **没有 [摘要级]/[未核验] 数值进入 Table 1**。被排除的：El-Nahas 2016/2017（摘要无 G、正文付费）、Guan 2022（摘要无 G、付费）、Manrique 2015（全文读取失败）、Dutta 2008（只有 I–V 图、无绝对 G）、Thygesen & Jacobsen 2005（k 点方法论文，无绝对 G）。
- **Brandbyge 2002** 只作为 transiesta/NEGF 方法源头引用，不含分子结 G；**Pauly et al. 2008** 作为 Bilan 的有限团簇 NEGF 实现引用。
- **Gulyaev 2026 的 1.3×10⁻² 是“PBE 透射 + PBE0 修正”的拟合值**，不是纯 PBE；**Li 2007 的 Gc/β 是系列拟合**，不是单一结 G；**Xia 2014 是实验**。三者在表中均已标注，不与“纯 PBE 单结计算”混列。
- **实验值本身有分布**：BDT 文献 10⁻⁴–10⁻¹、BDA 半宽 47%、C60 室温/10 K 差 200×；不要把某一个实验值当成精确靶。
- **smj 三套是验证样本**：按 2026-11-15 口径物理判定 `passed`（近 E_F 峰已按族归因）；Table 3 的数值是“本工作/未发表”，不是生产结论，也不代表方法已经验证通过。
- **没有编数**：所有 G 值、比值、方法参数均可回溯到上述 DOI/locator；比值标“≈”，是本文按原文数值做的算术，不是原文直接给出的（原文直接给出的比值在备注中注明，如 Quek 的 7×、Markussen 的 226×/0.75×/5.6×）。

## 7. 参考文献

1. Brandbyge, M.; Mozos, J.-L.; Ordejón, P.; Taylor, J.; Stokbro, K. Density-functional method for nonequilibrium electron transport. *Phys. Rev. B* **2002**, *65*, 165401. DOI: 10.1103/PhysRevB.65.165401.
2. Strange, M.; Kristensen, I. S.; Thygesen, K. S.; Jacobsen, K. W. Benchmark density functional theory calculations for nanoscale conductance. *J. Chem. Phys.* **2008**, *128*, 114714. DOI: 10.1063/1.2839275.
3. Strange, M.; Rostgaard, C.; Häkkinen, H.; Thygesen, K. S. Self-consistent GW calculations of electronic transport in thiol- and amine-linked molecular junctions. *Phys. Rev. B* **2011**, *83*, 115108. DOI: 10.1103/PhysRevB.83.115108.
4. Quek, S. Y.; Venkataraman, L.; Choi, H. J.; Louie, S. G.; Hybertsen, M. S.; Neaton, J. B. Amine−Gold Linked Single-Molecule Circuits: Experiment and Theory. *Nano Lett.* **2007**, *7*, 3477–3482. DOI: 10.1021/nl072058i.
5. Toher, C.; Sanvito, S. Effects of self-interaction corrections on the transport properties of phenyl-based molecular junctions. *Phys. Rev. B* **2008**, *77*, 155402. DOI: 10.1103/PhysRevB.77.155402.
6. Darancet, P.; Widawsky, J. R.; Choi, H. J.; Venkataraman, L.; Neaton, J. B. Quantitative Current–Voltage Characteristics in Molecular Junctions from First Principles. *Nano Lett.* **2012**, *12*, 6250–6254. DOI: 10.1021/nl3033137.
7. Markussen, T.; Jin, C.; Thygesen, K. S. Quantitatively accurate calculations of conductance and thermopower of molecular junctions. *Phys. Status Solidi B* **2013**, *250*, 2394–2402. DOI: 10.1002/pssb.201349217.
8. Gulyaev, A.; Hazarika, J.; Liu, Z.-F.; Venkataraman, L. A Computationally Efficient and Accurate Method for Predicting Conductance of Single-Molecule Junctions. *Nano Lett.* **2026**, *26*, 7429–7434. DOI: 10.1021/acs.nanolett.6c01462.
9. Bilan, S.; Zotti, L. A.; Pauly, F.; Cuevas, J. C. Theoretical study of the charge transport through C60-based single-molecule junctions. *Phys. Rev. B* **2012**, *85*, 205403. DOI: 10.1103/PhysRevB.85.205403.
10. Li, C.; Pobelov, I.; Wandlowski, T.; Bagrets, A.; Arnold, A.; Evers, F. Charge Transport in Single Au | Alkanedithiol | Au Junctions: Coordination Geometries and Conformational Degrees of Freedom. *J. Am. Chem. Soc.* **2007**, *130*, 318–326. DOI: 10.1021/ja0762386.
11. Xia, J.; Capozzi, B.; Wei, S.; Strange, M.; Batra, A.; Moreno, J. R.; Amir, R. J.; Amir, E.; Solomon, G. C.; Venkataraman, L.; Campos, L. M. Breakdown of Interference Rules in Azulene, a Nonalternant Hydrocarbon. *Nano Lett.* **2014**, *14*, 2941–2945. DOI: 10.1021/nl5010702.
12. Pauly, F.; Viljas, J. K.; Huniar, U.; Häfner, M.; Wohlthat, S.; Bürkle, M.; Cuevas, J. C.; Schön, G. Cluster-based density-functional approach to quantum transport through molecular and atomic contacts. *New J. Phys.* **2008**, *10*, 125019. DOI: 10.1088/1367-2630/10/12/125019.
13. Evers, F.; Korytár, R.; Tewari, S.; van Ruitenbeek, J. M. Advances and challenges in single-molecule electron transport. *Rev. Mod. Phys.* **2020**, *92*, 035001. DOI: 10.1103/RevModPhys.92.035001.
14. Thygesen, K. S.; Jacobsen, K. W. Interference and k-point sampling in the supercell approach to phase-coherent transport. *Phys. Rev. B* **2005**, *72*, 033401. DOI: 10.1103/PhysRevB.72.033401.

## 附录 A：换算与公式

- G0 = 2e²/h = 77.4809 μS（CODATA）。
- 1 Ha = 27.211386 eV；σ=0.05 Ha=1.3606 eV，0.10 Ha=2.7211 eV。
- 300 K 的热能 k_BT = 0.025852 eV。
- Landauer 零偏压：G = G0·T(E_F)；有限温度：G = G0·∫T(E)(−∂f/∂E)dE = G0·L0。
- 有效数字：Table 1 保留原文；Table 3 保留 smj 的 9 位小数（L0）与 6 位小数（T(E_F)）。

## 附录 B：与 smj 结果对齐的检查清单

在把 smj 的 G 与任何文献值并列之前，逐项打勾：

- [ ] 泛函相同（PBE vs PBE；不是 PBE vs GW/DFT+Σ/SIC）。
- [ ] 锚定相同（硫醇盐 vs 硫醇盐；硫醚 vs 硫醚；不是硫醇盐 vs 硫醚/腈基）。
- [ ] 电极模型相同/可解释（团簇 vs 周期；层数；尖端/平面；接触位点）。
- [ ] k∥/团簇尺寸收敛（Γ-only 会差 50%；37 Au vs 63 Au 会差）。
- [ ] 能量点相同（T(E_F) vs L0；T=0 vs 300 K）。
- [ ] 展宽相同/可解释（smj σ=0.05/0.10 Ha vs 文献 η=meV）。
- [ ] 几何相同（Au–S 距离、是否去质子、构象）。
- [ ] G 的定义相同（零偏压 T(E_F) vs 有限偏压微分电导 vs 系列拟合）。
- [ ] 实验值有分布（半宽/多几何），不是精确靶。

九项中任何一项不满足，比较结论最多写到“同量级”；不要写“验证通过”。
