# 单分子结电导实验 benchmark 表（一手文献）

> 范围：单分子结/单分子器件在 Au（及少量 Cu）电极上的零偏压/低偏压电导实验值，用于 AutoSMJ 的外部有效性（L4）对照。  
> 证据分级：**【一手全文】** = 已读论文全文/预印本并核对图表；**【一手摘要】** = 仅读摘要/题录，未读全文；**【换算】** = 由原文电阻/电导换算，公式见 §1；**【项目】** = 本仓库 `runs/` 的 DFT 值，仅用于 §4 对照，**不进入 benchmark 表**。  
> 检索源：arXiv、OpenAlex、Crossref、Europe PMC/PMC、DOI registry；检索日期 2026-11-15。  
> 关联文档：[SMJ 科学性要求与判定标准调研](smj-科学性要求与判定标准调研.md)（其 L4 明确指出“仅 BDT 理论带，无实验对照”）、[计算流程梳理](smj-计算流程梳理.md)。

---

## 0. 结论速览

- 本表收录 **9 个分子/体系、9 类结构**：芳香二硫醇、饱和烷烃二硫醇、吡啶锚定、胺锚定、非交替芳烃/薁、亚胺桥连、富勒烯、OPE、异腈/碳锚定。
- 推荐对照主值（G0）：BDT **0.011–0.020**；ODT **2.53×10^-4**；BPY **9.93×10^-3**；BDA **6.4×10^-3**；2,6Az **3.2×10^-4**（thiochroman 锚定，非自由二硫醇）；2CN **无同分子值**（Miao 的 5.3×10^-4 属双亚胺分子，≠ 本项目二腈）；C60 **1×10^-2（Au dumbbell）/0.25（Cu）**；OPE3 **2.8×10^-4**；DICB **3×10^-3**。
- **关键更正 1**：Xia 2014 的 “2,6Az” 是 **thiochroman 环内硫醚锚定**的薁，不是 “2,6-薁二硫醇”（不是自由 –SH）。
- **关键更正 2**：Miao 2022 的 “2CN” 是 **双亚胺（C=N）桥 + SMe 端基** 的另一分子，不是本项目 `2cn`（C30H20N2S2，两个 terminal –C≡N，1.147 Å，SMe 硫醚锚）。两者只是缩写撞名；5.3×10^-4 不能当本项目 2CN 的基准。
- 证据强度：BDA、2,6Az、Miao 的 2CN、OPE3、DICB、Leary 的 BDT 及 Néel 的 C60 已读全文；Xiao 的 BDT、Xu & Tao 的 ODT/BPY、Gillemot 的 C60 为摘要级；所有数值均为一手论文，未使用综述二手数字。
- 与项目 DFT 的关系：`runs/` 的 BDT 0.2789 G0、AZ 0.3806 G0 等是**计算值**，不能与本表实验值混用；数量级对照见 §4，仅作诊断提示。

---

## 1. 口径与换算

### 1.1 G0 的定义

电导量子取自旋简并的

\[
G_0=\frac{2e^2}{h}=7.748091729\times10^{-5}\ \mathrm{S},
\qquad
R_0=\frac{1}{G_0}=12.906\ \mathrm{k\Omega}.
\]

- 本表 G（G0）均以 `G/G0` 表示；G（S）由 `G/G0 × G0` 换算得到，保留 2–3 位有效数字。
- 不要把 `G0` 与 `e^2/h = 3.874×10^-5 S` 混淆；本表统一用 `2e^2/h`。
- 原文给电阻时，先 `G(S)=1/R`，再除以 `G0`；此类数值在表中标“由 … MΩ 换算”。

### 1.2 数值选取与标注规则

1. **优先一手实验论文**，且优先能读到全文/预印本、能定位到图表的值；只读到摘要的标“摘要级”。
2. **同一分子有多个可信值时，给主值 + 备注中的备选值**，不把不同接触构型的值平均成一个“真值”。
3. **不确定度/分布宽度按原文**：有的给 ±，有的给 σ（log-normal），有的只给“峰宽未给数值”；不替原文补误差棒。
4. **单结 vs 多结**：优先选明确声明“单分子/单结统计”的峰；crack-defined break junction 可能形成并联多结的，只作交叉参考，不作主值。
5. **换算值**保留原文有效数字；例如 51±5 MΩ 换算为 2.53×10^-4 G0，但原文只有 2 位有效数字（51），正文同时给出换算过程。

---

## 2. Benchmark 表

| # | 分子（中/英/缩写） | 结构类型 | 锚定基/电极 | G（G0） | G（S） | 不确定度/分布 | 方法 | 温度/环境 | 一手文献 DOI + locator | 备注 |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 苯-1,4-二硫醇<br>benzene-1,4-dithiol<br>BDT | 芳香二硫醇（交替芳烃，苯环） | –SH；Au（Au–S） | **0.020**（真 Au–S–Ph–S–Au 结；Leary 2018）<br>**0.011**（经典值；Xiao 2004） | 1.55×10^-6（0.020）<br>8.52×10^-7（0.011） | Leary：log(G/G0) = −1.7 ± 0.4（≈×2.5）；Xiao：峰宽未给数值；Kiguchi：0.004 ± 0.001 G0 | BJ（Leary；Vbias = 0.2 V）<br>STM-BJ（Xiao；摘要级）<br>STM-BJ（Kiguchi；Au tip/Au(111)，20 mV） | Leary：温度/环境未逐字核（正文讨论室温 Au 接触）；Xiao：未核；Kiguchi：室温、tetraglyme 1 mM | Leary 2018, J. Phys. Chem. C 122(6), 3211–3218, DOI 10.1021/acs.jpcc.7b11104；arXiv:1803.07411，全文/预印本，§Discussion（group 1）<br>Xiao 2004, Nano Lett. 4(2), 267–271, DOI 10.1021/nl035000m；摘要<br>Kiguchi 2006, Appl. Phys. Lett. 89(21), 203103, DOI 10.1063/1.2392816；Fig. 2e/f+正文 | 单结（Leary group 1 = 真结；group 2 = Au(SR)2 寡聚接触，约 3×10^-4 G0）；BDT 文献跨度约 4×10^-4–2×10^-2 G0，接触构型是主因；不要用单一 0.011 当“真值” |
| 2 | 辛二硫醇<br>octane-1,8-dithiol<br>ODT（C8） | 饱和烷烃二硫醇 | –SH；Au（Au–S） | **2.53×10^-4**（由 51 ± 5 MΩ 换算） | 1.96×10^-8 | ±5 MΩ → ±9.8%（G ±2.5×10^-5 G0；S ±1.9×10^-9） | STM-BJ（Tao 组；摘要级） | 未核（摘要未给；全文付费） | Xu & Tao 2003, Science 301(5637), 1221–1223, DOI 10.1126/science.1087481；摘要 | 单结；histogram 整数倍峰；零偏压电阻。同文 C6 10.5 ± 0.5 MΩ → 1.23×10^-3 G0，C10 630 ± 50 MΩ → 2.05×10^-5 G0；βN = 1.0 ± 0.1/C；G 为本表换算，原文给 R |
| 3 | 4,4'-联吡啶<br>4,4'-bipyridine<br>BPY | 吡啶锚定（联吡啶） | 吡啶 N；Au（Au–N） | **9.93×10^-3**（由 1.3 ± 0.1 MΩ 换算） | 7.69×10^-7 | ±0.1 MΩ → ±7.7%（G ±0.76×10^-3 G0；S ±0.59×10^-7） | STM-BJ（Tao 组；摘要级） | 未核（摘要未给；全文付费） | Xu & Tao 2003, Science 301(5637), 1221–1223, DOI 10.1126/science.1087481；摘要 | 单结；零偏压；pyridine N 锚定（Au–N），不是 –SH；与烷烃二硫醇同文同统计口径 |
| 4 | 1,4-苯二胺<br>benzene-1,4-diamine<br>BDA（1,4-diaminobenzene） | 胺锚定（芳香二胺） | –NH2；Au（Au–N） | **0.0064 ± 0.0004** | 4.96×10^-7 ± 0.31×10^-7 | ±0.0004 G0（Gaussian fit） | STM-BJ（Au point contacts；25 mV；>3000 traces） | 室温、ambient、TCB（1 mM） | Venkataraman 2006, Nano Lett. 6(3), 458–462, DOI 10.1021/nl052373+；arXiv:cond-mat/0603281，Fig. 2b inset，p. 4 | 单结；Gaussian fit to 1D histogram peak；amine 系列 β = 0.91 ± 0.03/CH2；BDA 是 amine 锚定的常用基准 |
| 5 | 2,6-薁<br>2,6-azulene<br>2,6Az（thiochroman 锚定） | 非交替芳烃/薁（5-7 稠环） | thiochroman 环内硫醚 S；Au（Au–S） | **3.2×10^-4**（32×10^-5） | 2.48×10^-8 | Table 1 平均电导，未给数值误差；直方图峰宽见 Fig. 3b | STM-BJ（Au；TCB；250 mV；5000 traces） | ambient、TCB 溶液 | Xia 2014, Nano Lett. 14(5), 2941–2945, DOI 10.1021/nl5010702；arXiv:1411.5463，Table 1，Fig. 3 caption（p. 15） | 单结；Table 1 为 2D histogram 平均电导；原文为 thiochroman-linked azulene，S 在 thiochroman 环内，不是自由 –SH；不要叫“2,6-薁二硫醇”。同表 1,3Az 3.2×10^-4、4,7Az 8×10^-5、5,7Az 2×10^-5、1,3PyAz 9×10^-5 |
| 6 | 2CN（Miao 2022）<br>2CN（双亚胺桥连、SMe 端基；≠ 本项目二腈） | 亚胺（C=N）桥连；硫醚锚定（此分子非腈基） | –SMe；Au（Au–S） | **5.3×10^-4** | 4.11×10^-8 | Gaussian 拟合，未给数值误差 | STM-BJ（Au；TCB；air；100 mV；5000–10,000 traces） | ambient/air、TCB 溶液 | Miao 2022, Nano Lett. 22(20), 8331–8338, DOI 10.1021/acs.nanolett.2c03288；PMC9614958，Fig. 2c+正文，Fig. 3b caption | 单结；“2” = two imine C=N bridge atoms，不是两个 CN；端基 SMe；对照 CN 4.3×10^-4、BN 2.4×10^-4、BO 8.7×10^-5、2CN-m 9.7×10^-5；与本项目 `2cn`（二腈取代 + SMe 锚，C30H20N2S2）不同分子，不作其基准 |
| 7 | C60 富勒烯<br>C60 fullerene<br>C60（dumbbell；对照 C60/Cu(100)） | 富勒烯/范德华接触（无传统锚定基） | 无；Au–C60（dumbbell）；对照 Cu–C60/Cu(100) | **~1×10^-2**（dumbbell/Au；Gillemot 2013）<br>**0.25**（C60/Cu(100)；Néel 2007） | 7.75×10^-7（1×10^-2）<br>1.94×10^-5（0.25） | Gillemot：“close to 10^-2 G0”，未给数值；Néel：~0.25，未给数值 | ambient STM（dumbbell；摘要级）<br>LT-STM 8 K（C60/Cu(100)） | ambient（dumbbell）；8 K、UHV（Néel） | Gillemot 2013, Small 9(22), 3812–3822, DOI 10.1002/smll.201300310；摘要<br>Néel 2007, Phys. Rev. Lett. 98, 065502, DOI 10.1103/PhysRevLett.98.065502；arXiv:cond-mat/0608476，摘要+正文 | 单分子；C60 无传统锚定基，接触几何/电极材料影响大；dumbbell/Au 值更接近 Au–S 体系；0.25 是 Cu(100) 上 tip–C60 接触，不宜直接与 Au–S 比 |
| 8 | OPE3<br>oligo(phenylene ethynylene)3<br>OPE3-SAc | OPE（共轭刚性棒） | –SAc（TBAH 原位脱保护 → –S–Au）；Au | **2.8×10^-4**（log-normal 众数） | 2.17×10^-8 | σ(ln G) = 0.84（几何 SD ≈×2.3）；>2000 traces | MCBJ（Au；空气；0.1 V；>2000 traces） | 室温、空气 | Frisenda 2015, Beilstein J. Nanotechnol. 6, 1558–1567, DOI 10.3762/bjnano.6.159；PMC4578406，Table 1，§Results | 单结；log-normal 拟合峰，mode = exp(μ−σ²)；同表 SMe 5.9×10^-5、Py 2.2×10^-6、NH2 7.0×10^-5；DuBois 2018 crack-defined BJ 1–3×10^-4 交叉一致但可能多结，未作主值 |
| 9 | 1,4-二异腈苯<br>1,4-diisocyanobenzene<br>DICB | 异腈（isonitrile）/碳锚定芳烃 | –N≡C（C 端与 Au 成键；Au–C）；Au | **3×10^-3 ± 1×10^-3** | 2.32×10^-7 ± 0.77×10^-7 | ±1×10^-3 G0（直方图峰） | STM-BJ（Au tip/Au(111)；20 mV；>1000 traces） | 室温、tetraglyme（1 mM；toluene 对照） | Kiguchi 2006, Appl. Phys. Lett. 89(21), 203103, DOI 10.1063/1.2392816；arXiv:cond-mat/0610388，Fig. 2a/b+正文 | 单结；同文 BDT 4×10^-3 ± 1×10^-3；1,4-二氰苯（Au–N≡C，N 锚定）无明确峰，G < 1×10^-3；真腈基没有可靠单结值，故本表用异腈/碳锚定代表“碳–Au”类 |

### 2.1 表注

- **主值加粗**；同一行给出两个值时，按“更接近理想单结/更常用”排序，并在备注说明差异。G（S）是换算值，原文的 G0 值是权威值。
- **“由 … MΩ 换算”**：Xu & Tao 2003 原文给的是零偏压电阻，不是 G0；换算过程见 §1.1，表中同时保留原文电阻与换算后的 G0。
- **“摘要级”**：只读到摘要/题录，未读全文；这类值可用，但方法细节、温度/环境、误差口径可能未核。
- **BDT 行**：Leary 的 0.020 是“真 Au–S–Ph–S–Au 结”的窄分布中心，更接近理想二硫醇结；Xiao 的 0.011 是经典 histogram 基峰，文献引用最广；两者差约 2 倍，落在 BDT 的接触构型分散内。
- **2,6Az 行**：Xia 的 Table 1 是 **2D histogram 的平均电导**，不是 1D 峰的众数；1D 峰见 Fig. 3b。两者不要混用。
- **C60 行**：Au–C60 dumbbell 与 Cu–C60/Cu(100) 是两种完全不同的接触/电极，不能把 0.25 G0 当成“C60 在 Au 上的电导”。

---

## 3. 关键更正与易混点

### 3.1 2,6Az 不是“2,6-薁二硫醇”

Xia et al. 2014 的四个分子是 **thiochroman-linked azulene derivatives**：硫原子在 thiochroman（苯并噻吡喃）环内，通过环内 S 与 Au 成键，不是两个自由 –SH。原文 Table 1 的 2,6Az 是 32×10^-5 G0；1,3Az 也是 32×10^-5 G0；4,7Az 8×10^-5；5,7Az 2×10^-5；1,3PyAz（吡啶端基）9×10^-5。把 2,6Az 叫成“2,6-薁二硫醇”会误导锚定化学与质子化状态。项目文档中的 “AZ” 若指另一个薁衍生物，不能仅凭“薁”这个名字套用 3.2×10^-4 G0。

### 3.2 2CN：同名不同物，Miao 的值不能当本项目基准

Miao et al. 2022 的 2CN 是 **两个 imine（C=N）桥连**的分子；“2” 表示 bridge atoms 的数目（two imine bridges），不是两个氰基。端基是 **aurophilic thioether（–SMe）**，与 Au 形成 Au–S 键。原文：G(2CN) = 5.3×10^-4 G0；G(CN) = 4.3×10^-4；G(BN) = 2.4×10^-4；G(BO) = 8.7×10^-5；2CN-m = 9.7×10^-5。原文用 STM-BJ、TCB、空气、100 mV、10,000 traces，Gaussian fit 到 1D histogram 峰。

**本项目的 `2cn` 是另一个分子**：`demo/2CN-656-4-thiomethylphenyl/molecule.xyz` 为 C30H20N2S2、54 原子，两个 **terminal 腈基 –C≡N（1.147 Å；N 只连一个 C）**，端基为 **–S–CH3（硫醚）** 锚定。项目里 “2CN” 指“二腈取代 + 甲硫基锚”的共轭分子，不是 Miao 的双亚胺桥分子，两者只是缩写撞名。**因此 Miao 的 5.3×10^-4 不能作为本项目 2CN 的实验基准**，§4 对照里这一行留空（或标“不同分子”）。

真腈基（1,4-二氰苯，Au–N≡C）在 Kiguchi 2006 中没有明确电导峰（G < 1×10^-3 G0）；同一文的 1,4-二异腈苯（Au–C≡N，C 锚定）给出 3×10^-3 ± 1×10^-3 G0。因此本表用 **异腈/碳锚定** 代表“碳–Au”类，而不伪造一个“腈基 BDT”值；本项目 2CN 的锚是 SMe，C≡N 只是取代基，不参与 Au 成键。

### 3.3 BDT 的 0.011 不是唯一值

BDT 是单分子结里最“著名”也最分散的 benchmark：

| 来源 | G（G0） | 条件/口径 |
|---|---|---|
| Xiao 2004 | 0.011 | 经典 histogram 基峰；摘要级 |
| Leary 2018 | 0.020 ± 0.4 dec | 真 Au–S–Ph–S–Au 结，Vbias = 0.2 V；全文/预印本 |
| Kiguchi 2006 | 0.004 ± 0.001 | 同文 STM-BJ 对照，20 mV；全文/预印本 |
| Tsutsui 2006 | ~0.01 与 ~0.1 两组态 | MCBJ；高电导态稳定到 1.6 V；摘要级 |

Kiguchi 2006 还转引了 BDT 文献范围约 **4.1×10^-4–1.1×10^-2 G0**；Leary 的 0.020 略高于该范围，是因为它专门分离出“真结”而不是低电导的 Au(SR)2 寡聚接触。**结论：BDT 应给范围 0.004–0.020 G0，并按接触构型说明；不要用单一 0.011 作为“真值”。**

### 3.4 C60 的电极依赖

- **Au–C60 dumbbell**：Gillemot 2013，ambient STM，单分子先成像确认，高电导态“close to 10^-2 G0”；这是 Au 电极、有分子桥的体系。
- **C60/Cu(100)**：Néel 2007，LT-STM，8 K，UHV，Cu(100)，W tip（Cu-covered）；从隧穿进入接触时 G ≈ 0.25 G0，继续压至 G ≈ G0 跳变；C60 无传统锚定基，靠范德华/轨道杂化接触。

两者电极、温度、接触几何都不同，不能直接比较；若项目用 Au 电极，优先参考 dumbbell 的 1×10^-2 G0，并把 0.25 G0 仅作为“C60–Cu 接触”的对照。

---

## 4. 与项目 `runs/` 的对照（仅数量级提示）

> **重要**：下表左列是项目 `runs/` 的 **DFT 值**，不是实验值；不同泛函、电极模型、接触几何、分子结构都可能不同。比值只用于判断“是否差一个数量级以上”，不是误差评估，也不构成对计算代码的否证。

| 项目分子 | `runs/` DFT（G0） | 本表实验基准（G0） | DFT/exp | 备注 |
|---|---:|---:|---:|---|
| BDT | 0.2789 | 0.011–0.020 | 14–25× | DFT 高估；理想 Au–S–Ph–S–Au 几何/接触与实验“真结”仍有差距 |
| HDT（C6） | 1.578×10^-3 | 1.229×10^-3 | 1.28× | 同系物，吻合较好 |
| ODT（C8） | 2.936×10^-4 | 2.531×10^-4 | 1.16× | 同系物，吻合较好 |
| DDT（C10） | 6.095×10^-5 | 2.049×10^-5 | 2.97× | 可接受，但已见偏差 |
| BPY | 1.159×10^-3 | 9.928×10^-3 | 0.117× | DFT 低 8.6×；pyridine N 锚定/能级对齐需检查 |
| 2CN | 1.325×10^-3 | 无同分子值 | — | Miao 的 2CN 是双亚胺桥分子（5.3×10^-4），与本项目二腈不同，不能相比；只能与同锚定 SMe 的长度序列对照 |
| AZ | 0.3806 | 3.2×10^-4 | ~1.2×10^3 | 结构/锚定未必一致；仅提示“差三个数量级”，不可直接归因于代码 |
| SMe | 0.05637 | — | — | 无同分子实验值；不要与 OPE3-SMe 的 5.9×10^-5 直接比 |

**读法**：HDT/ODT/DDT 的 1.2–3× 偏差在“接触几何 + 泛函 + 电极模型”的合理分散内；BPY 的 8.6× 低估值和 BDT 的 14–25× 高估值更可能是锚定/能级对齐/接触构型问题；AZ 的 10^3 差异先要确认项目分子与 Xia 的 2,6Az 是否同一结构，不能直接下结论。

---

## 5. 使用注意

### 5.1 方法、电极、锚定与统计口径的可比性

- **方法**：STM-BJ（Au tip/Au(111) 或 Au point contact）、MCBJ（机械可控断裂结）、crack-defined break junction、LT-STM 接触。不同方法的电极形貌、接触历史、拉伸速率不同；同一分子可差数倍。
- **电极**：Au（大多数）、Cu(100)（C60/Néel）。Au 的 (111) 表面、蒸镀膜、tip 原子构型都会改变接触；Cu 与 Au 的功函数/轨道杂化不同，不能混比。
- **锚定**：–SH（Au–S）、–SMe/–SAc（Au–S）、thiochroman 环内 S（Au–S）、–NH2（Au–N）、吡啶 N（Au–N）、–N≡C（Au–C）、C60（范德华/杂化）。锚定的质子化/去保护状态（SAC + TBAH）会改变结合与电导。
- **统计口径**：histogram 基峰、Gaussian 峰中心、log-normal 众数、2D histogram 平均电导、最后平台电导——这些不是同一个估计量。表中已逐行标注；跨行比较时先看“估计量是否一致”。
- **单结 vs 多结**：多数工作用“整数倍峰/基峰”排除多分子结；crack-defined BJ 可能形成并联结；有成像确认的（C60 dumbbell）更可靠。不要把“平均电导”与“单结基峰”直接相除。

### 5.2 G0 定义与单位换算

- `G0 = 2e^2/h = 7.748091729×10^-5 S`；`R0 = 12.906 kΩ`。
- `G(S) = (G/G0) × G0`；`G/G0 = G(S)/G0`。
- 原文给电阻 `R` 时：`G/G0 = R0/R`。例：51 MΩ → `12.906 kΩ / 51 MΩ = 2.53×10^-4`。
- 不要把 `G0` 与 `e^2/h` 或 `h/e^2` 的倒数混淆；本表全部用 `2e^2/h`。
- 表中 G（S）保留 2–3 位有效数字；需要更高精度时用 `G/G0` 乘 `7.748091729×10^-5`。

### 5.3 G 与 T(E_F) 的关系

- 低偏压、零温极限下，Landauer 公式：`G = G0 T(E_F) = G0 Σ_i τ_i(E_F)`。实验的 `G/G0` 是 `T(E_F)` 的**系综/最可几**估计，不是某一个静态几何的 `T(E_F)`。
- 实验峰对应“接触几何 + 分子构象 + 电极形貌”的分布；DFT 单几何的 `T(E_F)` 只是该分布中的一个点。把两者直接相等，会把接触差异误判为泛函误差。
- 远离共振的隧穿区，`G ∝ exp(−βL)`；烷烃二硫醇 `βN ≈ 1.0 ± 0.1/C`（Xu & Tao 2003），胺系列 `β = 0.91 ± 0.03/CH2`（Venkataraman 2006）。这是比单个绝对值更稳健的一致性检查。
- 近共振时，`G` 对 `E_F − ε0` 和耦合 `Γ` 极敏感；PBE 常把分子共振放在离 `E_F` 过近的位置，系统性高估线性响应电导（参见[科学性调研](smj-科学性要求与判定标准调研.md) §3.2）。BDT 的 DFT 0.2789 G0 与实验 0.011–0.020 G0 的差距，应先查接触几何/能级对齐，而不是直接“调参到 0.011”。
- `G/G0` 也不等于“透射谱在 E_F 的单一本征通道”；多通道、非相干效应、电极态都会让 `T(E_F)` 的解释复杂化。

### 5.4 温度、溶剂与偏压

- **温度**：C60/Cu(100) 在 8 K、UHV；大多数溶液结在室温、ambient。热展宽约 `3.5 kBT`：300 K ≈ 25.9 meV，8 K ≈ 0.69 meV。温度会影响能级对齐、峰宽和结稳定性；8 K 的 0.25 G0 不能直接搬到 300 K 的 Au 体系。
- **溶剂**：TCB（BDA、2,6Az、Miao 的 2CN）、tetraglyme（DICB）、二氯甲烷（OPE3 沉积）、空气（OPE3 MCBJ）。溶剂可稳定结、屏蔽电荷、改变 Au–S 结合与分子构象；真空 DFT 与溶液实验比较时要写明环境差异。
- **偏压**：20 mV（Kiguchi）、25 mV（Venkataraman）、100 mV（Miao）、250 mV（Xia）、0.1 V（Frisenda）、0.2 V（Leary）。有限偏压把 `T(E)` 在窗口内平均；越接近零偏压，越接近 `T(E_F)`。跨行比较时先看偏压是否同量级。
- **空气/氧/水**：硫醇易氧化；SAc 需 TBAH 原位脱保护；thiochroman/thioether 通常比自由硫醇稳定。样品的保存与制备历史会影响结果。

### 5.5 单结 vs 多结、峰选取

- **整数倍峰法**（Xiao、Xu & Tao、Kiguchi）：histogram 在 `G0`、`2G0`、`3G0`… 出现峰，取最低的基峰为单结；高次峰可能是 2–3 个分子并联。
- **Gaussian 峰法**（Venkataraman、Miao）：对 1D histogram 的峰做 Gaussian 拟合，取中心 ± 拟合误差；需要说明 bin 大小、拟合区间。
- **log-normal 法**（Frisenda）：对 `ln G` 分布拟合，取 mode `exp(μ−σ²)`；σ 是分布宽度，不是标准误。
- **2D histogram 平均**（Xia）：用 `G–z` 二维直方图确定平台，取平均电导；与 1D 峰可能不同，不能混用。
- **平台法**（Néel、Gillemot）：取接触/拉伸过程中的最后平台；受 tip 形貌和接触历史影响大。
- 报告时至少写清：**估计量、bin 大小、是否剔除空曲线、单结判据、traces 数量**；没有这些，跨论文的 ± 不可比。

### 5.6 用于 AutoSMJ 的 L4 判定

- 本表是**外部实验基准**，不是拟合目标。建议在报告中写“同一分子/锚定/电极/几何下的实验范围 + 本文值 + 偏差倍数 + 可能来源”，而不是只写“通过/未通过”。
- 若 DFT 值落在实验范围内（如 HDT/ODT），可标“数量级一致”；若差 3× 以上（BPY、BDT），标“需解释”，不要直接宣称“正确”。
- 若差 10^3（AZ），先核对项目分子与实验分子是否同一结构/锚定；结构不一致时，这个比值没有意义。
- 对摘要级数值（Xiao BDT、Xu & Tao ODT/BPY、Gillemot C60），在报告中标注“摘要级/未核全文”，不要给出比原文更小的误差棒。
- 本表只补了“实验值”这一半；L4 还需要不确定度合成、方法间散布、适用域声明。建议与[科学性调研](smj-科学性要求与判定标准调研.md) §3.3 的“基准表 + 不确定度块 + 适用域”一起使用。

---

## 6. 参考文献

1. Xiao, B.; Xu, B.; Tao, N. J. Measurement of Single Molecule Conductance: Benzenedithiol and Benzenedimethanethiol. *Nano Lett.* **2004**, *4* (2), 267–271. DOI: [10.1021/nl035000m](https://doi.org/10.1021/nl035000m).
2. Xu, B.; Tao, N. J. Measurement of Single-Molecule Resistance by Repeated Formation of Molecular Junctions. *Science* **2003**, *301* (5637), 1221–1223. DOI: [10.1126/science.1087481](https://doi.org/10.1126/science.1087481).
3. Venkataraman, L.; Klare, J. E.; Tam, I. W.; Nuckolls, C.; Hybertsen, M. S.; Steigerwald, M. L. Single-Molecule Circuits with Well-Defined Molecular Conductance. *Nano Lett.* **2006**, *6* (3), 458–462. DOI: [10.1021/nl052373+](https://doi.org/10.1021/nl052373+).
4. Xia, J.; Capozzi, B.; Wei, S.; Strange, M.; Batra, A.; Moreno, J. R.; Amir, R. J.; Amir, E.; Solomon, G. C.; Venkataraman, L.; Campos, L. M. Breakdown of Interference Rules in Azulene, a Nonalternant Hydrocarbon. *Nano Lett.* **2014**, *14* (5), 2941–2945. DOI: [10.1021/nl5010702](https://doi.org/10.1021/nl5010702).
5. Miao, Z.; Quainoo, T.; Czyszczon-Burton, T. M.; Rotthowe, N.; Parr, J. M.; Liu, Z.-F.; Inkpen, M. S. Charge Transport Across Dynamic Covalent Chemical Bridges. *Nano Lett.* **2022**, *22* (20), 8331–8338. DOI: [10.1021/acs.nanolett.2c03288](https://doi.org/10.1021/acs.nanolett.2c03288).
6. Néel, N.; Kröger, J.; Limot, L.; Frederiksen, T.; Brandbyge, M.; Berndt, R. Controlled Contact to a C60 Molecule. *Phys. Rev. Lett.* **2007**, *98*, 065502. DOI: [10.1103/PhysRevLett.98.065502](https://doi.org/10.1103/PhysRevLett.98.065502).
7. Gillemot, K.; Evangeli, C.; Leary, E.; La Rosa, A.; González, M. T.; Filippone, S.; Grace, I.; Rubio-Bollinger, G.; Ferrer, J.; Martín, N.; Lambert, C. J.; Agraït, N. A Detailed Experimental and Theoretical Study into the Properties of C60 Dumbbell Junctions. *Small* **2013**, *9* (22), 3812–3822. DOI: [10.1002/smll.201300310](https://doi.org/10.1002/smll.201300310).
8. Frisenda, R.; Tarkuç, S.; Galán, E.; Perrin, M. L.; Eelkema, R.; Grozema, F. C.; van der Zant, H. S. J. Electrical Properties and Mechanical Stability of Anchoring Groups for Single-Molecule Electronics. *Beilstein J. Nanotechnol.* **2015**, *6*, 1558–1567. DOI: [10.3762/bjnano.6.159](https://doi.org/10.3762/bjnano.6.159).
9. Kiguchi, M.; Miura, S.; Hara, K.; Sawamura, M.; Murakoshi, K. Conductance of a Single Molecule Anchored by an Isocyanide Substituent to Gold Electrodes. *Appl. Phys. Lett.* **2006**, *89* (21), 203103. DOI: [10.1063/1.2392816](https://doi.org/10.1063/1.2392816).
10. Leary, E.; Zotti, L. A.; Miguel, D.; Márquez, I. R.; Palomino-Ruiz, L.; Cuerva, J. M.; Rubio-Bollinger, G.; González, M. T.; Agraït, N. The Role of Oligomeric Gold–Thiolate Units in Single-Molecule Junctions of Thiol-Anchored Molecules. *J. Phys. Chem. C* **2018**, *122* (6), 3211–3218. DOI: [10.1021/acs.jpcc.7b11104](https://doi.org/10.1021/acs.jpcc.7b11104).
11. Tsutsui, M.; Teramae, Y.; Kurokawa, S.; Sakai, A. High-Conductance States of Single Benzenedithiol Molecules. *Appl. Phys. Lett.* **2006**, *89* (16), 163112. DOI: [10.1063/1.2363995](https://doi.org/10.1063/1.2363995).
12. Dubois, V.; Raja, S. N.; Gehring, P.; Caneva, S.; van der Zant, H. S. J.; Niklaus, F.; Stemme, G. Massively Parallel Fabrication of Crack-Defined Gold Break Junctions Featuring Sub-3 nm Gaps for Molecular Devices. *Nat. Commun.* **2018**, *9*, 3433. DOI: [10.1038/s41467-018-05785-2](https://doi.org/10.1038/s41467-018-05785-2).
