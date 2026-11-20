# 判据来源核验：Thygesen 2005 / Jackson 2009 与"分子族 → 预期共振"

> 目的：核验 `src/smj/junction/criteria.py` 中 `t_ef_magnitude`、`g_in_theory_band`、
> `au_s_bond`、`t_ef_no_peak_near_ef` 四条判据的 source{doi,locator} 是否成立，并为按族
> 登记 `applies_to` 提供一手依据。
>
> 纪律：本文所有数字与英文引文均来自**本次实际读到**的正文/摘要；读不到的写明"仅摘要级"或"未读到"，
> 不以记忆补数。日期 2026-11-15。

---

## 0. 可获取性（先说清证据等级）

| 文献 | 可获取性 | 证据等级 |
|---|---|---|
| Thygesen & Jacobsen 2005, Chem. Phys. **319**, 111–125, DOI 10.1016/j.chemphys.2005.05.032 | 开放版 **arXiv:cond-mat/0501238 全文可读**（本次按节读了 §V.2 全节、§V.1、参考文献表、摘要、结论） | **原文核验（全节）** |
| Jackson et al. 2009, Surf. Sci. **603**, 807–813, DOI 10.1016/j.susc.2009.01.022 | OpenAlex / Semantic Scholar 均无开放全文，出版社付费墙；Warwick WRAP 记录页本次抓取超时 | **仅摘要级**（摘要文本经 WRAP 记录与 infona 记录两处一致） |
| 各分子族一手文献 | 见 §3 逐行标注 | 摘要级 / 全文节选级 |

**页码问题（重要）**：registry 写 `locator='p.12'` / `'p.13'`。开放版是 arXiv 预印本（HTML，无期刊页码），
因此**页码本身无法用开放版核验**；可核验的是"该内容确实在 §V.2，且 §V.2 位于正文后段"。
建议 locator 改写为可核验的形式：`§V.2 (Benzene-dithiol between Au(111) surfaces), Fig. 9/10`。

---

## 1. 核验结论（判断先行）

| # | registry 原写法 | 原文实际内容 | 判定 |
|---|---|---|---|
| 1 | `au_s_bond` band **[2.40, 2.49]** Å，源 Thygesen p.12 `d(Au-S)=2.40 Å` | 原文确有 2.40 Å，但它是**输入几何参数**："in accordance with those used in Ref. 48"，即采用 Stokbro 2003 的值，**不是本文的计算或测量结果** | ⚠ **部分成立**：数值可引，但"源"错位（一手应为 Stokbro 2003）；且 **2.49 无源**（见 #5） |
| 2 | `t_ef_magnitude` 典型 ~0.1 G0 | 原文："the predicted conductance is around **0.1G₀**" | ✅ **成立，可直接引用** |
| 3 | 峰位 −3.0/+2.5 eV 宽峰、−1.5/−1.0 eV 窄峰 | 原文整句给出，四个峰位一字不差 | ✅ **成立，可直接引用** |
| 4 | band **[0.03, 0.40]** G0，源 Thygesen p.13 | 原文**从未给出任何区间**；只给"0.1 G₀"与两个**比例**（比 Ref.48 小 4 倍、比 Ref.15 大 3 倍） | ❌ **不成立（派生量）**：0.1×4=0.40，0.1/3=0.0333→registry 手取 0.03 |
| 5 | `au_s_bond` 上端 2.49 Å | 两个锚点合起来只到 2.42+0.02=**2.44** Å | ❌ **无源** |
| 6 | Jackson `S-Au = 2.42±0.02 Å` | 摘要原文一致（NIXSW + S 1s scanned-energy-mode photoelectron diffraction，Au(111)(√3×√3)R30°–CH₃S，atop 位） | ✅ **成立**，但**仅摘要级**，且是**表面吸附几何、甲硫醇盐**，不是电极–分子–电极结 |

### 对 `t_ef_magnitude` / `g_in_theory_band` 的处置建议

按用户规则"来源不支持 [0.03,0.40] 就给替代一手来源，没有则建议删除"——**本次检索未找到任何
"原文直接给出 Au-BDT-Au 的 T(E_F) 区间"的一手文献**（找到的都是单点值）。故给三档建议，按从严到宽：

1. **（从严，符合现行断言纪律）删除该区间判据**：把 `t_ef_magnitude` / `g_in_theory_band` 的 band 取消，
   只保留"0.1 G₀ 量级"这一条**单点**判据（可在 note 里写"±一档"作提示，但不作闸门）。
2. **（推荐）保留数值、改证据等级、降闸门级别**：把 `source` 从"Thygesen p.13 给出 0.03–0.40"改为
   `evidence: 'derived'`（派生包络），locator 写 `§V.2 + 其 Ref.15/Ref.48 的 ×3/×4 关系`，
   并把 `gate` 从 `error` 降为 `warn`。理由：**派生包络不足以构成"否证 H0"的强证据**，用 error 级
   等于把"复现某个二手区间"当成物理真伪判据。
3. **（不可取）维持现状**：即继续声称"源为该文"，与原文不符——这正是本次核验要消除的过度声明。

另有一条原文事实对判据定位很关键，建议写进 note：

> "In general, the calculations for a BDT molecule covalently bonded to the electrodes finds a
> conductance which is **2-3 orders of magnitude higher than the experimental value**."
> — Thygesen 2005 §V.2 首段

即该文的 0.1 G₀ 是**理论族内部的值**，与 BDT 实验电导相差 2–3 个数量级。用 [0.03, 0.40] 作 error 闸门，
实际检验的是"是否落在 DFT-NEGF 理论族内"，**不是**"是否复现实验"。这一点必须在判据 note 里写明。

---

## 2. 逐条引文（原文英文 + locator）

### 2.1 Thygesen & Jacobsen 2005（arXiv:cond-mat/0501238；期刊版 Chem. Phys. 319, 111–125）

体系与方法：$3\times3$ Au(111) 超胞 + 3 原子金字塔接触，金固定为 fcc（$a=4.18$ Å）；
Kohn–Sham DFT（平面波 + 超软赝势）→ 部分占据 Wannier 函数基组 → Landauer–Büttiker 透射。
泛函名称在本次读到的段落中未出现，**故本文不声称其泛函**（不写 PBE/LDA/RPBE）。

#### (a) 是否原文直接给出 Au–S = 2.40 Å？—— 给出，但**是输入参数**

> "We use the bond lengths $d_{\mathrm{Au-S}}=2.40$ Å, $d_{\mathrm{S-C}}=2.75$ Å,
> $d_{\mathrm{C-C}}=1.40$ Å and $d_{\mathrm{C-H}}=1.08$ Å **in accordance with those used in Ref. 48**."
> — §V.2，正文（对应 Fig. 9(a) 的散射区；期刊版约 p.12）

判定：**2.40 Å 是"采用值"而非该文结论**。Ref. 48 = Stokbro, Taylor, Brandbyge, Mozos, Ordejón,
*Comput. Mater. Sci.* **27**, 151 (2003)。若要把 `au_s_bond` 的一手源立住，应引 Stokbro 2003 的几何，
把 Thygesen 2005 降为"二手采用/一致性佐证"。

#### (b) 是否原文给出 T(E_F) ≈ 0.1 G₀？—— 给出

> "Finally we note, that the **predicted conductance is around $0.1G_0$**. This is in very good
> agreement with the results obtained in Refs. 47, 49, while it is **a factor of 4 smaller than the
> result of Ref. 48** and **a factor of 3 larger than the result of Ref. 15**."
> — §V.2 末段（对应 Fig. 10）

判定：**成立，可直接引用**。注意原文说的是 "predicted conductance"（图 10 在 $E_F$ 处的透射平台），
**没有**写 "T(E_F)" 字样，也没有给出不确定度。

#### (c) 峰位 −3.0/+2.5 eV（宽）与 −1.5/−1.0 eV（窄）？—— 一字不差

> "The calculated transmission shows **two broad peaks at around $-3.0$ eV and $2.5$ eV, and two more
> narrow peaks around $-1.5$ eV and $-1.0$ eV** (all energies are with respect to the Fermi level of
> the leads)."
> — §V.2 正文（对应 Fig. 10）

归属（同节后续原文）：

> "...the broad peak at $-3$ eV is exclusively due to orbital 1, while the peak at $2.5$ eV is
> exclusively due to orbital 3 ... the narrow peak at $-1.5$ eV is mostly due to orbital 3 while the
> narrow peak at $-1.0$ eV is mostly due to orbital 1"
> — §V.2（对应 Fig. 11）

> "...the transmission peaks at $2.5$ eV and $-1.5$ eV are due to orbital 3, while the peaks at
> $-3.0$ eV and $-1.0$ eV are due to orbital 1."
> — §V.2 末（对应 Fig. 12）

判定：**成立**。补充一条对本项目更有用的信息：该文的两个窄峰都在 **$E_F$ 以下**且最近的一个在
**−1.0 eV**，即"$E_F$ 近旁 1.0 eV 内无峰"这条判据的余量是 **0**（恰好边界）。

#### (d) 0.03 / 0.40 是原文还是推得？—— **推得**，原文无此二数

原文只给比例（见 (b) 引文）：比 Ref.48 小 4 倍、比 Ref.15 大 3 倍。故：

- 0.40 = 0.1 × 4（Ref.48 的值 ⇒ ≈0.4 G₀）
- 0.0333 = 0.1 ÷ 3（Ref.15 的值 ⇒ ≈0.033 G₀）；registry 取 **0.03**，比算术值更保守一点

Ref. 15 = M. Di Ventra, S. T. Pantelides, N. D. Lang, *Phys. Rev. Lett.* **84**, 979 (2000)；
Ref. 48 = K. Stokbro et al., *Comput. Mater. Sci.* **27**, 151 (2003)（均据该文参考文献表 [15]/[48]）。

判定：**[0.03, 0.40] 是 registry 自建的派生包络，不是 Thygesen 2005 的原文结论**；registry 注释里
"上下端点由该文 Ref.15(x3)/Ref.48(x4) 关系推得"**准确**。并且由于我**没有读 Ref.15/48 原文**，
这两个端点的绝对值在本轮仍是**二手**（经 Thygesen 的比例转述）。

### 2.2 Jackson et al. 2009（Surf. Sci. 603, 807–813）

**仅摘要级**（全文未获取）。摘要原文（WRAP / infona 两处记录一致）：

> "Measurements of the local adsorption geometry of the S head-group atom in the
> Au(111)(√3×√3)R30°–CH₃S surface have been made using **normal incidence X-ray standing waves
> (NIXSW)** and **S 1s scanned-energy mode photoelectron diffraction** on the same surface
> preparations. The results confirm that the local adsorption site is **atop an Au atom in a
> bulk-continuation site with a S–Au bondlength of 2.42 ± 0.02 Å**, and that there can be no
> significant fraction of coadsorbed bridging species..."

判定：**2.42±0.02 Å 成立，可直接引用**。三条必须写进 note 的适用边界：

1. 体系是 **Au(111)(√3×√3)R30°–CH₃S 表面吸附相**（甲硫醇盐在表面 atop 位），**不是**电极–分子–电极结；
2. 方法为 NIXSW + S 1s PhD 实验（非 DFT）；
3. 区间只到 **2.44 Å**（2.42+0.02），**不能**支撑 band 上端 2.49 Å。

### 2.3 `au_s_bond` 的 2.49 Å 上端：无源

两锚点并集 = [2.40, 2.44] Å。2.49 Å 在本轮读到的任何文献中都没有出处。
建议：band 改 **[2.40, 2.44]**（或保留 2.49 但在 note 显式标 `derived/unverified`）。
（注：registry 现有 note 已承认"仅方向性"、且硫醚 S–Au 可到 2.6–2.7 Å 而不适用——那个方向是对的，
只是 2.49 这个具体数字缺源。）

---

## 3. 分子族 → 预期共振表

列说明：**文献预期**=一手文献给出的通道与窗口；**本管线实测**=lead 提供的 `runs/*/TE.dat`
（`significant_peaks`，突出度 0.1）与 Seebeck 实算值，用于对照，**不作为判据源**。

| 族 | 主导通道 | 文献预期共振窗口（相对 E_F） | 允许符号 | 一手依据（DOI + locator） | 本管线实测 | 状态 |
|---|---|---|---|---|---|---|
| **BDT** Au-BDT-Au（硫醇盐） | **HOMO**（p 型） | HOMO 在 E_F 下 **≈1.2 eV**；透射窄峰 −1.0/−1.5 eV、宽峰 −3.0/+2.5 eV | S **> 0** | Reddy 2007, DOI 10.1126/science.1137149，**摘要**："the Au Fermi level position for Au-BDT-Au junctions was identified to be **1.2 eV above the highest occupied molecular orbital level of BDT**"，S=**+8.7±2.1 µV/K**，"positive sign unambiguously indicates **p-type (hole) conduction**"；Thygesen 2005 §V.2 Fig. 10（峰位，PBE 级，辅助） | 峰 −0.339 / −0.687 eV（T≈1.01/1.06）；S=**+34.5** | ✅ 有源 |
| **BPY** Au-(4,4′-bipyridine)-Au（吡啶 N） | **LUMO**（n 型） | LUMO 在 E_F **上 ≈1.5 eV**（1.47 eV 计算 / 1.53 eV 由 S、G 反推） | S **< 0** | Widawsky 2012, DOI 10.1021/nl203634m（开放全文 arXiv:1201.1837），**正文 §"resonance positions" + Table 1**：mol.4 = 4,4′-bipyridine，$S_{exp}$=**−9.5 µV/K**(HWHM 4.3)、$G_{exp}$=0.68×10⁻³ G₀；正文"pyridine-terminated molecular junctions have a **negative** Seebeck coefficient and conduct through the **LUMO**" | ±1 eV 内**无显著峰**；S=**−24.7** | ✅ 有源 |
| **AZ** Au-(2,6-azulenedithiol)-Au | **HOMO 侧**（非交替烃，谱不对称；E_F 落在 gap 中点以下） | GW 模型取 **E_F = −1.5 eV** 时给出与实验一致的连接位序（E 轴的零点定义原文未明说；按"gap 中点"理解属本人推断，见 §3.1-3） 1,3 > 2,6 ≈ 4,7 > 5,7 | S **> 0** | Xia 2014, DOI 10.1021/nl5010702（开放全文 arXiv:1411.5463），**正文 + Fig. 4(GW)**："we can infer that **E_F < 0 for the thiochroman linkers**, while E_F > 0 for the pyridine linker. Considering **E_F = −1.5 eV**, the conductance trends ... 1,3Az > 2,6Az ~ 4,7Az > 5,7Az"；"azulene ... transmission ... **non-zero at E = 0**"，1,3Az/5,7Az 在 gap 内有**两个透射零点但不在 E_F** | 峰 −0.177 / −0.406 eV（T=1.00/1.18）；S=**+54.1** | ✅ 有源（窗口为 GW 模型设定，非直接测得的 HOMO 位置） |
| **SME** Au-(1,4-bis(methylthio)benzene)-Au（硫醚） | 未定（硫醚接触；本族**无一手共振文献**） | **无源** —— 只有接触强度排序，无共振窗口 | 不作符号判断 | 仅能引：Czyszczon-Burton 2025, DOI 10.1021/acs.nanolett.5c01892，**摘要**（硫醚连接烷烃/寡聚苯，电导 Au > Ag ≥ Cu，耦合以 Au 最强）——**与共振位置无关** | 无显著峰；S=**+1.5**（近零） | ❌ 无源 → 保持 fail-safe（按未知族处理） |
| **烷烃 C4–C10**（HDT/ODT/DDT） | 无近 E_F 共振，宽 gap 隧穿 | HOMO/LUMO 均远离 E_F；长度标度 **β_N = 1.0±0.1 /碳** | 不判符号（E_F 近 gap 中心） | Xu & Tao 2003, DOI 10.1126/science.1087481，**摘要**：R = **10.5±0.5 (C6) / 51±5 (C8) / 630±50 (C10) MΩ**，β_N=1.0±0.1；Seebeck 无一手实验值，最接近的同类锚点为 Widawsky 2012 Table 1 mol.3（二苯基膦端烷烃）$S_{exp}$=**+1.1 µV/K**，"magnitude is small enough to conclude that ... E_F is very close to **mid-gap**" | 无显著峰；S=**+14 ~ +17** | ⚠ 通道有源、Seebeck 无源 |
| **2CN** Au-(2CN-656-4-thiomethylphenyl) | **推断为 LUMO 侧**（无直接文献） | **无源** | 不作符号判断 | 化合物名在 arXiv/OpenAlex/Crossref/PubMed 中**无对应文献**。只能引旁证：Venkataraman 2007, DOI 10.1021/nl062923j，**摘要**（取代基效应，苯二胺系列）："**Electron-donating** substituents ... increase the junction conductance, while **electron-withdrawing** substituents have the opposite effect ... the occupied states are closest to the gold Fermi energy ... **hole tunneling**" → 吸电子 CN 压低占据轨道、削弱 HOMO 通道；但"该化合物走 LUMO"是**推断而非文献结论** | 峰 **+0.241 eV**（T=0.871）；S=**−33.8** | ❌ 无源 → 保持 fail-safe |

### 3.1 表外需要记录的两处数字冲突（供判据 note 使用）

1. **BPY 的电导，两篇实验差 ~15×**：Xu & Tao 2003 摘要给 4,4′-bipyridine **1.3±0.1 MΩ**
   （≈9.9×10⁻³ G₀），而 Widawsky 2012 Table 1 给 **0.68×10⁻³ G₀**（≈19 MΩ）。
   两者都是 STM-BJ，环境/接触条件不同。→ BPY 的 G 只能当量级参考，**不可作闸门数值**。
2. **Seebeck 量级系统性偏差（符号稳、幅值不稳）**：BDT 实验 +8.7 vs 本管线 **+34.5**（约 4×）；
   BPY 实验 −9.5 vs 本管线 **−24.7**（约 2.6×）。两族**符号全对**。→ 这支持"`seebeck_sign` 只判符号"的现有口径，
   建议在 note 里把这条对照写成显式结论（不要用幅值判据）。
3. **AZ 的峰位与 GW 不一致**：GW 模型把 E_F 放在 gap 中点以下 ~1.5 eV，而本管线 PBE 在 −0.18/−0.41 eV
   就出现峰（即 HOMO 几乎贴上 E_F）。这与 AZ 行"PBE 低估 gap、把 HOMO 推向 E_F"一致，
   也解释了为何 S 被放大到 +54.1 µV/K。

### 3.2 `t_ef_no_peak_near_ef`（E_F 近旁 ≥1.0 eV 无峰）的源问题

- 源（Thygesen §V.2）**支持**该判据的字面形式：最近峰在 **−1.0 eV**，即恰好在 1.0 eV 边界上，余量为 0。
- 但该判据实际检验的是**泛函是否低估 gap**：本管线 PBE 对 BDT 给出 −0.339 eV 的峰，
  必然触发 warn。这与分子本身无关，是泛函/几何把 HOMO 推近 E_F 的系统偏差。
- 建议：保留 `warn`（不要升级为 error），并在 note 里写明"该判据反映泛函 gap 偏差，不是 H0 否证"；
  若要作闸门，阈值应按泛函分档（实验/GW 用 1.0 eV，纯 PBE 用更小值），否则恒触发、失去区分力。

---

## 4. 对 registry 的具体改动建议（汇总）

> **执行记录（2026-11-15）**：本表第 1–3 行的建议已落地——`t_ef_magnitude` 按「无一手区间则取消」删除；
> 派生包络并入 `g_in_theory_band` 并降为 `warn`（`evidence=derived`）；`au_s_bond` 带改为 [2.40,2.44]；
> `t_ef_no_peak_near_ef` 改名为 `t_ef_resonance_attribution`（按 `EXPECTED_RESONANCE` 族表归因，
> 未登记族按 Seebeck 符号自洽回退）；`seebeck_sign` 保留。registry 当前共 10 条判据。

| 判据 | 现状 | 建议 |
|---|---|---|
| `t_ef_magnitude` | band [0.03,0.4]，source=Thygesen p.13，gate=error | band 标注 `derived`（0.1 G₀ × Ref.15/48 的 ×3/×4），locator 改 `§V.2 + Fig.10`；**gate 降为 warn**（或按用户规则删除区间、只留 0.1 G₀ 单点） |
| `g_in_theory_band` | 同上 | 与 `t_ef_magnitude` 同处置（二者同带同源） |
| `au_s_bond` | band [2.40,2.49]，source=Thygesen p.12 | band 改 **[2.40,2.44]**；note 写明 2.40 是 Thygesen **采用 Stokbro 2003** 的输入参数；2.42±0.02 来自 Jackson 2009 **摘要**、体系为 CH₃S/Au(111) 表面相 |
| `t_ef_no_peak_near_ef` | band [1.0,None]，warn | 保留 warn；note 加"该文最近峰恰在 −1.0 eV，余量为 0；PBE 会系统性把峰推向 E_F，属泛函偏差" |
| `seebeck_sign` | gate=error，applies_to BDT | 保留；note 补"实验 vs 本管线：BDT +8.7/+34.5、BPY −9.5/−24.7，符号一致、幅值差 2.6–4×，故只判符号" |
| （新）族登记 | BDT 有源；AZ/BPY 有源但口径受限；SME/2CN 无源；烷烃通道有源、Seebeck 无源 | 见 §3"状态"列：只有 **BDT / BPY / AZ** 三族有可引的一手依据；**SME / 2CN 保持 fail-safe**；烷烃只登记"无近 E_F 共振"这一条 |

---

## 5. 引用清单（均为本次实际读到的条目）

- Thygesen, K. S., & Jacobsen, K. W. (2005). Molecular transport calculations with Wannier functions.
  *Chemical Physics*, 319(1–3), 111–125. https://doi.org/10.1016/j.chemphys.2005.05.032
  （开放全文：https://arxiv.org/abs/cond-mat/0501238 ；本文引文均出自其 §V.2 与参考文献表）
- Jackson, D. C., Chaudhuri, A., Lerotholi, T. J., Woodruff, D. P., Jones, R. G., & Dhanak, V. R. (2009).
  The local adsorption site of methylthiolate on Au(111): Bridge or atop?
  *Surface Science*, 603(5), 807–813. https://doi.org/10.1016/j.susc.2009.01.022 （仅摘要级）
- Reddy, P., Jang, S.-Y., Segalman, R. A., & Majumdar, A. (2007). Thermoelectricity in molecular junctions.
  *Science*, 315(5818), 1568–1571. https://doi.org/10.1126/science.1137149 （摘要）
- Widawsky, J. R., Darancet, P., Neaton, J. B., & Venkataraman, L. Simultaneous determination of
  conductance and thermopower of single-molecule junctions. *Nano Letters*, 12(1), 354–358.
  https://doi.org/10.1021/nl203634m （Crossref 记录年份 2011＝在线首发；开放全文 arXiv:1201.1837，
  本次读其全文，Table 1 与 p.8–9 正文）
- Xia, J., Capozzi, B., Wei, S., Strange, M., Batra, A., Moreno, J. R., Amir, R. J., Amir, E., Solomon, G. C.,
  Venkataraman, L., & Campos, L. M. (2014). Breakdown of interference rules in azulene, a non-alternant
  hydrocarbon. *Nano Letters*, 14(5), 2941–2945. https://doi.org/10.1021/nl5010702
  （题名按 Crossref 记录为 "Nonalternant"；开放全文 arXiv:1411.5463 作 "Non-Alternant"。
  本次读其全文，Fig. 4 与正文）
- Xu, B., & Tao, N. J. (2003). Measurement of single-molecule resistance by repeated formation of
  molecular junctions. *Science*, 301(5637), 1221–1223. https://doi.org/10.1126/science.1087481 （摘要）
- Venkataraman, L., Park, Y. S., Whalley, A. C., Nuckolls, C., Hybertsen, M. S., & Steigerwald, M. L. (2007).
  Electronics and chemistry: varying single-molecule junction conductance using chemical substituents.
  *Nano Letters*, 7(2), 502–506. https://doi.org/10.1021/nl062923j （摘要；全文 open access 为 arXiv:cond-mat/0703446）
- Czyszczon-Burton, T. M., et al. (2025). Thioether-linked single-molecule junctions with coinage metal
  contacts. *Nano Letters*. https://doi.org/10.1021/acs.nanolett.5c01892 （摘要，仅接触强度排序）
- El-Nahas, A. M., Staykov, A., & Yoshizawa, K. (2016). First-principles calculations of electron transport
  through azulene. *The Journal of Physical Chemistry C*, 120(17), 9043–9052.
  https://doi.org/10.1021/acs.jpcc.6b00767
  （**仅摘要级**，付费墙；摘要只给"2,6- 等四个连接位的电流较高""2,7- 整流比 ~4"，
  **未给共振峰相对 E_F 的位置**，故未用于 §3 的窗口列）

### 5.1 未能核验/未获取的项（明示）

- Jackson 2009 全文（付费墙；Warwick WRAP 记录页本次 HTTP 504 超时）→ 2.42±0.02 只有摘要级证据。
- Thygesen 2005 的 Ref.15（Di Ventra PRL 84, 979 (2000)）与 Ref.48（Stokbro CMS 27, 151 (2003)）**原文未读**
  → 0.03/0.40 的绝对值仍为二手（仅由 Thygesen 的比例转述）。
- Thygesen 2005 期刊版页码（p.12/p.13）无法用开放版核验 → 建议 locator 改用节号/图号。
- 1,4-bis(methylthio)benzene（SME）的共振通道：本次检索未找到一手输运文献 → 无源。
- 烷烃二硫醇的 Seebeck 实验值：本次检索未找到一手文献 → 无源（仅有用膦端烷烃的间接锚点）。
- "2CN-656-4-thiomethylphenyl" 该化合物本身：无文献身份 → 无源。
