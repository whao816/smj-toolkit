# AutoSMJ 自动化流程中的 DFT 方法

*AutoSMJ / smj-toolkit · 方法说明*

**泛函、基组档位与程序实际跑过的算例**

- **文档版本**：v1.1
- **编制日期**：2026-11-22
- **适用版本**：smj-toolkit 0.2.x（主链 FHI-aims + Aitranss）
- **证据基线**：docs/validation-0.2.0.md、docs/validation-0.2.2rc1.md、demo/、runs/（59 条历史 run + 2026-11-22 四臂）

*本报告只依据仓库中的源码默认值、验证记录与算例产物；未标注“已跑”的内容均为默认配置或建议，不当作已验证结论。*

> **摘要**
>
> 本报告以仓库源码默认值、验证文档、demo 与 runs/ 产物为依据，梳理 AutoSMJ 主链（FHI-aims + Aitranss）实际使用的 DFT 方法，并单独标出程序自己跑过的算例与尚未验证的边界。
>
> 核心结论：主链的电子结构统一为 FHI-aims 221103、PBE 泛函、light 档 species（H/C/N/S/Au）、Tkatchenko–Scheffler Hirshfeld vdW 与标量 ZORA；几何优化用 BFGS、力判据 1e-2 eV/Å；输运用 Aitranss NEGF、σ = 0.05/0.10 Ha、300 K。
>
> 0.2.0 版本在真实集群上完成了 BDT、AZ、2CN 三套体系的整链验证（86/92/128 原子，301 点），runs/ 中另有 40 条历史整链算例（最大 218 原子），2026-11-22 又补了 4 条 light/tight 臂（见 9.6）；所有进入 SCF 的算例均为 PBE。tight 档位原先只完成分子优化，2026-11-22 补做了 AZ 与 2CN 两条整链输运对照（见 9.6），但仍只有两个体系、一个能量窗，不构成基组收敛；主链也未做泛函、vdW 开关、k 点或电极层数收敛对照。

## 1 文档目的与范围

这份文档回答一个具体问题：AutoSMJ 自动化程序里那套 DFT 到底怎么设的——用什么泛函、什么基组、其余关键参数是什么，以及这些设置中哪些是程序自己真正跑过的，而不只是写在配置文件里的默认值。

范围限定为主链：FHI-aims 负责分子与结的电子结构、几何优化，Aitranss 负责 NEGF 输运。Gaussian 只作为对照分支在第 8 节单列，不计入主链结论。

**“测试过”的判定：** 在 demo/ 或 runs/ 中留下阶段产物（control.in、input_summary.json、results.json 等）且状态为 collected 的算例；只有配置、提交记录或计划、没有产物的，记为“已尝试/未完成”，不写成“跑过”。

**依据：** ① src/smj/config.py、src/smj/junction/tcontrol.py 与 control_scf.in/control_relax.in 的默认值；② docs/validation-0.2.0.md 与 docs/validation-0.2.2rc1.md 的验证记录；③ runs/ 下 59 条 run 的 state.json、run_spec.json 与各阶段 input_summary.json；④ demo/ 下的产物与 Gaussian 日志。

**术语：** 本文说的“基组”指 FHI-aims 的 species 档位（light/tight）。FHI-aims 用数值原子轨道，基函数、嵌入势和价电子配置都写在 species 文件里，和 Gaussian 的 6-31G(d) 不是同一套记号；文中出现 6-31G(d) 时，指的是第 8 节的 Gaussian 对照分支。

## 2 结论速览

- **生产口径：** FHI-aims 221103，PBE（GGA），light species（H/C/N/S/Au），vdw_correction_hirshfeld（TS），relativistic atomic_zora scalar；BFGS 力判据 1e-2 eV/Å；Aitranss NEGF，σ = 0.05/0.10 Ha，300 K。
- **已跑的硬证据：** 0.2.0 真实集群三体系整链（86/92/128 原子，301 点）；runs/ 40/59 条整链完成，最大 218 原子；9 条自旋极化、5 条 +1 电荷算例；所有 SCF 均为 PBE；除 2026-11-22 的两条 tight 整链（9.6）外均为 light。
- **尚未验证：** 2026-11-22 前 tight 只到 molecule_opt；现已有 AZ/2CN 两条 tight 整链（9.6），但只有两个体系、一个能量窗，仍不构成基组收敛；没有泛函对比（无 PBE0/HSE/SCAN）；没有 vdW 开关、k 点、电极层数或 SOC 对照；FHI-aims 分支不做虚频检查。
- **一句话：** PBE/light + TS-vdW + scalar ZORA 的 FHI-aims，配 0.05/0.10 Ha 的 Aitranss NEGF；这套组合在 86–218 原子的金结上跑通过；tight 档只在 AZ/2CN 上补了两条对照链。

## 3 计算流程与引擎分工

流程按固定阶段推进，每一步的状态与产物写入该 run 的 state.json；失败或未提交的阶段保持 pending/waiting，不会进入下游。主链的阶段、引擎与产物如下。

**表 1 主链阶段、引擎与产物**

| 阶段 | 引擎（主链） | 作用 | 主要产物 |
| --- | --- | --- | --- |
| molecule_opt | FHI-aims（可切 Gaussian） | 孤立分子几何优化 | geometry.in.next_step、control.in |
| aims_tip_relax | FHI-aims | Au tip 与分子接触弛豫 | tip 几何、接触距离 |
| electrode_build | 内置 fcc111 库 | 构建 37 Au/侧电极 | 结几何、表面原子索引 |
| review | 规则检查 + 人工确认 | 几何/接触/索引复核 | A2 放行记录 |
| aims_scf | FHI-aims | 结的静态自洽 | omat.aims、Mulliken、restart |
| aitranss_te | Aitranss | NEGF 透射谱 | T(E)、能级对齐 |
| junction_analyze | 后处理 | L0、S、κe、ZT | results.json |

电极来自内置库 fcc111.au.clusters.37_atoms，每侧 37 个 Au；review 阶段在验证记录中是 A2 人工确认（manual approve）。分子优化阶段可切 Gaussian（opt_engine=gaussian），但主链默认 aims；demo 中两条历史链即走 Gaussian，见第 8 节。

## 4 电子结构设置：泛函、基组与 SCF

主链的电子结构参数集中在 control_scf.in（静态 SCF）与 control_relax.in（弛豫）两份模板，由 src/smj/junction/tcontrol.py 的 make_control 填充。默认值与实际使用值汇总如下。

**表 2 主链电子结构参数（生产默认）**

| 参数 | 默认/主链取值 | 说明与来源 |
| --- | --- | --- |
| 程序与版本 | FHI-aims 221103 | 验证记录固定版本；二进制路径见 aitranss_te/input_summary.json |
| 交换关联泛函 | PBE（GGA） | xc pbe；主链唯一实际使用的泛函 |
| 基组 / species | light 档（H/C/N/S/Au） | 由 species_dir 决定；opt_basis 字段仅记录 |
| 色散校正 | vdw_correction_hirshfeld | Tkatchenko–Scheffler；后自洽能量/力修正 |
| 相对论 | atomic_zora scalar | 标量 ZORA；不含自旋轨道耦合 |
| 自旋 | none（默认）；collinear（按需） | 自旋极化算例见 9.5 |
| 电荷 | 0（默认）；+1（ttf 阳离子） | charge 字段 |
| 占据展宽 | gaussian 0.01 eV | 金属结实测 0.05 eV（见 9.5） |
| 混合器 | pulay；charge_mix 0.2；n_max_pulay 10 | 金属结实测 0.05 / 20（见 9.5） |
| SCF 阈值 | rho 1E-5；eev 1E-3；etot 1E-6 | control_scf.in |
| 力阈值 | sc_accuracy_forces 5E-3 eV/Å | 仅弛豫阶段；为几何判据的一半 |
| 迭代上限 | 500（静态 SCF）；1000（分子优化） | control_scf.in 与 molecule_opt 摘要 |
| 本征值求解 | KS_method serial | 保证 omat.aims 三角布局可被 Aitranss 读取 |
| 输出 | aitranss；mulliken；restart aims.restart | 静态 SCF |

泛函只有 PBE。仓库里没有 PBE0、HSE、SCAN 等的主链算例；Gaussian 的 B3LYP 只用于第 8 节的分子几何/频率对照，不进入输运。

vdW 是后自洽修正：它进入能量和力，因而改变弛豫几何；对固定几何的单点 SCF，它不改变哈密顿量，也就不改变 T(E)。这一点在解读“开/关 vdW”的差异时很重要。

ZORA 是标量相对论；Au 的 5d/6s 收缩被包含，但没有自旋轨道耦合，因此不适用于需要 SOC 的重元素效应。

KS_method serial 是 Aitranss 的接口约束，不是精度选择：默认的 ELPA/ELSI 路径会改用打包矩阵布局，Aitranss 读 omat.aims 时会出现负重叠本征值。control.in 里保留了这段注释。

## 5 基组档位与元素覆盖

FHI-aims 的基组不是单一的高斯收缩函数，而是每个元素一份 species 文件：里面同时给出价电子、嵌入势（如需要）以及从 minimal 到 tier 2 的数值原子轨道。所谓 light/tight 是这些轨道的两套推荐档位；主链用 light。

**表 3 主链 light 档 species 与价电子**

| 元素 | species 文件 | 价电子数 | 备注 |
| :---: | --- | :---: | --- |
| H | 01_H_default | 1 | light 档 |
| C | 06_C_default | 4 | light 档 |
| N | 07_N_default | 5 | light 档 |
| S | 16_S_default | 6 | light 档 2020 版加入 f 函数（文件头注释） |
| Au | 79_Au_default | 31 | 6s1 5p6 5d10 4f14；4f 在价层 |

Au 的 31 价电子值得单独说明：4f 被放进价层，意味着 f 态参与自洽和输运；这也是 100–200 原子金结 SCF 成本的主要来源之一。light 档是在这个规模下能承受的取舍，而不是基组收敛性研究的结果。

tight 档位：compare-2026-11-15 建了 species-tight，2cn、az-thiol、az-thiolate 三条链只跑到 molecule_opt。2026-11-22 另建四臂（AZ-dithiol、2CN 各 light/tight）跑到 junction_analyze，见 9.6；即便如此，tight 仍只有两个体系、一个能量窗，不能用来声称基组收敛。

Au 的 species 文件带 4f 是 light 档的默认；如果换成不含 4f 的 Au 势，价电子数和结果都会变，仓库没有做过这个对照。

## 6 几何优化与界面构建

**表 4 几何优化、tip 弛豫与电极参数**

| 参数 | 取值 | 说明 |
| --- | --- | --- |
| 分子优化引擎 | FHI-aims（opt_engine=aims） | 默认；Gaussian 为对照分支 |
| 分子优化泛函/基组 | PBE / light | 与静态 SCF 一致 |
| 分子优化力判据 | relax_geometry bfgs 1.e-2 eV/Å | 判据为最大力分量 |
| 弛豫 SCF | sc_accuracy_forces 5E-3 eV/Å；sc_iter_limit 1000 | control_relax.in |
| 能量判据 | aggregated_energy_tolerance 1.0e-1 | 防止 BFGS 在能量平台反复 |
| tip 初始 Au–S | 2.35 Å | ttf 算例用 2.6 Å |
| 接触约束 | pin_contact false；strip_thiol_h null | 主链默认；算例按需显式设置 |
| tip 优化 | 全原子 BFGS；电极成结后冻结 | relax=false |
| 电极 | fcc111.au.clusters.37_atoms；每侧 37 Au；nlayers 3 | 库文件建议 3 层 |
| 结规模 | 86–218 原子（实测） | 见 9.4 |
| 表面原子索引 | 每个结 3+3 个 | review 阶段人工复核 |

力判据用的是最大力分量，不是 RMS；因此 1e-2 eV/Å 比很多 Gaussian 默认的 RMS 判据更严。

tip 弛豫可能失败：AZ 的 tip 作业在第 80 步报 trust radius 过小（原文见 9.2），程序按 geometry.in.next_step + init_hess（diag 97）重启，20 步收敛。这是流程里已知的稳健性处理，不是方法参数。

strip_thiol_h 主链默认 null；实际算例中硫醇盐（bdt、az、ddt、hdt、odt 等）显式设 true，硫醚（2cn、sme 等）显式设 false。电极成结后冻结，因此 tip 弛豫只动分子和 tip 的 Au，不重塑电极。

## 7 输运计算参数

**表 5 Aitranss 输运参数**

| 参数 | 配置默认 | 验证/生产实际 |
| --- | --- | --- |
| 程序与方法 | Aitranss 071813；NEGF-Landauer | 零偏压 T(E) |
| 自能展宽 | σ_near 0.05 Ha；σ_far 0.10 Ha | 近/远电极 |
| 能量窗（默认） | (-2.0, 0.005, 2.0) eV 相对 E_F；801 点 | config.py |
| 能量窗（0.2.0 验证） | (-0.22, 0.0002, -0.16) Ha；301 点 | BDT/AZ/2CN 三体系 |
| 电子温度 | 300 K | 只进入费米分布/输运积分 |
| 输出 | T(E)、L0、log10(G/G0)、S、κe、ZT | results.json |

两套能量窗都实际跑过：一套是相对 E_F 的宽窗（默认，801 点），一套是围绕 Au 费米能级的绝对 Ha 窄窗（验证，301 点）。窄窗覆盖 E_F 上下约 30 meV，并保留两侧共振；这是精度与成本的折中，不是两套物理。

300 K 是输运积分里的电子温度，不是 FHI-aims 的电子展宽；后者是 0.01 eV（金属结 0.05 eV）的数值展宽，两者不要混为一谈。

L0 是 300 K 费米窗内积分得到的零阶输运矩（电导，单位 G0 = 2e²/h），T(E_F) 是费米能处的透射值，二者接近但不相等；S、κe、ZT 分别由 L1、L2 组合得到。

网格敏感性：BDT 的 estep 从 2e-4 Ha 加密到 5e-5 Ha（301→1201 点），L0/S/κe/ZT 相对差 &lt; 1e-8；T(E_F) 的差别只来自最近网格点。固定 S–Au = 2.42 Å 的几何重复了同一对照。

## 8 Gaussian 对照分支（非主链）

仓库保留了 Gaussian 分支，用于分子几何优化与频率检查；它不属于本次统计的主链口径，但属于“程序涉及过的方法”，因此单列。

**表 6 Gaussian 分支的实际运行记录**

| 用途 | 方法 | 结果 |
| --- | --- | --- |
| 2CN demo 分子优化 | B3LYP/6-31G(d)，opt+freq，G16 | Normal termination；156 模，0 虚频；E = -2061.10992343 Ha |
| 2CN fullrelax demo | 同上 | Normal termination；156 模，0 虚频；E = -2061.10992343 Ha |
| Azulene demo 分子优化 | B3LYP/def2-TZVP，opt+freq，G16 | Normal termination；54 模，0 虚频；E = -1182.42130005 Ha |
| g16tip 的 PBE 精修尝试 | PBE/def2-SVP（Gaussian） | 未正常结束（无 Normal termination）；最后 SCF 能约 -2329.5568 Ha；弃用 |
| 配置默认（仅 opt_engine=gaussian） | method b3lyp；basis 6-31g(d) | 主链不触发 |

这些 Gaussian 结果只用于分子几何和频率对照；所有 tip 弛豫、静态 SCF 与输运仍是 FHI-aims PBE/light。因此它们不能用来回答“主链泛函是什么”。

B3LYP 两条链都做了 freq，0 虚频，说明分子几何是极小点；FHI-aims 分支不做虚频检查（README 明确）。如果论文需要“无虚频”证据，目前只能引用这两条 Gaussian 对照。

g16tip 的失败值得记录：用 Gaussian PBE/def2-SVP 精修 Aitranss 前的 tip 团簇时，作业没有正常结束，最后一条 SCF 能量约 -2329.5568 Ha，未产生可用几何；该 demo 停留在轮询状态。这说明该分支不适合大团簇，主链继续用 FHI-aims。

## 9 程序实际跑过的算例

### 9.1 真实集群整链验证（0.2.0）

验证节点为 gpu-outside（8 CPU 核 + 1 GPU），FHI-aims 221103、Aitranss 071813；PBE/light，37 Au/侧、nlayers 3，301 点、300 K，pin_contact false，A2 人工确认。三套体系均整链完成，结果见下。

**表 7 0.2.0 验证体系与界面**

| 体系 | 分子 | 锚定化学 | 分子原子数 | 结原子数 | Au–S / Å | 表面索引（左/右） |
| :---: | --- | --- | :---: | :---: | :---: | --- |
| bdt | 1,4-BDT（C6H6S2） | 硫醇盐（去 H） | 14 | 86 | ≈2.25 | 40,41,44 / 76,77,80 |
| az | azulene-2,6-dithiol（C10H8S2） | 硫醇盐（去 H） | 20 | 92 | ≈2.25 | 46,47,50 / 82,83,86 |
| 2cn | 2CN-656-4-thiomethylphenyl | 硫醚（不去 H） | 54 | 128 | 2.42（配位） | 82,83,86 / 118,119,122 |

**表 8 0.2.0 三体系输运结果（301 点，300 K）**

| 体系 | E_F / Ha | T(E_F) | L0 / G0 | log10(G/G0) | S / (µV/K) | κe / (pW/K) | ZT | 最大 T（能量） |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| bdt | -0.188142 | 0.272587 | 0.278854 | -0.5546 | +34.49 | 165.58 | 0.04656 | 1.0545（-0.687 eV） |
| az | -0.188087 | 0.349526 | 0.380578 | -0.4196 | +54.07 | 223.99 | 0.11548 | 1.1795（-0.406 eV） |
| 2cn | -0.179465 | 0.001253 | 0.001325 | -2.8778 | -33.84 | 1.24 | 0.02837 | 0.8710（+0.241 eV） |

SCF 用时分别为 bdt 1564 s、az 1963 s、2cn 2896 s；TE 阶段 219–251 s（8 核节点）。这些数字说明 128 原子、31 价电子 Au 的 light 档计算在单节点上是可承受的，也是选择 light 的现实依据。

### 9.2 事件与数值对照

（1）AZ tip 弛豫的 trust radius 中止。作业 26119 在第 80 步报：

> The trust radius has fallen below its allowed minimum value, 0.100000E-03 Angstrom

处理方式是按 geometry.in.next_step + init_hess（diag 97）重启为作业 26120，20 步收敛。该事件只影响 AZ 的 tip 几何生成，不改变方法参数。

（2）BDT 能量网格加密。estep 0.0002 → 0.00005 Ha（301 → 1201 点），L0/S/κe/ZT 相对差 &lt; 1e-8；T(E_F) 的微小差别来自最近网格点。固定 S–Au = 2.42 Å 的几何上重复了同一对照，结论一致。这说明 301 点对这套体系已经足够。

（3）2CN 的低 T(E_F)。2CN 的 T(E_F) 仅 1.25e-3，原因归于硫醚配位接触（Au–S 2.42 Å）使前线轨道与 Au 耦合弱；谱在 +0.24 eV 处有 0.87 的峰。这不是计算失败，而是接触化学决定的物理结果，也解释了为何 2CN 的 S 与 bdt/az 反号。

0.2.0 的工程判定为通过；物理判定当时为 inconclusive，后续按 t_ef_resonance_attribution 复核为 passed。方法层面的结论不受这一判定变化影响。

### 9.3 纯 FHI-aims demo 链（2CN-aimsopt）

demo/2CN-656-4-thiomethylphenyl-aimsopt 是唯一分子优化也用 FHI-aims 的 demo；另外两条 2CN demo 的 molecule_opt 是 Gaussian（见第 8 节）。

**表 9 aimsopt demo 的方法与结果**

| 项目 | 数值/设置 |
| --- | --- |
| 体系与规模 | 2CN-656-4-thiomethylphenyl；54 原子分子 → 128 原子结；301 点 |
| 分子优化 | FHI-aims PBE/light；BFGS 1e-2；收敛最大力 0.9289e-2 eV/Å；HOMO-LUMO gap 1.0008 eV |
| tip 弛豫 | 2.35 Å 初值；pin_contact false；收敛最大力 0.9976e-2 eV/Å |
| 静态 SCF | 128 原子；xc pbe；spin none；charge 0；occupation 0.01；sc_iter_limit 500；vdw hirshfeld |
| TE 输入 | nsaos 4742；valence_electrons 2456；E_F -0.179360798048 Ha；301 点 |
| 输运结果 | T(E_F)=0.00124463；L0=0.001320209；log10=-2.8794；S=-35.1335 µV/K；κe=1.2664 pW/K；ZT=0.02991；max T=0.9992 @ +0.240 eV |

与 0.2.0 集群验证的 2CN 相比，同一分子、不同分子几何（aims vs Gaussian）给出 T(E_F)=1.24e-3 vs 1.25e-3、S=-35.13 vs -33.83 µV/K；数量级一致，但说明分子几何来源会影响热电量。

demo 目录下另有 fullrelax 与 base 两条链，其 molecule_opt 用 Gaussian B3LYP/6-31G(d)，aims 各阶段仍为 PBE/light；fullrelax 的 L0=0.001284826、S=-29.09 µV/K，base（钉扎接触）的 L0=0.000477053、S=+2.80 µV/K。它们不是纯 FHI-aims 分子优化，只作为对照，不进入主链统计。

### 9.4 runs/ 语料覆盖

runs/ 下共 59 条 run（不含 _archive），40 条所有阶段 collected、19 条部分完成。覆盖情况如下。

**表 10 runs/ 语料覆盖矩阵**

| 维度 | 覆盖情况 |
| --- | --- |
| run 总数 | 59（不含 _archive）；完整 40，部分 19 |
| 引擎 | FHI-aims 59/59；Gaussian 仅 demo 的分子优化 |
| 泛函 | PBE 59/59；进入 SCF 的 51 条均为 pbe |
| species | 留下 molecule_opt 摘要的 47 条中：light 44，tight 3（仅 molecule_opt）；另有 2026-11-22 的 2 条 tight 整链（见 9.6） |
| 自旋 | 51 条 SCF 摘要：none 42，collinear 9 |
| 电荷 | 0（多数）；+1（ttf 5 条） |
| 结规模 | 86–218 原子 |
| 锚定 | S–Au 多数；N–Au（bpy、ppd、di-fcfr/di-pcpr）；按 Au 索引（ttf、ifs、et7h） |
| 能量窗 | 绝对 Ha 301 点（验证体系）与 eV 相对 E_F（默认窗口）两套 |
| 最大体系 | 218 原子（et7h-bsm2）；214（di-fcfr）；210（di-pcpr） |

部分完成的 19 条分四类：44dpa 的 8 条参数扫描（d225–d285、tw15–tw90）；bdt 的 3 条几何/网格对照（bdt-dense、bdt-anchor24、bdt-anchor24-dense）；a2-z-pinned、et7h-bs2 各 1 条；compare-2026-11-15 的 6 条 light/tight 对照。

完整链的 40 条见附录 A；其中 ttf-lone-pair 与 ttf-lone-pair-release 状态完整但 junction_analyze 没有数值结果，附录中以“—”标出。

### 9.5 非默认参数的实际使用

下列设置不是默认值，但确实出现在 smj.yaml 中并被程序执行过；它们说明程序对金属结 SCF、接触几何、自旋/电荷和网格加密都至少跑过一次，但不构成系统收敛研究。

**表 11 非默认设置的实际算例**

| 类别 | 实测设置 | 算例 | 状态/结果 |
| --- | --- | --- | --- |
| 金属结 SCF | occupation 0.05；mix 0.05；pulay 20 | ifs-bs、ifs-triplet、di-fcfr、di-pcpr | 整链完成；152–214 原子 |
| 接触钉扎 | pin_contact true；2.413–2.55 Å | a2-e、ppd-pin、sme-axis（完整）；bdt-anchor24（部分） | 用于接触几何对照 |
| 自旋极化 | collinear；M = 1 或 2 | ttf-* 5 条、di-* 2 条、ifs-* 2 条 | 完成（ttf-lone-pair* 无结果数值） |
| 电荷 +1 | charge 1；fixed moment 1 | ttf-cation/fixed/release 等 | 完成 |
| 网格加密 | estep 5e-5 Ha（1201 点） | bdt-dense、bdt-anchor24-dense | 与 2e-4 一致（&lt;1e-8）；bdt-dense 仅部分阶段 |
| tip 初距 | 2.6 Å | ttf-cation 等 | 完成 |
| tight species | species-tight | compare-2026-11-15 的 3 个体系（仅 molecule_opt）；2026-11-22 的 AZ/2CN（整链，见 9.6） | 前者无输运；后者整链完成 |

### 9.6 light vs tight 整链输运对照（2026-11-22）

为回答“换成 tight 档会不会改变输运结果”，2026-11-22 补做了四臂整链：AZ-dithiol 与 2CN 各做 light/tight 两臂；同分子、同生产协议（不钉扎、tip 初距 2.35 Å、电极成结后冻结）与同能量窗（[-0.22, 0.0002, -0.16] Ha，301 点，300 K），唯一变量是 species 档位。四臂均跑到 junction_analyze，各阶段 collected。

AZ 用中性二硫醇（C10H8S2，20 原子；tip 阶段 strip 掉 S–H），与 0.2.0 的 AZ 构造一致。最初按 compare-2026-11-15 的 −2 游离阴离子起的两臂，在 tip 阶段的 SCF 跑满 500 步未收敛（Δρ 在 0.1–0.8 间振荡，反复出现 “VBM is higher than the internal energy zero”），已弃用：−2 分子配中性 Au tip 不是生产接触模型，也无法进入建结。

**表 12 light/tight 整链输运结果（301 点，300 K）**

| 体系 | 档位 | 基函数 | E_F (Ha) | L0 (2e²/h) | S (µV/K) | κe (pW/K) | ZT | T(E_F) | max T |
| --- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| AZ-dithiol（20 原子） | light | 4364 | -0.188093 | 0.373290 | +53.35 | 220.640 | 0.11194 | 0.342941 | 1.2406 |
| AZ-dithiol（20 原子） | tight | 5368 | -0.190031 | 0.486402 | +65.64 | 242.818 | 0.20060 | 0.444570 | 0.9832 |
| 2CN（54 原子） | light | 4742 | -0.179492 | 0.001313 | -33.53 | 1.239 | 0.02769 | 0.001244 | 0.9997 |
| 2CN（54 原子） | tight | 6436 | -0.181567 | 0.001361 | -3.53 | 0.853 | 0.00046 | 0.001338 | 0.8913 |

差值（tight − light）：AZ 的 L0 从 0.373290 变到 0.486402（+30.3%），S 从 +53.35 到 +65.64 µV/K，κe 从 220.640 到 242.818 pW/K，ZT 从 0.11194 到 0.20060；2CN 的 L0 从 0.001313 变到 0.001361（+3.7%），S 从 -33.53 到 -3.53 µV/K，κe 从 1.239 到 0.853 pW/K，ZT 从 0.02769 到 0.00046。

解读：AZ 的差异温和（L0 +30%、ZT +79%），2CN 的 S/ZT 则大幅塌缩（ZT −98%）。原因在于 S ∝ dT/dE 在 E_F 处对谱形极敏感：2CN 的 T(E_F) 只变 +7.5%，但斜率几乎被抹平（S 从 −33.5 到 −3.5 µV/K）；AZ 则因 E_F 附近的窄共振，tight 下 L0 与 ZT 反而更高。注意两档的分子与 tip 几何都在各自档位下重新优化，差异包含几何响应，不是同一几何的单点基组对照。两条链各只做了一次、一个能量窗，且没有做基函数逐档收敛，因此不能把 tight 当成“收敛值”；但 tight 档已有整链输运数据，旧结论“tight 无输运”不再成立，定量报 ZT 时也不应无条件只用 light。

## 10 未测试内容与适用边界

以下内容在仓库中没有完整算例，写方法部分时不能当作已验证：

- **基组：** light 是主链；tight 已有 AZ/2CN 两条整链（见 9.6），但只有两个体系、一个能量窗，不构成基组收敛。
- **泛函：** 主链只有 PBE；没有 PBE0、HSE、SCAN 的输运对比；Gaussian B3LYP 只到分子几何。
- **vdW：** 主链始终 hirshfeld，没有开/关 vdW 的对照。
- **周期性/k 点：** 团簇 + NEGF，Gamma 点；没有 k 网格、周期性电极或层数收敛。
- **相对论：** scalar ZORA；没有自旋轨道耦合。
- **偏压：** 只有零偏压 T(E)；没有有限偏压或非平衡计算。
- **温度：** 只跑 300 K；没有温度依赖。
- **环境：** 气相；没有溶剂或环境模型。
- **电极：** 37 Au/侧、3 层；没有 55 原子电极或 4/5 层的收敛研究。
- **频率：** FHI-aims 分支不做虚频；只有 Gaussian demo 有 freq。
- **接触扫描：** 44dpa 距离/扭转扫描与 bdt-anchor24 都只完成部分阶段，未形成完整曲线。
- **体系规模：** 最大 218 原子；更大体系未测。

如果要把这套方法写进论文或审稿回复，建议至少补第三/第四个体系的 tight 整链、一组 PBE0 单点、一组 vdW 开关，以及 4/5 层电极对照。以上是建议，仓库中没有这些算例，本报告不把它们当作已跑结果。

## 11 复现与核查路径

**表 13 参数与产物位置**

| 内容 | 位置 |
| --- | --- |
| 默认参数 | src/smj/config.py（分子/电极 129–200 行；输运 225–306 行） |
| control 模板 | src/smj/junction/templates/control_scf.in、control_relax.in |
| 生成逻辑 | src/smj/junction/tcontrol.py（make_control 的默认值） |
| 流程与参数表 | docs/smj-计算流程梳理.md |
| 集群验证 | docs/validation-0.2.0.md（BDT/AZ/2CN 整链） |
| demo 基线 | docs/validation-0.2.2rc1.md（4 份谱基线） |
| 单算例产物 | runs/&lt;组&gt;/&lt;run&gt;/{molecule_opt,aims_tip_relax,aims_scf,aitranss_te,junction_analyze} |
| tight species | runs/compare-2026-11-15/species-tight |
| light/tight 整链对照 | docs/validation-light-vs-tight-2026-11-22.md；runs/light-tight-2026-11-22/ |

核查顺序：state.json（阶段是否 collected）→ 各阶段 input_summary.json（本轮实际参数）→ control.in（交给 FHI-aims 的原文）→ junction_analyze/results.json（数值结果）。demo/ 下 state.json 的阶段名 gaussian_opt 是历史命名，对应 molecule_opt，不代表本轮一定跑了 Gaussian；以 molecule_opt/input_summary.json 的 engine 字段为准。

## 附录 A 40 条整链完成 run 清单

数据来自 runs/&lt;组&gt;/&lt;run&gt;/ 下的 state.json、run_spec.json、aims_scf/input_summary.json 与 junction_analyze/results.json。run id 即末级目录名；锚定为左/右锚定原子（S/N 或 Au 索引）；pin 为接触是否钉扎；col = collinear；L0 单位 G0；S 单位 µV/K；“—”表示状态完整但无结果数值。本表是 2026-11-22 之前的 40 条历史整链，不含当日新增的 light/tight 四臂（见 9.6）。

**表 A1 整链完成 run（按 run id 排序）**

| # | run | 原子数 | 锚定 | pin | 自旋 | q | L0 / G0 | S / (µV/K) |
| :---: | --- | :---: | :---: | :---: | :---: | :---: | ---: | ---: |
| 1 | 2cn | 128 | S/S | 否 | none | 0 | 0.0013 | -33.8 |
| 2 | 33dpa | 106 | S/S | 否 | none | 0 | 5.916e-06 | -20.6 |
| 3 | 33dpa-syn | 106 | S/S | 否 | none | 0 | 9.112e-06 | -24.6 |
| 4 | 34dpa | 106 | S/S | 否 | none | 0 | 3.787e-04 | -10.0 |
| 5 | 44dpa | 106 | S/S | 否 | none | 0 | 0.0138 | +6.2 |
| 6 | a2-e | 106 | S/S | 是 | none | 0 | 0.0091 | +9.2 |
| 7 | a2-z | 106 | S/S | 否 | none | 0 | 0.8938 | +17.8 |
| 8 | az | 92 | S/S | 否 | none | 0 | 0.3806 | +54.1 |
| 9 | az-sac | 104 | S/S | 否 | none | 0 | 0.0360 | -12.2 |
| 10 | az-sh | 92 | S/S | 否 | none | 0 | 0.4176 | +74.3 |
| 11 | az-sme | 100 | S/S | 否 | none | 0 | 0.0807 | -3.5 |
| 12 | azo-e | 98 | S/S | 否 | none | 0 | 0.1303 | +60.2 |
| 13 | azo-z | 98 | S/S | 否 | none | 0 | 0.0621 | +58.2 |
| 14 | bdt | 86 | S/S | 否 | none | 0 | 0.2789 | +34.5 |
| 15 | bpy | 94 | N/N | 否 | none | 0 | 0.0012 | -24.7 |
| 16 | ddt | 106 | S/S | 否 | none | 0 | 6.095e-05 | +17.0 |
| 17 | di-fcfr | 214 | N/N | 否 | col | 0 | 0.0139 | -119.3 |
| 18 | di-pcpr | 210 | N/N | 否 | col | 0 | 0.0064 | -119.8 |
| 19 | et7h-bsm2 | 218 | 135/136 | 否 | none | 0 | 2.875e-05 | +229.1 |
| 20 | hdt | 94 | S/S | 否 | none | 0 | 0.0016 | +13.6 |
| 21 | ifs-bs | 152 | 39/41 | 否 | col | 0 | 0.0371 | +6.2 |
| 22 | ifs-triplet | 152 | 39/41 | 否 | col | 0 | 0.0249 | +27.6 |
| 23 | ind-sac | 103 | S/S | 否 | none | 0 | 0.0142 | +4.1 |
| 24 | ind-sh | 91 | S/S | 否 | none | 0 | 0.3170 | +39.8 |
| 25 | ind-sme | 99 | S/S | 否 | none | 0 | 0.0512 | +4.3 |
| 26 | nap-sac | 104 | S/S | 否 | none | 0 | 0.0081 | +0.6 |
| 27 | nap-sh | 92 | S/S | 否 | none | 0 | 0.1903 | +42.1 |
| 28 | nap-sme | 100 | S/S | 否 | none | 0 | 0.0297 | +2.8 |
| 29 | odt | 100 | S/S | 否 | none | 0 | 2.936e-04 | +15.1 |
| 30 | ppd | 90 | N/N | 否 | none | 0 | 0.0452 | +10.2 |
| 31 | ppd-pin | 90 | N/N | 是 | none | 0 | 0.0254 | +33.6 |
| 32 | sme | 94 | S/S | 否 | none | 0 | 0.0564 | +1.5 |
| 33 | sme-axis | 94 | S/S | 是 | none | 0 | 0.0039 | +2.2 |
| 34 | sme13 | 94 | S/S | 否 | none | 0 | 0.0021 | -12.9 |
| 35 | sme14 | 94 | S/S | 否 | none | 0 | 0.0561 | +1.7 |
| 36 | ttf-cation | 88 | 6/10 | 否 | col | 1 | 0.2337 | +10.7 |
| 37 | ttf-fixed | 88 | 6/10 | 否 | col | 1 | 0.2337 | +10.7 |
| 38 | ttf-lone-pair | 88 | 6/10 | 否 | col | 1 | — | — |
| 39 | ttf-lone-pair-release | 88 | 6/10 | 否 | col | 1 | — | — |
| 40 | ttf-release | 88 | 6/10 | 否 | col | 1 | 0.2337 | +10.7 |

## 附录 B 术语表

**表 B1 文中缩写与术语**

| 术语 | 含义 |
| --- | --- |
| FHI-aims | 全电子数值原子轨道 DFT 程序；主链电子结构引擎 |
| species / light / tight | FHI-aims 的逐元素基组档位；light 为主链，tight 已有 AZ/2CN 两条整链（见 9.6） |
| PBE | Perdew–Burke–Ernzerhof GGA 交换关联泛函 |
| TS-vdW / hirshfeld | Tkatchenko–Scheffler 色散校正；后自洽能量/力修正 |
| ZORA | 零阶正则近似相对论（scalar，标量） |
| NEGF | 非平衡格林函数；Aitranss 用它算 T(E) |
| T(E_F) | 费米能处的透射 |
| L0 | 300 K 费米窗内积分得到的零阶输运矩（电导，单位 G0） |
| G0 | 电导量子 2e²/h |
| S / κe / ZT | 塞贝克系数 / 电子热导 / 热电优值 |
