# smj — 单分子结输运计算链（独立包）

> 从**分子结构**到**透射谱**的端到端自动化：FHI-aims 分子优化 → FHI-aims 加针尖/弛豫 →
> 构建电极（本地三角岛或 FHI-aims 电极库；默认 Au(111)）→ FHI-aims 静态 SCF → aitranss 透射 → Landauer 分析。
> 整链可一把跑完，也可**逐阶段独立运行**；每个阶段的输入/产物都落在
> `<run.dir>/<run.name>/<stage>/`，状态机在 `state.json` 里可查、可恢复。
>
> 本包是 AgenticDFT 主仓 SMJ 链的独立打包（发行名 `smj-toolkit`，
> import 名 `smj`，命令行 `smj`）。自 0.2.2 起 junction 内核
> （`src/smj/junction/`）由**本仓维护**，与主仓的关系从「代码同源」改为
> 「结果同契约」：来源与本地变更史记在 `PROVENANCE.json`，数值一致性由
> `demo/demo-reference.json` 与 `scripts/crosscheck_agenticdft.py` 守。

## 0. 定位与边界

| | |
|---|---|
| 输入 | 分子 `.xyz`（锚定原子默认硫醇 -SH 的 S，也支持 `-SCH3` 等硫醚，或 N/Se 等任意元素符号 / 1-based 序号指定锚）+ 电极/锚点/输运窗口参数（`smj.yaml`） |
| 输出 | `TE.dat`（透射谱）、`results.json`（`log10_G_over_G0_like`、`S_uV_per_K`、`k_e_pW_per_K`、`ZT_electronic`）；`smj report` 另生成 `manifest.json`、`report.md` 与判定（`manifest.verdict`） |
| 适用 | 单分子结 / 分子电子学：fcc(111) 电极（默认且已完整验证 Au(111)）间硫醇/硫醚及杂原子（N/Se 等）锚定分子的电输运与热电 |
| 不承诺 | 自带 Gaussian/FHI-aims/aitranss；自动替人选方法/基组/锚点；替代人工结构确认；把 20 Au/侧与 Au60 口径混比；**自动加对照臂/敏感性计算（pinned、变锚定、变团簇、空结等）——这些只在用户明确要求时跑** |

六阶段链（顺序固定，但可单独执行任意阶段）：

```text
molecule.xyz
  └─ molecule_opt      FHI-aims PBE/light 分子优化      → opt.xyz
      └─ aims_tip_relax  加电极金属 tip（默认 Au）+ 全原子弛豫 → geometry.in.next_step
          └─ electrode_build  构建电极：把每个 tip 扩展为完整 fcc(111) 团簇（默认 Au(111)）→ geometry.in + electrode_metadata.json
              └─ aims_scf     FHI-aims 静态 SCF      → aims.out / mos.aims
                  └─ aitranss_te  aitranss 透射      → TE.dat
                      └─ junction_analyze  Landauer  → results.json
```

`au_relax` 已从主链路移除：结几何在 `electrode_build` 后即冻结，后续 SCF/输运不再弛豫。
该阶段模块仍保留，仅供旧 run 或显式 `smj stage au_relax ...` 兼容使用。

首段规范名是 `molecule_opt`（2026-11-13 起；旧名 `gaussian_opt` 仍可作为命令别名）。
旧 run 的目录/状态键若仍是 `gaussian_opt/`，`RunLayout` 与 `RunState` 会回退读取，
`smj report` 仍按 `molecule_opt` 展示；新 run 一律写 `molecule_opt/`。

建结阶段规范名是 `electrode_build`（2026-11-13 起；旧名 `au60_build` 仍可作为
命令别名；更早的 `au20_build`/`auN_build` 也会归一化到该阶段）。旧 run 的
`au60_build/` 或 `auN_build/` 目录由 `RunLayout` 回退读取，`smj report` 按
`electrode_build` 展示；新 run 一律写 `electrode_build/`。阶段名不再绑定
“每侧 60 个 Au”的具体电极尺寸。

`molecule_opt` 的引擎由 `molecule.opt_engine` 选择：默认 `aims`（FHI-aims
PBE/light，与 tip/SCF 同一软件栈）；`gaussian` 保留旧口径（B3LYP/6-31G(d)，
需 Gaussian 16 许可）。两种引擎的产物契约相同：`opt.xyz` 供 `aims_tip_relax`
消费。collect 是 fail-closed 的：aims 分支要求日志出现
`Present geometry is converged`（FHI-aims 弛豫达到步数上限时也会正常退出并写
`geometry.in.next_step`，所以"有 Have a nice day"不等于收敛）；两种引擎都要过
`smj.checks` 的几何合理性检查（最近原子对 ≥0.7 Å；重原子键消失、孤立原子拒绝；
含 H 键消失（质子迁移）只警告），优化后的最大力、HOMO/LUMO/gap 写进
`results.json`。弛豫的 SCF 精度也显式写进 `control.in`（`sc_accuracy_rho 1E-5`、
`sc_accuracy_eev 1E-3`、`sc_accuracy_etot 1E-6`、`sc_accuracy_forces 5E-3` eV/Å、
`sc_iter_limit 1000`）——旧模板靠 FHI-aims 自动默认，力的精度是 `not_checked`；
开壳层写 `spin collinear` + `fixed_spin_moment = M−1`（=2S，不是多重度）固定总自旋，
可用 `molecule.initial_moments` 按 1-based 原子序号给初始磁矩（FHI-aims 在 collinear
下要求初始自旋密度；不给自己会在 notes 里警告）；结级 `aims_tip_relax`/`aims_scf`
会把 `molecule.initial_moments` 映射到结几何的原子前缀（tips/Au 追加在后，去质子只删
锚定 S–H 的氢，序号自动左移），存在逐原子初猜时不再写统一的 `default_initial_moment`；
`molecule.fixed_spin_moment: false` 则**不写**约束行（只留 `spin collinear` + 初猜），
用于"释放约束后的解"对照。非零 `molecule.charge` 时把总电荷
均分成 per-atom `initial_charge` 写进 `geometry.in`（FHI-aims 要求 `charge` 等于
该行之和；集群上的 221103 不认 `override_initial_charge_check`，实测会 parse_control
停算，所以采用满足一致性检查的写法）。日志里
出现 SCF 未收敛直接拒绝。
FHI-aims 分支仍不做虚频检查（无现成一键 freq）。

柔性分子先用 `smj conformers --smiles "Sc1ccc(S)cc1"` 做 ETKDGv3 + MMFF94s
构象粗筛（按能量排序 + 重原子骨架 RMSD 去重），把 top-k 各跑一次下游链再比较；
MMFF94s 只负责粗筛，不是电子结构判据。

本阶段优化的是中性分子，而器件里锚定原子与 Au 成键。`smj variants --xyz molecule.xyz`
可派生两条化学态对照臂：`thiolate`（去 S–H 的 H、电荷 −1/个，模拟去质子化）与
`aucap`（S–H 的 H 换成沿原方向的 Au，模拟已与电极成键），连同
`variants.json`（含 `config_hint`、`source_sha256` 与口径 caveat）。两条臂各跑一次
`molecule_opt`，比较 C–S 键长、gap 与 `aims_tip_relax` 后的锚–Au 接触几何——
同一 2CN-656 硫醚在两条 demo 里 S–Au 分别是 2.69 Å/141° 与 2.43 Å/111°，
这就是初猜化学态的影响量级。

基组口径现在可核对：`molecule.opt_basis` 只是标签，`smj species --dir <species_dir>`
会读 species 文件头声明的档位与基组规模、记录逐文件 sha256；每次 `molecule_opt`
prepare 也把同样的 `species.json` 写进阶段目录，collect 搬进 `results.json.species`，
标签与实际档位不一致时给警告。柔性分子用 `smj ensemble --smiles ... --dir ensemble/`
把 top-k 构象排成批处理（每个构象一个隔离 `--run-dir`），跑完
`smj ensemble --dir ensemble/ --report` 给 gap / 最大力的极差——这就是构象贡献的量级。

## 0.1 已验证体系（0.2.0）

三套真实体系已在 Slurm 集群上用本包跑通全链（FHI-aims + aitranss），数值与耗时见
[`docs/validation-0.2.0.md`](docs/validation-0.2.0.md)：

| 体系 | 锚定化学 | 结原子数 | 每侧 Au | G/G0 (L0) | 物理判定 |
|---|---|---|---|---|---|
| 1,4-benzenedithiol (BDT) | 硫醇盐 | 86 | 37 | 0.279 | passed（HOMO 峰 −0.34 eV 已归因） |
| azulene-2,6-dithiol (AZ) | 硫醇盐 | 92 | 37 | 0.381 | passed（HOMO 峰 −0.18 eV 已归因） |
| 2CN-656-4-thiomethylphenyl | 硫醚 | 128 | 37 | 0.00132 | passed（LUMO 峰 +0.24 eV 自洽归因） |

三套的**工程判定均通过**（SCF 收敛、TE 数据完整、阶段哈希闭环、连接核对有效）。
2026-11-15 起，E_F±1.0 eV 内的透射峰不再直接判 `inconclusive`，而是按族预期通道归因
（`t_ef_resonance_attribution`）：峰在预期侧且 0.15≤|ΔE|≤1.0 eV 记「物理信号已归因」；
2CN 未登记文献通道，按实测 Seebeck 负号自洽归入 LUMO 侧。BDT 族判据只对标定族生效，
AZ/2CN 不做该族比对。

## 1. 安装

要求 **Python ≥ 3.10**；运行时依赖只有 **PyYAML** 与 **numpy**。发布包已在 Ubuntu 26.04 / Python 3.14.4 上完成 wheel 与 sdist 的用户视角安装冒烟（`smj verify` / `init` / `plan`）；CI 矩阵覆盖 Python 3.10–3.14（Ubuntu）与 3.12–3.14（Windows）。

```bash
# 源码可编辑安装（在仓库根执行；pyproject.toml 在根，源码在 src/）
python -m pip install -e .

# 或构建 wheel 后安装（产物在 dist/）
python scripts/gen_junction_provenance.py           # 先刷新内核台账（PROVENANCE.json 随 wheel 分发）
python scripts/build_wheel.py --check       # 只用标准库，不需要 setuptools/build；--check 按 RECORD 校验
python -m pip install dist/smj_toolkit-0.2.2rc2-py3-none-any.whl
```

可选依赖：

```bash
python -m pip install -e ".[plot]"   # matplotlib，用于 T(E) 出图与连接核对图
python -m pip install -e ".[dev]"    # pytest，跑本包测试
python -m pip install -e ".[plot,dev]"    # plot + dev
```

安装后验证（在仓库根执行）：

```bash
smj --help
smj verify      # Python/PyYAML/numpy、junction 内核、模板、阶段模块（按当前 STAGE_ORDER；legacy au_relax 只作历史 run 兼容）
```

仓库级一致性检查（CI 同款，改动内核/demo 前后各跑一次）：

```bash
python scripts/gen_junction_provenance.py --check     # PROVENANCE.json 的 sha256 台账（repo_commit 不参与比对）
python scripts/gen_demo_reference.py --check  # demo 谱的数值基线（G/G0、S、k_e、ZT）
```

## 2. 5 分钟快速开始

```bash
# 0) 建一个工作目录，把 molecule.xyz 放进去（init 只写配置，不复制/生成结构）
mkdir smj-demo && cd smj-demo

# 1) 生成配置骨架（默认 dry-run，不提交任何作业）
smj init . --name demo --xyz molecule.xyz --mode dry-run

# 2) 编辑 ./smj.yaml：分子路径、电极参数；真跑时补 cluster 段
#    （见 §4；只想先看输入装配的话保持 run.mode: dry-run 即可）
#    注意：默认 electrode.source: library，需要 FHI-aims 的 electrodes.library
#    （用 electrode.library_dir 指向它，或写 source: builtin 回退到自建电极）

# 3) 看六阶段计划、依赖与当前状态
smj plan

# 4) 只生成第一阶段输入（人工检查 mol.gjf / sub_g16.sh）
smj stage molecule_opt --prepare

# 5) 干跑整链：逐阶段 prepare，不提交（预期停在缺上游产物的下一阶段，退出码 1）
smj run --dry-run
```

配置在 `./smj.yaml`，run 目录是 `./demo/`（`run.dir=.`、`run.name=demo`，
即 `<run.dir>/<run.name>/`）。`dry-run` 的终点是 `prepared`：
它**不会**伪造 `mol.log`/`opt.xyz` 让整链看起来跑通，所以整链干跑会停在需要上游产物的
下一阶段（该阶段标 `failed`，阻塞原因写进 note，`smj run` 退出码 1）——这是刻意的
fail-closed 语义，不是缺陷。

真跑把 `run.mode` 改成 `slurm-ssh`，补齐 `cluster.host/remote_root/partition` 与
`cluster.env`，然后（仍在 `smj-demo/` 下）：

```bash
smj preflight          # 提交前检查
smj run                # 整链
# 或从某阶段续跑 / 跑到某阶段：
smj run --from aims_scf --to aitranss_te
```

## 3. 逐阶段独立运行

阶段之间只通过文件耦合：每个阶段读上游 stage 目录里的约定文件，覆盖参数优先。
因此任意阶段都可以单独 `--prepare` / `--submit` / `--wait` / `--collect`。

下表示例省略 `-c <path>`（默认读当前目录的 `smj.yaml`）；在别处执行请加
`--config my-run/smj.yaml`。`--run` 在 `dry-run` 模式下只 prepare（见 §2）。

| 阶段 | 默认输入 | 主要产物 | 需提交？ | 示例命令 |
|---|---|---|---|---|
| `molecule_opt` | `molecule.xyz`（或 `--xyz`） | 默认（aims）：`geometry.in`、`control.in`、`sub_aims.sh`、`input_summary.json`，collect 后 `aims.out`、`opt.xyz`；gaussian 引擎：`mol.gjf`、`sub_g16.sh`，collect 后 `mol.log`、`opt.xyz`、`gaussian_summary.json` | 是（FHI-aims 或 Gaussian） | `smj stage molecule_opt --run --xyz molecule.xyz` |
| `aims_tip_relax` | `molecule_opt/opt.xyz`（或 `--geometry`） | `geometry.in`、`control.in`、`sub_aims.sh`；collect 后 `geometry.in.next_step` | 是（FHI-aims） | `smj stage aims_tip_relax --run --geometry molecule_opt/opt.xyz` |
| `electrode_build` | `aims_tip_relax/geometry.in.next_step`（或 `--geometry`） | `geometry.in`、`electrode_metadata.json`（含电极元素/构建方式/背层索引与库团簇来源）；tip 被扩展为完整电极团簇，S–金属距离/角度沿用弛豫结果 | **否**（本地） | `smj stage electrode_build --run --geometry aims_tip_relax/geometry.in.next_step` |
| `aims_scf` | `electrode_build/geometry.in`（或 `--geometry`） | `control.in`、`sub_aims.sh`；collect 后 `aims.out`、`mos.aims` | 是（FHI-aims） | `smj stage aims_scf --run` |
| `aitranss_te` | `aims_scf/`（或 `--scf-dir`） | `tcontrol`、`sub_te.sh`；collect 后 `TE.dat`、`te.out` | 是（aitranss） | `smj stage aitranss_te --run --scf-dir /path/to/scf` |
| `junction_analyze` | `aitranss_te/TE.dat`（或 `--te`） | `results.json` + 统一版式 `T_E.png`/`T_E.pdf`（无 matplotlib 时跳过图） | **否**（本地） | `smj stage junction_analyze --run --te TE.dat` |

FHI-aims 阶段（`aims_tip_relax`/`aims_scf`）还需要 `cluster.species_dir` 或
`paths.species_defaults` 指向 species 目录；缺失时 `prepare` 直接 `StageBlocked`
（`smj preflight` 也会报 error）。`--geometry` 在三个消费几何的阶段
（`aims_tip_relax`/`electrode_build`/`aims_scf`）都可用，但语义不同：`aims_tip_relax`
要**已优化分子的 XYZ**（`opt.xyz`），不是 FHI-aims 的 `geometry.in`；`electrode_build`
与 `aims_scf` 要的才是 FHI-aims 几何。`aitranss_te` 除 SCF 目录外还要
`electrode_build/electrode_metadata.json`（六个表面索引与左右 tip 索引），可用
`--metadata` 覆盖。

`electrode_build` 是**「构建电极」**阶段：**怎么构建**与**构建什么电极**都是本阶段参数。

- **构建方式** `electrode.source`：`library`（默认）用 FHI-aims 自带的 `electrodes.library`
  预定义团簇贴装电极（手册 §5.3.1 的要求），由 `electrode.library_dir`、
  `electrode.library_cluster` 控制（默认
  `fcc111.au.clusters.37_atoms/fcc111.cluster.natoms37.x1.aims`）；库文件不进 wheel，只在
  metadata 里记 `sha256` 与相对路径。`electrode.nlayers: null` 表示取库文件头部建议的层数。
  `builtin` 用纯 Python 在本地构建 fcc(111) 三角岛（旧自建路径，`Electrode111`，非 Au 走
  `ElectrodeFCC`）。
- **电极类型** `electrode.element`（默认 `Au`）+ `electrode.a_lattice`（默认 `null`）：
  fcc(111) 单晶面金属。已登记晶格常数 Au / Ag / Cu / Ni / Pt / Pd / Al；其它元素必须显式
  给出 `a_lattice`（2.0–6.0 Å）。`aims_tip_relax` 生成的 2 个 tip 也用同一元素，因此 tip–S
  距离/角度记录的是该金属的真实接触。
- **一致性**：`source: library` 时库团簇的元素必须与 `electrode.element` 一致，否则阶段
  fail-closed（提示改 `library_cluster` 或 `element`）；`builtin` 时每侧原子数仍由
  `electrode.au_per_electrode` 控制（`library` 模式忽略该键）。
- **产物**：metadata 里非 Au 电极会写 `electrode_element` / `electrode_a_lattice`
  （默认 Au 保持旧 metadata 逐字节不变；下游与旧 run 缺该键按 `Au` 处理），SCF/透射/自能/
  校验/连接核对出图都从 metadata 读取电极元素。

默认协议是 **tip 全原子弛豫 → 扩展成完整团簇 → 几何冻结**：`electrode.pin_contact: false`
（默认）让 `aims_tip_relax` 的全部原子（分子 + 2 个电极金属 tip，Au 体系即 56 原子）自由
弛豫，不用 `anchor_au_s_target` 限制金属–S，这一步给出的就是可靠的金属–S 键长与连接角度；
`electrode_build` 只沿弛豫后的 tip→linker 方向把每个 tip 扩展为完整团簇，不改 S–金属
距离/角度。扩展完成后结几何完全冻结，`aims_scf`/`aitranss_te` 不再做任何结构弛豫。

`electrode.relax` / `relax_freeze` 字段仅为旧 run 或显式重跑 `au_relax` 保留；
主链路已不消费它们。
默认 `relax: false` 时这个开关不生效。

`transport.energy_window: [起点, 步长, 终点]` 默认是 **eV 且相对费米能级**：
默认 `[-2.0, 0.005, 2.0]` 就是 `E − E_F ∈ [−2, +2] eV`、步长 5 meV（801 点）。
`aitranss_te` 在 prepare 时读 `aims_scf/aims.out` 的 efermi，把它平移成 `tcontrol`
需要的 Ha 绝对 `$ener/$estep/$eend`；缺 efermi 直接 `StageBlocked`，不猜。
换金属/换分子后 E_F 会移动，但相对窗口不跟着漂——这就是 0.2.x 起默认 ±2 eV 的理由。
旧配置的 Ha 绝对三元组（如 `[-0.21, 0.0002, -0.17]`）仍可读：跨度 ≤0.1 Ha 且量级 <2 Ha
时按旧口径解释并给迁移提示，也可显式写 `energy_window_unit: Hartree`；
`transport.ener/estep/eend` 仍是 **Ha 绝对**的单分量覆盖（优先级最高）。
`TE.dat` 第一列是绝对 Hartree，第二列才是
`E−E_F`（eV）。`transport.require_te_out`（默认 `true`）要求 `aitranss_te/te.out`
也在场——它是 aitranss 主输出，缺了就无法复核这次透射计算；外部/历史 run 才显式设
`false` 豁免。

`molecule.family` 声明判据标定族（`AZ` = Au-(2,6-azulenedithiol)-Au、`BDT` 等）。
BDT 的绿带（T(E_F)/G ∈ [0.03, 0.40]、Seebeck 为正）只对 `Au-BDT-Au` 闭壳层作闸门：
**族未登记**或**开壳层**（`multiplicity > 1`）时两条判据只记录、不作闸门，物理判定回落
`inconclusive`（`no_applicable_gate`），不会把别的族的包络套到当前体系上——参见
`runs/ttf-release` 的复核记录（G=0.2337 G0、S=+10.72 µV/K 落在 BDT 带内，但那不是
该自由基结的证据）。

动作开关（每个阶段通用）：

| 开关 | 含义 | 什么时候用 |
|---|---|---|
| `--prepare` | 只生成输入/作业脚本，不提交 | 人工检查 `control.in`、资源申请 |
| `--submit` | 清理残留 → 上传 → `sbatch`，记录 job id | 想让排队与准备分离 |
| `--wait` | 轮询到结束；超时不判失败 | 长作业先退出，稍后回来再等 |
| `--collect` | 回取产物并解析（可单独执行） | 手工跑完/从别处拷回产物后只解析 |
| `--run` | prepare → submit → wait → collect（默认） | 一把梭；本地阶段 = prepare → collect |

覆盖参数只改变**本阶段读哪个文件**，不会改动上游 `state.json`；独立跑出的产物不会自动
把上游阶段标成 `collected`。常用组合：

```bash
# 只做输入装配，拷到集群手工提交
smj stage aims_scf --prepare --geometry /data/run/electrode_build/geometry.in --config smj.yaml

# 外部已有 TE.dat，只做分析（--te 也可给 aitranss_te 输出目录）
smj stage junction_analyze --run --te /data/te/TE.dat --config smj.yaml
```

## 4. Slurm 配置

`cluster` 段只在 `run.mode: slurm-ssh` 时必需；`dry-run` / `local` 下可留空。
`slurm-ssh` 至少要求 `host`、`remote_root`、`partition`，缺一直接报配置错误。

```yaml
run:
  name: bdt-trans
  dir: .                 # run 根目录；最终产物在 <dir>/<name>/
  mode: slurm-ssh        # dry-run | local | slurm-ssh
  clean: false           # true 时 prepare 前清空本阶段目录（仅限 run 目录内）
  keep_going: false      # true 时整链遇到失败继续跑后续阶段

cluster:
  host: cluster.example.edu
  user: yourname
  port: 22
  key: ~/.ssh/id_ed25519
  remote_root: /data/yourname/smj
  partition: main
  ntasks: 4
  gres: ""
  account: ""
  qos: ""
  species_dir: /data/software/fhi-aims/species_defaults/defaults_2020/light
  time_limit: ""
  mem: ""
  poll_seconds: 60
  timeout_hours: 6
  env:
    g16root: /data/software/gaussian/g16
    aims_bin: /data/software/fhi-aims/bin/aims.221103.scalapack.mpi.x
    aitranss_bin: /data/software/aitranss/aitranss
    mkl_lib: /data/software/intel/mkl/lib/intel64
```

`env` 的键与内核 `junction/inputs.py` 的 `ENV_SPECS` 同源（**区分大小写，请写小写**），
常用项：

| 键 | 用途 | 默认 |
|---|---|---|
| `g16root` / `g16_module` | Gaussian 安装根 / module 名 | 空（`g16root` 或 `g16_module` 二选一） |
| `g16_scratch` / `g16_mem` / `g16_walltime` / `g16_gres` | Gaussian 临时目录/内存/墙钟/GRES | `16G` / `01:00:00` / 空 |
| `aims_prefix` / `aims_module` / `aims_bin` / `aims_ld_module` | FHI-aims 定位方式 | 空（三者择一） |
| `aims_mem` / `aims_walltime` / `aims_gpus` | FHI-aims 资源 | `16G` / `02:00:00` / 空 |
| `aitranss_bin` | aitranss 可执行文件 | 空 |
| `mkl_lib` | aitranss 依赖的 MKL 库路径 | 空 |
| `te_mem` / `te_walltime` / `te_gres` | aitranss 资源 | `8G` / `02:00:00` / 空 |

缺省值不静默成功：输入渲染器在缺 `g16root`/`g16_scratch` 等键时，会把
`UNCONFIGURED` 写进作业脚本（例如 `export g16root=UNCONFIGURED`），这样的脚本在集群上会
直接失败。`smj preflight` 在**非 `dry-run`** 模式下检查 Gaussian 许可、FHI-aims 与
aitranss 可执行文件，缺一项即 error（退出码 1）；`molecule.opt_engine: aims` 时分子优化
走 FHI-aims，Gaussian 路径/许可检查被跳过（记 info `gaussian_not_required`）。
`smj run` 不会自动替你跑 preflight，
真跑前请先跑一次。`dry-run` 不触网、也不核对二进制，所以干跑通过不等于集群能跑。
0.1.1 起提交前会再扫一遍作业脚本：仍含 `UNCONFIGURED` 时 `--submit` / `smj run` 直接
拒绝（`StageError`，退出 1），不再让它到集群上才失败。

### 4.1 GRES 与节点独占（全项目通用）

**纯 CPU 作业也必须显式写 `--gres=gpu:0`。** 这条不是本仓的偏好，是 Slurm 的行为：

- 节点声明了 GRES（`/etc/slurm/slurm.conf` 里 `NodeName=... Gres=gpu:4`），而分区没有
  为 GRES 设默认值；Slurm 的默认分配方式是**把节点全部 GRES 给该作业**。
- 后果：作业本身只用了 8 核，却把 4 张卡一起占住，**同一节点上第二个作业一律
  PENDING(Resources)** —— 即使 `squeue` 显示 `CPUAlloc=8/CPUTot=64`、内存也空着。
  判断方法：`sbatch --test-only` 给出的预计开始时间正好等于正在跑那个作业的时限终点。
- `--gres=none` **无效**（Slurm 不接受该写法）；`--gres=gpu:0` 才是"不要卡"。
- 实测（2026-11-24，gpu1，64 核 / 4×RTX 4090）：不写 GRES 时 12 h 的作业把后续
  24 个数组任务全部堵到时限之后；加上 `--gres=gpu:0` 后 8 个数组任务 + 2 个弛豫作业
  同时跑满 64 核。

规则实现在 `src/smj/slurm.py`（叶子模块，全项目共用，Gaussian / FHI-aims / aitranss /
热导链的作业脚本都走它）：

| 配置 | 渲染结果 |
|---|---|
| 不写（`gres` / `g16_gres` / `te_gres` / `aims_gres_sched` 留空） | `#SBATCH --gres=gpu:0` |
| 写 `gpu:1` | `#SBATCH --gres=gpu:1` |
| 写 `off` / `none` / `-` / `no` | **整行不渲染**（给没有 GRES 的集群用） |

注意 `aims_gpus` 是另一回事：它表示**真正的 GPU 计算模式**（1 GPU ↔ 1 rank），
走 `#SBATCH --gpus=N`，与调度用的 `aims_gres_sched` 互不影响。热导链的 FHI-aims
是纯 CPU（`mpirun -np N`），只用后者。

## 5. 三种运行模式与离线演练

| 模式 | 行为 | 适用 |
|---|---|---|
| `dry-run` | 只 `prepare`，产出输入/脚本；不建传输层、不提交 | 本机装配输入、CI、评审 |
| `local` | 在本机 `bash <script>` 同步执行，哨兵 `.dsh_done` 记退出码 | 有 Gaussian/FHI-aims 的工作站、调试 |
| `slurm-ssh` | 上传 → `sbatch` → 轮询 → 回取 | 真实集群 |

`local` 需要 `bash`（Windows 上装 Git-Bash/WSL 并加入 PATH）；没有 bash 时
`LocalTransport` 会明确报错并建议改用 `dry-run`。

不想连集群也想演练远端链路时，用 `smj.testing.FakeTransport`（内存实现
`run/put/get/mkdirs`，记录每条命令；`simulate_done(remote_dir, rc)` 可造哨兵）：

```python
from smj.config import load_config
from smj.runner import StageRunner
from smj.testing import FakeTransport
from smj.transport import CommandResult

transport = FakeTransport()

def on_run(command, cwd, timeout):
    if "sbatch" in command:                    # 假调度器：给出可解析的 job id
        return CommandResult(0, "Submitted batch job 4242")
    return None                                 # 其余命令按默认成功处理

transport.on_run = on_run
config = load_config("smj.yaml")               # 在 run 配置所在目录执行
# 演练远端链路要求 run.mode 为 slurm-ssh（local 走本机 bash，不演练上传/sbatch）；
# 不想改文件时可在内存里改（下面的 host 等只用于通过字段检查）：
config.run.mode = "slurm-ssh"
config.cluster.host = "fake-cluster"
config.cluster.remote_root = "/tmp/smj-fake"
config.cluster.partition = "main"
runner = StageRunner(config, transport=transport, dry_run=False)
runner.prepare("molecule_opt")
sub = runner.submit("molecule_opt")            # 只走内存：清理 → 上传 → 记录 job_id=4242
print(sub["job_id"], sub["cleaned"], len(sub["uploaded"]))
print(transport.commands[-2:])                 # 看实际发出的命令
```

`--wait` 的轮询探针本身是一条 shell 命令，FakeTransport 需要按
`@@SENT/@@QRC/@@ID` 协议模拟输出才能驱动（可参考 `tests/test_state_machine.py`
的 `_SbatchTransport`）；只验证提交/清理逻辑时到 `submit` 为止即可。FakeTransport
不模拟排队与真实程序行为。

## 6. run 目录布局

```text
<run.dir>/<run.name>/            # smj.yaml 的 run.dir 是 run 根，默认 "." → ./<name>/
├── state.json                   # 阶段状态机（status/job_id/时间戳/资源/哈希/results）
├── run_spec.json                # 运行参数快照（完整解析配置、协议指纹、软件版本）
├── run_spec.previous.json       # 仅 --force 重绑时保留的上一份快照
├── manifest.json                # 产物清单（sha256、大小、阶段、回取模式）
├── report.md                    # 人读报告：概览/阶段状态/输入摘要/资源/输运/失败/产物
├── saved_structures/            # 结构归档（已 gitignore）：xyz/geometry.in 原样副本 + 结果快照 + manifest.json
│   ├── 0001_molecule_opt_input.xyz
│   ├── 0002_molecule_opt_optimized.xyz
│   ├── 0003_aims_tip_relax_molecule_tips.in
│   ├── 0004_electrode_build_junction.in
│   ├── 0005_au_relax_junction_relaxed.in
│   └── 0005_au_relax_junction_relaxed.results.json
├── molecule_opt/                # 每阶段一个目录，输入+产物都在里面
│   ├── mol.gjf  sub_g16.sh  input_summary.json
│   └── mol.log  opt.xyz  gaussian_summary.json
├── aims_tip_relax/
│   ├── geometry.in  control.in  sub_aims.sh
│   └── geometry.in.next_step  aims.out
├── electrode_build/
│   └── geometry.in  electrode_metadata.json
├── aims_scf/
│   ├── control.in  sub_aims.sh  geometry.in
│   └── aims.out  mos.aims  basis-indices.out
├── aitranss_te/
│   ├── tcontrol  sub_te.sh  self.energy.in
│   └── TE.dat  te.out
└── junction_analyze/
    └── results.json  T_E.png  T_E.pdf
```

`state.json` 是唯一状态来源：`smj status` 读它；`smj report` 由它 + 产物目录生成
`manifest.json`/`report.md`（默认不覆盖已存在的报告，`--force` 才重写）。

`saved_structures/manifest.json` 把每份结构与产出它的计算绑在一起：`stage`/`kind`/
`format`（`xyz` / `aims`）/`engine`/`sha256`/`natoms`/`source`/`created_at`；带结果的条目还记
结果快照文件名、快照 `sha256` 与标量摘要。归档点覆盖结构链的每一版几何：

- `molecule_opt`：`input`（prepare 时的分子结构）、`optimized`（collect 的 `opt.xyz`）；
- `aims_tip_relax`：`molecule_tips`（分子 + 两枚 tip，collect 的最终几何）；
- `electrode_build`：`junction`（补全两侧电极后的完整分子结，collect 校验通过后）；
- `au_relax`：`junction_relaxed`（贴装后弛豫的新几何；`electrode.relax=false` 的直通
  分支不重复归档，直接复用 `junction` 条目）。

结构文件与源文件逐字节一致（xyz 仍是 `.xyz`，FHI-aims 的近邻表仍是原生 `.in`），
可直接用 `sha256` 对账；归档目录已在 `.gitignore` 中，重复 prepare/collect 会追加新条目
而不是覆盖历史。归档失败只记阶段 note，不阻断计算。

## 6.1 声子热导（thermal 链）

给定 Au–分子–Au 结，算声子热导 `G_ph`（pW/K）与声子透射谱 `T_ph(ħω)`。**它不是 κ**：
单分子结没有确定的截面积，块体的傅里叶定律 / Wiedemann–Franz / 动力学理论都不适用，
能报的只有热导 `G`（W/K）本身。

**前置条件比电子链少**：只需要 FHI-aims（任何带 `compute_forces` 的版本即可，有限位移
不依赖 DFPT），电极力常数随包发布，不用自己先算一份体相：

```bash
$ smj lead-fc --show          # 资产身份：Au / 4.078230 Å / 3×3 / 2 层 / sha256
$ smj lead-fc --export my_lead   # 导出后可用 thermal.lead_fc 指向自己的版本
```

`smj.yaml` 里打开（其余键都有可用默认值）：

```yaml
thermal:
  enabled: true
  electrode: {in_plane: [3, 3], layers: 4, device_layers: 2, vacuum: 20.0}
  fc:        {displacement: 0.01, cutoff: 10.0, array_parallel: 8, k_grid: [1, 1, 1]}
  negf:      {eta_meV: 0.05, energy_step_meV: 0.05, band_top_meV: 19.5}
```

热导链**独立于电子六阶段链**（D12）：它不出现在 `smj status`/`smj plan` 的默认输出里，
也不接受 `run --from/--to` 区间，逐阶段用 `smj thermal --stage X --run`。

| 阶段 | 做什么 | 产物 | 上机 |
|---|---|---|---|
| `thermal_device` | 建 Au(111) 3×3 周期板（`layers` 层 × 两侧）+ 去质子化分子，标 `layer_index`/`side`/`role`；`S–Au = gap/2 + au_s_distance − S–S/2` | `geometry.in`、`device_metadata.json` | 本地 |
| `thermal_relax` | 弛豫（冻结最外层 `freeze_outer_layers` 层），模板带 `compute_forces .true.`、`preconditioner kerker 1.5`、`k_grid K K 1` | `geometry.in.next_step`、`aims.out` | 1 个作业 |
| `thermal_fc` | 对器件区每个自由度做 ±0.01 Å 中心差分（6N 个点）→ 读力 → 力常数 | `force_constants.npy`、`displacements.json`、`forces_all.json` | 数组作业 |
| `thermal_analyze` | 拼 lead 力常数 → 动力学矩阵 → ASR 诊断 → 表面自能 → `T_ph` → Landauer 积分 | `thermal.json`、`thermal_spectrum.*`、`thermal_conductance_vs_temperature.*` | 本地 |

nap-sh 实测（arXiv 级标定，2026-11）：90 原子（Au 72 + 分子 18）、器件区 54、**324 个位移**、
8 路并发 **约 9 h**（单点 20–35 min）；`G_ph(300 K) = 13.72 pW/K`、`G_Q = 283.93 pW/K`、
有效通道数 `N_eff = 0.048`，300 K 已达高温平台的 98.9%。

判定走 `smj review` 的 `thermal` 块，十条判据：`ph_transmission_present`、`ph_relax_converged`、
`ph_band_top_zero`、`ph_quantum_bound`、`ph_acoustic_limit`、`ph_asr_residual`、`ph_fc_symmetry`、
`ph_imaginary_modes`、`ph_lead_cell_match`、`ph_grid_independence`。error 门未过 → 顶层
`blocked`；warn 门未过或**缺证据** → `thermal.verdict = inconclusive`（不会判 pass），
顶层 `status` 仍是 `passed/warning/blocked` 三值，不改变既有退出码语义。

两条容易踩的硬约束：

1. **热导器件与电子链器件不是同一体系**。热导用 Au(111) 周期板（有主层才能做半无限表面
   自能），电子链用金字塔团簇；`G_ph` 与 `G/G0` **不可并列**，`smj report` 的
   `## 5. 热导（声子）` 里有强制提示。
2. **η 是数值正则化参数，不是物理线宽**。换 η 等于换方法：本体系 η 由 0.5 降到 0.05 meV，
   `G_ph` 由 10.9 升到 13.7（+25%）。协议值冻结在 0.05 meV / 0.05 meV 步长（方案 D15），
   与文献横比必须同 η，或显式声明为「η→0 外推值」（外推约 13.8 pW/K）。

已知限制：简谐近似（非谐修正文献估计 ~15%）；线性响应给不出热整流；`ZT_total` 要求同一
几何的电子透射，本链不提供（`thermal.json` 里 `zt.available = false` 并写明原因）。

### 6.1.1 扫参与跨 run 比较

`G_ph` 对正则化参数（η、能量步长）极其敏感，所以报告数字前先把「这个数有多依赖 η」
量化出来，再谈分子之间的差别：

```bash
$ smj thermal --sensitivity -c smj.yaml     # 不碰集群：复用已有 force_constants.npy
  η\步长           0.100      0.050      0.025
  0.300         12.099     12.372          —
  0.200         12.799     13.122     12.961
  0.100         13.355     13.614     13.403
  0.050         13.494     13.717     13.449   ← 生产配置
  η→0 外推：13.819 pW/K
  步长引起的相对变化最大 2.46%（梯形积分底噪，不是收敛趋势）
```

写 `thermal_analyze/sensitivity.json`。η→0 只用**生产步长下最小的两个 η** 做渐近
线性外推——`G(η)` 是饱和的，全局拟合会被低端陡坡拖出比最小 η 还小的值。

多个分子放一张表时用 `smj thermal-compare`，它**先比协议再比数**：

```bash
$ smj thermal-compare runs/            # 输出 markdown 表；--out 写文件，--json 给机器
$ smj thermal-compare runs/ --strict   # 协议不一致直接返回非零
```

协议项为温度、η、能量步长、`band_top_meV`、`lead_fc_sha256`；任一不一致都会在表前
点名（跨协议的 `G_ph` 互比等于换方法，见 D15）。

数组作业掉点时不必重跑整批：

```bash
$ smj thermal --stage thermal_fc --retry missing            # 或 --retry 17,203
```

`--retry` 复用已 prepare 的脚本，只把 `--array=` 换成缺的那几个槽位；`fc_summary.json`
里若有 SCF **未收敛**的槽位会直接报错（未收敛的 SCF 照样写出形状正确的力，靠条数
和形状是查不出来的）。

### 6.1.2 环境清单与收敛轴

`thermal_analyze` 除 `thermal.json` 外还写两份东西，用来回答「换个团队能不能复现」
和「这个数扫过什么」：

- **`environment.json`**（`smj.environment.v1`）—— 决定数值的输入与工具：smj/Python/
  numpy 版本、FHI-aims 二进制 sha256、species 目录内容摘要、lead 力常数 sha256、
  热导配置快照。`identity()` 只取六项（版本、numpy、aims、species、lead FC），
  **不含路径**，所以两台机器装在不同目录仍是同一套工具；`digest()` 剔除
  `generated_at` 便于比对。`thermal-compare` 会单独报「环境不一致」，与「协议不一致」
  分开——前者是换了工具，后者是换了方法。
- **`thermal.json` 的 `convergence`** —— 三条收敛轴。力常数截断半径零机时可得
  （切已组装矩阵的 lead 块后重算 NEGF）；lead 层数与器件区大小必须重跑位移集，
  所以只登记状态与所需作业。

实测（nap-sh，η=0.05 meV）：

```text
截断 4/5/6/7/8/10/12 Å → 10.55 / 10.72 / 12.42 / 12.20 / 11.77 / 13.12 / 13.72 pW/K
相对散布 23.0%；该散布在 η=0.05–0.3 meV 都是 22–23%，所以是真实截断效应而非采样噪声
名义 fc_cutoff_A = 10 Å，但组装矩阵里最远的非零元素在 11.91 Å → 按名义值切会丢 4.3%
```

最后一行是这套轴最有用的产物：它说明 `10 Å` 的 lead 力常数**没有收敛**——2 层块的
平移对称映射把体相里 6 Å 内的一对摆到了器件坐标的 11.9 Å 处，「名义截断」和「实际
生效」不是一回事。报告会把它写成显式警告并给修复动作（把 `fc_cutoff_A` 提到 ≥12 Å
重算电极力常数，54 个 SCF，约半天）。

## 7. 附带功能

| 命令 | 作用 | 产物/输出 |
|---|---|---|
| `smj preflight` | 提交前检查：分子可读/原子重叠、species 覆盖、建结几何闸门、资源上限；非 `dry-run` 时再加集群/SSH 连通、远端目录可写、二进制/许可 | 结构化 issues（error/warning/info），有 error 退出 1 |
| `smj status` | 读 `state.json`，打印阶段/状态/job id/远端目录/耗时/备注 | 终端表格（`--json` 机器可读） |
| `smj report` | 由 run 目录生成 `manifest.json` + `report.md`（9 节：概览/判定/结构复核/阶段状态/输入摘要/资源/输运/失败原因/产物清单）；缺 `TE.dat` 时物理判定为「无法判定」而不是通过 | 报告文件；默认跳过已存在文件，`--force` 重写，`--no-checks` 不算判定，`--strict` 在判定失败/阻断时退出 1 |
| `smj analyze` | 对 `aitranss_te/TE.dat` 做 Landauer 分析（`log10_G_over_G0_like`、`S_uV_per_K`、`k_e_pW_per_K`、`ZT_electronic`），写 `junction_analyze/analysis.json`；`--plot` 可顺手按统一版式出 `T_E_plot.png`（`[plot]` extra） | `analysis.json`（含 `independent` 独立复核结果） |
| `smj plot` | 按**统一版式**画 T(E)（四边框/刻度朝内/E_F 虚线/log 纵轴/固定配色，见 `docs/plot_style.md`）；可给 TE.dat、run 目录或逗号标签叠加多条 | `T_E.png` + `T_E.pdf`（默认 PNG+PDF 同源导出，600 dpi） |
| `smj pdos` | 分组 PDOS：读 `aims_scf` 的 `mos.aims`（开壳层 `alpha.aims`/`beta.aims`，按 `$uhfmo_*` 头行判通道）与 `basis-indices.out`，按 `--groups molecule=1-20,left=21-52` 把原子分组投影到 α/β，占据数取 `aims.out` 的 `State Occupation` 表（闭壳层只取最后一张；开壳层取末两张分别给 α/β，行数允许不同，但与通道波函数态数不符时显式报错而不是按索引截断），落在 `TE.dat` 的 E − E_F 网格上（与 T(E) 同一能量轴）；权重按每态 `Σc²` 归一（NAO 非正交，真正 Mulliken 需要 `omat.aims` 重叠矩阵，本地只有表头）；`--plot` 出上下两栏共享 x 轴的图。默认回取只抓 `mos.aims` 表头，本地没有完整本征向量时会**显式失败**并提示用 `--fetch` 按清单整取 | 终端摘要或 `--json`；`-o` 写 `smj.pdos.v1` JSON；缺完整本征向量时退出 1 |
| `smj review` | 结构复核门：建结几何 + 连接核对（原 A2：指纹/表面索引偏差）。默认由管线在 `electrode_build` 之后自动登记（`who=auto`，审计不过则拒绝登记并阻断）；也可人工签 `--approve --who <名字>`（审计或出图不通过时**拒绝签署**，退出 2），`--no-figure` 只在你已在别处看过图时使用，`--connection-check/--no-connection-check` 临时覆盖 `review.connection_check` | `review.json`、`electrode_build/connection_check.json`（旧 `a2_confirmation.json` / `surface_check.json` 仍可读）；`blocked` 时退出 1 |
| `smj spin` | 开壳层自旋汇总（只读）：自旋声明、双通道 Tα/Tβ(E_F)、P(E_F)/P(L0)、Gα/Gβ/Gtotal、Sα/Sβ/S_total、P_S；旧 run 的 spin 块缺分通道 S 时直接从 `TE.dat` 重算；Hirshfeld 分区与分子内 +/− 自旋矩、review 已有的自旋 code 与口径提醒 | 终端摘要或 `--json`（`run_id`/`molecule`/`open_shell`/`spin`/`partition`/`findings`/`review_codes`/`caveats`）；不写任何文件 |
| `smj audit` | 本地一致性检查（不连远端）：state 与阶段目录/产物/job id 对账 | issues 列表；有 error 退出 1 |
| `smj verify` | 安装自检：Python/PyYAML/numpy、junction 内核、模板、阶段模块 | issues 列表；有 error 退出 1（不需要 `-c`） |
| `smj plan` | 打印六阶段计划、依赖与当前状态，不落盘、不提交 | 终端表格 |

退出码：`0` 成功/作业已结束；`1` 失败或配置/输入错误；`2` 未完成（`--wait` 超时、
作业可能仍在跑，或 `run` 有阶段停在 `waiting`），以及 **连接核对签署被拒绝**
（审计不过 / 出图失败 / 缺 metadata）。argparse 用法错误也返回 2。

### 7.1 判定与连接核对

- `smj report` 默认算物理判据并写进 `manifest.verdict`；缺 `TE.dat` 时登记
  `te_present`（engineering/error）且结论是 `inconclusive`（**证据缺失 ≠ 通过，
  也不等于否证 H0**），`report.md` 的「## 0.1 判定」会列出未过/未评估的判据。
- **默认自动**：`review.auto_connection_check: true`（默认）时，`electrode_build`
  收集完成即由管线自动跑一次表面审计并登记（`who="auto"`），不会因为"没人看图"
  把整条链卡住；审计不通过时**拒绝登记**，`smj review` 仍会 `blocked`。
- **闸门**：`review.connection_check: true`（默认）→ 没有有效记录（缺文件或指纹与当前
  `geometry.in` 不一致）时 `smj review` 直接 `blocked`。人工复看后可用
  `smj review --approve --who <名字>` 覆盖（图落在 `electrode_build/connection_check_viz.png`，
  审计结果落在 `electrode_build/connection_check_audit.json`）；`--no-figure` 表示已在
  别处看过图。想恢复"必须人工签"：设 `auto_connection_check: false`。
- **命令行临时开关**：`smj review --no-connection-check` 让本次复核不要求连接核对
  （缺记录/指纹失效只降为 info/warning），`--connection-check` 在 yaml 写 `false` 时
  强制开启；缺省沿用 `review.connection_check`。开关只作用于当次调用，不写回配置，
  也不进 `run_spec` 的协议哈希。
- **旧名兼容**：`review.require_a2` / `review.a2_deviation` / `review.surface_check` /
  `review.surface_check_deviation` / `review.auto_surface_check` 仍被接受（加载时提示改名），
  旧 run 目录里的 `a2_confirmation.json` / `a2_surface_audit.json` / `surface_check.json` /
  `surface_check_audit.json` 仍被读取；
  配置模板已改用 `connection_check` / `connection_check_deviation` / `auto_connection_check`。

### 7.2 运行溯源（run_spec / 哈希闭环）

- **第一次 `prepare`** 会在 run 根写 `run_spec.json`：完整解析后的配置（含模板默认值）、
  分子与电极库的 sha256、`protocol_sha256`（只覆盖会改变物理结果的分子/电极/输运子集）、
  SMJ/Python 版本、git commit/dirty、软件路径与本地可读二进制的 sha256。
- **配置漂移默认拒绝**：下一阶段 `prepare` 时若当前配置的 `protocol_sha256` 与
  `run_spec.json` 不同，直接报错；确认要重绑才加 `--force`，旧快照保留为
  `run_spec.previous.json`。
- **阶段哈希**：`state.json` 每阶段记录 `inputs_sha256` / `outputs_sha256` /
  `parent_stage` / `parent_fingerprint`；`aims_scf`、`aitranss_te`
  在 `prepare` 时校验上游 `geometry.in` 与电极 metadata 的 `surface_fingerprint`，
  不一致直接 `StageBlocked`（不会拿被换过的几何继续算）。
- **报告绑定**：`manifest.json` 带 `run_spec` 与 `provenance`（run_spec 文件 sha256、
  protocol、是否 reconstructed）；`report.md` 概览显示参数快照状态。旧 run 没有
  `run_spec.json` 时允许用当前配置重建，但标记 `reconstructed: true` 并告警。
- **复核门**：`smj review` 对缺 run_spec/重建只告警，对协议漂移记
  `provenance_drift` 并 `blocked`——结果不能按原参数解读时不会静默放行。
  两个例外都不改物理参数：**只差 `molecule.family`**（判据标定族，不参与物理计算）记
  `provenance_family_relabel`；**只差协议 schema 形状**（`run_spec.json` 记
  `protocol_schema`；新增/改名的记录键使哈希变化，但两侧都有的键取值全同、
  `PROTOCOL_CRITICAL` 一个不缺）记 `provenance_schema_evolved` warning——代码升级不必
  把历史 run 集体打成"不可解读"。物理参数一变或关键键缺失仍是 error。
  另一条独立的闸门：**分子几何内容**先按 `run_spec.json` 里的 `molecule.sha256`
  与当前文件哈希比对，同一路径下坐标一旦改动就直接记 `provenance_drift`
  （detail 带两个 sha 前缀），不会因为 `xyz`→`xyz_sha256` 两代键名不同而被误判成
  schema 演进。

### 7.3 收敛扫描（`scripts/scan_convergence.py`）

把"改一个参数重跑一遍"变成可复现的矩阵 + 对照表，**默认只写计划、不提交作业**：

```bash
# 1) 展开变体：每个组合一个目录、一份独立 smj.yaml（run.dir 固定为 run）
python scripts/scan_convergence.py plan smj.yaml \
    --axis electrode.nlayers=3,4,5 \
    --axis cluster.species_dir=/data/species/light,/data/species/tight \
    --out scan/
# 2) 真跑（必须 --yes；遵守 AGENTS.md「只跑用户点名的计算」）
python scripts/scan_convergence.py run scan/scan.json --yes
# 3) 汇总：读各变体 junction_analyze/results.json，逐轴、逐切片判是否收敛
python scripts/scan_convergence.py compare scan/scan.json --tolerance 0.05
```

- 轴写成 `section.key=v1,v2`，值会按 bool/int/float/None 还原；轴键先按 `smj.config`
  的字段表校验，不认识就直接报错（不会静默写进配置）；路径值在计划里绝对化，
  `scan.json` 记每个变体的真实 `run_dir`，写盘后还会回读自检。
- 基组"档位"不是配置键，用 `cluster.species_dir`（或 `paths.species_defaults`）这类
  目录轴表达；`paths.species_dir` 只是兼容别名，两个键同时给会报错。
- `compare` 读 `<run_dir>/junction_analyze/results.json`（旧数据回退 `analysis.json`），
  沿轴的取值顺序比较相邻两档的 `S_uV_per_K` / `log10_G_over_G0_like` /
  `ZT_electronic` / `k_e_pW_per_K`；多轴时只配对"其余轴取值相同"的行，逐切片比较，
  最大相对变化小于阈值才算 `converged`；缺结果、结果不可用或指标全空都记
  `incomplete`，**不把缺结果当 0**。
- 展宽 σ 属分析侧参数（`smj pdos --sigma`），不需要重跑 SCF/TE，所以不在扫描矩阵里。

## 8. 与上游 AgenticDFT 的关系

- **本仓自维护内核**：`src/smj/junction/` 是 junction 内核（含模板），自 0.2.2 起由本仓
  直接维护。它源自 `agenticdft/agenticdft/tools/junction/`（A 方案之前是逐字节副本），
  来源与本地变更史记在 `src/smj/junction/PROVENANCE.json`，但**不再与上游逐字节同步**：
  `scripts/sync_core.py` 已退役，不要再用它覆盖内核文件。
- **唯一的一致性契约是结果**：`demo/demo-reference.json` 记录 4 份示例谱的 G/G0、S、
  κe、ZT 与峰位；`python scripts/gen_demo_reference.py --check` 与
  `tests/test_demo_regression.py` 守住数值。若与 AgenticDFT 对同一 TE.dat 的结果分叉，
  用 `python scripts/crosscheck_agenticdft.py <TE.dat>` 对账（可选、不进 CI，CI 里没有主仓）。
- **台账**：`PROVENANCE.json` 的 `files` 是内核文件的 sha256（改过内核必须重跑
  `python scripts/gen_junction_provenance.py`，否则 `smj verify`/CI 报哈希过期）；
  `local_patches` 记录历史本地变更（id / commit / files / what / why / tests）；
  `repo_commit` 是生成时所在提交，**不参与**校验——否则每次 commit 都会让它变红。
- 主仓的 VASP 运行根（`agenticdft/work/`）与本包无关；本包的 run 产物落在
  `run.dir/run.name/`（默认是 `smj.yaml` 所在目录下的 `./<name>/`），不写主仓 `results.db`。

## 9. 已知限制

- **不自带 DFT 软件**：Gaussian 16、FHI-aims、aitranss 的二进制、许可证与 species 基组
  都不随包分发，需要集群/工作站已有并按 §4 配好路径；没有它们时只能 `dry-run` 或跑两个
  本地阶段。
- **Windows 本机执行需要 bash**：`run.mode: local` 依赖 `bash`；`slurm-ssh` 还需要可用的
  `ssh`/`scp`（Git-Bash 自带 OpenSSH）与目标集群免密登录。沙箱内网络受限时，请在宿主
  shell 里执行真跑。
- **结构确认默认自动、但可退回人工**：`electrode_build` 之后的连接核对默认由管线
  自动登记（`auto_connection_check: true`，`who=auto`）：审计（表面三元组/输运轴/接触）
  不通过时拒绝登记，`smj review` 仍 `blocked`；想强制人工看图再签，设
  `auto_connection_check: false`，再用 `smj review --approve --who <名字>` 签。签名绑定当前 `geometry.in` 的 sha256，并记录左右表面原子坐标；几何一改（哪怕只平移
  0.05 Å、表面索引集合没变）旧签名即失效，`smj review` 默认 `blocked`，只有显式
  `review.accept_deviation` 才降为 warning 并标注“人工豁免”。`review.connection_check`
  默认 `true`，没有有效签名时 `smj review` 为 `blocked`。
- **口径不可混**：Au60 与 20 Au/侧的绝对值不可直接比较；电子 `ZT_e` 不是总 ZT。
  热导链能测 `G_ph`，但它是**另一个器件**（周期板，见 §6.1），与电子 `G/G0` 不同构，
  `thermal.json` 的 `zt.available` 因此恒为 `false` 并写明原因——不要拿两个器件的
  分子结论互相拼成总 ZT。
- **电极厚度口径**：`electrode.au_per_electrode` 只在 `source: builtin` 生效
  （`null` = 完整 60 Au/侧；整数按完整原子层截断，30 → 35/侧；上限 30）。默认
  `source: library` 时该键被完全忽略，团簇大小由 `electrodes.library` 文件决定。
- **电极元素口径**：`electrode.element` 改变的是电极金属（tip 同元素），不只是标签；
  已登记 fcc 晶格常数的金属可直接换（Au/Ag/Cu/Ni/Pt/Pd/Al），其它元素必须显式
  `a_lattice`。`aims_scf`/`aitranss_te`/`self_energy`/`validate`/连接核对出图都从
  `electrode_metadata.json` 的 `electrode_element` 取元素（缺键按 Au），
  不要再用「Au」当所有体系的默认假设。
- **锚定化学不可混**：`2.40–2.44 Å` 是去质子**硫醇盐**的 Au–S **判据带**
  （Jackson 2009 一手锚点 S–Au = 2.42 ± 0.02 Å），且只按 Au 标定；
  `2.40–2.49 Å` 只是 `contacts.py` 的**建结目标/搜索带**（初猜排序用），不是判据带，
  两者不要混用；
  非 Au 电极只报告接触距离、不做该带有源判定。`-SCH3` 等**硫醚**的 S–Au 是配位接触
  （2.6–2.7 Å 也可能正常），`smj review` 按 S 的成键环境（S–C 数 / S–H）分流，不把该带
  套到硫醚上。**分子里没有 S 也可以建结**：`left_anchor`/`right_anchor` 写 N/Se 等任意
  元素符号或 1-based 序号；`strip_thiol_h: null`（默认）只对 S 锚自动去质子，其它锚保留 H
  （`true` = 强制去、`false` = 不去）。非 S 锚只报告锚–金属接触距离，不套任何键长带。
  默认 `anchor_au_s_target: null`、`pin_contact: false`：只用
  `au_tip_distance` 放初始几何，最终锚–金属键长/角度由 `aims_tip_relax` 的全原子弛豫决定；
  不要把 2.40 Å 当作目标键长或约束。
- **开壳层结（自由基）已接线，但口径有限**：`molecule.multiplicity > 1` 会一路带到
  `aims_tip_relax`/`aims_scf`（`spin collinear` + `fixed_spin_moment = M−1` +
  逐原子 `initial_moment` 初猜，`molecule.initial_moments` 会按结几何的分子前缀映射；
  不指定时退回 `default_initial_moment 0.1`；`molecule.fixed_spin_moment: false` 可
  释放总自旋约束做对照；非零 `charge` 把总电荷均分成 per-atom
  `initial_charge` 写进 `geometry.in`），
  SCF 产出 `alpha.aims`/`beta.aims`，`aitranss_te` 用 `$uhfmo_alpha`/`$uhfmo_beta`
  读双通道（远端回取时 `alpha.aims`/`beta.aims` 与 `mos.aims` 一样只取表头行）。
  `junction_analyze` 在 `results.json` 里另给 `spin` 块（每个通道的
  G/G0、E_F 处与积分后的自旋极化率）。**但**：总 G/S/ke/ZT 仍按
  `(Tα+Tβ)/2` 的双通道均值口径（与内核 Landauer 一致），且未考虑自旋轨道
  耦合、非共线磁构型、SOMO 与电极的强关联近藤效应；这些仍需人工判断，不要
  把平均透射当自旋分辨结论。另注意 `molecule.charge` 是整个模型的电荷（不是
  分子片段电荷），`fixed_spin_moment` 固定的是整个结——源文档建议把"释放约束
  后的解"作为对照，现在可用 `molecule.fixed_spin_moment: false` 直接跑；逐条对照
  与缺口清单见
  [自由基单分子结：实现对照与开发规划](docs/自由基单分子结-实现对照与开发规划.md)。
  展示层（2026-11-20）：`smj report` 在 `results.json` 有 `spin` 时输出 §4.1
  自旋通道小节（Tα/Tβ(E_F)、P(E_F)、P(L0)，并注明总 G/S 仍是双通道均值）；
  `smj review` 对 `multiplicity>1` 记 `open_shell_uncalibrated` 告警并在
  `review.json` 里给 `spin` 块——这两处都不改数值判定，只保证"这是开壳层、
  自旋占多少、判据未校准"不需要翻 results.json 才能看见。
  第一手证据（2026-11-20，默认关闭）：`transport.extra_outputs: [hirshfeld,
  spin_density]` 会让 `aims_scf` 的 `control.in` 写出 `output hirshfeld` /
  `output cube spin_density`，collect 核验 `spin_density.cube` 非空、
  `aims.out` 里确有 Hirshfeld 段——请求了却没产出会拒绝 `collected`。闭壳层
  要 `spin_density` 直接报错（没有共线自旋密度可写）。该关键字组合在 FHI-aims
  221103 上尚未实测，首次使用请在小系统确认；Gaussian 分支可用
  `molecule.stable_check: true` 加 `Stable=Opt`，`mol.log` 的
  `S**2 before annihilation … after …` 与 SCF 收敛状态会进
  `results.gaussian`（⟨S²⟩ 偏离 S(S+1) 超过 0.1 会在 notes 里点名自旋污染）。
- **自旋到底分布在哪（2026-11-21，S1）**：请求 `hirshfeld` 时，`aims_scf` 会把
  `aims.out` 的 `Hirshfeld spin moment` 逐原子解析到 `aims_scf/spin_partition.json`，
  `smj report` 的 §4.1 加"自旋分布（Hirshfeld，μB）"表：分子 / 接触区（距锚原子
  ≤3.0 Å 的电极原子）/ 电极的份额与占比。总磁矩守恒**不能**说明自由基还在分子上：
  `runs/ttf-cation` 实测总矩 1.0000 μB 里分子只占 0.0056（0.56%），`smj review` 会记
  `spin_localized_off_molecule`（info）。分区只陈述数字，不设阈值判定；Hirshfeld 段名
  在但解析不出 `| Atom N:` 行会拒绝 `collected`。
- **单列 TE 默认阻断（2026-11-21，S3）**：开壳层 TE.dat 只剩一列时，aitranss 已经把
  α/β 合并成均值，自旋极化不可恢复——`aitranss_te.collect` 现在直接 `StageError`，
  并提示检查 `$uhfmo_alpha/$uhfmo_beta`、SCF 是否真开壳层、nsaos 是否非零。外部/历史
  数据确知只能单列时显式设 `transport.allow_single_channel: true` 放行，结果里记
  `single_channel_allowed: true` 与 note。
- **自旋分辨热输运（2026-11-21，S2）**：`analyze_te_spin` 现在每通道也积 L1/L2，
  `results.json` 的 `spin` 块多出 `S_alpha_uV_per_K` / `S_beta_uV_per_K` /
  `S_total_uV_per_K` / `k_e_*` / `ZT_*` / `spin_polarization_S`；`smj report` §4.1
  在 P 行下加 Sα/Sβ 与 P_S。分通道公式与 `analyze_te` 共用
  `landauer.transport_moments()`（不复制公式），因此"总 S = 按 L0 加权的 Sα/Sβ"
  可被机械核对：engineering 判据 `spin_channel_consistency` 校验
  `G_total = 0.5(L0α+L0β)`（相对 1e-9）与 `S_total = (L0αSα+L0βSβ)/(L0α+L0β)`
  （相对 1e-6）；旧 `spin` 块（只有 L0）只标"不适用"，不误判故障。
  `runs/ifs-bs` 实测 Sα=+6.1560、Sβ=+6.1589 μV/K（相对差 4.8e-4），P_S=−2.4e-4，
  与 P(E_F)=−1.8e-5 同号且可忽略。
- **破缺对称（BS）单重态的证据（2026-11-21，S4）**：FHI-aims 单行列式 BS 近似不给
  ⟨S²⟩/NOON，能拿到的第一手证据是 Hirshfeld 局部自旋矩。`spin_partition.summarize`
  现在给出分子组内的正/负计数与净矩，正负共存、单个原子 ≥0.05 μB 且
  `|net| < 0.05·max_abs` 时标 `localized_antiparallel: true`（0.05 μB 下限挡住数值
  噪声）；`smj report` 的 §4.1 在 BS run 下加"破缺对称证据：分子内 +/− 自旋矩各
  N₁/N₂ 个，净矩 X μB"一行，`smj review` 记 info `broken_symmetry_local_moments`
  （只报数，不改任何判据结论）。要严格验证反平行请另跑 Gaussian 模型分子的 ⟨S²⟩
  交叉验证——默认不自动执行，也不自动跑闭壳层对照。
- **端到端集成回归（2026-11-21，S6）**：`tests/integration/test_radical_runs.py` 用真实
  `runs/ttf-cation` 资产重算 Hirshfeld 分区与双通道 S，断言六阶段齐全、通道数 2、
  |P(E_F)|<0.05、分子自旋占比 <1%，并在 run 副本上渲染报告、重跑 `review`（应出现
  `open_shell_uncalibrated`）。默认测试集用 `addopts = -m "not integration"` 排除它；
  要跑：`python -m pytest -m integration tests/integration -q`（本机 2 passed）。
- **开壳层协议自检与判据口径（2026-11-21，S5）**：`validate` 新增 engineering 检查
  `open_shell_protocol`（warn）：`multiplicity>1` 或 BS 单重态时核对 `spin` 块存在、
  `channels==2`、`transport.extra_outputs` 含 `hirshfeld`/`spin_density`，缺项逐条写进
  detail；闭壳层记"不适用"。该检查只进工程判定，不参与 `physics_verdict`。开壳层
  被判 `inconclusive` 且原因是 `no_applicable_gate` 时，`wording` 会明说"开壳层判据
  未校准；如需物理解读，需先用开壳层参考体系标定 `g_in_theory_band` / `seebeck_sign`"，
  避免把"没有闸门"读成"数据不足"。
- **`smj spin` 汇总子命令（2026-11-21，S8）**：`smj spin [RUN] [--json] [--temperature K]`
  只读汇总一条 run 的自旋通道（α/β 点数、Tα/Tβ(E_F)、P(E_F)/P(L0)、Gα/Gβ/Gtotal、
  Sα/Sβ/S_total、P_S）、Hirshfeld 分区与分子内 +/− 自旋矩、review 已有的自旋 code，
  并给出口径提醒。旧 run 的 spin 块没有分通道 S 时直接从 `TE.dat` 重算
  （`spin_source` 标明来源），分区侧车缺失就写"自旋分布: 未采集"。不写任何文件；
  需要新结论仍跑 `smj review`。

- **超时不是失败**：`--wait` 超时后状态留在 `waiting`，作业可能仍在跑；稍后重跑
  `smj stage <stage> --wait` 即可。
- **资源估算是申请下限**：132 原子锚点外推，体系明显偏离时应按实测回填，不要当保证值。

## 10. 故障排查（最短路径）

| 现象 | 先看 |
|---|---|
| `smj init` 说配置已存在 | 加 `--force`，或换目录 |
| 阶段卡在 `pending` | 先 `--prepare`；`smj status` 看 note |
| `--submit` 报"尚未 prepare" | 先 `smj stage <stage> --prepare` |
| 远端提交失败 | `smj preflight`；看 `cluster.env` 的 `g16root`/`aims_bin`/`aitranss_bin` |
| 轮询一直 `waiting` | 正常排队或超时；`smj stage <stage> --wait` 可重复执行 |
| `failed` 且 exit_code 非 0 | 看 `<stage>/slurm-*.out`、`aims.out`/`mol.log` 尾部 |
| 报告没更新 | `smj report --force`（默认跳过已存在文件） |
| 独立复核与管线不一致 | `smj analyze` 的 `independent_check` 结果，先查单位/口径 |

## 11. 许可证

MIT License，见 [LICENSE](LICENSE)。Copyright (c) 2026 AgenticDFT contributors。

