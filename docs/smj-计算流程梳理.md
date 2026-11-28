# SMJ 计算流程梳理（smj-toolkit 0.2.x）

> 依据 0.2.0rc1 工作树逐文件核对：`src/smj/config.py`、`state.py`、`runner.py`、`resources.py`、`preflight.py`、`review.py`、`report.py`、`provenance.py`，以及 `src/smj/stages/*` 与 `src/smj/junction/*`（inputs/electrode/electrode_library/tcontrol/self_energy/landauer/validate/criteria/analysis）。
> 2026-11-15 增补：步骤 4 由「本地构建 Au(111) 电极」改为「构建电极」——构建方式（`source: library|builtin`）与电极类型（`electrode.element` / `electrode.a_lattice`，fcc(111)）都是本阶段参数。
> 默认值一栏写的是 **不写 smj.yaml 时生效的值**；单位与量纲单独注明。

---

## 0. 一句话架构

SMJ 是一条「分子结构 → 透射谱」的端到端流水线：六阶段固定顺序、阶段之间只通过 `<run>/<stage>/` 下的**文件**耦合，状态全部记在 `<run>/state.json`；每个阶段都是 `prepare → submit → wait → collect` 四动作（本地阶段只有 `prepare → collect`），上游产物缺失时**拒绝猜测**（fail-closed）。

```
molecule.xyz
   │ ① molecule_opt      FHI-aims PBE/light 分子优化（或 Gaussian B3LYP/6-31G(d)）
   ▼ opt.xyz
   │ ② aims_tip_relax    分子 + 2 个电极金属 tip（默认 Au），全原子 BFGS 弛豫
   ▼ geometry.in.next_step
   │ ③ electrode_build   构建电极：tip 扩展为 fcc(111) 团簇（默认 Au(111)；library / builtin，本地阶段）
   ▼ geometry.in + electrode_metadata.json
   │ ④ review（默认自动记录；可人工复核）  结构闸门 + 连接核对（原 A2）
   ▼
   │ ⑤ aims_scf          PBE 静态 SCF（KS_method serial）
   ▼ aims.out / omat.aims / mos.aims / basis-indices.out
   │ ⑥ aitranss_te       NEGF 透射（自能 σ_near=0.05 / σ_far=0.10 Ha）
   ▼ TE.dat（+ te.out）
   │ ⑦ junction_analyze  Landauer 积分（本地阶段）
   ▼ results.json
   │ ⑧ smj report        manifest.json + report.md（物理判据 verdict）
   ▼ 报告 + 判定
```

- **远端阶段**（`needs_submit=True`）：molecule_opt、aims_tip_relax、aims_scf、aitranss_te
- **本地阶段**：electrode_build、junction_analyze（`needs_submit=False`，prepare 即产出，collect 只校验/汇总）
- 旧阶段 `au_relax` 已于 2026-11-13 移出主链路，只供历史 run 或显式 `smj stage au_relax` 使用。

---

## 1. 步骤总览（先看这张）

| # | 步骤 | 命令 | 软件 | 目的 | 关键默认参数 | 产物 |
|---|------|------|------|------|--------------|------|
| 0 | 初始化配置 | `smj init` | — | 生成 smj.yaml 全量默认配置 | `mode=dry-run`，`--xyz molecule.xyz` | `smj.yaml` |
| 1 | 提交前预检 | `smj preflight` | — | 结构/物种/资源/集群/二进制五类检查，error 级阻断提交 | 不触网；SSH 超时 15 s | 终端/`--json` 报告 |
| 2 | 分子几何优化 | `smj stage molecule_opt` | FHI-aims（默认）或 G16 | 得到单分子平衡几何 | PBE/light，charge 0，multiplicity 1，力阈值 `1.e-2` eV/Å，vdW=Hirshfeld | `opt.xyz` |
| 3 | 加金属 tip 并弛豫 | `smj stage aims_tip_relax` | FHI-aims | 每个锚定原子接一个电极金属 tip（默认 Au），得到 锚–金属 接触几何；锚定原子可以是 S（默认）、N、Se 等 | 放置距离 2.35 Å，`strip_thiol_h=null`（仅 S 锚自动去质子），`pin_contact=false`，全原子自由弛豫；tip 元素 = `electrode.element` | `geometry.in.next_step` |
| 4 | 构建电极 | `smj stage electrode_build` | 纯 Python（本地） | 按 `source` 选构建方式、按 `element` 选电极金属，把 tip 扩展为完整 fcc(111) 电极团簇，冻结结几何 | 构建方式 `source=library`（默认）/`builtin`；电极类型 `element=Au`、`a_lattice=null` | `geometry.in` + `electrode_metadata.json` |
| 4.5 | 结构复核 + 连接核对 | `smj review [--approve --who] [--connection-check\|--no-connection-check]` | — | 确认建结几何没被换错：结构闸门 + 连接核对（原 A2，默认自动登记） | `review.connection_check=true`、`connection_check_deviation=0.15`、`auto_connection_check=true`；CLI 可临时覆盖 | `review.json`、`connection_check.json`（+ `connection_check_audit.json`、`connection_check_viz.png`） |
| 5 | 静态 SCF | `smj stage aims_scf` | FHI-aims | 为 aitranss 准备冻结几何的单点波函数 | PBE，`KS_method serial`，sc_accuracy_rho 1e-5 / eev 1e-3 / etot 1e-6，sc_iter_limit 500；开壳层加 `spin collinear` + `fixed_spin_moment=M−1` | `aims.out`、`omat.aims`、`mos.aims`（闭壳层）或 `alpha.aims`+`beta.aims`（开壳层）、`basis-indices.out` |
| 6 | 透射谱 | `smj stage aitranss_te` | Aitranss | NEGF + 自能，算 T(E) | 能量窗默认 `[-2.0, 0.005, 2.0]` eV 相对 E_F（即 E_F ±2 eV，801 点）；prepare 时按 SCF 的 efermi 换成 Ha 绝对值；σ_near 0.05 / σ_far 0.10 Ha | `TE.dat`（+`te.out`） |
| 7 | Landauer 后处理 | `smj stage junction_analyze` | 纯 Python（本地） | 积分得 G/G0、S、k_e、ZT_e | `temperature=300 K` | `results.json` |
| 8 | 报告与判定 | `smj report` | — | 汇总 + 物理判据 verdict | 缺 TE.dat 判 `inconclusive` | `manifest.json`、`report.md` |
| 附 | 独立复算 | `smj analyze` | — | 独立解析器+Simpson 复核 + 可选出图 | 容差：log10G 0.1、S/ZT_e/k_e 10%、T(E_F) 5% | `junction_analyze/analysis.json` |

全链一条命令：`smj run --from molecule_opt --to junction_analyze`（`--run` 为 stage 子命令默认动作）。

---

## 2. 逐步详述

### 步骤 0：`smj init` —— 生成配置

- **命令**：`smj init [PATH] [--name NAME] [--xyz PATH] [--mode dry-run|local|slurm-ssh] [--force]`
- **目的**：把一份「默认就是模板」的 `smj.yaml` 写到目标目录，用户只改物理参数，不改结构。
- **默认**：`path="."`；`--xyz molecule.xyz`；`--mode dry-run`（注意：`RunConfig` 数据类里 `run.mode` 默认是 `"slurm-ssh"`，但 `init` 模板写的是 `dry-run`）；`--name` 缺省取目录名的小写 slug；同名文件存在时拒绝，`--force` 覆盖。
- **写出的关键默认**：`run.dir "."`、`run.clean false`、`run.keep_going false`；`molecule.opt_engine aims`、`opt_functional pbe`、`opt_basis light`、`opt_force_threshold 1.e-2`；`electrode.au_per_electrode 30`、`source library`、`element Au`、`a_lattice null`、`au_tip_distance 2.35`、`pin_contact false`、`strip_thiol_h null`、`vdw true`、`nlayers null`；`transport.temperature 300.0`、`energy_window [-2.0, 0.005, 2.0]` + `energy_window_unit eV`（相对 E_F 的 ±2 eV，步长 5 meV）；`review.connection_check true`、`connection_check_deviation 0.15`、`auto_connection_check true`。

> 💡 能量窗口径（2026-11-15 起）：`transport.energy_window` 默认 `(-2.0, 0.005, 2.0)` **eV 且相对费米能级**（`E − E_F ∈ [−2, +2] eV`，801 点）。`aitranss_te` 在 prepare 时从 `aims_scf/aims.out` 读 efermi，平移到 `tcontrol` 需要的 Ha 绝对 `$ener/$estep/$eend`；缺 efermi 直接 `StageBlocked`（不猜）。换金属/换分子后 E_F 会移动，相对窗口不跟着漂。旧 Ha 绝对三元组（跨度 ≤0.1 Ha、量级 <2 Ha）仍按旧口径解释并给迁移提示，也可显式写 `energy_window_unit: Hartree`；`transport.ener/estep/eend` 恒为 Ha 绝对值覆盖。

---

### 步骤 1：`smj preflight` —— 提交前预检

- **命令**：`smj preflight [--stage STAGE] [--no-remote]`
- **目的**：在花钱排队之前，把「结构错、缺物种、资源超限、集群连不上、二进制/许可不在」尽量暴露在本地。`error` 级 issue 阻断提交，`warning/info` 不阻断。
- **默认行为**：`transport=None` 时**完全不触网**；slurm-ssh 模式只记 warning「未执行远端检查」（`require_remote=True` 时升为 error）；`dry-run` 不核对二进制与许可（只提示「提交前在目标机器再跑一次非 dry-run 预检」）；SSH/远端探测超时 15 s。

检查项与阈值：

| 检查 | 级别 | 规则/默认 |
|------|------|-----------|
| 分子结构 | error | 路径未配置/文件不存在/不可解析 → error；可读时记 info（原子数、元素） |
| 分子重叠 | error | 元素对键长下限 = **0.7 × 共价半径和**，3 Å 预筛（`PAIR_CUTOFF`）；>2000 原子跳过 O(N²) 并与结构闸门兜底 |
| species 目录 | error | `cluster.species_dir` 或 `paths.species_defaults` 必须存在；逐元素核对 species 文件齐备（`*_<El>_default` 等 5 种命名） |
| 建结几何闸门 | error/warning | 若 `electrode_build/geometry.in` 已存在，跑 `validate_junction_geometry` 的 7 项检查（见 §6.1） |
| 资源估算 | error | 四个远端阶段逐一估算；超 72 h / 128 GiB 封顶报 `resource_capped`，超预算报 `resource_budget`；electrode_build 记「本地阶段不占预算」 |
| 集群配置 | error | slurm-ssh 必须有 `host / partition，远端工作根为 ~/smj-runs` |
| SSH 连通 / 远端目录 | error | `ssh true`、`test -d`、`test -w` 三步 |
| 二进制/许可 | error | FHI-aims 可执行、aitranss 二进制；Gaussian 许可仅在 `opt_engine=gaussian` 时检查（dry-run 全跳过） |

---

### 步骤 2：`molecule_opt` —— 分子几何优化

- **阶段名**：`molecule_opt`（旧别名 `gaussian_opt` / `g16` / `mol` / `opt`）；依赖：无（首段）。
- **目的**：把用户给的 `.xyz` 优化到平衡几何；**去质子化不在这里做**，留给 `aims_tip_relax`（`strip_thiol_h`）。
- **执行**：远端（`sub_aims.sh` 或 `sub_g16.sh`）。
- **默认参数**：

| 配置键 | 默认 | 说明 |
|--------|------|------|
| `molecule.opt_engine` | `aims` | `aims`（PBE/light，与 tip/SCF 同栈）或 `gaussian`（B3LYP/6-31G(d)，需 G16 许可） |
| `molecule.charge` | `0` | 非零时把总电荷均分成 per-atom `initial_charge` 写进 geometry.in：FHI-aims 要求 `charge` 等于这些行之和，而 221103 不认 `override_initial_charge_check`（实测 parse_control 停算） |
| `molecule.multiplicity` | `1` | `>1` 时写 FHI-aims `spin collinear`（`spin` 只接受 `none`/`collinear`，不接受数字）并用 `fixed_spin_moment = multiplicity-1`（=2S，μB；与多重度不是同一个数）把总自旋在弛豫中固定住 |
| `molecule.initial_moments` | `null` | 开壳层初始磁矩 `{1-based 原子序号: μB}`，写成 geometry.in 的 per-atom `initial_moment` 子 tag（独占一行）；FHI-aims 在 collinear 下要求初始自旋密度，未设只给警告。闭壳层设置会报错。结级 `aims_tip_relax`/`aims_scf` 会把分子坐标系的序号映射到结几何（去质子后自动左移） |
| `molecule.fixed_spin_moment` | `null`（依自旋态选择） | 普通开壳层的结级计算默认固定 (N_α−N_β=M−1)，BS 默认释放。BS 显式设为 `true` 时，FHI-aims 分子优化、尖端弛豫和整结 SCF 均固定净磁矩为零，同时允许两个中心保留反号局域矩。 |
| `molecule.opt_functional` / `opt_basis` | `pbe` / `light` | `opt_basis` 只记录；实际基组由 species 目录决定，`smj species` / prepare 写出的 `species.json` 记录真实档位与逐文件 sha256 |
| `molecule.opt_force_threshold` | `1.e-2` | 写进 `relax_geometry bfgs`；单位 **eV/Å**，判据是**最大力分量**（日志 `Maximum force component ... eV/A`），不是 Hartree/Bohr；合法范围 `(0, 0.1]` |
| `electrode.vdw` | `true` | → `vdw_correction_hirshfeld`；只改能量/力（几何），不改固定几何的 T(E) |

- **写出的输入**（aims 分支）：`geometry.in`、`control.in`（`control_relax.in` 模板：`xc pbe`、`spin none|collinear`、`charge`（≠0 时另有 per-atom `initial_charge`）、开壳层另有 `fixed_spin_moment`、`relativistic atomic_zora scalar`、**显式 SCF 精度块**（`occupation_type gaussian 0.01`、`mixer pulay`、`n_max_pulay 10`、`charge_mix_param 0.2`、`sc_accuracy_rho 1E-5`、`sc_accuracy_eev 1E-3`、`sc_accuracy_etot 1E-6`、`sc_accuracy_forces 5E-3`、`sc_iter_limit 1000`）、`relax_geometry bfgs 1.e-2`、`aggregated_energy_tolerance 1.0e-1`、`output cube eigenstate HOMO-1/HOMO/LUMO/LUMO+1`、`output dipole`、species 块）、`sub_aims.sh`、`input_summary.json`。
- **SCF 精度为什么必须显式写**：旧模板只写几何相关的键，FHI-aims 会打印 `No s.c.f. convergence criteria (sc_accuracy_*) were provided`，自行挑 `sc_accuracy_rho`（实测 `0.300000E-05`，随版本/体系变）并把 `sc_accuracy_forces` 设为 `not_checked`——**BFGS 吃的力没有独立精度保证**。现在与 `control_scf.in` 同源给出，`sc_accuracy_forces 5E-3 eV/Å` 取几何阈值（`1.e-2`）的一半，`sc_iter_limit` 保留 1000（FHI-aims 自带默认，不收紧成 SCF 模板的 500）。这些值同时写进 `input_summary.json` 的 `scf` 段，collect 再从 `control.in` 读回、写进 `results.json`。
- **产物**：`opt.xyz`（下游统一契约）、`aims.out`、`geometry.in.next_step`（aims 还会写出 4 个 `cube_*_eigenstate_*.cube`）；Gaussian 分支为 `mol.gjf` → `mol.log` → `opt.xyz` + `gaussian_summary.json`。
- **验收（collect 闸门，fail-closed）**：
  - aims 收敛：`aims.out` 必须含 `Have a nice day`，**并且**必须出现 `Present geometry is converged`。FHI-aims 达到弛豫步数上限时仍会正常退出（第 1 步的 `Have a nice day` 与 `geometry.in.next_step` 都照写），所以"正常结束"不等于"收敛"，必须显式核对，否则 StageError（消息带最后一步最大力与阈值）。
  - SCF 收敛（aims）：日志里若出现 "self-consistency cycle … not … converged"，说明那一步的力没有自洽保证，直接 StageError；并核对 `Self-consistency cycle converged.` 的出现次数不少于弛豫步数（demo 实测 106 / 105），不足只记 note（日志可能被截断），不阻断。
  - 电荷/多重度奇偶校验（prepare 阶段，两种引擎共用）：电子数 = ΣZ − charge；偶数电子必须奇多重度、奇数电子必须偶多重度，否则直接 StageError。写错时 SCF 会为"不存在的电子态"收敛出漂亮结果，链上没有别处能发现。
  - 几何合理性（两种引擎共用 `smj.checks.check_geometry`）：最近原子对 ≥ 0.7 Å；输入几何 vs 优化后几何的**键图对比**——重原子键消失/原子解离 → StageError，含 H 的键消失（质子迁移）→ note，元素组成变化 → StageError；未知共价半径的元素只做重叠检查。检查在**落盘 `opt.xyz` 之前**，不通过就不产出。
  - gaussian：必须 `Normal termination`，否则先删 `opt.xyz` 再报错；虚频只记 note（需人工确认）。
  - 锚定轨道权重（aims，`smj.cube`）：若阶段目录有 cube，按 `control.in` 的请求顺序（HOMO-1/HOMO/LUMO/LUMO+1）解析，给出**半径 1.6 Å 球内 |ψ|² 占全箱的比例**——权重高 = 该前线轨道在锚定原子上有振幅（有耦合通道），接近 0 = 轨道不经过锚定原子。只算占比，因此不受 FHI-aims "内容单位 1/Å³、体素单位原子单位" 警告影响；文件有负值按带符号 ψ 平方处理，全非负按密度处理。没有 cube（旧 run / gaussian）只记一条 note，不报错。
  - 结果写 `results.json` 的 `molecule_opt` 段：`engine`、`normal_termination`、`atom_count`、`converged`、`max_force_ev_per_A`、`relax_steps`、`electron_count`、`geometry_check`（problems / warnings / min_distance_A / rmsd_A / lost_bonds / bond_count…）、`frontier_orbitals`（`homo_ev` / `lumo_ev` / `gap_ev`，eV；分数占据时 `fractional_occupation=true`）、`frontier_orbital_weights`（每项 `orbit` / `eigenstate` / `anchor_fraction` / `anchors[]`）、`scf`（从 `control.in` 读回的 `sc_accuracy_*` / `sc_iter_limit`，加 `converged_steps` / `relax_steps`；旧 run 无这些行时不伪造）、`species`（`actual_tier` / 逐元素 `tier` / `basis_functions` / `sha256`；缺 `species.json` 的旧 run 只记 note）；aims 分支另写 `software` 段。
- **资源默认**：借用 `aims_tip_relax` 线性模型（60 s + 10 s/atom 是 resources 表里的分子项，但实际操作按 tip_relax 标定）；墙钟下限 900 s，`--time = 预期 × 1.9`，`ntasks 1`。

> 备注：FHI-aims 分支**没有虚频检查**（未提供一键 freq）。几何可靠性现在由四层把关：①`Present geometry is converged`（弛豫收敛，collect 强制）；②SCF 精度显式 + SCF 未收敛拒绝（collect 强制）；③`smj.checks` 几何合理性检查（重叠 / 断键 / 解离，collect 强制）；④`smj review` 人工确认。Gaussian 分支的虚频只提示、不阻断。

- **化学态对照臂（`smj variants`，`smj.variants`）**：器件里的锚定原子是 S–Au 成键的硫醇盐/硫醚，而本阶段优化的是中性分子（硫醇的 S–H 还在，去质子化在 `aims_tip_relax`）。同一条 2CN-656 硫醚在两条 demo 里的 S–Au 接触分别是 2.6877 Å / 141.5°（tip 落在硫醚上）与 2.4267 Å / 111.4°（硫醇盐），差 0.26 Å / 30°——tip 落点对初猜很敏感，只靠 `tip_relax` 局部弛豫未必能回到真实接触几何。因此提供两条**对照臂**做量化：

  ```bash
  smj variants --xyz molecule.xyz --dir variants/      # 生成 thiolate.xyz / aucap.xyz + variants.json
  ```

  - `thiolate`：去掉锚定元素上的 H（S–H ≤ 1.70 Å），`charge -= 1/个`，电子数与多重度不变（去的是 H⁺）；
  - `aucap`：把 S–H 的 H 换成沿原 S→H 方向的 Au（默认 2.40 Å），电荷不变——局部模拟"已与电极成键"的化学态；
  - `variants.json` 记 `source_sha256`、每条臂的 `atom_count / electron_count / removed_h_indices / note / config_hint`（charge/multiplicity 直接照抄进配置），以及 `caveat`：气相 / 0 K / 无溶剂、裸阴离子在 PBE/light 下 HOMO 可能高于真空能级、单原子 Au 不是真实电极；
  - 两条臂各跑一遍 `molecule_opt`，比较 C–S 键长、gap 以及 `aims_tip_relax` 后的锚–Au 接触几何；差异就是"质子化态 / 初猜"对结果的影响量级。**只生成输入与清单，不做判定**。

- **species 审计（`smj species`，写进每次 prepare）**：`molecule.opt_basis` 只是配置标签，真正生效的是 `cluster.species_dir` 里的 species 文件。prepare 会逐元素记录**文件头声明的档位**（`Suggested "light" defaults`）、生效基组行数、最高角动量与**逐文件 sha256**，写到阶段目录的 `species.json`；collect 原样搬进 `results.json` 的 `species` 段，标签与实际档位不一致时在 notes 里给警告（"以 species 文件为准"）。这样"这次算的到底是 light 还是 tight"是可核对的，而不是靠配置里的一个词。

  ```bash
  smj species --dir $SMJ_SPECIES --xyz molecule.xyz --label light --out species.json
  ```

- **基组敏感性怎么测**（离线也能做，不需要新代码）：
  1. `smj species --dir light/ --xyz molecule.xyz --label light` 与 `--dir tight/ --label tight`，确认两次的 `sha256` 不同、档位一致；
  2. 同一分子、同一 `control.in`（除 species 块外一致）各跑一次 `molecule_opt`，`--run-dir` 分开；
  3. 比较 `results.json` 的 `frontier_orbitals.gap_ev` 与 `opt.xyz` 的 C–S / 锚–锚键长（`smj checks` 的 `anchor_contacts` 也能用）；
  4. 把差值记成"基组 + 积分网格"的合并不确定度——**注意 tier 换档同时改基组与网格，不能把二者分开归因**，除非另做网格收敛测试。

### 步骤 2.5：构象自动分支（`smj ensemble`，不属于阶段链）

柔性分子（烷烃二硫醇 HDT/ODT/DDT、2CN-656 的甲硫基、azulene 取代基）换一个构象，锚–锚距离可能差 0.3–1 Å，gap 与 T(E) 都会变；只优化"手画的那一个构象"等于把构象误差当成方法误差。

```bash
smj ensemble --smiles "Sc1cccc(S)c1" --dir ensemble/ --keep 3   # 生成 conf_01..03.xyz + ensemble.json + run_ensemble.sh
sh ensemble/run_ensemble.sh                                       # 每个构象一个隔离 --run-dir（run/conf_XX）
smj ensemble --dir ensemble/ --report                             # 汇总 gap / 最大力 / 构象能量的极差
```

- 生成阶段：ETKDGv3 → MMFF94s → 按能量升序 + 重原子骨架 RMSD 去重（与 `smj conformers` 同一实现）；`run_ensemble.sh` 里路径相对脚本目录，便于整目录拷到集群；`--run-dir` 覆盖 run 目录，各臂的 `run_spec.json` 与产物互不覆盖。
- 汇总阶段：读各臂 `molecule_opt/results.json`，给 `gap_ev` 与 `max_force_ev_per_A` 的**极差**；`gap` 极差 <0.05 eV 可认为该分子对构象不敏感，>0.1 eV 则必须把构象分布（或至少多个极小）纳入结论。
- 口径：MMFF94s 只决定"哪些构象值得算"，不是能量判据；构象能量差在 kcal/mol 量级时应给 Boltzmann 权重，而不是只取最低能。

---

### 步骤 3：`aims_tip_relax` —— 加电极金属 tip 并全原子弛豫

- **阶段名**：`aims_tip_relax`（别名 `tip` / `relax`）；依赖：`molecule_opt`；远端（`sub_aims.sh`）。
- **目的**：把分子两端锚定原子（S/N/Se 任意元素）各接一个**电极金属** tip（默认 Au，随 `electrode.element`），**全原子第一性弛豫**，让 锚–金属 接触由力场（DFT）自己找位置，而不是人为固定距离。
- **输入**：上游 `molecule_opt/opt.xyz`（可 `--geometry` 覆盖）。
- **默认参数**：

| 配置键 | 默认 | 说明 |
|--------|------|------|
| `electrode.element` / `a_lattice` | `Au` / `null` | tip 与最终电极用同一元素；`a_lattice=null` 用内置 fcc 晶格常数（见 §3.3） |
| `electrode.left_anchor` / `right_anchor` | `S` / `S` | 任意元素符号（分子里出现**恰好 2 次**，第一个=左，第二个=右；无 S 分子就写 N/Se 等），也可写 1-based 序号 |
| `electrode.au_tip_distance` | `2.35` Å | tip 初始放置距离；校验 1.0–5.0 |
| `electrode.anchor_au_s_target` | `null` | 非空时直接把 tip 放到该 锚–金属 距离（1.5–4.0 Å）；生产一般留空，让弛豫决定 |
| `electrode.pin_contact` | `false` | `false` = 分子+2 tip **全部自由**；`true` = 冻结两锚定原子 + 两 tip（距离扫描/固定接触口径） |
| `electrode.strip_thiol_h` | `null` | 三态：`null` = 仅当锚定原子是 S 才移除 S–H（<1.70 Å，必须恰好 0/2 个），N/Se 锚保留 H；`true` = 强制去（每个锚原子恰好 1 个 H，-NH2 会 fail-closed）；`false` = 不去；去质子前原几何存 `opt_neutral.xyz` |
| `electrode.vdw` | `true` | Hirshfeld 色散 |
| `molecule.charge` / `multiplicity` | `0` / `1` | |
| `opt_force_threshold` | `1.e-2` | 同上（eV/Å） |

- **放置规则**：默认沿「左锚–右锚」连线向外（`add_tip_au` 默认方向 = 本锚减对侧锚），距离 = `au_tip_distance` 或 `anchor_au_s_target`；落盘后**立即回测** 金属–S 距离，偏差 >1e-6 Å 直接 RuntimeError（拒绝悄悄错位）。
- **产物**：`geometry.in.next_step`（BFGS 最终几何）、`aims.out`、`geometry.in`、`input_summary.json`（含 `tip_element` 与 `anchor_element`）。
- **验收**：优先 `geometry.in.next_step`，缺失才回退 `geometry.in` 并记 note；collect 计算**弛豫后真实 金属–S 距离** `relaxed_au_s_A`（tip 定义为文件末尾 2 个原子，元素 = `tip_element`）。
- **资源默认**：180 s + 30 s/atom；内存 8 GiB + 0.10 GiB/atom；`ntasks 4`；墙钟下限 900 s。
- **备注**：2.35 Å 只是初始放置值，最终 锚–金属 由弛豫得到；`preflight/review` 的 `au_s_bond` 判据只对**S 为锚的去质子硫醇盐 + Au 电极**用 2.40–2.49 Å 带（非 Au 电极或非 S 锚只报告接触距离），硫醚（-SCH3）是配位接触（2.6–2.7 Å），不适用该带。

---

### 步骤 4：`electrode_build` —— 构建电极（本地）

- **阶段名**：`electrode_build`（旧名 `au60_build`；更早 `au20_build`/`auN_build` 归一化）；依赖：`aims_tip_relax`；**本地阶段**（`needs_submit=False`，prepare 即产出）。
- **目的**：把 2 个 tip 沿 tip→linker 方向扩展成完整的 **fcc(111) 电极团簇**；**只扩展、不改 金属–S 距离/角度**，此后结几何冻结。**怎么构建**（`source`）与**构建什么电极**（`element`/`a_lattice`）都是本阶段参数。
- **默认参数**：

| 配置键 | 默认 | 说明 |
|--------|------|------|
| `electrode.source` | `library` | **构建方式**：`library` = FHI-aims `electrodes.library` 预定义团簇；`builtin` = 本地构建 fcc(111) 三角岛（`Electrode111`；非 Au 走 `ElectrodeFCC`） |
| `electrode.element` | `Au` | **电极金属**（fcc(111) 单晶面）：Au / Ag / Cu / Ni / Pt / Pd / Al 可直接用；其它元素必须显式给 `a_lattice`。tip 同元素 |
| `electrode.a_lattice` | `null` | fcc 晶格常数（Å，2.0–6.0）；`null` = 用内置登记值（Au 4.07823、Ag 4.08530、Cu 3.61496、Ni 3.52360、Pt 3.92360、Pd 3.89070、Al 4.04950）；未登记元素必须显式给出 |
| `electrode.library_cluster` | `fcc111.au.clusters.37_atoms/fcc111.cluster.natoms37.x1.aims` | 37 原子，头部 `-nlayers=3`（手册 §5.3.1 同族示例）；库团簇元素必须与 `element` 一致，否则 fail-closed |
| `electrode.library_dir` | `null` | 也可用环境变量 `SMJ_ELECTRODE_LIBRARY`；否则依次搜 `/data/apps/src/fhi-aims*/src/aitranss/electrodes.library`、`/opt/fhi-aims*`、`~/fhi-aims*`、`~/software/fhi-aims*`；slurm-ssh 下还会 `ssh cat` 远端取回并缓存到 `<run>/electrodes/` |
| `electrode.nlayers` | `null` | `null` = 取库文件头部建议值；与库建议不一致会 note 警告 |
| `electrode.au_per_electrode` | `30` | **仅 builtin 生效**；默认 `source=library` 时**完全忽略**（显式写会 warning）。上限 30（超过报错）。`null` → 全量 60 原子/侧（body 56 + 4 extras）；数值 N → 逐层累加到 ≥N 的最完整层数（N=30 → 35 原子/5 层），不残缺表面层 |
| 内置团簇晶格 | `a` 由 `element`/`a_lattice` 决定 | 三角层 1+3+6+10+15+21 = 56，+4 = 60 原子/侧；`Δz=a/√3`、`Δy=Δz/√2`、`Δx=Δz·√3/(2√2)`（Au：Δz=2.35455、Δy=1.66492、Δx=1.44186 Å） |

- **产物**：`geometry.in`（结几何）、`electrode_metadata.json`（左右 tip 索引、linker 索引、六表面索引 `surc/lsurx/lsury + rsurc/rsurx/rsury`、左右背层索引、左右电极索引集、`final_atom_count`、`surface_rule="det-v1"`、`surface_fingerprint`、库团簇 provenance；**非 Au 电极另写 `electrode_element` / `electrode_a_lattice`**，Au 默认保持旧 metadata 逐字节不变）。
- **闸门**：结构里必须**恰好 2 个** tip 元素（= `electrode.element`）；linker = 距 tip 最近的**非电极元素、非 H** 原子；`source=library` 时库团簇元素必须与配置一致（不一致 StageError）；collect 时 metadata 的 `final_atom_count` 必须等于实际原子数；缺 metadata 指纹会 note「审阅确认无法绑定几何」。
- **资源**：本地，不占集群预算。
- **备注**：库文件**不进 wheel**，run_spec 记录其 sha256 + 相对路径；库电极（37 原子）与内置 Au60、与旧的 20 Au/侧口径**不可直接横比**。换电极金属会改变 `protocol_sha256`（`element`/`a_lattice` 属于物理输入），旧 run 的协议漂移检查会报 `provenance_drift`。

---

### 步骤 4.5：`smj review` —— 结构复核 + 连接核对（原 A2，默认自动登记）

- **命令**：`smj review [--approve --who NAME] [--force] [--stage DIR] [--no-figure] [--accept-deviation FRAC] [--connection-check|--no-connection-check]`
- **目的**：三段式复核，确保错误结构不进入 SCF：① 结构几何闸门（`validate_junction_geometry` 7 项）；② 连接核对（原 A2，默认开启、默认自动登记）；③ 溯源一致性。
- **锚定接触报告**（只报告、不判定）：结构闸门之外，`review` 还对每个贴着金属的锚定原子（S/N/Se）给出 `M–X` 距离与 `C–X–M` 角，写 `review.json` 的 `structure.anchors` 并打一条 info。硫醇盐在 Au(111) 上通常 ≈2.4–2.5 Å、≈105–115°；明显更长（≳2.6 Å）或接近线形（≳140°）说明锚定轨道与表面耦合弱（2CN-656 的硫醚锚实测 2.69 Å / 141°，与 G/G0=0.00132 一致）。合格区间随锚定化学与位点变化，所以不进闸门，交人工判读。
- **默认（自动）**：`review.auto_connection_check=true` 时，`electrode_build` collect 完即由管线自动跑一次表面审计并登记 `electrode_build/connection_check.json`（`who="auto"`）；审计不通过则**拒绝登记**，`smj review` 仍会 `blocked`。想恢复"必须人工看图再签"，设 `auto_connection_check: false`。
- **闸门**：`review.connection_check=true`（默认）→ 没有有效记录（缺文件/指纹与当前 `geometry.in` 不符）时 **blocked**（退出码 1）；`connection_check_deviation=0.15`（表面索引偏差容差）；`--approve --who` 才人工写 `electrode_build/connection_check.json`（绑定 `geometry.in` 的 sha256 与签署人；已有有效记录默认幂等不覆盖，`--force` 重签）。CLI 可临时覆盖：`--no-connection-check` 关闭该项要求（缺记录/指纹失效降为 info/warning），`--connection-check` 强制开启；只作用于本次调用，不写回配置、不进协议哈希。
- **判定**：`passed | warning | blocked`，`ok = status != "blocked"`；`review.json` 落 run 根目录。
- **签署流程**（`--approve`）：先跑 `a2_viz.surface_audit`（表面索引/输运轴/接触审计，写 `electrode_build/connection_check_audit.json`）→ 不通过直接拒绝；默认渲染 `electrode_build/connection_check_viz.png`（**无图拒绝签署**，`--no-figure` 表示已在别处看过图）；再写 `electrode_build/connection_check.json`（who、when、geometry sha256、表面索引）。自动登记的失败回退：matplotlib 不可用时改 `--no-figure` 路径重试，仍失败则记 warning 交给闸门。
- **偏差算法**：`deviation = 1 - Σ|交集| / Σ max(|集合|)`（左右两侧加权，0..1）。≤0.15 → info；≤`--accept-deviation`（显式豁免）→ warning 且 `accepted=True`（报告须标注为人工豁免）；否则 error。
- **旧名兼容**：配置键 `review.require_a2` / `review.a2_deviation` / `review.surface_check*` 仍接受（加载时提示改名）；旧 run 的 `a2_confirmation.json` / `a2_surface_audit.json` / `surface_check.json` / `surface_check_audit.json` 仍被读取；CLI 旧开关 `--surface-check/--no-surface-check`、`--a2/--no-a2` 仍生效。
- **溯源检查**：`run_spec.json` 缺失 → warning；`reconstructed=true` → warning；`protocol_sha256` 与当前配置不符 → **error provenance_drift**（run 开始后物理相关配置变了，结果不能按原参数解读）。唯一例外：差异只在 `molecule.family`（判据标定族，不参与物理计算）→ `warning provenance_family_relabel`，按当前标注重算判据，不阻断。

---

### 步骤 5：`aims_scf` —— 静态 SCF

- **阶段名**：`aims_scf`（别名 `scf` / `aims`）；依赖：`electrode_build`；远端（`sub_aims.sh`）。
- **目的**：在**冻结的结几何**上做一次静态 SCF，产出 aitranss 需要的波函数、重叠矩阵、本征值；不做任何弛豫。
- **默认参数**（`control_scf.in`）：

| 条目 | 默认 | 目的/说明 |
|------|------|-----------|
| `xc` | `pbe` | + `vdw_correction_hirshfeld`（若 `electrode.vdw=true`） |
| `spin` / `fixed_spin_moment` | `none` / `0.0`；开壳层 `collinear` / `M−1` | 闭壳层非极化不变。`molecule.multiplicity > 1` 时 prepare 写 `spin collinear` + `fixed_spin_moment = M−1` + 初猜自旋密度（`molecule.initial_moments` 映射成 per-atom `initial_moment`；未指定时退回 `default_initial_moment 0.1`），非零 `charge` 把总电荷均分成 per-atom `initial_charge` 写进 geometry.in；`molecule.fixed_spin_moment: false` 则不写约束行（release 对照）；SCF 改产 `alpha.aims` + `beta.aims`（不再写 `mos.aims`） |
| `charge` | `0` | |
| `relativistic` | `atomic_zora scalar` | Au 必须（模板对所有元素统一开启；Ag/Cu/Ni/Pt/Pd/Al 同样适用） |
| `occupation_type` | `gaussian 0.01` | 展宽，帮助收敛 |
| `mixer` | `pulay`，`n_max_pulay 10`，`charge_mix_param 0.2` | |
| `sc_accuracy_rho` / `eev` / `etot` | `1E-5` / `1E-3` / `1E-6` | |
| `sc_iter_limit` | `500` | |
| `KS_method` | `serial` | **关键**：默认 ELPA/ELSI 走 ScaLAPACK 打包矩阵（`PM_index`），`omat.aims` 布局 aitranss 读不了，会出现负重叠本征值；serial LAPACK 保持三角重叠布局 |
| `output aitranss` | 开 | 写 `omat.aims` / `mos.aims` / `basis-indices.out` |
| `output mulliken` | 开 | |
| `restart aims.restart` | 开 | 可续算 |

- **产物**：`aims.out`、`omat.aims`、`mos.aims`、`basis-indices.out`（四件缺一不可）。
- **闸门**：`aims.out` 非空且含 `Have a nice day`；其余三件必须存在（tcontrol 的 `$basis/$read_omat/$scfmo` 直接读，提前 fail-closed）；prepare 时若 metadata 在，校验父几何指纹 `verify_parent_geometry`。
- **资源默认**：132 原子锚点 4320 s / 13.3 GiB，墙钟幂律指数 3.0、内存指数 2.0；`ntasks 4`；墙钟下限 1800 s；上限 72 h / 128 GiB。

---

### 步骤 6：`aitranss_te` —— NEGF 透射谱

- **阶段名**：`aitranss_te`（别名 `te` / `aitranss` / `transport`）；依赖：`aims_scf`；远端（`sub_te.sh`）。
- **目的**：用 Aitranss 在 SCF 结果上做自能 + Landauer 透射，输出 `TE.dat`（第 1 列绝对能量 Hartree，第 2 列 `E−E_F` eV，第 3/4 列自旋通道）。
- **输入**：`aims_scf/{aims.out, mos.aims, geometry.in}`（闭壳层）或 `aims_scf/{aims.out, alpha.aims, beta.aims, geometry.in}`（开壳层；FHI-aims 在 `spin collinear` 下不写 `mos.aims`，nsaos 从 `alpha.aims` 解析）+ `electrode_build/electrode_metadata.json`（可 `--metadata` 覆盖）。
- **默认参数**：

| 配置键 | 默认 | 说明 |
|--------|------|------|
| `transport.temperature` | `300.0` K | 后处理积分温度（不写进 tcontrol） |
| `transport.energy_window` | `(-2.0, 0.005, 2.0)` eV（代码默认与 init 模板一致） | `[起点, 步长, 终点]`；单位由 `transport.energy_window_unit` 决定（默认 `eV` = **相对 E_F**，即 E_F ±2 eV、801 点）；prepare 时按 `aims_scf` 的 efermi 换算成 aitranss 的 Ha 绝对值，缺 efermi 直接 `StageBlocked` |
| `transport.energy_window_unit` | `eV` | `eV` = 相对费米能级；`Hartree` = 旧口径（Ha 绝对三元组，逐字节兼容 0.2.0 行为）。缺省且三元组跨度 ≤0.1 Ha、量级 <2 Ha 时按旧 Ha 口径识别并给迁移提示 |
| `transport.ener/estep/eend` | `null` | 单独覆盖窗口分量，**恒按 Ha 绝对能量**解释（优先级最高） |
| `transport.require_te_out` | `true` | 要求 `te.out` 在场（外部/历史 run 才设 false；不进 protocol 指纹） |
| `electrode.nlayers` | `null` | 解析顺序：显式值 → metadata 库建议 → 回退 3 |
| 自能 `σ_near` / `σ_far` | `0.05` / `0.10` Ha | **模型常数，不是收敛参数**：离 tip 最近的 3 个电极原子 0.05、其余 6 个 0.10；tip 与分子为 `empty`(0.0)；生产必须原样使用，否则不再与本仓冻结基准可比 |

- **写出输入**：`tcontrol`（`$coord/geom`、`$basis=basis-indices.out`、`$read_omat=omat.aims`、`$scfmo=mos.aims`（闭壳层）或 `$uhfmo_alpha file=alpha.aims` + `$uhfmo_beta file=beta.aims`（开壳层）、`$nsaos`、六个表面索引、`$nlayers`、`$s1i 0.1/$s2i 0.05/$s3i 0.025`、`$ener/$estep/$eend`、`$output file=TE.dat`、`$testing off`、`$ecp off`、`$valence_electrons`（按几何里出现的元素自动求和：Au 4f 计入价层=31，Ag/Cu=11，Pt/Pd/Ni=10，Al=3；开壳层带电结按 `分子电荷` 扣减/增加）、`$efermi`）；`self.energy.in`（逐原子虚部自能，Fortran D 格式）；`sub_te.sh`（`LD_LIBRARY_PATH=$mkl_lib`，`te_cpus` 默认 16）。
- **闸门**：tcontrol 的 `$coord` 几何与 `self.energy.in` 用的几何必须**同一份**（sha256 相等，否则 StageBlocked）；六个表面索引必须恰好 3+3；电极金属数 ≥4 且为偶数（元素取 `electrode_metadata.json` 的 `electrode_element`，缺键按 Au）；`TE.dat` 可解析；`te.out` 缺失且 `require_te_out=true` → StageError。
- **资源默认**：132 原子锚点 2552 s / 8.0 GiB，墙钟指数 3.0、内存指数 1.5；`ntasks 1`、`cpus 1`（该构建 Aitranss 不并行，`-c 16` 实测反而慢 ~10%）；墙钟下限 1800 s。

---

### 步骤 7：`junction_analyze` —— Landauer 积分（本地）

- **阶段名**：`junction_analyze`（别名 `analyze` / `analysis` / `landauer`）；依赖：`aitranss_te`；本地阶段。
- **目的**：把 `TE.dat` 积分为可报告的输运量，写 `results.json`：
  - `L_n = |∫ (E−E_F)^n · T(E) · (−∂f/∂E) dE|`（梯形积分，E 用 eV 列）；
  - `log10_G_over_G0_like = log10(L0)`（L0 数值即 G/G0，以 `2e²/h` 每自旋为单位）；
  - `S_uV_per_K = 1e6 · (−L1/(L0·T))`；
  - `k_e_pW_per_K = 1.602e-7 · (L2 − L1²/L0) · 2/(h·T)`；
  - `ZT_electronic = 1 / ((L0·L2)/L1² − 1)`（**只是电子项，不是总 ZT**，未测晶格热导）。
- **默认参数**：`temperature = transport.temperature`（默认 300 K）；常数 `k_B=8.6173324e-5 eV/K`、`h=4.135667662e-15 eV·s`；自旋极化时主指标仍按两列取平均（与 vendored Landauer 口径一致）。
- **产物**：`junction_analyze/results.json`（含 `fermi_nearest`、`max_transmission`；开壳层另含 `spin` 块：`L0_alpha/L0_beta`、per-channel `g_over_g0`、`T_alpha_at_ef/T_beta_at_ef`、`spin_polarization_at_ef`、`spin_polarization_L0`、`max_transmission_alpha/beta`）。
- **免责**：点数 ≤7 时 note「7 点样本免责：仅验证管线，非生产科学结论（`_like` 语义）」。
- **可选**：`smj analyze` 生成 `analysis.json`（schema `smj.analysis.v1`），用**独立解析器 + Simpson 求积 + CODATA 常数**复核主结果，并可 `plot=True` 出 `T_E_plot.png`（matplotlib 属可选 extra `plot`）。图**纵轴为对数**：下限取最小正值的一半（`T_LOG_FLOOR=1e-12` 兜底），`T≤0` 的点在轴上截断而非丢弃，否则 T(E) 的 1e-3 量级深谷（反共振/干涉谷）在线性轴上会被压成一条贴底的线。

---

### 步骤 8：`smj report` —— 汇总与物理判定

- **命令**：`smj report [--force] [--no-checks]`
- **目的**：把状态、输入、资源、输运结果、失败原因、产物清单汇总成 `report.md`，并调用 `validate_junction_run` 生成物理判定 `verdict`（`physics_verdict`、`wording`、`allowed_action`、`failed/inconclusive/missing_evidence_criteria`）。
- **默认**：缺 `TE.dat` → `inconclusive`；`--force` 覆盖已有报告；`--no-checks` 跳过判据。
- **报告 7 节**：`0. 概览`、`1. 阶段状态`、`2. 输入摘要`、`3. 资源估算`、`4. 输运结果`、`5. 失败原因`、`6. 产物清单`。
- **判定措辞**：`passed` = 支持 H0、`failed` = 否证 H0、`inconclusive` = 无法判定；空 checks **不等于通过**（F2）。

---

## 3. 默认参数总表（smj.yaml）

### 3.1 run

| 键 | 默认 | 说明 |
|----|------|------|
| `name` | 必填 | 正则 `^[A-Za-z0-9][A-Za-z0-9._-]*$` |
| `dir` | `.` | 相对 yaml 所在目录；最终 run 目录 = `<dir>/<name>` |
| `mode` | `slurm-ssh`（数据类）/ `dry-run`（init 模板） | `dry-run` / `local` / `slurm-ssh` |
| `clean` | `false` | prepare 前清空本 run 的阶段目录 |
| `keep_going` | `false` | 阶段失败后是否继续（`--keep-going` 覆盖） |

### 3.2 molecule

| 键 | 默认 | 说明 |
|----|------|------|
| `xyz` | `null`（init 写 `molecule.xyz`） | 输入结构 |
| `charge` | `0` | |
| `multiplicity` | `1` | |
| `method` / `basis` | `b3lyp` / `6-31g(d)` | 仅 Gaussian 引擎 |
| `family` | `null` | 判据标定族；`null`/未登记 → 全部判据适用（fail-safe） |
| `opt_engine` | `aims` | `aims` / `gaussian` |
| `opt_functional` / `opt_basis` | `pbe` / `light` | |
| `opt_force_threshold` | `1.e-2` | eV/Å |

### 3.3 electrode

| 键 | 默认 | 说明 |
|----|------|------|
| `au_per_electrode` | `30` | 仅 builtin，上限 30；`null` = 60 原子/侧，`30` = 35 原子/侧（按层取整）；`source=library` 时完全忽略（显式写会 warning） |
| `element` | `Au` | 电极金属（fcc(111)）：Au/Ag/Cu/Ni/Pt/Pd/Al 有内置晶格常数；其它元素必须显式 `a_lattice`；tip 同元素 |
| `a_lattice` | `null` | fcc 晶格常数 Å（2.0–6.0）；`null` = 内置值（Au 4.07823、Ag 4.08530、Cu 3.61496、Ni 3.52360、Pt 3.92360、Pd 3.89070、Al 4.04950） |
| `left_anchor` / `right_anchor` | `S` / `S` | 任意元素符号或 1-based 序号（无 S 分子写 N/Se 等） |
| `au_tip_distance` | `2.35` Å | 1.0–5.0 |
| `anchor_au_s_target` | `null` | 1.5–4.0 Å |
| `pin_contact` | `false` | true = 冻结锚定原子 + 两 tip |
| `strip_thiol_h` | `null` | `null` = 仅 S 锚自动去质子；`true`/`false` = 强制/不去 |
| `vdw` | `true` | Hirshfeld |
| `nlayers` | `null` | 1–8；null = 库建议 |
| `source` | `library` | **构建方式**：`library`（FHI-aims electrodes.library）/ `builtin`（本地三角岛） |
| `library_dir` / `library_cluster` | `null` / `fcc111.au.clusters.37_atoms/...x1.aims` | |
| `relax` | `false` | legacy `au_relax` 口径 |
| `relax_force_threshold` | `1.e-2` | legacy |
| `relax_freeze` | `backplane` | legacy，`backplane` / `electrodes` |

### 3.4 transport

| 键 | 默认 | 说明 |
|----|------|------|
| `temperature` | `300.0` K | |
| `energy_window` | `(-2.0, 0.005, 2.0)` eV（代码默认 = init 模板） | `[起点, 步长, 终点]`；默认相对 E_F（E_F ±2 eV） |
| `energy_window_unit` | `eV` | `eV`（相对 E_F）或 `Hartree`（旧 Ha 绝对口径） |
| `ener` / `estep` / `eend` | `null` | 分量覆盖，恒按 Ha 绝对能量 |
| `require_te_out` | `true` | |

### 3.5 cluster

`host/user/key/remote_root/partition/ntasks/gres/account/qos/species_dir/time_limit/mem` 默认全 `null`；`port 22`；`poll_seconds 60.0`；`timeout_hours 6.0`；`env {}`（键小写、区分大小写；值做 `~`/`${VAR}` 展开）。

跨段校验：`slurm-ssh` 必须有 `host/partition`；`local` 且无 `species_dir` 也无 `paths.species_defaults` → warning（本地真跑 AIMS 会失败）。

### 3.6 paths / review

| 键 | 默认 | 说明 |
|----|------|------|
| `paths.templates` / `paths.species_defaults` | `null` / `null` | 缺省用包内 vendored 模板 |
| `review.connection_check` | `true` | 复核闸门：没有有效连接核对记录 → blocked（旧键 `surface_check` / `require_a2` 仍生效） |
| `review.connection_check_deviation` | `0.15` | 表面索引偏差容差（旧键 `a2_deviation` 仍生效） |
| `review.auto_connection_check` | `true` | `electrode_build` collect 后自动审计并登记（`who=auto`）；设 `false` 恢复人工签署 |
| `review.accept_deviation` | `null` | CLI `--accept-deviation` 覆盖门槛（显式人工豁免，报告须标注） |
| （CLI）`review --connection-check/--no-connection-check` | 不设 | 单次覆盖 `review.connection_check`，不写回配置（旧 `--a2/--no-a2` 仍可用） |

---

## 4. 执行模型：状态机与动作

每个阶段 `prepare → submit → wait → collect`：

| 动作 | 行为 | 默认/边界 |
|------|------|----------|
| `prepare` | 把上游产物 + 配置翻译成本地输入文件 | 拒绝在 submitted/waiting/finished/collected 上重跑（除 `--force`）；先 `ensure_run_spec`，协议漂移默认拒绝 |
| `submit` | 清理残留 → 上传 → `sbatch` 记 job id | 本地阶段报错；`dry-run` 报错；脚本含 `UNCONFIGURED` 占位即 StageError |
| `wait` | 轮询作业状态 | `poll_seconds` 默认 60 s；**超时保留 waiting，不判失败**（作业可能仍在跑） |
| `collect` | 回取产物并解析 | 缺件必须 StageError、**不得标 collected** |

状态：`pending → prepared → submitted → waiting → finished → collected`；另有 `failed`、`skipped`；`dry-run` 只到 `prepared`。

运行模式：

| 模式 | 行为 |
|------|------|
| `dry-run` | 只 `prepare`，不建传输层、不提交；整链干跑会停在缺上游产物的阶段（fail-closed） |
| `local` | 本机 `bash` 同步执行；哨兵 `.dsh_done` 记退出码；超时 = `timeout_hours × 3600` |
| `slurm-ssh` | 上传 → `sbatch` → 轮询 → 回取；需 `host/partition` |

退出码：`0` 成功/作业已结束；`1` 失败或配置/输入错误；`2` 未完成（`--wait` 超时、run 有阶段 waiting）、连接核对签署被拒绝、argparse 用法错误。

阶段跳过规则：`smj run` 对已 `collected/skipped` 的阶段默认跳过（`--force` 重跑）；`finished` 只补 `collect`；遇到 `waiting` 一律中断链条。

---

## 5. 溯源与可复现

- 第一次 `prepare` 写 `run_spec.json`：完整解析配置 + 模板默认值 + 分子/电极库 sha256 + `protocol_sha256` + SMJ/Python 版本 + git commit/dirty + 软件路径 + 本地可读二进制 sha256。
- `protocol_sha256` **只覆盖会改变物理结果的子集**：`molecule`（含 `xyz_sha256`）、`electrode`（去掉 `library_dir`，改记 `library_file_sha256`）、`transport`（`resolved_window()`，去掉 `require_te_out`）；路径/资源/报告设置不算。
- 协议漂移默认拒绝，`--force` 重绑并把旧快照存 `run_spec.previous.json`。
- 每阶段记 `inputs_sha256 / outputs_sha256 / parent_stage / parent_fingerprint`；`aims_scf`、`aitranss_te` prepare 时校验上游 geometry 与 metadata 的 `surface_fingerprint`，不一致 StageBlocked。
- `manifest.json` 带 `run_spec` + `provenance`（`reconstructed` 标记）。

---

## 6. 验收判据（建结几何 + 物理）

### 6.1 状态与结构闸门

`validate_junction_geometry` 7 项检查（`smj review`、`smj preflight`、`smj report` 共用）：

| 检查 | 级别 | 规则 |
|------|------|------|
| `junction_geom_inputs` | error | `geometry.in` 必须存在 |
| `junction_geom_metadata` | error | `electrode_metadata.json` 必须存在 |
| `junction_geom_surface_valid` | error | 左右表面索引非空、在 1..N、互不相交、全为电极元素（否则 aitranss 报「左右电极共享原子」） |
| `junction_geom_surface_side` | error | 输运轴取左右 linker 连线主导坐标轴；左索引在中点左侧、右索引在右侧 |
| `junction_geom_contacts` | error | 元素对键长下限 = 0.7×共价半径和，3 Å 预筛；短于下限视为原子重叠 |
| `junction_geom_anchor_aus` | warn/info | thiolate 落在 2.40–2.44 Å → ok（Jackson 2009 实验锚点 2.42±0.02 Å）；硫醚/硫醇按配位接触（≈2.6–2.7 Å vdW 和）；未识别化学 → info |
| `junction_geom_electrode_extent` | info | 仅报告电极沿输运轴延伸量 |
| `junction_connection_check` | error/warn | 连接核对有效 ok；有记录但指纹不符 error；无记录 warn（旧 check 名 `junction_surface_check` / `junction_a2_confirmation` 兼容） |

表面三角选择（det-v1）：每侧 3 个原子，按「与 tip 同元素、排除对侧 tip、按输运轴中点分侧」取外层最远层，排序键 =（最长边、到 tip 投影距离、面积）各按 **0.05 Å 分桶** + 顶点坐标字典序，结果与原子顺序无关、可复现；同桶多解记 `tie_count`。旧实现用绝对 1e-3 Å + 遍历顺序，同一几何不同次运行会给出不同索引（三种优化设置各一套），因此改为确定性规则。

### 6.2 物理判据（10 条）

- physical 判据必须带 `source{doi,locator}`（无源不得作闸门）；engineering 判据免源但须标 `kind`。
- `gate=error` 未过 → `failed`；`gate=warn` 未过 → `inconclusive`；`missing=True`（证据缺失）→ `inconclusive`，**证据缺失 ≠ 否证**。
- 有源 ≠ 有阈值：标 `evidence=derived` 的带是文献散布/派生包络，只用于提示，不足以否证 H0。

| 判据 | 类型 | 闸门 | 默认带 | 适用 | 来源（physical） |
|------|------|------|--------|------|------------------|
| `gaussian_no_imaginary` | physical | error | 无虚频 | Gaussian 引擎 | 项目既定口径（任务卡 §3）；低频阈值仍无源，见 `KNOWN_GAPS` |
| `scf_converged` | engineering | error | — | 全链 | — |
| `te_present` | engineering | error | `TE.dat` 必须存在且可解析 | 全链 | — |
| `grid_independence` | physical | warn | 抽稀后 L0/L1/L2 相对变化 <1% | 全链 | 任务卡 §5.1.3（实测 ≤0.09%） |
| `electrodes_disjoint` | engineering | error | 左右电极无共享原子 | 全链 | — |
| `au_s_bond` | physical | warn | 2.40–2.44 Å | **仅 S 锚去质子硫醇盐；硫醚、非 S 锚只报接触距离** | Jackson 2009 摘要：Au(111)(√3×√3)R30–CH₃S atop，NIXSW+S1s PhD，S–Au=2.42±0.02 Å（表面吸附相，仅方向性；Thygesen 2005 §V.2 的 2.40 Å 是采用 Stokbro 2003 的输入） |
| `t_ef_resonance_attribution` | physical | warn | 近 E_F 峰落在族预期通道侧且 0.15≤\|ΔE\|≤1.0 eV | 全链（族表；未登记族按实测 Seebeck 符号自洽回退） | Thygesen 2005 §V.2：窄峰 −1.5/−1.0 eV、宽峰 −3.0/+2.5 eV；0.15 eV 下限按 PBE 半局域能级误差 0.1–0.3 eV（Gulyaev 2026） |
| `g_in_theory_band` | physical | warn | 0.03–0.40 G0（**派生包络**，evidence=derived） | BDT | Thygesen 2005 §V.2：around 0.1 G0 的 1/3–4 倍（Ref.15 ×1/3 … Ref.48 ×4；原文无此区间） |
| `seebeck_sign` | physical | error | S > 0 | BDT | Reddy 2007 +8.7±2.1 µV/K、Lee 2014 +7.4±0.5 µV/K（均 Au–BDT–Au 实测） |
| `k_e_zt_order` | physical | warn | 未实现（来源待补） | 全链 | —（`KNOWN_GAPS`） |

**共振判据的语义**（2026-11-15 改）：E_F±1.0 eV 内的显著峰不再是「错误信号」，而是**需解释的物理信号**。

- 无显著峰 → `pass`（含「TE 窗口截断导致看不到峰」的情形，判据不要求峰必须存在）；
- 峰在预期通道侧且 0.15≤\|ΔE\|≤1.0 eV → `pass`，detail 记「物理信号已归因」；
- 侧别不符 / 距 E_F <0.15 eV / 通道未知 → `warn`，detail 说明需要解释的方向；
- 0.15 eV 下限是 PBE 能级误差量级的**工程容差**，不是「允许的物理共振距离」；PBE 会把峰系统性推向 E_F（本管线 BDT 最近峰 −0.339 eV），因此不要求复现文献峰位。

判据族（`molecule.family`）与预期通道（`EXPECTED_RESONANCE`，均带 `source{doi,locator}`）：

| 族 | 预期通道 | 依据 |
|----|----------|------|
| `Au-BDT-Au` | HOMO（占据侧） | Thygesen 2005 §V.2 窄峰在 −1.0 eV；Reddy 2007 S=+8.7±2.1 µV/K |
| `Au-(4,4'-bipyridine)-Au` | LUMO（空侧） | Widawsky 2012 Table 1：S=−9.5 µV/K，LUMO 在 E_F 上方 1.47–1.53 eV |
| `Au-(2,6-azulenedithiol)-Au` | HOMO（侧） | Xia 2014 GW；E_F 绝对位置依赖模型 |
| `Au-alkane-dithiol-Au` | 无（宽 gap） | Xu & Tao 2003：β_N≈1.0±0.1，无近 E_F 共振 |
| 其它（2CN、SME…） | 未登记 | 按实测 Seebeck 符号自洽归因（T(E) 与 S 是同一组数据的不同矩），不自称文献先验 |

未登记/`null` 族仍全部适用（fail-safe）；`t_ef_magnitude`（0.03–0.40 G0）已于 2026-11-15 按「无一手区间则取消」规则删除，其派生包络并入 `g_in_theory_band` 并降为 warn。0.2.0 三体系在新口径下：BDT/AZ/2CN 均 `passed`（AZ/2CN 原先因「E_F 附近有峰」被记 inconclusive）；`bdt-anchor24` 因 S=−22.2 µV/K 对 BDT 族判 `seebeck_sign` failed。

---

## 7. 目录布局与产物

```
<run.dir>/<run.name>/
├── state.json                 # 唯一状态源
├── run_spec.json              # 首次 prepare 的不可变配置快照
├── run_spec.previous.json     # --force 重绑时的旧快照
├── manifest.json              # 报告 + verdict
├── report.md                  # 7 节报告
├── review.json                # 复核结果
├── electrodes/                # 远端取回的库团簇缓存
├── logs/
├── molecule_opt/              # geometry.in control.in sub_aims.sh input_summary.json aims.out opt.xyz
├── aims_tip_relax/            # geometry.in control.in sub_aims.sh aims.out geometry.in.next_step
├── electrode_build/           # geometry.in electrode_metadata.json
│                              # + connection_check.json / connection_check_audit.json / connection_check_viz.png（登记后；旧 a2_* 兼容读）
├── aims_scf/                  # control.in sub_aims.sh geometry.in aims.out mos.aims omat.aims basis-indices.out
├── aitranss_te/               # tcontrol sub_te.sh self.energy.in TE.dat te.out
└── junction_analyze/          # results.json
```

旧目录名（`gaussian_opt`、`au60_build`、`auN_build`）有回退读取，报告按新名展示。

---

## 8. 物理口径与单位

- 能量：`transport.energy_window` 默认是 **eV 且相对费米能级**（`E − E_F`，默认 ±2 eV、步长 5 meV）；`aitranss_te.prepare` 从 `aims_scf/aims.out` 取 efermi，把它换成 tcontrol `$ener/$estep/$eend` 与 `TE.dat` 第 1 列所需的 **Ha 绝对能量**（`ener = efermi + 起点_eV/Ha2eV`、`estep = 步长_eV/Ha2eV`，Ha2eV = 27.211386245988）。缺 efermi 直接 `StageBlocked`。旧 Ha 绝对三元组（0.2.x 默认 `[-0.21, 0.0002, -0.17]`）仍可读，可用 `energy_window_unit: Hartree` 显式声明；`transport.ener/estep/eend` 恒为 Ha 绝对值覆盖。换金属/换分子后 E_F 会移动，但相对窗口不跟着漂。
- `TE.dat` 第 2 列是 `E − E_F`（eV）。
- 透射单位：aitranss 以 `G0 = 2e²/h` **每自旋**为单位；非极化单列，极化两列取平均；`L0` 数值即 `G/G0`。
- 后处理常数：`k_B = 8.6173324e-5 eV/K`、`h = 4.135667662e-15 eV·s`。
- `log10_G_over_G0_like` 的 `_like` 后缀表示项目近似命名，7 点样本不得当生产科学结论。
- fcc(111) 内置晶格：`Δz = a/√3`、`Δy = Δz/√2 = a/√6`、`Δx = Δz·√3/(2√2) = a/(2√2)`。登记值 `a`（Å）：Au 4.07823、Ag 4.08530、Cu 3.61496、Ni 3.52360、Pt 3.92360、Pd 3.89070、Al 4.04950；Au 即 Δz=2.35455、Δy=1.66492、Δx=1.44186。其它元素必须显式给 `electrode.a_lattice`。

---

## 9. 资源估算模型（默认）

| 阶段 | 模型 | 墙钟锚点 | 内存锚点 | 并行 | 墙钟下限 |
|------|------|----------|----------|------|----------|
| `molecule_opt` | 线性 | 60 s + 10 s/atom | 4 GiB + 0.05 GiB/atom | 1 | 900 s |
| `aims_tip_relax` | 线性 | 180 s + 30 s/atom | 8 GiB + 0.10 GiB/atom | 4 | 900 s |
| `au_relax`（legacy） | 线性 | 300 s + 40 s/atom | 10 GiB + 0.10 GiB/atom | 4 | 1800 s |
| `aims_scf` | 幂律（锚点 132 原子） | 4320 s × (N/132)^3.0 | 13.3 GiB × (N/132)^2.0 | 4 | 1800 s |
| `aitranss_te` | 幂律（锚点 132 原子） | 2552 s × (N/132)^3.0 | 8.0 GiB × (N/132)^1.5 | 1（cpus 1） | 1800 s |
| `electrode_build` | 本地 | — | — | — | — |

- 余量：墙钟 ×1.25、内存 ×1.20；Slurm `--time` = 预期墙钟 × **1.9**（SCF 收敛后写产物还要 ~25 min，132 原子曾被 `--time` 杀掉）。
- 硬上限：72 h / 128 GiB，超限则封顶并报 `resource_capped`。
- 估算值只是**申请下限**，按 132 原子锚点外推。
- 注：`molecule.opt_engine=aims` 时，`molecule_opt` 与 `preflight` 实际借用 `aims_tip_relax` 的线性行估算（表中 60 s + 10 s/atom 一行仅 Gaussian 口径）。

---

## 10. 已知限制与常见坑

1. **能量窗口径**（2026-11-15 起）：默认 `energy_window` 是 eV 且相对 E_F（`[-2.0, 0.005, 2.0]`，即 ±2 eV、801 点），`aitranss_te` 按 SCF 的 efermi 换成 Ha 绝对窗口；缺 efermi 会 `StageBlocked`。旧 Ha 绝对三元组仍可读（可用 `energy_window_unit: Hartree` 显式声明），但新 run 建议用 eV 相对口径，别照搬 0.2.0 的 `[-0.21, 0.0002, -0.17]`。
2. `transport.require_te_out=true`：缺 `te.out` 直接 StageError；只有外部/历史 run 才设 false。
3. `pin_contact=false` 是默认：分子+2 tip 全自由弛豫；库电极（37 原子）与内置 Au60、旧 20 Au/侧**不可横比**。
4. 7 个能量点只是管线验证样本，`_like` 字段不得当生产结论。
5. `ZT_electronic` 只是电子项；`k_e_zt_order` 判据尚未实现。
6. `species_dir` 必须显式配置（`cluster.species_dir` 或 `paths.species_defaults`），否则远端才失败、白排一次队。
7. 缺部署路径写 `UNCONFIGURED` 占位，提交前扫到即拒（StageError）。
8. `--wait` 超时不是失败（保留 waiting）；别把超时当作业挂了。
9. 连接核对（原 A2）默认自动登记：`auto_connection_check=true` 时管线在 `electrode_build` collect 后自动审计并写 `who=auto` 记录；审计不过则拒绝登记、`smj review` 仍 blocked。`review.connection_check=true` 是闸门本身；`smj review --no-connection-check` 只临时放宽本次复核，`--accept-deviation` 属人工豁免（报告须标注）。想强制人工：`auto_connection_check: false`。
10. Gaussian 与 FHI-aims 产物契约相同（`opt.xyz`），但 aims 无虚频检查，几何靠力阈值 + 人工复核。
11. 库电极文件不进 wheel；run_spec 记 sha256，换库文件 = 换协议。
12. `dry-run` 整链会停在缺上游产物的阶段（fail-closed），这是预期行为。
13. `smj stage au_relax` 只用于历史 run，不在主链。
14. `KS_method serial` 不能改：ELPA/ELSI 打包矩阵会让 aitranss 读到错误重叠布局（负本征值）。
15. `au_per_electrode` **只在 `source=builtin` 时生效**；默认 `source=library` 时被完全忽略（团簇大小由库文件决定，显式写会 warning）。builtin 口径：`null` = 60 原子/侧，`30` 会按层截到 **35 原子/侧（5 层）**（0.2.0rc1 已修正注释与模板说明）。生产建议用库电极，避免这个歧义。
16. `electrode.element`/`a_lattice` 是**物理参数**不是标签：tip 元素、电极团簇晶格、自能/校验/连接核对出图的元素判定全跟着变；库团簇元素与 `element` 不一致会 StageError。已登记晶格常数的只有 Au/Ag/Cu/Ni/Pt/Pd/Al，其它元素必须显式 `a_lattice`。非 Au 电极的 `au_s_bond` 2.40–2.49 Å 带不适用（只报告接触距离）。
17. 锚定原子不限于 S：`left_anchor`/`right_anchor` 可写任意元素符号（结构里恰好 2 个）或 1-based 序号。去质子三态 `strip_thiol_h: null|true|false`：`null`（默认）只对 S 锚自动去 H，N/Se 等锚保留 H；`true` 强制去，但每个锚原子必须恰好 1 个 H（-NH2 带 2 个 H 会 fail-closed，不猜）。非 S 锚的 `junction_geom_anchor_aus` 只报接触距离（info），不套 2.40–2.49 Å 带。

---

## 附：命令速查

```
smj init [DIR] --name myrun --xyz molecule.xyz --mode dry-run
smj preflight [--stage aims_scf] [--no-remote]
smj plan
smj run [--from molecule_opt] [--to junction_analyze] [--dry-run] [--force] [--keep-going]
smj stage molecule_opt --run          # 单阶段全流程
smj stage aims_scf --prepare          # 只生成输入
smj stage aims_scf --submit --wait --collect
smj status [--json]
smj review --approve --who "张三"     # 人工连接核对（想强制人工签：auto_connection_check: false）
smj review --no-connection-check          # 本次复核不要求连接核对（临时，不写回配置）
smj report [--force] [--no-checks]
smj analyze [--force]                 # 独立复核 + 出图
smj audit                             # 本地对账
smj verify                            # 安装自检
```
