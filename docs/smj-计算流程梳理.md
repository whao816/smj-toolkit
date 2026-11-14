# SMJ 计算流程梳理（smj-toolkit 0.2.0rc1）

> 依据 0.2.0rc1 工作树逐文件核对：`src/smj/config.py`、`state.py`、`runner.py`、`resources.py`、`preflight.py`、`review.py`、`report.py`、`provenance.py`，以及 `src/smj/stages/*` 与 `src/smj/_vendor/junction/*`（inputs/electrode/electrode_library/tcontrol/self_energy/landauer/validate/criteria/analysis）。
> 默认值一栏写的是 **不写 smj.yaml 时生效的值**；单位与量纲单独注明。

---

## 0. 一句话架构

SMJ 是一条「分子结构 → 透射谱」的端到端流水线：六阶段固定顺序、阶段之间只通过 `<run>/<stage>/` 下的**文件**耦合，状态全部记在 `<run>/state.json`；每个阶段都是 `prepare → submit → wait → collect` 四动作（本地阶段只有 `prepare → collect`），上游产物缺失时**拒绝猜测**（fail-closed）。

```
molecule.xyz
   │ ① molecule_opt      FHI-aims PBE/light 分子优化（或 Gaussian B3LYP/6-31G(d)）
   ▼ opt.xyz
   │ ② aims_tip_relax    分子 + 2 个 Au tip，全原子 BFGS 弛豫
   ▼ geometry.in.next_step
   │ ③ electrode_build   2 个 tip 扩展为 Au(111) 金字塔团簇（本地阶段）
   ▼ geometry.in + electrode_metadata.json
   │ ④ review（人工）    结构闸门 + A2 红线签署
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
| 2 | 分子几何优化 | `smj stage molecule_opt` | FHI-aims（默认）或 G16 | 得到单分子平衡几何 | PBE/light，charge 0，multiplicity 1，力阈值 `1.e-2` Ha/Bohr，vdW=Hirshfeld | `opt.xyz` |
| 3 | 加 Au tip 并弛豫 | `smj stage aims_tip_relax` | FHI-aims | 去质子硫醇盐 + 2 Au tip，得到 S–Au 接触几何 | 放置距离 2.35 Å，`strip_thiol_h=true`，`pin_contact=false`，全原子自由弛豫 | `geometry.in.next_step` |
| 4 | 构建电极团簇 | `smj stage electrode_build` | 纯 Python（本地） | tip 扩展为 Au(111) 金字塔，冻结结几何 | 库电极 `fcc111.au.clusters.37_atoms`（37 原子，nlayers 3） | `geometry.in` + `electrode_metadata.json` |
| 4.5 | 结构复核 + A2 签署 | `smj review [--approve --who] [--a2\|--no-a2]` | — | 人工确认几何，A2 红线 | `review.require_a2=true`，`a2_deviation=0.15`；CLI 可临时覆盖 | `review.json`、`a2_confirmation.json` |
| 5 | 静态 SCF | `smj stage aims_scf` | FHI-aims | 为 aitranss 准备冻结几何的单点波函数 | PBE，`KS_method serial`，sc_accuracy_rho 1e-5 / eev 1e-3 / etot 1e-6，sc_iter_limit 500 | `aims.out`、`omat.aims`、`mos.aims`、`basis-indices.out` |
| 6 | 透射谱 | `smj stage aitranss_te` | Aitranss | NEGF + 自能，算 T(E) | 能量窗 `[-0.21, 0.0002, -0.17]` Ha（默认，围绕 Au Fermi），σ_near 0.05 / σ_far 0.10 Ha | `TE.dat`（+`te.out`） |
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
- **写出的关键默认**：`run.dir "."`、`run.clean false`、`run.keep_going false`；`molecule.opt_engine aims`、`opt_functional pbe`、`opt_basis light`、`opt_force_threshold 1.e-2`；`electrode.au_per_electrode 30`、`source library`、`au_tip_distance 2.35`、`pin_contact false`、`strip_thiol_h true`、`vdw true`、`nlayers null`；`transport.temperature 300.0`、`energy_window [-0.21, 0.0002, -0.17]`（Ha，围绕 Au Fermi 的 ±0.02 Ha 口径）；`review.require_a2 true`、`a2_deviation 0.15`。

> 💡 能量窗默认值：0.2.0rc1 起 `transport.energy_window` 默认 `(-0.21, 0.0002, -0.17)` Ha（201 点，步长约 5.4 meV，围绕 Au Fermi ≈ −0.19 Ha），与 `smj init` 模板一致。生产仍须在 `aims_scf` 后按实际 E_F 调整窗口；不要用 5 Ha 全窗（0.02 Ha 步长 = 0.54 eV，既浪费又会让 E_F 附近采样过粗）。

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
| 集群配置 | error | slurm-ssh 必须有 `host / remote_root / partition` |
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
| `molecule.charge` | `0` | |
| `molecule.multiplicity` | `1` | `>1` 时转成 FHI-aims `spin = multiplicity-1` |
| `molecule.opt_functional` / `opt_basis` | `pbe` / `light` | `opt_basis` 只记录；实际基组由 species 目录决定 |
| `molecule.opt_force_threshold` | `1.e-2` | 写进 `relax_geometry bfgs`，单位 Ha/Bohr |
| `electrode.vdw` | `true` | → `vdw_correction_hirshfeld`；只改能量/力（几何），不改固定几何的 T(E) |

- **写出的输入**（aims 分支）：`geometry.in`、`control.in`（`control_relax.in` 模板：`xc pbe`、`spin`、`charge`、`relativistic atomic_zora scalar`、`relax_geometry bfgs 1.e-2`、`aggregated_energy_tolerance 1.0e-1`、`output cube eigenstate HOMO-1/HOMO/LUMO/LUMO+1`、`output dipole`、species 块）、`sub_aims.sh`、`input_summary.json`。
- **产物**：`opt.xyz`（下游统一契约）、`aims.out`、`geometry.in.next_step`；Gaussian 分支为 `mol.gjf` → `mol.log` → `opt.xyz` + `gaussian_summary.json`。
- **验收（collect 闸门）**：
  - aims：`aims.out` 必须含 `Have a nice day`，`geometry.in.next_step` 必须存在，否则 StageError（不把残缺几何交给下游）；
  - gaussian：必须 `Normal termination`，否则先删 `opt.xyz` 再报错；虚频只记 note（需人工确认）。
  - 结果写 `results.json` 的 `molecule_opt` 段（engine、normal_termination、atom_count）与 `software` 段。
- **资源默认**：借用 `aims_tip_relax` 线性模型（60 s + 10 s/atom 是 resources 表里的分子项，但实际操作按 tip_relax 标定）；墙钟下限 900 s，`--time = 预期 × 1.9`，`ntasks 1`。

> 备注：FHI-aims 分支**没有虚频检查**（未提供一键 freq），几何质量靠 BFGS 力阈值 + `smj review` 人工确认。Gaussian 分支的虚频只提示、不阻断。

---

### 步骤 3：`aims_tip_relax` —— 加 Au tip 并全原子弛豫

- **阶段名**：`aims_tip_relax`（别名 `tip` / `relax`）；依赖：`molecule_opt`；远端（`sub_aims.sh`）。
- **目的**：把分子两端锚定原子各接一个 Au tip，**全原子第一性弛豫**，让 S–Au 接触由力场（DFT）自己找位置，而不是人为固定距离。
- **输入**：上游 `molecule_opt/opt.xyz`（可 `--geometry` 覆盖）。
- **默认参数**：

| 配置键 | 默认 | 说明 |
|--------|------|------|
| `electrode.left_anchor` / `right_anchor` | `S` / `S` | 元素符号必须在结构中出现**恰好 2 次**（第一个=左，第二个=右），也可写 1-based 序号 |
| `electrode.au_tip_distance` | `2.35` Å | tip 初始放置距离；校验 1.0–5.0 |
| `electrode.anchor_au_s_target` | `null` | 非空时直接把 tip 放到该 Au–S 距离（1.5–4.0 Å）；生产一般留空，让弛豫决定 |
| `electrode.pin_contact` | `false` | `false` = 分子+2 tip **全部自由**；`true` = 冻结两锚定原子 + 两 tip（距离扫描/固定接触口径） |
| `electrode.strip_thiol_h` | `true` | 移除 S–H < 1.70 Å 的氢（去质子硫醇盐），必须恰好 0 或 2 个；原几何存 `opt_neutral.xyz` |
| `electrode.vdw` | `true` | Hirshfeld 色散 |
| `molecule.charge` / `multiplicity` | `0` / `1` | |
| `opt_force_threshold` | `1.e-2` | 同上（Ha/Bohr） |

- **放置规则**：默认沿「左锚–右锚」连线向外（`add_tip_au` 默认方向 = 本锚减对侧锚），距离 = `au_tip_distance` 或 `anchor_au_s_target`；落盘后**立即回测** Au–S 距离，偏差 >1e-6 Å 直接 RuntimeError（拒绝悄悄错位）。
- **产物**：`geometry.in.next_step`（BFGS 最终几何）、`aims.out`、`geometry.in`、`input_summary.json`。
- **验收**：优先 `geometry.in.next_step`，缺失才回退 `geometry.in` 并记 note；collect 计算**弛豫后真实 Au–S 距离** `relaxed_au_s_A`（tip 定义为文件末尾 2 个原子）。
- **资源默认**：180 s + 30 s/atom；内存 8 GiB + 0.10 GiB/atom；`ntasks 4`；墙钟下限 900 s。
- **备注**：2.35 Å 只是初始放置值，最终 Au–S 由弛豫得到；`preflight/review` 的 `au_s_bond` 判据只对**去质子硫醇盐**用 2.40–2.49 Å 带，硫醚（-SCH3）是配位接触（2.6–2.7 Å），不适用该带。

---

### 步骤 4：`electrode_build` —— 构建 Au(111) 电极团簇（本地）

- **阶段名**：`electrode_build`（旧名 `au60_build`；更早 `au20_build`/`auN_build` 归一化）；依赖：`aims_tip_relax`；**本地阶段**（`needs_submit=False`，prepare 即产出）。
- **目的**：把 2 个 tip 沿 tip→linker 方向扩展成完整 Au(111) 金字塔团簇；**只扩展、不改 S–Au 距离/角度**，此后结几何冻结。
- **默认参数**：

| 配置键 | 默认 | 说明 |
|--------|------|------|
| `electrode.source` | `library` | `library` = FHI-aims 预定义团簇；`builtin` = 内置 `Electrode111` |
| `electrode.library_cluster` | `fcc111.au.clusters.37_atoms/fcc111.cluster.natoms37.x1.aims` | 37 原子，头部 `-nlayers=3`（手册 §5.3.1 同族示例） |
| `electrode.library_dir` | `null` | 也可用环境变量 `SMJ_ELECTRODE_LIBRARY`；否则依次搜 `/data/apps/src/fhi-aims*/src/aitranss/electrodes.library`、`/opt/fhi-aims*`、`~/fhi-aims*`、`~/software/fhi-aims*`；slurm-ssh 下还会 `ssh cat` 远端取回并缓存到 `<run>/electrodes/` |
| `electrode.nlayers` | `null` | `null` = 取库文件头部建议值；与库建议不一致会 note 警告 |
| `electrode.au_per_electrode` | `30` | **仅 builtin 生效**；默认 `source=library` 时**完全忽略**（显式写会 warning）。上限 30（超过报错）。`null` → 全量 60 Au/侧（body 56 + 4 extras）；数值 N → 逐层累加到 ≥N 的最完整层数（N=30 → 35 原子/5 层），不残缺表面层 |
| 内置团簇 | 6 层 + 4 个额外原子 | 三角层 1+3+6+10+15+21 = 56，+4 = 60 Au/侧；`a=4.07823 Å`，`Δz=2.35455`、`Δy=1.66492`、`Δx=1.44186` Å |

- **产物**：`geometry.in`（结几何）、`electrode_metadata.json`（左右 tip 索引、linker 索引、六表面索引 `surc/lsurx/lsury + rsurc/rsurx/rsury`、左右背层索引、左右电极索引集、`final_atom_count`、`surface_rule="det-v1"`、`surface_fingerprint`、库团簇 provenance）。
- **闸门**：结构里必须**恰好 2 个** tip 元素（默认 Au）；linker = 距 tip 最近的**非 Au、非 H** 原子；collect 时 metadata 的 `final_atom_count` 必须等于实际原子数；缺 metadata 指纹会 note「审阅确认无法绑定几何」。
- **资源**：本地，不占集群预算。
- **备注**：库文件**不进 wheel**，run_spec 记录其 sha256 + 相对路径；库电极（37 原子）与内置 Au60、与旧的 20 Au/侧口径**不可直接横比**。

---

### 步骤 4.5：`smj review` —— 结构复核 + A2 红线（人工签署）

- **命令**：`smj review [--approve --who NAME] [--force] [--stage DIR] [--no-figure] [--accept-deviation FRAC] [--a2|--no-a2]`
- **目的**：三段式复核，确保错误结构不进入 SCF：① 结构几何闸门（`validate_junction_geometry` 7 项）；② A2 人工确认（红线，默认开启）；③ 溯源一致性。
- **默认**：`review.require_a2=true` → 没有 `a2_confirmation.json` 时直接 **blocked**（退出码 1）；`a2_deviation=0.15`（表面索引偏差容差）；`--approve` 才写 `electrode_build/a2_confirmation.json`（绑定 `geometry.in` 的 sha256 与签署人）。CLI 可临时覆盖：`--no-a2` 关闭 A2 要求（等价 `require_a2=false`，缺记录/指纹失效降为 info/warning），`--a2` 强制开启；只作用于本次调用，不写回配置、不进协议哈希。
- **判定**：`passed | warning | blocked`，`ok = status != "blocked"`；`review.json` 落 run 根目录。
- **签署流程**（`--approve`）：先跑 `a2_viz.surface_audit`（表面索引/输运轴/接触审计，写 `electrode_build/a2_surface_audit.json`）→ 不通过直接拒绝；默认渲染 `electrode_build/a2_surface_viz.png`（**无图拒绝签署**，`--no-figure` 表示已在别处看过图）；再写 `electrode_build/a2_confirmation.json`（who、when、geometry sha256、表面索引）。已有仍有效的确认默认幂等不覆盖，`--force` 重签；签署被拒返回退出码 2。
- **A2 偏差算法**：`deviation = 1 - Σ|交集| / Σ max(|集合|)`（左右两侧加权，0..1）。≤0.15 → info；≤`--accept-deviation`（显式豁免）→ warning 且 `accepted=True`（报告须标注为人工豁免）；否则 error。
- **溯源检查**：`run_spec.json` 缺失 → warning；`reconstructed=true` → warning；`protocol_sha256` 与当前配置不符 → **error provenance_drift**（run 开始后物理相关配置变了，结果不能按原参数解读）。

---

### 步骤 5：`aims_scf` —— 静态 SCF

- **阶段名**：`aims_scf`（别名 `scf` / `aims`）；依赖：`electrode_build`；远端（`sub_aims.sh`）。
- **目的**：在**冻结的结几何**上做一次静态 SCF，产出 aitranss 需要的波函数、重叠矩阵、本征值；不做任何弛豫。
- **默认参数**（`control_scf.in`）：

| 条目 | 默认 | 目的/说明 |
|------|------|-----------|
| `xc` | `pbe` | + `vdw_correction_hirshfeld`（若 `electrode.vdw=true`） |
| `spin` / `fixed_spin_moment` | `none` / `0.0` | 非极化 |
| `charge` | `0` | |
| `relativistic` | `atomic_zora scalar` | Au 必须 |
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
- **输入**：`aims_scf/{aims.out, mos.aims, geometry.in}`（可 `--scf-dir` 覆盖）+ `electrode_build/electrode_metadata.json`（可 `--metadata` 覆盖）。
- **默认参数**：

| 配置键 | 默认 | 说明 |
|--------|------|------|
| `transport.temperature` | `300.0` K | 后处理积分温度（不写进 tcontrol） |
| `transport.energy_window` | `(-0.21, 0.0002, -0.17)` Ha（代码默认与 init 模板一致） | `[ener, estep, eend]`，**绝对 Hartree**；要求 `estep>0`、`eend>ener`；默认围绕 Au Fermi ±0.02 Ha，生产按实际 E_F 调整 |
| `transport.ener/estep/eend` | `null` | 单独覆盖窗口分量 |
| `transport.require_te_out` | `true` | 要求 `te.out` 在场（外部/历史 run 才设 false；不进 protocol 指纹） |
| `electrode.nlayers` | `null` | 解析顺序：显式值 → metadata 库建议 → 回退 3 |
| 自能 `σ_near` / `σ_far` | `0.05` / `0.10` Ha | **模型常数，不是收敛参数**：离 tip 最近的 3 个电极原子 0.05、其余 6 个 0.10；tip 与分子为 `empty`(0.0)；生产必须原样使用，否则不再与本仓冻结基准可比 |

- **写出输入**：`tcontrol`（`$coord/geom`、`$basis=basis-indices.out`、`$read_omat=omat.aims`、`$scfmo=mos.aims`、`$nsaos`、六个表面索引、`$nlayers`、`$s1i 0.1/$s2i 0.05/$s3i 0.025`、`$ener/$estep/$eend`、`$output file=TE.dat`、`$testing off`、`$ecp off`、`$valence_electrons`（Au 4f 计入价层=31）、`$efermi`）；`self.energy.in`（逐原子虚部自能，Fortran D 格式）；`sub_te.sh`（`LD_LIBRARY_PATH=$mkl_lib`，`te_cpus` 默认 16）。
- **闸门**：tcontrol 的 `$coord` 几何与 `self.energy.in` 用的几何必须**同一份**（sha256 相等，否则 StageBlocked）；六个表面索引必须恰好 3+3；Au 数 ≥4 且为偶数；`TE.dat` 可解析；`te.out` 缺失且 `require_te_out=true` → StageError。
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
- **默认参数**：`temperature = transport.temperature`（默认 300 K）；常数 `k_B=8.6173324e-5 eV/K`、`h=4.135667662e-15 eV·s`；自旋极化时两列取平均。
- **产物**：`junction_analyze/results.json`（含 `fermi_nearest`、`max_transmission`）。
- **免责**：点数 ≤7 时 note「7 点样本免责：仅验证管线，非生产科学结论（`_like` 语义）」。
- **可选**：`smj analyze` 生成 `analysis.json`（schema `smj.analysis.v1`），用**独立解析器 + Simpson 求积 + CODATA 常数**复核主结果，并可 `plot=True` 出 `T_E_plot.png`（matplotlib 属可选 extra `plot`）。

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
| `opt_force_threshold` | `1.e-2` | Ha/Bohr |

### 3.3 electrode

| 键 | 默认 | 说明 |
|----|------|------|
| `au_per_electrode` | `30` | 仅 builtin，上限 30；`null` = 60 Au/侧，`30` = 35 原子/侧（按层取整）；`source=library` 时完全忽略（显式写会 warning） |
| `left_anchor` / `right_anchor` | `S` / `S` | 元素符号或 1-based 序号 |
| `au_tip_distance` | `2.35` Å | 1.0–5.0 |
| `anchor_au_s_target` | `null` | 1.5–4.0 Å |
| `pin_contact` | `false` | true = 冻结锚定原子 + 两 tip |
| `strip_thiol_h` | `true` | S–H < 1.70 Å 去质子 |
| `vdw` | `true` | Hirshfeld |
| `nlayers` | `null` | 1–8；null = 库建议 |
| `source` | `library` | `library` / `builtin` |
| `library_dir` / `library_cluster` | `null` / `fcc111.au.clusters.37_atoms/...x1.aims` | |
| `relax` | `false` | legacy `au_relax` 口径 |
| `relax_force_threshold` | `1.e-2` | legacy |
| `relax_freeze` | `backplane` | legacy，`backplane` / `electrodes` |

### 3.4 transport

| 键 | 默认 | 说明 |
|----|------|------|
| `temperature` | `300.0` K | |
| `energy_window` | `(-0.21, 0.0002, -0.17)` Ha（代码默认 = init 模板） | `[ener, estep, eend]`，绝对 Ha |
| `ener` / `estep` / `eend` | `null` | 分量覆盖 |
| `require_te_out` | `true` | |

### 3.5 cluster

`host/user/key/remote_root/partition/ntasks/gres/account/qos/species_dir/time_limit/mem` 默认全 `null`；`port 22`；`poll_seconds 60.0`；`timeout_hours 6.0`；`env {}`（键小写、区分大小写；值做 `~`/`${VAR}` 展开）。

跨段校验：`slurm-ssh` 必须有 `host/remote_root/partition`；`local` 且无 `species_dir` 也无 `paths.species_defaults` → warning（本地真跑 AIMS 会失败）。

### 3.6 paths / review

| 键 | 默认 | 说明 |
|----|------|------|
| `paths.templates` / `paths.species_defaults` | `null` / `null` | 缺省用包内 vendored 模板 |
| `review.require_a2` | `true` | F4 红线 |
| `review.a2_deviation` | `0.15` | A2 表面索引偏差容差 |
| `review.accept_deviation` | `null` | CLI `--accept-deviation` 覆盖门槛 |
| （CLI）`review --a2/--no-a2` | 不设 | 单次覆盖 `require_a2`，不写回配置 |

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
| `slurm-ssh` | 上传 → `sbatch` → 轮询 → 回取；需 `host/remote_root/partition` |

退出码：`0` 成功/作业已结束；`1` 失败或配置/输入错误；`2` 未完成（`--wait` 超时、run 有阶段 waiting）、A2 签署被拒绝、argparse 用法错误。

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
| `junction_geom_surface_valid` | error | 左右表面索引非空、在 1..N、互不相交、全为 Au（否则 aitranss 报「左右电极共享原子」） |
| `junction_geom_surface_side` | error | 输运轴取左右 linker 连线主导坐标轴；左索引在中点左侧、右索引在右侧 |
| `junction_geom_contacts` | error | 元素对键长下限 = 0.7×共价半径和，3 Å 预筛；短于下限视为原子重叠 |
| `junction_geom_anchor_aus` | warn/info | thiolate 落在 2.40–2.49 Å → ok；硫醚/硫醇按配位接触（≈2.6–2.7 Å vdW 和）；未识别化学 → info |
| `junction_geom_electrode_extent` | info | 仅报告电极沿输运轴延伸量 |
| `junction_a2_confirmation` | error/warn | A2 有效 ok；有记录但指纹不符 error；无记录 warn |

表面三角选择（det-v1）：每侧 3 个原子，按「与 tip 同元素、排除对侧 tip、按输运轴中点分侧」取外层最远层，排序键 =（最长边、到 tip 投影距离、面积）各按 **0.05 Å 分桶** + 顶点坐标字典序，结果与原子顺序无关、可复现；同桶多解记 `tie_count`。旧实现用绝对 1e-3 Å + 遍历顺序，同一几何不同次运行会给出不同索引（三种优化设置各一套），因此改为确定性规则。

### 6.2 物理判据（11 条）

- physical 判据必须带 `source{doi,locator}`（无源不得作闸门）；engineering 判据免源但须标 `kind`。
- `gate=error` 未过 → `failed`；`gate=warn` 未过 → `inconclusive`；`missing=True`（证据缺失）→ `inconclusive`，**证据缺失 ≠ 否证**。

| 判据 | 类型 | 闸门 | 默认带 | 适用 |
|------|------|------|--------|------|
| `gaussian_no_imaginary` | physical | error | 无虚频 | Gaussian 引擎 |
| `scf_converged` | engineering | error | — | 全链 |
| `te_present` | engineering | error | `TE.dat` 必须存在 | 全链 |
| `grid_independence` | physical | warn | 抽稀后 L0/L1/L2 相对变化 <1% | 全链 |
| `electrodes_disjoint` | engineering | error | 左右电极无共享原子 | 全链 |
| `au_s_bond` | physical | warn | 2.40–2.49 Å | **仅去质子硫醇盐** |
| `t_ef_magnitude` | physical | error | 0.03–0.40 G0 | BDT |
| `t_ef_no_peak_near_ef` | physical | warn | E_F ±1 eV 无显著峰 | 全链 |
| `g_in_theory_band` | physical | error | 0.03–0.40 G0 | BDT |
| `seebeck_sign` | physical | error | S > 0 | BDT |
| `k_e_zt_order` | physical | warn | 未实现（来源待补） | 全链 |

判据族（`molecule.family`）：`BDT`、`2CN`、`C4`、`AZ`；未登记/`null` → 全部判据适用（fail-safe）。

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
│                              # + a2_confirmation.json / a2_surface_audit.json / a2_surface_viz.png（签署后）
├── aims_scf/                  # control.in sub_aims.sh geometry.in aims.out mos.aims omat.aims basis-indices.out
├── aitranss_te/               # tcontrol sub_te.sh self.energy.in TE.dat te.out
└── junction_analyze/          # results.json
```

旧目录名（`gaussian_opt`、`au60_build`、`auN_build`）有回退读取，报告按新名展示。

---

## 8. 物理口径与单位

- 能量：`transport.energy_window`、tcontrol `$ener/$estep/$eend`、`TE.dat` 第 1 列均为 **Hartree 绝对能量**；Au 体系 Fermi ≈ **−0.19 Ha**，推荐窗口 ±0.02 Ha（init 模板即 `[-0.21, 0.0002, -0.17]`：起点 −0.21 Ha、步长 0.0002 Ha ≈ 5.4 meV、终点 −0.17 Ha）。
- `TE.dat` 第 2 列是 `E − E_F`（eV）。
- 透射单位：aitranss 以 `G0 = 2e²/h` **每自旋**为单位；非极化单列，极化两列取平均；`L0` 数值即 `G/G0`。
- 后处理常数：`k_B = 8.6173324e-5 eV/K`、`h = 4.135667662e-15 eV·s`。
- `log10_G_over_G0_like` 的 `_like` 后缀表示项目近似命名，7 点样本不得当生产科学结论。
- Au(111) 内置晶格：`a = 4.07823 Å`，`Δz = a/√3 = 2.35455`、`Δy = Δz/√2 = 1.66492`、`Δx = Δz·√3/(2√2) = 1.44186`。

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

1. **能量窗默认值已统一为围绕 Fermi**（0.2.0rc1）：代码默认与 init 模板都是 `(-0.21, 0.0002, -0.17)` Ha；生产仍应在 SCF 后按实际 E_F 调整，别用 5 Ha 全窗。
2. `transport.require_te_out=true`：缺 `te.out` 直接 StageError；只有外部/历史 run 才设 false。
3. `pin_contact=false` 是默认：分子+2 tip 全自由弛豫；库电极（37 原子）与内置 Au60、旧 20 Au/侧**不可横比**。
4. 7 个能量点只是管线验证样本，`_like` 字段不得当生产结论。
5. `ZT_electronic` 只是电子项；`k_e_zt_order` 判据尚未实现。
6. `species_dir` 必须显式配置（`cluster.species_dir` 或 `paths.species_defaults`），否则远端才失败、白排一次队。
7. 缺部署路径写 `UNCONFIGURED` 占位，提交前扫到即拒（StageError）。
8. `--wait` 超时不是失败（保留 waiting）；别把超时当作业挂了。
9. A2 是红线：`review.require_a2=true` 时没有签署 → blocked；`smj review --no-a2` 只临时放宽本次复核，`--accept-deviation` 属人工豁免，报告须标注。
10. Gaussian 与 FHI-aims 产物契约相同（`opt.xyz`），但 aims 无虚频检查，几何靠力阈值 + 人工复核。
11. 库电极文件不进 wheel；run_spec 记 sha256，换库文件 = 换协议。
12. `dry-run` 整链会停在缺上游产物的阶段（fail-closed），这是预期行为。
13. `smj stage au_relax` 只用于历史 run，不在主链。
14. `KS_method serial` 不能改：ELPA/ELSI 打包矩阵会让 aitranss 读到错误重叠布局（负本征值）。
15. `au_per_electrode` **只在 `source=builtin` 时生效**；默认 `source=library` 时被完全忽略（团簇大小由库文件决定，显式写会 warning）。builtin 口径：`null` = 60 Au/侧，`30` 会按层截到 **35 Au/侧（5 层）**（0.2.0rc1 已修正注释与模板说明）。生产建议用库电极，避免这个歧义。

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
smj review --approve --who "张三"     # A2 签署
smj review --no-a2                    # 本次复核不要求 A2（临时，不写回配置）
smj report [--force] [--no-checks]
smj analyze [--force]                 # 独立复核 + 出图
smj audit                             # 本地对账
smj verify                            # 安装自检
```
