# Changelog

本项目遵循语义化版本（0.x 阶段：次版本可含不兼容变更；rc 视为发布候选）。

## 未发布

### 引擎插件层（2026-11-18）

- 新增 `smj.engines` 注册表：`Engine` 协议声明脚本名、产物、资源/环境锚点、预检探测项、
  远端回取清单；`aims` 与 `gaussian` 两个插件分别负责写输入与解析产物，
  `molecule_opt` 只保留电子态校验、`check_geometry` 与 dispatch。
- 收拢 6 处引擎分派：`molecule.opt_engine` 白名单由 `engine_names()` 派生；
  `plan` / `preflight` / `runner` 分别从插件取脚本名、资源锚点与预检项、远端回取名。
  未知引擎在 config 解析阶段即报错，不再等到 prepare。
- 行为不变：默认仍是 FHI-aims PBE/light；旧 `gaussian_opt` 目录仍回退 `sub_g16.sh`。
  新增 `tests/test_engine_registry.py`（8 例）与 Gaussian prepare 用例；全量 276 passed。

## 0.2.1rc1 — 2026-11-18

0.2.0 之后的累积改进（rc 版，未重跑真实集群全链；三体系验证沿用 2026-11-15 的
记录，见 [docs/validation-0.2.0.md](docs/validation-0.2.0.md) 与各 demo report）。

### tip 放置方向 `electrode.tip_direction`（2026-11-18）

- 新增 `electrode.tip_direction: axis | lone_pair`，默认 `axis`——旧行为逐字节不变。
  `lone_pair` 沿锚原子的孤对方向外推 tip（S 为例：C–S–C 角平分线的反向），用于
  S···S 连线与 S–C(Me) 近共线的间位取代硫醚：沿 `axis` 外推会把金属放到甲基一侧，
  而不是锚原子的孤对一侧。
- `aims_tip_relax` 在写输入前校验：锚原子必须恰有 2 个 1.95 Å 内近邻，否则
  `StageError`（fail-closed，不猜方向）；两键共线导致方向退化也直接拒绝。
  `input_summary.json` 记录 `left_direction` / `right_direction`（`axis` 时为 `auto`），
  prepare 的 note 同时给出本次 `direction=` 口径。
- `provenance` 的 protocol 只在非默认值（`lone_pair`）时写入 `tip_direction`，
  历史 run 不会因新增默认键集体判 `provenance_drift`。
- 新增 `tests/test_tip_direction.py`（7 例：默认值/大小写归一化/非法值、`axis` 不动几何、
  `lone_pair` 方向与单位长、近邻数不足与未知 mode 均 fail-closed）；全量 267 passed。

### T(E) 绘图纵轴改为对数（2026-11-16）

- `smj.analysis.plot_transmission`：纵轴 `log`。新增 `T_LOG_FLOOR = 1e-12` 与 `_log_floor()`，
  下限取「最小正值 × 0.5」再夹到 `T_LOG_FLOOR`，`T≤0` 的点在下限处截断而不是丢点/炸轴；
  `ylim = (floor, max(max(T) × 1.15, floor × 10))`。返回 dict 增加 `log_scale`/`floor`，
  `analysis.json` 的 `plot` 键随之带出。纵轴标签统一为 `T(E) (G0)`，对数信息写在标题里
  （`transmission TE.dat (log scale)`）——标签保持 ASCII：DejaVu Sans 无 CJK 字形，
  中文标签会出方框（`Glyph 65288 missing`）。
- 同口径其余出图：`scripts/compare_controls.py` 两个 T(E) 面板（网格收敛差值面板保持线性）、
  `demo/azulene-26-dithiol/plot_results.py` 上方面板（−dT/dE 面板保持线性）、
  `runs/azo_compare.py` 的 E/Z 叠加图（图例由 `azo-e (E)` / `azo-z (Z)` 改为 `E (trans)` /
  `Z (cis)`，去掉与 run 名重复的异构体标注）。
- 顺带修 `demo/azulene-26-dithiol/plot_results.py`：`analysis.json` v1 的 `max_transmission`
  已是 dict（`{E_Hartree,E_minus_Ef_eV,T}`），旧代码 `float(dict)` 直接 `TypeError`，
  现按 `T` 取值。该 script 自 schema 变更后一直跑不起来。
- 新增 `tests/test_analysis_plot.py`（4 例：`_log_floor` 退化输入、log 轴元数据、全零 T(E)、
  `analyze_te(plot=True)` 落盘），全量 260 passed。

### 判据口径修订：共振归因、删 `t_ef_magnitude`、收窄 Au–S 带（2026-11-15）

- `t_ef_no_peak_near_ef` → `t_ef_resonance_attribution`：E_F±1.0 eV 内的显著峰不再直接判
  `inconclusive`，改为按族归因——无峰 `pass`；峰在预期通道侧（HOMO 占据侧 / LUMO 空侧）
  且 `0.15≤|ΔE|≤1.0 eV` → `pass`（detail 记「物理信号已归因」）；侧别不符 / 距 E_F 过近 /
  通道未知 → `warn`（需解释）。新增 `EXPECTED_RESONANCE` 族表（BDT=HOMO、BPY=LUMO、
  AZ=HOMO 侧、烷烃=无近 E_F 共振，均带 `source{doi,locator}`）；未登记族（2CN、SME…）
  按实测 Seebeck 符号自洽回退，不自称文献先验；0.15 eV 是 PBE 能级误差量级的工程容差。
- 删除 `t_ef_magnitude`：0.03–0.40 G0 并非 Thygesen 2005 原文区间（原文只给 around 0.1 G0
  与「比 Ref.48 小 4 倍、比 Ref.15 大 3 倍」），按「无一手源则取消」规则删除；派生包络并入
  `g_in_theory_band`，其 gate 由 `error` 降为 `warn` 并标 `evidence=derived`。
- `au_s_bond` 带 2.40–2.49 → **2.40–2.44 Å**：上端 2.49 无一手源；实验锚点为 Jackson 2009
  摘要 S–Au=2.42±0.02 Å（Au(111)(√3×√3)R30–CH₃S atop，NIXSW+S1s PhD）；Thygesen 2005 §V.2
  的 2.40 Å 是采用 Stokbro 2003 的模型输入，note 已写明；`validate.py` 的 `AU_S_BAND` 与
  4 处 detail 文案同步。
- `te_present`：`TE.dat` 存在但不可解析时也记 `missing`（原先只在文件缺失时记）。
- 证据与派生过程见新增的 [`docs/criteria-source-research.md`](docs/criteria-source-research.md)
  （Thygesen 2005 §V.2 与 Jackson 2009 原文核验、各族通道一手来源）；新增
  `tests/test_resonance_attribution.py`（9 个用例），全量 256 passed。
- 0.2.0 三体系在新口径下：BDT / AZ / 2CN 均 `passed`；`bdt-anchor24` 因 S=−22.2 µV/K 对
  BDT 族判 `seebeck_sign` failed。

### 分子弛豫的 SCF 精度与自旋固定

- `control_relax.in` 之前只写几何相关的键，FHI-aims 会打印
  `No s.c.f. convergence criteria (sc_accuracy_*) were provided` 并自行挑
  `sc_accuracy_rho`（实测 `0.300000E-05`，随版本/体系漂移）、`sc_accuracy_forces`
  置 `not_checked` —— **BFGS 用的力没有独立精度保证**。现在与 `control_scf.in` 同源
  显式写出 `occupation_type gaussian 0.01`、`mixer pulay`、`n_max_pulay 10`、
  `charge_mix_param 0.2`、`sc_accuracy_rho 1E-5`、`sc_accuracy_eev 1E-3`、
  `sc_accuracy_etot 1E-6`、`sc_accuracy_forces 5E-3`（eV/Å，取几何阈值 1.e-2 的一半）、
  `sc_iter_limit 1000`（保留 FHI-aims 默认，不收紧成 SCF 模板的 500）。
- 开壳层（`multiplicity > 1`）新增 `fixed_spin_moment = multiplicity-1`（μB），把总自旋
  固定在弛豫过程中，避免 BFGS 步间自旋态漂移；闭壳层不写该行，行为与旧版一致。
- `smj.aims_out` 新增 `scf_settings(control_text)`（从 `control.in` 读回实际生效的
  `sc_accuracy_*`）与 `relax_summary` 的 `scf_converged_steps` / `scf_not_converged`；
  collect 遇到 "self-consistency cycle … not … converged" 直接 `StageError`，
  收敛标记次数少于弛豫步数只记 note。`results.json` 的 `molecule_opt` 段新增 `scf`，
  `input_summary.json` 也记录 SCF 精度与 `spin_moment`。

### 化学态对照臂 `smj variants`

- 本阶段优化的是中性硫醇/硫醚，而器件里锚定 S 与 Au 成键；两个 2CN-656 硫醚 demo 的
  S–Au 接触分别是 2.6877 Å / 141.5°（tip 落在硫醚）与 2.4267 Å / 111.4°（硫醇盐），
  差 0.26 Å / 30°，说明 tip 落点对初猜化学态很敏感。新增 `smj.variants` 与
  `smj variants --xyz molecule.xyz` 派生两条对照臂：`thiolate`（去 S–H 的 H，
  `charge -= 1/个`，电子数与多重度不变）与 `aucap`（S–H 的 H 换成沿原方向的 Au，
  默认 2.40 Å，电荷不变），写 `<name>.xyz` + `variants.json`（源文件 sha256、
  `config_hint`、`removed_h_indices`、`electron_count` 与口径 caveat：气相/0 K/无溶剂、
  裸阴离子 HOMO 可能高于真空能级、单原子 Au 不是真实电极）。只生成输入与清单，
  判定仍由下游链给出。

### species 审计与构象分支

- 新增 `smj species --dir <species_dir> [--xyz molecule.xyz] [--label light] [--out species.json]`
  与 `smj.species`：读 `*_default`（或 `<element>.in`）文件头里的
  `Suggested "<tier>" defaults`、统计生效基组行数与最高角动量、逐文件 sha256。
  `molecule_opt` prepare 逐元素写 `species.json`，collect 搬进 `results.json.species`；
  配置标签 `opt_basis` 与实际档位不一致时 notes 明确警告"以 species 文件为准"——
  `opt_basis` 从"一个词"变成可核对的口径。文档补了 light/tight 敏感性对照的做法
  （注意换档同时改基组与积分网格，不能分开归因）。
- 新增 `smj ensemble --smiles ... --dir ensemble/` 与 `smj.ensemble`：构象从"手动前端"
  变成批处理——生成 `conf_XX.xyz` + `ensemble.json` + `run_ensemble.sh`（脚本内
  `cd "$(dirname "$0")"`、路径相对脚本目录，便于整目录拷到集群），每个构象一条
  `smj stage molecule_opt --run --run-dir run/conf_XX`；`smj ensemble --dir ... --report`
  读各臂 `results.json` 给出 `gap_ev` 与 `max_force_ev_per_A` 的极差与缺失清单。
- `smj stage` 新增 `--run-dir DIR` 覆盖 run 目录（`StageRunner(run_dir=...)` 早已支持，
  现在从 CLI 暴露），构象/对照臂并行跑时各臂的 `run_spec.json` 与产物互不覆盖。

### 分子优化加收敛门与几何合理性闸门

- `molecule_opt` collect 之前只核"进程正常结束 + `geometry.in.next_step` 存在"，而
  FHI-aims 弛豫达到步数上限时仍会正常退出并写出 next_step —— 未收敛几何会以"已优化"
  身份进入 `aims_tip_relax`。现在 aims 分支强制要求 `aims.out` 出现
  `Present geometry is converged`，未收敛直接 `StageError`（消息带最后一步最大力、
  阈值、步数，并说明超步数也会正常退出的原因）。
- 新增 `smj.checks.check_geometry`（纯 python，无第三方依赖）：原子重叠（最近对
  <0.7 Å）、键图对比（重原子键消失 = 断键/解离、孤立原子、元素组成变化 → 拒绝；
  含 H 键消失 = 质子迁移 → 警告）、平移 RMSD 只报告。两种引擎共用，且**在落盘
  `opt.xyz` 之前**执行，不通过就不产出（Gaussian 分支会先删中间几何）。
- 新增 `smj.aims_out`：`relax_summary`（收敛 / 最大力 / 步数）与
  `frontier_orbitals`（最后一张 Kohn-Sham 特征值表的 HOMO/LUMO/gap，分数占据时报
  flag）。`results.json` 的 `molecule_opt` 段新增 `converged` / `max_force_ev_per_A`
  / `relax_steps` / `geometry_check` / `frontier_orbitals`，notes 里的 gap 明确标注
  "PBE 气相值，只作输运对照，不是器件能隙"。
- **口径修正**：`molecule.opt_force_threshold` 的单位是 **eV/Å**（FHI-aims
  `relax_geometry` 的判据是最大力分量），不是 Hartree/Bohr；默认 `1.e-2 eV/Å`
  ≈ 1.9e-4 Ha/Bohr 已经偏紧，不改默认值。加载配置时校验 `(0, 0.1]`，非数字/超范围
  直接 `ConfigError`；文档里所有 "Ha/Bohr" 字样同步改为 "eV/Å"。
- 新增 `smj conformers` 子命令与 `smj.conformers`：SMILES → ETKDGv3 生成候选 →
  MMFF94s 优化 → 按能量升序 + 重原子骨架 RMSD(0.25 Å) 去重 → top-k 写
  `conformer_<rank>.xyz` + `conformers.json`（第 1 个同时写 `--out`，可直接交给
  `smj init --xyz`）。默认 seed 42：RDKit 2026.09.1 里 `randomSeed=0` 会让所有构象
  用同一初猜（两两 RMS=0），故拒绝 0。MMFF 只作粗筛，manifest 的 `caveat` 写明
  "top-k 应各跑一次下游链再比较"。
- 新增锚定轨道权重（`smj.cube`）：解析 `output cube eigenstate` 写出的
  `cube_<n>_eigenstate_<i>_spin_<s>.cube`，按 `control.in` 的请求顺序标注
  HOMO-1/HOMO/LUMO/LUMO+1，给出**半径 1.6 Å 球内 |ψ|² 占全箱的比例**，写进
  `results.json` 的 `frontier_orbital_weights` 并打一条 note。权重高 = 该前线
  轨道在锚定原子上有振幅（有耦合通道）。只算占比，绕开 FHI-aims "内容单位
  1/Å³、体素单位原子单位" 的口径坑；有负值按带符号 ψ 平方、全非负按密度处理；
  没有 cube 只记 note 不报错。
- 新增电荷/多重度奇偶校验（prepare，两引擎共用）：电子数 = ΣZ − charge，偶数
  电子必须奇多重度、奇数电子必须偶多重度，否则 `StageError`。写错时 SCF 会为
  "不存在的电子态"收敛出漂亮结果，链上没有别处能发现；`results.json` 记
  `electron_count`。
- 修 `_vendor/junction/molecules.py::molecule_from_smiles` 的多构象 bug：
  `n_conformers>1` 时所有构象都返回构象 0（循环外只取了一次 `mol.GetConformer()`）。
- 新增锚定接触报告（只报告、不判定）：`smj.checks.anchor_contacts` 给出每个贴着
  金属的锚定原子（S/N/Se）的 `M–X` 距离与 `C–X–M` 角，`smj review` 写进
  `review.json` 的 `structure.anchors` 并打一条 info。硫醇盐在 Au(111) 上通常
  ≈2.4–2.5 Å / 105–115°；明显更长（≳2.6 Å）或接近线形（≳140°）说明锚定轨道与
  表面耦合弱（2CN-656 硫醚锚实测 2.69 Å / 141°，与 G/G0=0.00132 一致）。合格区间
  随锚定化学与位点变化，故不进闸门、交人工判读。
- 测试：`tests/test_aims_molecule_opt.py` 夹具补收敛横幅/最大力/特征值表，新增
  未收敛、断键、重叠、缺 `geometry.in` 四条 + config 阈值三条；新增
  `tests/test_geometry_checks.py` 与 `tests/test_conformers.py`（RDKit 缺失时 skip）。
- 文档：README 分子优化段、`docs/smj-计算流程梳理.md` 步骤 2 的闸门与单位口径同步。

### 输运能量窗改 eV 相对费米能级（默认 E_F ±2 eV）

- `transport.energy_window` 的默认语义从 **Ha 绝对能量** 改为 **eV 且相对费米能级**：
  默认 `[-2.0, 0.005, 2.0]` = `E − E_F ∈ [−2, +2] eV`、步长 5 meV（801 点）。
  换金属/换分子后 E_F 移动，相对窗口不跟着漂，不再需要手动按 Au Fermi 估窗口。
- `aitranss_te.prepare` 现在先读 `aims_scf/aims.out` 的 efermi（`parse_efermi`），
  把 eV 相对窗口平移到 `tcontrol` 需要的 Ha 绝对 `$ener/$estep/$eend`
  （`HARTREE_EV = 27.211386245988`）；缺 efermi 直接 `StageBlocked`，不用占位值。
- 新增 `transport.energy_window_unit`（`eV` 默认 / `Hartree` 旧口径）。旧 yaml 里
  跨度 ≤0.1 Ha、量级 <2 Ha 的三元组按旧 Hartree 口径识别并给迁移 warning，
  行为与 0.2.0 逐字节一致；要显式保留旧口径可写 `energy_window_unit: Hartree`。
- `transport.ener/estep/eend` 仍是 **Ha 绝对能量** 的单分量覆盖（优先级最高）。
- `protocol_sha256`：eV 口径把用户写的三元组 + unit 入协议；旧 Hartree 口径不再新增
  key，历史 run_spec 不会仅因本改动被判漂移（无需重建 run）。
- 新增 `tests/test_energy_window.py`（12 条）：默认相对窗口、efermi 平移、旧 Ha 探测、
  unit 显式优先、分量覆盖、protocol 指纹形状、`aitranss_te.prepare` 真的把 eV 口径写进
  tcontrol、缺 efermi 阻断。
- 文档：README 输运段、`docs/smj-计算流程梳理.md`（总览/步骤 6/§3.4/§8/§10）、
  demo/azulene README（历史备注）同步；`smj init` 模板改 `[-2.0, 0.005, 2.0]` +
  `energy_window_unit: eV`。

### 连接核对（原 A2 / 表面核对）更名 + 默认自动登记

- 「A2 审核」→「表面核对」→**连接核对**（connection check）：它核对的是电极与分子之间的
  连接/表面三元组是否与结构一致，不再用任务卡内部代号。用户可见措辞（CLI help、提示、报告、文档）
  全部换成「连接核对」；`表面核对` 与 `A2` 都降为旧名兼容。
- **默认自动登记**：新增 `review.auto_connection_check`（默认 `true`）。`electrode_build` collect 之后，
  管线自动跑表面审计并写 `electrode_build/connection_check.json`（`who="auto"`）；
  审计不通过则**拒绝登记**，`smj review` 仍会 blocked。想要「必须人工看图再签」，
  设 `auto_connection_check: false`，仍可用 `smj review --approve --who NAME`。
- 配置键更名（旧键仍生效，加载时提示改名）：`review.require_a2` / `review.surface_check` →
  `review.connection_check`，`review.a2_deviation` / `review.surface_check_deviation` →
  `review.connection_check_deviation`，`review.auto_surface_check` → `review.auto_connection_check`。
- 产物文件更名：`a2_confirmation.json` / `surface_check.json` → `connection_check.json`，
  `a2_surface_audit.json` / `surface_check_audit.json` → `connection_check_audit.json`，
  `a2_surface_viz.png` / `surface_check_viz.png` → `connection_check_viz.png`；
  旧文件仍会被读取（签署记录、审计、review 闸门均兼容）。
- 复核码与 check 名更名：`a2_*` / `surface_check_*` → `connection_check_*`；`validate` 的
  `junction_a2_confirmation` / `junction_surface_check` → `junction_connection_check`（旧名兼容跳过）。
- CLI：`smj review --connection-check/--no-connection-check` 取代 `--surface-check` / `--a2`
  （旧开关保留为隐藏别名，只覆盖本次调用、不写回配置）。
- 新增 `tests/test_connection_check_approve.py`、`tests/test_cli_connection_check.py`、
  `tests/test_runner_connection_check.py`，并删除 `tests/test_cli_a2_switch.py`。

### 事后补登记判据族不再阻断复核（`family` 漂移容差）

- 背景：`molecule.family` 只决定判据适用族，但计入 `protocol_sha256`；run 完成后再补登记族
  （例如给 BPY 标 `family: BPY`）会让 `smj review` 报 `provenance_drift` 并 blocked。
- 现在复核侧做精确容差：把当前配置的 `family` 改成 run 时值重算 `protocol_sha256`，
  若与 run_spec 一致 → 记 `warning`（`provenance_family_relabel`），按当前标注重算判据，不阻断；
  其余任何物理参数差异仍为 `error`（`provenance_drift`）。`run_spec.json` 不重写，审计链不变。
- `tests/test_review_family_drift.py` 覆盖：仅 family 不同放行、tip 距离变化仍阻断、
  缺 config 快照 fail-closed。

### 判据族 `BPY`

- `SYSTEM_FAMILIES` 登记 `BPY: "Au-(4,4'-bipyridine)-Au"`（吡啶 N 锚定，非硫醇盐）。
  `t_ef_magnitude` / `g_in_theory_band` / `seebeck_sign` 的标定域是 Au–BDT–Au，
  对 BPY 记 **not applicable**（不再 fail-safe 误判「否证 H0」），只作参考对照
  （实验 1.3±0.1 MΩ ≈ 9.93e-3 G0）。
- `runs/bpy` 重新生成报告后：工程判定通过、物理判定「支持 H0（passed）」、可进入发布/预测卡。
- 测试：`tests/test_family.py` 新增族登记与「BDT 带不适用」用例；全套 `pytest` 162 passed。

### 电极构建参数化（`electrode_build` → 「构建电极」）

- 阶段 4 不再写死「本地构建 Au(111) 电极」：**怎么构建**与**构建什么电极**都是本阶段参数。
  - `electrode.source`：`library`（FHI-aims `electrodes.library` 预定义团簇，默认）/ `builtin`（本地三角岛）。
  - `electrode.element`：电极金属（fcc(111) 单晶表面），登记 Au / Ag / Cu / Ni / Pt / Pd / Al；
    **tip 与电极同元素**，`aims_tip_relax` 由该键决定加哪种金属 tip。
  - `electrode.a_lattice`：fcc 晶格常数（Å），`null` = 用内置登记值（Au 4.07823 … Al 4.04950）；
    未登记元素必须显式给出（2.0–6.0 Å）。
- `source: library` 时若库团簇元素与 `electrode.element` 不一致，`electrode_build` **fail-closed**
  （`StageError`），避免「Ag 电极配 Au 团簇」这类静默错配。
- `electrode_metadata.json` 新增 `electrode_element` / `electrode_a_lattice`，仅在非 Au 或缺省晶格
  被显式覆盖时写入；**Au 默认路径的 metadata 仍逐字节不变**，旧 run 缺键时下游按 `Au` 处理
  （`self.energy.in`、`tcontrol` 表面索引、`smj review` 表面审计、`validate` 表面/锚定/电极范围检查）。
- 规范检查与报告：`validate` 的表面、锚定、电极范围检查改读 metadata 的电极元素；
  2.40–2.49 Å 的 S–Au 有源带**只按 Au 标定**，非 Au 电极只报接触距离、不作偏离判定（severity info）。
  `smj report` 新增「电极构建方式 / 电极元素」行。
- 兼容性：`electrode.element`/`a_lattice` 进入 `protocol_sha256`，旧 run 的 provenance 对比会
  报 `provenance_drift`（预期；Au 默认值下重建输入才可复现旧哈希）。
- 测试：新增 `tests/test_electrode_metals.py`（配置校验、非 Au tip、builtin Ag/自定义晶格装配、
  library 元素一致性、Au metadata 字节兼容）；全套 `pytest` 162 passed。
- 文档：README、`docs/smj-计算流程梳理.md`（步骤 3/4 与配置表）、`docs/README.md` 同步改为
  「构建电极」口径。

### 无硫分子锚定（任意元素 / 去质子三态）

- 锚定原子不再限定硫醇盐的 S：`electrode.left_anchor` / `right_anchor` 可写任意元素符号
  （结构里恰好出现 2 次）或 1-based 序号，分子里没有 S 也能建结（N/Se/… 锚）。
  `aims_tip_relax` 不再报「需恰好 2 个 S 锚原子」。
- `electrode.strip_thiol_h` 改为**三态**：`null`（新默认）= 仅当锚定原子是 S 时才自动移除
  S–H（旧 S 分子的行为不变）；`true` = 强制去质子；`false` = 不去。非 S 锚默认保留 H。
- 去质子 fail-closed：每个锚原子必须恰好 0/1 个 H、总计 0/2 个才移除；-NH2 这类一个锚原子
  带 2 个 H 的情况直接报错（不猜着只去一个）。
- `aims_tip_relax/input_summary.json` 新增 `anchor_element`（锚定元素），`strip_thiol_h` 记录
  实际去质子信息（`anchor_element`/`anchor_indices`/`route` 等）；注释与错误信息里的
  「Au–S」推广为「锚–金属」，tip/锚元素分别来自 `electrode.element` 与结构。
- `validate` 的 `junction_geom_anchor_aus` 改用 metadata 的 linker 索引识别锚定原子（任意元素）：
  S 锚保持原有硫醇盐带/硫醚分流逻辑；非 S 锚 severity info，只报告锚–金属接触距离，
  2.40–2.49 Å 带不适用。（旧 metadata 无 linker 索引时退回按 S 找锚。）
- 兼容性：`strip_thiol_h` 默认值由 `true` 变为 `null`，进入 `protocol_sha256`，旧 run 的
  provenance 对比会报 `provenance_drift`（预期；S 分子的实际去质子行为不变）。
- 测试：新增 `tests/test_no_sulfur_anchor.py`（三态解析、无 S 锚按元素/序号解析、默认保留 H、
  无 S 结校验给 info、S 分子仍自动去质子、-NH2 强制去质子 fail-closed、CLI prepare 端到端）；
  全套 `pytest` 147 passed。
- 文档：README（输入/适用/口径）、`docs/smj-计算流程梳理.md`（步骤 3 参数表、§3.3、判据表、
  §10 第 17 条）同步说明任意元素锚定与三态去质子。

## 0.2.0 — 2026-11-15

与 0.2.0rc1 相同的代码；本版补齐发布记录与真实体系验证。

### 验证（真实集群，FHI-aims 221103 + aitranss 071813）

- 三套体系跑通全链（molecule_opt → tip_relax → electrode_build → SCF → TE → Landauer）：
  BDT（86 原子，L0=0.279）、azulene-2,6-dithiol（92 原子，L0=0.381）、
  2CN-656-4-thiomethylphenyl（硫醚，128 原子，L0=0.00132）；工程判定全部通过，
  物理判定均为 inconclusive（`t_ef_no_peak_near_ef`）。详见 `docs/validation-0.2.0.md`。
- 网格收敛对照：BDT 与「冻结 S–Au=2.42 Å」两套几何上，estep 0.0002→0.00005 Ha
  （301→1201 点）后 L0/S/k_e/ZT_e 相对差 < 1e-8；`runs/_compare/BDT-对照实验报告.md`。
- 新增 `RELEASING.md`（版本 bump → 测试 → wheel/sdist → 安装冒烟 → tag/CI → Release）；
  README 增加「已验证体系」小节。

## 0.2.0rc1 — 2026-11-14

### 打包与发布工程

- 新增 MIT LICENSE，并写入 `pyproject.toml` 与 wheel METADATA（PEP 639 `License-Expression`）。
- 补全分发元数据：`readme` 长描述、authors、keywords、classifiers（Python 3.10–3.12 / Linux / Science）。
- `scripts/build_wheel.py` 支持 `--sdist`：同时产出可复现的 wheel 与 sdist（固定时间戳、排序、RECORD sha256 自检）。
- 版本号升至 0.2.0rc1；`smj verify` 在从 wheel 解包的环境下通过。
- 新增 GitHub Actions CI（`.github/workflows/ci.yml`）：pytest 矩阵 Python 3.10–3.14（Ubuntu）与 3.12（Windows），
  以及 wheel/sdist 构建 + 干净 venv 安装冒烟（`smj verify` / `smj init` / `smj plan`）。

### 验证

- 用户视角安装验证（Ubuntu 26.04 / Python 3.14.4，超出原 CI 矩阵）：wheel 与 sdist 分别装入干净 venv，
  `smj verify` 六项通过、`smj init` + `smj plan` 输出六阶段；`pip show` 正确解析 METADATA 2.4
  （`License-Expression: MIT`），`pip check` 无破损依赖；PyYAML 6.0.3 / numpy 2.5.4 均有 cp314 wheel，
  安装过程未触发本地编译。
- 本地全量 `pytest`（125 passed）在 Python 3.14.4 上通过；package classifiers 已加入 3.13/3.14，
  CI 矩阵同时覆盖这两个版本（3.13 未在本机单独跑过，由 CI 观察项兜底）。

### 默认行为（可能影响既有配置）

- `transport.energy_window` 的代码默认由 `(-2.5, 0.02, 2.5)` Ha 改为
  `(-0.21, 0.0002, -0.17)` Ha（围绕 Au Fermi ≈ −0.19 Ha，±0.02 Ha，201 点）。
  与 `smj init` 模板一致；生产仍应在 `aims_scf` 后按实际 E_F 调整窗口。
- `electrode.au_per_electrode` 在 `source: library`（默认）下被完全忽略，团簇大小由
  `electrodes.library` 文件决定；若显式写了该键会给出 warning。`source: builtin` 时
  `null` = 完整 60 Au/侧，整数按完整原子层截断（30 → 35/侧），上限 30。
- `smj preflight` 在 `molecule.opt_engine: aims` 时跳过 Gaussian 路径/许可检查
  （此前 G16ROOT/G16_MODULE 为空会误报 error）；改为 info 提示 `gaussian_not_required`。

### 自 0.1.2 以来的功能与修复

- 阶段规范化命名：`molecule_opt`（旧 `gaussian_opt` 别名）、`electrode_build`（旧
  `au60_build`/`auN_build` 归一化）；旧 run 目录仍可回退读取。
- 分子优化默认改用 FHI-aims（`opt_engine: aims`，PBE/light），与 tip/SCF 同一软件栈；
  Gaussian B3LYP/6-31G(d) 分支保留兼容。
- 默认使用 FHI-aims `electrodes.library` 预定义团簇（`source: library`），
  `electrode.nlayers: null` 取库文件头部建议；`source: builtin` 保留自建三角岛。
- 新增 `smj review`（A2 结构确认）：表面审计、渲染图、`geometry.in` sha256 指纹绑定；
  `review.require_a2` 默认 true，`--accept-deviation` 单次豁免，`--no-figure` 跳过出图。
- 新增 `smj review --a2 / --no-a2` 命令行临时开关（不写回配置、不进 `run_spec`）。
- 锚定化学判定：2.40–2.49 Å 只对去质子硫醇盐适用；硫醚/硫醇按 vdW 接触上限判定。
- 表面三元组确定性规则 det-v1（0.05 Å 分桶、坐标字典序），同一几何不再因运行次序换索引。
- 可复现闭环：`run_spec.json` 记录完整解析配置与协议哈希 `protocol_sha256`，漂移默认拒绝；
  各阶段记录输入/输出 sha256 与父阶段指纹。
- fail-closed 判据：缺 `TE.dat` 记为 `te_present(missing=True)` 并判 inconclusive，不再静默通过。
- 报告 `manifest.verdict` 区分 passed / failed / inconclusive，并给出 `allowed_action`。
- `smj analyze` 增加独立解析器交叉校验（Simpson + CODATA），相对差 ~1e-8。
- 修复 `sub_g16.sh` 不建 scratch 目录、Gaussian 坐标块缺空行、SCF 不收敛兜底
  （`scf=(xqc,maxcycle=512) int=ultrafine`）。
- 修复 mos/omat 回取清单只对 `aims_scf` 生效（`AIMS_HEAD_STAGES`）。

## 0.1.2 — 2026-11-11

- 库电极（`electrodes.library`）+ 背层固定弛豫（`au_relax`，主链路可选）。
- `smj report` / `smj analyze` / `smj audit` / `smj verify` 等命令完善。

## 0.1.1 — 2026-11-11

- 独立 wheel 打包、`smj verify` 安装自检、vendored junction 本地补丁（F1/F2/F4）。

## 0.1.0 — 2026-11-10

- 首个独立包：六阶段编排、状态机、Slurm/SSH 传输层、Landauer 后处理。
