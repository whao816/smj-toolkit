# Changelog

本项目遵循语义化版本（0.x 阶段：次版本可含不兼容变更；rc 视为发布候选）。

## 未发布

### 新功能:单分子结热导（声子）计算链(2026-11-23)

- 新增独立的四阶段链 `thermal_device` → `thermal_relax` → `thermal_fc` →
  `thermal_analyze`，与电子六阶段链(D12)完全隔离：smj run/plan/status
  的默认输出零变化，stage_range 仍只认六阶段。
- `thermal_fc` 是一条 Slurm 数组作业：中心差分 6N 个位移（3x3 周期板、
  器件区 54 原子时 324 个），每任务独立目录 + `apply_displacement.py`
  改几何 + `extract_forces.py` 取力块，最后一个完成的任务汇总写哨兵。
- 数值内核 src/smj/phonon/：单位换算（EV_PER_SQRT=0.0646542、
  THz=15.6333，与 phonopy VaspToTHz 相对差 1.4e-7）、López-Sancho
  表面格林函数、Caroli 透射（含 1/N_k 归一化）、Landauer 线性响应积分、
  动力学矩阵与虚频。63 条解析基准（完美链 T=1、带边 0.5、n 通道给 n·G_Q、
  G∝T、L↔R 对称到 1e-12）。
- smj thermal：新增子命令；`thermal.json` 契约 smj.thermal.v1 由
  smj.analysis.read_thermal 校验。
- 十条热导判据进 CRITERIA（ph_*）：带顶归零、量子上限、声学极限、
  ASR 残差、FC 对称性、虚频、lead 胞匹配、网格独立性等；4 条 physical 各带
  一手源。smj review 对热导 run 出 passed/inconclusive/failed，
  **缺证据一律 inconclusive，不当作通过**。
- smj report 新增 ## 5. 热导（声子） 节，并显式标注"热导器件（周期板
  3x3，~90 原子）与电子链器件（金字塔团簇）不是同一体系，G_ph 与
  G/G0 不可并列"。
- smj report 的 ## 5. 同时嵌入两联谱图
  （thermal_analyze/thermal_spectrum.png）与热流贡献分位表（累积 G_ph 达
  25/50/75/90% 的能量，meV 与 cm⁻¹）：G_ph 是积分，分位才能说明热由哪段
  声子带走。nap-sh 实测分位 5.8/10.9/12.2/14.3 meV 全落在界面骨架模区间，
  分子内 C-C 伸缩（~190 meV）在 lead 带顶之上、弹性透射为零。
- 新增 G_ph(T) 曲线（thermal_analyze/thermal_conductance_vs_temperature.png，
  thermal.json 的 g_ph_by_temperature 列 14 个温度点）：谐波近似下 T_ph(ω)
  与温度无关，只换玻色因子，**不额外花机时**。它同时是"300 K 够不够收敛"
  与"低温端是否回到 G∝T 量子极限"的免费检验。nap-sh：300 K = 13.72 pW/K
  已是高温平台（13.87）的 98.9%，10→20 K 由 2.04 升到 3.83 pW/K（比值 1.88）。
- 新增软模归属 src/smj/phonon/modes.py（describe_soft_modes）：虚频只报波数
  没有意义，要报"谁在动"。给出 原子/区域(分子/器件 Au/lead)/键 三类权重、
  参与比、刚体占比、最大键伸缩。**本征矢在质量加权坐标下**，必须除以 √m
  才能当物理位移用（否则 √14≈3.7 的 H 振幅会让每个低频模都像 C-H 伸缩）。
  nap-sh 的两个虚模（-21.30、-18.23 cm⁻¹）分子占 99.5%/99.8%、lead 占 0%、
  刚体占比 0.14/0.04、最大键伸缩 ~0.01 Å，且 |λ| 比典型实模浅 30–500 倍
  → 分子骨架的浅弯曲/扭转，不是键拉伸、也不是 lead 截断假象。把全部负本征值
  投影掉后 G_ph 只变 +0.20%。报告里给出归属表与投影对照。
- 判据口径两处修正（方案 D13/D14）：ASR 残差改判**相对残差**（1e-3）且只判
  器件区行（lead 两块按构造不闭合，判它永远失败）；网格独立性门限由 1%
  重标到 5%（实测梯形积分底噪就有 1-2%），真正的收敛控制量是 η。
- gpu1 端到端跑通（90 原子 3x3 周期板、324 个位移、约 17 h 墙钟）：
  G_ph(300 K) = 13.72 pW/K、G_Q = 283.93 pW/K、N_eff = 0.048、
  ASR 相对残差 4e-15、带顶以上 max T = 1.2e-7。
- 调度经验提升为项目级：新增 smj.slurm.gres_line，纯 CPU 作业默认
  --gres=gpu:0（Slurm 默认把节点全部 GRES 分给单个作业，不写会把 64 核
  整节点锁死），README §4.1 记录实测依据。
- 唯一内核改动：write_aims_geometry(path, atoms, constrained=(), lattice=None)；
  lattice=None 时输出逐字节不变（有逐字节回归用例）。PROVENANCE.json
  已刷新，gen_demo_reference.py --check 最大相对偏差 0.00e+00。

### 可复现性：方法说明与验证数据入库（2026-11-23）

- 方法说明 md/docx 的生成脚本从 `.tmp/`（被 gitignore）迁入 `scripts/`：
  `build_method_doc_docx.py`（docx 内容源，固定 core props 与 zip 时间戳）、
  `build_method_doc_md.py`（docx→md）、`check_method_doc.py`（重建并逐字节校验）。
  docx 可复现构建：两次构建 sha256 相同（`B21C3381…D2DB24A`，62702 B）。
- 2026-11-22 两组验证的六臂数据从本机 `runs/` 导出为入库副本 `docs/validation-data/`
  （222 文件，1091736 B；含 state/run_spec/几何/TE.dat/results.json 与 `MANIFEST.json`，
  白名单导入，不含 `aims.out`）。导出/校验脚本：`scripts/export_validation_data.py [--check]`。
- 报告与图的生成脚本同样入库：`scripts/report_light_tight.py`、`scripts/report_vdw_onoff.py`
  （原先只在 `.tmp/`），重跑得到相同报告与图。
- `runs/light-tight-2026-11-22/` 的两条 SCF 不收敛弃用臂移入 `_abandoned/` 并加 README，
  避免被误读为已完成的对照臂；四臂对照数据不变。
- 本轮不改生产代码、不跑新计算；默认测试仍为 620 passed / 2 deselected。

## 0.2.2rc2 — 2026-11-23

相对 0.2.2rc1（`04739ba`，24 个提交）：① 自由基/开壳层支线 S1–S8——Hirshfeld 自旋分区、
自旋分辨 Seebeck 与分通道判据、破缺对称局部证据与端到端集成回归、开壳层协议自检、
`smj spin` 汇总子命令；② 11-21/11-22 三轮复审修复（判定链指纹、run_spec 执行面哈希、
原子写、CLI 契约、打包与 P2 八项）；③ light/tight 与 vdW 开关整链验证记录、DFT 方法
说明 v1.2 入库。默认全量 620 passed / 2 deselected（`-m integration` 2 passed）。

产物（2026-11-23 构建）：`smj_toolkit-0.2.2rc2-py3-none-any.whl`（465405 字节，
sha256 `F6B6B0026C4DB2B618679AB8F77A85E7DF8135740949351759C19F09AA3813E2`）、
`smj_toolkit-0.2.2rc2.tar.gz`（856433 字节，
sha256 `4E8788EA753169C769F96DE9A8BB4D122A4D4B4462C32E03EF39E59F1B1E0AB8`）。

### 验证与方法文档：light/tight、vdW 开关与方法说明 v1.2（2026-11-22）

本轮只入库验证记录与方法文档，不改代码、不跑新计算；测试基线不变（620 passed /
2 deselected，`-m integration` 2 passed）。数据目录 `runs/` 不入库，仓内保留的是
记录、数字与图。

- **light/tight 整链对照**（[docs/validation-light-vs-tight-2026-11-22.md](docs/validation-light-vs-tight-2026-11-22.md)；
  数据 `runs/light-tight-2026-11-22/`）：AZ-dithiol 与 2CN 各两条完整臂（301 点、300 K），
  同分子、同协议、同能量窗，只换 species 档位；两档几何各自优化，差异含几何响应。
  AZ：L0 0.373290→0.486402（+30.3%）、S +53.35→+65.64 µV/K、ZT 0.11194→0.20060
  （+79.2%）；2CN：L0 0.001313→0.001361（+3.7%），但 S −33.53→−3.53 µV/K、
  ZT 0.02769→0.00046（−98.3%）。结论：light 档对 2CN 的 S/ZT 不具备定量可靠性；
  两个体系、一个能量窗不构成基组收敛。同目录另有两条 az-thiolate 臂停在 tip/SCF
  阶段（SCF 500 步不收敛、无 `results.json`），文档已说明不能当作“已跑”。
- **vdW 开关对照**（[docs/validation-vdw-onoff-2026-11-22.md](docs/validation-vdw-onoff-2026-11-22.md)；
  数据 `runs/vdw-2026-11-22/`）：同一生产协议只切 `vdw_correction_hirshfeld`（TS）。
  AZ 的 ZT 0.11194→0.23279（+108.0%），2CN 的 ZT 0.02769→0.02360（−14.8%）；
  两臂 S–Au 接触只差 0.002–0.007 Å——TS 是后自洽修正、不改 KS 哈密顿量，量到的是
  弛豫几何的响应，方向/幅度体系相关，不能外推成“vdW 抬高/降低 ZT”。
- **方法说明 v1.2 入库**（[docs/AutoSMJ_DFT方法说明_v1.2.md](docs/AutoSMJ_DFT方法说明_v1.2.md)，
  含 `.docx`）：主链设置（FHI-aims 221103 / PBE / light / TS-vdW / scalar ZORA；
  Aitranss NEGF σ=0.05/0.10 Ha、300 K）与“程序实跑过的算例”分账，§9.6/§9.7 收录上述
  两组对照，§10 列出未验证边界。v1.0/v1.1 移入 `docs/archive/` 保留版本历史。
- **文档/图**：`docs/figures/light-vs-tight-TE.png`、`docs/figures/vdw-onoff-TE.png`；
  `docs/README.md` 索引补入三份文档（旧版方法说明列在 archive），`RELEASING.md`
  测试基线 607→620 并注明 `-m integration` 的用法。

### 修复：2026-11-22 全仓代码复审（批 1–4）

用户要求按复审建议顺序修复。四批都补了回归测试；默认全量 620 passed / 2 deselected。

1. **判定链**：`manifest` 记录判据指纹（`verdict.py` + `junction/criteria.py` +
   `junction/validate.py` 的源码 sha256），判据实现变化或旧 manifest 缺指纹时
   `smj report` 重建判定而不是复用旧 `passed`；`unregistered_checks` 不再在渲染时丢失，
   `no_applicable_gate` 不再写成“未过的 warn 级判据”；`spin_channel_consistency` 在
   `junction_analyze/results.json` 布局下真正可读；TE.dat 含 NaN/inf 视为证据损坏
   （`te_present` missing → inconclusive），不再穿透成 `value=nan` 被判 failed；
   负透射保留告警与原始峰值记录。
2. **协议/编排**：`run_spec` 新增执行面指纹（species defaults / templates / 可执行文件
   哈希），改这些输入默认阻断并提示 `--force` 重绑；旧 spec 原地补齐基线且明说无法追溯。
   `sbatch` 超时/无 job id 抛 `SubmitUncertain` 并保留 `submit_unknown`，不再自动重提交
   两个作业；collect 失败保留 `finished` + `collect_failed`；缺件哈希写进 state 而不是
   静默省略；`aims_scf` 缺 `electrode_metadata.json` 且无 `--geometry` 时 fail-closed；
   `state.json` / `manifest.json` 改 mkstemp+replace 原子写。
3. **连接核对/报告/CLI**：指纹失效默认 `blocked`，只有显式 `accept_deviation` 才降为
   warning；`geometry.in` 缺失不再用 metadata 自存指纹自比；`report.md` 新增
   §0.2 结构复核/连接核对；`--json` 的进度日志改走 stderr；顶层异常给一行错误而非
   traceback；`smj report --strict` 在判定失败时退出 1；`smj review --stage` 真正透传。
4. **打包/文档**：`build_wheel.py` 在 Python 3.10（无 tomllib）上可读版本号；补
   `Provides-Extra: all` 与 pyproject 的 `all` extra；sdist 带 `docs/` 与
   `demo/demo-reference.json`；README/RELEASING 的 CI 矩阵、报告节数、`independent`
   键名、verify 的 legacy 说明与测试基线同步更新。

### 修复：2026-11-22 复审 P2（八项遗留）

用户点名处理复审清单剩余的 P2。全部补了测试（`tests/test_p2_20261122.py` 13 例），
全量 620 passed / 2 deselected。

1. **P_S 近抵消**：`Sα+Sβ` 相对通道尺度 ≤1e-6 时 `spin_polarization_S` 记 `None`
   并写 `spin_polarization_S_note`，`smj spin` 打印原因，不再输出 `+199999900%`
   这种由数值噪声放大的巨值。
2. **Au–S 带单一来源**：`criteria.AU_S_BOND_BAND = (2.40, 2.44)` 成为判据带唯一来源
   （Jackson 2009 2.42±0.02）；`validate.AU_S_BAND` 改为引用它；`contacts.TARGET_BAND_A`
   （2.40–2.49）明确标注为建结搜索目标带、不是判据带；preflight/validate/README 文案同步。
3. **cube 负步长**：`read_cube` 把降序轴翻成升序（新原点=旧原点+(n-1)·step、体素翻转、
   步长取反），修复 `searchsorted` 空窗导致 `anchor_weight` 静默 `fraction=0`；
   步长为 0 时显式报错。同时修正原点换算 Å 的单位混用。
4. **plot 参数契约**：`--labels` 数量与输入不一致时拒绝作图，不再 `zip` 截断静默丢曲线；
   多输入叠加时显式非 `auto` 的 `--spin` 报错，不再悄悄丢掉自旋参数。
5. **independent 复核进判定**：`validate_junction_run` 新增未注册检查
   `independent_recheck`（独立解析 + Simpson 对 L0/L1/L2，相对差 1e-3）；一致记 info，
   不一致按 unknown 把物理判定压成 inconclusive；无 TE.dat / 点数 <50 / 无管线矩阵时
   applicable=False。
6. **junction/review.py 接线**：`aims_tip_relax`、`electrode_build` 的 prepare 写
   `review_tips.xyz` / `review_junction.xyz` 检查点；`smj review` 读取并报告
   未确认/已确认/确认记录失效（`structure_check_*`）。几何不可解析时降级成 note，
   不阻断阶段。
7. **verify legacy-stages**：`smj verify` 新增 warning 级检查，确认
   `state.LEGACY_STAGES`（`au_relax`）的模块仍可导入，旧 run 迁移不会因删模块而崩。
8. **manifest 单一实现**：`junction.manifest.build_manifest` 在有 `aitranss_te` 时委托
   `report.build_manifest` 并投影回旧形状（含判据指纹、completeness、unregistered），
   消除两套 inventory/verdict 逻辑；无 TE 阶段仍只做 inventory、不伪造 verdict。

### 新增：2026-11-21 自由基支线 S8（`smj spin` 汇总子命令）

计划见 `docs/plan-radical-2026-11-21.md`（S8，用户点名后实现）。此前开壳层的通道数、
Sα/Sβ、Hirshfeld 分区与 review 提示散在 `results.json`、侧车和 `review.json` 三处，
看一条 run 要翻三个文件。

- `smj spin [RUN] [--json] [--temperature K]`：只读汇总一条 run 的自旋声明、双通道
  积分（α/β 点数、Tα/Tβ(E_F)、P(E_F)/P(L0)、Gα/Gβ/Gtotal、Sα/Sβ/S_total、P_S）、
  Hirshfeld 分区与分子内 +/− 自旋矩、review 已有的自旋 code，并附口径提醒；不写任何
  文件（需要新结论仍跑 `smj review`）。
- 旧 run 的 `results.json` spin 块没有分通道 S 时，优先从 `TE.dat` 重算
  （`spin_source` 标明来源）；`aims_scf/spin_partition.json` 不存在时明说"自旋分布:
  未采集"，不补零。
- `--json` 输出同一 payload（`run_id`、`molecule`、`open_shell`、`spin`、`partition`、
  `findings`、`review_codes`、`caveats`），便于脚本消费。
- 测试：`tests/test_spin_cli.py` 5 条；默认全量 571 passed，`-m integration` 2 passed。

### 新增：2026-11-21 自由基支线批次 D（开壳层协议自检 + 判据口径）

计划见 `docs/plan-radical-2026-11-21.md`（批次 D = S5；S8 已在上一条单独实现）。此前开壳层 run
被判 `inconclusive` 时只写"无法判定"，读者容易当成数据不足；同时"该采的自旋证据有没有
采"没有任何工程检查。

- `physics_verdict(checks, *, open_shell=False)`：当归约结果为 `no_applicable_gate` 且
  `open_shell=True` 时，`wording` 追加"（开壳层判据未校准；如需物理解读，需先用开壳层
  参考体系标定 `g_in_theory_band` / `seebeck_sign`）"。默认调用（闭壳层/未传）文案不变。
- `validate_junction_run` 新增 `extra_outputs` 参数与 engineering 检查
  `open_shell_protocol`（warn）：开壳层（`multiplicity>1` 或 `broken_symmetry`）时核对
  `results.json` 的 `spin` 块存在、`channels==2`、`transport.extra_outputs` 含
  `hirshfeld` 或 `spin_density`；缺项逐条写进 detail，齐全记 info，闭壳层记"不适用"。
  该检查带 `physics: False`，只进工程判定/报告，不进 `physics_verdict` 归约——否则一条
  工程 warn 会把 `no_applicable_gate` 挤掉，把"未校准"伪装成"数据不足"。
- `report.build_manifest` 从 `run_spec.config.transport.extra_outputs` 取快照传给
  validate；不新增任何数值阈值，也不把未校准判据升级为 error。
- 测试：新增 `tests/test_open_shell_protocol.py` 6 条（齐全/缺 spin 块/缺 extra_outputs/
  闭壳层不适用/默认文案不含标定句/开壳层文案含标定句）；`test_verdict_gates.py`、
  `test_open_shell_gates.py`、`test_criteria_failopen.py` 保持绿。默认全量 566 passed /
  2 deselected；`-m integration` 2 passed；改过 `src/smj/junction/criteria.py`，
  `gen_junction_provenance.py` 已刷新（36 文件），`gen_demo_reference.py --check` 零偏差。

### 新增：2026-11-21 自由基支线批次 C（BS 局部证据 + 端到端集成回归）

计划见 `docs/plan-radical-2026-11-21.md`（批次 C = S4 + S6）。此前破缺对称（BS）单重态
只有"反号初猜 + 两条自旋通道"，没有局部自旋证据；真实开壳层 run 也没有任何端到端回归。

- `spin_partition.summarize` 新增全表 `net`/`max_abs` 与每分区 `n_positive`/`n_negative`/
  `net`/`max_abs`，并给出分子组的 `localized_antiparallel`（正负共存、单个原子至少
  0.05 μB 且 `|net| < 0.05·max_abs`；0.05 μB 下限挡住数值噪声）；无 `molecule` 分区时
  为 `None`，不把"缺证据"写成"未检出"。`runs/ifs-bs` 实测：总自旋矩 0.0017 μB、
  分子 78 原子净矩 +0.00003 μB，但局部矩仅 ~5e-5 μB → 判为 `False`（该 BS 路径实际
  收敛到近闭壳层解），报告如实写"未满足局部反平行判据"。
- `smj report` 的 §4.1 在 BS run 的自旋分布表下加"破缺对称证据：分子内 +/− 自旋矩各
  N₁/N₂ 个，净矩 X μB"一行，并加脚注说明 FHI-aims 单行列式不提供 ⟨S²⟩/NOON、如需严格
  验证要用 Gaussian 模型分子交叉验证（默认不自动执行，也不自动跑闭壳层对照）。
- `smj review` 对 BS run 记 info `broken_symmetry_local_moments`（带局部计数与净矩）；
  非 BS 开壳层不出现该代码，不改任何判据结论。
- 新增 `tests/integration/test_radical_runs.py`（`pytestmark = pytest.mark.integration`）：
  真实 `runs/ttf-cation` 资产断言六阶段齐全、`channels==2`、|P(E_F)|<0.05、分子自旋占比
  0.56%；在 run 副本上离线重算 `spin`/分区后渲染报告并重跑 review（出现
  `open_shell_uncalibrated`/`spin_localized_off_molecule`）。`pyproject.toml` 注册
  `integration` marker，并用 `addopts = ... -m "not integration"` 默认排除。
- 验收：默认全量 560 passed / 2 deselected；`python -m pytest -m integration tests/integration -q`
  → 2 passed（1.7 s）。

### 新增：2026-11-21 自由基支线批次 B（自旋分辨热输运）

计划见 `docs/plan-radical-2026-11-21.md`（批次 B = S2）。此前开壳层只有
`P(E_F)`/`P(L0)`：知道两个通道电导不同，但说不出 Sα 与 Sβ 各是多少、总 S 是否等于按
电导加权的平均。

- `src/smj/junction/landauer.py` 抽出 `transport_moments(l0, l1, l2, T)` 作为
  S/k_e/ZT 的唯一公式源；`analyze_te` 与 `landauer_spin.analyze_te_spin` 共用，杜绝
  "两份公式各自漂移"。
- `spin` 块新增 `L1_alpha/L1_beta/L2_alpha/L2_beta`、`S_alpha_uV_per_K` /
  `S_beta_uV_per_K` / `S_total_uV_per_K`、`k_e_*`、`ZT_*` 与
  `spin_polarization_S = (Sα−Sβ)/(Sα+Sβ)`（分母为 0 给 None）；`junction_analyze`
  自动带出，`smj report` §4.1 加 Sα/Sβ 与 P_S 两行。
- 新增 engineering 判据 `spin_channel_consistency`（error 级）：校验
  `G_total = 0.5(L0α+L0β)`（相对 1e-9）与 `S_total = (L0αSα+L0βSβ)/(L0α+L0β)`
  （相对 1e-6），不满足即工程失败（说明分通道数不可信，不写成否证 H0）。闭壳层与
  旧 `spin` 块（只有 L0）标 `applicable=False`，不误判。
- 实测 `runs/ifs-bs`：Sα=+6.15599、Sβ=+6.15894 μV/K（相对差 4.8e-4），
  S_total=+6.15746，P_S=−2.39e-4（与 P(E_F)=−1.84e-5 同号）。判据在真实 run 上
  `passed=True`。改 `src/smj/junction/` 后已刷新 provenance（36 文件/13 补丁）并
  `--check` 一致，demo 参考 4 份谱零偏差；全量 555 passed。

### 新增：2026-11-21 自由基支线批次 A（Hirshfeld 自旋分区 + 单列 TE 闸门）

计划见 `docs/plan-radical-2026-11-21.md`（批次 A = S1 + S3 + S7）。此前的开壳层链能保证
"总磁矩保持、双通道透射还在"，但回答不了"这 1 μB 里有多少真的在分子上"；单列 TE 也会
静默通过。本批次补上证据与闸门：

- **S1 Hirshfeld 自旋分区**：新模块 `src/smj/spin_partition.py` 解析 `aims.out` 的
  `Performing Hirshfeld analysis` 段（取最后一段；缺 `Hirshfeld spin moment` 的原子记
  `None` 而不是 0），按分子 / 接触区（距锚原子 ≤3.0 Å）/ 电极分区。`aims_scf.collect` 在
  `transport.extra_outputs` 含 `hirshfeld` 时写侧车 `aims_scf/spin_partition.json`，
  `results["aims_scf"]["spin_partition"]` 放摘要；`smj report` §4.1 新增"自旋分布
  （Hirshfeld，μB）"表；`smj review` 对分子份额 <1% 的开壳层 run 记
  `spin_localized_off_molecule`（info）。实测 `runs/ttf-cation`：总 1.0000 μB，分子 1–14
  0.0056、tip 15/16 0.0063、电极 17–88 0.9881，与既有 Mulliken 投影一致到 1e-3。请求了
  `hirshfeld` 但段里解析不出 `| Atom N:` 行 → `StageError`（没有证据不标 collected）。
- **S3 单列 TE 默认阻断**：`aitranss_te.collect` 在 `multiplicity>1`（或破缺对称）且
  TE.dat 只有一列时直接 `StageError`——α/β 已被 aitranss 合并成均值，自旋极化无法恢复。
  新增 `transport.allow_single_channel: false` 显式放行，结果记
  `single_channel_allowed` 与 note。既有 7 个开壳层 run 的 TE.dat 全 4 列，不受影响。
- 文档：`docs/自由基单分子结-实现对照与开发规划.md` 状态表与 §5.1、`README.md` §9、
  `docs/plan-radical-2026-11-21.md` 同步更新。全量 549 passed；未动 `src/smj/junction/`，
  provenance/demo 基线不变。

### 修复：2026-11-21 复审 8 项问题（批 1–4）

复审报告（针对 `6384ce0d`）列出 8 个问题（3 条叙述 + 5 个收敛扫描反例），附 8 条断言；
在本仓 HEAD `a7f2a45` 上逐条复现后 8 条全部失败，按 `docs/fix-plan-2026-11-21-r3.md` 的
R1–R7 分四批修完。

**批 1：费米窗求积充分性。** 12 点常数谱（步长 0.2 eV、300 K、T=0.25）的权重积分只有
精确覆盖的 57.7%、权重一阶矩 −1.96 kBT，因此虚报 S=+169.25 µV/K、G/G0=0.1443（解析值
0.25）；端点覆盖判据（f(E_min)−f(E_max)=1.0）照样判 passed。`landauer.py` 新增
`fermi_quadrature(rows, temperature)`（±4 kBT 窗内的权重积分/精确覆盖比值、权重一阶矩、
窗内点数与最小步长），`analyze_te` 输出 `fermi_quadrature`；`criteria.py` 新增
`te_fermi_quadrature`（physical/error，`MIN_QUADRATURE_RATIO=0.95`、
`MAX_FERMI_MOMENT_KBT=0.25`；窗口本身不完整时一阶矩只作诊断），不足时 `missing=True`
把物理判定归约为 inconclusive。53 份真实 TE.dat 的比值偏差 ≤5.2e-8、一阶矩 ≤1.7e-4 kBT，
历史 run 全部不受影响；报告输运表新增"费米窗求积"一行。

**批 2：变体配置自洽 + 逐切片收敛。** `scan_convergence` 的变体 yaml 现在把
`molecule.xyz`、`cluster.species_dir`、`paths.*` 绝对化，把 `paths.species_dir` 规范化到
`paths.species_defaults`，写盘后回读自检（未知轴键、轴值解析不符、输入文件缺失都在 plan
阶段报错），`scan.json` 记录每个变体的真实 `run_dir`；`variant_metrics` 改读
`<run_dir>/junction_analyze/results.json`（旧数据回退 `analysis.json`，状态分
ok/empty/unreadable/missing）；`compare_axis` 只配对"其余轴取值相同"的切片并逐片比较，
缺结果、结果不可用或指标全空记 `incomplete`——不再把"没可比数据"的 worst=0.0 当收敛。
`config.py` 把 `paths.species_dir` 做成兼容别名：单给时告警并入 `species_defaults`，
两个键同时给报错。

**批 3：结构指纹与占据表配对。** 复核此前只比较协议子集，快照 config 里是 `molecule.xyz`
路径键、当前是 `xyz_sha256`，于是同一路径下改坐标也被降级成 `provenance_schema_evolved`
warning；现在先比对 `run_spec.json` 的 `molecule.sha256` 与当前文件哈希，不同即记
`provenance_drift` error（detail 带两个 sha 前缀），快照缺哈希或文件不可读时仍回落到
原 family/schema 逻辑（旧的 schema 演进用例保持 warning）。PDOS 的占据表不再取"最后两张
等长表"：闭壳层只取最后一张（旧行为把上一轮迭代当 beta，+0.2 eV 的权重 0 → 39.894），
开壳层取末两张表逐通道配对（行数可不同），占据表行数与波函数态数不符抛 `PdosError`
而不是按索引静默截断。

**批 4：8 条断言入库 + 文档。** 评审附件 `test_smj_review_6384ce0.py` 的 8 条断言逐条
移植为 `tests/test_review_20261121_regressions.py`（`workdir` 夹具，文件头注明
"附件测试 → 修复项" 映射）；`tests/test_review_20261121.py` 保留轴键校验、别名冲突、
开壳层配对/态数校验等补充断言。README 补上几何内容哈希闸门、PDOS 占位数配对规则，
以及基组轴/别名/汇总读取口径。全量 **533 passed / 0 errors**，
`gen_junction_provenance.py --check`（36 文件）与 `gen_demo_reference.py --check`（4 份谱）
均无偏差；附件原版 8 条在 HEAD 上直接复跑 8 passed。

### 修复：2026-11-20 复审 7 项问题（批 1–4）

复审报告（针对 0.2.2rc1 / `04739baa`）列了 5 个 P1 + 2 个 P2，附 8 条断言；在本仓 HEAD
上逐条复现后 7 项全部成立，按 `docs/fix-plan-2026-11-20-r2.md` 的 R1–R7 分三批修完。

**批 1：输运积分只有一个能量序，费米窗覆盖成为闸门。**
`landauer.py` 新增 `sort_by_energy` / `read_te_ordered(path) -> (rows_asc, order)`，
`criteria.py`、`spectra.py`、`pdos.py`（`energy_grid` / `plot_pdos`）、`plotting.read_curve`、
`aitranss_te.collect` 全部改走它；`analyze_te_spin` 分通道积分前同样排序，并记录
`input_energy_order`。旧行为下同一份 TE.dat 打乱行序，`g_over_g0_total` 会从 0.2498 变成
0.5413、峰检测还会凭空报出 +0.025 eV 的伪峰并判 inconclusive——两个数字是同一根因（按文件
行序 trapz / 相邻行当能量相邻）。新增判据 `te_thermal_window_coverage`
（physical/error，`MIN_THERMAL_COVERAGE=0.95`，即窗口至少覆盖 ±3.66 kBT 的费米导数权重）：
覆盖不足时判据 `missing=True`，物理判定归约为 inconclusive，截窗积分仍照常输出；报告输运表
加"费米窗覆盖 / 能量窗"两行。实测窄窗样例（±0.025 eV、300 K、T=0.25）原来 G=0.0586
（解析 0.25）、S=+38.4 µV/K 还判 passed，现在如实标"证据不足"。52 份既有 run 的覆盖率都在
0.9999 以上，没有一份被翻转。

**批 2：波的"存在"要按通道逐个校验，归档序号要能续接。**
`aims_scf.collect` 先由 `spin_state(config).open_shell` 或已落盘的 alpha/beta 决定是否自旋
极化，再对 `("basis-indices.out", "omat.aims", *required_wf)` 逐个 `probe` + `problem`：
零字节、头部被截断、规模 ≤0 都直接 `StageError`/`StageBlocked`，`has_mos/has_alpha/has_beta`
与新增的 `wavefunction_files` / `wavefunction_bytes` 都来自观测而不是 `is_file()`；
`aitranss_te._spin_nsaos` 也补了同样的表头校验。`structure_archive` 新增
`_recover_entries`（清单缺失或损坏时从现存快照重建条目，坏清单仍留 `.bad` 副本）、
`_next_seq`（清单序号与现存文件名的 max+1）与 `open("xb")` 排他创建：清单写坏后再归档
不会从 0001 重新开始、更不会覆盖首份结构，结果快照重名时追加 `.N`。

**批 3：工程判定与物理判定分账，`finished` 不再算完成。**
`verdict.apply_completeness` 现在分别累计 `engineering_errors` / `physics_failures`，
`engineering_passed = not engineering_errors`（物理失败不再拖垮工程判定），`warnings` 同时
认 `warn` 与 `warning`（`validate.py` 产出的正是 `warn`，旧口径把它们全丢了）；
`completeness_check` / `apply_completeness` 接受 `awaiting=()`，新增
`stage_acceptance(statuses, expected) -> (awaiting, missing)`：`collected`/`skipped` 算验收，
`finished` 单列"已采集待验收"，只有 `pending`/`running` 才算缺。manifest 的
`verdict.completeness` 带上 `stages_awaiting`，报告渲染"已采集待验收 / 缺少"两段、
"判据分账"一行，并把 skipped 阶段的 `note` 打印成跳过原因（无 note 时写"未记录原因"）。
`state.DONE_STATUSES` 保持不动（runner 的跳过逻辑依赖它）。

**批 4：复审断言入库，测试不再用 `tmp_path`。**
评审附件 `test_smj_review_20261120.py` 的 8 条断言逐条移植为
`tests/test_review_20261120_regressions.py`（夹具换仓库 `workdir`），并在文件头注明每条
对应哪一批。`tests/test_gaussian_spin.py` 最后 4 处 `tmp_path` 也迁到 `workdir`——Python 3.12
把 POSIX `0700` 落成本机不可读的 Windows ACL，`tmp_path`/`tempfile.mkdtemp` 必然踩雷。
顺带把 `tests/test_spin_transport.py` 的假产物补成带 `$uhfmo_alpha.aims` / `$uhfmo_beta.aims`
头行的真实形状（旧占位被新校验正确拒绝）。全量 **505 passed / 0 errors**
（此前长期是 476 passed / 4 errors），`gen_junction_provenance.py --check` 与
`gen_demo_reference.py --check` 均无偏差。

### 修复：开壳层与未登记族不再借用 BDT 的判据包络（2026-11-20）

`criteria.py` 里只有两条判据声明了适用域（`applies_to={'families': ('Au-BDT-Au',)}`）：
`g_in_theory_band`（派生包络 G/G0 ∈ [0.03, 0.40]）与 `seebeck_sign`（"S 为正"铁律）。
旧口径在**族未知**时 fail-open（"宁可多判一条"），于是 `runs/ttf-release`——family 未登记、
多重度 2 的自由基结——因为 G=0.2337 G0 落进 BDT 包络、S=+10.72 µV/K 为正，被写成
`physics_verdict=passed / 支持 H0`。两条数字碰巧落在别人的带里，不是这个体系的证据。

现在的口径（`criteria.applies(..., open_shell=...)`）：

- **开壳层**（`multiplicity > 1`）：判据族只针对闭壳层标定过，两条判据一律降为记录，
  detail 尾注写"不适用(开壳层未校准, 仅记录)"；
- **族未知**：没有该族的标定依据，同样降为记录（"不适用(体系族未知, 无标定依据, 仅记录)"），
  fail-open 取消；
- **族已知且匹配**（Au-BDT-Au 闭壳层）：仍作闸门（回归不变）。

`evaluate_run_criteria` / `validate_junction_run` 新增 `multiplicity` 参数；
`report.build_manifest` 从 `run_spec` 取族、从输入摘要取多重度传下去。
实测 `runs/ttf-release`：`physics_verdict` 由 `passed` 变 `inconclusive`
（`inconclusive_criteria=["no_applicable_gate"]`，`not_applicable_criteria` 两条判据，
`allowed_action="停; 须补数据或降级声明"`），复核仍 `errors: 0`、只剩两条 warning。
新增 `tests/test_open_shell_gates.py`（12 例）并改写 `tests/test_family.py` 里断言旧
fail-open 语义的用例。全量 **475 passed**。

### 修复：输运轴改由电极质心定向、协议 schema 演进不再误判漂移（2026-11-20）

用 `runs/ttf-release` 复核时出现两个 error，都不是物理结论变坏，而是判据口径缺陷：

- **表面索引被判到"错误一侧"**：`validate_junction_geometry` 原来用两锚定原子连线的主方向
  定输运轴。分子相对电极倾斜时，锚定连线会被分子内几何主导——实测 S6→S10 主要沿 y，
  而两电极质心沿 z 分离 19.66 Å，于是合法的 82 号表面原子（越过中线仅 0.25 Å）被判 error。
  现在**优先用左右电极质心连线**定向，锚定连线只作兜底；check 的 `value` 多了
  `basis: electrode_centroids | linkers`，报告里写明轴是由什么定的。
- **代码升级把老 run 判成协议漂移**：`protocol_sha256` 会因新增/改名的记录键（`stable_check`、
  `xyz`→`xyz_sha256`、`library_dir`、`tip_direction` 等）变化，原实现一律报
  `provenance_drift` error。现在 `run_spec.json` 增 `protocol_schema`，复核时逐键比较：
  两侧都有的键取值全同、`PROTOCOL_CRITICAL` 一个不缺、且差异只来自单侧键 → 降为
  `warning: provenance_schema_evolved`（写明 shared / only_run / only_now 计数）；
  只要物理参数变了、或关键键缺失，仍是 error。
- 实测 `runs/ttf-release` 复核：两个 error 消失，`status: warning`，只剩
  `open_shell_uncalibrated` 与 `provenance_schema_evolved` 两条 warning。
- 新增 `tests/test_transport_axis.py`（3 例：质心定向、锚定兜底在倾斜分子上复现误判、
  两者都缺时 warn 而不误判）与 `tests/test_review_provenance.py`（4 例：形状演进只告警、
  改温度仍 error、关键键缺失仍 error、比对函数本身）。全量 **463 passed**。

### 修复：让 PDOS 读得懂真实的 FHI-aims 波函数文件（2026-11-20）

P2 的 PDOS 原先只对合成夹具正确，拿 `runs/ttf-release`（88 原子、4348 基函数）实测时
每通道积分权重 ~6e20——三个真实格式差异一起发作：

- **头行不是 `$scfmo.aims`**：自旋极化写 `$uhfmo_alpha.aims` / `$uhfmo_beta.aims`；
  两者态头字母**都是 `a`**，所以通道只能由头行决定，不能由态头字母决定。
- **系数是定宽字段**：`format(4d20.14)` 每行 4 个、每个 20 字符、无分隔符，负号会粘在
  上一字段尾部（`…D-08-.752…`）。按正则 `findall` 切会把指数尾当新数，得到 `-6.1e-81`
  这类假系数——现在按声明的宽度切，`format(…)` 里的 `4`（每行个数）与 `20`（宽度）分开解析。
- **占据数不在波函数文件里**：改从 `aims.out` 的 `State Occupation …` 表读（自旋极化时
  末尾 α/β 各一张，共 6 张是各迭代的中间结果），核验后取最后两张等长的表。
- **NAO 基非正交**：`Σc²` 并不等于 1，也没有 Mulliken 意义；现在按每个态的 `Σ_all c²`
  归一（`w_g = Σ_{ν∈g} c² / Σ_all c²`），并在 `notes` 里写明"真正的 Mulliken 需要
  `omat.aims` 的重叠矩阵，本地只有表头"。
- 解析改为**流式**（`_scan_states` 逐态累加，不保留系数）：`ttf-release` 的
  383 MB `alpha.aims` + 383 MB `beta.aims` 全程 23 s、内存平稳。
- 实测：`smj pdos runs/ttf-release/ttf-release --groups molecule=1-14,left=15-51,right=52-88 --fetch`
  → 每通道积分权重 `molecule≈52.94`、`left≈1461.09`、`right≈1460.47`，合计 2974.5，
  与 `aims.out` 里该自旋通道的占据数之和逐位一致；图见 `junction_analyze/PDOS.png`
  （上栏 T(E)、下栏 α 实线/β 虚线，共享 E − E_F 轴）。
- `tests/test_pdos.py` 增到 15 例：定宽 `$uhfmo_*` 夹具、`Σc²=9` 也必须归一、
  占据数来自 `aims.out`、截断文件报错必须指向 `smj pdos --fetch`。全量 **456 passed**。

### 新增：分组 PDOS 与收敛扫描（P2，2026-11-20）

P1 补了"有没有自旋证据"，P2 补"这些态在哪里"和"参数改动是否改变结论"：

- `smj pdos RUN --groups molecule=1-20,left=21-52,right=53-84`：读 `aims_scf/mos.aims`
  （开壳层读 `alpha.aims`/`beta.aims`）与 `basis-indices.out`，把每个基函数的系数
  `c²` 按 `basis-indices.out` 的 atom 列归到原子、再按 `--groups` 汇总、乘占据数后做
  高斯展宽，落在 `TE.dat` 的 E − E_F 网格上——**与 T(E) 同一能量轴**（`--sigma` 默认 0.1 eV）。
  输出 `smj.pdos.v1`：`groups`/`channels{alpha,beta}`/`integrated_weight`/`state_energies_eV`/
  `pdos{group: [...]}`；`--plot` 出上下两栏共享 x 轴的 T(E)+PDOS 图（`[plot]` extra）。
  缺完整本征向量时**显式失败**（本地默认只回取 `mos.aims` 表头，实测 demo 的
  `mos.aims` 只有 112 B），错误信息写明 `0/4364` 与"请整取"，`--fetch` 可按显式清单
  从远端全量取回 `mos.aims`/`alpha.aims`/`beta.aims`。
- `scripts/scan_convergence.py plan|run|compare`：把"改一个参数重跑一遍"做成矩阵。
  `plan` 按 `--axis section.key=v1,v2` 展开变体、每个变体一份独立 `smj.yaml`
  （`run.dir` 固定 `run`），只写计划**不提交作业**；`run` 必须 `--yes`
  （遵守 AGENTS.md「只跑用户点名的计算」）；`compare` 读各变体 `analysis.json`，
  沿轴按取值顺序比较相邻档的 S/G/ZT/κe，最大相对变化小于 `--tolerance`（默认 5%）
  才算 `converged`，缺结果记 `missing` 并判 `incomplete`，**不把缺结果当 0**。
- 新增 `tests/test_pdos.py`（12 例，含 CLI 与 `--fetch` 的配置检查）与
  `tests/test_scan_convergence.py`（7 例：轴解析、变体展开、收敛/未收敛/缺结果、
  `run` 无 `--yes` 不执行）。

### 新增：开壳层的第一手证据（P1，2026-11-20）

P0 让报告"说出了"自旋通道，但证据仍只有 `TE.dat` 的两列透射谱：没有分原子磁矩、
没有自旋密度，也没有 Gaussian 一侧的自旋纯度。本次补两条数据链，默认全部关闭：

- `transport.extra_outputs: [hirshfeld, spin_density]`（白名单校验，写错即 `ConfigError`）：
  `aims_scf` 的 `control.in` 追加 `output hirshfeld` / `output cube spin_density`；
  collect 核验 `spin_density.cube` 非空、`aims.out` 确有 Hirshfeld 分析段（扣掉
  `vdw_correction_hirshfeld` 的回显），请求了却没产出就拒绝 `collected`；
  `results.aims_scf` 记录 `extra_outputs`/`spin_density_cube_bytes`/`hirshfeld_mentions`。
  闭壳层（`multiplicity=1`）要 `spin_density` 直接 `StageError`——FHI-aims 没有
  共线自旋密度可写。该关键字组合在 221103 上**尚未实测**，首次使用请在小体系确认。
- `molecule.stable_check: true`（仅 Gaussian）：route 追加 `Stable=Opt`；
  `parse_gaussian_log` 解析最后的 `S**2 before annihilation … after …`、
  `SCF Done` 次数与 `Convergence failure`；`results.gaussian` 记录
  `scf_converged`/`multiplicity`/`s2_before_annihilation`/`s2_after_annihilation`/
  `s2_expected`/`spin_contamination`，⟨S²⟩ 偏离 S(S+1) 超过 0.1 时在 notes 里点名
  自旋污染，缺成功 SCF 时也会点名。
- 两项都**不进** `protocol_sha256`（只决定多存哪些证据文件，不改求解参数）；
  默认配置的 `control.in`/`.gjf` 逐字节不变。新增 `tests/test_extra_outputs.py`（18 例）
  与 `tests/test_gaussian_spin.py`（10 例）。

### 新增：开壳层结果进入报告与复核（P0 展示层，2026-11-20）

此前 `report.py`/`review.py` 里 0 处提及自旋：P(E)、α/β 通道与「判据族未校准开壳层」
只躺在 `junction_analyze/results.json` 和图里，评审者从 `report.md`/`review.json`
看不出"这是开壳层、自旋占多少、解是否可信"。本次只做展示层，不改任何数值判定：

- 新增 `src/smj/spin_summary.py`：合并 `run_spec` 的 molecule 块与 `results.json` 的
  `spin` 块，`report`/`review` 共用同一份字段与三条 caveat（判据族未校准、固定总自旋
  只证明约束被满足、单列 TE.dat 无法恢复自旋极化）。
- `report.md`：§2 新增「电荷 / 多重度」与开壳层「总自旋约束」；§4 末尾新增
  `### 4.1 自旋通道（开壳层）`（通道数、Tα/Tβ(E_F)、P(E_F)、P(L0)、Gα/Gβ/G_total），
  并显式写明总 G、S、k_e、ZT 仍是双通道均值。manifest 的 `transport.spin` 同步保留。
- `smj review`：新增 `open_shell_uncalibrated`（multiplicity>1，warning）、
  `spin_fixed_constraint`（info）、`spin_channels_missing`（开壳层却单列谱，warning）、
  `spin_channels_undeclared`（旧 run 有两列谱但无 multiplicity，info）；
  `review.json` 新增 `spin` 块。前两条 warning 把 `review.status` 降为 `warning`
  （不阻断，`ok` 仍为 true）。
- 新增 `tests/test_spin_surface.py`（6 例）。全量 `pytest -q` → 406 passed；
  `smj verify` 七项全 ok。

## 0.2.2rc1（2026-11-20）

### 修复：换行一致性导致 CI 变红（2026-11-20）

首次推送后 CI 有两类失败，都不是数值逻辑问题：

- **Windows 三个版本全红**：`actions/checkout` 在 Windows 上按 `core.autocrlf=true`
  检出，内核 `.py` 与 `TE.dat` 被转成 CRLF，`PROVENANCE.json` 与
  `demo-reference.json` 里的 sha256 对不上（6 个用例失败）。
- **`demo` job 红**：三处原因 —— ① `demo/demo-reference.json` 曾以 CRLF 提交，
  Linux 上 `gen_demo_reference.py --check` 直接比对文本 → 必然不一致；② 该 job 只装
  运行时依赖（`pip install -e .`）却调用 `python -m pytest`，缺 pytest → 退出码 1；
  ③ 修掉前两条后仍红：`--check` 拿整份文本逐字节比，而 `L0`/`ZT_electronic` 是按
  完整精度写入的浮点数，`np.exp`/`np.cosh` 在不同平台的 numpy 构建（以及不同 CPU 的
  SIMD 分派）上有 1 ulp 级差异 → 什么都没改也会红。

处理：新增 `.gitattributes`（`* text=auto eol=lf`，`*.dat -text`），全仓
`git add --renormalize`；`gen_demo_reference.py` 写文件固定 LF；`--check` 改为
**按字段比较**（哈希/点数/能量顺序精确比，浮点量按相对 `1e-9` 比），失败时逐条打印
「哪份谱的哪个量、期望/实际、相对差」，CI 再把这段文字写进 check-run annotation
（匿名也能从 annotations API 读到）；新增 CRLF 容忍与「1e-13 噪声不红、1e-6 漂移仍红」
两个用例。修复后本地与两种检出模拟（autocrlf=true/false）均 398 passed。

### A 方案：junction 内核改由本仓维护（2026-11-20）

复核后的结论是「运行时早就独立、维护模型还假装是镜像」：13 条本地补丁覆盖
17/36 个内核文件（47%），`upstream_commit` 为空，CI 里没有同步脚本，6 个文件
曾从上游未提交工作区同步。本次按 A 方案收口：

- **目录与命名**：`src/smj/_vendor/junction/` → `src/smj/junction/`；
  `VENDOR.json` → `PROVENANCE.json`；`scripts/gen_vendor_json.py` →
  `scripts/gen_junction_provenance.py`；`smj verify` 的检查项改为
  `junction-kernel` / `kernel-provenance`。
- **同步机制退役**：删除 `scripts/sync_core.py`。内核不再与 AgenticDFT 逐字节同步，
  算法改动直接改本仓（AGENTS.md 已改）。
- **一致性契约改为结果**：新增 `scripts/crosscheck_agenticdft.py`，把同一份 TE.dat
  交给两边内核对比 L0/S/κe/ZT（手动跑、不进 CI，因为 CI 里没有主仓）。
  实测薁示例两边逐位一致。
- **R12 收口**：`k_e_zt_order` 元数据 `implemented` 改回 `True`（评估器确实会 emit），
  删除 `PROVENANCE.json` note 与 README §9 的「已知偏离」登记。
- **测试**：`tests/test_vendor_manifest.py` → `tests/test_junction_provenance.py`；
  全量 397 passed；`smj verify` 7 项全 ok。


### 复核修复 R12/R13：产物定位去歧义、已知偏离登记（2026-11-20）

- **R13b 产物定位**：vendored `validate.py` 的 `_ARTIFACT_STAGE_HOME['aims.out']` 原为
  `('aims_scf', 'aims_tip_relax')`。弛豫日志每一步都写
  `Self-consistency cycle converged.`，于是"静态 SCF 阶段根本没跑"的 run 也能拿到一条
  `aims_scf=OK` 判据——它证明的是弛豫里的某一步 SCF。收紧为只认 `aims_scf` 阶段
  （扁平布局仍由 `run_dir/aims.out` 命中），以本地补丁
  `F-RECHECK-2026-11-20-aims-out-stage` 登记。
- **R13a 两级一致**：`tests/test_scf_artifacts.py` 新增一条用例，同一份未收敛日志
  必须同时让 stage 抛 `StageError` 与 validate 报 error——任何一边漏掉即红。
- **R12 已知偏离登记**：vendored `criteria.py` 里 `k_e_zt_order` 的元数据仍写
  `implemented=False`（实际会 emit），与主仓新注册表不一致。该判据不参与闸门，
  按计划**不改 vendored 判定语义**，只在 `PROVENANCE.json` 的 `note` 与 README §9 登记，
  待回迁主仓时一并对齐。
- **测试**：`tests/test_validate_artifacts.py`（4 例）+ `test_scf_artifacts.py` +1。


### 复核修复 R9/R10：demo 回归基线、VENDOR 校验语义（2026-11-20）

- **demo 回归（R9）**：新增 `demo/demo-reference.json` 与 `scripts/gen_demo_reference.py`
  （`--check` 可当 CI 门），把 4 份随发布走的 TE.dat 的
  `sha256 / 点数 / 温度 / log10(G/G0) / L0 / S / k_e / ZT / 峰值` 固化成基线；
  `tests/test_demo_regression.py`（16 例）用 `analysis.compare_reference` 做容差比对，
  并证明"改一个数 / 整体放大 50% / 换谱"都会红。CI 新增 `demo` job（只装运行时依赖，
  不需要 matplotlib）。
- **demo 文档纠错**：`demo/azulene-26-dithiol/README.md` 两处「193 点」改为实际的 190 点。
- **VENDOR 校验（R10）**：`scripts/gen_junction_provenance.py --check` 改为逐字段比较并
  **忽略 `repo_commit`**——旧实现比整份文本，而 `repo_commit` 记录生成时所在提交，
  每次 commit 都让它变红（正是复核指出的"绑 commit 必红"）。
  `tests/test_vendor_manifest.py`（6 例）守住该语义，并校验补丁台账里登记的文件真实存在
  （写错文件名会让 `sync_core` 静默覆盖补丁）。
- **README**：§1 增加三段仓库级一致性检查命令；§8 说明 `repo_commit` 不参与比对。

### 复核修复 R8：采样覆盖进协议指纹（2026-11-20）

- **问题**：`provenance._protocol_config` 无条件 `pop("ener"/"estep"/"eend")`，而
  `config.resolved_window` 里这些覆盖**确实生效**——eV 相对窗下把 `transport.estep`
  改成 `0.001` Ha，`protocol_sha256` 不变，`resolved_window` 却变了。
- **改法**：只保留 `require_te_out` 的 pop（验证策略，不进协议）；`ener/estep/eend`
  先 pop 再「显式写了才入协议」。默认配置不加键，0.2.x 的哈希形状逐字节不变
  （demo azulene 仍是 `635ae0f0f1b7…`）；`config_to_dict` 的 `resolved_window` 旁标注
  `resolved_window_unit: "Hartree"`。
- **测试**：`tests/test_provenance.py` +4（`estep` 覆盖改哈希、起终点覆盖改哈希、
  默认不加键、run_spec 标注单位）。

### 复核修复 R7：运行没走完就不能"通过"（2026-11-20）

- **问题**：`manifest.verdict` 只看已存在产物的判据；只放一份 TE.dat 的空 run
  也能 `passed: true`（复核实测 `engineering_passed=true`、`errors=0`）。
- **改法**：新增 `src/smj/verdict.py`——`completeness_check()` 生成
  `run_incomplete`（engineering/error）判据并按"工程判据"口径重算
  `errors/warnings/engineering_passed/passed`，`missing_stages()` 以
  `collected/finished/skipped` 为完成；`report.build_manifest` 在
  `validate_junction_run` 之后套一层，`verdict` 增加 `completeness` 字段，
  `render_report` 增「完工度」行。
- **测试**：`tests/test_report_verdict.py` +3——未完成时 `errors==1` 且
  `failed_criteria` 含 `run_incomplete`、`physics_verdict` 不被改写、只放 TE 的
  run 判定为未完成。


### 复核修复 R6：远端仍在跑的作业不再被"重跑"打成 failed（2026-11-20）

- **问题**：`state=waiting` + `job_id` 时执行 `run`/`run_stage`，`prepare` 抛
  "正在远端运行"的保护异常，被 `run_stage` 的兜底 `except SmjError` 写成
  `status=failed`；下一次 `run` 于是重新 `prepare`+`sbatch`，在集群上留下第二个
  作业（复核实测新作业 424242 覆盖 12345 的登记）。
- **改法**：新增 `ActiveJobError(StageError)`；`prepare` 抛它，`run_stage` 单独
  捕获——**不改状态**，返回 `StageOutcome(status="waiting", ok=False)`，保留
  `job_id`。`run()` 对 `submitted`/`waiting` 走新的 `_resume_active()`：先
  `wait()`，结束再 `collect()`；`finished` 仍只补 collect；`--force` 语义不变。
- **测试**：`tests/test_rerun_active_job.py`（4 例）——活跃作业保留状态且不重复
  提交、`run` 走 resume 分支、`--force` 仍重跑、等待超时保持 `waiting`。

### 复核修复 R4/R5：弛豫与 SCF 必须"真收敛、真产物"（2026-11-20）

- **问题**：`aims_tip_relax.collect` 只判 `geometry.in.next_step` 存在，从不打开
  `aims.out`；`au_relax.collect` 同样不解析收敛标记却在 metadata 里硬写
  `converged: True`。`aims_scf.collect` 只 `is_file()`，零字节
  `basis-indices.out`/`omat.aims`/`mos.aims` + 未收敛日志都能 collected，且
  `normal_termination`/`has_*` 是写死的常量。
- **新增** `src/smj/stages/_aims_log.py`：`relax_status()` 取**最后一次**
  `Present geometry (is|is not yet) converged.`、要求日志以 `Have a nice day.` 结束、
  识别步数上限（排除 "No maximum number of relaxation steps, defaulting to 1000"），
  大日志只读尾部 1 MB。
- **新增** `src/smj/stages/_aims_artifacts.py`：按 FHI-aims 头部（`$basis.indices.aims`、
  `$omat.aims`、`$scfmo.aims`）校验产物真实性，顺带解析 `nsaos`/基函数数量。
- **阶段行为**：`aims_tip_relax` / `au_relax` 未收敛或日志缺失 → `StageError`（不采集）；
  `aims_scf` 缺 SCF 收敛标记、零字节或头部不符 → `StageError`；`results` 改为观测值
  （`scf_converged`、`basis_functions`、各文件字节数等）；`electrode_build.prepare`
  与 `aitranss_te.prepare` 各加一道上游复核（`StageBlocked`）。
- **测试**：新增 `tests/test_relax_convergence.py`（7 例）、`tests/test_scf_artifacts.py`
  （9 例）、`tests/_aims_fixtures.py`（真实形状的最小产物片段）；既有测试的
  `placeholder` 产物改为走 fixtures。
- **实测**：4 份 demo 的 `aims_tip_relax/aims.out` 全部判 converged，4 套
  `aims_scf` 产物全部通过头部校验（无回归）。

### 复核修复 R2/R3：能量窗必须跨 E_F，"无适用闸门"不再判通过（2026-11-20）

- **问题**：`transport.energy_window`（eV 口径）只校验 `eend > ener`，不要求窗口跨越 E_F；
  整窗落在 E_F 一侧时 `criteria` 仍判 `physics_verdict=passed`（`g_in_theory_band` /
  `seebeck_sign` 都 `applicable=False`，等于没有任何闸门被真正评估）。
- **上游/本包口径**：新增 `te_window_covers_efermi`（physical / error）：`E−E_F` 必须
  `min <= 0 <= max`，否则 `failed`；`physics_verdict` 在所有带 `applies_to` 的族闸门都
  `applicable=False` 时返回 `inconclusive`（`inconclusive_criteria` 含 `no_applicable_gate`），
  不再落到 `passed`。
- **本包守卫**：`aitranss_te.prepare` 在解析出 `efermi` 后校验窗口含 E_F，不含则
  `StageBlocked`（不生成 `tcontrol`），错误信息给出相对窗与 `efermi`。
- **测试**：新增 `tests/test_verdict_gates.py`（7 例：窗口在 E_F 上方/下方/跨 0、AZ 常数谱
  → inconclusive、BDT 有适用闸门、S≈0 不判通过）；`tests/test_energy_window.py` 增
  `test_prepare_blocks_window_not_covering_fermi`。
- **VENDOR 维护**：`criteria.py` 新增本地补丁 `F-RECHECK-2026-11-20-window`
  （PROVENANCE.json 现 36 文件 / 12 条补丁）。

### 复核修复 R1：TE 谱行序不再影响 Seebeck 符号（2026-11-20）

- **问题**：`analyze_te` 按 `TE.dat` 行序直接梯形积分，`L1` 有符号而 `L0/L2` 取绝对值；
  能量列写成降序时 L1 反号 → Seebeck 静默变成反号（demo 薁谱 +25.94 ↔ −25.94 µV/K）。
- **上游修复**（`agenticdft/agenticdft/tools/junction/landauer.py`）：新增
  `energy_order(rows)`；`analyze_te` 发现非升序时按 `E−E_F` 排序后再积分，并返回
  `input_energy_order` / `reordered` / `warnings`。`independent_check.l0_l1_l2` 同口径排序。
- **本包**：`src/smj/analysis.py` 保留 vendored 的 `warnings`（旧实现直接清空，等于把
  "已被静默排序"这条唯一提示也吞掉）；`aitranss_te.collect` 记 `energy_order`，非升序时
  在 notes 里写明。新增 `tests/test_te_order.py`（5 例，含顺序无关参考积分对照）。
- **VENDOR 维护**：`scripts/sync_core.py` 落地（`--check`/`--list`/`--files`，
  默认跳过已登记补丁的文件）；`PROVENANCE.json` 新增 `F-UPSTREAM-PENDING-2026-11-20`
  条目，登记 7 个"副本领先于主仓 HEAD"的文件；CI 增加 `gen_junction_provenance.py --check`。

### 结级自旋初猜与释放约束（2026-11-19）

- `molecule.initial_moments` 现在会映射到结几何：`aims_tip_relax`/`aims_scf`
  把分子坐标系的 1-based 序号写到 `geometry.in` 对应 `atom` 行的 per-atom
  `initial_moment`（`strip_thiol_h` 删掉锚定 S–H 后自动左移；指向被删氢会报错）。
  有逐原子初猜时不再写全体统一的 `default_initial_moment 0.1`。
- 新增 `molecule.fixed_spin_moment: true|false`（默认 `true`）：`false` 时
  `control.in` 不写 `fixed_spin_moment`，只留 `spin collinear` + 初猜，用于
  "解除总自旋约束"的 release 对照；`input_summary.json` 记
  `spin_constraint: released`。`smj/spin.py` 的 `SpinState.fix_spin_moment` 统一
  两处 stage；`protocol_sha256` 对默认值不敏感（只在显式 `false` 时入协议）。
- 对照实测（`runs/ttf-fixed` vs `runs/ttf-release`，同一 88 原子几何）：两解都
  收敛；fixed（S=1/2）总能低 3.88 meV，但 1 μB 的 98.9% 在 Au 电极、分子 1–14
  仅 0.0055 μB；release 解收敛到 S≈0（分子 −0.0001 μB，gap 0.89 meV）。详见
  `runs/ttf-spin-comparison.md` 与开发规划文档 §5.1。
- 接触构型复核（2026-11-20）：`tip_direction: lone_pair` 重建（review 通过）后
  重跑 fixed/release，结论相同——ΔE 仍为 +3.88 meV、分子局域自旋仍 ≈0；
  四套 SCF 见 `runs/ttf-lone-pair*`。

### 自由基结实现对照文档（2026-11-19）

- 新增 `docs/自由基单分子结-实现对照与开发规划.md`：把外部调研文档
  （Gaussian/FHI-aims/AITRANSS 自由基输运核查整合版）翻译成 smj 口径，
  逐阶段对照自旋/电荷/输运链，列出与源文档一致处、口径差异（整结电荷 vs
  片段电荷、固定总自旋、initial_moment 粒度、OpenMP aitranss 等）、
  未实现项与后续优先级；`docs/README.md` 与根 README §9 已加链接。

### 开壳层分子输入修正（2026-11-19）

- 修 `molecule_opt` 的 FHI-aims `spin` 关键字：原先按多重度写 `spin 2` 这类数字，
  而 FHI-aims 的 `spin` 只接受 `none`/`collinear`（数字会因关键字取值非法而停止）。
  现在 `multiplicity > 1` 写 `spin collinear`，总自旋仍由
  `fixed_spin_moment = M−1`（=2S，不是多重度）固定。
- 新增 `molecule.initial_moments`（1-based 原子序号 -> μB）：写成 geometry.in 的
  per-atom `initial_moment` 子 tag（独占一行，不影响 `read_aims_geometry` 解析）；
  FHI-aims 在 `spin collinear` 下要求给出初始自旋密度。未给时 notes 警告；
  闭壳层设置或序号越界直接报错。
- `aims_tip_relax` / `aims_scf` 对 `multiplicity > 1` 先 fail-closed（当时结级自旋
  尚未接线；现已按下一节「结级自旋接线」打通）。
- 非零 `molecule.charge` 时把总电荷均分成 per-atom `initial_charge` 写进
  geometry.in：FHI-aims 要求 `charge` 等于这些行之和；项目使用的 221103 不认
  `override_initial_charge_check`（2026-11-19 集群实测 parse_control 停算），
  所以改用满足一致性检查的 per-atom 写法。写入值记入 `input_summary.json`
  的 `initial_charge_total/initial_charge_per_atom` 与 notes。
- `protocol_sha256` 只在 `initial_moments` 非空时纳入，历史 run_spec 不判漂移。
- 新增 `tests/test_open_shell_junction.py`、更新 `tests/test_aims_molecule_opt.py`；
  全量 306 passed。

### 结级自旋接线（2026-11-19）

- 新增 `src/smj/spin.py`：`SpinState` 统一分子/结的自旋口径（`multiplicity>1`
  → `spin collinear` + `fixed_spin_moment=M−1` + `default_initial_moment 0.1`
  初猜；`charge≠0` → 把总电荷均分成 per-atom `initial_charge` 写进 geometry.in），
  并负责把关键字插进 control.in 的 general 段、写入 `input_summary.json`。
- `aims_tip_relax` / `aims_scf` 不再 fail-closed：prepare 把 `spin`/`charge`/
  `fixed_spin_moment` 传进 FHI-aims 输入；开壳层 SCF 产出 `alpha.aims` +
  `beta.aims`（collect 接受两者替代 `mos.aims`，`results.json` 记
  `has_mos/has_alpha/has_beta/spin_unpolarized`）。
- `aitranss_te` 开壳层时从 `alpha.aims` 解析 `nsaos`（自旋极化无 `mos.aims`），
  把 tcontrol 的 `$scfmo` 换成 `$uhfmo_alpha file=alpha.aims` +
  `$uhfmo_beta file=beta.aims`；带电结按电荷修正 `$valence_electrons`；
  `input_summary.json` 记 `spin_polarized/wavefunction_files/valence_electrons`。
- 新增 `src/smj/landauer_spin.py`：解析 `TE.dat` 的 α/β 双列，给出 per-channel
  `L0`/`G/G0`、`E_F` 处与 L0 积分后的自旋极化率；`junction_analyze` 在开壳层
  时把 `spin` 块并入 `results.json`（总 G/S/ke/ZT 仍按双通道均值）。
- 新增 `tests/test_spin_transport.py`，重写 `tests/test_open_shell_junction.py`；
  全量 315 passed。口径限制（无 SOC、非共线、Kondo；平均透射不是自旋分辨结论）
  写入 README §9。
- 阶段间几何核对改用**结构指纹**（只哈希 atom/lattice_vector 行）：带电/开壳层
  SCF 会往 geometry.in 追加 `initial_charge`/`initial_moment` 子 tag，字节 sha256
  变化会让 aitranss_te 误判「几何与 metadata 不一致」。`electrode_build` 的
  metadata 新增 `surface_fingerprint_structural`；旧 run 缺该键时仍按字节
  sha256 回退（2026-11-19 TTF 实跑 prepare 时发现）。
- 远端回取补 `alpha.aims`/`beta.aims` 的表头（head，2 行）：vendored 台账只认
  `mos.aims`，否则自旋极化 SCF 在 slurm-ssh 模式下本地拿不到 `nsaos`；闭壳层
  run 里两者缺席记为 `optional_missing`，不算缺件。自旋极化 run 反之：取到
  α/β 表头时 `mos.aims` 缺席同样记 `optional_missing`，不再判回取缺件
  （2026-11-19 TTF 实跑 collect 时发现）。全量 317 passed。

### 统一 T(E) 绘图版式（2026-11-19）

- 新增 `src/smj/plotting.py`：单一版式实现（四边框、刻度朝内、`E_F` 虚线、
  log 纵轴、固定配色 `PALETTE`、单栏/双栏/宽幅尺寸 `FIG_SIZES`、600 dpi
  PNG+PDF 同源导出、panel 角标、`T(E_F)`/`G/G0` 注释）。版式常量集中在
  `STYLE` / `STYLE_VERSION`，顶层不 import matplotlib（渲染时才加载，且强制 Agg）。
- 新增 `smj plot` 子命令：输入可为 `TE.dat`、run 目录或阶段目录，支持单图与
  多 run 叠加（固定配色、带框图例）、`--scale log|linear`、`--figsize
  single|double|wide`、`--formats`、`--dpi`、`--xlim/--ylim`、`--panel`、
  `--no-annotate`；缺省输出到 run 根的 `T_E.png`（多图 `T_E_overlay.png`）。
- `junction_analyze` 阶段在写 `results.json` 的同时按同一版式出
  `T_E.png` + `T_E.pdf`（matplotlib 未安装时只记 note、不阻断数值结果）。
- `smj analyze --plot` 与 `smj.analysis.plot_transmission` 均委托新实现，
  旧返回键（`path`/`points`/`log_scale`/`floor`）保持不变。
- 版式依据用户 Zotero 中 Yi 2016 *Chin. Phys. B* Fig. 6 与 Guan 2022
  *Nanotechnology* Figs. 2–3；说明见 `docs/plot_style.md`，新增
  `tests/test_plotting.py`（含 CLI 与阶段出图）；全量 296 passed。
- ET7H-BS2 结果图按新版式重绘（`runs/et7h-bs2/ET7H-BS2_TE.png` +
  同源 `.pdf`），删除该 run 下各画各的临时脚本 `plot_te.py`。

### T(E) 自旋分辨自适应绘图（2026-11-19）

- `read_curve` 改用 `landauer_spin.read_te_channels`，返回 `t`（双通道均值，
  与旧口径一致）、`t_alpha`/`t_beta` 与 `channels`；旧调用不受影响。
- `plot_te` 新增 `spin="auto"|"mean"|"channels"|"polarization"`：`auto` 在
  TE.dat 有 α/β 两列时自动出双 panel——(a) `Tα`/`Tβ`（固定配色
  `PALETTE[0]`/`PALETTE[1]`），(b) 自旋极化率 `P(E)=(Tα−Tβ)/(Tα+Tβ)`；
  单列数据、或 `channels`/`polarization` 遇到单列时退化为旧版式 `mean`。
  双 panel 在 `figsize="single"` 下自动升到 `FIG_SIZES["double"]`，panel 角标
  自动 `(a)`/`(b)`；`meta` 增加 `spin`/`channels`。
- `smj plot` 新增 `--spin`；`smj analyze --plot` 与 `junction_analyze` 仍走
  默认 `auto`，因此自旋极化 run 会自动得到自旋分辨图。
- `docs/plot_style.md` 补自旋双 panel 说明；`tests/test_plotting.py` 增加
  四列 TE.dat 用例（通道拆分、自动双 panel、`mean`/`channels` 降级、非法值、
  `smj plot --spin`），全量 325 passed。

### 运行范围策略（2026-11-18）

- 新增仓库级 `AGENTS.md`：只跑用户明确点名的分子/异构体/协议；不默认加跑
  pinned 对照、敏感性扫描、替代锚定、团簇尺寸扫描、空结对照等控制计算；需要
  对照时先说明理由和代价，等用户明确同意。
- `electrode.pin_contact: true` 现在在配置加载时给出警告，提示它是 pinned 对照/
  距离扫描口径、不是生产默认；`smj plan` / `smj run` / `smj preflight` 等命令都会
  打印。README §0 "不承诺"同步写明不自动加对照臂/敏感性计算。
- 新增 `tests/test_full_relax_protocol.py::test_pinned_contact_is_flagged_as_control`；
  全量 286 passed。

### 脚本整理（2026-11-18）

- 原仓库根的一次性脚本 `_azo.py` 整理为 `scripts/azo_ez_ratio.py`：两个 run 的
  `TE.dat` 路径与统计窗口改为参数（支持 run 根目录、`runs/<name>/<name>` 布局或
  `TE.dat` 路径），新增 `--json` 输出；数值口径与原脚本逐项核对一致。
  `.gitignore` 不再忽略 `_azo.py`。

### 结构归档 `saved_structures/`（2026-11-18）

- `molecule_opt` 在 prepare 时归档输入分子结构、collect 成功后归档 `opt.xyz` 到
  `<run>/saved_structures/`；`manifest.json` 记录 `stage` / `kind` / `format` /
  `engine` / `sha256` / `natoms` / `source` / `created_at`，collect 条目另写
  `<seq>_<stage>_<kind>.results.json` 结果快照，并记快照 sha256 与标量摘要
  （能量/收敛/最大力/弛豫步数），结构文件与源文件逐字节一致。
- 分子结结构链同样归档：`aims_tip_relax` → `molecule_tips`（分子 + 两枚 tip），
  `electrode_build` → `junction`（补全两侧电极后的完整分子结，collect 校验通过后），
  `au_relax` → `junction_relaxed`（贴装后弛豫的新几何；`electrode.relax=false` 的
  直通分支不重复归档）。FHI-aims 几何按原生 `geometry.in`（`.in`）落盘，`format`
  字段区分 `xyz`/`aims`，`natoms` 由解析器核对。
- 新增 `safe_archive()`：归档写盘失败只返回一行 note，不阻断 prepare/collect。
- 归档目录已加入 `.gitignore`；新增 `tests/test_structure_archive.py`（9 例）。

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
- 修 `junction/molecules.py::molecule_from_smiles` 的多构象 bug：
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
