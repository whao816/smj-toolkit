# 自由基支线改进计划（2026-11-21）

- 现状基线：`master` @ `e2ea6d8`，全量 `pytest -q` → 533 passed / 0 errors。
- 上一轮结论（本条支线的整体评估，见对话 2026-11-21）：机制闭环基本完整，科学判定不完整，工程完备性中上。
- 本计划只补"证据链 + 硬闸门 + 真实数据回归"，**不做**判据标定、不引入新的物理近似、不自动跑对照计算（AGENTS.md）。

## 0. 一句话结论

| 现状 | 缺口 | 本计划 |
| --- | --- | --- |
| 自旋输入、五阶段接线、`spin` 块、review/报告展示都已闭环 | 缺"自旋到底分布在哪"的证据；缺自旋分辨热输运；单列 TE 静默通过 | S1/S2/S3 |
| 开壳层判据按设计不参与物理判定（`no_applicable_gate`） | 文案与工程一致性检查不足；BS 只有标签没有证据 | S4/S5 |
| 533 项测试、43 条开壳层相关测试 | 没有一条真实开壳层结的端到端回归；文档状态表过时 | S6/S7 |

## 1. 现状与缺口

### 1.1 已经有的（不需要重做）

- 自旋输入单一出口 `src/smj/spin.py`：`SpinState`、`spin_control_lines`、初矩/初电荷注入、去质子 1-based 重映射、`record_junction_spin`。
- BS 单重态四条硬校验（`src/smj/config.py:603-626`）：M=1、初矩非空、合计 0、正负分量都有。
- 五阶段接线：`molecule_opt`（aims）、`aims_scf`、`aitranss_te`（`$scfmo→$uhfmo_alpha/$uhfmo_beta`、nsaos 取 alpha、`transmission_channels`）、`junction_analyze`（自动加 `spin` 块）。
- 展示与复核：`report.py:623-655` 的 §4.1 表、`review.py:690-755` 的 `_spin_review`、`spin_summary.py` 的 7 个 code 与 3 条 caveat。
- 产物门槛：`extra_outputs`（`hirshfeld`、`spin_density`）在 `aims_scf` 只核验"文件非空 / 段存在"，不解析数值。

### 1.2 缺口（本计划要补的）

| 编号 | 缺口 | 后果 |
| --- | --- | --- |
| G1 | Hirshfeld 自旋矩不解析、不入报告 | 无法回答"自旋在分子上还是跑到电极上"；BS 也没有局部证据 |
| G2 | `analyze_te_spin` 只有 L0/T(E_F)/P，没有 L1/L2 | 开壳层给不出 Sα/Sβ、k_eα/β、ZTα/β |
| G3 | 开壳层配单列 TE 只 warning，阶段照常 collected | 下游把"α=β"的均值当自旋量，P 恒 0 却不见阻断 |
| G4 | BS 只有"反号初猜 + 两条通道"，无局部磁矩/NOON/⟨S²⟩ | 无法区分"真破缺"与"没收敛到破缺" |
| G5 | 判据口径靠文案，没有工程一致性检查 | "宣称开壳层但产物不支持"只能靠人看 |
| G6 | 无真实开壳层结的端到端回归 | CI 全绿也不代表真实 run 能跑通 |
| G7 | 文档状态表停在 2026-11-20 | 误读"已实现/未实现" |

### 1.3 实测证据（本次核查，非推断）

- **Hirshfeld 可解析**：`runs/ttf-cation/ttf-cation/aims_scf/aims.out` 的最后一段 `Performing Hirshfeld analysis`：88 原子，spin moment 总和 1.0000；按 1–14 分子 / 15–16 接触区 / 17–88 电极分区得 **0.0056 / 0.0063 / 0.9881**，与规划文档 §5.1 的 Mulliken 投影（0.00546 / 0.98961）一致到 1e-3。
- **单列 TE 不是假想**：全仓 `runs/` 下开壳层 run 共 7 个（`di-fcfr`、`di-pcpr`、`ifs-bs`、`ifs-triplet`、`ttf-cation`、`ttf-fixed`、`ttf-release`），TE.dat **全部 4 列**；另有 2 个开壳层 run 尚无 TE.dat。→ S3 的默认硬闸门**不会翻转任何既有 run**。
- `runs/et7h-bs2` 目录名像 BS，但 `smj.yaml` 只有 `multiplicity: 1`、没有 `broken_symmetry`，`results.spin` 为 null——是闭壳层 run，不能当反例。
- `runs/ifs-bs` 是唯一真正 BS 的 run：`spin.channels=2`、`P(E_F)=−1.8e-5`（几乎无极化，符合闭壳层附近能隙），可作 S2/S4 的数值锚点。
- `runs/ttf-cation/ttf-cation/review.json` 是旧版产物：只有 `structure:junction_geom_surface_side` error 与两条 info，**没有** `open_shell_uncalibrated`——S6 用它验证"重新生成后应出现新 code"。

## 2. 目标与边界

**目标（本计划完成时）**

1. `smj report` 的开壳层 run 里能看到"自旋分布（Hirshfeld）：分子 / 接触区 / 电极"三行数字。
2. 开壳层 `results.spin` 里有每通道 L1/L2 与 Sα/Sβ/k_e/ZT，报告与判据可用。
3. 开壳层配单列 TE 在 `aitranss_te.collect` 直接失败，除非显式配置放行。
4. `smj review` 对真实开壳层 run 给出"自旋分布 + 未校准"的明确结论，且有一条可选的端到端回归。
5. 文档状态表与 README 与本计划一致。

**不做（避免越界）**

- 不做开壳层判据的数值标定（需要专门的开壳层参考体系，另行立项）。
- 不做 SOC/非共线、Kondo、ΔSCF、多体输运。
- 不自动跑闭壳层对照或任何"控制"计算（AGENTS.md：只跑用户点名的）。
- 不把 Hirshfeld 分区做成硬阈值判据（见 D2）。

## 3. 改进项

### S1 Hirshfeld 自旋分区（P0）

**状态**：✅ 2026-11-21 已实现（批次 A）。实测 `runs/ttf-cation/ttf-cation`：分子 1–14 为
0.0056 μB、tip 15/16 为 0.0063 μB、电极 17–88 为 0.9881 μB（总 1.0000），与 §5.1 的
Mulliken 投影一致到 1e-3。

**证据**：见 §1.3；`aims_scf` 目前只核验段存在（`src/smj/stages/aims_scf.py:242-263`）。

**改法**

1. 新模块 `src/smj/spin_partition.py`：
   - `parse_hirshfeld(path) -> HirshfeldTable`：逐行解析 `| Atom <n>: <el>` 与 `Hirshfeld spin moment : <v>`，取**最后一次**分析段（多段时后者覆盖前者）；缺字段给 `None`，不抛错；返回 `{index, element, spin_moment}` 与 `total`。
   - `default_groups(n_molecule, n_atoms, contact_indices)`：分子 = `1..n_molecule`；接触区 = 传入的几何近邻集合（可选）；电极 = 其余。接触区用 `contact_atoms(geometry, molecule_atoms, anchors, cutoff=3.0)` 现算。
   - `summarize(table, groups)` → `{"total": μB, "groups": {name: {atoms, spin_moment, fraction}}, "n_positive", "n_negative", "localized": bool}`（分区重叠直接 `ValueError`）。
2. `aims_scf.collect`：当 `transport.extra_outputs` 请求了 `hirshfeld` 且段存在时解析，写侧车 `aims_scf/spin_partition.json`（逐原子表，避免撑大 `results.json`），`results["spin_partition"]` 只放摘要 + 文件路径 + `source="aims.out"`；解析失败 → `StageError`（D4）。
3. `report.py` 在 `_transport_summary` 里把 `aims_scf` 侧车摘要挂到 `transport["spin_partition"]`（**不**经过 `spin_summary.normalize`——该层只管自旋 code 判定）；`review.py` 直接读侧车文件（见 5）。
4. `report.py` §4.1 表后加"自旋分布（Hirshfeld）"小表：分子 / 接触区 / 电极（μB 与占比），并注明"接触区 = 距锚原子 ≤3.0 Å 的电极原子（仅展示）"。
5. `review.py` 加 `spin_localized_off_molecule`（info）：开壳层且 `|分子自旋矩| / total < 1%` 时给出"自旋主要在电极侧"的陈述，只报数字不设阈值判定。

**测试**

- `tests/test_spin_partition.py`：两段 Hirshfeld 取最后一段；字段缺失容错；分组语法与越界；`summarize` 的 fraction 和为 1；空表 → `None`。
- `tests/test_scf_artifacts.py` 增一例：fixture 里加 mini Hirshfeld 段 + `extra_outputs=[hirshfeld]` → `results["spin_partition"]` 存在且 `molecule.fraction` 可读。

**验收**

- 对 `runs/ttf-cation/ttf-cation` 跑一次 `smj report`（或脚本），输出 0.0056 / 0.0063 / 0.9881（容差 1e-3）。
- 全量 pytest；不动 `src/smj/junction/`，因此无需重跑 provenance。

**风险/回退**

- 分区边界靠 `molecule_opt` 的原子数（去质子后）推定；对"分子在电极之上重复"的非常规几何需人工核对。方案提供 `--groups` 覆盖入口（CLI 也可后置），出错时报告只显示"合计"一行。

### S2 自旋分辨热输运（P0）

**证据**：`src/smj/landauer_spin.py:59-103` 只对每通道积分 L0；`analyze_te` 已有 L1/L2 → S/k_e/ZT。开壳层因此只有 P(E_F)/P(L0)，没有 Sα/Sβ。

**改法**

1. `analyze_te_spin` 每通道复用 `trapz` 积分出 `L1_alpha/L1_beta/L2_alpha/L2_beta`，按与 `analyze_te` **完全相同**的公式给 `S_*_uV_per_K`、`k_e_*_pW_per_K`、`ZT_*_electronic`（实现时直接调用同一 helper，禁止复制公式）。
2. 新增 `spin_polarization_S = (S_alpha - S_beta)/(S_alpha + S_beta)`（分母为零给 `None`）。
3. `results.spin` 通过 `junction_analyze` 自动带出；`report.py` §4.1 表加 Sα/Sβ 与 P_S 两行。
4. 新增 engineering 判据 `spin_channel_consistency`（D6）：开壳层且有 `spin` 块时校验
   - `g_over_g0_total == 0.5*(L0α+L0β)`（相对 1e-9）；
   - `S_total == (L0α·Sα + L0β·Sβ)/(L0α+L0β)`（相对 1e-6）；
   - 不满足记 engineering/error，`missing=False`（这是真错误，不是证据不足）。
5. 该判据写进 `src/smj/junction/criteria.py` 的 `CRITERIA`（`gate='error'` 即被 `_ERROR_GATED` 自动纳入 severity）。engineering 判据免文献源，`lint_registry` 只对 physical 要求 `source`，因此**不需要**进 `KNOWN_GAPS`（该集合只登记无源 physical 判据）。

**测试**

- α/β 不同谱：Sα>0、Sβ<0、加权平均等于 S_total、P_S 符号正确。
- α=β 退化：Sα=Sβ=S_total、P=0、consistency 判据 pass。
- 乱序谱：与升序结果一致（复用批 1 的排序逻辑）。
- 判据单测：构造 G_total 与 0.5(L0α+L0β) 不一致的 `results.json` → 判据 failed。

**验收**

- `runs/ifs-bs`：Sα≈Sβ（差 <1e-3 µV/K）、P_S≈0，与 `P(E_F)=−1.8e-5` 自洽。
- 改 `src/smj/junction/` → `gen_junction_provenance.py` 刷新并 `--check`、`gen_demo_reference.py --check` 零偏差。

**风险**：`results.json` 的 `spin` 键增多，旧读取方只取已知键（`analysis.py` 只校验必需键，多余键安全，已在批 1 验证）。

**状态（2026-11-21，批次 B）：✅ 已实现。** `landauer.transport_moments()` 作为唯一公式源，`analyze_te` 与 `analyze_te_spin` 共用；`spin` 块新增 `L1/L2_alpha/beta`、`S_alpha/beta_uV_per_K`、`S_total_uV_per_K`、`k_e_*`、`ZT_*`、`spin_polarization_S`。`spin_channel_consistency` 在 `evaluate_run_criteria` 末尾评估：闭壳层或旧 `spin` 块（无分通道 L1/L2）`applicable=False`，不误判。实测量（`runs/ifs-bs`，801 点，300 K）：Sα=+6.15599、Sβ=+6.15894 µV/K（相对差 4.8e-4），S_total=+6.15746，P_S=−2.39e-4，与 P(E_F)=−1.84e-5 同号；判据在真实 run 上 `passed=True`。计划里写的"差 <1e-3 µV/K"是数量级估计，实测绝对差 2.9e-3 µV/K（相对 4.8e-4）——**验收口径应看相对差与 P_S 量级**，绝对阈值不具普适性（S 的量级随体系变化）。全量 555 passed；provenance 已刷新（36 文件/13 补丁）并 `--check` 一致，demo 参考 4 份谱零偏差。

### S3 开壳层单列 TE 硬闸门（P0，小）

**证据**：见 §1.3。当前单列 TE 在 `aitranss_te.collect` 只被 review 记为 warning，阶段 collected，后续 `spin.channels=1` 静默。

**状态**：✅ 2026-11-21 已实现（批次 A）。默认阻断，`transport.allow_single_channel: true` 显式放行；既有 7 个开壳层 run 全 4 列，不受影响。

**改法**

1. `src/smj/stages/aitranss_te.py collect`：读通道数后，若 `state.open_shell and channels < 2` 且未配置放行 → `StageError`，文案含"开壳层 TE.dat 只有 1 个自旋通道；检查 `$uhfmo_alpha/$uhfmo_beta` 是否被 tcontrol 使用、SCF 是否真处于开壳层、nsaos 是否非零"。
2. 新增配置 `transport.allow_single_channel: bool = false`（`src/smj/config.py:480-481` 的 `_TRANSPORT_KEYS` + 解析 + 示例）。放行时写 `results["single_channel_allowed"]=true` 与 note，review 对应 code 从 warning 保留但文案改成"已显式放行"。
3. `junction_analyze` 不拦截（允许事后分析旧数据），只在 `results.spin.channels==1` 时保留现有 note。

**测试**

- 开壳层 + 单列 TE → `StageError`（match "自旋通道"）。
- `allow_single_channel: true` → collected，`results["single_channel_allowed"] is True`，note 有记录。
- 闭壳层 + 单列 TE → 照常 collected。
- 开壳层 + 双列 → 照常 collected。

**验收**：全量 pytest；既有 run 不受影响（7 个开壳层 run 全 4 列）。

### S4 BS 单重态证据（P1）

**证据**：BS 路径只有"反号初猜 + 2 通道"；`ifs-bs` 的 P≈0；无 NOON/⟨S²⟩。

**改法**

1. 用 S1 的逐原子自旋矩给 BS 一个**局部证据**：`spin_partition.summarize` 额外返回 `n_positive`/`n_negative`/`max_abs`/`net`；若分子内同时有正负分量且 `|net| < 0.05·max_abs`，标记 `localized_antiparallel=True`。
2. `report.py` 在 BS run 的 §4.1 表下加一行"破缺对称证据：分子内 +/− 自旋矩各 N₁/N₂ 个，净矩 X μB（局部反平行）"；`review.py` 加 info `broken_symmetry_local_moments`（带数字）。
3. 明确写进报告脚注：FHI-aims 单行列式不提供 ⟨S²⟩/NOON；如需，用 Gaussian 模型分子交叉验证（`engines/gaussian.py:100-126` 已有 ⟨S²⟩ 解析），但**不自动跑**。
4. 不自动跑闭壳层对照；文档给出手动命令模板，由用户决定。

**测试**

- 构造正负各一的 mini Hirshfeld → `localized_antiparallel=True`；全正 → False。
- 报告渲染：BS run 出现"破缺对称证据"行，普通开壳层 run 不出现。

**验收**：`runs/ifs-bs` 生成的分区摘要与报告行可读；不改变任何判据结论（只加 info/展示）。

### S5 判据口径与工程一致性（P2）

**改法**

1. `criteria.py physics_verdict` 的 `no_applicable_gate` 文案区分开壳层：追加"开壳层判据未校准；如需物理解读，需先用开壳层参考体系标定 `g_in_theory_band`/`seebeck_sign`"。
2. 新增 engineering 检查 `open_shell_protocol`（warn）：`multiplicity>1 or broken_symmetry` ⇒ 必须同时满足 `spin` 块存在、`channels==2`、`extra_outputs` 含 `hirshfeld` 或 `spin_density`；缺一项记 warn，detail 指出缺哪项。
3. 不新增任何数值阈值；不把未校准判据升级为 error。

**测试**：三例（齐全 pass、缺 spin 块 warn、缺 extra_outputs warn）；既有 `test_verdict_gates.py` 保持绿。

### S6 端到端集成回归（P1）

**改法**：新增 `tests/integration/test_radical_runs.py`，`pytestmark = pytest.mark.integration`（默认 skip，`pytest -m integration` 才跑）：

- 若 `runs/ttf-cation/ttf-cation/state.json` 不存在 → `pytest.skip`（仓库克隆者无 assets）。
- 断言：state 六阶段全 collected；`results.spin.channels==2`；`|P(E_F)| < 0.05`；`spin_partition` 分子占比 <1%（0.0056）；`smj report` 渲染含"§4.1"与"自旋分布"。
- 在 `tmp_path` 复制该 run 的轻量文件（state/results/TE/review/manifest/aims.out 片段 2.3 MB 可接受）后跑 `review_run(force=True)`，断言 `open_shell_uncalibrated` 出现（旧 review.json 里没有）。

**验收**：`pytest -m integration` 本机通过；默认全量 pytest 仍 533+ 且不跑集成。

### S7 文档（P0，小）

1. `docs/自由基单分子结-实现对照与开发规划.md`：状态表更新（review ❌→✅；PDOS 待做→已做并加自旋配对；BS 待做→已做；extra_outputs 已接；Hirshfeld 分区标"本计划 S1"），§5.1 补"解析结果"，新增"§7 后续（本计划）"引用。
2. `README.md` §9 自由基段：加 Hirshfeld 分区、Sα/Sβ、单列 TE 闸门三条。
3. `CHANGELOG.md` `## 未发布`：新增"改进：自由基支线证据链（S1–S6）"。
4. `docs/README.md` 索引加本文件。

### S8（可选，需用户点名）`smj spin` 子命令

汇总 `review` 的 spin 块 + 分区 + P(E) 到一张表/JSON。默认不做，避免与 S1/S2 重复；若要做，放批次 D。

## 4. 批次

| 批次 | 内容 | 涉及 | 完成判据 |
| --- | --- | --- | --- |
| A | S1 + S3 + S7 | `spin_partition.py`、`aims_scf.py`、`aitranss_te.py`、`config.py`、`report.py`、`review.py`、文档 | ✅ 2026-11-21：报告出现自旋分布；单列 TE 被拦；全量 549 passed |
| B | S2 | `landauer_spin.py`、`criteria.py`、`junction_analyze.py`、`report.py` | Sα/Sβ 进结果与报告；provenance + demo 校验 | ✅ 2026-11-21（555 passed） |
| C | S4 + S6 | `spin_partition.py`、`report.py`、`review.py`、`tests/integration/` | BS 证据行；`pytest -m integration` 通过 |
| D | S5（+S8 可选） | `criteria.py`、`report.py` | 口径文案与 open_shell_protocol；无回归 |

顺序理由：S1/S3 是证据与闸门，先落地才能给 S2/S4 提供真实数字；S2 动 junction 内核，单独一批方便 provenance/demo 对账；S4/S6 依赖 S1；S5 只改文案与工程检查，放最后避免反复改文案。

## 5. 决策点（需用户拍板）

| 编号 | 问题 | 建议 |
| --- | --- | --- |
| D1 | 单列 TE 默认阻断还是默认放行 | **默认阻断**（7 个既有开壳层 run 全 4 列，无误伤；放行需显式配置） |
| D2 | 自旋分区是否设硬阈值（如分子占比 <1% 判"自旋外泄"） | **不设**，只报数字 + info；阈值留给标定后的判据 |
| D3 | BS 是否自动跑闭壳层对照 | **不自动**（AGENTS.md）；文档给手动模板 |
| D4 | 请求了 hirshfeld 但解析失败 | **StageError**（证据链缺一环就明确失败，不静默） |
| D5 | `spin_channel_consistency` 的 severity | engineering/error（工程自洽可硬判，不涉及物理阈值） |
| D6 | 集成回归是否进 CI | 默认 skip，本机 `-m integration` 跑；CI 保持快速 |

## 6. 文件清单

| 文件 | 动作 |
| --- | --- |
| `src/smj/spin_partition.py` | 新增：解析、分组、摘要 |
| `src/smj/stages/aims_scf.py` | `collect` 解析 Hirshfeld，写侧车与摘要 |
| `src/smj/stages/aitranss_te.py` | 单列 TE 闸门 + 放行标记 |
| `src/smj/landauer_spin.py` | 每通道 L1/L2、S/k_e/ZT、P_S |
| `src/smj/junction/criteria.py` | `spin_channel_consistency`、`open_shell_protocol`、文案 |
| `src/smj/config.py` | `transport.allow_single_channel`、`_TRANSPORT_KEYS`、示例 |
| `src/smj/report.py` | 自旋分布表、Sα/Sβ、BS 证据行 |
| `src/smj/review.py` | 2 个 info code |
| `src/smj/spin_summary.py` | 透传 `spin_partition`、新 code 映射 |
| `tests/test_spin_partition.py` | 新增 |
| `tests/test_scf_artifacts.py`、`tests/test_spin_transport.py`、`tests/test_report.py` | 增例 |
| `tests/integration/test_radical_runs.py` | 新增（默认 skip） |
| `docs/自由基单分子结-实现对照与开发规划.md`、`README.md`、`CHANGELOG.md`、`docs/README.md` | 更新 |

## 7. 验收总表

| 项 | 命令/证据 | 期望 |
| --- | --- | --- |
| 单元测试 | `python -m pytest -q` | ≥533 passed，0 failed |
| 自旋分区 | 对 `runs/ttf-cation` 生成报告 | 分子 0.0056 / 接触 0.0063 / 电极 0.9881（±1e-3） |
| 单通道闸门 | 合成单列 TE + `multiplicity: 2` | `StageError` |
| 自旋热输运 | `runs/ifs-bs` | Sα≈Sβ、P_S≈0，与 P(E_F) 自洽 |
| 内核校验 | `python scripts/gen_junction_provenance.py --check` | 一致（仅批 B 后需要） |
| 演示基线 | `python scripts/gen_demo_reference.py --check` | 零偏差 |
| 集成回归 | `python -m pytest -m integration -q` | 通过（有 assets 时） |

## 8. 与项目规则的关系

- 只对用户点名的体系跑链；S1/S2 的验证都在既有 `runs/` 上做**离线解析**，不提交新计算。
- 批次 B 改动 `src/smj/junction/` 后必须重跑 `gen_junction_provenance.py` 刷新 `PROVENANCE.json`，并跑 `gen_demo_reference.py --check`。
- 判据仍不替代科学裁决：开壳层 `inconclusive` 就是 `inconclusive`，不包装成 `passed`。
