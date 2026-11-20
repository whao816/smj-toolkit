# 对《smj-toolkit 评审报告》的复核（2026-11-20）

本文复核 `D:/codex_sessions/2026-11-19/https-github-com-whao816-smj-toolkit/outputs/smj-toolkit-review.md`
提出的 7 条问题、3 条后续完善建议与「验证范围」声明，逐条给出**成立 / 部分成立 / 不成立 / 无法复核**的判定，
并记录复核中发现的、原报告未提的新问题。

**不改原报告、不改仓库代码**：本次只读代码、跑测试与复现脚本，结论落在本文。

## 0. 方法与口径

- 被复核对象：提交 `0f7b3c4a3ce3106029c168d17f47ae1e29772f9b`（release 0.2.1rc1）。
- 复核环境：本仓库 `D:\dsh_projects\AutoSMJ`，HEAD `14bc480`（比被评审提交多 5 个提交），**工作区含未提交改动**
  （`aims_scf.py`、`aims_tip_relax.py`、`provenance.py`、`config.py` 等）。因此每条都要区分
  「被评审提交的行为」与「当前工作区的行为」。
- 旧提交代码的取得方式：`git archive 0f7b3c4 tests src pyproject.toml` 导出到 `.tmp/verify/scope/tree-0f7b3c4/`，
  **不 checkout、不切换分支、不碰工作区**。行号引用以该副本为准，与当前工作区不同处单独标注。
- 四类证据：
  1. **只读代码对照**：把原报告引用的 12 处 `file:line` 逐条打印到被评审提交的源码上；
  2. **数值复现**：用仓库自带真实谱（薁 190 点）与合成谱复算 `L0/L1/L2/S/k_e/ZT`；
  3. **构造输入**：合成「未收敛日志 + 已写 `geometry.in.next_step`」、零字节产物、只放 TE.dat 的空 run、
     `waiting`2 状态的重跑 —— 均为 `.tmp/verify/*/work/` 下的一次性目录，未进仓库；
  4. **测试与 CLI 实测**：`pytest`、`smj verify --json`、`smj report --json`。
- 结论口径：**成立**（现象与数字均可复现）、**部分成立**（机制真实但原报告的表述/范围/前提需修正）、
  **不成立**（与实测冲突）、**无法复核**（无产物可审计）。另标「说轻」「说重」两种偏差方向。
- 严重度只按「后果 × 可达性」排，不按代码行数。

## 1. 总览

| # | 原报告主张 | 复核结论 | 与原报告的差异 |
|---|---|---|---|
| P1-1 | 倒转谱序改变 Seebeck 符号 | **成立**（逐位复现） | 无；补一句：L0/L2/k_e/ZT 不变，**只有 S 的符号被翻转** |
| P1-2 | 能量窗完全偏离 E_F 仍 `physics_verdict=passed` | **成立**（逐位复现） | 补：该 fail-open **可由通过校验的配置直接触发**（原报告未证明可达性）；行号引用范围偏窄 |
| P1-3 | 针尖弛豫缺收敛验收，未收敛几何流入建电极 | **成立** | **说重**：「日志明写未收敛」不能用随仓库发布的 demo 举证（四份 demo 日志最终都收敛）；真实危险路径是「日志未收敛但写出了 `next_step`」 |
| P1-4 | 静态 SCF 可采集空矩阵 | **成立** | 补：`validate` 能拦未收敛标记、`collect` 不拦，两级判定不一致；**零字节产物两级都不拦**；`results` 里的 `has_*` 是常量不是观测 |
| P1-5 | 重跑会把活跃作业改标 failed，可能重复提交 | **成立** | **说轻**：不是「可能」，第二次 `run` **必然**重新 prepare + sbatch 覆盖作业号 |
| P1-6 | 仅有透射谱的目录仍 `engineering_passed=true`、`errors=0` | **成立** | **说轻**：顶层 `verdict.passed` **也为 true**；另有 `warnings=2` 被原报告省略；「报告缺完成度」**说重**（报告已印阶段数/未完成阶段） |
| P2-1 | 改步长协议指纹不变 | **部分成立** | 仅在 **eV 相对窗 + `estep/ener/eend` 覆盖键**下成立；Hartree 窗改 `estep` 哈希会变；原报告把「默认 eV 窗步长 0.000183747 Ha」当成 demo 实际步长（demo 是 0.0002 Ha） |
| 后续① | demo 谱做成发布回归样本 | **部分成立** | 方向对；但「哈希、温度全无」不成立：温度 4/4 已记录，3/4 已有与实际一致的 TE.dat sha256；真缺口是**参考值、容差、统一清单、CI 重算** |
| 后续② | 增加结构与机理证据 | **部分成立** | 结构归档**已实现并在跑**（`structure_archive.py` + 每阶段 `saved_structures/manifest.json`）；缺口是机理证据字段与跨 run 对照，不是「没有结构证据」 |
| 后续③ | 统一 VENDOR 维护说明 | **成立** | 补：`scripts/gen_junction_provenance.py --check` **是存在的可执行核验**，当前 exit=1 只因 `repo_commit` 绑死旧 HEAD（36 文件全一致）；原报告说「本仓库没有该脚本」属实，但不能推出「无法自动核验」 |
| 验证范围 | Python 3.12.14；267 项全通过 18.40 s；`smj verify` 六项全过 | **部分成立** | 3.12.14 **不成立**（本机 venv 3.12.10，无 3.12.14）；267 项与 `0f7b3c4` 吻合（实测 19.08–20.06 s），但**未绑定提交**（当前工作区 330 项）；`smj verify` 六项**逐字成立**；「装包后全通过」未说明 matplotlib 只是 extra、rdkit 完全未声明、CI 实际 248 passed + 3 模块级 skip |

**行号引用质量（对被评审提交）**：12 处引用**全部落在其声称的代码上**，无错位；唯一瑕疵是 P1-2 的
`criteria.py:412-443` 只覆盖到「无显著峰 → 记 pass」，真正放行的两处（`seebeck_sign` 的 `applicable=False`
在 :476-482、`physics_verdict` 归约在 :273-323）在引用范围之外。

## 2. P1-1 倒转谱序 → Seebeck 反号（成立）

**原报告**：薁示例 190 点，升/降序 `L0=0.74614719` 相同，`S=+25.93621 / -25.93621 µV/K`。

**复核**：`.tmp/verify/lead/repro_claims_1_2.py` 读真实谱
`demo/azulene-26-dithiol/azulene-26-dithiol/aitranss_te/TE.dat`（10400 B，194 行 = 4 行注释 + 190 数据行），
分别按文件原序（升序）与反转（降序）复核，另做一次「先按 E−E_F 排序再积分」的顺序无关参考：

| 序列 | L0 | L1 | L2 | S (µV/K) | k_e (pW/K) | ZT |
|---|---|---|---|---|---|---|
| 升序（原文件） | 0.74614719 | −5.80566880e−03 | 1.54749514e−03 | **+25.93621** | 388.004467 | 0.030069 |
| 降序（行序反转） | 0.74614719 | +5.80566880e−03 | 1.54749514e−03 | **−25.93621** | 388.004467 | 0.030069 |
| 排序后积分（物理正确） | 0.7461471866885822 | −0.005805668801406427 | 0.0015474951380115456 | **+25.936208945458937** | — | — |

- `L1_desc / L1_asc = -1.000`，`S_asc + S_desc = 0.000e+00`，`L0_equal = True`。
- 机制：`landauer.py:70` 与 `:72` 对偶数阶矩取 `abs(...)`，`:71` 的 `l1` 保留符号；降序使梯形积分整体变号，
  绝对值把 L0/L2 的变号吸收掉，于是**只有 S 反号**，其余派生量照旧。
- 排序后积分与升序结果在 8 位有效数字上一致 → **升序才是物理正确值**，降序是静默错误。
- 全仓 TE→积分路径**无排序/单调性校验**：`landauer.read_te`、`landauer_spin.read_te_channels` 都按行序 append；
  调用方 `analysis.py:266`、`stages/junction_analyze.py:43`、`stages/aitranss_te.py:148`（只数点数）均不检查顺序。

**结论：成立，数字与原文一致。** 原文未写、但值得补的是：该缺陷**只翻转 S 的符号**，不污染 G/κe/ZT，
所以它能穿过所有「量级在带内」类判据，只在 `seebeck_sign` 上现形 —— 而 `seebeck_sign` 对非 BDT 族是
`applicable=False`（见 P1-2），两条缺陷因此**互相掩护**。

## 3. P1-2 窗口偏离 E_F 仍物理通过（成立）

**复核做法**：构造恒定 T=0.25 的合成谱，200 点，`E−E_F ∈ [+0.250, +0.449] eV`（整窗在 E_F 之上），
`temperature=300`，`system="AZ"`（族 `Au-(2,6-azulenedithiol)-Au`），走
`evaluate_run_criteria()` + `physics_verdict()`：

| 量 | 实测 |
|---|---|
| L0 / L1 | 1.5774441100123155e−05 / 4.349896753482421e−06 |
| S | **−919.1866610620011 µV/K**（原报告 −919.19 ✓） |
| k_e / ZT | 0.002650266303394581 / 116.89288932777318 |
| `physics_verdict` | **passed**，`wording="支持 H0"`，`allowed_action="可进入发布(A5)/预测卡(A4)"` |
| 未通过的判据 | `g_in_theory_band`（warn）、`seebeck_sign`（error）—— **两条都 `applicable=False`** |
| 通过的判据 | `t_ef_resonance_attribution`（"E_F±1.0 eV 内无显著峰"）、`k_e_zt_order` |
| `not_applicable_criteria` | `["g_in_theory_band", "seebeck_sign"]` |

**机制**：`g_in_theory_band` 与 `seebeck_sign` 的注册项带 `applies_to={'families': ('Au-BDT-Au',)}`，
对 AZ 记 `applicable=False`；`physics_verdict` 在 `criteria.py:292-294` 直接 `continue`，不计闸门；
剩下的 warn 级判据「无显著峰 → pass」反而被记为通过，于是 fail-open 闭环。

**原报告未证明的一点（本次补上）**：该场景**不是**只能靠手搓 API 触发，而是可由通过 `ConfigError` 校验的配置
直接产生。`config.py:641 _parse_transport` 只校验 `estep > 0` 与 `eend > ener`（:649-650），**不要求窗口含 E_F**；
实测 `energy_window=[0.25, 0.000499, 0.449]` + `energy_window_unit="eV"` 被接受、零 warning，
`resolved_window(-0.188452434085)` 得 `(-0.17926510354108627, 1.833791176565184e-05, -0.17195198842813092)` Ha，
整窗位于 E_F 之上 —— 这正是原报告构造的场景，且走的是**正常配置路径**。

**结论：成立。** 偏差：引用范围 `criteria.py:412-443` 没有覆盖真正放行的 `seebeck_sign` 适用性（:476-482）
与 `physics_verdict` 归约（:273-323），读者按引用逐个看代码会找不到「为什么 passed」。

## 4. P1-3 针尖弛豫缺收敛验收（成立，但举证方式说重）

**复核做法**：`aims_tip_relax` 的 `collect` 走 `StageRunner`，四组输入 + 一组 0f7b3c4 原文件对照。

| 用例 | 输入 | 结果 |
|---|---|---|
| 合成「未收敛」 | `aims.out` = 全程 `not yet converged` + 步数上限 + `Have a nice day`；**有** `geometry.in.next_step` | **`collected`**，`source="geometry.in.next_step"`，`relaxed_au_s_A = 2.4000000018`（左右） |
| 无日志 | 删掉 `aims.out`，保留 `next_step` | **`collected`**，同 payload |
| 无 `next_step` | 删掉 `next_step` | `collected`，`source="geometry.in"`，note「…按初始 geometry.in 记录」 |
| 读文件跟踪 | 记录 `collect` 期间的 `open()` | 只读 `geometry.in.next_step` 与 `input_summary.json`，**`aims_out_read = false`** |
| 被评审提交原文件 | 同上未收敛输入 | **`collected`**，payload **与当前实现完全相同**（`same_payload_as_current: true`） |

**流入下游**：把上述未收敛 `next_step` 交给 `electrode_build`，`prepare` 成功
（`parent_fingerprint_equals_tip_next_step: true`），`collect` 得 `final_atom_count=26`、4 Au/电极 →
**未收敛的针尖几何确实被当作「弛豫后」结构用于建电极**，原报告此处判断正确。

**原报告说重的地方**：「日志明确写着几何尚未收敛」**无法用随仓库发布的样本举证**。四份 `demo/` 下的
`aims_tip_relax/aims.out` 最终都收敛：

| demo | 步数 | `not yet converged` 出现次数 | `Present geometry is converged.` | 结尾 |
|---|---|---|---|---|
| azulene-26-dithiol | 54 | 54（:12035…:213869，中间每步一行） | 1 次，:216550 | `Have a nice day.` :216778 |
| 2CN-656-4-thiomethylphenyl | 87 | 87 | 1 次，:395489 | 同 |
| 2CN …-aimsopt | 75 | 75 | 1 次 | 同 |
| 2CN …-fullrelax | 73 | 73 | 1 次 | 同 |

即：`not yet converged` 是**中间过程行**，最后一行是 converged。（日志里的 `maximum number of relaxation …`
（azulene :740 / 2CN :921）是 `control.in` 的步数上限设置行，不是失败标记 —— 按关键词检索时容易误判。）
正确表述是：**只要 `geometry.in.next_step` 存在就采集，与日志是否收敛无关**（本复核用合成日志证明），
而真实世界里 BFGS 走到步数上限时**同样会写出 `next_step`** —— 这才是危险路径。

**补充（原报告未提）**：`no_next_step` 回落分支不会静默流到下游 —— `collect` 登记的 artifacts 是
`["geometry.in"]`，而 `electrode_build.prepare` 按文件名 `geometry.in.next_step` 找上一阶段的**同名文件**
（`context.py:102 upstream` → `require_stage_file` → `layout.artifact`，路径式而非登记表），因此回落分支在下游会
`StageBlocked`。由此推论：修复的验收测试应构造「**有 `next_step` + 日志未收敛**」，构造「缺 `next_step`」测不到缺陷。

## 5. P1-4 静态 SCF 采集空矩阵（成立）

**复核做法**：6 组合成输入直接调 `aims_scf.collect`。

| 用例 | `aims.out` | 三个产物 | 结果 |
|---|---|---|---|
| A 零字节 | 330 B，含未收敛信息 + `Have a nice day` | `basis-indices.out`/`omat.aims`/`mos.aims` **各 0 B** | **`collected`**，`has_basis_indices/has_omat/has_mos = true` |
| B 有内容 | 同上日志 | 各 49–52 B | `collected`，同 payload |
| C 零字节 + 干净日志 | 226 B，无失败标记 | 各 0 B | **`collected`**，`has_* = true` |
| D 无退出语 | 无 `Have a nice day` | 正常 | `StageError: FHI-aims 未正常结束…` |
| E 空日志 | 5 B | 正常 | `StageError: aims_scf/aims.out 为空…` |
| F 缺 mos | 正常 | 无 `mos.aims` | `StageError: …SCF 产物不完整…` |

**关键点**：`results.aims_scf` 里 `normal_termination / has_basis_indices / has_omat / has_mos` 是**写死的常量**
（当前 `aims_scf.py:183-185`），不是对文件的观测 —— 零字节时它们照样为 `true`。

**两级判定不一致（原报告未提）**：同一 run 交给报告期 `validate`（`junction/validate.py:96-101`）：

| 用例 | `validate` 的 `aims_scf` 检查 | 结论 |
|---|---|---|
| A（未收敛日志 + 零字节） | `passed=false`，error「FHI-aims 输出含失败标记(SCF 未收敛/异常): 不得静默降级」 | 报错 2 条 |
| B（未收敛日志 + 有内容） | 同上 | 报错 2 条 |
| C（干净日志 + 零字节） | **`passed=true`**，「aims.out 存在且无失败标记(SCF 完成)」 | 只报 1 条（别的原因） |

→ **日志层**的问题 `validate` 能拦、`collect` 不拦（阶段已经 `collected` 并落盘 artifacts，报告期才变红）；
**零字节层**的问题两级都不拦。另：`_find_artifact` 按「run 根 → 阶段子目录首个命中」找 `aims.out`，
而 `aims_tip_relax` 也产 `aims.out`，报告期这条检查究竟读到哪一份并不确定（代码推断，未构造复现）。

## 6. P1-5 重跑把活跃作业改标 failed（成立，原报告说轻）

**复核做法**：构造 `state.json`：`aims_scf` 阶段 `status="waiting"`、`job_id=12345`、有 `remote_dir`。

| 步骤 | 实测 |
|---|---|
| 调 `run_stage`（等价一次 `smj run`） | prepare 抛「正在远端运行」保护异常 → **`status` 被写成 `failed`**；`job_id`、`remote_dir`、时间戳**全部保留**；未生成输入、未 sbatch |
| 再调 `run` | `prepare=1`、`sbatch=1`：**新作业 424242 覆盖 12345**，`clean_stale` 还删掉了远端 `.dsh_done` 等标记 |

- 代码位置：`runner.py:273-280`（保护异常）、`608-616`（`except SmjError` 无条件写 `failed`）；
  被评审提交为 `270-277` / `608-616`，行为相同。`run()` 里**只有 `finished` 的补 collect 分支**（:644-659），
  **没有 `waiting` 分支**，所以 `waiting` 必然掉进 `run_stage`。
- 原报告写「再跑可能重复提交」—— 实测是**必然重复提交**，且 `failed` 与旧 `job_id` 并存，
  manifest 读者会看到「失败 + 有个作业号」，与实际「远端仍在跑」冲突。

## 7. P1-6 仅有透射谱也工程通过（成立，原报告说轻）

**复核做法**：空 run 目录只放 `demo/azulene-26-dithiol/.../aitranss_te/TE.dat` 一份，跑
`generate_report(force=True, include_checks=True)`，并用 CLI `python -m smj.cli report -c <smj.yaml> --force --json`
对照。

| 字段 | 实测 |
|---|---|
| `state.is_done` | **false** |
| 6 个阶段 | 全 `pending`，`first_incomplete = molecule_opt` |
| `stages_inventory` | `["aitranss_te"]`（`totals.stages = 1`） |
| `verdict.engineering_passed` | **true** |
| `verdict.passed` | **true**（原报告只提 engineering_passed） |
| `verdict.errors` / `verdict.warnings` | **0** / **2**（原报告未提 warnings） |
| `verdict.physics_verdict` | `inconclusive` |
| `checks` | 仅 6 条：`aitranss_te`(info,190 点)、`efermi_source`(info)、`t_ef_resonance_attribution`、`g_in_theory_band`、`seebeck_sign`、`k_e_zt_order` |

CLI 与库路径**逐字段一致**，exit=0。

**机制**：`validate.py` 的检查表只由「已存在的产物」驱动（:68/:89/:91/:104 全是 `if ... is not None:`），
缺失阶段**不产生任何 check**，`_eng` 为空即 `engineering_passed=true`（:155-156），`passed = (errs == 0)`（:158）。
现状只堵住「TE.dat 自身缺失」（空目录对照 `errors=1`）。

**原报告说重的地方**：「报告缺完成度」—— 报告概览其实已印「阶段数 1」、6 个 pending、首个未完成阶段。
真缺口是：`is_done` 只写进 `manifest.json`（`report.py:463`）**正文不渲染**，且 `_failure_summary`（:361-366）
只收 `failed` 阶段 → 打印「无失败阶段」，与 `is_done=false` 并存。

## 8. P2-1 协议指纹不变（部分成立，需限定）

**复核做法**：对 `provenance.py` 的 `protocol_sha256` 输入逐项扰动。

| 场景 | `resolved_window` | `protocol_sha256` |
|---|---|---|
| eV 默认窗（未改） | `(-0.2619143201, 0.0001837466, -0.1149170314)` Ha | `4ead928b…` |
| 同上 + `transport.estep = 0.001` Ha | 步长 `0.0001837466 → 0.001` | **`4ead928b…`（不变）** |
| azulene demo（Hartree 窗 `-0.2075 / 0.0002 / -0.1695`）+ `estep=0.001` | 变 | `635ae0… → f8f3ef…`（**变了**） |
| 改 `energy_window` 起点 / 终点 | 变 | `367bbc0a…` / `58bad4ab…`（**变了**） |

**机制**：`provenance.py:217` 把 `energy_window_unit` 取走后，eV 口径只写入 `relative_window()` 三元组
（:220，即用户写的相对窗），随后 :226-228 把 `ener`/`estep`/`eend` 三个覆盖键 **`pop` 掉**；
而真正生效的是 `config.py:256-257` 的 `if self.estep is not None: estep = float(self.estep)`。
→ **覆盖键不进指纹**，但 `energy_window` 的起终点仍进指纹。

**原文偏差两处**：
1. 「改步长协议指纹不变」只在 **eV 相对窗 + 使用 `estep/ener/eend` 覆盖**时成立；Hartree 窗改 `estep` 会让哈希变化。
   原报告没写这个前提，读者会以为普遍成立。
2. 「原实际步长 ≈ 0.000183747 Ha」是 **eV 默认窗**解析出的步长，而 `demo/azulene` 用的是 Hartree 窗、步长
   `0.0002 Ha`。数字本身没错，但把它挂在 demo 上会误导。
3. 行号为被评审提交；当前工作区已漂到 `provenance.py:250-267`、`config.py:242-269`。

**结论：部分成立。** 缺陷真实（覆盖键不进协议指纹、`resolved_window` 只进 `run_spec` 不进 `protocol_sha256`），
但作用面比原文窄。现有测试均未覆盖该场景。

## 9. 三条「后续完善」的核实

### ① demo 谱做成发布回归样本 —— 部分成立

| 事实 | 实测 |
|---|---|
| 随仓库 TE.dat | 4 份：薁 **190 数据行**（+4 行注释，sha256 `01d89887fe02…`）；2CN 三条完成路线各 **301 行**，sha256 互异（`ce273be9…` / `b67442e9…` / `327136be…`） |
| 半成品 | `demo/2CN-…-g16tip` **无 TE.dat**，`gaussian_opt` 仍 `waiting`（job_id 26072），无 manifest |
| 温度 | 4/4 run 记录 300.0 K（`manifest.json` / `run_spec.json`） |
| TE.dat 哈希 | **3/4** 已写入 `state.json` 的 `outputs_sha256["TE.dat"]` 且与实际文件一致；薁只有 `README.md` 散文里的 `01D89887…` |
| 参考值与容差 | **完全没有**机器可读字段；薁的参考值仅 `README.md:61-65` 散文表（T=0.746、G=57.8 µS、S=+25.94 µV/K、κe=388.0 pW/K、ZT=0.0301），无容差 |
| 已有但未用 | `analysis.py:161 compare_reference()` + `:52 DEFAULT_TOLERANCES` 自 `0f7b3c4` 起就在，但 `junction_analyze.py:43`、`report.py:349` 调 `analyze_te` 均**不传 `reference`**，CLI 无入口，测试零覆盖 |
| CI | 只有 `tests` / `dist` 两个 job，**不读 `demo/`、不重算 G/G0/S/κe、不校验哈希**（`.github/workflows/ci.yml:14-37`） |

→ 方向成立，但**「保留哈希、温度」不是从零开始**：温度全有、3/4 哈希已有且一致。真正缺的是
参考值、容差、统一清单（含「以 run + 谱文件关联」）、CI 重算，以及把已有哈希纳入校验。
附带发现：`demo/azulene-26-dithiol/README.md:17,95` 写「193 点」，实际 190 点（`manifest.json` 与
`analysis.json` 均记 190）—— 若把点数写进回归基线，需以 190 为准并修 README。

### ② 增加结构与机理证据 —— 部分成立

结构归档**已经实现且在每个阶段跑**：`src/smj/structure_archive.py`，产物形如
`saved_structures/0001_aims_tip_relax_<name>.in` + `.results.json` + `manifest.json`
（复核中每次 `collect` 都生成了它们）。机理侧也已有判据（共振归因、Seebeck 符号、`g_in_theory_band`）。
因此「增加结构与机理证据」若被读作「目前没有结构证据」，前提不成立；真正缺的是
(a) 机理结论的**可判定字段**（而非文字 detail）、(b) **跨 run / 跨协议** 的对照证据。
（本条只核到「归档机制存在且在工作流中生效」，未逐一审阅归档内容是否完整可复现。）

### ③ 统一 VENDOR 维护说明 —— 成立

- `README.md:8-11`（两版本逐字相同）要求 `python scripts/sync_core.py --check` 守漂移；
  `README.md §8`（`0f7b3c4` 的 :457-464）却明说「本仓库没有该脚本」；`scripts/` 下 7 个 `.py` 确无
  `sync_core.py`。`AGENTS.md:18` 也引用了这个不存在的脚本，并指向第 8 节 —— 三处自相矛盾。
- `PROVENANCE.json`（与 `0f7b3c4` 相同）：36 条目 / 36 匹配 / 0 缺失 / 0 多余，28 个 py + 8 个模板，
  `local_patches = 10`。
- **原报告未提的关键事实**：仓库里**已有一个可执行的核验** —— `python scripts/gen_junction_provenance.py --check`
  （:187-201）。它当前 **exit=1**，唯一原因是 `repo_commit` 字段（:179 `_git_head()`）绑死在 `8eaa8fb`，
  而 36 个文件的哈希全部一致。也就是说这个检查**每次 commit 后必然变红**，红的原因与漂移无关 ——
  这比「没有核验」更麻烦：它训练人忽略红灯。
- CI 无任何 VENDOR 指纹核对（`ci.yml` 搜 `VENDOR` 零命中）；`smj verify` 只试导入 vendored 模块 + 查 3 个模板
  （`cli.py:1062-1073`），不校哈希。
- 顺带：`docs/runs-reverdict-2026-11-16.md:19`「`--check` 通过」只在 `8eaa8fb` 成立，现已过期。

## 10. 「验证范围」声明核实

| 陈述 | 结论 | 实测 |
|---|---|---|
| Windows / **Python 3.12.14** | **不成立** | `.venv` = **3.12.10**（`pyvenv.cfg:3`）；`py -0p` 无 3.12.14。若原报告在别的机器跑，应注明 |
| 267 项 pytest 全通过，18.40 s | **对 `0f7b3c4` 成立**（计时不同），**对当前工作区不成立** | `git archive` 副本：**267 passed / 19.08 s**（另一次 20.06 s）；当前工作区：**330 passed / 32.72 s** |
| `smj verify` 六项全通过 | **成立** | 实测恰好 6 项全 `ok: true`：python 3.12.10、pyyaml 6.0.3、numpy 2.5.1、junction-kernel、templates、stages（7 个阶段模块） |
| 「装 matplotlib + RDKit 后全通过」 | **部分成立（说过头）** | matplotlib 只是 `plot` extra（`pyproject.toml:31`），**rdkit 完全未声明**；CI 只装 `.[dev]`（`ci.yml:35`）→ 屏蔽两包后 `0f7b3c4` 为 **248 passed + 3 skipped**、当前树 **295 passed + 4 skipped**；`importorskip` 按模块计 1 skip，**19（当前 35）个用例根本没被收集**，CI 绿并不覆盖它们 |
| 「额外复核：原始谱重算/行序对照/合成谱/日志与作业状态」 | **无法复核** | 仓库无报告原文与这四项的执行记录；本文独立重做了这四项（§2–§7），结论与之吻合，但原报告的自述本身不可审计 |

**+63 项的去向（当前 vs `0f7b3c4`）**：未跟踪新测试文件 37（`test_open_shell_junction` 13、`test_plotting` 16、
`test_spin_transport` 8）+ HEAD 已跟踪新文件 17（`test_engine_registry` 8、`test_structure_archive` 9）
+ 已跟踪文件新增 9 —— 不是「只是多了未跟踪文件」。

## 11. 复核新增发现（原报告未提）

1. **`validate` 与阶段 `collect` 判定标准不一致（SCF）**：未收敛日志被 `validate` 判 error，却已被 `collect`
   收入并落盘（§5）。同一份证据在流水线下游才变红，属于「阶段门形同虚设」。
2. **零字节产物无人拦**：`collect` 只 `is_file()`；`validate` 只扫日志标记 → 干净日志 + 零字节矩阵可以一路绿灯。
3. **`results` 里的布尔是常量**：`normal_termination`、`has_basis_indices`、`has_omat`、`has_mos`
   写死为 `true`，会污染 `manifest.json` 的可审计性。
4. **`au_relax` 同族缺陷更重**：`src/smj/stages/au_relax.py:183` 在 metadata 里硬写 `"converged": True`，
   同样不解析任何收敛标记（虽然它额外校验了退出语、原子数、冻结原子未位移）。
5. **判据注册表元数据不自洽**：`k_e_zt_order` 标 `implemented=False`（"无闸门"），但评估器实际会 emit 结果。
6. **`_find_artifact` 的目标不确定**：报告期对 `aims.out` 的检查可能命中 `aims_tip_relax` 的那份
   （`_ARTIFACT_STAGE_HOME` 首个命中），使 `aims_scf` 判据的语义依赖目录布局（代码推断）。
7. **`electrode_build` 的回落分支会 `StageBlocked`**（§4）：`artifacts=["geometry.in"]` 与下游按名索取
   `geometry.in.next_step` 不匹配 —— 意味着「缺 `next_step`」不是危险路径，危险路径是「有 `next_step` + 未收敛」。
8. **`docs/` 与数据不一致**：薁 README 写 193 点，实际 190 点。

## 12. 修订建议（只针对报告文字，不改代码）

1. **P1-3 换举证**：删掉「日志明确写着几何尚未收敛」的 demo 依据，改为
   「`collect` 从不打开 `aims.out`（读文件跟踪实测 `aims_out_read=false`）；合成未收敛日志 + 存在的 `next_step`
   仍 `collected`，并以其建电极」。
2. **P1-4 分开写两件事**：阶段门不查收敛标记、`results` 布尔是常量；并补「报告期 `validate` 能拦日志标记、
   拦不住零字节」这一双层不一致。
3. **P1-5 改成「必然」**：实测第二次 `run` 必 `prepare+sbatch` 覆盖作业号，并保留旧 `job_id` 与 `failed` 并存。
4. **P1-6 补三个数**：顶层 `passed=true`、`warnings=2`、`totals.stages=1`；把「报告缺完成度」改为
   「报告印了阶段数与 pending，但 `is_done` 不渲染、`_failure_summary` 只收 failed 故印『无失败阶段』」。
5. **P2-1 加限定**：仅 eV 相对窗 + `estep/ener/eend` 覆盖键；Hartree 窗改 `estep` 哈希会变；
   删掉把 0.000183747 Ha 说成 demo 实际步长的表述。
6. **P1-2 补引用范围**：把 `criteria.py:476-482`（`seebeck_sign` 的 `applies_to`）与 `:273-323`
   （`physics_verdict` 归约）一并列入，否则引用链不自洽。
7. **验证范围绑定提交**：写成「`0f7b3c4`：267 passed（本机 19.08–20.06 s，原环境 18.40 s）」，
   另注明当前工作区 330 passed；环境改 3.12.10 或注明机器；依赖说明补 matplotlib 为 extra、rdkit 未声明、
   CI 覆盖为 248 passed + 3 模块级 skip。
8. **后续①②③ 按 §9 改写**：①补「温度与 3/4 哈希已有」；②改为「机理证据字段与跨 run 对照」；
   ③补「`gen_junction_provenance.py --check` 存在但绑死 `repo_commit`，需只比 files 并进 CI」。

## 13. 复现与证据

所有证据在 `.tmp/verify/`（未纳入版本库）：

| 路径 | 内容 |
|---|---|
| `.tmp/verify/lead/repro_claims_1_2.py`、`result.json` | P1-1 升/降序 + 顺序无关参考积分；P1-2 合成谱 → 判据/verdict |
| `.tmp/verify/lead/citations_0f7b3c4.txt` | 12 处 `file:line` 在被评审提交源码上的逐行打印 |
| `.tmp/verify/stages/repro_tip_collect.py`、`evidence_tip.json`、`log_tip.txt` | P1-3：未收敛/无日志/无 next_step/读文件跟踪/建电极/0f7b3c4 对照 |
| `.tmp/verify/stages/repro_scf_collect.py`、`evidence_scf.json`、`log_scf.txt` | P1-4：零字节/有内容/干净日志/三种拒绝 + `validate` 对照 |
| `.tmp/verify/runner/repro_a_waiting.py`、`repro_b_protocol.py`、`report.md` | P1-5 `waiting` 重跑；P2-1 哈希扰动矩阵 |
| `.tmp/verify/verdict/repro_only_te.py`、`keyfields_A_only_te.json`、`check_vendor_hashes.py` | P1-6 空 run；VENDOR 36 文件哈希 |
| `.tmp/verify/scope/report.md`、`pytest_full.txt`、`pytest_0f7b3c4.txt`、`tree-0f7b3c4/` | 验证范围：267 vs 330、缺包 248+3skip、`smj verify` 输出 |

最小复现片段：

```powershell
cd D:\dsh_projects\AutoSMJ
$env:PYTHONPATH = 'src'

# P1-1：同一份谱，升序 vs 降序（TE.dat 三列：E_Ha / E-E_F(eV) / T）
@'
from pathlib import Path
from smj.junction.landauer import analyze_te
src = Path("demo/azulene-26-dithiol/azulene-26-dithiol/aitranss_te/TE.dat")
lines = src.read_text().splitlines()
asc = Path(".tmp/x_asc.dat"); desc = Path(".tmp/x_desc.dat")
asc.write_text("\n".join(lines) + "\n")
desc.write_text("\n".join(reversed(lines)) + "\n")
print(analyze_te(asc)["S_uV_per_K"])    #  25.936208945458937
print(analyze_te(desc)["S_uV_per_K"])   # -25.936208945458937
'@ | .\.venv\Scripts\python.exe -

# P1-2：整窗在 E_F 之上 → physics_verdict=passed（判据在同一目录找 TE.dat）
@'
from pathlib import Path
from smj.junction.criteria import evaluate_run_criteria, physics_verdict
run = Path(".tmp/synth/aitranss_te"); run.mkdir(parents=True, exist_ok=True)
rows = "\n".join(f"{-0.20 + 0.001*i:.9f} {0.250 + 0.001*i:.6f} 0.2500000" for i in range(200))
(run / "TE.dat").write_text("#synthetic: constant T, window fully above E_F\n" + rows + "\n")
checks = evaluate_run_criteria(run, temperature=300.0, system="AZ")
print({c["name"]: (c["passed"], c.get("applicable")) for c in checks})
v = physics_verdict(checks)
print(v["physics_verdict"], v["wording"], v["not_applicable_criteria"])   # passed 支持 H0 [...]
'@ | .\.venv\Scripts\python.exe -
```

## 14. 本次复核的边界

- 不重跑任何 DFT / AIMS / 输运计算；所有「未收敛日志」「零字节产物」「空 run」都是**构造输入**，
  用于证明「验收逻辑放行」，不代表仓库里有真实坏 run。
- P1-5 用构造的 `state.json` + 假 `sbatch`，未接真实集群：证明的是**状态机与控制流**，不是 Slurm 侧后果。
- `_find_artifact` 命中哪份 `aims.out`（§11.6）为**代码推断**，未构造多份 `aims.out` 的实验。
- 结构归档只验证「机制存在并生效」，未逐份审阅归档文件的内容完整性。
- 未核对原报告是否还有其他未在此列出的主张；本复核只覆盖它自述的 7 条 + 3 条后续 + 验证范围。
