# smj-toolkit 0.1.0 项目分析报告

> 分析对象：`D:\dsh_projects\AutoSMJ`（工作区 3 个文件）
> 分析方式：解包 wheel 后逐模块精读（54 个 `.py`，约 11,175 行）+ 在本机临时目录跑通全部 11 个 CLI 子命令（dry-run 全链路）+ 关键路径用 Python 直接调用验证
> 分析日期：2026-11-11

---

## 0. 一页结论

1. **这不是一个源码仓库，而是一份"发行物 + 文档"**。工作区只有 `README.md`、`smj-standalone-package.md` 和一个 228 KB 的纯 Python wheel。源码树（`smj/`、`scripts/`、`tests/`、`pyproject.toml`）、`LICENSE`、`VENDOR.json` **全部不在工作区，也不在 wheel 里**。
2. **代码质量明显高于文档所示的水位**。fail-closed 设计贯穿始终（缺上游产物就阻塞、日志不完整就拒绝产出、写文件用原子替换、远端命令全走模板），Landauer 物理公式经我独立推导与 Wiedemann–Franz 极限校验**正确**，资源估算的 12 项常量与设计文档逐条一致。
3. **它现在跑不完整条链**。默认 `slurm-ssh` 模式下，`aims_scf` 的 `mos.aims` 永远不会被回取（`head` 计划写了但没人消费），而 `aitranss_te` 生成 `tcontrol` 时必须读本地 `aims_scf/mos.aims` —— 第 5 阶段必被 `StageBlocked` 卡死。这是**唯一一处能让默认主路径整体不可用**的缺陷（已实证）。
4. **"科学判定"这一层在 CLI 里是断开的**。`smj report` 的 `include_checks` 默认 `False` 且无 CLI 开关，`manifest.json` 的 `verdict` 实测恒为 `null`；而一旦有人把它打开，一个**只含空目录的 run 会被判定为 `passed` / "支持 H0" / "可进入发布(A5)"**（已实证）。
5. 其余问题以**状态机语义与文档不符**和**打包缺件**为主，属于可以按清单修的程度，不是设计崩坏。

一句话定位：**一个工程实现扎实、但交付物裁剪过度、且主链路有一处断点的科研工作流工具包**。

---

## 1. 交付物形态（最先要看清的事实）

### 1.1 工作区内容

| 文件 | 大小 | 说明 |
|---|---|---|
| `README.md` | 18,367 B / 339 行 | 面向使用者的说明 |
| `smj-standalone-package.md` | 18,659 B / 224 行 | 面向维护者的设计文档（架构、契约、边界） |
| `smj_toolkit-0.1.0-py3-none-any.whl` | 233,741 B | 纯 Python wheel，`setuptools 84.2.1` 构建 |

无 `.git`、无子目录、无构建产物、无测试。

### 1.2 wheel 内实有

- 54 个 `.py`（571.1 KB），约 11,175 行：
  - 核心编排层 16 文件 / 5,057 行（`cli.py` 860、`report.py` 647、`runner.py` 556、`preflight.py` 469、`config.py` 468、`review.py` 406、`remote.py` ~380、`analysis.py` 350、`state.py` 311、`context.py` 180、`transport.py` 189、`resources.py` 157）
  - 六阶段适配器 9 文件 / 847 行（`stages/`）
  - vendored `junction/` 29 文件 / 5,271 行（只读副本）
- 8 个模板：`templates/control_{orbital,relax,scf}.in`、`sub_{aims,g16,te}.sh`、`gaussian_opt_freq.gjf`、`tcontrol.template`
- `dist-info`：`METADATA`（Name `smj-toolkit`，Version `0.1.0`，Requires-Python `>=3.10`，依赖仅 `PyYAML<7,>=6`、`numpy>=1.24`；extras `plot`/`dev`/`all`）、`entry_points.txt`（`smj = smj.cli:main`）、`RECORD`（67 行）

### 1.3 wheel 内**没有**（这是本次分析最重要的发现之一）

| 缺失项 | 被谁引用 | 后果 |
|---|---|---|
| `pyproject.toml` | `README.md:40` 的 `pip install -e smj` | 只装 wheel 就无法做可编辑安装 / 重建 |
| `scripts/sync_core.py` | `README.md:10-11` | vendor 漂移门不存在 |
| `smj/_vendor/junction/VENDOR.json` | `README.md:307`、`md:46` 的 sha256+commit 台账 | 即使有脚本也无台账可比对 |
| `scripts/smoke_install.py` | `README.md:59` | 安装冒烟脚本不存在 |
| `tests/` | `README.md:249` 引用的 `smj/tests/test_remote_semantics.py` | 全部单元测试不在交付物内 |
| `LICENSE`、`License` 元数据 | 无 | 再分发授权不明（METADATA 只有 `Author:`） |
| `CHANGELOG`、`dist/`、sdist | `md` 有提及 | 无法审计版本历史 |

**结论**：README 头部的"主仓为权威源、由 `sync_core.py --check` 守漂移"这套治理机制，在**当前交付物中完全不可执行**。它描述的是一个源码仓库的使用方式，而这个工作区只有 wheel。

---

## 2. 架构与运行模型

### 2.1 六阶段固定链

```
gaussian_opt    (Slurm)  mol.gjf → mol.log/opt.xyz
   ↓
aims_tip_relax  (Slurm)  Au 针尖 + 分子弛豫 → geometry.in.next_step
   ↓
au60_build      (本地)   拼 Au(111) 电极（每侧 30 Au）→ geometry.in + electrode_metadata.json
   ↓
aims_scf        (Slurm)  FHI-aims 静态 SCF → aims.out / mos.aims / basis-indices.out
   ↓
aitranss_te     (Slurm)  tcontrol → TE.dat（透射谱）
   ↓
junction_analyze(本地)   Landauer 积分 → results.json
```

四阶段提交集群、两阶段纯本地（`au60_build`、`junction_analyze`）。`run.mode` 三态：`dry-run`（只生成输入）、`local`（本地 bash 执行）、`slurm-ssh`（上传 → `sbatch` → 轮询 → 回取）。

### 2.2 分层

| 层 | 模块 | 职责 |
|---|---|---|
| 入口 | `cli.py` | 11 个子命令、统一退出码、aux 模块动态发现 |
| 配置 | `config.py` | `smj.yaml` 解析、路径解析、校验（未知键只 warning） |
| 状态 | `state.py` | `state.json` 唯一状态源、8 状态机定义、原子写 |
| 编排 | `runner.py` | prepare/submit/wait/collect/run 五动作 |
| 适配 | `stages/*.py` | 六阶段各自的 `prepare`/`collect`，返回 `Prepared`/`Collected` |
| 执行 | `transport.py`、`remote.py` | 本地子进程 / SSH；上传过滤、清理、回取计划、轮询 |
| 适配缝 | `compat.py`（45 行） | 让 vendored 代码脱离主仓：3 个 re-export + `_LocalBackend` 占位 |
| 科学 | `_vendor/junction/{landauer,tcontrol,self_energy,criteria,validate,…}` | tcontrol 渲染、Landauer 积分、判据注册表、审计 |
| 报告 | `report.py`、`review.py`、`analysis.py` | manifest/report.md、A2 复核门、独立复核 |

**vendored 缝非常干净**：`from smj.compat import` 恰好 3 处（`_vendor/junction/preflight.py:171 / :292 / :360`），无任何 `from agenticdft` 残留；仅有注释与 usage 字符串还指向主仓（见 §4.4）。

### 2.3 状态机

8 状态：`pending → prepared → submitted → waiting → finished → collected`，另有 `failed`、`skipped`。`state.json` 是唯一事实来源，`smj status`/`report`/`audit` 全部从这里读。

---

## 3. 文档 ↔ 代码 逐条核查

### 3.1 成立的声明

| 文档声明 | 实测结果 |
|---|---|
| 11 个子命令 `init/plan/preflight/stage/run/status/report/analyze/review/audit/verify` | ✅ `smj --help` 逐个确认 |
| fail-closed：dry-run 不伪造上游产物，链路停在缺失处 | ✅ `smj run --dry-run` 在 `aims_tip_relax` 报"gaussian_opt 阶段产物缺失: opt.xyz"，退出 1 |
| `smj analyze` 写 `<run>/junction_analyze/analysis.json` | ✅ 实测生成，含 `independent_check` |
| 六阶段顺序/依赖/kind 列 | ✅ `smj plan` 表格一致 |
| 资源估算 12 项常量（60+10·N s、4+0.05·N GiB、180+30·N、8+0.10·N、4320/13.3@132、2552/8.0@132、×1.25/×1.20/×1.9、900/1800 下限、72h/128GiB、ntasks 1/4/4/1、N³+N²、N³+N^1.5） | ✅ 与 `resources.py:19-56` 逐条对上；实算 SCF@132 → `01:30:00`/`02:51:00`/`16G` |
| 装完 `smj verify` 自检 | ✅ 全通过退出 0（Python/PyYAML/numpy/vendored import/模板/6 阶段模块） |
| README §2 产物树、§5 目录结构 | ✅ 除两处（见 §3.2）与实际一致 |
| `md:160` "队列缺席仅在 squeue 成功且连续两次时判定" | ✅ 代码 `jobwait.py:166-177` 严格如此 |

### 3.2 不成立 / 失实的声明

| # | 文档位置 | 文档说 | 实际 |
|---|---|---|---|
| D1 | `README.md:11`、`md:46` | `VENDOR.json` 记录每文件 sha256 与主仓 commit，`sync_core.py --check` 守漂移 | wheel 内**两者都没有** |
| D2 | `README.md:249` | 引用 `smj/tests/test_remote_semantics.py` | 测试目录不在交付物内 |
| D3 | `README.md:40/59` | `pip install -e smj`、`smj/scripts/smoke_install.py` | 无 `pyproject.toml`、无脚本 |
| D4 | `README.md:269`、`md:73` | `aims_scf` collect 后本地有 `mos.aims` | 回取 `must` 不含它，`head` 计划从未消费（**导致 §4.1 的 P0**） |
| D5 | `README.md:272` | `aitranss_te` 产出 `orbital_match.json` | 代码从不生成、不回取；`aims_tip_relax` 的 `aims.out`/`input_summary.json`、`aims_scf` 的 `geometry.in`/`input_summary.json` 等反而未在文档产物树里列出 |
| D6 | `README.md:18` | `results.json` 键为 `G/G0`、`S`、`k_e`、`ZT_e` | 实际键：`log10_G_over_G0_like`、`S_uV_per_K`、`k_e_pW_per_K`、`ZT_electronic`（`md:75` 才是对的，README 自相矛盾） |
| D7 | `README.md:284-291` §7 表 | 只列 8 个命令 | 实际 11 个（少的 3 个在 `md:39`） |
| D8 | `md:167` §判据 | `report.md` 区分 `engineering_passed` 与 `physics_verdict` | `smj report` 的 `include_checks` 默认 `False` 且无 CLI 开关 → `verdict` **恒为 null**，`report.md` 无判定行（实测两次） |
| D9 | `README.md:186` | env 键大小写不敏感 | `config.py:379` 原样保存键名，不折叠大小写 |
| D10 | `README.md:189-197` | `cluster.mem` / `cluster.time_limit` 有默认值可用 | 4 个远端阶段的 `estimate()` 恒给出值，这两个默认值是死分支 |
| D11 | `README.md:174`、`md:198` | `species_dir` 指向"集群既有目录" | `config.py:471-472` 用本地 `Path.resolve()`，`preflight.py:272-274` 用本地 `is_dir()` 校验；但 `tcontrol.py:45-50` 是**读取内容内联进 control.in**（不是写路径）→ 实际要求**本地**有 species 副本，与文档措辞相反（此处比 exec 报告的结论温和：不存在"把 Windows 路径写进远端 control.in"的问题） |
| D12 | `README.md:307` | A2 人工确认是红线，`review.json` 会 blocked | `review.require_a2` **默认 `False`**（`config.py:213`），默认不 block；且**没有任何 CLI 能签署**（见 §4.5） |

---

## 4. 缺陷清单

严重度定义：**P0** 阻断默认主路径或产出错误科学结论；**P1** 语义错误/可导致状态或数据不一致；**P2** 健壮性与体验问题；**P3** 卫生问题。

### 4.1 P0 — `slurm-ssh` 主链在第 5 阶段必断（`mos.aims` 永不回取）

**机理**（三处代码串联）：

1. `_vendor/junction/artifacts.py:23-25` — `aims_scf` 的 `must` 清单是 `geometry.in, control.in, aims.out, basis-indices.out, sub_aims.sh, .dsh_done`，**不含 `mos.aims`**；
2. `artifacts.py:33,40-55` — `mos.aims` 被放进 `head` 计划（只取表头 2 行），注释明确写"自动归入 head，不因阶段表未列出而漏掉"；
3. `remote.py:363-366` — `download_plan()` **只返回 `must`**；`download()`（`:368-380`）逐个全量 `transport.get()`，**全程没有消费 `head`**。

于是 `aims_scf` collect 后本地目录永远没有 `mos.aims`。而 `stages/aitranss_te.py:49-57` 调用 `inputs.make_aitranss_input(...)` → `tcontrol.make_tcontrol()` → `tcontrol.py:194-195` 在 `nsaos is None` 时读 `scf/mos.aims`。

**实证**：在只有 `geometry.in`、没有 `mos.aims` 的目录上直接调用 `make_aitranss_input`：

```
raised: FileNotFoundError -> [Errno 2] No such file or directory: '...\probe_mos\mos.aims'
```

`aitranss_te.prepare` 的 `except (FileNotFoundError, KeyError, ValueError)` 会把它转成 `StageBlocked`（`aitranss_te.py:54-57`，错误文案已经写了"mos.aims / efermi / 表面索引"）。

**影响**：默认配置（`run.mode: slurm-ssh`）下六阶段链无法自动跑完，必须人工 `scp` 一份 `mos.aims` 进 `aims_scf/`。README/`md` 里描述的"一键六阶段"在集群模式下不成立。

**修法（二选一）**：① `download()` 支持 `head` 模式（`ssh head -2` 或 `tail`/`sed` 取表头），把 `retrieval_plan()['head']` 真正接上；② 或把 `mos.aims` 加入 `aims_scf.must`（代价：单次多拉约 0.96 GB，正是 P1-7 想省的）。推荐 ①，并在 `aims_scf.collect` 时校验 `mos.aims` 存在，否则 `StageError`。

### 4.2 P0 — 物理判定层：CLI 断线 + 开启后 fail-open（空格判"支持 H0"）

**两段独立事实**：

**(a) CLI 从不产出判定。** `report.py:347` `include_checks: bool = False`；`cli.py:620-631` `cmd_report` 调 `_invoke(func, config, force=args.force)`，**不传 `include_checks`**；`report.py:402` 的整块判定逻辑被跳过。实测两次（demo run 与新建空 run）`manifest.json` 都是 `"verdict": null`，`report.md` 里 grep 不到任何判定行。`md:167` 承诺的 `engineering_passed` / `physics_verdict` 分离在交付形态中**不存在**。

**(b) 一旦开启就 fail-open。** 如果调用 `build_manifest(..., include_checks=True)` 或直接用 `validate_junction_run`，而 run 目录里**只有一个空的 `aitranss_te/` 目录**、连 `TE.dat` 都没有：

```
passed: true, engineering_passed: true, errors: 0, warnings: 0
physics_verdict: "passed", wording: "支持 H0", allowed_action: "可进入发布(A5)/预测卡(A4)"
checks: []
```

路径：`report.py:403-406` 只检查 `te_dir.is_dir()` → `validate.py:140` → `criteria.py:286-287` 无 `TE.dat` 时 `return []` → `physics_verdict([])` 在 `criteria.py:197-203` 落到 `else: v = 'passed'`。"零判据 = 通过"。

**叠加风险**：`runner.collect()`（`runner.py:406-439`）对 download 缺件只写 `note`（`:433-434`）仍标 `collected` —— 一次 `scp` 抖动就能造出"已 collected、无 TE.dat"的 run。二者组合 = 报告宣称"支持 H0"。

**修法**：`criteria.py` 在无 `TE.dat` 时返回一条 error 级"缺透射谱，判据不可评估"的 check（而非 `[]`）；`physics_verdict` 对 `checks == []` 显式返回 `inconclusive`；`cli.py report` 加 `--checks` 开关并在文档里说明判定何时产出。

### 4.3 P1 — 物理判据的三个科学表述问题

| 问题 | 位置 | 说明 |
|---|---|---|
| 判据温度硬编码 300 K | `criteria.py:276` `temperature: float = 300.0`；`validate.py:140-141` 只传 `system`；`report.py:407-408`、`manifest.py:107` 同样不传 | `transport.temperature` 配置被忽略；同一份产物在 `results.json`（用真实温度）与 `manifest.verdict`（用 300 K）之间可能自相矛盾 |
| 解析失败被表述为"否证 H0" | `criteria.py:291-293` 抛出的异常被包成 `t_ef_magnitude`（error 级）失败 | 实现事故（文件损坏/格式异常）在措辞上与"物理被否证"不可区分 |
| 3 条判据注册即"已实现" | `criteria.py:31-41` `grid_independence` / `electrodes_disjoint` / `scf_converged` `implemented=True` 但从不产出 | 实际最多评估 9 条，文档口径 11 条 |

### 4.4 P1 — 状态机语义与设计文档不符（5 处）

设计文档 `md:94-106` 给了一张状态图；代码是三处散落 `if` 白名单，实际行为有 5 处偏差：

| # | 位置 | 现象 |
|---|---|---|
| S1 | `runner.py:350-354` | `wait()` 在 `slurm-ssh` 下**不读当前状态**就 `set(status="waiting")`。对已 `collected` 的阶段执行 `smj stage X --wait`（`README.md:325/336` 明确推荐）→ 状态退回 `waiting`，探测到哨兵后落到 `finished` → 丢掉 `collected`；`smj run` 的"已 collected 跳过"分支（`runner.py:498-505`）随之失效。若远端哨兵已被清理，还可能经"连续两次缺席"降为 `failed` |
| S2 | `runner.py:284-293` | `submit()` 的拦截白名单漏了 `failed`：失败阶段可以直接 `--submit` 而不 `--force`，绕过 `prepare`，提交的可能是上一轮的脚本 |
| S3 | `runner.py:498-519` vs `state.py:208-218` | `skipped` 在 `is_done()`/`first_incomplete()` 里算"已完成"，`run()` 却照跑不误 |
| S4 | `runner.py:406-439` | `collect()` 只拦 `collected`，从 `pending`/`prepared`/`failed` 都能直接进；且**无论缺件与否都写 `collected` + `results`**，状态与产物脱节 |
| S5 | `runner.py:256-266` | `prepare()` 清了 `job_id`/各类时间戳，但**不清 `results`** 与 `extra.download_*` → `--force` 重跑期间 `status`/`report` 仍展示上一轮结果 |

### 4.5 P1 — A2 人工确认门：默认不拦、且无签署入口

- `review.require_a2` 默认 `False`（`config.py:213`），A2 缺失只记 `warning`（`review.py:335-378`），`review.json` 状态仍是 `warning`、退出 0 —— 与 `README` 把 A2 称作"红线"的表述不符。
- **独立包的 CLI 没有任何签署命令**。`smj review` 只读 `a2_confirmation.json`；写入函数 `record_confirmation` 只存在于 `_vendor/junction/fingerprint.py:42` 与 `a2_viz.py:456-460`。
- 真实签署路径是未文档化的：`python -m smj._vendor.junction.a2_viz --stage <au60_build> --sign-who <name>`（且必须审计 PASS、且**无图时拒绝签署**返回 2，`a2_viz.py:450-459`）。
- 更糟的是 `_vendor/junction/review.py:230-237` 的 `pending_message` 让用户运行 `agenticdft smj-review <run> --approve` —— 这个命令在独立包里**不存在**（主仓 CLI 的残留文案）。

### 4.6 P1 — 阶段适配器的不对称

| 问题 | 位置 | 说明 |
|---|---|---|
| `aims_scf.collect` 不 fail-closed | `stages/aims_scf.py:59-85` | `aims.out` 缺 `"Have a nice day"` 只写进 `results.normal_termination=False` 与 note，**不抛错**；对比 `gaussian_opt.py:65-72` 会在同一情形下删掉 `opt.xyz` 并 `StageError`。SCF 崩了也会被标成 `collected` |
| 陈旧产物冒充新结果 | `aims_scf.py:65-66`、`:71-74` | `basis-indices.out` / `mos.aims` 只按"文件是否存在"判定，不校验时间戳/内容归属 |
| 本地模式不清理旧哨兵 | `runner.py:315-324` | 读 `.dsh_done` 优先于进程返回码；上一轮遗留的哨兵可能让超时被杀的作业"假成功" |
| `tip_relax` 缺 `next_step` 回退，但下游硬依赖 | `aims_tip_relax.py:75-81` vs `au60_build.py:31-34` | tip 阶段 `geometry.in` 可用时仍可能没有 `geometry.in.next_step`，而 `au60_build` 只认后者 |

### 4.7 P1 — 等待/缺席判定依赖 `squeue` 返回码

`jobwait.py:63-65` 在探测脚本里打印 `squeue` 的真实退出码，`parse_probe`（`:95-96`）只在 `@@QRC==0` 时置 `queue_known=True`；缺席累计（`:166-174`）严格要求 `queue_known`。**代码与 `md:160` 的规则完全一致**（查询失败不计缺席，避免网络抖动误杀）——所以这不是"代码违反文档"。真正的工程风险在于：若目标集群的 `squeue -j <已完成 id>` 返回非 0，`queue_known` 永远为假，**被杀/未起来"的作业就只能靠 `timeout_hours`（默认 6 h）兜底**，`md:94/106` 承诺的"连续两次缺席 → failed"事实上不会触发。

建议：把 `squeue` 的 `rc==1`（无效 job id）与 `rc==255`/超时区分开，前者恰恰是"作业已不在队列"的强证据，可用于缺席判定。

### 4.8 P1 — preflight 的"未检查即通过"

- `require_remote` 默认 `False` 且 CLI 不传（`cli.py:449-451`）→ SSH/远端目录/二进制三项只在 `slurm-ssh` 真连时检查；建不起连接时只记 `warning`，`ok=True`、退出 0（`preflight.py:437-463`）。
- `UNCONFIGURED` 只是模板占位（如 `inputs.py:355` `aitranss_bin=e['aitranss_bin'] or 'UNCONFIGURED'`），`preflight` 不扫描脚本内容 → 占位作业可以被提交，然后在集群上失败。
- 提交路径**没有**接 `validate_junction_run`（只有 `report.py:406` 与 `manifest.py:104` 调它）→ "物理判定"与"能不能提交"在代码里毫无耦合。

### 4.9 P2 — 数值与解析边界

| 问题 | 位置 | 说明 |
|---|---|---|
| `ZT_electronic` 可除零 | `landauer.py:83` | `1.0 / ((l0*l2)/(l1*l1) - 1.0)` 只在 `if l1` 判空；`l0*l2 == l1*l1`（平带 + 精确对称）时 `ZeroDivisionError` |
| 能量列不做 `D→E` 替换 | `landauer.py:40-44` | Fortran 风格 `1.23D+00` 只在通道列被替换，能量列 `float()` 直接失败 → 该行被静默丢弃（`except ValueError: pass`） |
| `walltime_seconds` 的 `MM:SS` 分支按 `HH:MM` 算 | `resources.py:67-77` | docstring 说支持 `MM:SS`，实现 `parts[0]*3600 + parts[1]*60` → `"30:00"` = 108000 s（30 小时）而非 1800 s |
| 预算检查在正常路径上是空操作 | `resources.py:112-120` vs `:139-157` | `estimate()` 先把值钳到 72 h/128 GiB，`check_smj_budget()` 再拿同一上限比较 → 永远 `ok`；超限信号只剩 `capped` 字段 |
| `junction_analyze` 未列入资源表 | `resources.py:89-110` | 只特判了 `au60_build`，`junction_analyze` 会落到 `{'unknown': True}`（虽然当前调用方都跳过本地阶段） |
| CODATA 常数偏旧 | `landauer.py` | `8.6173324e-5` 为 2010 版；与 `independent_check` 用的 2018 值不同源（量级无影响） |

### 4.10 P2 — 报告层的静默行为

- `report.py:704-706`：`manifest.json` 已存在且未加 `--force` 时，用**旧** manifest 渲染 `report.md` → 新完成的阶段静默不出现在报告里（`README.md:278` 只说了"默认不覆盖"）。
- 原子写使用**固定 `.tmp` 文件名**（同目录并发会撞）。
- `smj report`/`analyze`/`review` 走 `cli.py` 的 `_optional_module` + `_find_callable` 动态发现；若模块缺判定键，`_truthy_ok()`（`cli.py:144`）在四类键全缺时**默认返回 True** → aux 模块返回无键 dict 会被判"通过"。

### 4.11 P3 — 卫生问题

- `self_energy.py:187`、`a2_viz.py:387` 的 usage 仍写 `python -m agenticdft.tools.junction.*`。
- `molecules.py:53` 提示 `pip install "agenticdft[dft]"`（独立包不适用）。
- `METADATA` 无 `License` / `Project-URL` / `Home-page`；无 `README` 项目链接。
- `dist-info` RECORD 仅 67 行，无 `INSTALLER` 之外的溯源信息（无 git commit、无 build 时间戳）。
- `README.md:307` 第 307 行的 A2 描述与默认配置不一致（见 §4.5）。

---

## 5. 设计亮点（不应被缺陷清单淹没的部分）

1. **fail-closed 贯彻到底**：`gaussian_opt.collect` 在没有 `"Normal termination"` 时删掉 `opt.xyz` 并报错（"宁可不产出，也不静默交给下游"）；`context.require_stage_file` 在上游缺失时直接 `StageBlocked`；`parse_efermi` 三处都缺时 fail-closed。
2. **vendored 缝只有 3 处**：`smj/compat.py` 45 行撑起整个脱离主仓的适配层（3 个 re-export + `_LocalBackend` 占位）。这是本次审阅里最"干净"的工程决策。
3. **远端命令全部模板化 + `str.format`**（`tcontrol.py:19-26`），无手写字符串拼接；`shlex.quote` 到位；远端路径含空白直接拒绝（`transport.py:169-173`）；删除范围不越出 run 目录。
4. **完成判定不依赖 Slurm 退出码**：远端脚本末尾写 `.dsh_done`（首列退出码），轮询以哨兵为准（`md:159` 的"Slurm 的 shell 退出码不可靠"落在了代码里）。
5. **回取策略有真实标定**：`artifacts.py:1-6` 记录了"单 run 集群侧 >2.2 GB（mos.aims 0.96 GB + omat.aims 0.48 GB + restart 0.76 GB）"的实测来源，头表回取的设计动机写得很清楚（只可惜 `head` 没接上，即 §4.1）。
6. **资源估算有来源**：`resources.py:6-11` 的幂律锚点（132 原子 4320 s / 13.3 GiB）与 `TIME_LIMIT_FACTOR = 1.9` 的理由（"收敛后写产物还要 ~25 min，132 原子实测曾因 `--time` 卡在产物写完之前被 Slurm 杀掉"）都写在代码里。
7. **独立复核**：`analysis.py` 的 `independent_check` 用独立的解析器 + Simpson 法 + CODATA 2018 重算 L0/L1/L2/S/ZT，实测相对差 ~1e-8。虽非完全独立（非均匀网格会回退到与主算相同的梯形法，且只校验 L0/L1/L2），但这个"第二实现"的存在说明作者对数值正确性有自觉。
8. **Landauer 物理正确**：`-∂f/∂E` 解析式（`landauer.py:54-59`）、Mott 公式（`:81`）、`k_e` 的 `1.60217733e-7` 因子（`:82`）、`ZT_e = 1/(L0L2/L1² − 1)`（`:83`）均正确。我独立核算：平带 `T=1`、300 K 时 `k_e = 567.9 pW/K`，与 Wiedemann–Franz 极限 `L·G0·T = 5.673e-10 W/K` 相差 0.1%（CODATA 差异），公式自洽。
9. **文档本身质量高**：`md` 有明确的状态图、语义约定表（`md:155-167`）、"不承诺"边界（`README.md:20`），甚至把设计债写进了文档（如退出码 2 的双重语义）。

---

## 6. 实测记录（dry-run 全链路，临时目录，未污染工作区）

| 命令 | 结果 | 退出码 |
|---|---|---|
| `smj init . --name demo --xyz molecule.xyz --mode dry-run` | 写 `smj.yaml`（run/molecule/electrode/transport/cluster/paths/review 七段） | 0 |
| `smj --help` | 11 个子命令 | 0 |
| `smj plan` | 六阶段表，`kind` = slurm/slurm/local/slurm/slurm/local | 0 |
| `smj stage gaussian_opt --prepare` | 写 `mol.gjf`（b3lyp/6-31g(d)）、`sub_g16.sh`、`input_summary.json`，状态 `prepared` | 0 |
| `smj run --dry-run` | `gaussian_opt` 成功 → `aims_tip_relax` 被阻塞（缺 `opt.xyz`）→ `failed` | **1** |
| `smj status` / `--json` | state 表格 / JSON 均正常 | 0 |
| `smj preflight` | `[error] 未配置 species_dir`，结论未通过；其余为 info/warning | 1 |
| `smj review` | `status=warning`、`a2.present=false` | 0 |
| `smj audit` | 通过（error=0, warning=0） | 0 |
| `smj report` | 写 `manifest.json`（`smj.report.v1`）+ `report.md`（7 节） | 0 |
| `smj analyze --te <合成 TE.dat>` | `analysis.json`：L0=0.0225137、S=−34.97 µV/K、ZT_e=0.03787、`independent_check.ok=true` | 0 |
| `smj stage junction_analyze --run` | dry-run 只 prepare，打印"未提交" | 0 |
| `smj verify` | Python 3.12.10 / PyYAML 6.0.3 / numpy 2.5.1 / vendored import / 模板 / 6 阶段模块全部通过 | 0 |
| 未知配置段 `bogus_section` | 仅打印"未知配置段，已忽略"，继续 | 0 |
| 空 run（只有空 `aitranss_te/`）+ `validate_junction_run` | `passed: true`、`支持 H0`、`可进入发布` | — |
| 无 `mos.aims` 的 `make_aitranss_input` | `FileNotFoundError`（→ `StageBlocked`） | — |

环境说明：本机无 FHI-aims / Gaussian / Slurm，因此**集群真跑与真实量化计算未验证**；`matplotlib` 未安装，A2 出图路径未执行。全部测试在 `%TEMP%` 下，工作区文件未改动（本报告除外）。

---

## 7. 修复优先级建议

| 顺序 | 动作 | 对应 | 工作量 |
|---|---|---|---|
| 1 | `download()` 接上 `retrieval_plan()['head']`（或把 `mos.aims` 移入 `must`），并在 `aims_scf.collect` 缺 `mos.aims` 时 `StageError` | §4.1 | 小 |
| 2 | `criteria.evaluate_run_criteria` 无 `TE.dat` 时返回 error 级 check；`physics_verdict([])` → `inconclusive`；`cli report` 加 `--checks` 并落盘判定 | §4.2 | 小 |
| 3 | `wait()` 先读状态，只在合法前态（`submitted`/`waiting`）写 `waiting`；`submit` 白名单补 `failed`；`collect` 增加前置状态与产物校验；`prepare` 清 `results` | §4.4 | 中 |
| 4 | 判据温度从 `transport.temperature` 透传（`validate`/`report`/`manifest` 三处） | §4.3 | 小 |
| 5 | 给 `require_a2` 一个 CLI 签署入口（`smj review --approve --who`），修掉 `review.py:234` 的 `agenticdft smj-review` 文案 | §4.5 | 小 |
| 6 | `aims_scf.collect` 与 `gaussian_opt` 对齐（未正常结束即 `StageError`） | §4.6 | 小 |
| 7 | `squeue` rc 细分（1 = 已出队，可计缺席；255/超时 = 未知） | §4.7 | 中 |
| 8 | 把 `pyproject.toml`、`scripts/`、`tests/`、`VENDOR.json`、`LICENSE` 一并放进 sdist；wheel 至少保留 `VENDOR.json` 与 `LICENSE` | §1.3 | 中 |
| 9 | 文档修订：D1–D12 全部（尤其 `README:18` 键名、`README:269` 产物表、`README:307` A2 红线表述） | §3.2 | 小 |
| 10 | `ZT` 除零保护、能量列 `D→E`、`walltime_seconds` `MM:SS` 修正、固定 `.tmp` 改名 | §4.9/§4.10 | 小 |

---

## 8. 本次分析未覆盖的边界

- 真实集群：`sbatch`/`squeue`/`scp` 的实际返回码与 `sshpass`/`ProxyJump` 行为（代码逻辑已读，行为未测）。
- 真实计算：Gaussian 与 FHI-aims 的产物格式变体（如不同版本文本的 `nsaos`、`efermi` 写法）。
- A2 出图：`matplotlib` 未安装，`render_figure` 与"无图拒绝签署"路径未执行。
- 与主仓 `agenticdft/agenticdft/tools/junction/` 的逐字节漂移：`VENDOR.json` 缺失，无法比对。
- 性能与并发：全部为单进程 dry-run，未做长作业/多 run 并发测试。

---

## 9. 补遗：0.1.1 的修复记录（2026-11-11）

本报告分析的是 0.1.0（wheel 还原）。分析结束后，工作区已改造为可编辑源码仓库
（`src/` 布局 + `pyproject.toml` + `tests/`），并按 `smj-修复计划.md` 修掉了四个优先项。
以下结论不影响本报告对 0.1.0 的判断，但引用时请区分版本。

| 编号 | 0.1.0 的问题 | 0.1.1 的修复 | 提交 | 测试 |
|---|---|---|---|---|
| F2（P0） | 缺 `TE.dat` 时 `validate_junction_run` 返回 `passed=true` / "支持 H0" | 登记 `te_present`（engineering，gate=error，`missing=True`）；证据缺失归入 `inconclusive`（新字段 `missing_evidence_criteria`），不再算否证 | `bc7febf` | `tests/test_criteria_failopen.py`（7） |
| F1（P0） | `retrieval_plan` 的 `head` 从不被消费，mos.aims 不落地，第 5 阶段必断 | `remote.download()` 按 must 全量取 + head 用 `head -n 2` 取表头；`collect` 缺件直接 `StageError`（不再标 collected） | `f7d3970` | `tests/test_remote_mos_head.py`（4） |
| F3（P1） | 状态机 5 处偏差（wait 回退、submit 黑名单漏 failed、prepare 不清旧结果、collect 缺件、skipped 语义） | `DONE_STATUSES` 统一判定；wait 拒绝已结束状态；submit 白名单 + UNCONFIGURED 拦截；prepare 清 `results`/`download_*` | `fd52550` | `tests/test_state_machine.py`（8） |
| F4（P1） | 独立包无 A2 签署入口，`require_a2` 默认 False | 新增 `smj/a2.py` 与 `smj review --approve --who <名字>`（审计不过/出图失败退出码 2）；`require_a2` 三处默认改 True；修 vendored `pending_message` 文案 | `fe48c29` | `tests/test_a2_approve.py`（13） |

打包与可追溯性：新增 `scripts/gen_vendor_json.py`（`VENDOR.json`：36 文件 sha256 + 3 条本地补丁）
与 `scripts/build_wheel.py`（纯标准库 PEP 427 构建，可复现）；产物
`dist/smj_toolkit-0.1.1-py3-none-any.whl`（238193 B，69 条目，RECORD 自检通过），
从安装树跑 `smj verify` 六项全 ok，`pytest` 32 passed。

**仍未修复**（计划中列为 P2/P3，0.1.1 未动）：资源估算先钳到 72h/128GiB 再比较导致
`check_smj_budget` 恒 ok、`jobwait` 仅在 `queue_known` 时累计缺席、`preflight` 的
`require_remote` 默认 False、`ZT` 的 `l1=0` 除零、能量列 `D→E` 只覆盖通道列、
README 之外的 `smj-standalone-package.md` 未同步修订。

## 10. 0.1.2 补遗（库电极 + 背层固定弛豫）

0.1.2 依据 FHI-aims 手册 §5.3.1 补上了两处与官方流程的差距，两者都在
commit `00294d9` 里：

- **库电极**：`electrode.source: library`（默认）直接读 `electrodes.library` 里的
  `fcc111.au.clusters.37_atoms/fcc111.cluster.natoms37.x1.aims`，用 `mount_electrode`
  的四元数贴装；库文件不进 wheel，metadata 记 `library_cluster.{file,relpath,sha256,
  natoms,nlayers_suggested}` 作为 provenance。`source: builtin` 保留 0.1.1 的自建三角岛。
  库文件缺失时的报错给出三条出路（`library_dir` / `SMJ_ELECTRODE_LIBRARY` / `builtin`）。
- **背层固定弛豫**：新增第 4 阶段 `au_relax`（aims BFGS），`prepare` 给两侧背层整层加
  `constrain_relaxation .true.`，`collect` 要求 `Have a nice day` + `geometry.in.next_step`、
  校验冻结原子位移 ≤1e-4 Å、且用 det-v1 重算的表面三元组必须与 `au60_build` metadata 一致，
  然后才允许 `aims_scf`/`aitranss_te` 消费；A2 审计对象也随之改为弛豫后的几何
  （`find_build_dir` 优先 `au_relax`）。`electrode.relax: false` 时该阶段本地直通并写
  `relax_skipped.json`。
- 测试从 32 条扩到 **48 条**（`tests/test_library_au_relax.py` 16 条：解析硬规则、
  au60_build 库装配、au_relax 的四条 fail-closed、relax=false 直通、nlayers auto）；
  `smj verify` 的 stages 项改为动态检查 7 个阶段。
- 产物：`dist/smj_toolkit-0.1.2-py3-none-any.whl`（251694 B，71 条目，RECORD 自检通过；
  解包树 `smj verify` 全 ok），`VENDOR.json` 记 5 条本地补丁（新增 F-LIB、F-RELAX）。
