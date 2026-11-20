# smj-toolkit 0.2.2rc1 验证记录

> 目的：记录 0.2.2rc1 的**可验证状态**与证据来源。
> 与 `validation-0.2.0.md` 不同，本轮**没有新增集群计算**：AGENTS.md 规定只跑用户点名的
> 体系，0.2.2rc1 的数值证据来自随仓库发布的 demo 谱、外部评审复核与两边对账。

## 1. 版本与范围

| 项 | 值 |
|---|---|
| 版本 | smj-toolkit 0.2.2rc1（tag `v0.2.2rc1`） |
| 平台 | Windows 11 / Python 3.12.10 |
| 依赖 | PyYAML 6.0.3、numpy 2.5.1（运行时）；matplotlib 3.11.1、rdkit 2026.09.1（可选，测试时已装） |
| 本版主题 | 外部评审 7 条问题 + 3 条后续的修复；junction 内核改由本仓维护（A 方案） |

## 2. 修了什么（对应 `docs/fix-plan-2026-11-20.md`）

| 编号 | 问题 | 现在 |
|---|---|---|
| R1 | 倒序 TE.dat 静默翻转 Seebeck 符号 | 按 E−E_F 排序后积分，结果带 `reordered` / `warnings`；倒序与升序同值 |
| R2/R3 | 能量窗完全偏离 E_F 仍物理通过 | `te_window_covers_efermi`（error 闸门）；全部判别性闸门不适用时判 `inconclusive`；`prepare` 直接拒绝不含 E_F 的窗口 |
| R4 | 针尖弛豫未验收敛 | 解析 `aims.out` 尾部，未收敛 / 缺日志 / 缺 next_step 一律拒绝采集 |
| R5 | 静态 SCF 可采集空矩阵 | 要求 `Self-consistency cycle converged` + 产物头部/规模校验；`results` 改为观测值 |
| R6 | 重跑把活跃作业误标 failed | 新增 `ActiveJobError`；`submitted/waiting` 走 wait→collect，超时保持 waiting |
| R7 | 只有透射谱也工程通过 | 合成 `run_incomplete` 判据，`passed` 现在要求整链完成 |
| R8 | 改步长协议指纹不变 | 显式 `ener/estep/eend` 进协议哈希；默认配置哈希形状不变 |
| R9 | demo 谱无回归基线 | `demo/demo-reference.json` + 16 例回归测试（哈希 + 容差 + 篡改检测） |
| R10 | VENDOR 维护说明自相矛盾 | 台账改 `PROVENANCE.json`，`--check` 不绑 commit；删 `sync_core.py` |
| R12/R13 | 判据元数据不自洽、`aims.out` 定位歧义、两级判定不一致 | 元数据改 `implemented=True`；`aims.out` 只认 `aims_scf`；新增两级一致性用例 |

## 3. 验证证据

| 检查 | 命令 | 结果 |
|---|---|---|
| 单元/集成测试 | `python -m pytest -q -p no:cacheprovider` | **397 passed**（≈19 s） |
| 安装自检 | `smj verify` | **7 项全 ok**（`kernel-provenance`：36 文件 / 13 条历史补丁） |
| 内核台账 | `python scripts/gen_junction_provenance.py --check` | 一致（`repo_commit` 不参与比对） |
| demo 数值基线 | `python scripts/gen_demo_reference.py --check` | 4 份谱一致 |
| 与上游对账 | `python scripts/crosscheck_agenticdft.py <TE.dat>` | 薁谱 L0/S/κe/ZT 两边逐位一致 |
| 构建产物 | `python scripts/build_wheel.py --check --sdist` | wheel RECORD 自检通过 |

demo 4 份谱的基准值（T = 300 K，全部 `input_energy_order=ascending`）：

| 谱 | 点数 | L0 | S (µV/K) | κe (pW/K) | ZT_e |
|---|---|---|---|---|---|
| azulene-26-dithiol | 190 | 0.746147187 | +25.936209 | 388.004467 | 0.030068861 |
| 2CN-656-4-thiomethylphenyl | 301 | 0.000477053 | +2.796263 | — | 0.000301646 |
| 2CN-…-aimsopt | 301 | 0.001320209 | −35.133531 | 1.266416885 | 0.029910651 |
| 2CN-…-fullrelax | 301 | 0.001284826 | −29.091841 | 1.127745509 | 0.022412585 |

## 4. 边界

- **没有新增集群验证**：0.2.2rc1 不改变 DFT / aitranss 端的计算内容，改动集中在判定、
  验收与记账层；`validation-0.2.0.md` 的三套体系结论仍代表整链的科学证据。
- demo 谱是示例数据，不是生产级结论；`S` 的符号本版已不再依赖文件行序。
- `aims_tip_relax` 的收敛验收只读日志尾部 1 MiB：超长日志的早期内容不参与判定（判定只看最后一步）。
- 与 AgenticDFT 不再共享代码；若两边对同一 TE.dat 给出不同数值，以
  `scripts/crosscheck_agenticdft.py` 对账为入口排查。
- `rdkit` 未列入运行时依赖，只是 `tests/` 里 `importorskip` 的可选路径；CI 只装 `.[dev]`。

## 5. 复现

```powershell
cd D:\dsh_projects\AutoSMJ
.venv\Scripts\python.exe -m pytest -q -p no:cacheprovider
.venv\Scripts\python.exe -m smj.cli verify
.venv\Scripts\python.exe scripts\gen_junction_provenance.py --check
.venv\Scripts\python.exe scripts\gen_demo_reference.py --check
.venv\Scripts\python.exe scripts\build_wheel.py --check --sdist
```

## 6. 发布产物

| 文件 | 字节 | SHA-256 |
|---|---|---|
| `dist/smj_toolkit-0.2.2rc1-py3-none-any.whl` | 395818 | `B7436AC34D52A1FDAB9733056F49509541763605FE8D0CD56477023D40143BC2` |
| `dist/smj_toolkit-0.2.2rc1.tar.gz` | 474907 | `9486E99B422AD130B42EFE3EE64E09FA5473BE935FBA0A4ECCFFC0F1857D5E87` |

wheel 的 RECORD 自检 91 条；sdist 内 kernel/PROVENANCE.json 与工作树一致，可由源码独立重建。
