# smj-toolkit 0.1.2 修复计划：库电极 + 背层固定弛豫

> 触发：2026-11-11 用户确认（m01166）「默认用 33/37；要一起加“背层固定弛豫”阶段，按 0.1.2 的计划落实执行」。
> 基线：0.1.1，commit `2c0f86e`。本文档先定方案，再逐里程碑实施；每个里程碑完成即更新本文档状态。

## 1. 目标与依据

FHI-aims 手册 §5.3.1（https://fhi-aims.org/uploads/manual/Ch5/S3.html）的要求：

1. 电极要用 `electrodes.library/` 里的**预定义 Au 团簇**，不要自造；
2. 分子—电极装配后，**固定两个背层做一次几何弛豫**，防止电极在后续计算里塌陷；
3. `tcontrol` 的 `-nlayers` 取库文件头部建议值。

我们当前的差距：

| 现状 | 手册要求 | 0.1.2 动作 |
| --- | --- | --- |
| `Electrode111` 自造三角岛（默认 30/侧，截断到 35） | 用库里的紧凑壳层团簇 | 新增 `electrode.source: library`，默认 37_atoms |
| 无弛豫，直接 SCF | 背层固定后弛豫一次 | 新增第 4 阶段 `au_relax`（7 阶段链） |
| `nlayers` 固定 3 | 取库头部建议值 | 读库头 `-nlayers`，与配置不一致时告警并记录 |

## 2. 关键设计决策

- **D1 库文件不进 wheel**：AIMS 分发许可 + 版本漂移。只记录来源路径与 sha256（metadata `library_cluster`），运行时从 `electrode.library_dir` 读；本地不存在但 `mode=slurm-ssh` 时可经 ssh `cat` 取回并缓存到 `<run>/electrodes/`。
- **D2 默认 `source: library` + 37_atoms/x1**：与用户决定一致；老行为用 `source: builtin` 恢复。没有库时 `au60_build` 失败并给出明确指引（`SMJ_ELECTRODE_LIBRARY` / `library_dir` / `source: builtin` 三条出路），不静默回退。
- **D3 弛豫阶段独立成阶段**（而不是塞进 aims_scf）：`au_relax` 有独立产物（弛豫后 `geometry.in`、`aims.out`、`geometry.in.next_step`），A2 指纹与签署对象随之为**弛豫后**的几何；失败/未收敛可单独重跑，不必重跑 SCF。
- **D4 背层冻结**：`write_aims_geometry(..., constrained=...)` 给两侧背层全部原子加 `constrain_relaxation .true.`；背层索引由 `au60_build` 在 metadata 里写好（`left_backplane_indices` / `right_backplane_indices`，1-based，整层而非 3 个表面原子），弛豫阶段只消费、不重算。
- **D5 弛豫后一致性断言（fail-closed）**：`au_relax.collect` 用 det-v1 在弛豫后几何上重算左右表面三元组，必须与 metadata 完全一致——背层被冻结时这是恒等式；不一致说明约束失效或原子迁移，直接 StageError，不放行到 SCF。
- **D6 `electrode.relax: false` 是本地直通**：prepare 只拷贝几何与 metadata、不写作业脚本（`Prepared.local_only=True`），runner 不提交、不下载；A2 直接在 `au_relax/` 上签署。用于算力紧张时保留旧行为。
- **D7 A2 目标迁移**：`find_build_dir` 优先返回“含 geometry.in + electrode_metadata.json 的 `au_relax`”，否则回退 `auN_build`；`report`/`review`/`a2` 共用这一规则，报告中同时打印签署目录。
- **D8 不改公开 CLI 语义**：只新增 `--library-dir/--library-cluster/--source`（init/stage 覆盖）与 `smj stage au_relax`；`smj plan` 自动多一行。

## 3. 配置与数据结构

```yaml
electrode:
  source: library                  # library | builtin
  library_dir: null                # 本地目录/文件，或集群上的绝对路径（slurm-ssh 时经 ssh cat）
  library_cluster: fcc111.au.clusters.37_atoms/fcc111.cluster.natoms37.x1.aims
  au_per_electrode: 30             # 仅 builtin 生效
  nlayers: 3                       # 与库头建议值不一致时告警
  relax: true
  relax_force_threshold: "1.e-2"   # eV/Å，写进 control.in 的 relax_geometry bfgs
```

metadata 新增键：

- `electrode_source`: `library` / `builtin`
- `library_cluster`: `{file, relpath, sha256, natoms, nlayers_suggested, fetched_from}`
- `left_backplane_indices` / `right_backplane_indices`：背层全部原子（1-based）
- `relax`: `{enabled, converged, steps, force_threshold, input_fingerprint}`（由 `au_relax` 写入）

阶段链（7）：

`gaussian_opt → aims_tip_relax → au60_build → au_relax → aims_scf → aitranss_te → junction_analyze`

## 4. 里程碑

| # | 交付 | 关键文件 | 测试 |
| --- | --- | --- | --- |
| M1 | 库电极解析与装配 | `electrode_library.py`(新)、`electrode.py`、`au60_build.py`、`config.py` | `test_electrode_library.py` |
| M2 | au_relax 阶段 | `stages/au_relax.py`(新)、`inputs.py`、`state.py`、`stages/__init__.py`、`resources.py`、`preflight.py`、`remote.py`、`artifacts.py` | `test_au_relax_stage.py` |
| M3 | 下游接线与 A2 迁移 | `aims_scf.py`、`aitranss_te.py`、`review.py`、`report.py`、`a2.py`、`runner.py`(local_only) | 既有 32 条 + 新增 |
| M4 | 文档与自检 | `README.md`、`cli.py`(verify 7)、`scripts/gen_vendor_json.py` | `smj verify` |
| M5 | 0.1.2 打包 | `pyproject.toml`、`__init__.py`、wheel、`VENDOR.json` | `build_wheel.py --check` + 安装树 verify |
| M6 | 收尾提交 | 计划/分析报告补遗、git commit | 全量 pytest |

## 5. 测试矩阵

1. 库文件解析：层数/头部 nlayers/元素一致/共面/背层 ≥3/接触原子唯一；坏文件（缺 atom、混元素、两层缺失）全部拒绝。
2. 装配：37_atoms 默认 → 每侧 37、总原子 = 分子 + 74；metadata 的 provenance、背层索引、surface 三元组齐全；`source: builtin` 与 0.1.1 逐字节一致（golden）。
3. au_relax prepare：geometry.in 的背层行后带 `constrain_relaxation .true.`；control.in 含 `relax_geometry bfgs`；`relax: false` 时 script=None、local_only=True。
4. au_relax collect：next_step 覆盖 geometry.in、指纹更新、三元组一致；缺 "Have a nice day" / 缺 next_step / 原子数变化 / 三元组变化 → StageError。
5. 链路：`smj plan` 7 行；aims_scf 上游指向 au_relax；aitranss 优先读 au_relax 的 geometry+metadata；A2 默认签署目录是 au_relax。
6. 回归：`smj verify`、`smj run --dry-run`（无库时 au_relax 之前不报错）、0.1.1 的 32 条测试全绿（必要的 6→7 断言同步）。

## 6. 验收 DoD

- `python -m pytest` 全绿（≥ 40 条）；
- `smj verify` 7 个阶段模块、模板、vendored junction 全 ok；
- `dist/smj_toolkit-0.1.2-py3-none-any.whl` 构建通过、RECORD 自检、解包树 verify 通过；
- 用真实库文件（.tmp/aims_lib 的 37_atoms x1）跑通 `au60_build prepare` 与 `au_relax prepare`（不提交集群）；
- 文档三处一致：README、本计划、分析报告补遗。

## 7. 风险与回滚

| 风险 | 处置 |
| --- | --- |
| 弛豫把分子—电极接触拉坏 | 背层冻结 + 接触几何仍由 au60_build 决定；`relax: false` 可整体跳过 |
| det-v1 三元组在弛豫后变化 | collect 直接拒绝（D5），用户可 `--force` 重跑或回退 builtin |
| 集群上没有可读库路径 | 报错列出 `SMJ_ELECTRODE_LIBRARY`、`library_dir`、`source: builtin` 三条出路 |
| 7 阶段破坏旧 run 的 state | `STAGE_ORDER` 增加但旧 state 的 6 阶段记录仍可读；缺失的 au_relax 由 `--run` 从该阶段开始补 |
| 回滚 | 0.1.1 wheel 与 commit `2c0f86e` 保留；`electrode.source: builtin` + `relax: false` 等价旧行为 |

## 8. 明确不做

- 不做电极库文件的自动下载/版本管理（只记录 sha256）；
- 不做 fcc001 (001) 取向与 Ag/Cu 等其它金属的库支持（解析器接受任意元素，但默认只验证 Au）；
- 不做弛豫的收敛参数自动扫描（只透传 `relax_force_threshold`）；
- 不改 `aims_scf` 之后的所有物理判据。

## 9. 状态

| 里程碑 | 状态 | commit |
| --- | --- | --- |
| M1 库电极 | ✅ 代码+测试（tests/test_library_au_relax.py） | — |
| M2 au_relax | ✅ 代码+测试 | — |
| M3 下游/A2 | ✅ aims_scf/aitranss_te/review/report 已接线 | — |
| M4 文档 | ✅ README 七阶段/库电极/背层弛豫 | — |
| M5 wheel | ✅ 251694 B/71 条目，解包树 verify 全 ok | `87e8cd0` |
| M6 收尾 | ✅ README/分析报告补遗，全量 pytest 48 passed | `87e8cd0` |
