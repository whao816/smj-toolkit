# SMJ 独立打包与逐阶段独立运行设计

> 本文回答两个问题：**单分子结（SMJ）计算链能否单独打包？六个阶段能否各自独立运行？**
> 面向要复用/审计这条链的研究者，只写契约、状态与边界，不写营销话术。
> 截至 2026-11-11，`smj/` 包内的配置、状态机、远端语义、CLI、六个阶段适配器、报告/分析/
> 预检与 vendored 核心均已落地；本文的命令与参数已按 `smj --help` / `smj/smj/cli.py`
> 核对，产物名与覆盖参数已按 `smj/smj/stages/` 源码核对。**本文与代码冲突时以代码为准**。
> 主仓仍是非 VASP 运行产物与结果的权威落点（`work/<name>/`）。

## 0. 结论（30 秒版）

| 问题 | 结论 | 边界 |
|---|---|---|
| 能否单独打包？ | **能**。发行名 `smj-toolkit`，import 名 `smj`，console script `smj`；`pip install -e smj` 或 `python -m build` 出 wheel | 只打包**代码与输入模板**：Gaussian / FHI-aims / aitranss 二进制、许可证、species 基组不随包分发；`smj._vendor.junction` 是主仓的只读副本，由 `scripts/sync_core.py --check` 守漂移 |
| 每阶段能否独立运行？ | **能**。阶段之间只通过 `<run>/<stage>/` 里的文件耦合；可单独 `--prepare`、单独 `--submit/--wait/--collect`，也可用 `--xyz/--geometry/--scf-dir/--te/--metadata` 绕过上游直接跑单阶段 | 远端阶段（4 个）仍要集群与对应程序；`au60_build`、`junction_analyze` 纯本地；`dry-run` 只到 `prepared` |
| 能否不装主仓？ | **能**。运行时依赖只有 PyYAML + numpy（+ 可选 matplotlib），vendored 代码只改 3 处 import | 想同步 vendored 代码或跑主仓 CI 时仍需要主仓 checkout |

一句话：**打包是「可移植的执行壳」，不是「自带 DFT 软件的一体机」**；独立运行是
「文件级解耦」，不是「阶段之间没有物理依赖」——每个阶段的输入契约仍然是硬约束。

## 1. 包结构与职责边界

```text
AgenticDFT 主仓（权威源，只读）
  agenticdft/agenticdft/tools/junction/      # 单分子结算法 + 模板 + inputs.py
        │  python smj/scripts/sync_core.py [--check]
        ▼
smj/  （dist: smj-toolkit, import: smj, console: smj）
  pyproject.toml
  smj/
    config.py        # smj.yaml 解析、RunInfo/Molecule/Electrode/Transport/Cluster/Paths/Review
    state.py         # STAGE_ORDER、状态机、RunState、RunLayout
    context.py       # StageContext：路径解析、上游文件、覆盖参数、资源 env
    resources.py     # 资源估算（纯函数，无 IO）
    compat.py        # vendored 代码唯一适配缝（get_backend / 资源函数）
    transport.py     # Transport 四原语；LocalTransport / SshTransport
    remote.py        # Slurm 提交、残留清理、哨兵轮询、产物回取
    runner.py        # StageRunner：prepare/submit/wait/collect/run
    cli.py           # 11 个子命令（argparse，console script 入口）
    stages/          # 六个阶段适配器（base.py 定义契约）
    preflight.py     # 外部程序/species/许可证/资源预检
    report.py        # state.json + 产物 -> manifest.json + report.md
    analysis.py      # TE.dat -> G / S / k_e / ZT（薄封装 landauer）
    review.py        # 结构检查点（tips / electrodes）签署
    testing.py       # FakeTransport：内存提交/清理演练（轮询探针需子类模拟）
    _vendor/junction/  # 主仓副本（含 VENDOR.json 的 sha256 台账）
  scripts/sync_core.py
  tests/
```

三条边界：

1. **主仓是权威源**。vendored 目录不手工改；改算法去主仓 `agenticdft/agenticdft/tools/junction/`，
   再 `python smj/scripts/sync_core.py` 同步并提交。`--check` 比对内容 sha256、文件清单与
   `VENDOR.json` 记录的 `git_head`，不一致返回 1——这是 CI 可用的漂移门。
2. **兼容缝只有 `smj/compat.py`**。vendored 代码里对主仓 backends / resources 的 3 处 import
   被 `sync_core.py` 确定性地改写为 `from smj.compat import ...`，每处必须命中一次，否则同步
   直接失败（防止主仓重构后静默漂移）。
3. **包内不复制科学结论**。`report.md` / `manifest.json` 由 run 目录事实生成，不内置基准值；
   示例数据（`tests/data/bdt/`）只用于回归，不能当生产结论。

## 2. 六阶段契约

所有阶段都实现 `stages/base.py` 的同一接口：`prepare(ctx, stage_dir) -> Prepared`、
`collect(ctx, stage_dir) -> Collected`。`Prepared` 声明 `inputs/resources/script/artifacts/
notes/commands`，`Collected` 声明 `artifacts/results/notes`。

| 阶段 | 输入（默认来源） | 产物 | 执行位置 | 资源模型 |
|---|---|---|---|---|
| `gaussian_opt` | `molecule.xyz`（配置或 `--xyz`） | `mol.gjf`、`sub_g16.sh`、`input_summary.json`；collect 后 `mol.log`、`opt.xyz`、`gaussian_summary.json` | Slurm（Gaussian 16） | `LINEAR`：`60 s + 10 s/atom`，`4 GiB + 0.05 GiB/atom` |
| `aims_tip_relax` | `gaussian_opt/opt.xyz`（XYZ 格式，或 `--geometry`） | `geometry.in`、`control.in`、`sub_aims.sh`、`input_summary.json`；collect 后 `geometry.in.next_step`、`aims.out` | Slurm（FHI-aims） | `LINEAR`：`180 s + 30 s/atom`，`8 GiB + 0.10 GiB/atom` |
| `au60_build` | `aims_tip_relax/geometry.in.next_step`（或 `--geometry`） | `geometry.in`、`electrode_metadata.json` | **本地**（纯 Python，无提交） | 无（`{local: True}`） |
| `aims_scf` | `au60_build/geometry.in`（或 `--geometry`） | `geometry.in`、`control.in`、`sub_aims.sh`、`input_summary.json`；collect 后 `aims.out`、`mos.aims`、`basis-indices.out` | Slurm（FHI-aims） | `ANCHORS`：132 原子 `4320 s / 13.3 GiB`，幂律外推 |
| `aitranss_te` | `aims_scf/` 目录（或 `--scf-dir`） | `tcontrol`、`sub_te.sh`、`self.energy.in`；collect 后 `TE.dat`、`te.out`、`orbital_match.json`、图 | Slurm（aitranss，**复用 `aims_scf` 远端目录**） | `ANCHORS`：132 原子 `2552 s / 8.0 GiB`，幂律外推；`ntasks=1` |
| `junction_analyze` | `aitranss_te/TE.dat`（或 `--te`） | `results.json`（`log10_G_over_G0_like`、`S_uV_per_K`、`ZT_electronic` 等） | **本地**（Landauer 积分） | 无（`{local: True}`） |

要点：

- **输入只认文件**。阶段不 import 上游阶段，也不读上游内存对象；`RunLayout` 给出默认
  路径（`gaussian_xyz()`、`tip_geometry()`、`electrode_geometry()`、`scf_dir()`、`te_dat()`），
  覆盖参数优先于默认路径，覆盖文件不存在直接 `StageBlocked`（fail-closed，不猜）。
- **本地阶段无需提交**。`au60_build`、`junction_analyze` 的 `needs_submit=False`；`--run`
  会 prepare 后直接 collect，`--submit` 会明确报「本地阶段，无需提交」。
- **aitranss 的目录复用是刻意的**。`tcontrol` 里写的是相对 species 路径，aitranss 又直接读
  SCF 的 `mos.aims`/`aims.out`，所以它**在 `aims_scf` 的远端目录里执行**，只上传
  `tcontrol`、`sub_te.sh`、`self.energy.in` 三个文件；不要把它挪进独立目录。

## 3. 状态机

```text
pending ──prepare──▶ prepared ──submit──▶ submitted ──wait──▶ waiting ──wait──▶ finished ──collect──▶ collected
   ▲                    │                    │                    │                │
   │                    └──── dry-run 到此为止 │                    │                └── force 重跑 ──▶ prepared
   │                                          │                    └── 连续两次缺席 ──▶ failed
   │                                          └── 本地模式同步执行完 ──▶ finished/failed（无 waiting）
   └──────────────  force 重来（从 finished/collected 回 prepared）──────────────────┘
```

| 状态 | 谁写入 | 含义 |
|---|---|---|
| `pending` | `state.py` 初始化 | 尚未生成输入 |
| `prepared` | `runner.prepare` | 输入/脚本/资源已落盘；dry-run 的终点 |
| `submitted` | `runner.submit` | 远端已有 job id；本地模式表示脚本已同步跑完 |
| `waiting` | `runner.wait` | 轮询中；**超时也留在这里**，不判失败 |
| `finished` | `runner.wait` | 哨兵 `.dsh_done` 存在，退出码已记录 |
| `failed` | `runner.wait` / `run_stage` 捕获 `SmjError` | 退出码非 0，或连续两次队列缺席，或阶段异常 |
| `collected` | `runner.collect` | 产物已回取并解析，`results` 已写入 `state.json` |
| `skipped` | 预留状态（`RunState.is_done` 认可） | 供外部编排/手工标记；当前 CLI 跳过已完成阶段时**不改**状态 |

配套事实：

- `state.json` 在 run 根（`<run.dir>/<run.name>/state.json`），原子写（临时文件 + `os.replace`）。
- `prepare` 会清空上一轮的 `job_id/remote_dir/时间戳/exit_code/note/error`，避免「新输入配旧作业」。
- 远端重提交前必清 `.dsh_done`、`slurm-*.out` 与本阶段旧产物（`mol.log`/`aims.out`/`TE.dat`/`te.out`）；
  **清理失败直接拒绝提交**——陈旧哨兵会让轮询把上一轮当成这一轮完成。
- `finished` 与 `collected` 是两件事：作业结束不等于产物已解析；`smj run` 遇到 `finished`
  只补 `collect`，不重跑。

## 4. 逐阶段独立运行矩阵

「独立运行」= 不依赖 `smj run` 的顺序推进，只要求该阶段自己的输入存在。

| 阶段 | 单独跑的最小上游 | 覆盖方式 | 示例命令 |
|---|---|---|---|
| `gaussian_opt` | 一个 `.xyz` | `--xyz <file>`（默认 `molecule.xyz`） | `smj stage gaussian_opt --prepare --xyz mol.xyz` → `--run` |
| `aims_tip_relax` | Gaussian 优化后的 `opt.xyz` + species 目录 | `--geometry <file>`（默认上一阶段 `opt.xyz`；**XYZ 格式**，不是 `geometry.in`） | `smj stage aims_tip_relax --prepare --geometry gaussian_opt/opt.xyz` → `--run` |
| `au60_build` | tip 弛豫后的 `geometry.in.next_step` | `--geometry <file>` | `smj stage au60_build --run --geometry aims_tip_relax/geometry.in.next_step` |
| `aims_scf` | 带电极的 `geometry.in` + species 目录 | `--geometry <file>`（默认 `au60_build/geometry.in`） | `smj stage aims_scf --prepare --geometry au60_build/geometry.in` → `--run` |
| `aitranss_te` | 一份 SCF 目录（`aims.out`、`mos.aims`、`geometry.in`、`control.in`）**加** `au60_build/electrode_metadata.json` | `--scf-dir <dir>`、`--metadata <json>` | `smj stage aitranss_te --prepare --scf-dir /path/to/scf --metadata electrode_metadata.json` |
| `junction_analyze` | 一个 `TE.dat` | `--te <file>`（默认 `aitranss_te/TE.dat`） | `smj stage junction_analyze --run --te TE.dat` |

分步控制（每个远端阶段都支持）：

| 动作 | 含义 | 典型用途 |
|---|---|---|
| `--prepare` | 只生成输入/作业脚本，不提交 | 人工检查 `control.in`、`sub_*.sh`、资源申请 |
| `--submit` | 清理残留 → 上传 → `sbatch`，记 job id | 集群排队策略与 prepare 分离 |
| `--wait` | 轮询到 finished/failed；超时不判失败 | 长作业可先退出、稍后回来再 `--wait` |
| `--collect` | 回取产物并解析（可单独执行） | 手工跑完/从别处拷回产物后只做解析 |
| `--run` | prepare → submit → wait → collect（默认） | 一把梭；本地阶段 = prepare → collect |

两条实用路径：

1. **只做输入装配**：`smj stage <stage> --prepare --dry-run`（`run.mode: dry-run` 下等价），
   产出可直接拷到集群手工提交的脚本与输入。
2. **只做结果解析**：把外部跑出的 `TE.dat` 放到任意路径，
   `smj stage junction_analyze --run --te /path/TE.dat`，不碰上游状态机。

注意：覆盖参数只改变**本阶段读哪个文件**，不改上游 `state.json`；独立跑出来的产物
不会自动把上游阶段标成 `collected`。

## 5. 与主仓引擎语义的一致性

独立包不是重写引擎，而是把主仓 `tools/submit.py` 里 SMJ 的提交语义**复刻到 `Transport`
四原语上**（`run/put/get/mkdirs`），因此可以换成本地 shell、SSH 或 `FakeTransport`。

| 语义 | 独立包实现 | 与主仓一致的点 |
|---|---|---|
| 完成判定 | 远端 `.dsh_done` 首字段为退出码 | 哨兵文件名/格式一致；Slurm 的 shell 退出码不可靠，不拿它当判据 |
| 队列缺席 | 仅当 `squeue` **成功**且连续两次看不到 job id 才判 `failed` | 查询失败（非零退出/断连）不计缺席，避免网络抖动误杀 |
| 超时 | `finished=False, failed=False`，状态留在 `waiting` | 作业可能仍在跑；可稍后 `--wait`，不产生假失败 |
| 残留清理 | 提交前 `rm -f .dsh_done slurm-*.out <本阶段旧产物>`，失败即拒绝提交 | 与主仓 `SMJ_STALE_UPLOADS` 同表 |
| 上传过滤 | 跳过旧回取产物（`mol.log`/`aims.out`/`TE.dat`/`te.out`/`.dsh_done`） | 不把上一轮结果当输入传回去 |
| aitranss 目录 | 在 `aims_scf` 远端目录执行，只传 3 个 TE 输入 | 相对 species 路径与 SCF 产物依赖 |
| job id | 从 `sbatch` stdout/stderr 解析 `Submitted batch job <id>` | 解析不到就报错并要求先 `squeue` 确认，不盲重提 |
| 环境/资源 | `cluster.env` 覆盖 → 阶段默认值；`G16ROOT`/`AIMS_BIN`/`AITRANSS_BIN`/`MKL_LIB` 等键与 vendored `ENV_SPECS` 同源 | 缺省值不静默成功：渲染器写 `UNCONFIGURED` 占位，脚本在集群上会失败；`preflight` 非 dry-run 时逐项检查二进制/许可 |
| 判据 | 工程判据（`engineering_passed`）与物理判据（`physics_verdict`）分离 | `report.md` 是机械单一来源；电子 ZT ≠ 总 ZT |

`FakeTransport`（`smj.testing`）在内存里记录 `run/put/get/mkdirs` 并可用
`simulate_done(remote_dir, rc)` 造哨兵；`sbatch` 输出与轮询探针（`@@SENT/@@QRC/@@ID`
协议）需要在子类 `on_run` 里模拟，参考 `smj/tests/test_remote_semantics.py` 的
`_probe`。它验证状态机与命令序列，不模拟真实调度器。

## 6. 资源估算数值来源

`smj/resources.py` 是纯函数模块（无 IO），输出 `walltime/time_limit/mem/ntasks/cpus`，
供 `prepare` 写进作业脚本与 `state.json`。数值不是拍的：

| 模型 | 数值 | 来源/依据 |
|---|---|---|
| `aims_scf` 锚点 | 132 原子 → `4320 s`、峰值 `13.3 GiB` | 主仓真跑 SCF 的墙钟/内存记录（含 `basis-indices.out` 等尾部写入） |
| `aitranss_te` 锚点 | 132 原子 → `2552 s`、峰值 `8.0 GiB` | 主仓真跑 TE 记录；`-c 16` 比单核慢约 10%（2026-11-05 A/B），故 `ntasks=1` |
| 线性模型 | Gaussian `60 s + 10 s/atom`、tip 弛豫 `180 s + 30 s/atom` | 小分子/中等体系的保守拟合，按原子数线性放大 |
| 安全系数 | 墙钟 ×1.25、内存 ×1.20；`time_limit = 1.9 ×` 估算墙钟 | SCF 收敛后仍会写产物约 25 分钟，`time_limit` 必须留尾 |
| 下限 | Gaussian/tip `900 s`，SCF/aitranss `1800 s` | 避免小体系申请过短被调度器秒杀 |
| 上限 | 墙钟 `72 h`、内存 `128 GiB`；超限标 `capped` 并给 `check_smj_budget` 原因 | 集群单作业上限；超限应先降体系或拆分 |
| 并行度 | `gaussian_opt=1`、`aims_tip_relax=4`、`aims_scf=4`、`aitranss_te=1` | 与主仓阶段默认一致；aitranss 不做 MPI |

外推是**幂律**而非线性（`aims_scf`：`wall ∝ n^3`、`mem ∝ n^2`；`aitranss_te`：
`wall ∝ n^3`、`mem ∝ n^1.5`），132 原子是锚点；体系明显偏离锚点时
估算只作申请下限，真跑后应把实测值回填到 `state.json`/台账，而不是继续信外推。

## 7. 风险与未验证项

| 项 | 影响 | 现状/缓解 |
|---|---|---|
| 二进制与许可证不随包分发 | 无 Gaussian/FHI-aims/aitranss 时只能 dry-run 或本地阶段 | 文档明确「不自带 DFT 软件」；`preflight` 逐项检查可执行文件/许可证 |
| species 基组不在包内 | `aims_*` 缺 species 直接失败 | `paths.species_defaults` / `cluster.species_dir` 指向集群既有目录；preflight 检查元素覆盖 |
| Windows 本机执行需 bash 与 ssh | `run.mode: local` 无 bash 不可用；`slurm-ssh` 需可用的 ssh/scp 与免密登录 | `LocalTransport` 明确报错并提示装 Git-Bash/WSL 或改 `dry-run`；真跑建议在宿主 shell |
| 真实集群未在本次交付中验收 | 提交/轮询/回取只在 `FakeTransport` 与本地 fixture 上验证过 | 首次真跑按 §8 路径：先 `smj preflight`，再小体系 `smj run --from gaussian_opt --to au60_build` 分段验收 |
| Au60 与 20 Au/侧不是同一口径 | 绝对值不可直接比较 | 文档与判据表标注口径；不混用基准 |
| Gaussian/FHI-aims/aitranss 版本漂移 | 输入模板/解析可能不兼容 | `preflight` 记录版本；模板随 vendored 同步，不手改 |
| 超时/排队策略 | 长作业可能跨天 | `--wait` 可重复执行；超时不算失败，状态留 `waiting` |
| `dry-run` 只到 `prepared` | 不能产出整链假物理 | 刻意 fail-closed；演示/CI 只验证输入装配 |
| 六阶段适配器 | 已落地；离线用例覆盖「预置上游产物后单阶段 prepare/collect」与缺上游 `StageBlocked` | 本次未重跑真实 Gaussian/FHI-aims/aitranss；首次真跑按 §8 分段验收 |
| FakeTransport 不是集群 | 离线演练不覆盖真实调度器行为 | 只验证提交/轮询/回取的状态机；真跑仍需集群验收 |
| 物理判据依赖参考值 | 新分子没有基准时 `physics_verdict` 可能不适用 | 判据表区分「不适用/不确定/失败」，不把缺基准当通过 |

## 8. 最小验收路径（研究者视角）

```bash
python -m pip install -e smj                 # 或 python -m build smj && pip install smj/dist/*.whl

mkdir smj-demo && cd smj-demo                # molecule.xyz 放这里
smj init . --name demo --xyz molecule.xyz --mode dry-run
smj plan                                     # 六阶段计划、依赖、当前状态
smj stage gaussian_opt --prepare             # 只生成输入
smj run --dry-run                            # 干跑到缺上游产物处停住（预期边界，退出码 1）
cd .. && python smj/scripts/smoke_install.py # 打包/导入/vendored/同步门 冒烟（在仓库根执行）
```

真跑时把 `run.mode` 改成 `slurm-ssh`，补 `cluster.host/remote_root/partition` 与
`cluster.env`（`G16ROOT`、`AIMS_BIN`、`AITRANSS_BIN`、`MKL_LIB` 等），先 `smj preflight`
再 `smj run`。判读看 `report.md` 的「概览/阶段状态/输运结果/失败原因」，不要只看退出码。
