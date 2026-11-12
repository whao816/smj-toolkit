# smj — 单分子结输运计算链（独立包）

> 从**分子结构**到**透射谱**的端到端自动化：Gaussian 优化 → FHI-aims 加针尖/弛豫 →
> 本地构建 Au(111) 电极 → FHI-aims 静态 SCF → aitranss 透射 → Landauer 分析。
> 整链可一把跑完，也可**逐阶段独立运行**；每个阶段的输入/产物都落在
> `<run.dir>/<run.name>/<stage>/`，状态机在 `state.json` 里可查、可恢复。
>
> 本包是 AgenticDFT 主仓 SMJ 链的独立打包（发行名 `smj-toolkit`，
> import 名 `smj`，命令行 `smj`）。**主仓仍是权威源**：算法以
> `agenticdft/agenticdft/tools/junction/` 为准，本包内是只读副本，由
> `scripts/sync_core.py --check` 守漂移。

## 0. 定位与边界

| | |
|---|---|
| 输入 | 分子 `.xyz`（默认硫醇 -SH 锚定）+ 电极/锚点/输运窗口参数（`smj.yaml`） |
| 输出 | `TE.dat`（透射谱）、`results.json`（`log10_G_over_G0_like`、`S_uV_per_K`、`k_e_pW_per_K`、`ZT_electronic`）；`smj report` 另生成 `manifest.json`、`report.md` 与判定（`manifest.verdict`） |
| 适用 | 单分子结 / 分子电子学：Au(111) 电极间硫醇锚定分子的电输运与热电 |
| 不承诺 | 自带 Gaussian/FHI-aims/aitranss；自动替人选方法/基组/锚点；替代人工结构确认；把 20 Au/侧与 Au60 口径混比 |

七阶段链（顺序固定，但可单独执行任意阶段）：

```text
molecule.xyz
  └─ gaussian_opt      Gaussian 16 opt+freq            → opt.xyz
      └─ aims_tip_relax  FHI-aims 加 Au 针尖 + 弛豫      → geometry.in.next_step
          └─ au60_build  本地贴装 Au(111) 电极（库团簇）  → geometry.in + electrode_metadata.json
              └─ au_relax    FHI-aims 背层固定弛豫        → 弛豫后的 geometry.in（背层冻结）
                  └─ aims_scf     FHI-aims 静态 SCF      → aims.out / mos.aims
                      └─ aitranss_te  aitranss 透射      → TE.dat
                          └─ junction_analyze  Landauer  → results.json
```

## 1. 安装

要求 **Python ≥ 3.10**；运行时依赖只有 **PyYAML** 与 **numpy**。

```bash
# 源码可编辑安装（在仓库根执行；pyproject.toml 在根，源码在 src/）
python -m pip install -e .

# 或构建 wheel 后安装（产物在 dist/）
python scripts/gen_vendor_json.py           # 先刷新 vendored 补丁台账（VENDOR.json 随 wheel 分发）
python scripts/build_wheel.py --check       # 只用标准库，不需要 setuptools/build；--check 按 RECORD 校验
python -m pip install dist/smj_toolkit-0.1.1-py3-none-any.whl
```

可选依赖：

```bash
python -m pip install -e ".[plot]"   # matplotlib，用于 T(E) 出图与 A2 核对图
python -m pip install -e ".[dev]"    # pytest，跑本包测试
python -m pip install -e ".[all]"    # plot + dev
```

安装后验证（在仓库根执行）：

```bash
smj --help
smj verify      # Python/PyYAML/numpy、vendored junction、模板、7 个阶段模块
```

## 2. 5 分钟快速开始

```bash
# 0) 建一个工作目录，把 molecule.xyz 放进去（init 只写配置，不复制/生成结构）
mkdir smj-demo && cd smj-demo

# 1) 生成配置骨架（默认 dry-run，不提交任何作业）
smj init . --name demo --xyz molecule.xyz --mode dry-run

# 2) 编辑 ./smj.yaml：分子路径、电极参数；真跑时补 cluster 段
#    （见 §4；只想先看输入装配的话保持 run.mode: dry-run 即可）
#    注意：0.1.2 默认 electrode.source: library，需要 FHI-aims 的 electrodes.library
#    （用 electrode.library_dir 指向它，或写 source: builtin 回退到自建电极）

# 3) 看七阶段计划、依赖与当前状态
smj plan

# 4) 只生成第一阶段输入（人工检查 mol.gjf / sub_g16.sh）
smj stage gaussian_opt --prepare

# 5) 干跑整链：逐阶段 prepare，不提交（预期停在缺上游产物的下一阶段，退出码 1）
smj run --dry-run
```

配置在 `./smj.yaml`，run 目录是 `./demo/`（`run.dir=.`、`run.name=demo`，
即 `<run.dir>/<run.name>/`）。`dry-run` 的终点是 `prepared`：
它**不会**伪造 `mol.log`/`opt.xyz` 让整链看起来跑通，所以整链干跑会停在需要上游产物的
下一阶段（该阶段标 `failed`，阻塞原因写进 note，`smj run` 退出码 1）——这是刻意的
fail-closed 语义，不是缺陷。

真跑把 `run.mode` 改成 `slurm-ssh`，补齐 `cluster.host/remote_root/partition` 与
`cluster.env`，然后（仍在 `smj-demo/` 下）：

```bash
smj preflight          # 提交前检查
smj run                # 整链
# 或从某阶段续跑 / 跑到某阶段：
smj run --from aims_scf --to aitranss_te
```

## 3. 逐阶段独立运行

阶段之间只通过文件耦合：每个阶段读上游 stage 目录里的约定文件，覆盖参数优先。
因此任意阶段都可以单独 `--prepare` / `--submit` / `--wait` / `--collect`。

下表示例省略 `-c <path>`（默认读当前目录的 `smj.yaml`）；在别处执行请加
`--config my-run/smj.yaml`。`--run` 在 `dry-run` 模式下只 prepare（见 §2）。

| 阶段 | 默认输入 | 主要产物 | 需提交？ | 示例命令 |
|---|---|---|---|---|
| `gaussian_opt` | `molecule.xyz`（或 `--xyz`） | `mol.gjf`、`sub_g16.sh`、`input_summary.json`；collect 后 `mol.log`、`opt.xyz`、`gaussian_summary.json` | 是（Gaussian） | `smj stage gaussian_opt --run --xyz molecule.xyz` |
| `aims_tip_relax` | `gaussian_opt/opt.xyz`（或 `--geometry`） | `geometry.in`、`control.in`、`sub_aims.sh`；collect 后 `geometry.in.next_step` | 是（FHI-aims） | `smj stage aims_tip_relax --run --geometry gaussian_opt/opt.xyz` |
| `au60_build` | `aims_tip_relax/geometry.in.next_step`（或 `--geometry`） | `geometry.in`、`electrode_metadata.json`（含背层索引与库团簇来源） | **否**（本地） | `smj stage au60_build --run --geometry aims_tip_relax/geometry.in.next_step` |
| `au_relax` | `au60_build/geometry.in` + metadata | `geometry.in`（两侧背层全部 `constrain_relaxation`）、`control.in`、`sub_aims.sh`；collect 后 `geometry.in.next_step` | 是（FHI-aims） | `smj stage au_relax --run` |
| `aims_scf` | `au_relax/geometry.in`（回退 `au60_build/geometry.in`，或 `--geometry`） | `control.in`、`sub_aims.sh`；collect 后 `aims.out`、`mos.aims` | 是（FHI-aims） | `smj stage aims_scf --run` |
| `aitranss_te` | `aims_scf/`（或 `--scf-dir`） | `tcontrol`、`sub_te.sh`；collect 后 `TE.dat`、`te.out` | 是（aitranss） | `smj stage aitranss_te --run --scf-dir /path/to/scf` |
| `junction_analyze` | `aitranss_te/TE.dat`（或 `--te`） | `results.json` | **否**（本地） | `smj stage junction_analyze --run --te TE.dat` |

FHI-aims 阶段（`aims_tip_relax`/`au_relax`/`aims_scf`）还需要 `cluster.species_dir` 或
`paths.species_defaults` 指向 species 目录；缺失时 `prepare` 直接 `StageBlocked`
（`smj preflight` 也会报 error）。`--geometry` 在三个消费几何的阶段
（`aims_tip_relax`/`au60_build`/`aims_scf`）都可用，但语义不同：`aims_tip_relax`
要**已优化分子的 XYZ**（`opt.xyz`），不是 FHI-aims 的 `geometry.in`；`au60_build`
与 `aims_scf` 要的才是 FHI-aims 几何。`aitranss_te` 除 SCF 目录外还要
`au60_build/electrode_metadata.json`（六个表面索引与左右 tip 索引），可用
`--metadata` 覆盖。

`au60_build` 默认用 FHI-aims 自带的 `electrodes.library` 预定义团簇贴装电极（手册
§5.3.1 的要求），由 `electrode.source`、`electrode.library_dir`、`electrode.library_cluster`
控制（默认 `fcc111.au.clusters.37_atoms/fcc111.cluster.natoms37.x1.aims`）；库文件不进
wheel，只在 metadata 里记 `sha256` 与相对路径。`electrode.nlayers: null` 表示取库文件头部
建议的层数。`electrode.source: builtin` 回到旧的自建三角岛。`electrode.relax: true`
（默认）让 `au_relax` 在贴装后冻住两侧背层整层、做一次小力阈值弛豫（手册 §5.3.1 的
"固定两个背层做一次弛豫"）；置 `false` 时该阶段本地直通（写 `relax_skipped.json`）。

`electrode.relax_freeze: electrodes` 则冻结全部 Au、只让分子弛豫——理想晶格电极本来
就已弛豫好，把 60+ 个 Au 一起交给 BFGS 会让最大力在 0.4–0.9 eV/Å 之间来回；`backplane`
（默认）只冻两侧背层，适用于电极本身还要微调的场合。

`transport.energy_window: [ener, estep, eend]` 的单位是 **Hartree（绝对能量）**，与
`tcontrol` 的 `$ener/$estep/$eend` 一致，要围绕 `aims_scf/aims.out` 的 Fermi 能级取值
（Au 体系约 −0.19 Ha 附近取 ±0.02 Ha）。`TE.dat` 第一列是绝对 Hartree，第二列才是
`E−E_F`（eV）。

`molecule.family` 声明判据标定族（`AZ` = Au-(2,6-azulenedithiol)-Au、`BDT` 等）；
未登记体系按"全部判据适用"处理（fail-safe）。BDT 的绿带（T(E_F)/G ∈ [0.03, 0.40]、
Seebeck 为正）不该给薁作闸门，声明 `family: AZ` 后这三条记"不适用"，判定不再误判 failed。

动作开关（每个阶段通用）：

| 开关 | 含义 | 什么时候用 |
|---|---|---|
| `--prepare` | 只生成输入/作业脚本，不提交 | 人工检查 `control.in`、资源申请 |
| `--submit` | 清理残留 → 上传 → `sbatch`，记录 job id | 想让排队与准备分离 |
| `--wait` | 轮询到结束；超时不判失败 | 长作业先退出，稍后回来再等 |
| `--collect` | 回取产物并解析（可单独执行） | 手工跑完/从别处拷回产物后只解析 |
| `--run` | prepare → submit → wait → collect（默认） | 一把梭；本地阶段 = prepare → collect |

覆盖参数只改变**本阶段读哪个文件**，不会改动上游 `state.json`；独立跑出的产物不会自动
把上游阶段标成 `collected`。常用组合：

```bash
# 只做输入装配，拷到集群手工提交
smj stage aims_scf --prepare --geometry /data/run/au60_build/geometry.in --config smj.yaml

# 外部已有 TE.dat，只做分析（--te 也可给 aitranss_te 输出目录）
smj stage junction_analyze --run --te /data/te/TE.dat --config smj.yaml
```

## 4. Slurm 配置

`cluster` 段只在 `run.mode: slurm-ssh` 时必需；`dry-run` / `local` 下可留空。
`slurm-ssh` 至少要求 `host`、`remote_root`、`partition`，缺一直接报配置错误。

```yaml
run:
  name: bdt-trans
  dir: .                 # run 根目录；最终产物在 <dir>/<name>/
  mode: slurm-ssh        # dry-run | local | slurm-ssh
  clean: false           # true 时 prepare 前清空本阶段目录（仅限 run 目录内）
  keep_going: false      # true 时整链遇到失败继续跑后续阶段

cluster:
  host: cluster.example.edu
  user: yourname
  port: 22
  key: ~/.ssh/id_ed25519
  remote_root: /data/yourname/smj
  partition: main
  ntasks: 4
  gres: ""
  account: ""
  qos: ""
  species_dir: /data/software/fhi-aims/species_defaults/defaults_2020/light
  time_limit: ""
  mem: ""
  poll_seconds: 60
  timeout_hours: 6
  env:
    g16root: /data/software/gaussian/g16
    aims_bin: /data/software/fhi-aims/bin/aims.221103.scalapack.mpi.x
    aitranss_bin: /data/software/aitranss/aitranss
    mkl_lib: /data/software/intel/mkl/lib/intel64
```

`env` 的键与 vendored `inputs.py` 的 `ENV_SPECS` 同源（**区分大小写，请写小写**），
常用项：

| 键 | 用途 | 默认 |
|---|---|---|
| `g16root` / `g16_module` | Gaussian 安装根 / module 名 | 空（`g16root` 或 `g16_module` 二选一） |
| `g16_scratch` / `g16_mem` / `g16_walltime` / `g16_gres` | Gaussian 临时目录/内存/墙钟/GRES | `16G` / `01:00:00` / 空 |
| `aims_prefix` / `aims_module` / `aims_bin` / `aims_ld_module` | FHI-aims 定位方式 | 空（三者择一） |
| `aims_mem` / `aims_walltime` / `aims_gpus` | FHI-aims 资源 | `16G` / `02:00:00` / 空 |
| `aitranss_bin` | aitranss 可执行文件 | 空 |
| `mkl_lib` | aitranss 依赖的 MKL 库路径 | 空 |
| `te_mem` / `te_walltime` / `te_gres` | aitranss 资源 | `8G` / `02:00:00` / 空 |

缺省值不静默成功：输入渲染器在缺 `g16root`/`g16_scratch` 等键时，会把
`UNCONFIGURED` 写进作业脚本（例如 `export g16root=UNCONFIGURED`），这样的脚本在集群上会
直接失败。`smj preflight` 在**非 `dry-run`** 模式下检查 Gaussian 许可、FHI-aims 与
aitranss 可执行文件，缺一项即 error（退出码 1）；`smj run` 不会自动替你跑 preflight，
真跑前请先跑一次。`dry-run` 不触网、也不核对二进制，所以干跑通过不等于集群能跑。
0.1.1 起提交前会再扫一遍作业脚本：仍含 `UNCONFIGURED` 时 `--submit` / `smj run` 直接
拒绝（`StageError`，退出 1），不再让它到集群上才失败。

## 5. 三种运行模式与离线演练

| 模式 | 行为 | 适用 |
|---|---|---|
| `dry-run` | 只 `prepare`，产出输入/脚本；不建传输层、不提交 | 本机装配输入、CI、评审 |
| `local` | 在本机 `bash <script>` 同步执行，哨兵 `.dsh_done` 记退出码 | 有 Gaussian/FHI-aims 的工作站、调试 |
| `slurm-ssh` | 上传 → `sbatch` → 轮询 → 回取 | 真实集群 |

`local` 需要 `bash`（Windows 上装 Git-Bash/WSL 并加入 PATH）；没有 bash 时
`LocalTransport` 会明确报错并建议改用 `dry-run`。

不想连集群也想演练远端链路时，用 `smj.testing.FakeTransport`（内存实现
`run/put/get/mkdirs`，记录每条命令；`simulate_done(remote_dir, rc)` 可造哨兵）：

```python
from smj.config import load_config
from smj.runner import StageRunner
from smj.testing import FakeTransport
from smj.transport import CommandResult

transport = FakeTransport()

def on_run(command, cwd, timeout):
    if "sbatch" in command:                    # 假调度器：给出可解析的 job id
        return CommandResult(0, "Submitted batch job 4242")
    return None                                 # 其余命令按默认成功处理

transport.on_run = on_run
config = load_config("smj.yaml")               # 在 run 配置所在目录执行
# 演练远端链路要求 run.mode 为 slurm-ssh（local 走本机 bash，不演练上传/sbatch）；
# 不想改文件时可在内存里改（下面的 host 等只用于通过字段检查）：
config.run.mode = "slurm-ssh"
config.cluster.host = "fake-cluster"
config.cluster.remote_root = "/tmp/smj-fake"
config.cluster.partition = "main"
runner = StageRunner(config, transport=transport, dry_run=False)
runner.prepare("gaussian_opt")
sub = runner.submit("gaussian_opt")            # 只走内存：清理 → 上传 → 记录 job_id=4242
print(sub["job_id"], sub["cleaned"], len(sub["uploaded"]))
print(transport.commands[-2:])                 # 看实际发出的命令
```

`--wait` 的轮询探针本身是一条 shell 命令，FakeTransport 需要按
`@@SENT/@@QRC/@@ID` 协议模拟输出才能驱动（可参考 `tests/test_state_machine.py`
的 `_SbatchTransport`）；只验证提交/清理逻辑时到 `submit` 为止即可。FakeTransport
不模拟排队与真实程序行为。

## 6. run 目录布局

```text
<run.dir>/<run.name>/            # smj.yaml 的 run.dir 是 run 根，默认 "." → ./<name>/
├── state.json                   # 阶段状态机（status/job_id/时间戳/资源/哈希/results）
├── run_spec.json                # 运行参数快照（完整解析配置、协议指纹、软件版本）
├── run_spec.previous.json       # 仅 --force 重绑时保留的上一份快照
├── manifest.json                # 产物清单（sha256、大小、阶段、回取模式）
├── report.md                    # 人读报告：概览/阶段状态/输入摘要/资源/输运/失败/产物
├── gaussian_opt/                # 每阶段一个目录，输入+产物都在里面
│   ├── mol.gjf  sub_g16.sh  input_summary.json
│   └── mol.log  opt.xyz  gaussian_summary.json
├── aims_tip_relax/
│   ├── geometry.in  control.in  sub_aims.sh
│   └── geometry.in.next_step  aims.out
├── au60_build/
│   └── geometry.in  electrode_metadata.json
├── au_relax/
│   ├── geometry.in  control.in  sub_aims.sh
│   └── aims.out  geometry.in.next_step  electrode_metadata.json
├── aims_scf/
│   ├── control.in  sub_aims.sh  geometry.in
│   └── aims.out  mos.aims  basis-indices.out
├── aitranss_te/
│   ├── tcontrol  sub_te.sh  self.energy.in
│   └── TE.dat  te.out
└── junction_analyze/
    └── results.json
```

`state.json` 是唯一状态来源：`smj status` 读它；`smj report` 由它 + 产物目录生成
`manifest.json`/`report.md`（默认不覆盖已存在的报告，`--force` 才重写）。

## 7. 附带功能

| 命令 | 作用 | 产物/输出 |
|---|---|---|
| `smj preflight` | 提交前检查：分子可读/原子重叠、species 覆盖、建结几何闸门、资源上限；非 `dry-run` 时再加集群/SSH 连通、远端目录可写、二进制/许可 | 结构化 issues（error/warning/info），有 error 退出 1 |
| `smj status` | 读 `state.json`，打印阶段/状态/job id/远端目录/耗时/备注 | 终端表格（`--json` 机器可读） |
| `smj report` | 由 run 目录生成 `manifest.json` + `report.md`（7 节：概览/判定/阶段状态/输入摘要/资源/输运/失败原因/产物清单）；缺 `TE.dat` 时物理判定为「无法判定」而不是通过 | 报告文件；默认跳过已存在文件，`--force` 重写，`--no-checks` 不算判定 |
| `smj analyze` | 对 `aitranss_te/TE.dat` 做 Landauer 分析（`log10_G_over_G0_like`、`S_uV_per_K`、`k_e_pW_per_K`、`ZT_electronic`），写 `junction_analyze/analysis.json`；出图与参考值比较用 Python API（`smj.analysis.analyze_run(..., plot=True, reference=...)`，需 `[plot]` extra） | `analysis.json`（含 `independent_check` 结果） |
| `smj review` | 结构复核门：建结几何 + A2 人工确认（指纹/表面索引偏差）；`--approve --who <名字>` 登记 A2 确认（审计或出图不通过时**拒绝签署**，退出 2），`--no-figure` 只在你已在别处看过图时使用 | `review.json`、`auN_build/a2_confirmation.json`；`blocked` 时退出 1 |
| `smj audit` | 本地一致性检查（不连远端）：state 与阶段目录/产物/job id 对账 | issues 列表；有 error 退出 1 |
| `smj verify` | 安装自检：Python/PyYAML/numpy、vendored junction、模板、阶段模块 | issues 列表；有 error 退出 1（不需要 `-c`） |
| `smj plan` | 打印七阶段计划、依赖与当前状态，不落盘、不提交 | 终端表格 |

退出码：`0` 成功/作业已结束；`1` 失败或配置/输入错误；`2` 未完成（`--wait` 超时、
作业可能仍在跑，或 `run` 有阶段停在 `waiting`），以及 **A2 签署被拒绝**
（审计不过 / 出图失败 / 缺 metadata）。argparse 用法错误也返回 2。

### 7.1 判定与 A2 红线

- `smj report` 默认算物理判据并写进 `manifest.verdict`；缺 `TE.dat` 时登记
  `te_present`（engineering/error）且结论是 `inconclusive`（**证据缺失 ≠ 通过，
  也不等于否证 H0**），`report.md` 的「## 0.1 判定」会列出未过/未评估的判据。
- `review.require_a2` 默认 `true`：没有 `a2_confirmation.json` 或指纹与当前
  `geometry.in` 不一致时，`smj review` 直接 `blocked`。签署用
  `smj review --approve --who <名字>`（A2 核对图落在 `auN_build/a2_surface_viz.png`，
  审计结果落在 `auN_build/a2_surface_audit.json`）。

### 7.2 运行溯源（run_spec / 哈希闭环）

- **第一次 `prepare`** 会在 run 根写 `run_spec.json`：完整解析后的配置（含模板默认值）、
  分子与电极库的 sha256、`protocol_sha256`（只覆盖会改变物理结果的分子/电极/输运子集）、
  SMJ/Python 版本、git commit/dirty、软件路径与本地可读二进制的 sha256。
- **配置漂移默认拒绝**：下一阶段 `prepare` 时若当前配置的 `protocol_sha256` 与
  `run_spec.json` 不同，直接报错；确认要重绑才加 `--force`，旧快照保留为
  `run_spec.previous.json`。
- **阶段哈希**：`state.json` 每阶段记录 `inputs_sha256` / `outputs_sha256` /
  `parent_stage` / `parent_fingerprint`；`au_relax`、`aims_scf`、`aitranss_te`
  在 `prepare` 时校验上游 `geometry.in` 与电极 metadata 的 `surface_fingerprint`，
  不一致直接 `StageBlocked`（不会拿被换过的几何继续算）。
- **报告绑定**：`manifest.json` 带 `run_spec` 与 `provenance`（run_spec 文件 sha256、
  protocol、是否 reconstructed）；`report.md` 概览显示参数快照状态。旧 run 没有
  `run_spec.json` 时允许用当前配置重建，但标记 `reconstructed: true` 并告警。
- **复核门**：`smj review` 对缺 run_spec/重建只告警，对协议漂移记
  `provenance_drift` 并 `blocked`——结果不能按原参数解读时不会静默放行。

## 8. 与主仓的关系

- **主仓是权威源**：`src/smj/_vendor/junction/` 是 `agenticdft/agenticdft/tools/junction/`
  的副本（含模板），算法改动一律先改主仓，再同步；不要随意手改 vendored 文件。
- **本仓库的来历**：由 0.1.0 wheel 还原成源码树（wheel 本身不含 `scripts/`、`tests/`、
  `VENDOR.json`），因此**没有** `sync_core.py` 可用；0.1.1 的本地补丁都登记在
  `src/smj/_vendor/junction/VENDOR.json` 的 `local_patch` 字段里（逐条写明文件、原因、
  过了哪条测试），便于日后回迁主仓。
- 主仓的 VASP 运行根（`agenticdft/work/`）与本包无关；本包的 run 产物落在
  `run.dir/run.name/`（默认是 `smj.yaml` 所在目录下的 `./<name>/`），不写主仓 `results.db`。

## 9. 已知限制

- **不自带 DFT 软件**：Gaussian 16、FHI-aims、aitranss 的二进制、许可证与 species 基组
  都不随包分发，需要集群/工作站已有并按 §4 配好路径；没有它们时只能 `dry-run` 或跑两个
  本地阶段。
- **Windows 本机执行需要 bash**：`run.mode: local` 依赖 `bash`；`slurm-ssh` 还需要可用的
  `ssh`/`scp`（Git-Bash 自带 OpenSSH）与目标集群免密登录。沙箱内网络受限时，请在宿主
  shell 里执行真跑。
- **结构确认是人的责任**：`au60_build` 之后的 A2 表面确认是红线，管线不会自动放行。
  签署是 `smj review --approve --who <名字>`：审计不过或出图失败直接拒绝（退出 2），
  签名绑定当前 `geometry.in` 的 sha256，几何一改旧签名即失效；`review.require_a2`
  默认 `true`，没有有效签名时 `smj review` 为 `blocked`。
- **口径不可混**：Au60 与 20 Au/侧的绝对值不可直接比较；电子 `ZT_e` 不是总 ZT
  （未测晶格热导就不下总 ZT 结论）。
- **超时不是失败**：`--wait` 超时后状态留在 `waiting`，作业可能仍在跑；稍后重跑
  `smj stage <stage> --wait` 即可。
- **资源估算是申请下限**：132 原子锚点外推，体系明显偏离时应按实测回填，不要当保证值。

## 10. 故障排查（最短路径）

| 现象 | 先看 |
|---|---|
| `smj init` 说配置已存在 | 加 `--force`，或换目录 |
| 阶段卡在 `pending` | 先 `--prepare`；`smj status` 看 note |
| `--submit` 报"尚未 prepare" | 先 `smj stage <stage> --prepare` |
| 远端提交失败 | `smj preflight`；看 `cluster.env` 的 `g16root`/`aims_bin`/`aitranss_bin` |
| 轮询一直 `waiting` | 正常排队或超时；`smj stage <stage> --wait` 可重复执行 |
| `failed` 且 exit_code 非 0 | 看 `<stage>/slurm-*.out`、`aims.out`/`mol.log` 尾部 |
| 报告没更新 | `smj report --force`（默认跳过已存在文件） |
| 独立复核与管线不一致 | `smj analyze` 的 `independent_check` 结果，先查单位/口径 |
