# SIESTA 零偏压路线验证记录（2026-11-27）

> 本记录回答一个问题：这条路线能否在真实集群上，从分子 XYZ 一路跑到可复核的透射谱、电导和报告。
>
> **结论先行：工程链路已贯通，两套体系都完成了六阶段并生成了可追溯的 `report.md`；但物理量仍是 `inconclusive`，不能当作材料电导结论，也不能用这两例比较分子趋势。**

## 1. 环境与方法

| 项 | 设置 |
|---|---|
| 代码状态 | 0.2.2rc2 工作区；SIESTA 路线当时尚未提交，本记录对应 2026-11-27 的工作区快照 |
| 集群 | 入口 gpu1；Slurm 侧 `cluster.host: gpu-outside`，partition `main`（站点上限 7 天），作业配置 `time_limit: 20:00:00`，4 个 MPI rank + 1 GPU，32 GB 内存 |
| 程序 | 模块 `siesta/5.4.2-gpu-netcdf`（SIESTA 5.4.2 / TranSIESTA / TBtrans，NetCDF 输出） |
| 泛函与基组 | PBE GGA；C、H、N、Au 均为 DZP；元素级 PSML 赝势 |
| 网格与展宽 | `MeshCutoff 300 Ry`；`ElectronicTemperature 0.025 eV`；`SCF.DM.Tolerance 1e-4`；最大 300 次 SCF |
| 电极 | Au(111) `3x3`，6 个 principal layer，8 个 slab layer，1 个可动 Au 层；电极 k 点 `3x3x12` |
| 器件 k 点 | device `3x3x1`；TBtrans `5x5x1` |
| 输运 | 温度 300 K；能量窗 `[-1.0, +1.0]` eV，步长 5 meV，共 401 点 |
| 接触 | 两侧 `flat/top`；bpy 初始高度 2.35 Å，ppd 2.55 Å；ppd 复用 bpy 的电极阶段 |
| 色散修正 | 当前 FDF 模板未写入 vdW/色散项；定量比较接触距离或结合能时需另行评估 |

`MeshCutoff`、`ElectronicTemperature`、`SCF.DM.Tolerance`、`MD.MaxForceTol 0.02 eV/Ang` 等来自 `src/smj/siesta/stages.py` 的 FDF 生成逻辑；其余来自两份运行配置。

## 2. 两套体系

| 运行名 | 分子 | 锚定 | 分子原子数 | 散射区原子数 | 接触高度 | 电极来源 |
|---|---|---|---:|---:|---:|---|
| `bpy-siesta-flat` | family `BPY` | N / N | 20 | 164 | 2.35 Å | 新建 `3x3x6` 电极 |
| `ppd-siesta-flat` | family `BDA` | N / N | 16 | 160 | 2.55 Å | 复用 `bpy-siesta-flat/siesta_electrode` |

两套体系都从运行目录的 `runs/<分子>/molecule.xyz` 出发，按 `strip_thiol_h: false` 处理（锚点是 N，不是硫醇 S），自旋 `unpolarized`，器件混合权重 0.015。

## 3. 六阶段工程证据

两例的 `state.json` 都显示六个阶段 `collected`、各集群阶段作业退出码 0。作业号可供追溯：

| 阶段 | bpy | ppd |
|---|---|---|
| `molecule_opt` | 26958 | 26959 |
| `siesta_electrode` | 26960 | 复用 bpy，无新作业 |
| `siesta_junction_relax` | 26963 | 26964 |
| `transiesta_scf` | 26965 | 26966 |
| `tbtrans_te` | 26967 | 26968 |
| `junction_analyze` | 本地 | 本地 |

逐阶段核对结果：

| 检查项 | bpy | ppd |
|---|---|---|
| 分子优化 SCF / 几何 | 收敛；4 次 SCF；最大力 0.012088 eV/Å | 收敛；4 次 SCF；最大力 0.017451 eV/Å |
| 电极 SCF / 主层检查 | 9 次 SCF；`principal_layer_checked=true` | 复用同一电极 |
| 结弛豫最大力 | 0.015574 eV/Å（阈值 0.02） | 0.018645 eV/Å（阈值 0.02） |
| 接触 | 左/右 `top`，配位 1，Au–N 2.3503 Å | 左/右 `top`，配位 1，Au–N 2.7612 Å |
| 电极匹配 | 左 `−a3`、右 `+a3`，`matched=true` | 同左 |
| 电极闭合误差 | 3.83×10⁻¹⁰ Å | 同左 |
| 器件 SCF | `open_boundary_converged=true`；16 次；Δq = −0.005790 e | 23 次；Δq = −0.007490 e |
| 左右化学势 | 相等 | 相等 |
| TBtrans 谱 | 401 点，`transmission.tsv` 为有限值 | 同左 |
| 屏蔽标记 | `unverified` | `unverified` |

“收敛”在这里只表示该阶段自己的判据通过（SCF 残差、力阈值、开放边界标志），**不表示**横向 k 点、电极厚度或屏蔽已经收敛。

## 4. 输运结果

| 运行名 | E_F（eV） | G/G₀ | G（S） |
|---|---:|---:|---:|
| `bpy-siesta-flat` | 0.1981398615 | 1.814203377×10⁻³ | 1.405661418×10⁻⁷ |
| `ppd-siesta-flat` | 0.1793218747 | 0.4439993888 | 3.440147992×10⁻⁵ |

两例都在 300 K、0.025 eV 展宽下按 Landauer 积分，费米窗积分覆盖率为 1.0（权重积分 0.9999999999999999），G 按 G₀ = 2e²/h 归一。数值与 `report.md`、`manifest.json` 一致。

bpy 与 ppd 相差约 245 倍，但两者分子不同、接触高度不同、ppd 还复用了电极，**不能把这一比值解释成化学趋势**。

## 5. 判定与边界

### 5.1 机械判定

- `smj report`：两例均为 `engineering_passed=true`、`physics_verdict=inconclusive`、`errors=0`、输运警告 0。
- `smj review`：两例均为 `status=warning`、`accepted=false`、`errors=0`；警告原文是“物理判定 inconclusive：横向采样、电极厚度与屏蔽的收敛证据待核实”。

`inconclusive` 不是遗漏，而是本次验证时的显式选择：当时 `src/smj/siesta/report.py` 在写出物理判定时固定为 `inconclusive`，避免工程通过被误读成物理通过。此后代码已改为由 `screening` 与 `kconv` 两项证据共同决定判定；本记录保留的是 3x3/ppd 那次运行的原始结果，不把后续代码能力倒填进旧数据。

### 5.2 收敛与对照状态

- 横向 k 点：本记录中的器件只用了 `3x3x1`，TBtrans 只用了 `5x5x1`。后续 [k 点验证记录](validation-siesta-kconv-2026-11-28.md) 完成了 bpy 的 3×3→5×5（+22.853%）与 5×5→7×7（+5.646%），两者都超过 1% 判据；ppd 的 3×3→5×5 为 −0.726%，通过。因此 bpy 的 3×3 数值不收敛，本记录中的 bpy `G/G0=1.814×10⁻³` 只能作为旧网格结果，不能当收敛值；ppd 的 3×3 结果在 1% 内稳定。
- 电极厚度：6 个 principal layer 是固定设置，没有做层数扫描。
- 屏蔽：本记录中的 3×3/ppd 运行时，`screening` 字段仍写死为 `unverified`；后续 5×5 运行已从逐层 Mulliken 电荷得到 `passed`（最外两层 |Δq| 在 10⁻³ e 以下，电极体相层 spread 约 10⁻⁵ e），旧运行不追溯获得该证据。
- 色散：FDF 未启用 vdW，接触距离/结合能的定量解释需另行评估。
- 偏压与温度：只做了零偏压、300 K 单点，没有有限偏压或温度依赖。
- 外部标定：没有和实验、其他代码或更高精度设置对账。

### 5.3 代码测试状态

- SIESTA 相关测试共 36 个测试函数（48 个用例）。本机 Windows 沙箱下 **29 passed、19 setup errors**；19 个错误全部是 pytest 的 `tmp_path` 在 `.tmp/pytest-of-whao` 下建目录时被 ACL 拒绝（`PermissionError: [WinError 5]`），不是断言失败。
- 完整套件同一时点为 **1040 passed、10 failed、19 errors**；10 个失败中 9 个是 Git Bash 在沙箱内无法创建 signal pipe（`Win32 error 5`），1 个是 `test_cluster_paths` 把 Windows `HOME` 传给 git bash。均定位到本机环境，但未在 Linux/CI 复跑。
- 因此本记录的依据是真实 run 链路与 `report.md`/`review.json`，不是“测试全绿”。

## 6. 复现

运行目录 `runs/siesta-flat-pair/` 不入库，下面是当时的命令与产物位置：

```bash
smj preflight -c runs/siesta-flat-pair/bpy.yaml
smj run       -c runs/siesta-flat-pair/bpy.yaml
smj report    -c runs/siesta-flat-pair/bpy.yaml
smj review    -c runs/siesta-flat-pair/bpy.yaml

# ppd 同理；它的 siesta.electrode_reuse 指向 bpy-siesta-flat/siesta_electrode
smj run       -c runs/siesta-flat-pair/ppd.yaml
```

关键产物：

- `runs/siesta-flat-pair/<run>/manifest.json`：协议、阶段、产物 sha256、判定。
- `runs/siesta-flat-pair/<run>/report.md`：人读报告。
- `runs/siesta-flat-pair/<run>/review.json`：review 结论与警告。
- `runs/siesta-flat-pair/<run>/tbtrans_te/transmission.tsv`：标准化透射谱。
- `runs/siesta-flat-pair/<run>/junction_analyze/results.json`：G、G/G₀、E_F 与自旋约定。

复现需要同一套 `siesta/5.4.2-gpu-netcdf` 模块、可读的 C/H/N/Au PSML 赝势，以及配置中声明的 Slurm 资源；换模块或赝势后，电极复用检查会拒绝旧电极。
