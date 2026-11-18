# smj run 报告: 2CN-656-4-thiomethylphenyl

> 由 `smj.report` 从 run 目录（state.json / 阶段子目录 / manifest.json）机械生成；
> 数值结论仍须人工复核，报告不含图与独立科学裁决。

## 0. 概览

| 项 | 值 |
|---|---|
| run_id | `2CN-656-4-thiomethylphenyl` |
| 根目录 | `demo\2CN-656-4-thiomethylphenyl\2CN-656-4-thiomethylphenyl` |
| 模式 | slurm-ssh |
| 参数快照 | run_spec.json（sha256 100ad476936e…）；protocol d5596b57baeb… |
| 阶段数 / 文件总量 | 6 / 0.025 GiB |
| 阶段状态 | molecule_opt=collected, aims_tip_relax=collected, electrode_build=collected, aims_scf=collected, aitranss_te=collected, junction_analyze=collected |
| 首个未完成阶段 | `—` |
| 生成时间 | 2026-11-16T10:27:32+08:00 |

## 0.1 判定

| 项 | 值 |
|---|---|
| 工程判定 | 通过 |
| 物理判定 | 无法判定（inconclusive） |
| 允许动作 | 停; 须补数据或降级声明 |

- 未过的 warn 级判据: `t_ef_resonance_attribution`, `g_in_theory_band`

## 1. 阶段状态

| 阶段 | 状态 | 目录 | 作业 | 退出码 | 耗时 (s) | 说明 |
|---|---|---|---|---|---|---|
| `molecule_opt` | `collected` | 有 | 26047 | 0 | 1587.5 |  |
| `aims_tip_relax` | `collected` | 有 | 26048 | 0 | 1334.7 |  |
| `electrode_build` | `collected` | 有 | — | — | — |  |
| `aims_scf` | `collected` | 有 | 26051 | 0 | 11539.0 |  |
| `aitranss_te` | `collected` | 有 | 26063 | 0 | 312.3 |  |
| `junction_analyze` | `collected` | 有 | — | — | — |  |

## 2. 输入摘要

| 项 | 值 |
|---|---|
| 分子结构 | D:\dsh_projects\AutoSMJ\demo\2CN-656-4-thiomethylphenyl\2CN-656-4-thiomethylphenyl\gaussian_opt\opt.xyz |
| 分子原子数 | 54 |
| 建结阶段 | electrode_build |
| 结原子数 | 128 |
| 电极构建方式 | library |
| 电极元素 | Au |
| 每侧电极原子数 | 37 |
| 输运温度 (K) | 300.0 |
| TE.dat | demo\2CN-656-4-thiomethylphenyl\2CN-656-4-thiomethylphenyl\aitranss_te\TE.dat |
| TE 数据点 | 301 |
| efermi (Ha) | -0.181943661113 |

## 3. 资源估算

| 阶段 | 墙钟 | --time | 内存 | rank/核 | 核时 | 来源 / 预算 |
|---|---|---|---|---|---|---|
| `molecule_opt` | 00:15:00 | 00:29:00 | 9G | 1 | 0.25 | state; 预算内 |
| `aims_tip_relax` | 00:39:00 | 01:14:00 | 17G | 4 | 2.58 | state; 预算内 |
| `aims_scf` | 01:23:00 | 02:36:00 | 16G | 4 | 5.47 | state; 预算内 |
| `aitranss_te` | 00:49:00 | 01:33:00 | 10G | 1 | 0.81 | state; 预算内 |

## 4. 输运结果

| 量 | 值 |
|---|---|
| T(E_F) [G0] | 0.000471 |
| G/G0 = L0 | 0.000477053 |
| log10(G/G0) | -3.321433 |
| S (µV/K) | 2.7963 |
| k_e (pW/K) | 0.2874 |
| ZT_electronic | 0.000302 |
| 数据点 | 301 |
| 来源 | `demo\2CN-656-4-thiomethylphenyl\2CN-656-4-thiomethylphenyl\aitranss_te\TE.dat` |

## 5. 失败原因

- 无失败阶段。

## 6. 产物清单

| 阶段 | 文件数 | 体积 (MiB) | 关键产物 sha256(前12) |
|---|---|---|---|
| `molecule_opt` | 7 | 1.19 | mol.log:55a5766f847f |
| `aims_tip_relax` | 7 | 23.00 | aims.out:fb64cfc78b4f, geometry.in:6c7f2d1ead60 |
| `electrode_build` | 2 | 0.01 | geometry.in:5ca654c64af8 |
| `aims_scf` | 9 | 1.69 | aims.out:f9b8d7cf5a8a, geometry.in:951152c12dd1 |
| `aitranss_te` | 7 | 0.06 | TE.dat:ce273be9d9ee |
| `junction_analyze` | 1 | 0.00 | — |
