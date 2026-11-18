# smj run 报告: 2CN-656-4-thiomethylphenyl-fr

> 由 `smj.report` 从 run 目录（state.json / 阶段子目录 / manifest.json）机械生成；
> 数值结论仍须人工复核，报告不含图与独立科学裁决。

## 0. 概览

| 项 | 值 |
|---|---|
| run_id | `2CN-656-4-thiomethylphenyl-fr` |
| 根目录 | `demo\2CN-656-4-thiomethylphenyl-fullrelax\2CN-656-4-thiomethylphenyl-fr` |
| 模式 | slurm-ssh |
| 参数快照 | run_spec.json（sha256 c1313aa82c49…）；protocol 58a7f3dc7638… |
| 阶段数 / 文件总量 | 6 / 0.021 GiB |
| 阶段状态 | molecule_opt=collected, aims_tip_relax=collected, electrode_build=collected, aims_scf=collected, aitranss_te=collected, junction_analyze=collected |
| 首个未完成阶段 | `—` |
| 生成时间 | 2026-11-16T10:27:33+08:00 |

## 0.1 判定

| 项 | 值 |
|---|---|
| 工程判定 | 通过 |
| 物理判定 | 否证 H0（failed） |
| 允许动作 | 停; 不得进入预测模式(任务卡 §7.3②) |

- 否证 H0 的 error 级判据: `seebeck_sign`
- 未过的 warn 级判据: `g_in_theory_band`

## 1. 阶段状态

| 阶段 | 状态 | 目录 | 作业 | 退出码 | 耗时 (s) | 说明 |
|---|---|---|---|---|---|---|
| `molecule_opt` | `collected` | 有 | 26065 | 0 | 1549.0 |  |
| `aims_tip_relax` | `collected` | 有 | 26067 | 0 | 1088.2 |  |
| `electrode_build` | `collected` | 有 | — | — | — |  |
| `aims_scf` | `collected` | 有 | 26068 | 0 | 3186.1 |  |
| `aitranss_te` | `collected` | 有 | 26070 | 0 | 279.7 |  |
| `junction_analyze` | `collected` | 有 | — | — | — |  |

## 2. 输入摘要

| 项 | 值 |
|---|---|
| 分子结构 | D:\dsh_projects\AutoSMJ\demo\2CN-656-4-thiomethylphenyl-fullrelax\2CN-656-4-thiomethylphenyl-fr\gaussian_opt\opt.xyz |
| 分子原子数 | 54 |
| 建结阶段 | electrode_build |
| 结原子数 | 128 |
| 电极构建方式 | library |
| 电极元素 | Au |
| 每侧电极原子数 | 37 |
| 输运温度 (K) | 300.0 |
| TE.dat | demo\2CN-656-4-thiomethylphenyl-fullrelax\2CN-656-4-thiomethylphenyl-fr\aitranss_te\TE.dat |
| TE 数据点 | 301 |
| efermi (Ha) | -0.17951925121 |

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
| T(E_F) [G0] | 0.001217 |
| G/G0 = L0 | 0.001284826 |
| log10(G/G0) | -2.891156 |
| S (µV/K) | -29.0918 |
| k_e (pW/K) | 1.1277 |
| ZT_electronic | 0.022413 |
| 数据点 | 301 |
| 来源 | `demo\2CN-656-4-thiomethylphenyl-fullrelax\2CN-656-4-thiomethylphenyl-fr\aitranss_te\TE.dat` |

## 5. 失败原因

- 无失败阶段。

## 6. 产物清单

| 阶段 | 文件数 | 体积 (MiB) | 关键产物 sha256(前12) |
|---|---|---|---|
| `molecule_opt` | 7 | 1.19 | mol.log:8308ebb9e101 |
| `aims_tip_relax` | 7 | 18.90 | aims.out:300718c28c22, geometry.in:eefa2b185db5 |
| `electrode_build` | 2 | 0.01 | geometry.in:a127578d12c0 |
| `aims_scf` | 9 | 1.75 | aims.out:b74d53c9e964, geometry.in:a127578d12c0 |
| `aitranss_te` | 7 | 0.06 | TE.dat:327136be2b06 |
| `junction_analyze` | 1 | 0.00 | — |
