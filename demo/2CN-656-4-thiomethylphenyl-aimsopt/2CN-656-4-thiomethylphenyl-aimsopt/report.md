# smj run 报告: 2CN-656-4-thiomethylphenyl-aimsopt

> 由 `smj.report` 从 run 目录（state.json / 阶段子目录 / manifest.json）机械生成；
> 数值结论仍须人工复核，报告不含图与独立科学裁决。

## 0. 概览

| 项 | 值 |
|---|---|
| run_id | `2CN-656-4-thiomethylphenyl-aimsopt` |
| 根目录 | `D:\dsh_projects\AutoSMJ\demo\2CN-656-4-thiomethylphenyl-aimsopt\2CN-656-4-thiomethylphenyl-aimsopt` |
| 模式 | slurm-ssh |
| 参数快照 | run_spec.json（sha256 3de72e06c870…）；protocol d0a6584b1c44… |
| 阶段数 / 文件总量 | 6 / 0.043 GiB |
| 阶段状态 | molecule_opt=collected, aims_tip_relax=collected, electrode_build=collected, aims_scf=collected, aitranss_te=collected, junction_analyze=collected |
| 首个未完成阶段 | `—` |
| 生成时间 | 2026-11-13T22:12:47+08:00 |

## 0.1 判定

| 项 | 值 |
|---|---|
| 工程判定 | 通过 |
| 物理判定 | 否证 H0（failed） |
| 允许动作 | 停; 不得进入预测模式(任务卡 §7.3②) |

- 否证 H0 的 error 级判据: `t_ef_magnitude`, `g_in_theory_band`, `seebeck_sign`
- 未过的 warn 级判据: `t_ef_no_peak_near_ef`

## 1. 阶段状态

| 阶段 | 状态 | 目录 | 作业 | 退出码 | 耗时 (s) | 说明 |
|---|---|---|---|---|---|---|
| `molecule_opt` | `collected` | 有 | 26082 | 0 | 964.0 |  |
| `aims_tip_relax` | `collected` | 有 | 26083 | 0 | 1212.8 |  |
| `electrode_build` | `collected` | 有 | — | — | — |  |
| `aims_scf` | `collected` | 有 | 26084 | 0 | 2878.5 |  |
| `aitranss_te` | `collected` | 有 | 26085 | 0 | 249.3 |  |
| `junction_analyze` | `collected` | 有 | — | — | — |  |

## 2. 输入摘要

| 项 | 值 |
|---|---|
| 分子结构 | D:\dsh_projects\AutoSMJ\demo\2CN-656-4-thiomethylphenyl-aimsopt\2CN-656-4-thiomethylphenyl-aimsopt\gaussian_opt\opt.xyz |
| 分子原子数 | 54 |
| 建结阶段 | electrode_build |
| 结原子数 | 128 |
| 每侧电极 Au | 37 |
| 输运温度 (K) | 300.0 |
| TE.dat | D:\dsh_projects\AutoSMJ\demo\2CN-656-4-thiomethylphenyl-aimsopt\2CN-656-4-thiomethylphenyl-aimsopt\aitranss_te\TE.dat |
| TE 数据点 | 301 |
| efermi (Ha) | -0.17942907342 |

## 3. 资源估算

| 阶段 | 墙钟 | --time | 内存 | rank/核 | 核时 | 来源 / 预算 |
|---|---|---|---|---|---|---|
| `molecule_opt` | 00:38:00 | 01:12:00 | 17G | 4 | 2.5 | state; 预算内 |
| `aims_tip_relax` | 00:39:00 | 01:14:00 | 17G | 4 | 2.58 | state; 预算内 |
| `aims_scf` | 01:23:00 | 02:36:00 | 16G | 4 | 5.47 | state; 预算内 |
| `aitranss_te` | 00:49:00 | 01:33:00 | 10G | 1 | 0.81 | state; 预算内 |

## 4. 输运结果

| 量 | 值 |
|---|---|
| T(E_F) [G0] | 0.001245 |
| G/G0 = L0 | 0.001320209 |
| log10(G/G0) | -2.879357 |
| S (µV/K) | -35.1335 |
| k_e (pW/K) | 1.2664 |
| ZT_electronic | 0.029911 |
| 数据点 | 301 |
| 来源 | `D:\dsh_projects\AutoSMJ\demo\2CN-656-4-thiomethylphenyl-aimsopt\2CN-656-4-thiomethylphenyl-aimsopt\aitranss_te\TE.dat` |

## 5. 失败原因

- 无失败阶段。

## 6. 产物清单

| 阶段 | 文件数 | 体积 (MiB) | 关键产物 sha256(前12) |
|---|---|---|---|
| `molecule_opt` | 8 | 22.19 | aims.out:未哈希, geometry.in:0ce581f33704 |
| `aims_tip_relax` | 7 | 19.97 | aims.out:1688caff2079, geometry.in:d00e00b6c383 |
| `electrode_build` | 5 | 0.24 | geometry.in:0638355fa2c6 |
| `aims_scf` | 9 | 1.73 | aims.out:5b895cd86c8c, geometry.in:0638355fa2c6 |
| `aitranss_te` | 7 | 0.06 | TE.dat:b67442e9fe44 |
| `junction_analyze` | 1 | 0.00 | — |
