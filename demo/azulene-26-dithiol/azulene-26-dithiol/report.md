# smj run 报告: azulene-26-dithiol

> 由 `smj.report` 从 run 目录（state.json / 阶段子目录 / manifest.json）机械生成；
> 数值结论仍须人工复核，报告不含图与独立科学裁决。

## 0. 概览

| 项 | 值 |
|---|---|
| run_id | `azulene-26-dithiol` |
| 根目录 | `D:\dsh_projects\AutoSMJ\demo\azulene-26-dithiol\azulene-26-dithiol` |
| 模式 | — |
| 阶段数 / 文件总量 | 6 / 0.014 GiB |
| 阶段状态 | gaussian_opt=collected, aims_tip_relax=collected, au60_build=collected, au_relax=collected, aims_scf=collected, aitranss_te=collected, junction_analyze=collected |
| 首个未完成阶段 | `—` |
| 生成时间 | 2026-11-12T14:49:05+08:00 |

## 0.1 判定

| 项 | 值 |
|---|---|
| 工程判定 | 通过 |
| 物理判定 | 无法判定（inconclusive） |
| 允许动作 | 停; 须补数据或降级声明 |

- 未过的 warn 级判据: `t_ef_no_peak_near_ef`

## 1. 阶段状态

| 阶段 | 状态 | 目录 | 作业 | 退出码 | 耗时 (s) | 说明 |
|---|---|---|---|---|---|---|
| `gaussian_opt` | `collected` | 有 | — | — | — | 复用主仓 B3LYP/def2-TZVP opt+freq（0 虚频）日志 |
| `aims_tip_relax` | `collected` | 有 | 26020 | 0 | 187.5 |  |
| `au60_build` | `collected` | 有 | — | — | — |  |
| `au_relax` | `collected` | — | 26029 | 0 | 4210.4 |  |
| `aims_scf` | `collected` | 有 | 26033 | 0 | 2230.9 |  |
| `aitranss_te` | `collected` | 有 | 26035 | 0 | 186.3 |  |
| `junction_analyze` | `collected` | 有 | — | — | — |  |

## 2. 输入摘要

| 项 | 值 |
|---|---|
| 分子结构 | D:\dsh_projects\AutoSMJ\demo\azulene-26-dithiol\azulene-26-dithiol\gaussian_opt\opt.xyz |
| 分子原子数 | 20 |
| 建结阶段 | au_relax |
| 结原子数 | 92 |
| 每侧电极 Au | 37 |
| 输运温度 (K) | — |
| TE.dat | D:\dsh_projects\AutoSMJ\demo\azulene-26-dithiol\azulene-26-dithiol\aitranss_te\TE.dat |
| TE 数据点 | 190 |
| efermi (Ha) | -0.188452434085 |

## 3. 资源估算

| 阶段 | 墙钟 | --time | 内存 | rank/核 | 核时 | 来源 / 预算 |
|---|---|---|---|---|---|---|
| `gaussian_opt` | 00:15:00 | 00:29:00 | 6G | 1 | 0.25 | state; 预算内 |
| `aims_tip_relax` | 00:18:00 | 00:34:00 | 13G | 4 | 1.17 | state; 预算内 |
| `au_relax` | 01:23:00 | 02:38:00 | 24G | 4 | 5.53 | state; 预算内 |
| `aims_scf` | 00:31:00 | 00:58:00 | 8G | 4 | 2.03 | state; 预算内 |
| `aitranss_te` | 00:30:00 | 00:57:00 | 6G | 1 | 0.5 | state; 预算内 |

## 4. 输运结果

| 量 | 值 |
|---|---|
| T(E_F) [G0] | 0.761846 |
| G/G0 = L0 | 0.746147187 |
| log10(G/G0) | -0.127175 |
| S (µV/K) | 25.9362 |
| k_e (pW/K) | 388.0045 |
| ZT_electronic | 0.030069 |
| 数据点 | 190 |
| 来源 | `D:\dsh_projects\AutoSMJ\demo\azulene-26-dithiol\azulene-26-dithiol\aitranss_te\TE.dat` |

## 5. 失败原因

- 无失败阶段。

## 6. 产物清单

| 阶段 | 文件数 | 体积 (MiB) | 关键产物 sha256(前12) |
|---|---|---|---|
| `gaussian_opt` | 6 | 0.40 | mol.log:2a5b13347cb5 |
| `aims_tip_relax` | 8 | 12.02 | aims.out:f625af15ad5a, geometry.in:55615e5c2f2a |
| `au60_build` | 5 | 0.19 | geometry.in:9b8d6dd4f538 |
| `aims_scf` | 9 | 1.60 | aims.out:275126b950dd, geometry.in:0a5af56c3eb7 |
| `aitranss_te` | 7 | 0.04 | TE.dat:01d89887fe02 |
| `junction_analyze` | 1 | 0.00 | — |
