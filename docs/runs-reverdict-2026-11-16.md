# 旧 run 按当前判据口径重评（2026-11-16）

本文记录把 2026-11-15 判据修订（见 [判据来源核验](criteria-source-research.md)、
[计算流程梳理 §6.2](smj-计算流程梳理.md)、CHANGELOG「判据口径修订」）应用到仓库里全部
**已有 TE.dat 的 run** 之后的结果。重评只重算 manifest / report，不重跑 DFT、不改
`run_spec.json`、不改阶段产物。

## 0. 方法与口径

- 入口：`smj.report.generate_report(Path(run_dir), force=True, include_checks=True, system=<分子族>)`。
- **必须显式传 `system`**：以 run 目录为 target 时 `config=None`，`build_manifest` 无法从配置推断
  分子族，默认 `system=None` 会让全部判据 fail-safe 适用（连 BDT 标定带、Seebeck 符号也套到
  非 BDT 体系上）。各族 run 的 `smj.yaml` 里都写了 `molecule.family`，逐个传入。
- 当前口径 = registry 10 条判据：`gaussian_no_imaginary`、`scf_converged`、`te_present`、
  `grid_independence`、`electrodes_disjoint`、`au_s_bond`（2.40–2.44 Å）、
  `t_ef_resonance_attribution`（0.15–1.0 eV 归因）、`g_in_theory_band`（派生包络，warn）、
  `seebeck_sign`、`k_e_zt_order`。`t_ef_magnitude` 已删除。
- 重评时点：2026-11-16 10:27（+08:00）；当时代码全量测试 256 passed，
  `scripts/gen_vendor_json.py --check` 通过（36 文件 / 10 条本地补丁）。

## 1. `runs/` 正式 run（12 个）

| run | family | 旧判定 | 新判定 | 变化原因 |
| --- | --- | --- | --- | --- |
| `runs/2cn/2cn` | 2CN | inconclusive（`t_ef_no_peak_near_ef`） | **passed** | 近 E_F 峰 +0.241 eV，S<0 → 自洽归入 LUMO 侧，距离落在 0.15–1.0 eV |
| `runs/az/az` | AZ | inconclusive（`t_ef_no_peak_near_ef`） | **passed** | 峰 −0.177/−0.406 eV，S>0 → HOMO 侧，归因通过 |
| `runs/bdt/bdt` | BDT | inconclusive（`t_ef_no_peak_near_ef`） | **passed** | 峰 −0.339/−0.687 eV 在 HOMO 侧；G=0.2789 落在 0.03–0.40 |
| `runs/bdt-dense/bdt-dense` | BDT | inconclusive（`t_ef_no_peak_near_ef`） | **passed** | 同上（密网格） |
| `runs/bdt-anchor24/bdt-anchor24` | BDT | failed（`t_ef_magnitude`、`g_in_theory_band`、`seebeck_sign`） | **failed**（`seebeck_sign`），inconclusive（归因、G） | `t_ef_magnitude` 删除后仍失败：S=−22.2 µV/K 与 BDT 的 HOMO 输运要求 S>0 冲突；峰 +0.048 eV 在 LUMO 侧且 <0.15 eV；G=0.8288 超带 |
| `runs/bdt-anchor24-dense/bdt-anchor24-dense` | BDT | （原无 manifest） | **failed**（`seebeck_sign`），inconclusive（归因、G） | 与上一条同几何、密网格；本次新建 manifest/report |
| `runs/bpy/bpy` | BPY | passed | **passed** | G/Seebeck 对 BPY 不适用（族表）；无近 E_F 峰 |
| `runs/ddt/ddt` | DDT | passed | **passed** | 烷烃族：无近 E_F 共振，G/Seebeck 不适用 |
| `runs/hdt/hdt` | HDT | passed | **passed** | 同上 |
| `runs/odt/odt` | ODT | passed | **passed** | 同上 |
| `runs/sme/sme` | SME | passed | **passed** | 无近 E_F 峰；G/Seebeck 不适用 |
| `runs/sme-axis/sme-axis` | SME | （原无 manifest） | **passed** | 无近 E_F 峰；本次新建 manifest/report |

唯一仍被硬闸门拦下的是 `bdt-anchor24*`：它把 Au–S 锚距人为放到 2.42 Å（`anchor24` 臂），
输运结果 Seebeck 变号（S<0），与 BDT 的 HOMO 输运图像不符。这是该臂本来就该暴露的物理问题，
不是新判据的误伤。

## 2. `demo/` 示例 run（4 个）

| run | family | 旧判定 | 新判定 | 说明 |
| --- | --- | --- | --- | --- |
| `demo/2CN-656-4-thiomethylphenyl/…` | 未登记（`family: null`） | failed（`t_ef_magnitude`、`g_in_theory_band`；inc 归因） | **inconclusive**（归因、G） | 配置显式声明「新体系，只做工程验收」→ fail-safe 全部判据适用；峰 +0.390 eV 与实测 S=+2.8 µV/K 的 HOMO 回退侧别不符 |
| `demo/2CN-656-4-thiomethylphenyl-aimsopt/…` | 未登记 | failed（`t_ef_magnitude`、`g_in_theory_band`、`seebeck_sign`） | **failed**（`seebeck_sign`），inconclusive（G） | fail-safe 下 S=−35.1 µV/K 触发 Seebeck 硬闸门 |
| `demo/2CN-656-4-thiomethylphenyl-fullrelax/…` | 未登记 | failed（同上） | **failed**（`seebeck_sign`），inconclusive（G） | fail-safe 下 S=−29.1 µV/K |
| `demo/azulene-26-dithiol/…` | AZ | inconclusive（`t_ef_no_peak_near_ef`） | **inconclusive**（归因） | 最近峰 −0.0938 eV，距 E_F 过近（<0.15 eV）；与 `runs/az`（−0.177 eV，passed）不是同一几何/取点 |

`demo/` 的 2CN 系列没有登记分子族，这是配置里的显式选择，不是遗漏。本次不替它编一个族名：
按当前 fail-safe 语义，未登记族的 run 会被 BDT 标定带和 Seebeck 符号检查覆盖，结论就是
「inconclusive / failed」——这正是「未登记体系不得声称物理通过」的设计意图。若要把 demo 2CN
作为已登记体系对待，应先在 `smj.yaml` 写 `family: 2CN` 再重评。

## 3. 产物与边界

- 每个 run 的 `manifest.json` 与 `report.md` 均已用当前版本重写（`force=True`）；
  `report.md` 的「生成时间」即重评时间，`manifest.json` 的 checks 段是当前 10 条判据的输出。
- `run_spec.json`、阶段目录、`review.json`、`state.json` 均未改动。
  **注意**：`review.json` 里的几何 anchor 明细（`junction_geom_anchor_aus`）是旧口径生成的历史
  阶段产物，detail 文案仍写 2.40–2.49 Å；判据评估不读它，但人工看 review 时请以
  [计算流程梳理 §6.1](smj-计算流程梳理.md) 的 2.40–2.44 Å 为准。
- `demo/azulene-26-dithiol` 的 run 目录没有 `run_spec.json`，重评时 manifest 带一条
  「run_spec.json 缺失」警告（不伪造参数快照）；其余 run 的 protocol sha256、输入 sha256 原样保留。
- 跳过的目录：`demo/2CN-656-4-thiomethylphenyl-g16tip`（无 TE.dat）、`runs/compare-2026-11-15/*`
  （只做 molecule_opt，无输运产物）。`runs/compare-2026-11-15/REPORT.md` 是分子优化对照，
  不涉及判据判定，未改。

## 4. 复现

```powershell
$env:PYTHONPATH='src'
@'
from pathlib import Path
from smj.report import generate_report
generate_report(Path(r"runs/bdt/bdt"), force=True, include_checks=True, system="BDT")
'@ | .\.venv\Scripts\python.exe -
```

逐 run 的 `(run, family)` 清单与新旧判定对照见本次重评脚本输出
（临时文件 `.tmp/reverdict-2026-11-16.json`，未纳入版本库）。
