# smj-toolkit 0.2.0 验证记录

> 目的：记录 0.2.0 发布前在真实集群上跑通的三套体系（BDT / AZ / 2CN）及它们的输运数值，
> 作为「整链可用」的证据；同时保留已知事件与口径限制，避免把验证样本当生产结论。

## 1. 环境与方法

| 项 | 值 |
|---|---|
| 版本 | smj-toolkit 0.2.0（commit/tag 见仓库 `v0.2.0`） |
| 集群 | `gpu-outside`（Slurm 20.11.9，partition `main`，8 核 + 1 GPU / 24G 每作业） |
| FHI-aims | `/data/apps/fhi-aims.221103-cuda12.8-oneapi2025.3/bin/aims.221103.scalapack.mpi.x` |
| aitranss | `aitranss.071813.openmp.ifx2025.3.x` |
| 交换关联 / 基组 | PBE / light（species 目录 `runs/<run>/species`：H、C、N、S、Au） |
| 电极 | `electrode.source: library`，`fcc111.au.clusters.37_atoms/fcc111.cluster.natoms37.x1.aims`（37 Au/侧，nlayers=3） |
| 分子优化 | `molecule.opt_engine: aims`，`relax_geometry bfgs`，力阈值 `1.e-2` eV/Å |
| 结几何 | tip 全原子弛豫 → 库团簇扩展 → 冻结；`pin_contact: false` |
| SCF | 静态，`KS_method serial`，`output aitranss`（保证 `omat.aims` 布局可被 aitranss 读取） |
| 输运 | 能量窗 `[-0.22, 0.0002, -0.16]` Ha（301 点），T = 300 K，自能 σ_near 0.05 / σ_far 0.10 Ha |
| A2 | 人工签署 `smj review --approve --who whao`，表面审计（15/37 最外层、面积 3.6009 Å²、法向-轴 0.000°）PASS |

## 2. 三套体系

| run | 分子 | 锚定化学 | 分子原子数 | 结原子数 | 表面三元组（左/右） | 松弛后 Au–S (Å) |
|---|---|---|---|---|---|---|
| `bdt` | 1,4-benzenedithiol（C6H6S2） | 硫醇盐（去质子，`strip_thiol_h: true`） | 14 | 86 | [40,41,44] / [76,77,80] | ≈2.253 |
| `az` | azulene-2,6-dithiol（C10H6S2） | 硫醇盐（去质子） | 20 | 92 | [46,47,50] / [82,83,86] | 2.2511 / 2.2549 |
| `2cn` | 2CN-656-4-thiomethylphenyl | 硫醚（`strip_thiol_h: false`） | 54 | 128 | [82,83,86] / [118,119,122] | 2.4248（配位接触） |

## 3. 输运结果（301 点，T = 300 K）

| run | efermi (Ha) | L0 = G/G0 | log10(G/G0)_like | S (µV/K) | k_e (pW/K) | ZT_e | T(E_F)（最近点） | max T（位置 E−E_F） |
|---|---|---|---|---|---|---|---|---|
| `bdt` | −0.188142 | 0.278853973 | −0.554623 | +34.488 | 165.5769 | 0.046561 | 0.272587 @ −1.59 meV | 1.0545 @ −0.687 eV |
| `az` | −0.188086748 | 0.380578051 | −0.419556 | +54.074 | 223.9924 | 0.115479 | 0.349526 @ +2.36 meV | 1.1795 @ −0.406 eV |
| `2cn` | −0.179464717 | 0.001324925 | −2.877808 | −33.835 | 1.2426 | 0.028373 | 0.001253 @ +1.76 meV | 0.8710 @ +0.241 eV |

判据（`smj report` → `manifest.verdict`）：

| run | 工程判定 | 物理判定（0.2.0 原记录） | 按 2026-11-15 口径重评 |
|---|---|---|---|
| `bdt` | 通过 | inconclusive（`t_ef_no_peak_near_ef`） | passed（HOMO 峰 −0.339 eV 已归因） |
| `az` | 通过 | inconclusive（`t_ef_no_peak_near_ef`） | passed（HOMO 侧峰 −0.177 eV 已归因） |
| `2cn` | 通过 | inconclusive（`t_ef_no_peak_near_ef`） | passed（LUMO 侧峰 +0.241 eV，按 S<0 自洽归因） |

三套都是「管线通过」：工程判定检查的是数据完整性与 fail-closed 行为
（SCF 收敛、TE.dat 在场、阶段哈希闭环）。0.2.0 记录里物理判定为 `inconclusive`，
原因是当时的 `t_ef_no_peak_near_ef` 把「E_F 1 eV 内存在透射峰」一律当错误信号；
2026-11-15 起该判据改为 `t_ef_resonance_attribution`（按族预期通道归因，见
[`docs/smj-计算流程梳理.md`](smj-计算流程梳理.md) §6.2 与
[`docs/criteria-source-research.md`](criteria-source-research.md)），三套重评均 `passed`。
BDT 标定族（`g_in_theory_band`、`seebeck_sign`）只对 `family: BDT` 生效，AZ/2CN 不做该族比对。

## 4. 阶段耗时（实测，单作业 8 核 + 1 GPU）

| run | molecule_opt | aims_tip_relax | aims_scf | aitranss_te | 合计 |
|---|---|---|---|---|---|
| `bdt` | 32.0 s（job 26101） | 94.2 s（job 26102） | 1563.9 s（job 26103） | 250.7 s（job 26107） | ≈32 min |
| `az` | 32.0 s（job 26118） | 311.5 s（首次中止于 26119；重启 26120 后 20 步收敛） | 1963.4 s（job 26122） | 218.6 s（job 26124） | ≈43 min |
| `2cn` | 963.3 s（job 26117） | 1274.6 s（job 26121） | 2896.4 s（job 26123） | 218.8 s（job 26125） | ≈89 min |

- 128 原子 SCF 实测 48.3 min，高于资源模型对同规模的外推（约 31 min）；`--time` 为
  1.5–1.9 倍余量，本次未触顶。大规模体系建议按「实测 + 50%」再给 `--time`。
- TE 阶段与体系大小关系不大（218–251 s），主要成本是 SCF。

## 5. 已记录事件

1. **AZ tip 弛豫 trust radius 中止**：job 26119 在第 80 步报
   `The trust radius has fallen below its allowed minimum value, 0.100000E-03 Angstrom`
   （伴随 `forrtl: severe (66): output statement overflows record`）。按 FHI-aims 手册
   从 `geometry.in.next_step` 重启、`init_hess diag 97`，job 26120 用 20 步收敛。
   这是 BFGS 在能量/力噪声下的常见失败，不是 SMJ 的编排错误；`state.json` 里
   `aims_tip_relax` 仍保留首次失败 job 的退出码 66，几何已由重启结果替换。
2. **陈旧 review.json**：2CN 的 `review.json` 在一次 A2 签署前生成（status `blocked`,
   `a2_missing`），签署后未重写；`smj review --force` 重算后为 `passed`。说明：
   `review` 默认不覆盖已存在的 `review.json`，改过 A2 记录后需 `--force` 刷新。
3. **2CN 的 T(E_F) 很低**（G/G0 ≈ 1.3×10⁻³）：硫醚（S–CH3）锚定是配位接触，
   与硫醇盐的直接 Au–S 化学键不同，弱耦合是预期行为；该体系不在 BDT 标定族内。
4. **网格收敛对照（BDT）**：同一 BDT 体系把 estep 从 0.0002 加密到 0.00005 Ha
   （301→1201 点），L0/S/k_e/ZT_e 逐位一致（相对差 < 1e-8）；`smj report` 的 T(E_F)
   数字差异来自最近格点位置（−1.59 meV vs −0.23 meV），插值到 E=0 一致。
   同一对照在「冻结 S–Au=2.42 Å」几何上也成立（runs/bdt-anchor24、runs/bdt-anchor24-dense）。
   接触效应的对照结论见 `runs/_compare/BDT-对照实验报告.md`。

## 6. 复现

```bash
# 每个 run 目录下（runs/<run>/smj.yaml），任选其一：
smj run     -c runs/bdt/smj.yaml --from molecule_opt --to electrode_build
smj review  -c runs/bdt/smj.yaml --approve --who <name>   # A2 人工确认（需看图）
smj run     -c runs/bdt/smj.yaml --from aims_scf --to junction_analyze
smj report  -c runs/bdt/smj.yaml
smj analyze -c runs/bdt/smj.yaml        # 独立解析器交叉校验 + T(E) 图（需 .[plot]）
```

`aims_tip_relax` 首次中止时，按 §5.1 的方式重启后再 `smj stage aims_tip_relax --collect`。
run 目录布局见 README §6；`state.json` 是唯一状态源。

## 7. 版本验证矩阵（0.2.0）

| 项 | 结果 |
|---|---|
| `pytest`（本机 Windows / Python 3.12.10） | 125 passed |
| `pytest`（WSL Ubuntu 26.04 / Python 3.14.4） | 125 passed |
| GitHub Actions `ci.yml`（run 43610992842，rc1 代码） | 9/9 success：Ubuntu 3.10–3.14、Windows 3.12–3.14、wheel/sdist 构建+安装冒烟 |
| wheel 安装（干净 venv，WSL 3.14） | `smj verify` 6/6、`init` + `plan` 正常、`pip check` 无破损 |
| sdist 安装（干净 venv，WSL 3.14） | 现场构建 wheel 后 `smj verify` 通过 |

> 说明：0.2.0 相对 0.2.0rc1 只增加本验证记录与发布流程文档（无代码路径变化）；
> CI 在 `v0.2.0` 推送后会再跑一遍同矩阵。
