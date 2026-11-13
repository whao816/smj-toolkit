# AIMS vs Gaussian 全原子 tip 弛豫对比（2CN-656-4-thiomethylphenyl）

## 目的

同一 56 原子 tip 初猜（54 原子分子 + 2 个 Au tip，S 为 1-based 45/46，
初始 Au–S = 2.35 Å）分别用 FHI-aims 与 Gaussian 做全原子弛豫，比较
Au–S 键长与 S–C–Au 连接角。

## AIMS（已完成）

- 阶段：`demo/2CN-656-4-thiomethylphenyl-fullrelax/.../aims_tip_relax`
- 作业 26067，1088 s，56 原子，`relax_geometry bfgs 1.e-2`，无约束。
- 结果：Au–S = **2.4267 / 2.4246 Å**；C–S–Au = **111.35 / 111.32°**；
  S–C–Au = **40.49 / 40.49°**；`contact_pinned=false`。
- 初猜 → 弛豫结构 RMSD 0.5265 Å，Au 位移 2.553 Å，S 位移 0.057 Å。

## Gaussian（未完全收敛）

模板最初为 B3LYP/def2-SVP；56 原子含 Au ECP，SCF 在 129 周期后
`Convergence failure`。改为 PBEPBE/def2-SVP、`scf=(qc,maxcycle=256)`、
32 核后：

1. 从初猜结构出发：42 min 后 Gaussian 优化器在第 1 步线搜索报
   `Convergence failure`（20h51 CPU）。
2. 从 AIMS 弛豫结构出发（pbe2）：56 min 到达第 2 步；随后从该
   第 2 步几何重启（pbe3），30 min 到达第 2 步，力已收敛
   （max force 0.0021，RMS 0.0002），但位移未收敛
   （max displacement 0.195，RMS 0.031）。因单步仍需 ~20 min，
   为避免无限等待而停止。

第 2 步几何（未完全收敛）与 AIMS 比较：

- Au–S = **2.5082 / 2.5082 Å**（比 AIMS 长约 0.08 Å）
- C–S–Au = **111.39 / 111.44°**；S–C–Au = **41.02 / 40.99°**
- Kabsch RMSD = **0.160 Å**，最大位移 0.553 Å（Au 0.40 Å）

## 结论与边界

- 连接角在两种方法间几乎一致（差 ≤0.5°），说明局部成键取向稳健。
- Gaussian 的 Au–S 比 AIMS 长约 0.08 Å，但 Gaussian 结构尚未完成
  位移收敛，不能据此断言方法间存在系统偏差。
- 全 56 原子簇的 Gaussian DFT 优化在当前 32 核/4 h 配额下成本过高；
  这不是流水线错误。若要得到收敛的 Gaussian 结构，建议：缩小到
  局部接触模型（S + Au + 近邻 C）、使用 ONIOM/QM-MM、或分配更长
  walltime。

## 文件

- `mol_pbe_from_aims.gjf` / `.log`：从 AIMS 结构出发的 Gaussian 输入与日志
- `mol_pbe_step2_restart.log`：第 2 步重启日志
- `last_standard_step2.xyz`：第 2 步标准取向几何（未完全收敛）
