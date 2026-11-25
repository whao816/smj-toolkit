# Au lead 映射与主层厚度审计（2026-11-25）

审计对象：`scripts/gen_au_lead_fc.py` 生成的体相 Au(111) lead 力常数、
`src/smj/phonon/fc.py` 的组装逻辑，以及 `src/smj/stages/thermal_analyze.py`
的切块口径。审计用 nap-sh 现有产物（`runs/nap-sh-thermal`、`runs/lead-fc`），
**没有提交新的 DFT 作业**。

## 结论先行

**不需要重算那 54 个 SCF。** 现有 6 层体相力常数（54 个位移点）已经包含
重建两层主层 + 两层 hopping 所需的全部信息；主层厚度问题可以纯后处理修复。

但当前生产口径（一层主层）**不能用**：修正 lead 层序后它给
`G_ph = 15.93 pW/K`，而两层主层给 `11.12 pW/K`，差 **−30.2%**。
两层主层与三层主层只差 **2.3%**，低于 D11 的 5% 门限。

## 缺陷 1：lead 层序反了（已修）

`lead_fc.npy` 按「参考层 0（次外层）→ 参考层 1（主层）」排列，前 54 维
是层 0。`thermal_analyze` 原来传 `left_pl + left_outer`（主层在前），
等于把 `K(0,1)` 当 `K(1,0)`。矩阵仍对称、ASR/形状检查全过，只错一次
转置：

| 口径 | `G_ph(300 K)` |
|---|---:|
| 修复前 | 13.7166 pW/K |
| 修复后（L=1） | **15.9264 pW/K** |

回归：`tests/phonon/test_lead_mapping.py::test_lead_block_is_ordered_outer_then_principal`
用非对称层间耦合先复现失败，再钉住 `K(principal,outer) = K(0,1)^T`。

## 缺陷 2：截断轴用了裸坐标距离（已修）

lead 块来自周期体相，组装矩阵里的 lead 对可能跨周期边界。`_atom_distances`
原来只算裸距离，导致：

- 原子对 (18,35)：裸距 **11.889 Å**，面内最小像 **3.115 Å**；
- 原子对 (54,71)：裸距 **11.907 Å**，面内最小像 **3.141 Å**。

所以之前 `cutoff_sensitivity` 报的 **23% 散布**、`nominal_shortfall=4.3%`、
「10 Å 未收敛」都是伪结论。改用面内最小像后：

| `r_c` (Å) | 保留对数 | `G_ph` (pW/K) |
|---:|---:|---:|
| 4 | 324 | 14.361 |
| 5 | 538 | 15.289 |
| 6 | 576 | 15.484 |
| 7 | 576 | 15.484 |
| 8 | 612 | 15.926 |
| 10 | 612 | 15.926 |
| 12 | 612 | 15.926 |

8/10/12 Å 三档逐位相同。剩余 9.8% 散布只来自把 `r_c` 刻意切到 6 Å 以下
（切掉真实近邻），**不是名义截断不足**。

回归：`test_atom_distances_wrap_in_plane_periodic_images`。

## 主层厚度：单层不够（需改协议，未重算 DFT）

从 `runs/lead-fc` 的 54 个原始力文件重建 6 层周期胞，物理层序（由器件向
外）`q = [1, 0, 5, 4, 3, 2]`。主层厚度 L 的含义：`d_PL` 取 L 层，
`V_PL` 取到下一组 L 层；器件侧耦合只覆盖实际存在的 layer 1 / layer 0，
虚拟层列补零。

| 主层厚度 | `G_ph(300 K)` | 带顶以上最大透射 |
|---|---:|---:|
| L=1（当前生产） | 15.9264 pW/K | 8.4e-8 |
| L=2 | **11.1211 pW/K** | 2.6e-7 |
| L=3 | 11.3759 pW/K | 4.1e-7 |

- L=1 → L=2：**−30.2%**，说明一层主层漏掉了次外层参与的自能以及更远层耦合；
- L=2 → L=3：**+2.3%**，低于 D11 的 5% 门限；
- 1D 最近邻链解析验证：L=1 与 L=2 的 decimation 都给出 `T=1`，
  说明多层调用方式本身正确，差异来自 Au 力常数的长程耦合。

L=2 用到的层间对最远 7.06 Å（1–4）、其余 4.71 Å，都在现有 10 Å 截断内，
**因此不需要为截断重算**。

## 失败过的口径（不要再用）

把体相 K 直接按 6/8/10/12 Å 整体置零后做 L=2，四档全部被带顶闸门拒绝
（报 `T_ph` 在 19.50 meV 以上仍非零，最大值落在 19.69–19.81 meV）。
原因是截断破坏 ASR/平移不变性，凭空造出带顶以上通道。**判断体相截断收敛
必须用保持平移不变性的构造，不能对整块 K 直接按距离置零。**

## 建议的落地顺序

1. 保留已修的层序与最小像逻辑（本次提交）。
2. 用现有 `runs/lead-fc` 重建 lead 资产：principal 2 层 + hopping 2 层，
   资产 schema 升到 v2（记录 `principal_layers` 与 hopping 文件）。
3. `thermal_analyze` 的 `d_L`/`V_L` 改用 2 层主层块与 hopping 块，
   器件侧 `v_l` 覆盖两层真实 lead。
4. 重跑 `thermal_analyze --collect --force`，预期 `G_ph ≈ 11.12 pW/K`；
   把 L=3 的 11.38 pW/K 写进收敛轴作为对照。
5. 此时才恢复 az/ind 同协议 run；ODT/D9 仍作为单独控制臂。

## 复现

```powershell
$env:PYTHONPATH="src"
python -m pytest tests/phonon/test_lead_mapping.py tests/phonon/test_environment.py -q
python -m pytest -q
```

诊断脚本（已在审计后删除）：`.tmp/lead_audit.py`、`.tmp/principal_layer.py`。
