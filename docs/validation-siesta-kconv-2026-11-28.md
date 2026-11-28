# SIESTA 横向 k 点与屏蔽验证记录（2026-11-28）

> 本记录回答的问题：当前 SIESTA/TranSIESTA 零偏压路线给出的 `G/G0`，在横向 k 点加密后是否还在 1% 以内；电极边界是否已经屏蔽到体相。两项证据齐全后，`report`/`review` 才允许把物理判定从 `inconclusive` 升级为 `passed`。
>
> **状态：已完成（2026-11-28）。** bpy 的 3×3→5×5 与 5×5→7×7 均未达到 1% 判据；ppd 的 3×3→5×5 达到。屏蔽两套体系均 `passed`。这是一次“把 inconclusive 变成可判定结论”的验证：bpy 的旧 3×3 结果被证明未收敛，ppd 的 5×5 结果在 1% 内稳定。

## 1. 判据

在 `docs/siesta-zero-bias-transport-design.md` §3 中预注册：

- **横向 k 收敛**：以 3×3 为基线，器件 SCF 加密到 5×5（必要时 7×7），比较各自 `E_F` 处的 `G/G0`；相邻网格相对变化 `< 1%` 记为通过。TBtrans 的横向 k 单独加密（5×5→7×7）再比较一次，避免把 SCF 密度收敛与输运积分收敛混为一谈。
- **屏蔽**：取电极中间体相层的 Mulliken 平均电荷为参考；器件两侧最外两层（面向电极的 Au 层）的 |Δq| 均 `≤ 0.05 e`，且电极体相层间 spread `≤ 0.05 e`，记为通过。
- 任一项缺失或不满足，物理判定保持 `inconclusive`（证据缺失）或 `failed`（明确越界）。

## 2. 计算设置

| 项 | 3×3 基线 | 5×5 加密 |
|---|---|---|
| 电极 | Au(111) 3×3，6 principal / 8 slab / 1 movable，`3x3x12` | 同结构，`5x5x12` |
| 器件 k 点 | `3x3x1` | `5x5x1` |
| TBtrans k 点 | `5x5x1` | `5x5x1`（另跑 7×7 对照） |
| 结几何 | 独立弛豫 | `siesta.reuse_relax_from` 复用，`reused_from` 记录在 `relax.json` |
| Mulliken | 未开 | `WriteMullikenPop 1`（电极与器件都开） |

其余设置同 `docs/validation-siesta-2026-11-27.md`：PBE、DZP、`MeshCutoff 300 Ry`、`ElectronicTemperature 0.025 eV`、300 K、能量窗 `[-1, +1] eV`、步长 5 meV。5×5 几何没有在新 k 网格上重新弛豫，这是有意的成本控制；报告会注明这一点。

## 3. bpy 5×5 SCF 与屏蔽

`runs/siesta-flat-pair/bpy-k5/transiesta_scf/scf.json`：

| 量 | 值 |
|---|---:|
| SCF 收敛 | true |
| 开放边界 SCF 迭代 | 18 |
| `E_F` | 0.2679214969 eV |
| 电荷偏差 | −0.005770 e |
| 目标电荷 | 2794.0 e |
| 轨道数 | 2932 |
| 屏蔽状态 | **passed** |

电极参考（5×5 电极，job 26970）：`fermi_ev = −1.0703260417 eV`，`mulliken_bulk_charge_e = −1.85×10⁻⁸ e`，层间 spread `1.3×10⁻⁵ e`。器件两侧的层分辨 Δq：

| 侧 | 最外层 | 第二层 | 靠近分子的最内层 |
|---|---:|---:|---:|
| 左 | −7.43×10⁻⁶ e | +6.07×10⁻⁴ e | −2.061×10⁻² e |
| 右 | −1.03×10⁻⁵ e | +6.10×10⁻⁴ e | −2.061×10⁻² e |

最外两层都在 10⁻³ e 以下，离 0.05 e 的判据很远；到最内层才升到约 −0.021 e，符合表面电荷向体相衰减的预期。左右两侧的层 profile 在 10⁻⁴ e 内对称，最内层偏差也同为 −0.0206 e，说明两侧电极与接触没有单侧伪影。因此 **屏蔽这一半已经从 `unverified` 变成有数据的 `passed`**。

注意：`bpy-k5` 的 `scf.json` 是在进程启动时加载的旧版代码采集的，侧级字段名仍是 `surface_delta_e`；新版代码已改名为 `innermost_delta_e`，数值含义不变。

### 3.1 3×3 与 5×5 的能级移动

| 量 | 3×3 | 5×5 | 差 |
|---|---:|---:|---:|
| 电极 `E_F` | −1.1187587182 eV | −1.0703260417 eV | +0.04843 eV |
| 器件 `E_F` | 0.1981398615 eV | 0.2679214969 eV | +0.06978 eV |
| 器件 `E_F` − 电极 `E_F` | 1.31689858 eV | 1.33824754 eV | +0.02135 eV |

也就是说，k 加密不仅改变 G，也把费米能级抬高了约 70 meV。最终比较必须同时报告 `G/G0` 相对变化和 `E_F` 偏移；只报 G 的百分比会掩盖这一系统移动。

## 4. bpy k 收敛

| 比较 | 基线 `G/G0` | 加密 `G/G0` | 相对变化 | 判据 |
|---|---:|---:|---:|---|
| 3×3 SCF + 5×5 TBtrans → 5×5 SCF + 5×5 TBtrans | 1.814203377×10⁻³ | 2.228801541×10⁻³ | **+22.853%** | <1%，**failed** |
| 5×5 SCF + 5×5 TBtrans → 5×5 SCF + 7×7 TBtrans | 2.228801541×10⁻³ | 2.354640452×10⁻³ | **+5.646%** | <1%，**failed** |

对应的 `G` 分别为 `1.405661418×10⁻⁷ S`、`1.726895879×10⁻⁷ S`、`1.824397022×10⁻⁷ S`。结论很明确：**bpy 的 3×3 结果没有收敛，且 5×5 仍不够。** 3×3→5×5 变化约 23%，5×5→7×7 TBtrans 又变化 5.6%；在 1% 判据下，bpy 需要 7×7 SCF 与至少 9×9 TBtrans 才可能过关，当前只能给出“旧值是低估”的定量结论，不能给出收敛值。

`E_F` 也随 k 点移动：3×3 器件 `E_F=0.1981398615 eV`，5×5 为 `0.2679214969 eV`（+69.78 meV）。7×7 TBtrans 复用 5×5 SCF，因此 `E_F` 不变，只改变透射积分。

## 5. ppd k 收敛

| 比较 | 基线 `G/G0` | 加密 `G/G0` | 相对变化 | 判据 |
|---|---:|---:|---:|---|
| 3×3 SCF + 5×5 TBtrans → 5×5 SCF + 5×5 TBtrans | 4.439993888×10⁻¹ | 4.407773409×10⁻¹ | **−0.726%** | <1%，**passed** |

`G` 从 `3.440147992×10⁻⁵ S` 降到 `3.415183270×10⁻⁵ S`；5×5 `E_F=0.2712939129 eV`，比 3×3 的 `0.1793218747 eV` 高 `91.97 meV`。`E_F` 移动比 bpy 还大，但 `G/G0` 只变了 0.7%，说明 ppd 的透射在费米能级附近更平缓，积分对 k 点不敏感。

ppd 的 5×5 屏蔽同样 `passed`：左最外两层 Δq=`−2.13×10⁻⁵`/`+5.80×10⁻⁴ e`，右=`−6.47×10⁻⁶`/`+6.15×10⁻⁴ e`，电极层间 spread=`9.0×10⁻⁶ e`，SCF 电荷偏差 `−0.00524 e`。因此 `ppd-k5` 的 `report`/`review` 判定为 **passed**（仅数值判据）；`bpy-k5` 与 `bpy-k5-tb7` 均为 **failed**。

**结论。** 这次验证没有把 bpy 变成“通过”，而是把它从“未知”变成了“可定量地不收敛”：3×3 的 `G/G0` 比 5×5 低约 23%，5×5 的 TBtrans 结果在 7×7 下又低约 5.6%，所以 bpy 的 3×3 旧值不能用于定量比较。ppd 的 `G/G0` 在 3×3→5×5 只变 0.7%，在 1% 判据内通过；但它的 `E_F` 移动了 92 meV，说明“G 稳定”不等于“能级稳定”，引用时必须同时报告两者。两套体系的屏蔽都已从写死的 `unverified` 变成逐层 Mulliken 证据，且最外两层在 10⁻³ e 以下，电极边界可按体相处理。

## 6. 复现

```bash
smj stage siesta_electrode      -c runs/siesta-flat-pair/bpy-k5.yaml --run --force
smj stage siesta_junction_relax -c runs/siesta-flat-pair/bpy-k5.yaml --run --force
smj stage transiesta_scf        -c runs/siesta-flat-pair/bpy-k5.yaml --run --force
smj stage tbtrans_te            -c runs/siesta-flat-pair/bpy-k5.yaml --run --force
smj stage junction_analyze      -c runs/siesta-flat-pair/bpy-k5.yaml --run --force
smj siesta-kconv runs/siesta-flat-pair/bpy-siesta-flat runs/siesta-flat-pair/bpy-k5 \
    --json-out runs/siesta-flat-pair/bpy-k5/kconv.json
```

7×7 TBtrans 现在用 `runs/siesta-flat-pair/bpy-k5-tb7-v2.yaml` 复现：该配置用
`siesta.electrode_reuse: bpy-k5/siesta_electrode`、`siesta.reuse_relax_from:
bpy-k5/siesta_junction_relax`、`siesta.reuse_scf_from: bpy-k5/transiesta_scf`
复用 5×5 链，只把 `siesta.tbtrans_kgrid` 设为 `[7,7,1]`：

```bash
smj stage siesta_electrode      -c runs/siesta-flat-pair/bpy-k5-tb7-v2.yaml --run --force
smj stage siesta_junction_relax -c runs/siesta-flat-pair/bpy-k5-tb7-v2.yaml --run --force
smj stage transiesta_scf        -c runs/siesta-flat-pair/bpy-k5-tb7-v2.yaml --run --force
smj stage tbtrans_te            -c runs/siesta-flat-pair/bpy-k5-tb7-v2.yaml --run --force
smj stage junction_analyze      -c runs/siesta-flat-pair/bpy-k5-tb7-v2.yaml --run --force
```

前三个阶段只做本地复用与摘要核对，不重跑 SIESTA；`tbtrans.fdf` 的 `TBT.k`
为 `7 7 1`，`scf.json` 的 `reused_from` 指向 `bpy-k5/transiesta_scf`。验证时
已用该配置只执行到 `tbtrans_te --prepare` 检查输入，未重跑 7×7 输运。

## 7. 范围

即使两项都通过，结论也只覆盖“数值判据”：它说明当前设置下的 `G/G0` 对横向 k 与电极屏蔽已经稳定，不说明 PBE/DZP、未含 vdW 的接触模型、或 300 K 弹性输运近似足以给出实验可比的电导。`report`/`review` 的 `passed` 必须与这段范围声明一起读。
