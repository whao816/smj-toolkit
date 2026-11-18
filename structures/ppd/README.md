# PPD (1,4-苯二胺 / BDA) 单分子结结构文件

来源: `runs/<arm>/<arm>/` 下的阶段产物, 坐标原样导出, 未重排未旋转。

| 文件 | 原子数 | 说明 |
|---|---|---|
| `ppd_junction.xyz` / `.in` | 90 | 自由弛豫臂, 输运计算实际几何 (16 分子 + 74 Au) |
| `ppd_pin_junction.xyz` / `.in` | 90 | 刚性 2.55 A 臂, 输运计算实际几何 |
| `ppd_tiprelax18.xyz` | 18 | 自由弛豫臂 aims_tip_relax 收敛几何 (分子 + 2 个 Au tip) |
| `ppd_pin_tiprelax18.xyz` | 18 | 刚性臂 aims_tip_relax 几何 (N 与 tip 被冻结) |

## 原子顺序

两种文件都保持 run 内 `geometry.in` 的顺序: 1-16 为分子
(1-6 环 C, 7-10 环 H, 11/12 为 N, 13-16 为胺 H), 17 起为 Au。
18 原子文件中 17/18 是两侧 tip Au。

## 关键几何

| | ppd (自由) | ppd-pin (刚性) |
|---|---|---|
| N11-Au17 / N12-Au18 | 2.466 / 2.465 A | 2.550 / 2.550 A |
| C-N-Au 角 | 122.5 deg | 177.6 deg |
| 两尖端间距 | 9.167 A | 10.767 A |
| N-N | 5.613 A | 5.667 A |

.xyz 可直接用 VESTA / Avogadro / ASE 打开; .in 是 FHI-aims 原始格式, 可直接丢回集群。
