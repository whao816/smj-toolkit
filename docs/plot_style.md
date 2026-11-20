# T(E) 统一出图版式（style 2026.11）

## 为什么固定一套

单分子结的 `T(E)` 图是每次跑完的"主图"。此前 `smj.analysis.plot_transmission`
自带一套尺寸/字号/配色，跑完后再临时写脚本补图又各画各的（单栏和双栏混用、
有的带标题有的不带、log 轴下限各不相同）——图放不进同一份稿件，也不便横向比较。

现在所有 T(E) 出图都走同一份实现：

- 代码：`src/smj/plotting.py`（版式常量 `STYLE`、尺寸 `FIG_SIZES`、配色 `PALETTE`）
- 命令：`smj plot`
- 兼容入口：`smj.analysis.plot_transmission`（内部委托 `plotting.plot_te`，
  保留原有返回键：`path` / `points` / `log_scale` / `floor`）
- 分析顺手出图：`smj analyze --plot`
- 跑链自动出图：`junction_analyze` 阶段（本地）在写 `results.json` 的同时，
  用同一版式出 `T_E.png` + `T_E.pdf`；没装 matplotlib 时只跳过出图，不影响数值

**版式与物理口径分离**：floor 截断只影响渲染（log 轴不能画 0），
`analysis.json` 里的数值不受影响。

## 固定骨架

| 项目 | 取值 |
|---|---|
| 字体 | Arial → Helvetica → DejaVu Sans；正文 8 pt，轴标题 9 pt，刻度 8 pt，图例 7 pt |
| 边框 | 四边全封闭，0.8 pt |
| 刻度 | 朝内，四边都有；主刻度 3.5 pt、次刻度 1.8 pt；log 轴次刻度 2–9 |
| 线 | 1.1 pt 实线，圆端点 |
| 单曲线颜色 | `#1f4e79`（深蓝） |
| 叠加配色 | `#1f4e79` / `#c0392b` / `#1e8449` / `#7d3c98` / `#b9770e` / `#117a8b`（按顺序固定） |
| x 轴 | `E − E_F`（eV）；`E_F` 处 0.8 pt 灰色虚线 + 顶部 `E_F` 标记 |
| y 轴 | `T(E)`（`G_0 = 2e^2/h`，per spin）；默认 log |
| log 下限 | 最小正值的 1/2，且不低于 `1e-12` |
| 网格 | 无 |
| 标题 | 默认不加（标题交给图注），需要时 `--title` |
| 图例 | 无边框基础风格；实际绘制用白底细框，保证压在曲线上也可读 |
| 尺寸 | single 3.35×2.55 in / double 7.10×3.00 in / wide 6.20×2.80 in |
| 导出 | PNG 600 dpi + PDF 矢量（`pdf.fonttype 42`），`bbox_inches=tight` |

panel 角标 `(a)`、能隙底色、`T(E_F)` / `G/G_0` 注释都是可选项，不改变骨架。

## 文献依据（用户 Zotero 中的同类图）

- Yi et al., *Chin. Phys. B* **25**, 128503 (2016), Fig. 6：四边封闭边框、刻度朝内、
  `E_F` 处竖直虚线、log 透射纵轴、panel 字母在左上角、图例在 panel 内且无边框。
- Guan et al., *Nanotechnology* **34**, 095401 (2023), Figs. 2–3：四边框、刻度朝内、
  panel 内小图例、`E_F` 参考虚线、紧凑的单栏/双栏尺寸。

骨架取了这两处的交集；颜色不照搬（文献多为黑/红/蓝按偏压分色），而是固定
本工具链的调色板，保证不同 run 之间叠加时颜色稳定、可复现。

## 用法

```bash
# 1) 从配置的 run 自动找 TE.dat（默认 single，PNG+PDF）
smj plot -c runs/et7h-bs2/smj.yaml

# 2) 直接给 run 目录 / TE.dat，指定输出
smj plot runs/et7h-bs2/et7h-bs2 -o fig/et7h_TE.png --figsize single

# 3) 线性纵轴 + 指定能量窗口
smj plot path/TE.dat --scale linear --xlim -1 1

# 4) 多条曲线叠加（固定配色 + 白底图例）
smj plot runs/a runs/b --labels "a,b" --figsize double -o cmp.png

# 5) 做 Landauer 分析时顺手出图
smj analyze --plot

# 6) 自旋极化 run：默认自动出双 panel（α/β 通道 + 自旋极化率）
smj plot runs/ttf-cation -o fig/ttf_TE.png            # auto：四列 TE.dat 自动双 panel
smj plot runs/ttf-cation --spin mean                  # 只画双通道均值（旧版式）
smj plot runs/ttf-cation --spin channels              # 单 panel 叠画 α/β
smj plot runs/ttf-cation --spin polarization          # 强制双 panel
```

自旋画法（`plot_te(..., spin=...)` / `--spin`）：

| 值 | 行为 |
| --- | --- |
| `auto`（默认） | TE.dat 有 α/β 两列时用 `polarization`，单列时用 `mean` |
| `mean` | 只画 `(Tα+Tβ)/2`，与旧版式逐像素一致 |
| `channels` | 单 panel 叠画 `Tα`/`Tβ`（深蓝/砖红，`PALETTE[0]`/`PALETTE[1]`） |
| `polarization` | 双 panel：(a) `Tα`,`Tβ`；(b) `P(E)=(Tα−Tβ)/(Tα+Tβ)`；单列数据退化为 `mean` |

双 panel 在 `figsize="single"` 下自动升到 `FIG_SIZES["double"]` 宽度；panel 角标若给了
`panel="a"`，右图自动记为 `(b)`。`meta` 里多返回 `"spin"`（实际画法）与 `"channels"`
（数据列数），便于报告里注明图是按均值还是分通道画的。

Python API：

```python
from smj.plotting import plot_te, plot_te_multi

meta = plot_te("runs/x/aitranss_te/TE.dat", "fig/te.png",
               panel="a", gap=(-1.2, -0.8), formats=("png", "pdf"))
multi = plot_te_multi(["runs/a", "runs/b"], "fig/cmp.png", labels=["a", "b"])
```

返回值含 `style_version` / `floor` / `points` / `paths`（叠加图另有 `curves`），
方便报告和 `RESULTS.md` 里标注图是用哪一版版式出的。

## 维护约定

- 改版式 = 改 `STYLE` / `FIG_SIZES` / `PALETTE` 常量，**同时** bump
  `STYLE_VERSION`，并更新本文件与 `tests/test_plotting.py`；
- 不要在调用方临时改 rcParams 或字号；需要新样式（例如双栏多 panel）优先在
  `plotting.py` 里加函数，而不是复制一份画图代码；
- 文献图只作版式参考，不复制其数据或配色。
