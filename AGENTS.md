# 项目级 Agent 规则（AutoSMJ / smj-toolkit）

## 运行范围：只跑用户明说的

- 只对用户明确点名的分子、异构体、体系与 protocol 跑链。
- **不要默认加跑 pinned 对照、敏感性扫描、替代锚定、团簇尺寸扫描、空结对照等“控制/对照”计算。**
- 生产默认就是工具链默认：`pin_contact: false`、`tip_direction: axis`、`anchor_au_s_target: null`；不要为了“保险”切到 pinned 或 `lone_pair`。
- 如果结果异常、确实需要对照才能解释，先说明理由和代价，等用户明确同意再跑；不要静默加跑。
- 最终回复先给用户要求的 production 结果；对照只能作为建议，不能把已经跑了的对照混进“默认结果”。

## 判定口径

- `smj review` / `smj report` 的判定是机械判据，不替代科学裁决；物理判定 `inconclusive` 时不要包装成 `passed`。
- 接触几何/协议不一致时，必须在报告里显式标注，不得把不同 protocol 的 G 直接当作纯异构体比值。

## 代码

- **junction 内核在本仓维护**：算法改动直接改 `src/smj/junction/`（A 方案，2026-11-20）。不要再从 `agenticdft/agenticdft/tools/junction/` 覆盖同步，`scripts/sync_core.py` 已退役。
- 改过内核文件必须重跑 `python scripts/gen_junction_provenance.py` 刷新 `PROVENANCE.json` 的 sha256，并跑 `python scripts/gen_demo_reference.py --check` 确认数值基线没漂；与主仓的分叉用 `scripts/crosscheck_agenticdft.py` 对账，不进 CI（见 README §8）。
