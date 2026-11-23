# 验证数据（2026-11-22）

本目录是下列两份验证记录所引用 run 的**精简导出**，用于在没有 `runs/` 的
机器上核对文档里的数字，也是这些数据进入版本库的唯一途径（`runs/` 被
`.gitignore` 忽略）。

- [light vs tight 整链对照](../validation-light-vs-tight-2026-11-22.md)
- [vdW 开关整链对照](../validation-vdw-onoff-2026-11-22.md)

## 内容与范围

每个臂保留：`state.json`、`run_spec.json`、各阶段的 `input_summary.json`、
`control.in`、输入/优化几何、`saved_structures/`、`aitranss_te/TE.dat` 与
`junction_analyze/results.json`。文档里的 L0 / S / κe / ZT 可直接由
`TE.dat` + `results.json` 复算。

**不包含**：`aims.out`（单臂 10–25 MB）、PDF/PNG 渲染和对照图（图在
[`../figures/`](../figures/)）。本目录也不含未跑完的
`az-thiolate-*` 臂——它们已移到 `runs/light-tight-2026-11-22/_abandoned/`，
状态是 SCF 不收敛、无 `results.json`，不属于对照数据。

## 导出与校验

```bash
python scripts/export_validation_data.py           # 需要本机有 runs/，重新导出
python scripts/export_validation_data.py --check   # 校验 MANIFEST 与入库副本
```

`--check` 先核对入库副本的 sha256/字节数，再在 `runs/` 存在时逐字节比对
源文件；两边一致才算通过。`MANIFEST.json` 记录每个文件的来源路径与哈希，
`docs/` 其余文件请见 [文档索引](../README.md)。
