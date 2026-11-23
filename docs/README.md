# 文档索引

## 当前参考

- [计算流程梳理](smj-计算流程梳理.md) —— 六阶段流程、逐步命令与默认参数、验收判据、资源模型、已知限制（0.2.x 口径；步骤 4 已参数化为「构建电极」：source / element / a_lattice；步骤 3 锚定支持任意元素，strip_thiol_h 三态；复核门含连接核对（原 A2），默认自动登记）
- [科学性要求与判定标准调研](smj-科学性要求与判定标准调研.md) —— 四层判据（实现/数值/物理/外部有效）、分阶段必要条件、量化阈值与文献证据
- [判据来源核验](criteria-source-research.md) —— `au_s_bond`、共振归因、`g_in_theory_band` 的一手来源核验与处置记录
- [计算 benchmark 表](benchmark-research-computational.md) —— 12 个文献计算行（BDT/BDA/BDCN/BP/HDA/C60/SMePh 等）+ smj 0.2.0 实测对齐口径
- [实验 benchmark 表](benchmark-research-experimental.md) —— 单分子结电导/热电势的实验值、分布与可比性说明
- [0.2.0 验证记录](validation-0.2.0.md) —— Slurm 三体系真实运行数据
- [0.2.2rc1 验证记录](validation-0.2.2rc1.md) —— 复核修复后的可验证状态：397 项测试、七项自检、demo 数值基线与产物 sha256
- [AutoSMJ 自动化流程中的 DFT 方法（v1.2）](AutoSMJ_DFT方法说明_v1.2.md) —— 主链 DFT 设置（FHI-aims 221103 / PBE / light / TS-vdW / scalar ZORA；Aitranss NEGF σ=0.05/0.10 Ha、300 K）与程序实跑过的算例（59 条历史 run + 2026-11-22 新臂）；§10 列出未验证边界
- [light vs tight 整链输运对照（2026-11-22）](validation-light-vs-tight-2026-11-22.md) —— AZ-dithiol/2CN 四臂整链：tight 相对 light 的 L0/S/ZT 差异；2CN 的 S/ZT 对档位高度敏感（light 档不可定量）
- [vdW 开关整链输运对照（2026-11-22）](validation-vdw-onoff-2026-11-22.md) —— 只切 TS-vdW 的两臂整链差异；差异来自弛豫几何，方向/幅度体系相关，不可外推
- [验证数据入库副本](validation-data/README.md) —— 上述两组对照六臂的精简导出（state/run_spec/几何/TE.dat/results.json + MANIFEST sha256；`runs/` 的入库快照，`aims.out` 不入库）
- [旧 run 重评记录（2026-11-16）](runs-reverdict-2026-11-16.md) —— 当前判据口径下全部既有 run 的 manifest/report 重算、新旧判定对照与 fail-safe 说明
- [外部评审报告的逐条复核（2026-11-20）](review-recheck-2026-11-20.md) —— 对 `0f7b3c4` 评审报告 7 条问题 / 3 条后续 / 验证范围声明的核实：成立与不成立、证据与复现、原报告说轻或说重之处，以及复核新增发现
- [复核问题修复计划（2026-11-20）](fix-plan-2026-11-20.md) —— 依据复核报告落地的 13 项修复（R1–R13）：分层（上游同步 / 本包 wrapper）、阶段与依赖、逐项验收标准、测试清单、兼容性与需要拍板的 5 个决策点
- [A 方案实施计划（2026-11-20）](split-plan-2026-11-20.md) —— junction 内核改由本仓维护：目录与命名去「副本」化、同步机制退役、结果契约（demo 基线 + 对账脚本）替代代码同源、R12 收口
- [自由基单分子结：实现对照与开发规划](自由基单分子结-实现对照与开发规划.md) —— 把外部调研文档翻译成 smj 口径：自旋/电荷/输运链逐条对照、口径差异、缺口清单与后续工作优先级
- [自由基支线改进计划（2026-11-21）](plan-radical-2026-11-21.md) —— 在既有自旋链上补证据与闸门：Hirshfeld 自旋分区（S1）、自旋分辨热输运（S2）、单列 TE 硬闸门（S3）、BS 局部证据（S4）、口径与一致性检查（S5）、真实 run 集成回归（S6）、文档（S7），含批次、决策点与验收

## 历史档案（archive/）

- [AutoSMJ_DFT方法说明 v1.0](archive/AutoSMJ_DFT方法说明_v1.0.md)、[v1.1](archive/AutoSMJ_DFT方法说明_v1.1.md) —— 方法说明的早期版本（v1.2 见上）
- [项目分析报告](archive/smj-项目分析报告.md)
- [修复计划](archive/smj-修复计划.md)
- [0.1.2 计划](archive/smj-0.1.2计划.md)
- [独立打包设计](archive/smj-standalone-package.md)

根目录保留：README.md（项目说明与安装）、CHANGELOG.md、LICENSE、pyproject.toml。
