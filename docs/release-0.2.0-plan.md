# 0.2.0 发布计划

> 状态：执行中（批准于 2026-11-14，方案 = 并行编排；渠道 = 只上 GitHub；A2 = 看图后签署）
> 基线：`v0.2.0rc1`（commit `efc9544`，wheel 288,263 B / sdist 278,890 B，pytest 125 passed）

## 1. 目标与完成标准（DoD）

1. 在 0.2.0rc1 代码上跑通 **2CN-656-4-thiomethylphenyl（54 原子）** 与 **azulene-26,2-dithiol（20 原子）** 全链（molecule_opt → aims_tip_relax → electrode_build → A2 → aims_scf → aitranss_te → junction_analyze），产出 `analysis.json`、`report.md` 与数值记录（`docs/validation-0.2.0.md`）。
2. `v0.2.0` tag 指向新 commit；GitHub CI（Ubuntu 3.10–3.14 + Windows 3.12）全绿。
3. `dist/` 产出 wheel+sdist，解包 `smj verify` 通过；WSL Ubuntu 3.14 用户视角安装复验通过。
4. `pytest` ≥125 passed；README 增加「已验证体系」；`RELEASING.md` 落地。
5. 发布渠道：**只上 GitHub**（tag + Release + 附件）；本轮不上 PyPI。

## 2. 编排（并行）

| 时间 | 事项 | 负责人/依赖 |
|---|---|---|
| Day 1 上午 | 2CN、AZ 配置与 preflight；并行提交 molecule_opt→electrode_build | 集群 `gpu-outside`，两个作业可并行排队 |
| Day 1 上午 | 推 GitHub（需 repo URL）→ CI 首跑 | 用户提供 repo URL |
| Day 1 下午 | 两条链 A2 出图 → 用户签署 → 提交 aims_scf→junction_analyze | 用户签名（whos = whao） |
| Day 1 晚 | 收集结果、`smj report`、记录数值；修暴露的 bug（预留半天） | — |
| Day 2 上午 | bump 0.2.0 → 构建 → RECORD/解包/`pip --root`/WSL 复验 → RELEASING.md | — |
| Day 2 下午 | tag `v0.2.0` → 推送 → CI 绿 → GitHub Release + 附件 | repo URL |

## 3. 验证细节

| 体系 | 原子数（分子/tip/全结） | 关键配置 | 预估墙钟 |
|---|---|---|---|
| 2CN-656-4-thiomethylphenyl | 54 / 56 / ~130 | family 2CN、S 锚、`strip_thiol_h:false`（硫醚）、`anchor_au_s_target:2.40`、库团簇 37 原子/侧 | SCF 1–1.5 h |
| azulene-26,2-dithiol | 20 / 22 / ~96 | family AZ、S 锚、`strip_thiol_h:true`、`anchor_au_s_target:2.40` | SCF 40–60 min |

统一：`mode: slurm-ssh`、`remote_root /home/haow/smj-runs`、`energy_window [-0.22, 0.0002, -0.16] Ha`、`require_te_out:true`、`review.require_a2:true`。

验收：两链均到 `junction_analyze` 且 `analysis.json` 有 `independent_check.ok`；`manifest.verdict` 可计算；无 fail-closed 中断。数值记入 `docs/validation-0.2.0.md`（T(E_F)、G/G0、S、k_e、ZT_e、峰位、Au–S、SCF 墙钟）。

## 4. 发布工程

1. 版本：`pyproject.toml` / `src/smj/__init__.py` → `0.2.0`；CHANGELOG 增 0.2.0 节（自 rc1 的验证与修复；若验证暴露并修复 bug，逐条列出）。
2. 构建：`python scripts/build_wheel.py --check --sdist`；记录 wheel/sdist 大小与 sha256。
3. 复验：解包 wheel → `smj verify`；`pip install --root`；WSL 3.14 venv 从 dist 安装 → verify/init/plan；`pip check`。
4. `RELEASING.md`：bump → pytest → build → verify → commit → tag → push → CI 绿 → GitHub Release（twine 步骤标注为"仅上 PyPI 时"）。
5. tag `v0.2.0`（annotated）→ push `master` + tags → CI 绿 → Release 附件上传 wheel/sdist。

## 5. 风险与回滚

- **Q 队列**：两作业并行占 8 核 + 1 GPU ×2；若排队久，AZ 先跑（96 原子便宜），2CN 后补。
- **资源上限**：130 原子 SCF 在 24G 下若 OOM，升 mem 到 48G（节点 124G）并重提交；`--time` 已按 ×1.9 余量。
- **bug**：任何 fail-closed 中断都先修代码+测试，再重跑受影响阶段；0.2.0 版本号复用（rc1 未公开）。
- **回滚**：PyPI 本轮不做，无需 yank；GitHub tag 在 push 前可移动；Release 可删除重发。
- **A2 偏差**：若表面索引与审计偏差超 0.15，签前讨论（`--accept-deviation` 仅作显式豁免并在报告中标注）。
