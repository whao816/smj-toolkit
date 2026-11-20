# 从半镜像到自有内核：A 方案实施计划（2026-11-20）

## 0. 背景与判据

复核修复（`docs/fix-plan-2026-11-20.md`）完成后，本包与 AgenticDFT 的关系已从
「镜像」退化为「半镜像」：

- 运行时本来就没有耦合：`src/smj` 里 `import agenticdft` 为 **0 处**；
- 但维护模型还假装是镜像：13 条本地补丁覆盖 **17/36 个文件（47%）**，
  `PROVENANCE.json` 的 `upstream_commit` 是 **None**，CI 里没有 `sync_core`，
  6 个文件曾从上游**未提交**工作区同步、`molecules.py` 反向落后。

用户决策（m01398）：**走 A** —— 本仓成为 junction 内核的主开发位置，
与 AgenticDFT 的关系从「代码同源」改为「结果同契约」。

## 1. 目标（可验收）

1. 目录与命名不再声称「副本」：`src/smj/junction/` → `src/smj/junction/`，
   `PROVENANCE.json` → `PROVENANCE.json`（记录来源与本地变更史，不再当镜像台账）。
2. 删除 `scripts/sync_core.py`（漂移守卫随之退役），README/AGENTS 的约定改为
   「内核在本仓维护」。
3. 与 AgenticDFT 的唯一联系是**可验证的结果契约**：
   `demo/demo-reference.json`（本仓回归基线）+ 可选对账脚本，不进 CI 的强依赖。
4. 消掉 R12「已知偏离」：`k_e_zt_order` 元数据改回 `implemented=True`（评估器实际会 emit），
   删除 PROVENANCE.json note 与 README §9 的登记条目。
5. 全量测试、`smj verify`、`gen_*_reference --check` 全绿。

## 2. 提交拆分

| 提交 | 内容 | 验收 |
|---|---|---|
| S1 | `docs/split-plan-2026-11-20.md`（本文件） | 无代码改动 |
| S2 | 目录迁移 + 全仓 import/路径引用改写（58 个文件）+ `_vendor` 包删除 | `pytest` 全绿；`smj verify` 全 ok |
| S3 | 清单改名与语义：`PROVENANCE.json`、`scripts/gen_junction_provenance.py`、`smj verify` 检查项改名；删除 `sync_core.py` | `gen_junction_provenance.py --check` exit 0 |
| S4 | 约定文档：AGENTS.md、README §1/§8/§9、CHANGELOG「0.2.2：内核自维护」 | 文档一致 |
| S5 | R12 收口：`k_e_zt_order` `implemented=True`，删除偏离登记 | 判据元数据自洽用例 |

## 3. 风险与回滚

- **风险 1：import 面广（58 个文件）**。用脚本做机械替换（`smj.junction` →
  `smj.junction`、`.junction` → `.junction`、路径字符串同步），
  以全量测试 + `smj verify` 兜底；S2 是单独提交，回滚只需 `git revert`。
- **风险 2：打包**。`pyproject` 用 `packages.find include=["smj*"]`，
  新目录 `smj.junction` 自动包含；模板随包数据目录一起搬，`artifacts`/`verify`
  读的是 `Path(__file__).parent / "junction" / "templates"`，随代码同步改。
- **风险 3：与 AgenticDFT 结果分叉**。用 `demo-reference.json` 与
  `scripts/crosscheck_agenticdft.py`（可选、手动跑）对账；CI 不依赖兄弟仓。
- 回滚：三个代码提交各自独立，`git revert` 任意一个都能回到上一状态。

## 4. 不做的事

- 不升级 `criteria.py` 到主仓的新注册表（那是行为变更，需单独评审）；
  R12 只改元数据，使「注册表声明」与「评估器行为」一致。
- 不把对账脚本接进 CI（CI 拿不到 AgenticDFT，接了就是假绿）。
- 不改 vendored 判定语义，不动 demo 数据。

## 5. 后续（本计划之外）

- R11（结构与机理证据）仍按原决策待办；
- 主仓若仍活跃，可反向 vendor `smj.junction` 或固定版本引用，
  由主仓自己决定，不再由本仓承担同步职责。
