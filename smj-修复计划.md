# smj-toolkit 0.1.0 修复计划（F1–F4）

> 依据：[smj-项目分析报告.md](smj-项目分析报告.md) 的 P0/P1 结论
> 目标版本：从 `0.1.0` 修复为本地 `0.1.1`（本地补丁，未回灌主仓前不对外发布）
> 状态：**待执行**（本文件先定稿，执行进度在 §7 回填）

---

## 0. 目标与范围

### 0.1 修复目标（4 项，按用户确认范围）

| 编号 | 问题 | 现在的后果 | 修好后的判定标准 |
|---|---|---|---|
| **F1** | `mos.aims` 从不回取（`head` 计划写了没人用） | 默认 `slurm-ssh` 六阶段链在第 5 阶段 `StageBlocked` | 集群跑一次 `scf → te` 全程自动通过；本地 `aims_scf/mos.aims` 存在且只需 2 行 |
| **F2** | 物理判定 CLI 断线 + 开启后 fail-open | 报告 `verdict` 恒为 `null`；空目录 run 被判"支持 H0" | 缺 `TE.dat` → `failed`（`te_present`）；`smj report` 必产出非 null `verdict`；温度取配置值 |
| **F3** | 状态机 5 处偏差（S1–S5） | `--wait` 抹掉 `collected`；`failed` 可绕过 `--force` 提交；缺件仍标 `collected` | 5 条回归测试全绿；`stage X --wait` 二次执行状态不变 |
| **F4** | A2 确认门默认不拦、CLI 无签署入口 | "红线"可被默认配置绕过；文档命令不存在 | `smj review --approve --who` 可签署；审计不过/无图拒绝（退出 2）；`require_a2` 默认 `true` |

### 0.2 明确不做（本轮范围外）

- 不修 D1–D12 的**全部**文档失实，只修被 F1–F4 直接影响的 4 条（D4 产物表、D8 判定说明、D9 的 A2 表述、README §7 命令表）。
- 不做 `jobwait` 的 `squeue` rc 细分（P1-4.7）——它需要真集群样本，留到下一轮。
- 不做 `landauer` 除零/`D→E` 等数值边界（P2）——单独一轮，避免与主逻辑混在一起。
- 不重构 vendored 代码架构；只做最小必要改动。

### 0.3 硬约束

1. **不放松 fail-closed**：任何"缺产物"路径只能从 `passed` 往 `failed/inconclusive` 收，不能反向放宽。
2. **不改公开接口语义**：`smj` 11 个子命令、`smj.yaml` 七段、`state.json` schema 不变；只新增可选参数。
3. **vendor 改动须可追溯**：`smj/_vendor/junction/*` 的每处改动都要在 `VENDOR.json` 记录"本地补丁 + 原文件 sha256"，否则回灌主仓时会静默丢补丁。
4. **Python 3.10+、仅标准库 + PyYAML + numpy**；A2 出图的可选依赖 matplotlib 不得变成硬依赖。
5. 所有工作在本工作区 `D:\dsh_projects\AutoSMJ` 内；不往系统 Python 或全局 `site-packages` 写。

---

## 1. 环境与基线（步骤 0，执行前一次性完成）

### 1.1 为什么要重建源码树

工作区只有 wheel（`smj_toolkit-0.1.0-py3-none-any.whl`），wheel 里**没有** `pyproject.toml`、`tests/`、`scripts/`、`VENDOR.json`。修复必须针对可编辑源码；直接装 wheel 只会得到 `site-packages` 里的只读副本。

### 1.2 步骤

```powershell
cd D:\dsh_projects\AutoSMJ

# 1) 从已解包并 compileall 通过的副本还原源码树
New-Item -ItemType Directory -Force src | Out-Null
Copy-Item -Recurse -Force "$env:TEMP\dsh-hxCm7F\smj_whl_analysis\src\smj" .\src\smj

# 2) 补 pyproject.toml（内容见 §1.3）、.gitignore、tests/ 骨架

# 3) 基线提交（任何修复都可 `git checkout` 回滚）
git init
git add -A
git commit -m "baseline: smj-toolkit 0.1.0 源码树（wheel 还原）+ pyproject + 分析报告"

# 4) 工作区 venv（复用系统已有 PyYAML 6.0.3 / numpy 2.5.1，不联网）
python -m venv --system-site-packages .venv
.\.venv\Scripts\python.exe -m pip install -e .\src --no-deps --no-build-isolation

# 5) 基线自检
.\.venv\Scripts\smj.exe verify
.\.venv\Scripts\smj.exe --help
```

### 1.3 `src/pyproject.toml`（唯一新增的打包文件）

> **执行时的偏差（2026-11-11）**：pyproject.toml 实际放在**仓库根**（不是 `src/`）：
> pytest 的 rootdir 与 `package-dir = {"" = "src"}` 都以根为基准更稳；`pip install -e .`
> 也改成在根执行。另外本机系统 Python 没有 setuptools、`python -m venv` 的 ensurepip
> 也失败，导致"可编辑安装"不可行，最终采用**等价方案**：改 `pyvenv.cfg` 的
> `include-system-site-packages = true`、写入 `.venv/Lib/site-packages/smj_src.pth`
> 指向 `D:\dsh_projects\AutoSMJ\src`、并用 `.venv/Scripts/smj.cmd` 作 console script
> 替身（`import smj.cli; main()`）。`smj verify` 通过、`smj.cmd plan` 正常。

```toml
[build-system]
requires = ["setuptools>=68"]
build-backend = "setuptools.build_meta"

[project]
name = "smj-toolkit"
version = "0.1.1"
description = "Single-molecule junction transport workflow (standalone)"
requires-python = ">=3.10"
dependencies = ["PyYAML>=6,<7", "numpy>=1.24"]

[project.optional-dependencies]
plot = ["matplotlib>=3.7"]
dev = ["pytest>=8"]

[project.scripts]
smj = "smj.cli:main"

[tool.setuptools]
package-dir = {"" = "src"}

[tool.setuptools.packages.find]
where = ["src"]
include = ["smj*"]

[tool.setuptools.package-data]
smj = ["templates/*"]

[tool.pytest.ini_options]
testpaths = ["tests"]
pythonpath = ["src"]
addopts = "-p no:cacheprovider"
```

### 1.4 基线冒烟脚本 `scripts/smoke.ps1`

把分析阶段跑通的 14 条命令固化（`init/plan/stage --prepare/run --dry-run/status/report/analyze/preflight/verify/...`），每条断言退出码与关键输出。**任何修复后必须重跑全套**，作为"没有回归"的证据。

---

## 2. F2 — 物理判定层（先做：改动最小、收益最大）

### 2.1 现状与根因

| 位置 | 问题 |
|---|---|
| `smj/_vendor/junction/criteria.py:286-287` | 无 `TE.dat` 返回 `[]` |
| `smj/_vendor/junction/criteria.py:197-203` | 空判据列表落到 `else: v = 'passed'` |
| `smj/report.py:347` | `include_checks: bool = False` |
| `smj/cli.py:620-631` | `cmd_report` 只传 `force`，无开关 |
| `smj/report.py:403-406` | 只判 `te_dir.is_dir()`，空目录算有谱 |
| `smj/_vendor/junction/validate.py:140-141`、`report.py:407-408`、`manifest.py:107` | 判据温度固定 300 K，不读配置 |

### 2.2 改动

**(a) `criteria.py` — 补一条"产物缺失"判据**

- 在 `CRITERIA`（`:31-41`）注册 `te_present`：`{'name': 'te_present', 'kind': 'engineering', 'gate': 'error'}`；
- `evaluate_run_criteria` 的 `if not te.exists(): return []` 改为
  ```python
  return [_check('te_present', None, False,
                 '缺 TE.dat：物理判据无法评估，不得视为通过')]
  ```
- `physics_verdict` 开头加防御：
  ```python
  if not checks:
      return {'physics_verdict': 'inconclusive', ...}
  ```

**(b) `report.py` + `cli.py` — 让判定真的产出**

- `cmd_report`：`report` 子命令新增 `--no-checks`（默认开启判定），`_invoke(func, config, force=args.force, include_checks=not args.no_checks)`；
- `report.py:403-406`：`te_dir.is_dir()` → `(te_dir / "TE.dat").is_file()`；
- `report.py` 组装 manifest 时把 `verdict` 写进 `manifest.json` 与 `report.md` 的固定小节（`md:167` 的设计意图）。

**(c) 温度透传**

- `report.py`、`manifest.py`、`validate.py` 三处调用 `validate_junction_run(...)` / `evaluate_run_criteria(...)` 时传 `temperature=config.transport.temperature`。
- `evaluate_run_criteria` 的默认值保留 300.0（兼容直接调用），但所有生产路径显式传参。

### 2.3 测试与验收

`tests/test_criteria_failopen.py`：

1. 空 run（只有空 `aitranss_te/`）→ `physics_verdict == 'failed'`，`failed_criteria == ['te_present']`；
2. 有 `TE.dat` 且三条阈值全过 → `passed`；
3. `TE.dat` 损坏 → `failed` 且 detail 含"解析失败"；
4. `smj report`（临时 run）→ `manifest.json["verdict"]` 非 `null`；
5. 温度：`transport.temperature=77` 时 `evaluate_run_criteria` 收到的温度是 77（用 monkeypatch 记录入参）。

**验收**：分析阶段那两个探针（空目录 → passed；`verdict` 恒 null）双双反转。

---

## 3. F1 — `mos.aims` 回取（打通主链）

### 3.1 根因链

`artifacts.py:23-25`（`must` 不含 `mos.aims`）→ `:33`（`AIMS_HEAD_ONLY`）→ `:51`（`head.update`）→ `remote.py:363-366`（`download_plan` 只返回 `must`）→ `download()` 从不消费 `head` → `tcontrol.py:194-195` 读不到 `mos.aims` → `aitranss_te.py:54-57` `StageBlocked`。

### 3.2 改动

**(a) `remote.py:368-392` — 实现表头回取**

```python
plan = _artifacts.retrieval_plan(stage)
heads = {} if names is not None else dict(plan.get("head") or {})
wanted = list(names) if names is not None else list(plan.get("must") or ())
# …现有 must 全量回取…
for name, n_lines in heads.items():
    res = self.transport.run(
        "head -n %d %s" % (int(n_lines), shlex.quote(remote_join(remote_dir, name))))
    if res.ok and res.out.strip():
        (target / name).write_text(res.out, encoding="utf-8")
        fetched.append(name)
    else:
        missing.append(name)
        errors[name] = (res.err or res.out).strip()[:300]
```

返回值增加 `"heads": list(heads)`，便于 `collect` 记录。

**(b) `stages/aims_scf.py:64-70` — 缺件即 fail-closed**

```python
if not has_mos:
    raise StageError("aims_scf 缺 mos.aims：aitranss_te 无法确定 $nsaos")
if not has_basis:
    raise StageError("aims_scf 缺 basis-indices.out：aitranss_te 的 $basis 不完整")
```

**(c) 头行数校验**：用真实 `mos.aims` 前若干行确认 `tcontrol` 解析所需行（疑似 nsaos 在第 1–3 行）；必要时把 `AIMS_HEAD_ONLY['mos.aims']` 从 2 调到实际值。

### 3.3 测试与验收

`tests/test_download_head.py`（FakeTransport 扩展一个 `on_run` 钩子返回 `head` 输出）：

1. `download('aims_scf', tmp)` → 本地生成 `mos.aims` 且仅 N 行；`fetched` 含 `mos.aims`；
2. 远端 `head` 失败 → `missing` 含它，且 `collect` 抛 `StageError`；
3. `names=[...]` 显式指定时不触发 head 逻辑（保持旧行为）。

**验收**：真集群（或本地 `local` 模式伪造目录）跑 `aims_scf → aitranss_te` 不再 `StageBlocked`；`tcontrol` 能从 2 行表头解析出 nsaos。

**风险**：若 `tcontrol` 需要更多行，增大 N（每行约几百字节，代价可忽略）。

---

## 4. F3 — 状态机 5 处（S1–S5）

### 4.1 改动（全部在 `smj/runner.py`，S3 可能涉及 `state.py`）

| 编号 | 位置 | 改法 |
|---|---|---|
| S1 | `runner.py:350-354` | `wait()` 入口：`collected`/`finished` → 幂等返回；仅 `submitted`/`waiting` 允许进入，其余 `StageError` |
| S2 | `runner.py:284-293` | 拦截集合补 `failed`：`("submitted","waiting","finished","collected","failed")` 均需 `--force` |
| S3 | `runner.py:498-519` / `state.py:208-218` | 统一 `skipped` 语义：`run()` 跳过 `skipped`（推荐），`is_done` 保持不变 |
| S4 | `runner.py:406-439` | 前置状态限 `{"submitted","finished"}`（`failed` 需 `--force`）；`missing` 非空 → `StageError`，不写 `collected` |
| S5 | `runner.py:256-266` | `prepare()` 一并清 `results` 与 `extra.download_*`；`--force` 用新记录整体替换 |

### 4.2 测试与验收

`tests/test_runner_states.py`（`FakeTransport` + 内存 state）：

1. S1：`collected` 阶段 `wait()` → 状态仍 `collected`，无异常；
2. S2：`failed` 阶段 `submit()` 无 `--force` → `StageError`；
3. S3：`skipped` 阶段 `run()` → 不执行；
4. S4：`download` 缺件 → `collect()` 抛 `StageError`；
5. S5：`prepare()` 后 `results == {}`。

**验收**：README:325 推荐的 `smj stage X --wait` 连跑两次，第二次不退状态；`smj run` 的"已 collected 跳过"分支恢复生效。

---

## 5. F4 — A2 确认门

### 5.1 改动

**(a) 新增 `smj/a2.py`**：封装"签署"——读 `review.json` → 校验工程审计 PASS、`a2.present`、附图（或显式 `--no-figure`）→ 调 `_vendor/junction/fingerprint.py:42` 的 `record_confirmation` 写 `a2_confirmation.json`；拒绝时返回 `(False, reason)`。

**(b) `cli.py:648-667` + `:947-948`**：`smj review` 增加 `--approve`、`--who NAME`、`--no-figure`；`--approve` 时不走只读路径，直接调 `a2.approve()`，失败退出 **2**（与 `a2_viz.py:450-459` 一致）。

**(c) 默认值**：`config.py:213` `require_a2` 默认改为 `True`；`smj init` 生成的模板显式写 `require_a2: true`，并保留 `false` 作为逃生门。

**(d) 文案**：`_vendor/junction/review.py:230-237` 的 `pending_message` 由 `agenticdft smj-review <run> --approve` 改为 `smj review --approve --who <name>`。

### 5.2 测试与验收

`tests/test_a2_approve.py`：

1. 审计 PASS + 有图 → 签署成功，`a2_confirmation.json` 存在且 `who` 正确；
2. 审计 FAIL → 拒绝，退出 2，提示不用 `--force`；
3. 无图且无 `--no-figure` → 拒绝，退出 2；
4. 默认 `require_a2=True` 时 `smj run --dry-run` 到 `au60_build` 前被 `StageBlocked`（回归"红线"语义）；
5. 旧配置显式 `require_a2: false` → 行为与 0.1.0 一致（兼容）。

---

## 6. 测试矩阵与新增文件

| 文件 | 覆盖 |
|---|---|
| `tests/conftest.py` | 临时 run 工厂、FakeTransport 增强（`on_run`） |
| `tests/test_criteria_failopen.py` | F2 的 5 条 |
| `tests/test_download_head.py` | F1 的 3 条 |
| `tests/test_runner_states.py` | F3 的 5 条 |
| `tests/test_a2_approve.py` | F4 的 5 条 |
| `tests/test_cli_contract.py` | 11 子命令存在、退出码、`--json` 稳定性（回归护栏） |

运行方式：`.\.venv\Scripts\python.exe -m pytest -q`（`pyproject` 已配 `pythonpath=src`）。

---

## 7. 执行顺序、里程碑与进度

| 阶段 | 内容 | 产出 | 状态 |
|---|---|---|---|
| M0 | 源码树 + pyproject + venv + git 基线 + 冒烟脚本 | 可运行基线，`smj verify` 通过 | ✅ 2026-11-11（`39a37b4`；pyproject 实际落在仓库根，见 §1.3 备注） |
| M1 | F2 判定层 | 5 测试绿，空目录不再 passed | ✅ `bc7febf`（7 条，`tests/test_criteria_failopen.py`） |
| M2 | F1 头回取 | 3 测试绿，链路打通（本地模拟） | ✅ `f7d3970`（4 条，`tests/test_remote_mos_head.py`） |
| M3 | F3 状态机 | 5 测试绿 | ✅ `fd52550`（8 条，`tests/test_state_machine.py`） |
| M4 | F4 A2 门 | 5 测试绿 | ✅ `fe48c29`（13 条，`tests/test_a2_approve.py`） |
| M5 | 文档修订 D4/D8/D9 + README §7 | 文档与代码一致 | ✅ 2026-11-11（README §1/§4/§6/§7/§8/§9；新增 §7.1 判定与 A2 红线） |
| M6 | 构建 0.1.1 wheel，干净 venv 安装，重跑冒烟 + 全部测试 | `dist/smj_toolkit-0.1.1-py3-none-any.whl` | ✅ 2026-11-11（238193 B / 69 条目；`pip install --target` 成功；从安装树 `python -m smj.cli verify` 六项全 ok；pytest 32 passed） |
| M7 | 生成 `VENDOR.json`（sha256 + 本地补丁清单） | 可追溯 | ✅ 2026-11-11（`scripts/gen_vendor_json.py`；36 文件、3 条补丁，随 wheel 一起分发） |

每完成一个里程碑：`git commit`，并在本表回填状态与验证证据（命令 + 关键输出）。

### 7.1 M6 的两次返工（供回迁主仓时避坑）

1. **模板路径错**：`smj verify` 的模板检查读的是 `smj/_vendor/junction/templates/`，
   但最初 pyproject 的 package-data 与 `build_wheel.py` 的 glob 都写成 `smj/templates/*`，
   结果第一次构建的 wheel 里**一个模板都没有**，从安装树跑 `verify` 报
   `templates FAIL 缺少: sub_g16.sh, sub_aims.sh, sub_te.sh`。修正：package-data 改为
   `"smj._vendor.junction" = ["templates/*", "VENDOR.json"]`，构建脚本加"没有
   templates/* 就报错"的守卫。源码树里 `verify` 先前能通过，只是因为读到了源码目录，
   掩盖了打包缺失 —— **打包验证必须从 wheel 解出的树里跑**。
2. **`__version__` 漂移**：`src/smj/__init__.py` 仍写 `0.1.0`，而 pyproject 已是
   `0.1.1`。构建脚本现在会比较两处版本，不一致直接拒绝构建。
3. **本机 pip 的沙箱限制**：`tempfile.mkdtemp()` 建出的目录（mode 0700）在 DSH
   Windows 沙箱里 ACL 不可访问（`WinError 5`，pytest 的 `tmp_path` 也因此不可用）。
   `pip install --target` 内部要用 mkdtemp 暂存，直接跑必失败；M6 的安装验证是在
   in-process 调用 pip、把 `tempfile.mkdtemp` 换成"默认 mode 的等价实现"后完成的。
   这是本机沙箱问题，不是 wheel 问题：普通机器上 `pip install dist\*.whl` 即可。

---

## 8. 验收标准（Definition of Done）

1. `pytest -q` 全绿（当前 32 条），`smj verify` 全通过；
2. 四个"反转探针"成立：
   - 空 `aitranss_te/` 的 run → `failed`（`te_present`），不再是"支持 H0"；
   - `smj report` 后 `manifest.json["verdict"]` 非 `null`；
   - 模拟 `aims_scf` 产物 → `aitranss_te.prepare` 不再 `StageBlocked`；
   - `smj review --approve` 在审计不过/无图时退出 2，通过时写出确认文件；
3. `smj stage X --wait` 对 `collected` 幂等；
4. 新 wheel 在干净 venv 中 `pip install` 后 `smj verify` 通过；
5. 文档：README §7 命令表补全 11 个、删 `orbital_match.json` 与 `mos.aims` 作为全量产物的表述、A2 表述与默认值一致。

---

## 9. 风险登记与回滚

| 风险 | 影响 | 缓解 |
|---|---|---|
| 直接改 `_vendor/junction/*` 与主仓漂移 | 下次 `sync_core` 覆盖补丁 | `VENDOR.json` 记录补丁 + sha256；条件允许时回灌主仓 |
| `require_a2` 默认改 `True` 破坏既有脚本 | 老 run 被 `StageBlocked` | `smj init` 模板写 `true`；配置显式 `false` 保持兼容；CHANGELOG 说明 |
| `head -n 2` 不足以解析 nsaos | `tcontrol` 仍失败 | 用真实文件校准行数；测试固定该行数 |
| 空 `TE.dat` 判据过严导致误报 | 部分体系被 `failed` | 判据只在"缺文件"时 error；有文件仍按原阈值 |
| 无主仓访问权 | 无法 `sync_core --check` | 本地 `VENDOR.json` 兜底，并在 README 声明"本地补丁分支" |

**回滚**：每个里程碑独立 commit；任一阶段失败 `git revert <commit>` 即可，M0 前的原始 wheel 保持不动。

---

## 10. 附：与主仓的同步约定

- 若后续拿到 `agenticdft` 主仓：把 `criteria.py`、`artifacts.py`、`remote.py`、`aims_scf.py`、`runner.py`、`review.py` 的补丁按原路径回灌，再跑 `scripts/sync_core.py` 重建 vendor 副本，最后 `sync_core.py --check` 应无差异。
- 回灌前 `VENDOR.json` 是唯一能证明"本地改过哪些文件"的台账，必须先生成。
