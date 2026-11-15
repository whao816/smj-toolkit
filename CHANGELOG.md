# Changelog

本项目遵循语义化版本（0.x 阶段：次版本可含不兼容变更；rc 视为发布候选）。

## 未发布

### 输运能量窗改 eV 相对费米能级（默认 E_F ±2 eV）

- `transport.energy_window` 的默认语义从 **Ha 绝对能量** 改为 **eV 且相对费米能级**：
  默认 `[-2.0, 0.005, 2.0]` = `E − E_F ∈ [−2, +2] eV`、步长 5 meV（801 点）。
  换金属/换分子后 E_F 移动，相对窗口不跟着漂，不再需要手动按 Au Fermi 估窗口。
- `aitranss_te.prepare` 现在先读 `aims_scf/aims.out` 的 efermi（`parse_efermi`），
  把 eV 相对窗口平移到 `tcontrol` 需要的 Ha 绝对 `$ener/$estep/$eend`
  （`HARTREE_EV = 27.211386245988`）；缺 efermi 直接 `StageBlocked`，不用占位值。
- 新增 `transport.energy_window_unit`（`eV` 默认 / `Hartree` 旧口径）。旧 yaml 里
  跨度 ≤0.1 Ha、量级 <2 Ha 的三元组按旧 Hartree 口径识别并给迁移 warning，
  行为与 0.2.0 逐字节一致；要显式保留旧口径可写 `energy_window_unit: Hartree`。
- `transport.ener/estep/eend` 仍是 **Ha 绝对能量** 的单分量覆盖（优先级最高）。
- `protocol_sha256`：eV 口径把用户写的三元组 + unit 入协议；旧 Hartree 口径不再新增
  key，历史 run_spec 不会仅因本改动被判漂移（无需重建 run）。
- 新增 `tests/test_energy_window.py`（12 条）：默认相对窗口、efermi 平移、旧 Ha 探测、
  unit 显式优先、分量覆盖、protocol 指纹形状、`aitranss_te.prepare` 真的把 eV 口径写进
  tcontrol、缺 efermi 阻断。
- 文档：README 输运段、`docs/smj-计算流程梳理.md`（总览/步骤 6/§3.4/§8/§10）、
  demo/azulene README（历史备注）同步；`smj init` 模板改 `[-2.0, 0.005, 2.0]` +
  `energy_window_unit: eV`。

### 连接核对（原 A2 / 表面核对）更名 + 默认自动登记

- 「A2 审核」→「表面核对」→**连接核对**（connection check）：它核对的是电极与分子之间的
  连接/表面三元组是否与结构一致，不再用任务卡内部代号。用户可见措辞（CLI help、提示、报告、文档）
  全部换成「连接核对」；`表面核对` 与 `A2` 都降为旧名兼容。
- **默认自动登记**：新增 `review.auto_connection_check`（默认 `true`）。`electrode_build` collect 之后，
  管线自动跑表面审计并写 `electrode_build/connection_check.json`（`who="auto"`）；
  审计不通过则**拒绝登记**，`smj review` 仍会 blocked。想要「必须人工看图再签」，
  设 `auto_connection_check: false`，仍可用 `smj review --approve --who NAME`。
- 配置键更名（旧键仍生效，加载时提示改名）：`review.require_a2` / `review.surface_check` →
  `review.connection_check`，`review.a2_deviation` / `review.surface_check_deviation` →
  `review.connection_check_deviation`，`review.auto_surface_check` → `review.auto_connection_check`。
- 产物文件更名：`a2_confirmation.json` / `surface_check.json` → `connection_check.json`，
  `a2_surface_audit.json` / `surface_check_audit.json` → `connection_check_audit.json`，
  `a2_surface_viz.png` / `surface_check_viz.png` → `connection_check_viz.png`；
  旧文件仍会被读取（签署记录、审计、review 闸门均兼容）。
- 复核码与 check 名更名：`a2_*` / `surface_check_*` → `connection_check_*`；`validate` 的
  `junction_a2_confirmation` / `junction_surface_check` → `junction_connection_check`（旧名兼容跳过）。
- CLI：`smj review --connection-check/--no-connection-check` 取代 `--surface-check` / `--a2`
  （旧开关保留为隐藏别名，只覆盖本次调用、不写回配置）。
- 新增 `tests/test_connection_check_approve.py`、`tests/test_cli_connection_check.py`、
  `tests/test_runner_connection_check.py`，并删除 `tests/test_cli_a2_switch.py`。

### 事后补登记判据族不再阻断复核（`family` 漂移容差）

- 背景：`molecule.family` 只决定判据适用族，但计入 `protocol_sha256`；run 完成后再补登记族
  （例如给 BPY 标 `family: BPY`）会让 `smj review` 报 `provenance_drift` 并 blocked。
- 现在复核侧做精确容差：把当前配置的 `family` 改成 run 时值重算 `protocol_sha256`，
  若与 run_spec 一致 → 记 `warning`（`provenance_family_relabel`），按当前标注重算判据，不阻断；
  其余任何物理参数差异仍为 `error`（`provenance_drift`）。`run_spec.json` 不重写，审计链不变。
- `tests/test_review_family_drift.py` 覆盖：仅 family 不同放行、tip 距离变化仍阻断、
  缺 config 快照 fail-closed。

### 判据族 `BPY`

- `SYSTEM_FAMILIES` 登记 `BPY: "Au-(4,4'-bipyridine)-Au"`（吡啶 N 锚定，非硫醇盐）。
  `t_ef_magnitude` / `g_in_theory_band` / `seebeck_sign` 的标定域是 Au–BDT–Au，
  对 BPY 记 **not applicable**（不再 fail-safe 误判「否证 H0」），只作参考对照
  （实验 1.3±0.1 MΩ ≈ 9.93e-3 G0）。
- `runs/bpy` 重新生成报告后：工程判定通过、物理判定「支持 H0（passed）」、可进入发布/预测卡。
- 测试：`tests/test_family.py` 新增族登记与「BDT 带不适用」用例；全套 `pytest` 162 passed。

### 电极构建参数化（`electrode_build` → 「构建电极」）

- 阶段 4 不再写死「本地构建 Au(111) 电极」：**怎么构建**与**构建什么电极**都是本阶段参数。
  - `electrode.source`：`library`（FHI-aims `electrodes.library` 预定义团簇，默认）/ `builtin`（本地三角岛）。
  - `electrode.element`：电极金属（fcc(111) 单晶表面），登记 Au / Ag / Cu / Ni / Pt / Pd / Al；
    **tip 与电极同元素**，`aims_tip_relax` 由该键决定加哪种金属 tip。
  - `electrode.a_lattice`：fcc 晶格常数（Å），`null` = 用内置登记值（Au 4.07823 … Al 4.04950）；
    未登记元素必须显式给出（2.0–6.0 Å）。
- `source: library` 时若库团簇元素与 `electrode.element` 不一致，`electrode_build` **fail-closed**
  （`StageError`），避免「Ag 电极配 Au 团簇」这类静默错配。
- `electrode_metadata.json` 新增 `electrode_element` / `electrode_a_lattice`，仅在非 Au 或缺省晶格
  被显式覆盖时写入；**Au 默认路径的 metadata 仍逐字节不变**，旧 run 缺键时下游按 `Au` 处理
  （`self.energy.in`、`tcontrol` 表面索引、`smj review` 表面审计、`validate` 表面/锚定/电极范围检查）。
- 规范检查与报告：`validate` 的表面、锚定、电极范围检查改读 metadata 的电极元素；
  2.40–2.49 Å 的 S–Au 有源带**只按 Au 标定**，非 Au 电极只报接触距离、不作偏离判定（severity info）。
  `smj report` 新增「电极构建方式 / 电极元素」行。
- 兼容性：`electrode.element`/`a_lattice` 进入 `protocol_sha256`，旧 run 的 provenance 对比会
  报 `provenance_drift`（预期；Au 默认值下重建输入才可复现旧哈希）。
- 测试：新增 `tests/test_electrode_metals.py`（配置校验、非 Au tip、builtin Ag/自定义晶格装配、
  library 元素一致性、Au metadata 字节兼容）；全套 `pytest` 162 passed。
- 文档：README、`docs/smj-计算流程梳理.md`（步骤 3/4 与配置表）、`docs/README.md` 同步改为
  「构建电极」口径。

### 无硫分子锚定（任意元素 / 去质子三态）

- 锚定原子不再限定硫醇盐的 S：`electrode.left_anchor` / `right_anchor` 可写任意元素符号
  （结构里恰好出现 2 次）或 1-based 序号，分子里没有 S 也能建结（N/Se/… 锚）。
  `aims_tip_relax` 不再报「需恰好 2 个 S 锚原子」。
- `electrode.strip_thiol_h` 改为**三态**：`null`（新默认）= 仅当锚定原子是 S 时才自动移除
  S–H（旧 S 分子的行为不变）；`true` = 强制去质子；`false` = 不去。非 S 锚默认保留 H。
- 去质子 fail-closed：每个锚原子必须恰好 0/1 个 H、总计 0/2 个才移除；-NH2 这类一个锚原子
  带 2 个 H 的情况直接报错（不猜着只去一个）。
- `aims_tip_relax/input_summary.json` 新增 `anchor_element`（锚定元素），`strip_thiol_h` 记录
  实际去质子信息（`anchor_element`/`anchor_indices`/`route` 等）；注释与错误信息里的
  「Au–S」推广为「锚–金属」，tip/锚元素分别来自 `electrode.element` 与结构。
- `validate` 的 `junction_geom_anchor_aus` 改用 metadata 的 linker 索引识别锚定原子（任意元素）：
  S 锚保持原有硫醇盐带/硫醚分流逻辑；非 S 锚 severity info，只报告锚–金属接触距离，
  2.40–2.49 Å 带不适用。（旧 metadata 无 linker 索引时退回按 S 找锚。）
- 兼容性：`strip_thiol_h` 默认值由 `true` 变为 `null`，进入 `protocol_sha256`，旧 run 的
  provenance 对比会报 `provenance_drift`（预期；S 分子的实际去质子行为不变）。
- 测试：新增 `tests/test_no_sulfur_anchor.py`（三态解析、无 S 锚按元素/序号解析、默认保留 H、
  无 S 结校验给 info、S 分子仍自动去质子、-NH2 强制去质子 fail-closed、CLI prepare 端到端）；
  全套 `pytest` 147 passed。
- 文档：README（输入/适用/口径）、`docs/smj-计算流程梳理.md`（步骤 3 参数表、§3.3、判据表、
  §10 第 17 条）同步说明任意元素锚定与三态去质子。

## 0.2.0 — 2026-11-15

与 0.2.0rc1 相同的代码；本版补齐发布记录与真实体系验证。

### 验证（真实集群，FHI-aims 221103 + aitranss 071813）

- 三套体系跑通全链（molecule_opt → tip_relax → electrode_build → SCF → TE → Landauer）：
  BDT（86 原子，L0=0.279）、azulene-2,6-dithiol（92 原子，L0=0.381）、
  2CN-656-4-thiomethylphenyl（硫醚，128 原子，L0=0.00132）；工程判定全部通过，
  物理判定均为 inconclusive（`t_ef_no_peak_near_ef`）。详见 `docs/validation-0.2.0.md`。
- 网格收敛对照：BDT 与「冻结 S–Au=2.42 Å」两套几何上，estep 0.0002→0.00005 Ha
  （301→1201 点）后 L0/S/k_e/ZT_e 相对差 < 1e-8；`runs/_compare/BDT-对照实验报告.md`。
- 新增 `RELEASING.md`（版本 bump → 测试 → wheel/sdist → 安装冒烟 → tag/CI → Release）；
  README 增加「已验证体系」小节。

## 0.2.0rc1 — 2026-11-14

### 打包与发布工程

- 新增 MIT LICENSE，并写入 `pyproject.toml` 与 wheel METADATA（PEP 639 `License-Expression`）。
- 补全分发元数据：`readme` 长描述、authors、keywords、classifiers（Python 3.10–3.12 / Linux / Science）。
- `scripts/build_wheel.py` 支持 `--sdist`：同时产出可复现的 wheel 与 sdist（固定时间戳、排序、RECORD sha256 自检）。
- 版本号升至 0.2.0rc1；`smj verify` 在从 wheel 解包的环境下通过。
- 新增 GitHub Actions CI（`.github/workflows/ci.yml`）：pytest 矩阵 Python 3.10–3.14（Ubuntu）与 3.12（Windows），
  以及 wheel/sdist 构建 + 干净 venv 安装冒烟（`smj verify` / `smj init` / `smj plan`）。

### 验证

- 用户视角安装验证（Ubuntu 26.04 / Python 3.14.4，超出原 CI 矩阵）：wheel 与 sdist 分别装入干净 venv，
  `smj verify` 六项通过、`smj init` + `smj plan` 输出六阶段；`pip show` 正确解析 METADATA 2.4
  （`License-Expression: MIT`），`pip check` 无破损依赖；PyYAML 6.0.3 / numpy 2.5.4 均有 cp314 wheel，
  安装过程未触发本地编译。
- 本地全量 `pytest`（125 passed）在 Python 3.14.4 上通过；package classifiers 已加入 3.13/3.14，
  CI 矩阵同时覆盖这两个版本（3.13 未在本机单独跑过，由 CI 观察项兜底）。

### 默认行为（可能影响既有配置）

- `transport.energy_window` 的代码默认由 `(-2.5, 0.02, 2.5)` Ha 改为
  `(-0.21, 0.0002, -0.17)` Ha（围绕 Au Fermi ≈ −0.19 Ha，±0.02 Ha，201 点）。
  与 `smj init` 模板一致；生产仍应在 `aims_scf` 后按实际 E_F 调整窗口。
- `electrode.au_per_electrode` 在 `source: library`（默认）下被完全忽略，团簇大小由
  `electrodes.library` 文件决定；若显式写了该键会给出 warning。`source: builtin` 时
  `null` = 完整 60 Au/侧，整数按完整原子层截断（30 → 35/侧），上限 30。
- `smj preflight` 在 `molecule.opt_engine: aims` 时跳过 Gaussian 路径/许可检查
  （此前 G16ROOT/G16_MODULE 为空会误报 error）；改为 info 提示 `gaussian_not_required`。

### 自 0.1.2 以来的功能与修复

- 阶段规范化命名：`molecule_opt`（旧 `gaussian_opt` 别名）、`electrode_build`（旧
  `au60_build`/`auN_build` 归一化）；旧 run 目录仍可回退读取。
- 分子优化默认改用 FHI-aims（`opt_engine: aims`，PBE/light），与 tip/SCF 同一软件栈；
  Gaussian B3LYP/6-31G(d) 分支保留兼容。
- 默认使用 FHI-aims `electrodes.library` 预定义团簇（`source: library`），
  `electrode.nlayers: null` 取库文件头部建议；`source: builtin` 保留自建三角岛。
- 新增 `smj review`（A2 结构确认）：表面审计、渲染图、`geometry.in` sha256 指纹绑定；
  `review.require_a2` 默认 true，`--accept-deviation` 单次豁免，`--no-figure` 跳过出图。
- 新增 `smj review --a2 / --no-a2` 命令行临时开关（不写回配置、不进 `run_spec`）。
- 锚定化学判定：2.40–2.49 Å 只对去质子硫醇盐适用；硫醚/硫醇按 vdW 接触上限判定。
- 表面三元组确定性规则 det-v1（0.05 Å 分桶、坐标字典序），同一几何不再因运行次序换索引。
- 可复现闭环：`run_spec.json` 记录完整解析配置与协议哈希 `protocol_sha256`，漂移默认拒绝；
  各阶段记录输入/输出 sha256 与父阶段指纹。
- fail-closed 判据：缺 `TE.dat` 记为 `te_present(missing=True)` 并判 inconclusive，不再静默通过。
- 报告 `manifest.verdict` 区分 passed / failed / inconclusive，并给出 `allowed_action`。
- `smj analyze` 增加独立解析器交叉校验（Simpson + CODATA），相对差 ~1e-8。
- 修复 `sub_g16.sh` 不建 scratch 目录、Gaussian 坐标块缺空行、SCF 不收敛兜底
  （`scf=(xqc,maxcycle=512) int=ultrafine`）。
- 修复 mos/omat 回取清单只对 `aims_scf` 生效（`AIMS_HEAD_STAGES`）。

## 0.1.2 — 2026-11-11

- 库电极（`electrodes.library`）+ 背层固定弛豫（`au_relax`，主链路可选）。
- `smj report` / `smj analyze` / `smj audit` / `smj verify` 等命令完善。

## 0.1.1 — 2026-11-11

- 独立 wheel 打包、`smj verify` 安装自检、vendored junction 本地补丁（F1/F2/F4）。

## 0.1.0 — 2026-11-10

- 首个独立包：六阶段编排、状态机、Slurm/SSH 传输层、Landauer 后处理。
