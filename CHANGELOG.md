# Changelog

本项目遵循语义化版本（0.x 阶段：次版本可含不兼容变更；rc 视为发布候选）。

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
