# 发布流程（RELEASING）

本文件是 `smj-toolkit` 的发布清单。发布前先确认 `docs/release-0.2.0-plan.md` 对应版本的
DoD 已满足；本文以 `vX.Y.Z` 代指目标版本。

## 0. 前置

- 工作树干净（`git status`）。
- Python ≥ 3.10 环境可跑测试：`python -m pytest -q`（当前基线 620 passed / 2 deselected；
  另有 2 条依赖 `runs/` 资产的集成回归用 `python -m pytest -m integration -q`）。
- 能访问 FHI-aims / aitranss 的验证环境（集群），或至少能完成 `smj verify` + `init` + `plan` 冒烟。
- 若从沙箱内操作：git-over-ssh 与 GitHub 凭证不可用，**push / tag 推送须在用户自己的终端执行**。

## 1. 版本号与变更日志

同时改三处，保持一致：

1. `pyproject.toml` → `[project] version`
2. `src/smj/__init__.py` → `__version__`
3. `CHANGELOG.md` → 新增/更新该版本段落（含默认行为变化与验证结果）

## 2. 测试

```bash
python -m pytest -q
```

要求全绿。若改动涉及 `smj/junction/` 内核，必须先跑
`python scripts/gen_junction_provenance.py` 刷新 `PROVENANCE.json`，并跑
`python scripts/gen_demo_reference.py --check` 确认数值基线没漂
（与 AgenticDFT 的对比用 `scripts/crosscheck_agenticdft.py`，不进 CI）。

## 3. 构建发布产物

```bash
python scripts/gen_junction_provenance.py          # 刷新 PROVENANCE.json（含 repo_commit）
python scripts/build_wheel.py --check --sdist
```

- 产物：`dist/smj_toolkit-<ver>-py3-none-any.whl`、`dist/smj_toolkit-<ver>.tar.gz`。
- `--check` 会按 RECORD 校验 sha256 与大小；sdist 使用固定时间戳（可复现）。

## 4. 安装冒烟

```bash
python -m venv .tmp/venv-release
# Windows: .tmp\venv-release\Scripts\activate ；Linux: source .tmp/venv-release/bin/activate
python -m pip install dist/smj_toolkit-<ver>-py3-none-any.whl
smj --version           # 应为 <ver>
smj verify              # 七项 ok
mkdir .tmp/smoke && cd .tmp/smoke
smj init --name smoke --xyz ../molecule.xyz
smj plan
```

再用 sdist 做一次同样冒烟（验证源码包可独立重建）：

```bash
python -m pip install dist/smj_toolkit-<ver>.tar.gz
```

> Windows + DSH 沙箱内 `pip --target` / `venv ensurepip` 可能因 WinError 5 失败，
> 这是环境限制；换 Linux/WSL 或 CI 完成该步骤。

## 5. 提交与标签

```bash
git add -A
git commit -m "release: <ver>"
git tag -a v<ver> -m "smj-toolkit <ver>"
git push origin master --tags      # 在用户终端执行
```

推送后确认 GitHub Actions `ci.yml` 全绿（pytest 3.10–3.14 × Ubuntu、3.12–3.14 × Windows、
wheel/sdist 构建 + 安装冒烟）。**CI 未绿不得发布。**

## 6. GitHub Release

在 GitHub → Releases → Draft a new release：

- Tag：`v<ver>`（已推送）
- 标题：`smj-toolkit <ver>`
- 附件：`dist/` 下的 wheel 与 sdist
- 说明：引用 `CHANGELOG.md` 对应段落 + `docs/validation-<ver>.md` 的验证结论

## 7. PyPI（可选，本轮默认不做）

```bash
python -m twine upload dist/smj_toolkit-<ver>-py3-none-any.whl dist/smj_toolkit-<ver>.tar.gz
```

- 需要 PyPI API token（`~/.pypirc` 或 `TWINE_USERNAME=__token__`）。
- rc/预发布：`twine upload --repository testpypi ...` 先走 TestPyPI。
- 上传顺序无关；PyPI 不可覆盖同名文件，发错版本只能 yank + 递增版本号。

## 8. 发布后

- 在 README 更新「已验证体系」链接与安装示例版本号（若示例里写死了版本）。
- 记录发布产物 sha256（`Get-FileHash` / `sha256sum`）到本文件或 CHANGELOG。
- 下一个版本的开发线：从 `master` 开 `feature/*` 分支；默认行为变更写进 CHANGELOG 的
  “默认行为（可能影响既有配置）”小节。
