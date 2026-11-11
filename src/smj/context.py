"""阶段运行上下文：配置 + 目录布局 + 状态 + 传输层 + 资源估算。

每个 :class:`~smj.stages.base.Stage` 在 ``prepare`` / ``collect`` 里只通过这个对象
访问外部世界，因此：

- 单元测试可以构造 ``StageContext``（``dry_run=True``，无传输层）直接调用阶段；
- 阶段独立运行时，上游产物可用 ``overrides`` 显式替换，不要求整条链重跑；
- ``stage_env()`` 把配置翻译成 vendored 输入写出器认识的 ``env`` 字典
  （``g16_mem`` / ``aims_bin`` / ``te_cpus`` …），数值语义与 AgenticDFT 引擎一致。
"""

from __future__ import annotations

import os
import shlex
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Sequence

from .config import RunConfig
from .errors import SmjError, StageBlocked
from .resources import estimate_smj_resources
from .state import (RunLayout, RunState, normalize_stage, previous_stage,
                    stage_index)
from .transport import CommandResult, Transport, remote_join


@dataclass
class StageContext:
    """一次阶段执行所需的全部外部依赖。"""

    config: RunConfig
    layout: RunLayout
    state: RunState
    transport: Transport | None = None
    dry_run: bool = False
    overrides: dict[str, Any] = field(default_factory=dict)
    logger: Callable[[str], None] = print
    extra: dict[str, Any] = field(default_factory=dict)

    # -- 基本 -------------------------------------------------------------- #
    def log(self, message: str) -> None:
        self.logger(str(message))

    def resolve(self, value: str | Path | None) -> Path | None:
        return self.config.resolve(value)

    @property
    def cluster_remote(self) -> bool:
        return self.config.mode == "slurm-ssh"

    def stage_dir(self, stage: str) -> Path:
        return self.layout.stage_dir(stage)

    def ensure_stage_dir(self, stage: str) -> Path:
        return self.layout.ensure_stage(stage)

    def save_state(self) -> None:
        self.state.save(self.layout.state_path)

    # -- 输入发现 ---------------------------------------------------------- #
    def molecule_xyz(self, *, required: bool = True) -> Path | None:
        """gaussian_opt 的分子结构：``--xyz`` 覆盖 > ``molecule.xyz``。"""
        override = self.overrides.get("xyz")
        if override:
            # CLI 覆盖路径按当前工作目录解析（配置文件里的路径才按 yaml 目录解析）。
            path = Path(str(override)).expanduser()
            path = path if path.is_absolute() else path.resolve()
        else:
            path = self.config.molecule_xyz()
        if path is None:
            if required:
                raise StageBlocked(
                    "gaussian_opt",
                    "未指定分子结构：请在 smj.yaml 设 molecule.xyz，或传 --xyz <文件>")
            return None
        if required and not path.exists():
            raise StageBlocked("gaussian_opt", "分子结构文件不存在: %s" % path, [str(path)])
        return path

    def artifact(self, stage: str, *parts: str) -> Path:
        return self.layout.artifact(stage, *parts)

    def require_stage_file(self, stage: str, filename: str, *,
                           target: str | None = None) -> Path:
        path = self.layout.artifact(stage, filename)
        if not path.exists():
            raise StageBlocked(target or stage,
                               "%s 阶段产物缺失: %s" % (stage, filename), [str(path)])
        return path

    def require_stage_files(self, stage: str, filenames: Sequence[str], *,
                            target: str | None = None) -> Path:
        missing = [str(self.layout.artifact(stage, name)) for name in filenames
                   if not self.layout.artifact(stage, name).exists()]
        if missing:
            raise StageBlocked(target or stage,
                               "%s 阶段产物不完整（需要 %s）" % (stage, ", ".join(filenames)),
                               missing)
        return self.layout.stage_dir(stage)

    def upstream(self, stage: str, filename: str, *,
                 override: str | Path | None = None) -> Path:
        """目标阶段默认消费**上一阶段**的同名产物；``override`` 可显式替换。

        例：``aims_scf`` 的 ``geometry.in`` 默认取 ``au60_build/geometry.in``；
        ``--geometry other/geometry.in`` 时直接用它（阶段独立运行的关键）。
        """
        if override is not None:
            path = self.resolve(override)
            if path is None or not path.exists():
                raise StageBlocked(
                    normalize_stage(stage),
                    "显式输入不存在: %s" % (override,),
                    [str(path)] if path else [])
            return path
        prev = previous_stage(stage)
        if prev is None:
            raise StageBlocked(normalize_stage(stage), "这是第一个阶段，没有上游产物")
        return self.require_stage_file(prev, filename, target=normalize_stage(stage))

    # -- 资源 / 环境 ------------------------------------------------------- #
    def estimate(self, stage: str, natoms: int) -> dict:
        return estimate_smj_resources(normalize_stage(stage), natoms)

    def stage_env(self, stage: str, resources: dict | None = None) -> dict:
        """配置 -> 输入写出器的 ``env``（键名与 vendored inputs.py 的 ENV_SPECS 对齐）。"""
        stage = normalize_stage(stage)
        cluster = self.config.cluster
        res = dict(resources or {})
        env: dict[str, Any] = {"partition": cluster.partition or "normal"}
        # 用户显式 env 先行：后面只用估算值补缺，不覆盖用户值。
        env.update({k: v for k, v in cluster.env.items() if v is not None})

        if stage == "gaussian_opt":
            env.setdefault("ntasks", cluster.ntasks or 1)
            env.setdefault("nproc", cluster.ntasks or 8)
            env.setdefault("g16_mem", cluster.mem or res.get("mem") or "16G")
            env.setdefault("g16_walltime", cluster.time_limit or res.get("time_limit") or "01:00:00")
            if cluster.gres and str(cluster.gres).startswith("gpu"):
                env.setdefault("g16_gres", cluster.gres)
            return env

        if stage in ("aims_tip_relax", "au_relax", "aims_scf"):
            env.setdefault("ntasks", cluster.ntasks or res.get("ntasks") or 4)
            env.setdefault("aims_mem", cluster.mem or res.get("mem") or "16G")
            env.setdefault("aims_walltime", cluster.time_limit or res.get("time_limit") or "02:00:00")
            if cluster.gres and str(cluster.gres).startswith("gpu"):
                env.setdefault("aims_gres_sched", cluster.gres)
            return env

        if stage == "aitranss_te":
            env.setdefault("ntasks", 1)
            env.setdefault("te_cpus", res.get("cpus") or 16)
            env.setdefault("te_mem", cluster.mem or res.get("mem") or "8G")
            env.setdefault("te_walltime", cluster.time_limit or res.get("time_limit") or "02:00:00")
            if cluster.gres and str(cluster.gres).startswith("gpu"):
                env.setdefault("te_gres", cluster.gres)
            return env

        return env

    # -- 远端 -------------------------------------------------------------- #
    def remote_stage_dir(self, stage: str) -> str:
        """远端阶段目录；``aitranss_te`` 与 ``aims_scf`` 共用（Aitranss 读 SCF 产物）。"""
        stage = normalize_stage(stage)
        root = self.config.cluster.remote_root or ""
        if not root:
            raise SmjError("run.mode=slurm-ssh 需要 cluster.remote_root")
        physical = "aims_scf" if stage == "aitranss_te" else stage
        return remote_join(root, self.config.run.name, physical)

    def remote_run(self, command: str, *, timeout: float | None = None,
                   cwd: str | Path | None = None) -> CommandResult:
        if self.transport is None:
            raise SmjError("当前上下文没有传输层（dry-run？）")
        return self.transport.run(command, timeout=timeout, cwd=cwd)

    def remote_upload(self, local: Path, remote_dir: str, name: str | None = None) -> str:
        if self.transport is None:
            raise SmjError("当前上下文没有传输层（dry-run？）")
        remote = remote_join(remote_dir, name or Path(local).name)
        self.transport.put(local, remote).expect_ok("上传 %s" % local)
        return remote

    def remote_download(self, remote_dir: str, name: str, local_dir: Path) -> Path:
        if self.transport is None:
            raise SmjError("当前上下文没有传输层（dry-run？）")
        local = local_dir / name
        self.transport.get(remote_join(remote_dir, name), local).expect_ok("下载 %s" % name)
        return local

    # -- 常用路径 ---------------------------------------------------------- #
    def species_dir(self) -> Path | None:
        raw = self.config.cluster.species_dir or self.config.paths.species_defaults
        if not raw:
            return None
        return self.resolve(raw)

    def templates_dir(self) -> Path | None:
        """用户覆盖的模板目录；None = 用 vendored 包内模板。"""
        raw = self.config.paths.templates
        return self.resolve(raw) if raw else None

    def override_path(self, key: str) -> Path | None:
        """CLI 覆盖路径：按当前工作目录解析（与 yaml 内的配置路径区分开）。"""
        value = self.overrides.get(key)
        if not value:
            return None
        path = Path(str(value)).expanduser()
        return path if path.is_absolute() else path.resolve()

    def describe(self) -> str:  # pragma: no cover - 日志/调试
        return "run=%s mode=%s dir=%s" % (self.config.name, self.config.mode, self.layout.run_dir)
