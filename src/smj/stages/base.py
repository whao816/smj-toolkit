"""阶段基类与注册表。

一个阶段只做两件事：

- :meth:`Stage.prepare` —— 把上游产物 + 配置翻译成**本地输入文件**（gjf / geometry.in /
  control.in / tcontrol / 作业脚本）。不提交、不等待、不解析产物。
- :meth:`Stage.collect` —— 解析该阶段目录里的产物，返回结构化结果。不猜测、缺文件就报错。

阶段之间**只通过文件耦合**（``<run>/<stage>/``），所以任意阶段都可以：

* 单独 ``prepare`` 出来人工检查（``--prepare --dry-run``）；
* 拿到别人给的输入文件后单独运行（``--xyz`` / ``--geometry`` / ``--scf-dir`` / ``--te``）；
* 单独 ``collect`` 已有产物（``--collect``），不重跑计算。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, ClassVar

from ..errors import StageBlocked


@dataclass
class Prepared:
    """``prepare()`` 的返回值：写了什么、要提交什么、预期产出什么。"""

    stage: str = ""
    inputs: list[Path] = field(default_factory=list)   # 生成的本地输入文件（含作业脚本）
    resources: dict[str, Any] = field(default_factory=dict)
    script: str | None = None                          # 远端作业脚本名（sub_*.sh）
    artifacts: list[str] = field(default_factory=list)  # 该阶段期望的产物（相对 stage dir）
    notes: list[str] = field(default_factory=list)      # 人读提示（CLI 打印）
    commands: list[str] = field(default_factory=list)   # dry-run 时建议用户手工执行的命令
    #: 2026-11-11 本地补丁 F-LIB：true = 该阶段被配置关闭（如 electrode.relax=false），
    #: prepare 只做本地直通，collect 直接解析本地文件，不提交、不回取。
    local_only: bool = False

    def to_dict(self) -> dict[str, Any]:
        return {
            "stage": self.stage,
            "inputs": [str(p) for p in self.inputs],
            "resources": dict(self.resources),
            "script": self.script,
            "artifacts": list(self.artifacts),
            "notes": list(self.notes),
            "commands": list(self.commands),
            "local_only": self.local_only,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Prepared":
        return cls(
            stage=str(data.get("stage") or ""),
            inputs=[Path(p) for p in data.get("inputs") or []],
            resources=dict(data.get("resources") or {}),
            script=data.get("script"),
            artifacts=list(data.get("artifacts") or []),
            notes=list(data.get("notes") or []),
            commands=list(data.get("commands") or []),
            local_only=bool(data.get("local_only")),
        )


@dataclass
class Collected:
    """``collect()`` 的返回值。"""

    stage: str = ""
    artifacts: list[str] = field(default_factory=list)  # 实际解析到的产物（相对 stage dir）
    results: dict[str, Any] = field(default_factory=dict)  # 结构化结果（写入 state/results.json）
    notes: list[str] = field(default_factory=list)
    output_dir: Path | None = None                       # 产物写入的目录（本地阶段）

    def to_dict(self) -> dict[str, Any]:
        return {
            "stage": self.stage,
            "artifacts": list(self.artifacts),
            "results": dict(self.results),
            "notes": list(self.notes),
            "output_dir": str(self.output_dir) if self.output_dir else None,
        }


class Stage:
    """阶段基类。子类只需声明元数据并实现 prepare/collect。"""

    #: 标准阶段名（见 :data:`smj.state.STAGE_ORDER`）
    name: ClassVar[str] = ""
    #: 直接依赖的阶段（用于 plan/状态检查；通常就是上一阶段）
    depends: ClassVar[tuple[str, ...]] = ()
    #: 是否需要提交远端 Slurm 作业
    needs_submit: ClassVar[bool] = True
    #: 提交用的作业脚本名（``None`` = 本地阶段）
    script: ClassVar[str | None] = None
    #: 人读输入提示（``smj plan`` 打印）
    input_hint: ClassVar[str] = ""
    #: 该阶段的标准产物（相对 stage dir）
    outputs: ClassVar[tuple[str, ...]] = ()

    # -- 子类实现 ---------------------------------------------------------- #
    def prepare(self, ctx, stage_dir: Path) -> Prepared:  # pragma: no cover - 抽象
        raise NotImplementedError

    def collect(self, ctx, stage_dir: Path) -> Collected:  # pragma: no cover - 抽象
        raise NotImplementedError

    # -- 公共工具 ---------------------------------------------------------- #
    def prepared(self, ctx, stage_dir: Path, **kwargs: Any) -> Prepared:
        """带上默认值的 :class:`Prepared` 构造器。"""
        kwargs.setdefault("stage", self.name)
        kwargs.setdefault("artifacts", list(self.outputs))
        kwargs.setdefault("script", self.script)
        return Prepared(**kwargs)

    def collected(self, ctx, stage_dir: Path, **kwargs: Any) -> Collected:
        kwargs.setdefault("stage", self.name)
        kwargs.setdefault("output_dir", Path(stage_dir))
        return Collected(**kwargs)

    # -- 规模估算辅助 ------------------------------------------------------ #
    @staticmethod
    def count_xyz_atoms(path: str | Path) -> int:
        """XYZ 原子数（跳过前两行注释/计数行）。文件缺失时抛 :class:`StageBlocked`。"""
        p = Path(path)
        if not p.exists():
            raise StageBlocked("", "结构文件不存在: %s" % p, [str(p)])
        lines = [ln for ln in p.read_text(encoding="utf-8", errors="replace").splitlines() if ln.strip()]
        if len(lines) < 2:
            raise StageBlocked("", "XYZ 文件为空或格式不完整: %s" % p, [str(p)])
        try:
            declared = int(lines[0].split()[0])
        except (ValueError, IndexError):
            raise StageBlocked("", "XYZ 第一行不是原子数: %s" % p, [str(p)]) from None
        atoms = 0
        for line in lines[2:]:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    float(parts[1]); float(parts[2]); float(parts[3])
                except ValueError:
                    continue
                atoms += 1
        return atoms or declared

    @staticmethod
    def count_aims_atoms(path: str | Path) -> int:
        """FHI-aims ``geometry.in`` 的原子数。"""
        p = Path(path)
        if not p.exists():
            raise StageBlocked("", "geometry.in 不存在: %s" % p, [str(p)])
        count = 0
        for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
            parts = line.split()
            if len(parts) >= 5 and parts[0] == "atom":
                count += 1
        return count


# --------------------------------------------------------------------------- #
# 注册表
# --------------------------------------------------------------------------- #
_REGISTRY: dict[str, Stage] = {}


def register(stage_cls: type[Stage]) -> type[Stage]:
    """类装饰器：把阶段实例登记到全局表（按标准名）。"""
    if not stage_cls.name:
        raise ValueError("阶段类 %s 未声明 name" % stage_cls.__name__)
    if stage_cls.name in _REGISTRY:
        raise ValueError("阶段名重复注册: %s" % stage_cls.name)
    _REGISTRY[stage_cls.name] = stage_cls()
    return stage_cls


def get_stage(name: str) -> Stage:
    from ..state import normalize_stage
    key = normalize_stage(name)
    if key not in _REGISTRY:
        raise KeyError("阶段 %r 未注册（已注册: %s）" % (name, ", ".join(sorted(_REGISTRY))))
    return _REGISTRY[key]


def all_stages() -> list[Stage]:
    """按注册顺序返回全部阶段实例（导入 ``smj.stages`` 后即为 6 阶段链）。"""
    return list(_REGISTRY.values())
