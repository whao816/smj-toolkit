r"""run 目录布局与 ``state.json``（阶段状态机）。

布局（与 AgenticDFT 引擎默认一致，因此可以直接对既有 run 做 report/analyze）::

    <run_dir>/
        state.json                 # 本模块维护
        manifest.json              # 产物清单（可选，report 生成）
        report.md                  # 人读报告（可选）
        logs/                      # 本包自己的日志
        gaussian_opt/              # 每个阶段一个目录，阶段之间只通过文件耦合
        aims_tip_relax/
        au60_build/
        aims_scf/
        aitranss_te/
        junction_analyze/

状态机::

    pending -> prepared -> submitted -> waiting -> finished -> collected
                     \-> failed -> (prepare 后重试) prepared

- ``dry-run`` 模式只走到 ``prepared``。
- ``local`` / ``slurm-ssh`` 模式下，``--run`` = prepare -> submit -> wait -> collect。
- ``force`` 允许从 collected/failed 重新 prepare。
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Mapping

from .errors import SmjError

STAGE_ORDER: tuple[str, ...] = (
    "gaussian_opt",
    "aims_tip_relax",
    "au60_build",
    "aims_scf",
    "aitranss_te",
    "junction_analyze",
)

#: 远程阶段（需要 Slurm 作业）；au60_build / junction_analyze 是纯本地阶段。
REMOTE_STAGES: tuple[str, ...] = ("gaussian_opt", "aims_tip_relax", "aims_scf", "aitranss_te")
LOCAL_STAGES: tuple[str, ...] = ("au60_build", "junction_analyze")

STATUSES: tuple[str, ...] = (
    "pending", "prepared", "submitted", "waiting", "finished", "failed", "collected", "skipped",
)

_ALIASES = {
    "gaussian": "gaussian_opt",
    "g16": "gaussian_opt",
    "opt": "gaussian_opt",
    "gaussian-opt": "gaussian_opt",
    "tip": "aims_tip_relax",
    "tip_relax": "aims_tip_relax",
    "tip-relax": "aims_tip_relax",
    "relax": "aims_tip_relax",
    "au60": "au60_build",
    "build": "au60_build",
    "au60-build": "au60_build",
    "electrode": "au60_build",
    "scf": "aims_scf",
    "aims": "aims_scf",
    "aims-scf": "aims_scf",
    "te": "aitranss_te",
    "aitranss": "aitranss_te",
    "aitranss-te": "aitranss_te",
    "transmission": "aitranss_te",
    "transport": "aitranss_te",
    "analyze": "junction_analyze",
    "analysis": "junction_analyze",
    "junction-analyze": "junction_analyze",
    "landauer": "junction_analyze",
}


def now_iso() -> str:
    """本地时区的 ISO 时间戳（人读友好，含偏移）。"""
    return datetime.now().astimezone().isoformat(timespec="seconds")


def normalize_stage(name: str) -> str:
    """把 ``Gaussian-Opt`` / ``te`` / ``scf`` 等写法归一为标准阶段名。"""
    if not name:
        raise SmjError("阶段名不能为空")
    key = str(name).strip().lower().replace(" ", "_")
    key = re.sub(r"_+", "_", key)
    if key in STAGE_ORDER:
        return key
    if key in _ALIASES:
        return _ALIASES[key]
    raise SmjError("未知阶段 %r；可选: %s" % (name, ", ".join(STAGE_ORDER)))


def stage_index(name: str) -> int:
    return STAGE_ORDER.index(normalize_stage(name))


def next_stage(name: str) -> str | None:
    idx = stage_index(name)
    return STAGE_ORDER[idx + 1] if idx + 1 < len(STAGE_ORDER) else None


def previous_stage(name: str) -> str | None:
    idx = stage_index(name)
    return STAGE_ORDER[idx - 1] if idx > 0 else None


def stage_range(start: str | None = None, end: str | None = None) -> list[str]:
    """``--from/--to`` -> 阶段列表（含端点，按执行顺序）。"""
    lo = 0 if start is None else stage_index(start)
    hi = len(STAGE_ORDER) - 1 if end is None else stage_index(end)
    if lo > hi:
        raise SmjError("--from %s 在 --to %s 之后，区间为空" % (start, end))
    return list(STAGE_ORDER[lo:hi + 1])


# --------------------------------------------------------------------------- #
# 单阶段状态
# --------------------------------------------------------------------------- #
@dataclass
class StageState:
    status: str = "pending"
    job_id: str | None = None
    remote_dir: str | None = None
    prepared_at: str | None = None
    submitted_at: str | None = None
    finished_at: str | None = None
    collected_at: str | None = None
    elapsed_s: float | None = None
    exit_code: int | None = None
    note: str | None = None
    error: str | None = None
    artifacts: list[str] = field(default_factory=list)
    resources: dict[str, Any] = field(default_factory=dict)
    results: dict[str, Any] = field(default_factory=dict)
    extra: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        data = {
            "status": self.status,
            "job_id": self.job_id,
            "remote_dir": self.remote_dir,
            "prepared_at": self.prepared_at,
            "submitted_at": self.submitted_at,
            "finished_at": self.finished_at,
            "collected_at": self.collected_at,
            "elapsed_s": self.elapsed_s,
            "exit_code": self.exit_code,
            "note": self.note,
            "error": self.error,
            "artifacts": list(self.artifacts),
            "resources": dict(self.resources),
            "results": dict(self.results),
        }
        for key, value in self.extra.items():
            data.setdefault(key, value)
        return data

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "StageState":
        known = {f for f in cls.__dataclass_fields__ if f != "extra"}
        kwargs = {k: data.get(k) for k in known if k in data}
        kwargs["status"] = str(kwargs.get("status") or "pending")
        kwargs["artifacts"] = list(kwargs.get("artifacts") or [])
        kwargs["resources"] = dict(kwargs.get("resources") or {})
        kwargs["results"] = dict(kwargs.get("results") or {})
        state = cls(**kwargs)
        state.extra = {k: v for k, v in data.items() if k not in known}
        return state


# --------------------------------------------------------------------------- #
# run 状态
# --------------------------------------------------------------------------- #
@dataclass
class RunState:
    run: str
    path: Path | None = None
    created_at: str = field(default_factory=now_iso)
    updated_at: str = field(default_factory=now_iso)
    stages: dict[str, StageState] = field(default_factory=dict)
    results: dict[str, Any] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    extra: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        for stage in STAGE_ORDER:
            self.stages.setdefault(stage, StageState())

    # -- 查询 -------------------------------------------------------------- #
    def stage(self, name: str) -> StageState:
        return self.stages[normalize_stage(name)]

    def status(self, name: str) -> str:
        return self.stage(name).status

    def statuses(self) -> dict[str, str]:
        return {name: self.stages[name].status for name in STAGE_ORDER}

    def is_done(self, name: str) -> bool:
        return self.status(name) in ("finished", "collected", "skipped")

    def is_failed(self, name: str) -> bool:
        return self.status(name) == "failed"

    def first_incomplete(self) -> str | None:
        for name in STAGE_ORDER:
            if self.status(name) not in ("finished", "collected", "skipped"):
                return name
        return None

    # -- 修改 -------------------------------------------------------------- #
    def touch(self) -> None:
        self.updated_at = now_iso()

    def set(self, name: str, **fields: Any) -> StageState:
        stage = self.stage(name)
        for key, value in fields.items():
            if key == "status":
                if value not in STATUSES:
                    raise SmjError("非法阶段状态 %r（可选 %s）" % (value, "/".join(STATUSES)))
                stage.status = value
            elif key == "results":
                stage.results.update(value or {})
            elif key == "artifacts":
                stage.artifacts = list(value or [])
            elif key in stage.__dataclass_fields__:
                setattr(stage, key, value)
            else:
                stage.extra[key] = value
        self.touch()
        return stage

    def set_result(self, key: str, value: Any) -> None:
        self.results[key] = value
        self.touch()

    # -- 持久化 ------------------------------------------------------------ #
    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": "smj.state.v1",
            "run": self.run,
            "created_at": self.created_at,
            "updated_at": self.updated_at,
            "stages": {name: self.stages[name].to_dict() for name in STAGE_ORDER},
            "results": dict(self.results),
            "warnings": list(self.warnings),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any], *, run: str = "", path: Path | None = None) -> "RunState":
        stages: dict[str, StageState] = {}
        raw_stages = data.get("stages") or {}
        if isinstance(raw_stages, Mapping):
            for name, payload in raw_stages.items():
                if not name or not isinstance(payload, Mapping):
                    continue
                try:
                    stages[normalize_stage(name)] = StageState.from_dict(payload)
                except SmjError:
                    continue
        known = {"schema", "run", "created_at", "updated_at", "stages", "results", "warnings"}
        state = cls(
            run=str(data.get("run") or run),
            path=path,
            created_at=str(data.get("created_at") or now_iso()),
            updated_at=str(data.get("updated_at") or now_iso()),
            stages=stages,
            results=dict(data.get("results") or {}),
            warnings=list(data.get("warnings") or []),
        )
        state.extra = {k: v for k, v in data.items() if k not in known}
        return state

    @classmethod
    def load(cls, path: str | Path, *, run: str = "") -> "RunState":
        p = Path(path)
        if not p.exists():
            return cls(run=run or p.parent.name, path=p)
        try:
            data = json.loads(p.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise SmjError("state.json 读取失败 (%s): %s" % (p, exc)) from exc
        if not isinstance(data, Mapping):
            raise SmjError("state.json 顶层应为对象: %s" % p)
        return cls.from_dict(data, run=run or p.parent.name, path=p)

    def save(self, path: str | Path | None = None) -> Path:
        target = Path(path) if path else self.path
        if target is None:
            raise SmjError("RunState 未绑定路径，无法保存")
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = json.dumps(self.to_dict(), indent=2, ensure_ascii=False) + "\n"
        tmp = target.with_name(target.name + ".tmp")
        tmp.write_text(payload, encoding="utf-8", newline="\n")
        os.replace(tmp, target)
        self.path = target
        return target


# --------------------------------------------------------------------------- #
# 目录布局
# --------------------------------------------------------------------------- #
@dataclass
class RunLayout:
    """``<run_dir>`` 下各阶段目录与元数据文件的位置。"""

    run_dir: Path

    def __post_init__(self) -> None:
        self.run_dir = Path(self.run_dir).resolve()

    @property
    def state_path(self) -> Path:
        return self.run_dir / "state.json"

    @property
    def manifest_path(self) -> Path:
        return self.run_dir / "manifest.json"

    @property
    def report_path(self) -> Path:
        return self.run_dir / "report.md"

    @property
    def logs_dir(self) -> Path:
        return self.run_dir / "logs"

    def stage_dir(self, name: str) -> Path:
        return self.run_dir / normalize_stage(name)

    def exists(self, name: str) -> bool:
        return self.stage_dir(name).is_dir()

    def artifact(self, stage: str, *parts: str) -> Path:
        return self.stage_dir(stage).joinpath(*parts)

    def outputs(self, stage: str, names: tuple[str, ...] | list[str]) -> list[Path]:
        return [self.artifact(stage, name) for name in names]

    def ensure(self) -> "RunLayout":
        self.run_dir.mkdir(parents=True, exist_ok=True)
        return self

    def ensure_stage(self, name: str) -> Path:
        path = self.stage_dir(name)
        path.mkdir(parents=True, exist_ok=True)
        return path

    # 常见上游产物的快捷方式（便于"阶段独立运行"时发现默认输入）
    def gaussian_xyz(self) -> Path:
        return self.artifact("gaussian_opt", "opt.xyz")

    def tip_geometry(self) -> Path:
        return self.artifact("aims_tip_relax", "geometry.in.next_step")

    def electrode_geometry(self) -> Path:
        return self.artifact("au60_build", "geometry.in")

    def electrode_metadata(self) -> Path:
        return self.artifact("au60_build", "electrode_metadata.json")

    def scf_dir(self) -> Path:
        return self.stage_dir("aims_scf")

    def te_dat(self) -> Path:
        return self.artifact("aitranss_te", "TE.dat")
