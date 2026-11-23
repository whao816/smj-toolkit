r"""run 目录布局与 ``state.json``（阶段状态机）。

布局（与 AgenticDFT 引擎默认一致，因此可以直接对既有 run 做 report/analyze）::

    <run_dir>/
        state.json                 # 本模块维护
        manifest.json              # 产物清单（可选，report 生成）
        report.md                  # 人读报告（可选）
        logs/                      # 本包自己的日志
        molecule_opt/              # 每个阶段一个目录，阶段之间只通过文件耦合
                                   # （旧 run 仍可能是 gaussian_opt/，RunLayout 回退读取）
        aims_tip_relax/
        electrode_build/           # 旧 run 仍可能是 au60_build/ 或 auN_build/
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
import tempfile
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Mapping

from .errors import SmjError

STAGE_ORDER: tuple[str, ...] = (
    "molecule_opt",
    "aims_tip_relax",
    "electrode_build",
    "aims_scf",
    "aitranss_te",
    "junction_analyze",
)

#: 2026-11-23：热导（声子）链。与六阶段电子链**并行独立**，只在显式
#: ``smj stage thermal_xxx`` / ``--run`` 时使用；按设计 DECISIONS D12，
#: 它不出现在 ``smj status`` / ``smj plan`` 的默认阶段序列里，因此不进
#: :data:`STAGE_ORDER`，:func:`stage_index` 对它显式报错。
THERMAL_STAGE_ORDER: tuple[str, ...] = (
    "thermal_device",
    "thermal_relax",
    "thermal_fc",
    "thermal_analyze",
)

#: 2026-11-13：``au_relax`` 已从主链路移除。保留该名字仅用于读取/操作旧 run
#: （``smj stage au_relax`` 仍可用），不再参与 ``run``、``plan``、report 的阶段序列。
LEGACY_STAGES: tuple[str, ...] = ("au_relax",)

#: 远程阶段（需要 Slurm 作业）；electrode_build / junction_analyze 是纯本地阶段。
REMOTE_STAGES: tuple[str, ...] = (
    "molecule_opt", "aims_tip_relax", "aims_scf", "aitranss_te",
)
LOCAL_STAGES: tuple[str, ...] = ("electrode_build", "junction_analyze")

STATUSES: tuple[str, ...] = (
    "pending", "prepared", "submitted", "submit_unknown", "waiting",
    "finished", "failed", "collected", "skipped",
)

#: "该阶段无需再做" 的状态集合。2026-11-11 本地补丁 F3：状态机里 is_done()、
#: ``--run`` 的跳过判定、``--wait`` 的拒绝条件都引用这一个常量，避免三处各自
#: 维护一份清单（此前 ``skipped`` 在 is_done 里算完成，``--run`` 却会重跑它）。
DONE_STATUSES: tuple[str, ...] = ("finished", "collected", "skipped")

_ALIASES = {
    # 2026-11-13：阶段规范名改为 molecule_opt；gaussian_opt 仅作旧 run/旧命令别名。
    "molecule": "molecule_opt",
    "mol": "molecule_opt",
    "molecule-opt": "molecule_opt",
    "molecule_optimization": "molecule_opt",
    "gaussian_opt": "molecule_opt",
    "gaussian-opt": "molecule_opt",
    "gaussian": "molecule_opt",
    "g16": "molecule_opt",
    "opt": "molecule_opt",
    "tip": "aims_tip_relax",
    "tip_relax": "aims_tip_relax",
    "tip-relax": "aims_tip_relax",
    "relax": "aims_tip_relax",
    "au60": "electrode_build",
    "au60_build": "electrode_build",
    "au60-build": "electrode_build",
    "build": "electrode_build",
    "electrode": "electrode_build",
    "electrode-build": "electrode_build",
    "junction_build": "electrode_build",
    "junction-build": "electrode_build",
    "au-relax": "au_relax",
    "aurelax": "au_relax",
    "relax_au": "au_relax",
    "backplane_relax": "au_relax",
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
    # 热导链（2026-11-23）
    "thermal-device": "thermal_device",
    "tdevice": "thermal_device",
    "t_device": "thermal_device",
    "thermal-relax": "thermal_relax",
    "trelax": "thermal_relax",
    "t_relax": "thermal_relax",
    "thermal-fc": "thermal_fc",
    "tfc": "thermal_fc",
    "t_fc": "thermal_fc",
    "finite_fc": "thermal_fc",
    "thermal-fc-array": "thermal_fc",
    "thermal-analyze": "thermal_analyze",
    "tanalyze": "thermal_analyze",
    "t_analyze": "thermal_analyze",
    "thermal-analysis": "thermal_analyze",
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
    if key in STAGE_ORDER or key in LEGACY_STAGES or key in THERMAL_STAGE_ORDER:
        return key
    if key in _ALIASES:
        return _ALIASES[key]
    # 旧 run 的目录/状态键可能叫 au20_build/au60_build/auN_build；统一到
    # electrode_build（具体目录由 RunLayout.stage_dir 回退定位）。
    if re.match(r"^au\d+_build$", key):
        return "electrode_build"
    raise SmjError("未知阶段 %r；可选: %s（legacy: %s；热导链: %s）"
                   % (name, ", ".join(STAGE_ORDER), ", ".join(LEGACY_STAGES),
                      ", ".join(THERMAL_STAGE_ORDER)))


def stage_index(name: str) -> int:
    normalized = normalize_stage(name)
    if normalized in THERMAL_STAGE_ORDER:
        raise SmjError(
            "%s 属于热导链，不在电子六阶段链上；请直接 "
            "`smj stage %s --run`（它不参与 --from/--to 区间）"
            % (normalized, normalized))
    if normalized not in STAGE_ORDER:
        raise SmjError("阶段 %s 已从主链路移除（仅兼容旧 run）" % normalized)
    return STAGE_ORDER.index(normalized)


def next_stage(name: str) -> str | None:
    normalized = normalize_stage(name)
    if normalized == "au_relax":
        return "aims_scf"
    idx = stage_index(normalized)
    return STAGE_ORDER[idx + 1] if idx + 1 < len(STAGE_ORDER) else None


def previous_stage(name: str) -> str | None:
    normalized = normalize_stage(name)
    # au_relax 已不在链路上；旧 run 单独操作它时，其上游仍是 electrode_build。
    if normalized == "au_relax":
        return "electrode_build"
    idx = stage_index(normalized)
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
    #: 2026-11 P0-3：prepare 时记录的输入文件哈希（名字 -> sha256）。
    inputs_sha256: dict[str, str] = field(default_factory=dict)
    #: collect 时记录的产物哈希（相对 stage dir 的名字 -> sha256）。
    outputs_sha256: dict[str, str] = field(default_factory=dict)
    #: 父阶段名与父几何指纹（"这一阶段的输入来自哪份上游文件"）。
    parent_stage: str | None = None
    parent_fingerprint: str | None = None
    #: prepare 时绑定的 run_spec 协议指纹（配置漂移后与 run_spec 对比）。
    protocol_sha256: str | None = None
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
            "inputs_sha256": dict(self.inputs_sha256),
            "outputs_sha256": dict(self.outputs_sha256),
            "parent_stage": self.parent_stage,
            "parent_fingerprint": self.parent_fingerprint,
            "protocol_sha256": self.protocol_sha256,
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
        kwargs["inputs_sha256"] = dict(kwargs.get("inputs_sha256") or {})
        kwargs["outputs_sha256"] = dict(kwargs.get("outputs_sha256") or {})
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
        key = normalize_stage(name)
        if key not in self.stages:
            # 2026-11-13：au_relax 已从主链路移除，但旧 run 的 state 和显式的
            # `smj stage au_relax ...` 仍需可读；按需创建 legacy 段，不污染新 run。
            if key in LEGACY_STAGES:
                self.stages[key] = StageState()
            else:
                raise KeyError(key)
        return self.stages[key]

    def status(self, name: str) -> str:
        return self.stage(name).status

    def statuses(self) -> dict[str, str]:
        return {name: self.stages[name].status for name in STAGE_ORDER}

    def is_done(self, name: str) -> bool:
        return self.status(name) in DONE_STATUSES

    def is_failed(self, name: str) -> bool:
        return self.status(name) == "failed"

    def first_incomplete(self) -> str | None:
        for name in STAGE_ORDER:
            if self.status(name) not in DONE_STATUSES:
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

    def clear_results(self, name: str) -> None:
        """清空某阶段上一次解析出来的 results（2026-11-11 F3：重 prepare 时用）。

        ``set(results=...)`` 是**合并**语义，无法把旧键去掉；重新 prepare 后
        旧 results 仍然躺在 state.json 里，status/report 会显示已失效的数值。
        """
        self.stage(name).results.clear()
        self.touch()

    # -- 持久化 ------------------------------------------------------------ #
    def to_dict(self) -> dict[str, Any]:
        ordered = list(STAGE_ORDER)
        ordered += [name for name in LEGACY_STAGES if name in self.stages]
        return {
            "schema": "smj.state.v1",
            "run": self.run,
            "created_at": self.created_at,
            "updated_at": self.updated_at,
            "stages": {name: self.stages[name].to_dict() for name in ordered},
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
        # 2026-11-22 复审 P2：固定 "<name>.tmp" 在两个进程同时 save 时会互相
        # 截断（后写者把自己的半成品替换进去）。用同目录唯一临时文件 + os.replace，
        # 至少保证每次落盘是完整 JSON；并发仍是"最后写者赢"，但不丢文件结构。
        fd, tmp_name = tempfile.mkstemp(
            prefix=target.name + ".", suffix=".tmp", dir=str(target.parent))
        try:
            with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as handle:
                handle.write(payload)
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(tmp_name, target)
        except BaseException:
            Path(tmp_name).unlink(missing_ok=True)
            raise
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
    def run_spec_path(self) -> Path:
        return self.run_dir / "run_spec.json"

    @property
    def run_spec_previous_path(self) -> Path:
        return self.run_dir / "run_spec.previous.json"

    @property
    def logs_dir(self) -> Path:
        return self.run_dir / "logs"

    def stage_dir(self, name: str) -> Path:
        canonical = normalize_stage(name)
        path = self.run_dir / canonical
        # 2026-11-13：molecule_opt 是规范名；旧 run 的目录仍叫 gaussian_opt。
        # 只在规范目录不存在时回退，避免新 run 误写旧目录。
        if canonical == "molecule_opt":
            legacy = self.run_dir / "gaussian_opt"
            if not path.exists() and legacy.is_dir():
                return legacy
        elif canonical == "electrode_build" and not path.exists():
            # 旧 run：au60_build 是默认名；更早的自定义 run 可能是 auN_build。
            legacy = self.run_dir / "au60_build"
            if legacy.is_dir():
                return legacy
            variants = sorted(
                (child for child in self.run_dir.glob("au*_build") if child.is_dir()),
                key=lambda child: (int(re.match(r"^au(\d+)_build$", child.name).group(1))
                                   if re.match(r"^au(\d+)_build$", child.name) else 10 ** 9,
                                   child.name))
            if variants:
                return variants[0]
        return path

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
    def molecule_xyz(self) -> Path:
        return self.artifact("molecule_opt", "opt.xyz")

    # 旧名字保留（第三方脚本兼容）；molecule_opt 是规范名。
    def gaussian_xyz(self) -> Path:
        return self.molecule_xyz()

    def tip_geometry(self) -> Path:
        return self.artifact("aims_tip_relax", "geometry.in.next_step")

    def electrode_geometry(self) -> Path:
        return self.artifact("electrode_build", "geometry.in")

    def electrode_metadata(self) -> Path:
        return self.artifact("electrode_build", "electrode_metadata.json")

    def scf_dir(self) -> Path:
        return self.stage_dir("aims_scf")

    def te_dat(self) -> Path:
        return self.artifact("aitranss_te", "TE.dat")
