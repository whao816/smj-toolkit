"""阶段编排：prepare -> submit -> wait -> collect，并维护 ``<run>/state.json``。

设计约束：

* **不 import 阶段实现**。阶段注册表 ``smj.stages`` 的 ``__init__`` 会 import 全部
  六个阶段模块；这里只在真正需要时（函数内）``from .stages.base import get_stage``，
  测试也可以注入 ``stage_getter`` 用假阶段跑通整条链路。
* 只依赖 :class:`Stage` 的两个方法（``prepare`` / ``collect``）与
  ``Prepared.script`` / ``Prepared.artifacts``，其余字段按 duck typing 读取。
* 状态机见 :data:`smj.state.STATUSES`；超时只把状态留在 ``waiting``，
  **不写 failed**（作业可能还在跑）。
"""

from __future__ import annotations

import shutil
import time
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Mapping, Sequence

from .config import RunConfig
from .context import StageContext
from .engines import engine_for_config
from .errors import ActiveJobError, SmjError, StageError, SubmitUncertain
from .provenance import ensure_run_spec, hash_stage_inputs, hash_stage_outputs
from .remote import (
    SMJ_SCRIPT_ARGS,
    RemoteSession,
    read_local_sentinel,
    run_script_local,
    stage_script as default_stage_script,
)
from .state import (
    DONE_STATUSES,
    RunLayout,
    RunState,
    normalize_stage,
    now_iso,
    stage_range,
)
from .transport import Transport, make_transport

if TYPE_CHECKING:  # pragma: no cover - 仅类型检查
    from .stages.base import Collected, Prepared, Stage


# --------------------------------------------------------------------------- #
# 延迟导入的默认阶段获取器
# --------------------------------------------------------------------------- #
def lazy_stage_getter(name: str) -> "Stage":
    """从 ``smj.stages.base`` 取阶段实例；模块未就绪时给出可读错误。"""
    try:
        from .stages.base import get_stage
    except ImportError as exc:  # 阶段模块尚未落地/安装不完整
        raise SmjError(
            "阶段模块不可用（smj.stages 未完整安装）: %s\n"
            "提示: 在仓库开发环境里确认 smj/stages/ 下六个阶段模块存在。"
            % exc) from exc
    try:
        return get_stage(name)
    except KeyError as exc:
        raise SmjError(str(exc).strip("'\"")) from exc


# --------------------------------------------------------------------------- #
# duck-typed 结果归一化
# --------------------------------------------------------------------------- #
def _prepared_payload(prepared: Any) -> dict[str, Any]:
    """把 Prepared / 假阶段返回值归一化成可 JSON 化的 dict。"""
    if prepared is None:
        return {}
    to_dict = getattr(prepared, "to_dict", None)
    if callable(to_dict):
        return dict(to_dict() or {})
    if isinstance(prepared, Mapping):
        return dict(prepared)
    return {
        "stage": getattr(prepared, "stage", "") or "",
        "inputs": [str(p) for p in (getattr(prepared, "inputs", ()) or ())],
        "resources": dict(getattr(prepared, "resources", {}) or {}),
        "script": getattr(prepared, "script", None),
        "artifacts": list(getattr(prepared, "artifacts", ()) or ()),
        "notes": list(getattr(prepared, "notes", ()) or ()),
        "commands": list(getattr(prepared, "commands", ()) or ()),
        "local_only": bool(getattr(prepared, "local_only", False)),
        "parent_stage": getattr(prepared, "parent_stage", None),
        "parent_fingerprint": getattr(prepared, "parent_fingerprint", None),
    }


def _collected_payload(collected: Any) -> dict[str, Any]:
    if collected is None:
        return {}
    to_dict = getattr(collected, "to_dict", None)
    if callable(to_dict):
        return dict(to_dict() or {})
    if isinstance(collected, Mapping):
        return dict(collected)
    output_dir = getattr(collected, "output_dir", None)
    return {
        "stage": getattr(collected, "stage", "") or "",
        "artifacts": list(getattr(collected, "artifacts", ()) or ()),
        "results": dict(getattr(collected, "results", {}) or {}),
        "notes": list(getattr(collected, "notes", ()) or ()),
        "output_dir": str(output_dir) if output_dir else None,
    }


# --------------------------------------------------------------------------- #
# 运行结果
# --------------------------------------------------------------------------- #
@dataclass
class StageOutcome:
    """``run_stage`` / ``run`` 的单阶段结果（CLI 打印与退出码依据）。"""

    stage: str
    action: str = "run"
    status: str = "pending"
    ok: bool = False
    message: str = ""
    prepared: dict[str, Any] | None = None
    submission: dict[str, Any] | None = None
    wait: dict[str, Any] | None = None
    collected: dict[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "stage": self.stage,
            "action": self.action,
            "status": self.status,
            "ok": self.ok,
            "message": self.message,
            "prepared": self.prepared,
            "submission": self.submission,
            "wait": self.wait,
            "collected": self.collected,
        }


# --------------------------------------------------------------------------- #
# 编排器
# --------------------------------------------------------------------------- #
class StageRunner:
    """按阶段推进一次 SMJ run；所有状态都落到 ``<run>/state.json``。"""

    def __init__(self, config: RunConfig, *, run_dir: str | Path | None = None,
                 transport: Transport | None = None,
                 logger: Callable[[str], Any] | None = None,
                 overrides: Mapping[str, Any] | None = None,
                 dry_run: bool | None = None,
                 stage_getter: Callable[[str], Any] | None = None,
                 force: bool = False,
                 poll_s: float | None = None,
                 timeout_s: float | None = None,
                 sleep: Callable[[float], Any] | None = None,
                 clock: Callable[[], float] | None = None):
        self.config = config
        self.run_dir = Path(run_dir or config.run_dir).resolve()
        self.layout = RunLayout(self.run_dir)
        self.layout.ensure()
        self.state = RunState.load(self.layout.state_path, run=config.name, backend=config.transport.backend)
        self.state.path = self.layout.state_path

        self.transport = transport
        self.dry_run = (config.mode == "dry-run") if dry_run is None else bool(dry_run)
        self.overrides = dict(overrides or {})
        self.logger = logger or (lambda message: None)
        if config.transport.backend == "transiesta" and stage_getter is None:
            from .siesta import get_stage
            stage_getter = get_stage
        self.stage_getter = stage_getter or lazy_stage_getter
        self.force = bool(force)
        self.poll_s = poll_s
        self.timeout_s = timeout_s
        self.sleep = sleep
        self.clock = clock

        self._remote_session: RemoteSession | None = None
        self.context = StageContext(
            config=config, layout=self.layout, state=self.state,
            transport=transport, dry_run=self.dry_run,
            overrides=self.overrides, logger=self.logger)

    # -- 基础设施 ---------------------------------------------------------- #
    def _log(self, message: str) -> None:
        try:
            self.logger(message)
        except Exception:  # noqa: BLE001 - 日志永远不能打断运行
            pass

    def _save(self) -> None:
        self.state.touch()
        self.state.save(self.layout.state_path)

    def get_stage(self, name: str) -> "Stage":
        key = normalize_stage(name)
        try:
            return self.stage_getter(key)
        except KeyError as exc:
            raise SmjError(str(exc).strip("'\"")) from exc
        except SmjError:
            raise
        except Exception as exc:  # noqa: BLE001 - 注入的假 getter 可能抛别的
            raise SmjError("取阶段 %r 失败: %s" % (name, exc)) from exc

    def _invalidate_siesta_children(self, name):
        from .state import StageState
        affected = {name}
        for child in self.state.order:
            if affected.intersection(self.get_stage(child).depends):
                affected.add(child)
        for child in affected:
            if self.state.status(child) in ("submitted", "waiting", "submit_unknown"):
                raise ActiveJobError(f"{child}: 作业仍活跃，请先核对调度状态")
        for child in affected:
            self.state.stages[child] = StageState()
        self._save()

    def _ensure_transport(self) -> Transport:
        if self.transport is None:
            self.transport = make_transport(self.config)
        self.context.transport = self.transport
        return self.transport

    def _session(self) -> RemoteSession:
        if self._remote_session is None:
            self._remote_session = RemoteSession(
                self.context, transport=self._ensure_transport())
        return self._remote_session

    def session(self) -> RemoteSession:
        """对外的远端会话（``smj pdos --fetch`` 这类显式回取要用）。"""
        return self._session()

    def _stage_script(self, name: str) -> str | None:
        """优先用 prepare 记下的脚本名（跨进程 --submit 也能工作）。"""
        state = self.state.stage(name)
        script = (state.extra or {}).get("script") if state else None
        if script:
            return str(script)
        try:
            stage = self.get_stage(name)
        except SmjError:
            stage = None
        attr = getattr(stage, "script", None) if stage is not None else None
        # 旧 run 没有 prepared.script：优先用兼容表的回退值（molecule_opt 是
        # sub_g16.sh，历史 gaussian_opt 目录），再退到阶段类的缺省脚本。
        return default_stage_script(name) or attr

    def _reject_unconfigured(self, name: str, script: str, stage_dir: str | Path) -> None:
        """脚本里仍有 ``UNCONFIGURED`` 占位时拒绝提交（2026-11-11 本地补丁 F3）。

        六个脚本模板把集群可执行路径写成 ``${AIMS_BIN:-UNCONFIGURED}`` 这类占位；
        未在 ``cluster.env`` 里配置时，sbatch 照样能提交，作业在远端立刻失败，
        用户要等一轮排队 + 轮询才知道。这里在本地就挡住，并把缺的键名报出来。
        """
        path = Path(stage_dir) / script
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            return  # 脚本还没落盘（例如 --force 前的探测）交给后续步骤报错
        if "UNCONFIGURED" not in text:
            return
        hits = [line.strip() for line in text.splitlines() if "UNCONFIGURED" in line]
        raise StageError(
            "阶段 %s 的作业脚本 %s 仍含未配置的占位 UNCONFIGURED：%s；"
            "请先在 smj.yaml 的 cluster.env 里填好对应路径（如 aims_bin / g16_bin / "
            "aitranss_bin）再提交" % (name, script, "; ".join(hits[:3])))

    def _clean_stage_dir(self, name: str) -> None:
        """``run.clean: true`` 时清空阶段目录（只允许删 run 目录下的子目录）。"""
        stage_dir = self.layout.stage_dir(name).resolve()
        run_dir = self.run_dir.resolve()
        try:
            stage_dir.relative_to(run_dir)
        except ValueError:
            raise StageError(
                "拒绝清空 run 目录之外的路径: %s（run_dir=%s）"
                % (stage_dir, run_dir)) from None
        if stage_dir.exists():
            self._log("[%s] run.clean=true：清空 %s" % (name, stage_dir))
            shutil.rmtree(stage_dir)
        stage_dir.mkdir(parents=True, exist_ok=True)

    # -- prepare ----------------------------------------------------------- #
    def prepare(self, name: str, *, force: bool | None = None) -> dict[str, Any]:
        """生成阶段输入并写 ``state.status = prepared``；dry-run 到此为止。"""
        key = normalize_stage(name)
        stage = self.get_stage(key)
        status = self.state.status(key)
        force = self.force if force is None else bool(force)
        if not force:
            if status in ("submitted", "waiting"):
                # 2026-11-20 复核 P1-5：用 ActiveJobError 区分"作业还在跑"与
                # "阶段失败"。run_stage 的兜底 except 不能再把它写成 failed。
                raise ActiveJobError(
                    "阶段 %s 正在远端运行（status=%s）；确认要重来请加 --force"
                    % (key, status))
            if status == "submit_unknown":
                raise StageError(
                    "阶段 %s 上一次 sbatch 结果未知（可能已有作业在跑）：先按远端目录 "
                    "`squeue`/`sacct` 核对，确认没有作业后再加 --force 重新 prepare"
                    % key)
            if status in ("finished", "collected"):
                raise StageError(
                    "阶段 %s 已完成（status=%s）；重跑请加 --force" % (key, status))

        siesta_record = {}
        if self.state.backend == "transiesta":
            from .siesta.workflow import parents, stage_protocol
            siesta_record = {"siesta_protocol": stage_protocol(self.config, key),
                             "siesta_parents": parents(self.context, stage)}

        # 2026-11 P0-1：run 级配置快照先落盘，再开始产生任何阶段输入。
        # 配置漂移默认在这里拒绝；--force 才重绑（旧快照留 run_spec.previous.json）。
        spec = ensure_run_spec(self.config, self.run_dir, force=force)
        if siesta_record:
            self._invalidate_siesta_children(key)

        if getattr(self.config.run, "clean", False):
            self._clean_stage_dir(key)
        stage_dir = self.layout.ensure_stage(key)

        # 2026-11-11 本地补丁 F3：重新 prepare 时清掉上一次的采集结果与回取台账。
        # 旧实现只把 status 改回 prepared，results / download_* 原样保留，
        # status/report 于是把上一轮的数值当成"本次已产出"展示。
        self.state.clear_results(key)
        for stale in ("download_fetched", "download_head_fetched", "download_missing"):
            self.state.stage(key).extra.pop(stale, None)

        # 0.1.2：slurm-ssh 下 prepare 也可能要访问集群（electrode_build 从远端
        # electrodes.library 取回库团簇）。__init__ 里传输层是惰性创建的，
        # 这里先建好，否则 ctx.transport 为 None，resolve_cluster 的远端分支
        # 永远走不到，默认 source: library 的本地阶段直接失败。
        if self.config.mode == "slurm-ssh" and not self.dry_run:
            self._ensure_transport()

        self._log("[%s] prepare ..." % key)
        prepared = stage.prepare(self.context, stage_dir)
        payload = _prepared_payload(prepared)
        artifacts = payload.get("artifacts") or list(getattr(stage, "outputs", ()) or ())
        local_only = bool(payload.get("local_only"))
        if self.state.backend == "transiesta":
            from .siesta.io import sha
            inputs_sha256 = {str(Path(p).resolve()): sha(p) for p in payload["inputs"]}
        else:
            inputs_sha256 = hash_stage_inputs(payload.get("inputs") or ())
        self.state.set(
            key, status="prepared", prepared_at=now_iso(),
            job_id=None, remote_dir=None, submitted_at=None, finished_at=None,
            collected_at=None, elapsed_s=None, exit_code=None, note=None, error=None,
            artifacts=list(artifacts),
            resources=dict(payload.get("resources") or {}),
            script=None if local_only else (payload.get("script") or getattr(stage, "script", None)),
            inputs=[str(p) for p in payload.get("inputs") or ()],
            notes=list(payload.get("notes") or ()),
            commands=list(payload.get("commands") or ()),
            local_only=local_only,
            inputs_sha256=inputs_sha256,
            outputs_sha256={},
            parent_stage=payload.get("parent_stage"),
            parent_fingerprint=payload.get("parent_fingerprint"),
            protocol_sha256=spec.get("protocol_sha256"), **siesta_record,
        )
        self._save()
        for note in payload.get("notes") or ():
            self._log("[%s] %s" % (key, note))
        return payload

    # -- submit ------------------------------------------------------------ #
    def submit(self, name: str, *, force: bool | None = None,
               sbatch_args: Sequence[str] = ()) -> dict[str, Any]:
        """提交阶段作业（slurm-ssh 用 sbatch；local 同步执行脚本）。

        ``sbatch_args`` 透传给 ``sbatch``，目前只有数组作业的 ``--array=17,203``
        （只重跑失败槽位）用得上；local 模式不支持。
        """
        key = normalize_stage(name)
        stage = self.get_stage(key)
        if self.state.stage(key).extra.get("local_only") or not getattr(stage, "needs_submit", True):
            raise StageError("阶段 %s 是本地阶段，无需提交（直接 --run 或 --collect）" % key)
        status = self.state.status(key)
        force = self.force if force is None else bool(force)
        if self.dry_run:
            raise StageError(
                "dry-run 模式只生成输入，不提交；去掉 --dry-run 或把 run.mode 改成 slurm-ssh")
        # 2026-11-11 本地补丁 F3：显式白名单。此前 failed 只是"没被任何分支拦住"
        # 而意外可提交，语义靠默认分支兜底；现在 prepared / failed 明确可提交
        # （failed 允许直接重提，方便临时错误后重试），其余状态一律要 --force。
        if not force and status not in ("prepared", "failed"):
            if status == "submit_unknown":
                raise StageError(
                    "阶段 %s 上一次 sbatch 结果未知（可能已有作业在跑）：先按远端目录 "
                    "`squeue`/`sacct` 核对，确认没有作业后再加 --force 重提" % key)
            if status in ("submitted", "waiting"):
                raise StageError(
                    "阶段 %s 已提交（status=%s, job_id=%s）；重提交请加 --force"
                    % (key, status, self.state.stage(key).job_id))
            if status in DONE_STATUSES:
                raise StageError("阶段 %s 已完成（status=%s）；重提交请加 --force" % (key, status))
            if status == "pending":
                raise StageError(
                    "阶段 %s 尚未 prepare；先跑 `smj stage %s --prepare` 或直接用 `smj run`"
                    % (key, key))
            raise StageError(
                "阶段 %s 状态为 %r，不能提交；先 --prepare 或加 --force" % (key, status))

        if self.state.backend == "transiesta":
            from .siesta.workflow import check_prepared
            check_prepared(self.context, stage)
        script = self._stage_script(key)
        if not script:
            raise StageError("阶段 %s 没有作业脚本，无法提交" % key)
        args = list(SMJ_SCRIPT_ARGS.get(script, ()))
        stage_dir = self.layout.ensure_stage(key)
        self._reject_unconfigured(key, script, stage_dir)

        if self.config.mode == "slurm-ssh":
            session = self._session()
            try:
                submission = session.submit(key, stage_dir, script=script,
                                            sbatch_args=tuple(sbatch_args))
            except SubmitUncertain as exc:
                # 2026-11-22 复审 P0：sbatch 超时/无 job id 不能当"提交失败"。
                # 作业可能已经在排队，状态必须是显式的“未知”，让下一次 run 停下来
                # 要求人工 squeue 核对；否则 --run 会自动重提，同一目录两个作业。
                self.state.set(
                    key, status="submit_unknown",
                    remote_dir=str(exc.remote_dir or session.stage_dir(key)),
                    script=script, job_id=None,
                    note=str(exc), error=str(exc))
                self._save()
                self._log("[%s] sbatch 结果未知：%s" % (key, str(exc).splitlines()[0]))
                raise
            self.state.set(
                key, status="submitted", submitted_at=now_iso(),
                job_id=submission.job_id, remote_dir=submission.remote_dir,
                script=script, note=None, error=None)
            self._save()
            self._log("[%s] 已提交 job %s -> %s"
                      % (key, submission.job_id, submission.remote_dir))
            return submission.to_dict()

        if self.config.mode != "local":
            raise StageError("run.mode=%r 不支持提交" % self.config.mode)

        # local：同步执行脚本，退出码从 .dsh_done / 进程返回码读取。
        transport = self._ensure_transport()
        timeout = self.timeout_s if self.timeout_s is not None else (
            float(getattr(self.config.cluster, "timeout_hours", 6.0) or 0.0) * 3600.0)
        self._log("[%s] local 执行 %s ..." % (key, script))
        started = (self.clock or time.time)()
        result = run_script_local(transport, stage_dir, script, args, timeout=timeout)
        elapsed = max(0.0, (self.clock or time.time)() - started)
        sentinel = read_local_sentinel(stage_dir)
        exit_code = sentinel if sentinel is not None else result.code
        note = None if exit_code == 0 else (
            (result.err or result.out).strip()[:1000] or "本地脚本返回 rc=%d" % exit_code)
        self.state.set(
            key, status="submitted", submitted_at=now_iso(), job_id=None,
            remote_dir=str(stage_dir), script=script, elapsed_s=elapsed,
            exit_code=exit_code, note=note, error=None)
        self._save()
        self._log("[%s] local 执行结束 rc=%s（%.1fs）" % (key, exit_code, elapsed))
        return {
            "stage": key, "mode": "local", "script": script,
            "remote_dir": str(stage_dir), "job_id": None,
            "exit_code": exit_code, "elapsed_s": elapsed,
            "stdout": result.out[-2000:], "stderr": result.err[-2000:],
        }

    # -- wait -------------------------------------------------------------- #
    def wait(self, name: str, *, poll_s: float | None = None,
             timeout_s: float | None = None,
             sleep: Callable[[float], Any] | None = None,
             clock: Callable[[], float] | None = None) -> dict[str, Any]:
        """等待作业结束；超时保持 ``waiting``（不算 failed）。"""
        key = normalize_stage(name)
        state = self.state.stage(key)
        status = state.status

        # 2026-11-11 本地补丁 F3：已结束/已采集的阶段不能再被 --wait 打回 waiting。
        # 旧实现在 slurm-ssh 分支里先 set(status="waiting") 再轮询，于是对一个
        # collected 阶段跑 `smj stage X --wait` 会把状态机回退（README 推荐的
        # `--prepare --submit --wait` 组合正好能触发这一点）。
        if status in DONE_STATUSES:
            raise StageError(
                "阶段 %s 已结束（status=%s）；--wait 只用于 submitted/waiting 的作业"
                "（已结束未采集请用 --collect，重跑请用 --force）" % (key, status))

        if status == "submit_unknown":
            raise StageError(
                "阶段 %s 上一次 sbatch 结果未知、没有可轮询的 job_id：先按远端目录 "
                "`squeue`/`sacct` 核对；确认没有作业后再 --submit --force" % key)

        if self.config.mode == "slurm-ssh":
            if not state.job_id:
                raise StageError("阶段 %s 没有 job_id；先 --submit 或 --run" % key)
            self.state.set(key, status="waiting")
            self._save()

            def _tick(ticks: int, pending: int) -> None:
                self.state.set(key, note="轮询中：第 %d 次，%d 个作业未结束" % (ticks, pending))
                self._save()
                self._log("[%s] 等待中（第 %d 次轮询，剩 %d）" % (key, ticks, pending))

            result = self._session().wait(
                key, str(state.job_id),
                poll_s=self.poll_s if poll_s is None else poll_s,
                timeout_s=self.timeout_s if timeout_s is None else timeout_s,
                sleep=self.sleep if sleep is None else sleep,
                clock=self.clock if clock is None else clock,
                on_tick=_tick)
        else:
            exit_code = read_local_sentinel(self.layout.stage_dir(key))
            if exit_code is None:
                exit_code = state.exit_code
            if exit_code is None:
                raise StageError(
                    "阶段 %s 没有本地执行记录；先 --submit / --run" % key)
            result = {
                "stage": key, "job_id": None, "finished": True,
                "failed": exit_code != 0, "exit_code": exit_code,
                "note": None if exit_code == 0 else "本地脚本退出码 %d" % exit_code,
                "elapsed_s": state.elapsed_s, "polls": 0,
            }

        if result.get("failed"):
            self.state.set(
                key, status="failed", finished_at=now_iso(),
                exit_code=result.get("exit_code"),
                elapsed_s=result.get("elapsed_s", state.elapsed_s),
                note=result.get("note") or "作业失败",
                error=result.get("note") or "作业失败")
        elif result.get("finished"):
            self.state.set(
                key, status="finished", finished_at=now_iso(),
                exit_code=result.get("exit_code"),
                elapsed_s=result.get("elapsed_s", state.elapsed_s),
                note=result.get("note"), error=None)
        else:
            # 超时：作业可能还在跑，保留 waiting，不判 failed。
            self.state.set(
                key, status="waiting",
                elapsed_s=result.get("elapsed_s", state.elapsed_s),
                note=result.get("note") or "等待超时，作业可能仍在运行（可稍后再 --wait）",
                error=None)
        self._save()
        return result

    # -- collect ----------------------------------------------------------- #
    def collect(self, name: str, *, force: bool | None = None) -> dict[str, Any]:
        """回取产物（远端模式）并解析，写 ``state.status = collected``。"""
        key = normalize_stage(name)
        stage = self.get_stage(key)
        status = self.state.status(key)
        force = self.force if force is None else bool(force)
        if status == "collected" and not force:
            raise StageError("阶段 %s 已 collected；重新解析请加 --force" % key)
        stage_dir = self.layout.ensure_stage(key)

        if self.state.backend == "transiesta":
            from .siesta.workflow import check_prepared
            check_prepared(self.context, stage)
        download: dict[str, Any] | None = None
        local_only = bool(self.state.stage(key).extra.get("local_only"))
        # 0.1.2：本地阶段（electrode_build / junction_analyze）没有远端产物，
        # slurm-ssh 模式下也不能走 download —— 它们不在 artifacts 台账里，
        # 旧实现会抛「未知 smj 阶段」。
        remote_stage = bool(getattr(stage, "needs_submit", True))
        if (self.config.mode == "slurm-ssh" and not self.dry_run
                and not local_only and remote_stage):
            self._log("[%s] 回取产物 ..." % key)
            names = None
            if self.config.transport.backend == "transiesta":
                from .siesta import remote_outputs
                names = remote_outputs(self.context, key)
            elif key == "molecule_opt":
                # 回取清单由引擎声明；None = 沿用 vendored 的 gaussian_opt 计划
                # （mol.gjf/mol.log/sub_g16.sh，Gaussian 口径）。
                remote_files = engine_for_config(self.config).remote_files
                if remote_files is not None:
                    names = list(remote_files)
            download = self._session().download(key, stage_dir, names=names)
            self.state.set(
                key, download_fetched=download["fetched"],
                download_head_fetched=download.get("head_fetched") or [],
                download_missing=download["missing"])
            if download["missing"]:
                # 2026-11-11 本地补丁 F1/F3: 必须产物缺席时**不得**标 collected。
                # 旧实现只记 note 仍标 collected，下游阶段随后在 StageBlocked 处
                # 才炸，且状态机显示"已回取"，误导用户。
                self._save()
                raise StageError(
                    "阶段 %s 回取缺件（远端产物不完整，拒绝标 collected）: %s"
                    % (key, ", ".join(download["missing"])))

        self._log("[%s] collect ..." % key)
        collected = stage.collect(self.context, stage_dir)
        payload = _collected_payload(collected)
        artifacts = payload.get("artifacts") or list(getattr(stage, "outputs", ()) or ())
        unhashed: list[str] = []
        if self.state.backend == "transiesta":
            from .siesta.io import sha
            out = Path(payload.get("output_dir") or stage_dir).resolve()
            outputs_sha256 = {str(p if out == stage_dir.resolve() else out / p): sha(out / p)
                              for p in artifacts}
        else:
            outputs_sha256 = hash_stage_outputs(stage_dir, artifacts, missing=unhashed)
        note = None if not download or not download["missing"] else (
            "未取到: " + ", ".join(download["missing"]))
        if unhashed:
            extra_note = "未记录哈希（缺件或超过大小上限）: " + ", ".join(unhashed)
            note = "%s；%s" % (note, extra_note) if note else extra_note
        self.state.set(
            key, status="collected", collected_at=now_iso(),
            artifacts=list(artifacts), results=dict(payload.get("results") or {}),
            outputs_sha256=outputs_sha256, output_hash_missing=unhashed,
            note=note, error=None)
        self._save()
        for note in payload.get("notes") or ():
            self._log("[%s] %s" % (key, note))
        if key == "electrode_build":
            self._auto_connection_check(stage_dir)
        return payload

    # -- 连接核对自动登记 --------------------------------------------------- #
    def _auto_connection_check(self, stage_dir: Path) -> None:
        """``electrode_build`` 采集后自动登记连接核对（``review.auto_connection_check``）。

        自动登记只做 ``approve_connection_check``：审计/出图不通过时**不**写记录，
        只记一条日志；真正的拦截在 ``smj review`` 的连接核对闸门（fail-closed）。
        """
        cfg = getattr(self.config, "review", None)
        if not bool(getattr(cfg, "auto_connection_check", True)):
            return
        if self.dry_run:
            return
        from .a2 import approve_connection_check  # 本地 import：避免循环依赖

        try:
            result = approve_connection_check(self.run_dir, who="auto", stage=stage_dir)
            if not result.ok and result.code == "figure_failed":
                # 无 matplotlib/无显示：退化为只写记录（审计已在上面通过）。
                result = approve_connection_check(self.run_dir, who="auto",
                                               stage=stage_dir, no_figure=True)
        except SmjError as exc:
            self._log("[electrode_build] 连接核对自动登记失败: %s" % exc)
            return
        if result.ok:
            self._log("[electrode_build] 连接核对: %s（签署人 %s）"
                      % (result.status, result.who or "auto"))
        else:
            self._log("[electrode_build] 连接核对未登记: %s"
                      % (result.reason or result.code))

    # -- run --------------------------------------------------------------- #
    def run_stage(self, name: str, *, force: bool | None = None) -> StageOutcome:
        """单阶段全链路；业务失败返回 ``ok=False``（不抛异常）。"""
        key = normalize_stage(name)
        outcome = StageOutcome(stage=key, action="run")
        try:
            stage = self.get_stage(key)
            force_now = self.force if force is None else bool(force)
            if self.state.status(key) == "submit_unknown" and not force_now:
                # 不能先 prepare：旧状态会被覆盖成 prepared，人工核对的前提就没了。
                outcome.ok = False
                outcome.status = "submit_unknown"
                outcome.message = (
                    "上一次 sbatch 结果未知：作业可能已经在集群排队。先按远端目录 "
                    "`squeue`/`sacct` 核对；确认没有作业后 `smj stage %s --submit "
                    "--force`。" % key)
                return outcome
            outcome.prepared = self.prepare(key, force=force)
            outcome.status = "prepared"
            outcome.ok = True
            if self.dry_run:
                outcome.message = "dry-run：已生成输入，未提交"
                return outcome

            if outcome.prepared.get("local_only") or not getattr(stage, "needs_submit", True):
                outcome.collected = self.collect(key)
                outcome.status = "collected"
                return outcome

            outcome.submission = self.submit(key, force=force)
            outcome.status = "submitted"
            outcome.wait = self.wait(key)
            if outcome.wait.get("failed"):
                outcome.ok = False
                outcome.status = "failed"
                outcome.message = outcome.wait.get("note") or "作业失败"
                return outcome
            if not outcome.wait.get("finished"):
                outcome.ok = False
                outcome.status = "waiting"
                outcome.message = outcome.wait.get("note") or "等待超时，作业可能仍在运行"
                return outcome

            return self._collect_after_job(key, outcome)
        except ActiveJobError as exc:
            # 2026-11-20 复核 P1-5：作业还在远端跑不是失败。保留
            # submitted/waiting 与 job_id，否则下一次 run 的保护失效、会重复提交。
            outcome.ok = False
            outcome.status = self.state.status(key)
            outcome.message = str(exc)
            self._log("[%s] 作业仍在远端运行（status=%s），未重来"
                      % (key, outcome.status))
            return outcome
        except SubmitUncertain as exc:
            # 状态已由 submit() 写成 submit_unknown；不要覆盖、也不要再写 failed。
            outcome.ok = False
            outcome.status = "submit_unknown"
            outcome.message = str(exc)
            self._log("[%s] 提交结果未知，等待人工核对: %s"
                      % (key, outcome.message.splitlines()[0]))
            return outcome
        except SmjError as exc:
            self.state.set(key, status="failed", error=str(exc),
                           note=str(exc).splitlines()[0][:300])
            self._save()
            outcome.ok = False
            outcome.status = self.state.status(key)
            outcome.message = str(exc)
            self._log("[%s] 失败: %s" % (key, outcome.message.splitlines()[0]))
            return outcome
        except (ValueError, OSError) as exc:
            # 2026-11-22 复审 P1：junction/inputs 的配置错误抛 ValueError（如缺
            # species_dir），CLI 顶层兜底只保证退出码，state.json 里看不到原因。
            # 这里统一记 failed + note，保持 "dry-run 失败会标 failed" 的承诺。
            message = "阶段输入/环境错误: %s" % exc
            self.state.set(key, status="failed", error=message,
                           note=message.splitlines()[0][:300])
            self._save()
            outcome.ok = False
            outcome.status = "failed"
            outcome.message = message
            self._log("[%s] 失败: %s" % (key, message.splitlines()[0]))
            return outcome

    def _collect_after_job(self, key: str, outcome: StageOutcome) -> StageOutcome:
        """远端作业已结束：只补 collect，失败时保留 ``finished``。

        2026-11-22 复审 P0：旧实现把 collect 的任意 ``SmjError`` 统一写成
        ``failed``，下一次 ``run`` 于是重新 ``prepare`` + 提交，而远端结果还在、
        作业已跑完——既浪费机时又可能覆盖/删除已完成产物。collect 失败不是作业
        失败：状态留在 ``finished``，``error`` 记原因，下次 run 走 ``_resume_active``
        只重试 collect。
        """
        try:
            outcome.collected = self.collect(key)
        except SmjError as exc:
            message = str(exc)
            self.state.set(
                key, status="finished", error=message, collect_failed=True,
                note="远端作业已结束，但 collect 失败（不重提作业；修复后再次 run "
                     "只重试 collect）: %s" % message.splitlines()[0][:300])
            self._save()
            outcome.ok = False
            outcome.status = "finished"
            outcome.message = message
            self._log("[%s] collect 失败（作业结果保留）: %s"
                      % (key, message.splitlines()[0]))
            return outcome
        outcome.status = "collected"
        outcome.ok = True
        return outcome

    def _resume_active(self, key: str, *, waited: bool) -> StageOutcome:
        """阶段的作业已在远端：只等待/补 collect，不重新提交。

        ``waited=False`` 用于 ``finished``（远端已结束，只差收尾）。
        """
        outcome = StageOutcome(stage=key, action="resume")
        try:
            if waited:
                outcome.wait = self.wait(key)
                if outcome.wait.get("failed"):
                    outcome.status = "failed"
                    outcome.message = outcome.wait.get("note") or "作业失败"
                    return outcome
                if not outcome.wait.get("finished"):
                    outcome.status = "waiting"
                    outcome.message = (outcome.wait.get("note")
                                       or "等待超时，作业可能仍在运行")
                    return outcome
            return self._collect_after_job(key, outcome)
        except SmjError as exc:
            self.state.set(key, status="failed", error=str(exc),
                           note=str(exc).splitlines()[0][:300])
            self._save()
            outcome.status = "failed"
            outcome.message = str(exc)
            return outcome

    def run(self, start: str | None = None, end: str | None = None, *,
            keep_going: bool | None = None,
            force: bool | None = None) -> list[StageOutcome]:
        """按 ``--from/--to`` 链式推进；已 collected 的阶段默认跳过。"""
        names = stage_range(start, end, backend=self.config.transport.backend)
        keep = bool(getattr(self.config.run, "keep_going", False)
                    if keep_going is None else keep_going)
        force = self.force if force is None else bool(force)
        outcomes: list[StageOutcome] = []

        for key in names:
            status = self.state.status(key)
            if not force and status == "collected":
                if self.state.backend == "transiesta":
                    from .siesta.workflow import check_prepared, check_outputs
                    check_prepared(self.context, self.get_stage(key))
                    check_outputs(self.layout, key, self.state.stage(key))
                outcomes.append(StageOutcome(
                    stage=key, status="collected", ok=True,
                    message="已完成，跳过（--force 可重跑）"))
                self._log("[%s] 已完成，跳过" % key)
                continue
            if not force and status == "skipped":
                # 2026-11-11 F3：skipped 在 is_done() 里算完成，这里也必须一致，
                # 否则链式 --run 会静默重跑一个用户显式跳过的阶段。
                outcomes.append(StageOutcome(
                    stage=key, status="skipped", ok=True,
                    message="已标记 skipped，跳过（--force 可重跑）"))
                self._log("[%s] 已跳过，跳过" % key)
                continue
            if not force and status == "submit_unknown":
                message = (
                    "上一次 sbatch 结果未知：作业可能已经在集群排队。先按远端目录 "
                    "`squeue`/`sacct` 核对；确认没有作业后 `smj stage %s --submit "
                    "--force`，不要直接重跑 --run。" % key)
                outcomes.append(StageOutcome(
                    stage=key, status="submit_unknown", ok=False, message=message))
                self._log("[%s] %s" % (key, message))
                break
            if not force and status in ("submitted", "waiting", "finished"):
                # 2026-11-20 复核 P1-5：远端已有作业（或已结束待收尾）时只
                # 等待/补 collect，绝不重新 prepare+submit。旧实现只处理
                # "finished"，"waiting" 会掉进 run_stage → prepare 抛错 →
                # 被标 failed → 下一次 run 保护失效、重复提交。
                out = self._resume_active(key, waited=status != "finished")
                outcomes.append(out)
                if not out.ok:
                    if out.status == "waiting":
                        break
                    if not keep:
                        break
                continue

            outcome = self.run_stage(key, force=force)
            outcomes.append(outcome)
            if not outcome.ok:
                if outcome.status == "waiting":
                    # 上游没结束，链式推进没有意义。
                    break
                if not keep:
                    break
        return outcomes

    # -- 展示辅助 ---------------------------------------------------------- #
    def status_table(self) -> list[dict[str, Any]]:
        """给 CLI status/plan 用的扁平表。"""
        rows = []
        for key in self.state.stages:
            state = self.state.stage(key)
            rows.append({
                "stage": key,
                "status": state.status,
                "job_id": state.job_id,
                "remote_dir": state.remote_dir,
                "elapsed_s": state.elapsed_s,
                "note": state.note,
                "artifacts": list(state.artifacts),
                "error": state.error,
            })
        return rows


__all__ = [
    "StageOutcome",
    "StageRunner",
    "lazy_stage_getter",
]
