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
from typing import TYPE_CHECKING, Any, Callable, Mapping

from .config import RunConfig
from .context import StageContext
from .errors import SmjError, StageError
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
        self.state = RunState.load(self.layout.state_path, run=config.name)
        self.state.path = self.layout.state_path

        self.transport = transport
        self.dry_run = (config.mode == "dry-run") if dry_run is None else bool(dry_run)
        self.overrides = dict(overrides or {})
        self.logger = logger or (lambda message: None)
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
        return attr or default_stage_script(name)

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
                raise StageError(
                    "阶段 %s 正在远端运行（status=%s）；确认要重来请加 --force"
                    % (key, status))
            if status in ("finished", "collected"):
                raise StageError(
                    "阶段 %s 已完成（status=%s）；重跑请加 --force" % (key, status))

        if getattr(self.config.run, "clean", False):
            self._clean_stage_dir(key)
        stage_dir = self.layout.ensure_stage(key)

        # 2026-11-11 本地补丁 F3：重新 prepare 时清掉上一次的采集结果与回取台账。
        # 旧实现只把 status 改回 prepared，results / download_* 原样保留，
        # status/report 于是把上一轮的数值当成"本次已产出"展示。
        self.state.clear_results(key)
        for stale in ("download_fetched", "download_head_fetched", "download_missing"):
            self.state.stage(key).extra.pop(stale, None)

        self._log("[%s] prepare ..." % key)
        prepared = stage.prepare(self.context, stage_dir)
        payload = _prepared_payload(prepared)
        artifacts = payload.get("artifacts") or list(getattr(stage, "outputs", ()) or ())
        local_only = bool(payload.get("local_only"))
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
        )
        self._save()
        for note in payload.get("notes") or ():
            self._log("[%s] %s" % (key, note))
        return payload

    # -- submit ------------------------------------------------------------ #
    def submit(self, name: str, *, force: bool | None = None) -> dict[str, Any]:
        """提交阶段作业（slurm-ssh 用 sbatch；local 同步执行脚本）。"""
        key = normalize_stage(name)
        stage = self.get_stage(key)
        if not getattr(stage, "needs_submit", True):
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

        script = self._stage_script(key)
        if not script:
            raise StageError("阶段 %s 没有作业脚本，无法提交" % key)
        args = list(SMJ_SCRIPT_ARGS.get(script, ()))
        stage_dir = self.layout.ensure_stage(key)
        self._reject_unconfigured(key, script, stage_dir)

        if self.config.mode == "slurm-ssh":
            submission = self._session().submit(key, stage_dir, script=script)
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

        download: dict[str, Any] | None = None
        local_only = bool(self.state.stage(key).extra.get("local_only"))
        if self.config.mode == "slurm-ssh" and not self.dry_run and not local_only:
            self._log("[%s] 回取产物 ..." % key)
            download = self._session().download(key, stage_dir)
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
        self.state.set(
            key, status="collected", collected_at=now_iso(),
            artifacts=list(artifacts), results=dict(payload.get("results") or {}),
            note=None if not download or not download["missing"] else (
                "未取到: " + ", ".join(download["missing"])),
            error=None)
        self._save()
        for note in payload.get("notes") or ():
            self._log("[%s] %s" % (key, note))
        return payload

    # -- run --------------------------------------------------------------- #
    def run_stage(self, name: str, *, force: bool | None = None) -> StageOutcome:
        """单阶段全链路；业务失败返回 ``ok=False``（不抛异常）。"""
        key = normalize_stage(name)
        outcome = StageOutcome(stage=key, action="run")
        try:
            stage = self.get_stage(key)
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

            outcome.collected = self.collect(key)
            outcome.status = "collected"
            outcome.ok = True
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

    def run(self, start: str | None = None, end: str | None = None, *,
            keep_going: bool | None = None,
            force: bool | None = None) -> list[StageOutcome]:
        """按 ``--from/--to`` 链式推进；已 collected 的阶段默认跳过。"""
        names = stage_range(start, end)
        keep = bool(getattr(self.config.run, "keep_going", False)
                    if keep_going is None else keep_going)
        force = self.force if force is None else bool(force)
        outcomes: list[StageOutcome] = []

        for key in names:
            status = self.state.status(key)
            if not force and status == "collected":
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
            if not force and status == "finished":
                # 远端已结束但还没 collect：只补收尾，不重跑。
                try:
                    payload = self.collect(key)
                except SmjError as exc:
                    self.state.set(key, status="failed", error=str(exc))
                    self._save()
                    outcomes.append(StageOutcome(
                        stage=key, status="failed", ok=False, message=str(exc)))
                    if not keep:
                        break
                    continue
                outcomes.append(StageOutcome(
                    stage=key, status="collected", ok=True,
                    message="作业已结束，补 collect 完成", collected=payload))
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
