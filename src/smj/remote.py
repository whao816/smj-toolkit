"""SMJ 远端执行层：Slurm 提交、残留清理、哨兵轮询、产物回取。

这里复刻主仓 ``agenticdft/agenticdft/tools/submit.py`` 里单分子结（SMJ）的提交
语义，但只依赖 :class:`smj.transport.Transport` 四个原语（run/put/get/mkdirs），
因此可以用 :class:`smj.testing.FakeTransport` 在离线环境里把整条链路跑一遍。

一次远端阶段执行固定四步，全部落在 ``<remote_root>/<run>/<stage>``：

1. ``mkdir -p`` 远端目录；
2. 清理上一轮残留：``.dsh_done`` / ``slurm-*.out`` / 本阶段旧产物。
   陈旧的 ``.dsh_done`` 会让轮询把上一轮当成这一轮完成（取到半截日志/中间几何），
   所以**清理失败直接拒绝提交**；
3. 上传输入：``aitranss_te`` 只传 ``tcontrol`` / ``sub_te.sh`` /
   ``self.energy.in`` 三个文件，并在 ``aims_scf`` 的远端目录里运行；
   其余阶段上传 stage dir 下所有顶层普通文件，但跳过旧回取产物
   （``mol.log`` / ``aims.out`` / ``TE.dat`` / ``te.out`` / ``.dsh_done``）；
4. ``cd <dir> && sbatch <script> [args...]``，从输出解析 ``Submitted batch job <id>``。

轮询交给 vendored 的 :func:`smj.junction.jobwait.wait_many`，本模块只提供
``Transport -> (code, out, err)`` 的适配器。判定语义（与主仓一致）：

* ``squeue`` 查询失败（非零退出/连不上）**不等于**作业缺席；
* 连续两次成功查询都看不到 job id、且没有 ``.dsh_done`` → ``failed=True``；
* 超时只表示「还没结束」：``finished=False`` / ``failed=False``，**不算失败**；
* 作业退出码从 ``.dsh_done`` 首字段读取（Slurm 的 shell 退出码不可靠）。
"""

from __future__ import annotations

import re
import shlex
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

from .junction import artifacts as _artifacts
from .junction import jobwait as _jobwait
from .errors import SmjError, StageError, SubmitUncertain, TransportError
from .state import normalize_stage
from .transport import CommandResult, Transport, remote_join

#: 阶段 -> 远端作业脚本（与主仓 SMJ_STAGE_SCRIPTS 一致）
SMJ_STAGE_SCRIPTS: dict[str, str] = {
    # 2026-11-13：规范名 molecule_opt；远端脚本由 prepared.script 决定。
    # 这里的值是**旧 run 的回退**（历史 gaussian_opt 目录没有 prepared.script）；
    # 新 run 一律按引擎写进 state 的脚本名（aims → sub_aims.sh）。
    "molecule_opt": "sub_g16.sh",
    "gaussian_opt": "sub_g16.sh",
    "aims_tip_relax": "sub_aims.sh",
    "au_relax": "sub_aims.sh",
    "aims_scf": "sub_aims.sh",
    "aitranss_te": "sub_te.sh",
}

#: 全部作业脚本（detect_stage_script 用）
SMJ_SCRIPTS: tuple[str, ...] = ("sub_g16.sh", "sub_aims.sh", "sub_te.sh")

#: 作业脚本 -> 额外命令行参数
SMJ_SCRIPT_ARGS: dict[str, tuple[str, ...]] = {
    "sub_g16.sh": ("mol.gjf",),
}

#: aitranss_te 只上传这三个输入（复用 aims_scf 的远端目录）
SMJ_AITRANSS_UPLOADS: tuple[str, ...] = ("tcontrol", "sub_te.sh", "self.energy.in")

#: 每个脚本在远端目录里的旧产物（重提交前删除，避免拿上一轮结果）
SMJ_STALE_UPLOADS: dict[str, tuple[str, ...]] = {
    "sub_g16.sh": ("mol.log",),
    "sub_aims.sh": ("aims.out", "geometry.in.next_step"),
    "sub_te.sh": ("TE.dat", "te.out"),
}

#: 完成哨兵（jobwait 与本模块共用）
DONE_SENTINEL: str = _jobwait.SENTINEL

#: 远端专有的旧输出（不属于上传文件，清理时要带上）
_REMOTE_STALE_GLOBS: tuple[str, ...] = ("slurm-*.out",)

_SBATCH_RE = re.compile(r"Submitted\s+batch\s+job\s+(\d+)", re.IGNORECASE)


# --------------------------------------------------------------------------- #
# 纯函数：脚本映射、清理清单、输出解析
# --------------------------------------------------------------------------- #
def stage_script(stage: str) -> str | None:
    """阶段名 -> 作业脚本名；本地阶段返回 ``None``。"""
    try:
        key = normalize_stage(stage)
    except SmjError:
        return None
    return SMJ_STAGE_SCRIPTS.get(key)


def stale_names(script: str | None) -> tuple[str, ...]:
    """重提交前要 ``rm -f`` 的远端残留（含 ``slurm-*.out`` 通配符）。"""
    names = (DONE_SENTINEL, *_REMOTE_STALE_GLOBS)
    if script:
        names += SMJ_STALE_UPLOADS.get(script, ())
    return names


def skip_uploads(script: str | None) -> frozenset[str]:
    """上传时要跳过的本地文件（旧回取产物；``slurm-*.out`` 是远端专有，不在本地）。"""
    names = {DONE_SENTINEL}
    if script:
        names.update(n for n in SMJ_STALE_UPLOADS.get(script, ()) if "*" not in n)
    return frozenset(names)


def detect_stage_script(run_dir: str | Path) -> str | None:
    """在一个 stage 目录里探测作业脚本（依次尝试 sub_g16.sh / sub_aims.sh / sub_te.sh）。"""
    directory = Path(run_dir)
    for script in SMJ_SCRIPTS:
        if (directory / script).is_file():
            return script
    return None


def parse_job_id(text: str | None) -> str | None:
    """从 sbatch 输出里解析 job id；解析不到返回 ``None``。"""
    if not text:
        return None
    match = _SBATCH_RE.search(text)
    return match.group(1) if match else None


def parse_sentinel(text: str | None) -> int | None:
    """解析 ``.dsh_done`` 首字段（退出码）；空/损坏返回 ``None``。"""
    parts = (text or "").split()
    if not parts:
        return None
    try:
        return int(parts[0])
    except (TypeError, ValueError):
        return None


def read_local_sentinel(stage_dir: str | Path) -> int | None:
    """读本地 ``<stage_dir>/.dsh_done`` 的退出码（local 模式收尾用）。"""
    path = Path(stage_dir) / DONE_SENTINEL
    if not path.is_file():
        return None
    return parse_sentinel(path.read_text(encoding="utf-8", errors="replace"))


# --------------------------------------------------------------------------- #
# 数据结构
# --------------------------------------------------------------------------- #
@dataclass
class Submission:
    """一次成功提交的完整记录（runner 写进 state.json）。"""

    stage: str
    mode: str
    script: str
    remote_dir: str
    job_id: str | None = None
    sbatch: str = ""
    uploaded: list[str] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)
    cleaned: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "stage": self.stage,
            "mode": self.mode,
            "script": self.script,
            "remote_dir": self.remote_dir,
            "job_id": self.job_id,
            "sbatch": self.sbatch,
            "uploaded": list(self.uploaded),
            "skipped": list(self.skipped),
            "cleaned": list(self.cleaned),
        }


class _ProbeAdapter:
    """把 :class:`Transport` 适配成 jobwait 期望的 ``rc.run(cmd) -> (code, out, err)``。

    ``jobwait.wait_many`` 内部按 ``code, out, _err = rc.run(cmd, timeout=120)`` 解包，
    直接传 Transport 会因为返回 :class:`CommandResult` 而报错。
    """

    def __init__(self, transport: Transport):
        self.transport = transport

    def run(self, command: str, timeout: float | None = None) -> tuple[int, str, str]:
        result = self.transport.run(command, timeout=timeout)
        return result.code, result.out, result.err


# --------------------------------------------------------------------------- #
# 远端会话
# --------------------------------------------------------------------------- #
class RemoteSession:
    """绑定一个 :class:`smj.context.StageContext` 的远端执行会话。"""

    def __init__(self, ctx, *, transport: Transport | None = None):
        self.ctx = ctx
        self.transport = transport if transport is not None else getattr(ctx, "transport", None)
        if self.transport is None:
            raise TransportError(
                "远端会话缺少 Transport：请传入 transport，或先设置 ctx.transport")

    # -- 路径 -------------------------------------------------------------- #
    def stage_dir(self, stage: str) -> str:
        """阶段的远端目录（``aitranss_te`` 复用 ``aims_scf`` 的目录）。"""
        return self.ctx.remote_stage_dir(stage)

    def remote_dir(self, stage: str) -> str:
        return self.stage_dir(stage)

    # -- 原语 -------------------------------------------------------------- #
    def _run(self, command: str, *, timeout: float = 60) -> CommandResult:
        return self.transport.run(command, timeout=timeout)

    def ensure_dir(self, remote_dir: str) -> str:
        result = self.transport.mkdirs(remote_dir)
        if not result.ok:
            raise TransportError(
                "创建远端目录失败: %s (rc=%d) %s"
                % (remote_dir, result.code, (result.err or result.out).strip()[:500]))
        return remote_dir

    # -- 清理 / 上传 ------------------------------------------------------- #
    def clean_stale(self, stage: str, script: str | None) -> list[str]:
        """删除上一轮残留；失败即拒绝提交（陈旧哨兵会污染完成判定）。"""
        remote_dir = self.stage_dir(stage)
        names = stale_names(script)
        command = "cd %s && rm -f %s" % (
            shlex.quote(remote_dir), " ".join(names))
        result = self._run(command)
        if not result.ok:
            raise TransportError(
                "清理远端残留失败 (%s, rc=%d): %s\n"
                "拒绝继续提交：陈旧的 %s 会让轮询把上一轮当成这一轮完成。"
                % (remote_dir, result.code,
                   (result.err or result.out).strip()[:800], DONE_SENTINEL))
        return list(names)

    def upload_inputs(self, stage: str, stage_dir: str | Path,
                      script: str | None) -> tuple[list[str], list[str]]:
        """上传阶段输入，返回 ``(uploaded, skipped)``。"""
        stage = normalize_stage(stage)
        directory = Path(stage_dir)
        if not directory.is_dir():
            raise StageError("阶段目录不存在: %s（先跑 prepare）" % directory)

        if stage == "aitranss_te":
            # 只上传三个 TE 输入，在 aims_scf 的远端目录里跑。
            uploaded = [n for n in SMJ_AITRANSS_UPLOADS if (directory / n).is_file()]
            skipped = [n for n in SMJ_AITRANSS_UPLOADS if n not in uploaded]
        else:
            skip = skip_uploads(script)
            uploaded = []
            skipped = []
            for entry in sorted(directory.iterdir()):
                if not entry.is_file():
                    continue
                if entry.name in skip:
                    skipped.append(entry.name)
                else:
                    uploaded.append(entry.name)

        if script and script not in uploaded:
            raise StageError(
                "阶段脚本 %s 不在 %s：先跑 `smj stage %s --prepare` 生成作业脚本"
                % (script, directory, stage))

        remote_dir = self.stage_dir(stage)
        for name in uploaded:
            result = self.transport.put(directory / name, remote_join(remote_dir, name))
            if not result.ok:
                raise TransportError(
                    "上传 %s 失败 (rc=%d): %s"
                    % (name, result.code, (result.err or result.out).strip()[:500]))
        return uploaded, skipped

    # -- 提交 -------------------------------------------------------------- #
    def submit(self, stage: str, stage_dir: str | Path | None = None,
               *, script: str | None = None,
               sbatch_args: Sequence[str] = ()) -> Submission:
        """清理 → 上传 → sbatch，返回 :class:`Submission`。"""
        stage = normalize_stage(stage)
        if getattr(self.ctx, "dry_run", False):
            raise StageError(
                "dry-run 模式不提交远端作业（只生成输入）；"
                "去掉 --dry-run 或把 run.mode 改成 slurm-ssh")
        if stage_dir is None:
            stage_dir = self.ctx.layout.stage_dir(stage)
        script = script or stage_script(stage)
        if not script:
            raise StageError("阶段 %s 没有作业脚本（本地阶段无需提交）" % stage)

        remote_dir = self.stage_dir(stage)
        self.ensure_dir(remote_dir)
        cleaned = self.clean_stale(stage, script)
        uploaded, skipped = self.upload_inputs(stage, stage_dir, script)

        args = list(SMJ_SCRIPT_ARGS.get(script, ()))
        # sbatch_args 用于"只重跑数组作业里失败的槽位":命令行 --array=... 会
        # 覆盖脚本里的 #SBATCH --array 行,所以不必重写脚本。
        extra = " ".join(shlex.quote(str(a)) for a in sbatch_args)
        command = "cd %s && sbatch%s %s" % (
            shlex.quote(remote_dir),
            (" " + extra) if extra else "",
            shlex.quote(script))
        if args:
            command += " " + " ".join(shlex.quote(str(a)) for a in args)
        result = self._run(command, timeout=120)
        if result.code == 124:
            # 超时 ≠ 提交失败：作业可能已经创建。必须让调用方保留“未知”状态，
            # 人工 squeue 核对后再决定重提，否则同一 remote_dir 会出现两个作业。
            raise SubmitUncertain(
                "sbatch 超时（120 s），无法确认作业是否已创建（远端目录 %s）。\n"
                "请先在集群上按该目录/提交时间 `squeue`/`sacct` 核对；确认没有作业后"
                "再 `smj stage <name> --submit --force`。" % remote_dir,
                remote_dir=remote_dir)
        if not result.ok:
            raise TransportError(
                "sbatch 提交失败 (rc=%d): %s"
                % (result.code, (result.err or result.out).strip()[:1000]))
        job_id = parse_job_id(result.out) or parse_job_id(result.err)
        if not job_id:
            raise SubmitUncertain(
                "sbatch 已执行但未解析到 job id（输出: %r）。\n"
                "作业可能已经排上：请先 `squeue` 确认；确认没有作业后再重提交。"
                % ((result.out or result.err or "").strip()[:500],),
                remote_dir=remote_dir)
        return Submission(
            stage=stage, mode="slurm-ssh", script=script, remote_dir=remote_dir,
            job_id=job_id, sbatch=command, uploaded=uploaded,
            skipped=skipped, cleaned=cleaned)

    # -- 轮询 -------------------------------------------------------------- #
    def wait(self, stage: str, job_id: str, *, poll_s: float | None = None,
             timeout_s: float | None = None,
             sleep: Callable[[float], Any] | None = None,
             clock: Callable[[], float] | None = None,
             on_tick: Callable[[int, int], Any] | None = None) -> dict[str, Any]:
        """轮询单个作业，返回 ``{finished, failed, exit_code, note, elapsed_s, polls}``。"""
        jobs = [{
            "stage": normalize_stage(stage),
            "job_id": str(job_id),
            "remote_dir": self.stage_dir(stage),
        }]
        return self.wait_many(jobs, poll_s=poll_s, timeout_s=timeout_s,
                              sleep=sleep, clock=clock, on_tick=on_tick)[0]

    def wait_many(self, jobs: Sequence[Mapping[str, Any]], *,
                  poll_s: float | None = None, timeout_s: float | None = None,
                  sleep: Callable[[float], Any] | None = None,
                  clock: Callable[[], float] | None = None,
                  on_tick: Callable[[int, int], Any] | None = None) -> list[dict[str, Any]]:
        """批量轮询（沿用 jobwait 的语义，见模块 docstring）。"""
        pending_jobs = [dict(job) for job in jobs]
        if not pending_jobs:
            return []

        cluster = getattr(getattr(self.ctx, "config", None), "cluster", None)
        default_poll = float(getattr(cluster, "poll_seconds", 60.0) or 60.0)
        default_timeout = float(getattr(cluster, "timeout_hours", 6.0) or 0.0) * 3600.0
        poll = default_poll if poll_s is None else float(poll_s)
        timeout = default_timeout if timeout_s is None else float(timeout_s)

        kwargs: dict[str, Any] = {"poll_s": poll, "timeout_s": timeout}
        if sleep is not None:
            kwargs["sleep"] = sleep
        if clock is not None:
            kwargs["clock"] = clock
        if on_tick is not None:
            kwargs["on_tick"] = on_tick

        results = _jobwait.wait_many(_ProbeAdapter(self.transport), pending_jobs, **kwargs)
        merged: list[dict[str, Any]] = []
        for job, result in zip(pending_jobs, results):
            item = dict(result)
            item.setdefault("stage", job.get("stage"))
            item.setdefault("remote_dir", job.get("remote_dir"))
            merged.append(item)
        return merged

    # -- 回取 -------------------------------------------------------------- #
    def download_plan(self, stage: str) -> list[str]:
        """该阶段应**全量**回取的文件清单（vendored retrieval_plan 的 ``must``）。"""
        plan = _artifacts.retrieval_plan(normalize_stage(stage))
        return list(plan.get("must") or ())

    def download_plan_detail(self, stage: str) -> dict[str, Any]:
        """完整回取计划：``must``（全量）+ ``head``（只取表头行数）+ ``skip``。

        ``optional_head`` 列闭壳层 run 不存在的候选（自旋极化的
        ``alpha.aims``/``beta.aims``），缺席不记入 ``missing``。

        2026-11-11 本地补丁 F1：``mos.aims`` 属 head —— 旧实现只消费 ``must``，
        ``head`` 里的 ``mos.aims`` 永远不回取，导致 aitranss_te 阶段读不到 nsaos
        而在第 5 阶段必断（P0）。这里保留完整计划供 :meth:`download` 消费。
        """
        canonical = normalize_stage(stage)
        # vendored 计划仍以 gaussian_opt 命名首段；molecule_opt 沿用它的
        # Gaussian 产物清单（AIMS 引擎在 runner 里显式传 names，不走到这里）。
        plan_name = "gaussian_opt" if canonical == "molecule_opt" else canonical
        plan = _artifacts.retrieval_plan(plan_name)
        if plan.get("unknown"):
            raise StageError(plan.get("reason") or "未知阶段: %s" % stage)
        head = {str(k): int(v) for k, v in (plan.get("head") or {}).items()}
        optional_head: list[str] = []
        if canonical == "aims_scf":
            # 2026-11-19 结级自旋：FHI-aims 在 spin collinear 下把 mos.aims 拆成
            # alpha.aims + beta.aims；vendored 清单只认 mos.aims，本地只需表头
            # 读 nsaos（aitranss_te.prepare 的 $nsaos），不整取 GB 级波函数。
            # 闭壳层 run 没有这两个文件，缺席不算缺件。
            for name in ("alpha.aims", "beta.aims"):
                if name not in head:
                    head[name] = 2
                    optional_head.append(name)
        return {
            "stage": normalize_stage(stage),
            "must": list(plan.get("must") or ()),
            "head": head,
            "optional_head": optional_head,
            "skip": list(plan.get("skip") or ()),
        }

    def download(self, stage: str, local_dir: str | Path, *,
                 names: Iterable[str] | None = None) -> dict[str, Any]:
        """把阶段产物从远端拷回 ``local_dir``（缺文件只记录，不抛错）。

        ``must`` 全量 ``get``；``head`` 只取前 N 行（``get_head``），避免为读一个
        ``nsaos`` 而拉回 0.9 GB 的 mos.aims。``names`` 显式指定时按全量处理。
        """
        stage = normalize_stage(stage)
        remote_dir = self.stage_dir(stage)
        optional: set[str] = set()
        if names is not None:
            full = list(names)
            head: dict[str, int] = {}
        else:
            detail = self.download_plan_detail(stage)
            full = detail["must"]
            head = detail["head"]
            optional = set(detail.get("optional_head") or ())
        target = Path(local_dir)
        target.mkdir(parents=True, exist_ok=True)
        fetched: list[str] = []
        head_fetched: list[str] = []
        missing: list[str] = []
        optional_missing: list[str] = []
        errors: dict[str, str] = {}
        for name in full:
            result = self.transport.get(remote_join(remote_dir, name), target / name)
            if result.ok:
                fetched.append(name)
            else:
                missing.append(name)
                errors[name] = (result.err or result.out).strip()[:300]
        for name, lines in head.items():
            result = self.transport.get_head(remote_join(remote_dir, name),
                                             target / name, int(lines))
            if result.ok:
                head_fetched.append(name)
            elif name in optional:
                # 闭壳层 run 没有 alpha.aims/beta.aims：缺席是正常状态，不算缺件。
                optional_missing.append(name)
            else:
                missing.append(name)
                errors[name] = (result.err or result.out).strip()[:300]
        if (stage == "aims_scf" and "mos.aims" in missing
                and {"alpha.aims", "beta.aims"}.issubset(set(head_fetched))):
            # 2026-11-19 结级自旋：spin collinear 的 SCF 不写 mos.aims，
            # 而是写 alpha.aims + beta.aims。两者都取到表头时，mos.aims
            # 缺席是自旋极化的正常状态，不该判成回取缺件。
            missing.remove("mos.aims")
            optional_missing.append("mos.aims")
            errors.pop("mos.aims", None)
        return {
            "stage": stage,
            "remote_dir": remote_dir,
            "fetched": fetched,
            "head_fetched": head_fetched,
            "missing": missing,
            "optional_missing": optional_missing,
            "errors": errors,
        }


# --------------------------------------------------------------------------- #
# local 模式执行
# --------------------------------------------------------------------------- #
def run_script_local(transport: Transport, stage_dir: str | Path, script: str,
                     args: Sequence[str] = (), *,
                     timeout: float | None = None) -> CommandResult:
    """在本地 stage dir 里执行作业脚本（``run.mode: local``）。"""
    directory = Path(stage_dir)
    if not (directory / script).is_file():
        raise StageError("阶段脚本缺失: %s（先跑 prepare）" % (directory / script))
    command = "bash %s" % shlex.quote(script)
    if args:
        command += " " + " ".join(shlex.quote(str(a)) for a in args)
    return transport.run(command, cwd=directory, timeout=timeout)


__all__ = [
    "DONE_SENTINEL",
    "SMJ_AITRANSS_UPLOADS",
    "SMJ_SCRIPTS",
    "SMJ_SCRIPT_ARGS",
    "SMJ_STAGE_SCRIPTS",
    "SMJ_STALE_UPLOADS",
    "RemoteSession",
    "Submission",
    "detect_stage_script",
    "parse_job_id",
    "parse_sentinel",
    "read_local_sentinel",
    "run_script_local",
    "skip_uploads",
    "stage_script",
    "stale_names",
]
