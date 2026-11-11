"""传输层：把命令、上传、下载抽象成两个实现。

- :class:`LocalTransport` —— 本机 ``bash -lc``（``run.mode: local``，AIMS/Gaussian
  装在当前机器上的情形）。
- :class:`SshTransport` —— OpenSSH 客户端（``run.mode: slurm-ssh``，真实计算在集群侧）；
  上传/下载用 ``scp``，同样的选项（端口、私钥、BatchMode）也用于 ``ssh``。

只暴露执行原语（run/put/get/mkdirs/exists），**不含 Slurm 语义**；Slurm 的提交、
状态轮询、完成判定在 :mod:`smj.remote` 里，便于用假传输层做离线测试。
"""

from __future__ import annotations

import posixpath
import shlex
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

from .errors import TransportError


@dataclass
class CommandResult:
    """一条命令的执行结果。``code`` 为 124 表示超时（与 timeout(1) 一致）。"""

    code: int
    out: str = ""
    err: str = ""
    command: str = ""

    @property
    def ok(self) -> bool:
        return self.code == 0

    def expect_ok(self, what: str = "命令") -> "CommandResult":
        if not self.ok:
            raise TransportError("%s失败 (rc=%d): %s\n%s" % (
                what, self.code, self.command or "", (self.err or self.out).strip()[:2000]))
        return self

    def __str__(self) -> str:  # pragma: no cover - 调试/日志
        return "rc=%d %s" % (self.code, (self.out or self.err).strip()[:200])


class Transport:
    """传输层接口。"""

    name = "transport"

    def run(self, command: str, *, timeout: float | None = None,
            cwd: str | Path | None = None) -> CommandResult:
        raise NotImplementedError

    def put(self, local: str | Path, remote: str) -> CommandResult:
        raise NotImplementedError

    def get(self, remote: str, local: str | Path) -> CommandResult:
        raise NotImplementedError

    def mkdirs(self, remote: str) -> CommandResult:
        return self.run("mkdir -p %s" % shlex.quote(remote))

    def exists(self, remote: str) -> bool:
        r = self.run("test -e %s" % shlex.quote(remote))
        return r.ok

    def put_many(self, files: Sequence[tuple[str | Path, str]]) -> None:
        for local, remote in files:
            self.put(local, remote).expect_ok("上传 %s" % local)

    def get_many(self, files: Sequence[tuple[str, str | Path]]) -> None:
        for remote, local in files:
            self.get(remote, local).expect_ok("下载 %s" % remote)


# --------------------------------------------------------------------------- #
# 本地
# --------------------------------------------------------------------------- #
class LocalTransport(Transport):
    """本机 bash 执行（Windows 上需要 WSL/Git-Bash 提供的 ``bash``）。"""

    name = "local"

    def __init__(self, bash: str | None = None):
        self.bash = bash or shutil.which("bash")
        if not self.bash:
            raise TransportError(
                "run.mode=local 需要 bash，但当前 PATH 中找不到: 请安装 Git-Bash/WSL，"
                "并把 bash 加入 PATH；只想检查输入时用 run.mode=dry-run")

    def _run(self, argv: list[str], *, timeout: float | None, cwd: str | Path | None) -> CommandResult:
        try:
            proc = subprocess.run(
                argv, cwd=str(cwd) if cwd else None, timeout=timeout,
                capture_output=True, text=True, errors="replace",
            )
        except subprocess.TimeoutExpired as exc:
            return CommandResult(124, "", "本地命令超时 (%ss): %s" % (timeout, " ".join(argv)),
                                 command=" ".join(argv))
        except OSError as exc:
            raise TransportError("本地命令启动失败: %s" % exc) from exc
        return CommandResult(proc.returncode, proc.stdout or "", proc.stderr or "", command=" ".join(argv))

    def run(self, command: str, *, timeout: float | None = None,
            cwd: str | Path | None = None) -> CommandResult:
        return self._run([self.bash, "-lc", command], timeout=timeout, cwd=cwd)

    def put(self, local: str | Path, remote: str) -> CommandResult:
        src, dst = Path(local), Path(remote)
        try:
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
        except OSError as exc:
            return CommandResult(1, "", str(exc), command="cp %s %s" % (src, dst))
        return CommandResult(0, "", "", command="cp %s %s" % (src, dst))

    def get(self, remote: str, local: str | Path) -> CommandResult:
        return self.put(remote, local)  # 本机拷贝是对称的


# --------------------------------------------------------------------------- #
# SSH / SCP
# --------------------------------------------------------------------------- #
class SshTransport(Transport):
    """OpenSSH 客户端；目标形如 ``user@host`` 或 ``~/.ssh/config`` 别名。"""

    name = "ssh"

    def __init__(self, target: str, *, port: int = 22, key: str | None = None,
                 connect_timeout: int = 15, ssh_bin: str | None = None,
                 scp_bin: str | None = None):
        if not target:
            raise TransportError("SshTransport 需要 cluster.host（或 user@host）")
        self.target = target
        self.port = port
        self.key = str(key) if key else None
        self.connect_timeout = int(connect_timeout)
        self.ssh_bin = ssh_bin or shutil.which("ssh")
        self.scp_bin = scp_bin or shutil.which("scp")
        if not self.ssh_bin:
            raise TransportError("找不到 ssh 命令: 请安装 OpenSSH 客户端并加入 PATH")

    # -- 基础命令 ---------------------------------------------------------- #
    def _opts(self) -> list[str]:
        opts = ["-o", "BatchMode=yes", "-o", "ConnectTimeout=%d" % self.connect_timeout]
        if self.port and int(self.port) != 22:
            opts += ["-P" if self._is_scp else "-p", str(int(self.port))]
        if self.key:
            opts += ["-i", self.key]
        return opts

    _is_scp = False

    def _ssh_argv(self, remote_command: str) -> list[str]:
        return [self.ssh_bin, *self._opts(), self.target, remote_command]  # type: ignore[list-item]

    def _scp_argv(self, *args: str) -> list[str]:
        opts = []
        if self.port and int(self.port) != 22:
            opts += ["-P", str(int(self.port))]
        if self.key:
            opts += ["-i", self.key]
        opts += ["-o", "BatchMode=yes", "-o", "ConnectTimeout=%d" % self.connect_timeout]
        return [self.scp_bin or "scp", *opts, *args]  # type: ignore[list-item]

    @staticmethod
    def _check_remote_path(remote: str) -> str:
        if not remote or any(ch.isspace() for ch in remote):
            raise TransportError("远端路径不能为空或含空白: %r" % (remote,))
        return remote

    @staticmethod
    def _remote_spec(target: str, remote: str) -> str:
        return "%s:%s" % (target, remote)

    def _exec(self, argv: list[str], *, timeout: float | None) -> CommandResult:
        try:
            proc = subprocess.run(argv, timeout=timeout, capture_output=True,
                                  text=True, errors="replace")
        except subprocess.TimeoutExpired:
            return CommandResult(124, "", "ssh/scp 超时 (%ss)" % timeout, command=" ".join(argv))
        except OSError as exc:
            raise TransportError("ssh/scp 启动失败: %s" % exc) from exc
        return CommandResult(proc.returncode, proc.stdout or "", proc.stderr or "", command=" ".join(argv))

    # -- 接口 -------------------------------------------------------------- #
    def run(self, command: str, *, timeout: float | None = None,
            cwd: str | Path | None = None) -> CommandResult:
        remote_cmd = command
        if cwd is not None:
            remote_cmd = "cd %s && %s" % (shlex.quote(str(cwd)), command)
        return self._exec(self._ssh_argv(remote_cmd), timeout=timeout)

    def put(self, local: str | Path, remote: str) -> CommandResult:
        remote = self._check_remote_path(remote)
        return self._exec(self._scp_argv("-q", str(local), self._remote_spec(self.target, remote)),
                          timeout=600)

    def get(self, remote: str, local: str | Path) -> CommandResult:
        remote = self._check_remote_path(remote)
        Path(local).parent.mkdir(parents=True, exist_ok=True)
        return self._exec(self._scp_argv("-q", self._remote_spec(self.target, remote), str(local)),
                          timeout=600)

    def check(self, timeout: float = 30) -> None:
        """一次真实往返：确认 SSH 与远程 shell 可用。"""
        r = self.run("true", timeout=timeout)
        if not r.ok:
            raise TransportError(
                "SSH 连接失败 (%s): %s\n提示: 先手工 `ssh %s true` 确认免密登录与别名"
                % (self.target, (r.err or r.out).strip(), self.target))


# --------------------------------------------------------------------------- #
# 工厂
# --------------------------------------------------------------------------- #
def make_transport(config) -> Transport:
    """按 ``run.mode`` 创建传输层（dry-run 不创建）。"""
    mode = getattr(config, "mode", None)
    if mode == "slurm-ssh":
        cluster = config.cluster
        return SshTransport(cluster.ssh_target, port=cluster.port, key=cluster.key)
    if mode == "local":
        return LocalTransport()
    raise TransportError("run.mode=%r 不需要传输层（dry-run 只生成输入）" % (mode,))


def remote_join(root: str, *parts: str) -> str:
    """拼接远端 POSIX 路径（忽略空段）。"""
    path = root
    for part in parts:
        if part:
            path = posixpath.join(path, str(part))
    return path
