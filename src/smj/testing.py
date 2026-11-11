"""测试/演示用的假传输层：不碰网络、不碰真实文件系统，只记录命令。

典型用法（用户也可以用它在本机演练配置，不连集群）::

    from smj.testing import FakeTransport

    class MyRun(FakeTransport):
        def on_run(self, command, cwd, timeout):
            if command.startswith("sbatch"):
                return CommandResult(0, "Submitted batch job 12345")
            return super().on_run(command, cwd, timeout)
"""

from __future__ import annotations

from pathlib import Path

from .transport import CommandResult, Transport


class FakeTransport(Transport):
    """内存传输层。

    - ``files``: 远端路径 -> 内容（文本；二进制文件按 bytes 保存）。
    - ``commands``: 依次记录 ``run/put/get``。
    - ``on_run(command, cwd, timeout)``: 钩子，返回 ``CommandResult`` 则采用，
      返回 ``None`` 走默认行为（成功、空输出）。
    """

    name = "fake"

    def __init__(self) -> None:
        self.files: dict[str, bytes] = {}
        self.commands: list[tuple[str, str]] = []
        self.cwd = "/fake/remote"
        self.default_rc = 0

    # -- 钩子 -------------------------------------------------------------- #
    def on_run(self, command: str, cwd: str | Path | None, timeout: float | None) -> CommandResult | None:
        return None

    # -- 接口 -------------------------------------------------------------- #
    def run(self, command: str, *, timeout: float | None = None,
            cwd: str | Path | None = None) -> CommandResult:
        self.commands.append(("run", command))
        result = self.on_run(command, cwd, timeout)
        if result is not None:
            return result
        return CommandResult(self.default_rc, "", "", command=command)

    def put(self, local: str | Path, remote: str) -> CommandResult:
        src = Path(local)
        self.commands.append(("put", "%s -> %s" % (src, remote)))
        try:
            self.files[_norm(remote)] = src.read_bytes()
        except OSError as exc:
            return CommandResult(1, "", str(exc), command="put %s" % src)
        return CommandResult(0, "", "", command="put %s" % src)

    def get(self, remote: str, local: str | Path) -> CommandResult:
        self.commands.append(("get", "%s -> %s" % (remote, local)))
        data = self.files.get(_norm(remote))
        if data is None:
            return CommandResult(1, "", "not found: %s" % remote, command="get %s" % remote)
        dst = Path(local)
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_bytes(data)
        return CommandResult(0, "", "", command="get %s" % remote)

    def get_head(self, remote: str, local: str | Path, lines: int = 2) -> CommandResult:
        """模拟 ``head -n N``：只写前 N 行，用来验证 GB 级文件不会被整份拉回。"""
        self.commands.append(("get_head", "%s -> %s (%d 行)" % (remote, local, lines)))
        data = self.files.get(_norm(remote))
        if data is None:
            return CommandResult(1, "", "not found: %s" % remote, command="head %s" % remote)
        text = data.decode("utf-8", "replace").splitlines(keepends=True)[:int(lines)]
        dst = Path(local)
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text("".join(text), encoding="utf-8")
        return CommandResult(0, "", "", command="head %s" % remote)

    def mkdirs(self, remote: str) -> CommandResult:
        self.commands.append(("mkdirs", remote))
        return CommandResult(0, "", "", command="mkdir -p %s" % remote)

    def exists(self, remote: str) -> bool:
        return _norm(remote) in self.files

    # -- 测试辅助 ---------------------------------------------------------- #
    def simulate_done(self, remote_dir: str, rc: int = 0) -> None:
        """写 ``.dsh_done`` 哨兵（模拟作业结束）。"""
        self.files[_norm(remote_dir + "/.dsh_done")] = ("%d 2026-11-11T00:00:00+08:00\n" % rc).encode()

    def add_text(self, remote: str, text: str) -> None:
        self.files[_norm(remote)] = text.encode("utf-8")

    def text(self, remote: str) -> str | None:
        data = self.files.get(_norm(remote))
        return data.decode("utf-8", "replace") if data is not None else None

    def names(self) -> list[str]:
        return sorted(self.files)


def _norm(remote: str) -> str:
    text = str(remote).replace("\\", "/")
    while "//" in text:
        text = text.replace("//", "/")
    return text.rstrip("/") if text != "/" else text
