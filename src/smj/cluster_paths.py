"""Cluster workspaces are rooted in the SSH login account's home directory."""
import posixpath
import re
import shlex

from .errors import ConfigError, TransportError

DEFAULT_ROOT = "~/smj-runs"
MAX_RUN_NAME = 40


def validate_run_name(name):
    if len(name) > MAX_RUN_NAME or not re.fullmatch(r"[A-Za-z0-9]+(?:-[A-Za-z0-9]+)*", name):
        raise ConfigError("集群 run.name 限 1–40 个英文字母、数字，可用单个连字符分词，如 tpm-b-rad")


def validate_root(root):
    if root in (None, DEFAULT_ROOT, "$HOME/smj-runs", "${HOME}/smj-runs"):
        return
    if (not isinstance(root, str) or not root.startswith("/") or
            not root.endswith("/smj-runs") or any(p in ("", ".", "..") for p in root[1:].split("/"))):
        raise ConfigError("cluster.remote_root 必须是远端登录用户的 ~/smj-runs")


def checked_directory(transport, root, run, stage=None, *, timeout=30):
    """Resolve on the login node before creating or consuming a job directory."""
    validate_root(root)
    validate_run_name(run)
    suffix = "smj-runs/" + run
    if stage:
        suffix += "/" + stage
    # Checking the lexical target against realpath -m also covers missing directories.
    command = 'home=$(cd -P -- "${HOME:?}" && pwd -P) || exit 1; '
    if root and root.startswith("/"):
        command += 'test %s = "$HOME/smj-runs" || test %s = "$home/smj-runs"' % (
            shlex.quote(root), shlex.quote(root))
        command += ' || exit 1; '
    command += 'dir="$home"/%s && test "$(realpath -m -- "$dir")" = "$dir" && ' % shlex.quote(suffix)
    if stage:
        command += ('if test -d "$dir"; then '
                    'link=$(find "$dir" -type l -print -quit) && test -z "$link"; fi && ')
    command += 'printf "%s\\n" "$dir"'
    result = transport.run(command, timeout=timeout)
    directory = result.out.rstrip("\r\n")
    if (not result.ok or not directory.startswith("/") or "\n" in directory or
            not directory.endswith("/" + suffix) or posixpath.normpath(directory) != directory):
        raise TransportError("集群工作目录检查失败：请核对登录用户的 ~/smj-runs、目录权限及符号链接")
    return directory
