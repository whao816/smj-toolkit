"""smj 的异常体系。

约定：
- 所有异常都继承 :class:`SmjError`，CLI 顶层捕获它 -> 一行可读错误 + 退出码 1
  （``--debug`` 时再抛完整 traceback）。
- :class:`StageBlocked` 表示"当前阶段的前置条件不满足"，属于**可解释的拒绝**
  而不是 bug：CLI 打印原因与缺失文件路径，建议用户先跑上游阶段或显式传参。
"""

from __future__ import annotations


class SmjError(Exception):
    """smj 错误的基类。"""


class ConfigError(SmjError):
    """smj.yaml 不存在 / YAML 语法错误 / 字段值非法。"""


class DependencyMissing(SmjError):
    """运行所需的第三方包或外部命令缺失（附安装/配置提示）。"""


class TransportError(SmjError):
    """SSH/SCP/Slurm 命令执行失败或远端路径不可用。"""


class StageError(SmjError):
    """阶段内部错误（输入生成失败、产物解析失败、作业失败等）。"""


class StageBlocked(StageError):
    """阶段无法开始：上游产物缺失或参数不足。

    Parameters
    ----------
    stage:
        被阻塞的阶段名（标准名，如 ``aims_scf``）。
    reason:
        人读原因。
    missing:
        缺失文件/目录的绝对路径（可选，用于 CLI 直接列出）。
    """

    def __init__(self, stage: str, reason: str, missing: list[str] | None = None):
        self.stage = stage
        self.reason = reason
        self.missing = list(missing or [])
        detail = reason
        if self.missing:
            detail += "\n缺失: " + ", ".join(self.missing)
        super().__init__("阶段 %s 被阻塞: %s" % (stage, detail))
