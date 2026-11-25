"""把数组作业每任务的 ``forces.json`` 汇总成一个文件（远端 stdlib 脚本）。

为什么要有这一步
----------------
一个 324 位移的数组作业会产生 324 个 ``fc/%04d/forces.json``。回取时逐个
scp 既慢又容易漏（漏一个，力常数矩阵就少一行，而 :func:`smj.phonon.fc.assemble_fc`
只会在"条数不对"时报错，报错信息离现场很远）。所以让数组作业的收尾步骤在
远端先拼成单个 ``forces_all.json``，回取计划里就只有一个文件。

顺带做两件在远端更划算的事（本地没有 ``aims.out``）：

1. 逐槽位判 SCF 是否收敛（``Self-consistency cycle converged`` 恰出现一次）。
   不收敛但照样产出了力的槽位，是"数看着有、其实不该用"的最坏情况；
2. 把缺槽位、未收敛槽位写进 ``fc_summary.json``，让 collect 的报错直接给出
   可执行的重试清单（``smj thermal --stage thermal_fc --retry 17,203``）。

输出格式::

    {"n_atoms": 90, "count": 324, "amplitude_A": 0.01,
     "forces": [[[fx, fy, fz], ...], ...]}   # [slot][atom][axis]

只依赖标准库: 该文件会被拷到集群上、用节点上的 ``python`` 直接跑。
"""

from __future__ import annotations

import json
import sys

USAGE = (
    "用法: python concat_forces.py <fc_dir> <count> <out.json> "
    "[--summary summary.json] [--allow-missing]"
)

CONVERGED = "Self-consistency cycle converged"
NOT_CONVERGED = "Self-consistency cycle NOT converged"


def scan_aims(text: str) -> dict:
    """从一份 ``aims.out`` 里读 SCF 收敛情况（纯文本，只认两个固定标记）。"""
    n_cycles = 0
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("Self-consistency cycle"):
            n_cycles += 1
    converged = (CONVERGED in text) and (NOT_CONVERGED not in text)
    return {"converged": bool(converged), "n_cycles": n_cycles}


def _read_summary(fc_dir: str, slot: int) -> dict:
    import os

    entry = {"slot": slot}
    aims = os.path.join(fc_dir, "%04d" % slot, "aims.out")
    if os.path.exists(aims):
        try:
            with open(aims, "r", encoding="utf-8", errors="replace") as handle:
                entry.update(scan_aims(handle.read()))
        except OSError as exc:  # pragma: no cover - 远端 I/O 失败
            entry.update({"converged": None, "n_cycles": 0, "error": str(exc)})
    else:
        entry.update({"converged": None, "n_cycles": 0, "error": "aims.out 缺失"})
    return entry


def concat(fc_dir: str, count: int, *, allow_missing: bool = False,
           with_summary: bool = False) -> dict:
    """按 slot 顺序读 ``<fc_dir>/%04d/forces.json`` 并拼成一个 payload。

    ``allow_missing=True`` 时保留缺槽位（以 ``None`` 占位）并把情况写进
    ``payload["summary"]``；默认缺一个就退出（早期行为）。
    """
    import os

    forces = []
    n_atoms = None
    missing = []
    slots = []
    for slot in range(count):
        path = os.path.join(fc_dir, "%04d" % slot, "forces.json")
        if with_summary or allow_missing:
            slots.append(_read_summary(fc_dir, slot))
        if not os.path.exists(path):
            if not allow_missing:
                raise SystemExit("数组作业没跑完: %s 不存在" % path)
            missing.append(slot)
            forces.append(None)
            continue
        with open(path, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
        block = payload["forces"]
        if n_atoms is None:
            n_atoms = int(payload.get("n_atoms") or len(block))
        elif int(payload.get("n_atoms") or len(block)) != n_atoms:
            raise SystemExit(
                "forces.json 原子数不一致: slot %d 是 %d,前面是 %d"
                % (slot, len(block), n_atoms))
        forces.append(block)
    out = {"n_atoms": n_atoms, "count": count, "forces": forces}
    if with_summary or allow_missing:
        not_converged = [s["slot"] for s in slots if s.get("converged") is False]
        out["summary"] = {
            "count": count,
            "missing": missing,
            "not_converged": not_converged,
            "slots": slots,
        }
    return out


def write_summary(payload: dict, target: str) -> None:
    summary = payload.get("summary") or {"count": payload.get("count")}
    with open(target, "w", encoding="utf-8") as handle:
        json.dump(summary, handle)


def main(argv=None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in ("-h", "--help"):
        print(USAGE)
        return 0 if args else 2
    summary_path = None
    allow_missing = False
    positional = []
    index = 0
    while index < len(args):
        item = args[index]
        if item == "--summary" and index + 1 < len(args):
            summary_path = args[index + 1]
            index += 2
            continue
        if item == "--allow-missing":
            allow_missing = True
            index += 1
            continue
        positional.append(item)
        index += 1
    if len(positional) != 3:
        print(USAGE, file=sys.stderr)
        return 2
    fc_dir, count_text, target = positional
    payload = concat(fc_dir, int(count_text),
                     allow_missing=allow_missing,
                     with_summary=bool(summary_path) or allow_missing)
    with open(target, "w", encoding="utf-8") as handle:
        json.dump(payload, handle)
    if summary_path:
        write_summary(payload, summary_path)
    summary = payload.get("summary")
    if summary and (summary.get("missing") or summary.get("not_converged")):
        print("缺失槽位: %s" % (summary.get("missing") or "-"), file=sys.stderr)
        print("未收敛槽位: %s" % (summary.get("not_converged") or "-"),
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":  # pragma: no cover - 远端入口
    sys.exit(main())
