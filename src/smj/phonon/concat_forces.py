"""把数组作业每任务的 ``forces.json`` 汇总成一个文件（远端 stdlib 脚本）。

为什么要有这一步
----------------
一个 324 位移的数组作业会产生 324 个 ``fc/%04d/forces.json``。回取时逐个
scp 既慢又容易漏（漏一个，力常数矩阵就少一行，而 :func:`smj.phonon.fc.assemble_fc`
只会在"条数不对"时报错，报错信息离现场很远）。所以让数组作业的收尾步骤在
远端先拼成单个 ``forces_all.json``，回取计划里就只有一个文件。

输出格式::

    {"n_atoms": 90, "count": 324, "amplitude_A": 0.01,
     "forces": [[[fx, fy, fz], ...], ...]}   # [slot][atom][axis]

只依赖标准库: 该文件会被拷到集群上、用节点上的 ``python`` 直接跑。
"""

from __future__ import annotations

import json
import sys

USAGE = "用法: python concat_forces.py <fc_dir> <count> <out.json>"


def concat(fc_dir: str, count: int) -> dict:
    """按 slot 顺序读 ``<fc_dir>/%04d/forces.json`` 并拼成一个 payload。"""
    import os

    forces = []
    n_atoms = None
    for slot in range(count):
        path = os.path.join(fc_dir, "%04d" % slot, "forces.json")
        if not os.path.exists(path):
            raise SystemExit("数组作业没跑完: %s 不存在" % path)
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
    return {"n_atoms": n_atoms, "count": count, "forces": forces}


def main(argv=None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in ("-h", "--help"):
        print(USAGE)
        return 0 if args else 2
    if len(args) != 3:
        print(USAGE, file=sys.stderr)
        return 2
    fc_dir, count_text, target = args
    payload = concat(fc_dir, int(count_text))
    with open(target, "w", encoding="utf-8") as handle:
        json.dump(payload, handle)
    return 0


if __name__ == "__main__":  # pragma: no cover - 远端入口
    sys.exit(main())
