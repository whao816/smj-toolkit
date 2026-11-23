#!/usr/bin/env python
"""给一个原子加一个位移,写出该位移点的 ``geometry.in``。

**只用标准库**:这个脚本要跟着作业一起传到集群,目标机上的 python 不一定有
numpy;而且它只做"改一行坐标"这一件事。

输入:

* ``--geometry``:thermal_device/thermal_relax 产出的 FHI-aims ``geometry.in``
  (带 ``lattice_vector`` 与若干 ``constrain_relaxation`` 行);
* ``--manifest``:``displacements.json``,元素形如
  ``{"slot": 0, "atom": 3, "axis": 0, "axis_name": "x", "sign": 1, "delta_A": 0.01}``;
* ``--index``:取哪一条(slot 与数组下标一致)。

输出:``--out`` 指向的 ``geometry.in``,除目标原子那一行外**逐字节抄写原文件**。
``lattice_vector`` 与 ``constrain_relaxation`` 行原样保留 —— 位移只改一个原子的
三个分量,不动任何约束。
"""
from __future__ import annotations

import argparse
import json
import sys


def load_manifest(path):
    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, list) or not data:
        raise ValueError("displacements.json 应为非空数组")
    return data


def displace(text, entry):
    """把 ``entry`` 描述的位移加到 ``text`` 里第 ``atom`` 个原子(0-based)上。"""
    atom = int(entry["atom"])
    axis = int(entry["axis"])
    delta = float(entry["sign"]) * float(entry["delta_A"])
    if not 0 <= axis <= 2:
        raise ValueError("axis 必须是 0/1/2,实得 %r" % (entry["axis"],))

    lines = text.splitlines()
    out = []
    seen = 0
    hit = False
    for line in lines:
        parts = line.split()
        if len(parts) >= 5 and parts[0] == "atom":
            if seen == atom:
                coords = [float(parts[1]), float(parts[2]), float(parts[3])]
                coords[axis] += delta
                rest = " ".join(parts[4:])
                out.append("atom %22.14f %22.14f %22.14f  %s"
                           % (coords[0], coords[1], coords[2], rest))
                hit = True
            else:
                out.append(line)
            seen += 1
        else:
            out.append(line)
    if not hit:
        raise ValueError("geometry 里找不到第 %d 个 atom 行(共 %d 个)" % (atom, seen))
    return "\n".join(out) + "\n", seen


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    manifest = load_manifest(args.manifest)
    if not 0 <= args.index < len(manifest):
        print("index %d 超出 0..%d" % (args.index, len(manifest) - 1),
              file=sys.stderr)
        return 1
    entry = manifest[args.index]
    with open(args.geometry, "r", encoding="utf-8") as handle:
        text = handle.read()
    shifted, n_atoms = displace(text, entry)
    expected = entry.get("n_atoms")
    if expected is not None and int(expected) != n_atoms:
        print("位移点记录的 n_atoms=%s 与 geometry 的 %d 不一致"
              % (expected, n_atoms), file=sys.stderr)
        return 1
    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write(shifted)
    return 0


if __name__ == "__main__":
    sys.exit(main())
