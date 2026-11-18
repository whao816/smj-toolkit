"""Azulene E/Z 异构体透射比复核（原仓库根 _azo.py 的整理版）。

用途
----
比较两个 run 的 ``aitranss_te/TE.dat``：在若干能量窗内计算
``T_Z(E) / T_E(E)`` 的均值/标准差/极值，给出 E_F 处的比值，并统计低透射
样本数（T < 阈值，作为"能隙"的粗略代理）。这是对 E/Z 结论的一次性数值复核，
不是流程阶段，也不进入 ``smj`` 包。

TE.dat 列约定
-------------
第 1 列 E [Hartree]，第 2 列 E−E_F [eV]，第 3 列 T(E) per spin；``#`` 开头为注释。

用法
----
    python scripts/azo_ez_ratio.py
    python scripts/azo_ez_ratio.py --e-run runs/azo-e --z-run runs/azo-z
    python scripts/azo_ez_ratio.py --window 0.4 --sub-windows 0.1,0.2 --json out.json

``--e-run`` / ``--z-run`` 既可以给 run 根目录（兼容 ``runs/<name>/<name>/`` 布局）、
``runs/<name>/<name>`` 目录，也可以直接给 ``TE.dat`` 路径。
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

DEFAULT_E_RUN = Path("runs/azo-e")
DEFAULT_Z_RUN = Path("runs/azo-z")
RATIO_FLOOR = 1e-12


def resolve_te(path: Path) -> Path:
    """接受 run 根目录、run.dir/run.name 目录或 TE.dat 路径，返回实际 TE.dat。"""
    path = Path(path)
    if path.is_file():
        return path
    candidates = [path / "aitranss_te" / "TE.dat",
                  path / path.name / "aitranss_te" / "TE.dat"]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise SystemExit("找不到 TE.dat: %s（试过 %s）"
                     % (path, " / ".join(str(c) for c in candidates)))


def load_te(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """返回 (E−E_F [eV], T(E) per spin)；注释行由 numpy 自动跳过。"""
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] < 3:
        raise SystemExit("TE.dat 至少需要 3 列（E、E−E_F、T(E)），实际形状 %r: %s"
                         % (data.shape, path))
    return data[:, 1], data[:, 2]


def ratio_stats(energy: np.ndarray, numerator: np.ndarray, denominator: np.ndarray,
                lo: float, hi: float) -> dict[str, Any]:
    """在 (lo, hi) 窗内统计 T_num / T_den。"""
    mask = (energy > lo) & (energy < hi)
    count = int(mask.sum())
    if count == 0:
        return {"window": [lo, hi], "n": 0}
    ratio = numerator[mask] / np.clip(denominator[mask], RATIO_FLOOR, None)
    return {
        "window": [lo, hi],
        "n": count,
        "mean": float(ratio.mean()),
        "sd": float(ratio.std()),
        "min": float(ratio.min()),
        "max": float(ratio.max()),
    }


def fermi_ratio(energy_e: np.ndarray, trans_e: np.ndarray,
                energy_z: np.ndarray, trans_z: np.ndarray) -> dict[str, float]:
    """E_F（E=0）处各支与比值的线性插值（各自用自己的能量网格）。"""
    at_e = float(np.interp(0.0, energy_e, trans_e))
    at_z = float(np.interp(0.0, energy_z, trans_z))
    return {"T_e": at_e, "T_z": at_z,
            "ratio": at_z / max(at_e, RATIO_FLOOR)}


def low_transmission(energy: np.ndarray, transmission: np.ndarray,
                     *, threshold: float, scan: float) -> int:
    mask = (energy > -scan) & (energy < scan) & (transmission < threshold)
    return int(mask.sum())


def parse_floats(text: str) -> list[float]:
    return [float(part) for part in text.split(",") if part.strip()]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Azulene E/Z 异构体透射比复核")
    parser.add_argument("--e-run", type=Path, default=DEFAULT_E_RUN,
                        help="E 异构体的 run 目录或 TE.dat 路径（默认 %(default)s）")
    parser.add_argument("--z-run", type=Path, default=DEFAULT_Z_RUN,
                        help="Z 异构体的 run 目录或 TE.dat 路径（默认 %(default)s）")
    parser.add_argument("--window", type=float, default=0.6,
                        help="主统计窗 ±window eV（默认 %(default)s）")
    parser.add_argument("--sub-windows", default="0.2,0.4,1.0",
                        help="附加统计窗（逗号分隔，单位 eV）")
    parser.add_argument("--low-t-threshold", type=float, default=0.02,
                        help="低透射阈值（默认 %(default)s）")
    parser.add_argument("--low-t-scan", type=float, default=1.2,
                        help="低透射统计范围 ±scan eV（默认 %(default)s）")
    parser.add_argument("--json", type=Path, default=None,
                        help="把结果同时写入 JSON 文件")
    args = parser.parse_args(argv)

    e_path = resolve_te(args.e_run)
    z_path = resolve_te(args.z_run)
    energy_e, trans_e = load_te(e_path)
    energy_z, trans_z = load_te(z_path)
    if energy_e.shape != energy_z.shape:
        raise SystemExit("两个 TE.dat 的能量网格不一致：%d vs %d 行"
                         % (energy_e.size, energy_z.size))

    windows = [args.window]
    for value in parse_floats(args.sub_windows):
        if value not in windows:
            windows.append(value)
    stats = [ratio_stats(energy_e, trans_z, trans_e, -value, value)
             for value in windows]
    main_stats = stats[0]
    fermi = fermi_ratio(energy_e, trans_e, energy_z, trans_z)
    low_t = {
        "e": low_transmission(energy_e, trans_e,
                              threshold=args.low_t_threshold, scan=args.low_t_scan),
        "z": low_transmission(energy_z, trans_z,
                              threshold=args.low_t_threshold, scan=args.low_t_scan),
    }

    print("E: %s" % e_path)
    print("Z: %s" % z_path)
    if main_stats.get("n"):
        print("ratio Z/E over ±%.2f eV: mean %.3f sd %.3f min %.3f max %.3f (n=%d)"
              % (args.window, main_stats["mean"], main_stats["sd"],
                 main_stats["min"], main_stats["max"], main_stats["n"]))
    else:
        print("ratio Z/E over ±%.2f eV: 窗内无数据" % args.window)
    for item in stats[1:]:
        if not item.get("n"):
            print("  ±%.2f eV: 窗内无数据" % item["window"][1])
            continue
        print("  ±%.2f eV: mean %.3f sd %.3f (n=%d)"
              % (item["window"][1], item["mean"], item["sd"], item["n"]))
    print("T(E_F): E %.4f Z %.4f ratio %.4f"
          % (fermi["T_e"], fermi["T_z"], fermi["ratio"]))
    print("low-T samples (T < %.3g) within ±%.2f eV: E=%d Z=%d"
          % (args.low_t_threshold, args.low_t_scan, low_t["e"], low_t["z"]))

    if args.json is not None:
        payload = {
            "e_te": str(e_path),
            "z_te": str(z_path),
            "main": main_stats,
            "windows": stats,
            "fermi": fermi,
            "low_t": {"threshold": args.low_t_threshold,
                      "scan_eV": args.low_t_scan, **low_t},
        }
        args.json.write_text(
            json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8")
        print("JSON: %s" % args.json)
    return 0


if __name__ == "__main__":
    sys.exit(main())
