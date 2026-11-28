"""画单个 run 的声子透射谱与累积热导，用于论文/报告图.

这是 ``smj.plotting.plot_thermal`` 的命令行封装：图的设计与 ``smj thermal``
阶段、``smj report`` 嵌入的那张完全一致（同一函数），避免"报告里的图和脚本
出的图不是一个风格"。额外把谱数据落成 CSV（表头带单位），便于外部复算。

输入 ``thermal.json``（``smj.thermal.v1``），输出：

* ``<stem>.png`` / ``<stem>.pdf`` —— 两联图：(a) ``T_ph(ω)``（lead 带顶以上
  阴影化），(b) 归一化累积 ``G_ph``（标 25/50/75/90% 分位）
* ``<stem>.csv`` —— ``energy_meV, transmission, cumulative_g_ph_pct``
* ``<stem>.meta.json`` —— 画图参数快照（G_ph、η、分位能量），保证图能追回数值

用法示例::

    python scripts/plot_thermal_spectrum.py \\
        --thermal runs/nap-sh-thermal/nap-sh-thermal/thermal_analyze/thermal.json \\
        --title "naphthalene-2,6-dithiol / Au(111) 3x3"
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--thermal", type=Path, required=True,
                   help="thermal.json 路径（smj.thermal.v1）")
    p.add_argument("--output", type=Path, default=None,
                   help="输出目录（默认与 thermal.json 同目录）")
    p.add_argument("--stem", type=str, default="thermal_spectrum",
                   help="输出文件名主干（默认 thermal_spectrum，与阶段产物同名）")
    p.add_argument("--title", type=str, default=None, help="图题；缺省不写标题")
    p.add_argument("--xmax-meV", type=float, default=22.0,
                   help="横轴上限（meV，默认 22）")
    p.add_argument("--formats", type=str, default="png,pdf",
                   help="输出格式，逗号分隔（默认 png,pdf）")
    p.add_argument("--verbose", action="store_true", help="输出调试日志")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    from smj.plotting import plot_thermal

    if not args.thermal.is_file():
        raise SystemExit("找不到 %s" % args.thermal)
    payload = json.loads(args.thermal.read_text(encoding="utf-8"))
    if payload.get("schema") != "smj.thermal.v1":
        raise SystemExit("schema 应为 smj.thermal.v1，得到 %r"
                         % payload.get("schema"))

    output_dir = args.output or args.thermal.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    formats = tuple(s.strip() for s in args.formats.split(",") if s.strip())

    result = plot_thermal(args.thermal, output_dir / ("%s.png" % args.stem),
                          formats=formats, title=args.title)
    if not result.get("ok"):
        raise SystemExit("出图失败：%s" % result.get("reason"))
    for path in result["written"]:
        logger.info("写入 %s", path)

    # CSV：用 thermal.json 自带的累积曲线，保证与阶段内数字同源
    spectrum = payload["spectrum"]
    csv_path = output_dir / ("%s.csv" % args.stem)
    rows = ["energy_meV,transmission,cumulative_g_ph_pw_per_k"]
    total = float(spectrum["cumulative_g_ph_pw_per_k"][-1]) or 1.0
    for energy, t_ph, cumulative in zip(spectrum["energy_meV"],
                                        spectrum["transmission"],
                                        spectrum["cumulative_g_ph_pw_per_k"]):
        if float(energy) > args.xmax_meV:
            break
        rows.append("%.6f,%.8e,%.8e"
                    % (float(energy), float(t_ph), float(cumulative)))
    csv_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    logger.info("写入 %s", csv_path)

    meta = {
        "source": str(args.thermal.resolve()),
        "title": args.title,
        "xmax_meV": args.xmax_meV,
        "n_atoms": payload.get("n_atoms"),
        "temperature_K": payload.get("temperature_K"),
        "g_ph_pw_per_k": payload.get("g_ph_pw_per_k"),
        "g_quantum_pw_per_k": payload.get("g_quantum_pw_per_k"),
        "n_channels_effective": payload.get("n_channels_effective"),
        "band_top_meV": payload.get("band_top_meV"),
        "eta_meV": payload.get("eta_meV"),
        "energy_step_meV": payload.get("energy_step_meV"),
        "quantile_energies_meV": result.get("quantiles"),
        "written": list(result["written"]) + [str(csv_path)],
    }
    meta_path = output_dir / ("%s.meta.json" % args.stem)
    meta_path.write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n",
                         encoding="utf-8")
    logger.info("写入 %s", meta_path)
    logger.info("分位: %s", result.get("quantiles"))
    logger.info("G_ph(%.0f K) = %.3f pW/K（累计 %.6e pW/K）",
                float(payload.get("temperature_K", 300.0)),
                float(payload.get("g_ph_pw_per_k", 0.0)), total)


if __name__ == "__main__":
    main()
