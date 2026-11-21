"""junction_analyze：TE.dat → Landauer 输运量（本地阶段，不提交）。

用 vendored ``landauer.analyze_te`` 把 ``TE.dat`` 积分成
log10(G/G0)、Seebeck 系数与电子 ZT，写 ``results.json``。
顺带用统一版式（``smj.plotting``）出 ``T_E.png`` / ``T_E.pdf``；
matplotlib 未安装时只跳过出图，不影响数值结果。
prepare 即产出（纯本地解析）；collect 汇总同一条结果，供 runner 落库。
"""

from __future__ import annotations

import json
from pathlib import Path

from ..junction.landauer import analyze_te
from .. import plotting
from ..errors import StageError
from ..landauer_spin import analyze_te_spin
from ..spin import spin_state
from .base import Collected, Prepared, Stage, register

#: 少于等于该点数时，结果只应被当作管线验证，不足以支撑科学结论。
MIN_PRODUCTION_POINTS = 7


@register
class JunctionAnalyzeStage(Stage):
    name = "junction_analyze"
    depends = ("aitranss_te",)
    needs_submit = False
    input_hint = "消费 aitranss_te/TE.dat（--te 可显式覆盖文件或目录）"
    outputs = ("results.json",)

    # -- prepare（本地阶段：直接产出） -------------------------------------- #
    def prepare(self, ctx, stage_dir: Path) -> Prepared:
        stage_dir = Path(stage_dir)
        te = self._te_path(ctx)
        out_dir = ctx.override_path("out") or stage_dir
        out_dir.mkdir(parents=True, exist_ok=True)
        results_path = out_dir / "results.json"
        temperature = ctx.config.transport.temperature
        try:
            result = analyze_te(te, temperature=temperature, output=results_path)
        except ValueError as exc:
            raise StageError("junction_analyze 无法解析 TE.dat: %s" % exc) from exc
        state = spin_state(ctx.config)
        if state.open_shell:
            # vendored analyze_te 对 α/β 取均值，只保总 G/S；这里保留分通道量。
            result["spin"] = analyze_te_spin(te, temperature=temperature)
            results_path.write_text(
                json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False),
                encoding="utf-8")
        points = int(result.get("points") or 0)
        ctx.log("junction_analyze: %d 个能量点 → %s" % (points, results_path.name))
        notes = [
            "T = %.1f K，%d 个能量点；log10(G/G0)_like=%s"
            % (temperature, points, self._fmt(result.get("log10_G_over_G0_like"))),
            "本地阶段：只做 Landauer 积分，无集群提交",
        ]
        spin_block = result.get("spin") if isinstance(result, dict) else None
        if isinstance(spin_block, dict):
            notes.append(
                "自旋通道：%d 列；E_F 处 P=(Tα−Tβ)/(Tα+Tβ)=%s，"
                "积分 P(L0)=%s（G、S 仍按双通道均值）"
                % (spin_block.get("channels"),
                   self._fmt(spin_block.get("spin_polarization_at_ef")),
                   self._fmt(spin_block.get("spin_polarization_L0"))))
        artifacts = ["results.json"]
        try:
            meta = plotting.plot_te(
                te, out=out_dir / "T_E.png", scale="log", figsize="single",
                formats=("png", "pdf"), annotate=True, result=result,
                temperature=temperature)
        except Exception as exc:  # matplotlib 缺失/渲染失败时不让阶段失败
            notes.append("统一版式出图跳过: %s（pip install -e .[plot]）" % exc)
        else:
            artifacts += [Path(p).name for p in meta.get("paths") or []]
            notes.append("统一版式图（style %s）: %s"
                         % (meta.get("style_version"), " + ".join(
                             Path(p).name for p in meta.get("paths") or [])))
        if points <= MIN_PRODUCTION_POINTS:
            notes.append("7 点样本免责: 仅验证管线, 非生产科学结论(_like 语义)")
        return self.prepared(
            ctx, stage_dir,
            inputs=[te], artifacts=artifacts,
            resources={"local": True, "points": points,
                       "temperature_K": result.get("temperature_K")},
            notes=notes,
            commands=[],
        )

    # -- collect ----------------------------------------------------------- #
    def collect(self, ctx, stage_dir: Path) -> Collected:
        stage_dir = Path(stage_dir)
        out_dir = ctx.override_path("out") or stage_dir
        results_path = Path(out_dir) / "results.json"
        if not results_path.is_file():
            raise StageError(
                "junction_analyze 未产出 results.json（本地阶段应先执行 prepare）: %s"
                % results_path)
        try:
            data = json.loads(results_path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise StageError("junction_analyze 无法解析 results.json: %s" % exc) from exc
        if not isinstance(data, dict):
            raise StageError("junction_analyze 的 results.json 不是 JSON 对象")
        notes: list[str] = []
        if int(data.get("points") or 0) <= MIN_PRODUCTION_POINTS:
            notes.append("7 点样本免责: 仅验证管线, 非生产科学结论(_like 语义)")
        return self.collected(
            ctx, stage_dir, artifacts=["results.json"],
            results={"junction": data}, notes=notes)

    # -- 内部 -------------------------------------------------------------- #
    def _te_path(self, ctx) -> Path:
        override = ctx.override_path("te")
        if override is not None and override.is_dir():
            # --te 既接受 TE.dat 本身，也接受 aitranss_te 输出目录。
            override = override / "TE.dat"
        return ctx.upstream(self.name, "TE.dat", override=override)

    @staticmethod
    def _fmt(value) -> str:
        try:
            return "%.6g" % float(value)
        except (TypeError, ValueError):
            return str(value)
