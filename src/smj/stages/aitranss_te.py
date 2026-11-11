"""aitranss_te：自洽场结果 → aitranss 透射谱（远端阶段）。

prepare 写 ``tcontrol``、``self.energy.in`` 与 ``sub_te.sh``：
- ``tcontrol`` 由 ``aims_scf`` 目录里的 ``mos.aims``/``geometry.in`` 与
  ``au60_build/electrode_metadata.json`` 的表面索引生成；
- ``self.energy.in`` 按几何里实际 Au 数分层生成（左右以 metadata 的 tip 索引为准）。
aitranss 跑完产出 ``TE.dat``（collect 只做轻量核对，真正的量在 junction_analyze）。
"""

from __future__ import annotations

import json
from pathlib import Path

from .._vendor.junction.geometry import read_aims_geometry
from .._vendor.junction.inputs import make_aitranss_input
from .._vendor.junction.landauer import read_te
from .._vendor.junction.self_energy import make_self_energy
from ..errors import StageBlocked, StageError
from ._helpers import est_record, submit_command
from .base import Collected, Prepared, Stage, register


@register
class AitranssTeStage(Stage):
    name = "aitranss_te"
    depends = ("aims_scf",)
    script = "sub_te.sh"
    input_hint = "消费 aims_scf/{aims.out,mos.aims,geometry.in}（--scf-dir 可显式覆盖）"
    outputs = ("TE.dat",)

    # -- prepare ----------------------------------------------------------- #
    def prepare(self, ctx, stage_dir: Path) -> Prepared:
        stage_dir = Path(stage_dir)
        scf_dir = self._scf_dir(ctx)
        scf_geometry = scf_dir / "geometry.in"
        if not scf_geometry.is_file():
            raise StageBlocked(
                self.name, "缺 aims_scf/geometry.in（tcontrol 需要结几何）", [str(scf_geometry)])
        metadata = self._metadata(ctx)
        geometry = ctx.artifact("au_relax", "geometry.in")
        if not geometry.is_file():
            geometry = ctx.artifact("au60_build", "geometry.in")
        if not geometry.is_file():
            # 独立运行（--scf-dir）时上游目录可能不在：SCF 目录里的几何是同一份。
            geometry = scf_geometry
        natoms = self.count_aims_atoms(scf_geometry)
        estimate = ctx.estimate(self.name, natoms)
        env = ctx.stage_env(self.name, estimate)
        ener, estep, eend = ctx.config.transport.resolved_window()
        nlayers = self._resolve_nlayers(ctx.config.electrode.nlayers, metadata)
        try:
            written = make_aitranss_input(
                scf_dir, stage_dir, metadata=metadata, nlayers=nlayers,
                ener=ener, estep=estep, eend=eend, env=env,
            )
        except (FileNotFoundError, KeyError, ValueError) as exc:
            raise StageBlocked(
                self.name, "aitranss_te 输入不完整（mos.aims / efermi / 表面索引）: %s" % exc,
                [str(scf_dir)]) from exc
        values = written.get("values") or {}
        # self.energy.in：左右电极划分以 metadata 的 tip 索引为准（不得按距离猜）。
        au_count = sum(1 for atom in read_aims_geometry(geometry) if atom.element == "Au")
        if au_count < 4 or au_count % 2:
            raise StageBlocked(
                self.name, "结内 Au 数异常(需两侧对称且 >=4): %d" % au_count, [str(geometry)])
        try:
            self_energy = make_self_energy(
                geometry, stage_dir / "self.energy.in",
                au_per_electrode=au_count // 2, metadata=metadata,
            )
        except ValueError as exc:
            raise StageError("aitranss_te 无法生成 self.energy.in: %s" % exc) from exc
        ctx.log("aitranss_te: 已写 tcontrol + self.energy.in（%d 原子，%d Au）"
                % (natoms, au_count))
        return self.prepared(
            ctx, stage_dir,
            inputs=[Path(written["tcontrol"]), stage_dir / "self.energy.in",
                    Path(written["submit_script"]), stage_dir / "input_summary.json"],
            resources=est_record(estimate),
            notes=[
                "能量窗 ener=%s / estep=%s / eend=%s，nlayers=%d"
                % (values.get("ener", ener), values.get("estep", estep),
                   values.get("eend", eend), values.get("nlayers", nlayers)),
                "表面索引 l=%s r=%s（来自 %s）"
                % ([values.get("lsurc"), values.get("lsurx"), values.get("lsury")],
                   [values.get("rsurc"), values.get("rsurx"), values.get("rsury")],
                   Path(metadata).name),
                "self.energy.in：每侧 %d 个电极 Au（含 tip，tip 记 empty）"
                % (au_count // 2),
                "prepare 只写输入与作业脚本，未提交任何远端命令",
            ],
            commands=submit_command(self.script),
        )

    # -- collect ----------------------------------------------------------- #
    def collect(self, ctx, stage_dir: Path) -> Collected:
        stage_dir = Path(stage_dir)
        te = ctx.require_stage_file(self.name, "TE.dat")
        try:
            rows = read_te(te)
        except ValueError as exc:
            raise StageError("aitranss_te 无法解析 TE.dat: %s" % exc) from exc
        notes: list[str] = []
        if not (stage_dir / "te.out").is_file():
            notes.append("缺 te.out（aitranss 主输出）：请确认作业日志已回传")
        if not (stage_dir / "self.energy.in").is_file():
            notes.append("缺 self.energy.in：结果不可复现")
        results = {
            "aitranss_te": {
                "points": len(rows),
                "te_bytes": te.stat().st_size,
                "max_transmission": max(row[2] for row in rows),
                "has_te_out": (stage_dir / "te.out").is_file(),
            }
        }
        return self.collected(
            ctx, stage_dir, artifacts=["TE.dat"], results=results, notes=notes)

    # -- 内部 -------------------------------------------------------------- #
    def _scf_dir(self, ctx) -> Path:
        override = ctx.override_path("scf_dir")
        if override is not None:
            if not override.is_dir():
                raise StageBlocked(
                    self.name, "显式 scf 目录不存在: %s" % override, [str(override)])
            return override
        return ctx.upstream(self.name, "aims.out").parent

    def _metadata(self, ctx) -> Path:
        override = ctx.override_path("metadata")
        if override is not None:
            if not override.is_file():
                raise StageBlocked(
                    self.name, "显式 metadata 不存在: %s" % override, [str(override)])
            return override
        # 0.1.2：au_relax 会重写 metadata（指纹 + relax 记录）；存在时以它为准，
        # 否则回退 au60_build（electrode.relax=false 时两者等价）。
        for stage in ("au_relax", "au60_build"):
            candidate = ctx.artifact(stage, "electrode_metadata.json")
            if candidate.is_file():
                return candidate
        raise StageBlocked(
            self.name,
            "缺 au_relax（或 au60_build）/electrode_metadata.json（tcontrol 需要六个表面索引，"
            "self.energy.in 需要左右 tip 索引）",
            [str(ctx.artifact("au60_build", "electrode_metadata.json"))],
        )

    @staticmethod
    def _resolve_nlayers(configured: int | None, metadata: Path) -> int:
        """``electrode.nlayers`` 为 null 时用库文件建议值（手册 §5.3.2），再回退 3。"""
        if configured is not None:
            return int(configured)
        try:
            data = json.loads(Path(metadata).read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return 3
        suggested = (data.get("library_cluster") or {}).get("nlayers_suggested")
        try:
            return int(suggested) if suggested else 3
        except (TypeError, ValueError):
            return 3
