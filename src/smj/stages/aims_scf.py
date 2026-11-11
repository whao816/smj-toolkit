"""aims_scf：Au60 结的 FHI-aims 自洽场计算（远端阶段）。

prepare 把 ``au_relax/geometry.in``（固定背层弛豫后的几何）拷进本阶段目录
（作业脚本里用相对路径），写 ``control.in`` + ``sub_aims.sh``；AIMS 跑完产出
``aims.out``、``basis-indices.out`` 与 ``mos.aims`` —— 后两者是 aitranss_te 的
必要输入。``electrode.relax=false`` 时 au_relax 为本地直通，其 geometry.in
就是 au60_build 的原几何，这里无需分支。
"""

from __future__ import annotations

import shutil
from pathlib import Path

from .._vendor.junction.inputs import make_aims_scf_input
from ..errors import StageBlocked
from ._helpers import est_record, submit_command, vdw_mode
from .base import Collected, Prepared, Stage, register


@register
class AimsScfStage(Stage):
    name = "aims_scf"
    depends = ("au_relax",)
    script = "sub_aims.sh"
    input_hint = "消费 au_relax/geometry.in（弛豫后的结几何）"
    outputs = ("aims.out",)

    # -- prepare ----------------------------------------------------------- #
    def prepare(self, ctx, stage_dir: Path) -> Prepared:
        stage_dir = Path(stage_dir)
        upstream = ctx.upstream(
            self.name, "geometry.in", override=ctx.override_path("geometry"))
        species = self._species_dir(ctx)
        geometry = stage_dir / "geometry.in"
        shutil.copy2(upstream, geometry)
        natoms = self.count_aims_atoms(geometry)
        estimate = ctx.estimate(self.name, natoms)
        env = ctx.stage_env(self.name, estimate)
        written = make_aims_scf_input(
            geometry, stage_dir, species_dir=species, env=env,
            vdw=vdw_mode(ctx.config.electrode.vdw),
        )
        ctx.log("aims_scf: 已写 control.in（%d 原子）" % natoms)
        return self.prepared(
            ctx, stage_dir,
            inputs=[geometry, Path(written["control"]), Path(written["submit_script"]),
                    stage_dir / "input_summary.json"],
            resources=est_record(estimate),
            notes=[
                "%d 原子，ntasks=%s，partition=%s"
                % (natoms, env.get("ntasks"), env.get("partition")),
                "prepare 只写输入与作业脚本，未提交任何远端命令",
            ],
            commands=submit_command(self.script),
        )

    # -- collect ----------------------------------------------------------- #
    def collect(self, ctx, stage_dir: Path) -> Collected:
        stage_dir = Path(stage_dir)
        output = ctx.require_stage_file(self.name, "aims.out")
        try:
            text = output.read_text(encoding="utf-8", errors="replace")
        except OSError:
            text = ""
        notes: list[str] = []
        has_basis = (stage_dir / "basis-indices.out").is_file()
        has_mos = (stage_dir / "mos.aims").is_file()
        if not has_basis:
            notes.append("缺 basis-indices.out：aitranss_te 的 $basis 输入不完整")
        if not has_mos:
            notes.append("缺 mos.aims：aitranss_te 无法自动确定 $nsaos")
        results = {
            "aims_scf": {
                "atom_count": self.count_aims_atoms(stage_dir / "geometry.in"),
                "normal_termination": "Have a nice day" in text,
                "has_basis_indices": has_basis,
                "has_mos": has_mos,
            }
        }
        return self.collected(
            ctx, stage_dir,
            artifacts=[name for name in ("aims.out", "control.in", "geometry.in",
                                         "basis-indices.out", "mos.aims")
                       if (stage_dir / name).is_file()],
            results=results, notes=notes,
        )

    # -- 内部 -------------------------------------------------------------- #
    def _species_dir(self, ctx) -> Path:
        species = ctx.species_dir()
        if species is None:
            raise StageBlocked(
                self.name,
                "缺少 FHI-aims species 目录：请在 smj.yaml 设 cluster.species_dir 或 paths.species_defaults")
        if not species.is_dir():
            raise StageBlocked(self.name, "FHI-aims species 目录不存在: %s" % species, [str(species)])
        return species
