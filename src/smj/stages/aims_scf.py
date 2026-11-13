"""aims_scf：Au60 结的 FHI-aims 自洽场计算（远端阶段）。

prepare 把 ``electrode_build/geometry.in``（贴装后的结几何；``au_relax`` 已从主链路
移除）拷进本阶段目录，写 ``control.in`` + ``sub_aims.sh``；AIMS 跑完产出
``aims.out``、``basis-indices.out`` 与 ``mos.aims`` —— 后两者是 aitranss_te 的
必要输入。结几何在 SCF/输运阶段保持冻结，不再弛豫。
"""

from __future__ import annotations

import shutil
from pathlib import Path

from .._vendor.junction.inputs import make_aims_scf_input
from ..errors import StageBlocked, StageError
from ..provenance import (
    parse_aims_versions,
    record_observed_software,
    sha256_file,
    verify_parent_geometry,
)
from ._helpers import est_record, submit_command, vdw_mode
from .base import Collected, Prepared, Stage, register


@register
class AimsScfStage(Stage):
    name = "aims_scf"
    depends = ("electrode_build",)
    script = "sub_aims.sh"
    input_hint = "消费 electrode_build/geometry.in（冻结的结几何）"
    # aitranss 的 tcontrol 直接消费后三个文件；缺任何一个都跑不了透射。
    outputs = ("aims.out", "basis-indices.out", "omat.aims", "mos.aims")

    # -- prepare ----------------------------------------------------------- #
    def prepare(self, ctx, stage_dir: Path) -> Prepared:
        stage_dir = Path(stage_dir)
        upstream = ctx.upstream(
            self.name, "geometry.in", override=ctx.override_path("geometry"))
        species = self._species_dir(ctx)
        geometry = stage_dir / "geometry.in"
        shutil.copy2(upstream, geometry)
        notes: list[str] = []
        # P0-3：SCF 必须用 electrode_build 产出并已由 A2 确认/锚定的那份几何。
        metadata_path = ctx.artifact("electrode_build", "electrode_metadata.json")
        parent_stage = "electrode_build"
        if metadata_path.is_file():
            verify_parent_geometry(self.name, upstream, metadata_path)
        else:
            notes.append("未找到 electrode_metadata.json：跳过父几何指纹校验（独立运行模式）")
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
                *notes,
                "prepare 只写输入与作业脚本，未提交任何远端命令",
            ],
            commands=submit_command(self.script),
            parent_stage=parent_stage,
            parent_fingerprint=sha256_file(upstream),
        )

    # -- collect ----------------------------------------------------------- #
    def collect(self, ctx, stage_dir: Path) -> Collected:
        stage_dir = Path(stage_dir)
        output = ctx.require_stage_file(self.name, "aims.out")
        try:
            text = output.read_text(encoding="utf-8", errors="replace")
        except OSError as exc:
            raise StageError("aims_scf 无法读取 aims.out: %s" % exc) from exc
        if not text.strip():
            raise StageError("aims_scf/aims.out 为空：FHI-aims 没有产生可用输出")
        if "Have a nice day" not in text:
            raise StageError(
                "FHI-aims 未正常结束（aims.out 无 'Have a nice day'）；"
                "拒绝 collect 残缺的 SCF 结果")

        # tcontrol 的 $basis/$read_omat/$scfmo 直接读这三个文件；缺任何一个，
        # 下游 aitranss_te 都会在 prepare 才炸——这里必须提前 fail-closed。
        required = ("basis-indices.out", "omat.aims", "mos.aims")
        missing = [name for name in required if not (stage_dir / name).is_file()]
        if missing:
            raise StageError(
                "FHI-aims SCF 产物不完整（缺 %s）；aitranss_te 的 tcontrol 需要 "
                "basis-indices.out + omat.aims + mos.aims" % "、".join(missing))

        notes: list[str] = []
        software = parse_aims_versions(text)
        if software:
            record_observed_software(ctx.layout.run_dir, self.name, software)
        else:
            notes.append("未从 aims.out 解析到 FHI-aims 版本/uuid：远端软件版本未入档")
        results = {
            "aims_scf": {
                "atom_count": self.count_aims_atoms(stage_dir / "geometry.in"),
                "normal_termination": True,
                "has_basis_indices": True,
                "has_omat": True,
                "has_mos": True,
                "software": software or None,
            }
        }
        return self.collected(
            ctx, stage_dir,
            artifacts=[name for name in ("aims.out", "control.in", "geometry.in",
                                         "basis-indices.out", "omat.aims", "mos.aims")
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
