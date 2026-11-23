"""thermal_relax（弛豫声子器件）：FHI-aims bfgs 弛豫周期板（远端阶段）。

与电子链的 ``aims_tip_relax`` 有三处不同,都是为了让力常数的有限位移有意义:

1. **带 k 网格**:``control_thermal_relax.in`` 写 ``k_grid``(二维板,kz=1)。
   团簇不需要 k 点,但周期板没有 k 点等于只算 Γ,金属的电子密度不收敛。
2. **力必须被打印**:``compute_forces .true.`` + ``final_forces_cleaned .true.``。
   电子链的 relax 模板不带这两行(它只用 ``sc_accuracy_forces`` 判收敛),
   所以 ``aims.out`` 里根本没有力块 —— 报告 §8 的 ``calculate_forces`` 也是错的,
   手册里的关键字是 ``compute_forces``(附录 A-3)。
3. **冻结最外层**:``constrain_relaxation`` 已在 ``thermal_device`` 写入
   geometry.in,弛豫时最外层不动,层间距不被整体拉伸改变。
"""

from __future__ import annotations

import shutil
from pathlib import Path

from ..errors import StageBlocked, StageError
from ..junction.tcontrol import make_control
from ..phonon.forces import ForceBlockError, read_forces
from ..provenance import sha256_file
from ._aims_log import describe, relax_status
from ._helpers import est_record, submit_command
from ._thermal_common import (
    THERMAL_OCCUPATION_WIDTH_EV,
    render_thermal_sub,
    upstream_file,
    write_script,
)
from .base import Collected, Prepared, Stage, register


@register
class ThermalRelaxStage(Stage):
    name = "thermal_relax"
    depends = ("thermal_device",)
    script = "sub_thermal_relax.sh"
    input_hint = "消费 thermal_device/geometry.in（周期板器件）"
    outputs = ("aims.out", "geometry.in.next_step")

    def prepare(self, ctx, stage_dir: Path) -> Prepared:
        stage_dir = Path(stage_dir)
        thermal = ctx.config.thermal
        source = upstream_file(
            ctx, self.name, "geometry.in", "thermal_device",
            override=ctx.override_path("geometry"))
        geometry = stage_dir / "geometry.in"
        shutil.copy2(source, geometry)
        species = ctx.species_dir()
        if species is None:
            raise StageBlocked(
                self.name,
                "缺少 FHI-aims species 目录：请设 cluster.species_dir 或 paths.species_defaults")
        fc = thermal.fc
        electrode = thermal.electrode
        control = make_control(
            geometry, stage_dir / "control.in", "thermal_relax", species,
            xc="pbe",
            occupation_width=THERMAL_OCCUPATION_WIDTH_EV,
            force_threshold=str(fc.relax_force_threshold),
            k_grid=tuple(fc.k_grid),
            sc_accuracy_forces=str(fc.sc_accuracy_forces),
            # 金属板 + ~20 Å 真空 + Γ 是 charge sloshing 的典型场合:300 K 下第一版
            # (sc_accuracy_rho 1E-5 / charge_mix_param 0.2 / n_max_pulay 10) 跑满
            # 175 圈还卡在 4.5e-3。放松 ρ 到 1E-4(力只需要这个量级)、加大混参与
            # Pulay 历史。见 thermal-conductance-devplan §2.4。
            sc_accuracy_rho="1E-4",
            charge_mix_param=0.35,
            n_max_pulay=20,
            vdw="hirshfeld" if getattr(ctx.config.electrode, "vdw", True) else None,
        )
        natoms = self.count_aims_atoms(geometry)
        estimate = ctx.estimate(self.name, natoms)
        env = ctx.stage_env(self.name, estimate)
        script = write_script(
            stage_dir / self.script,
            render_thermal_sub(self.script, env,
                               ntasks=int(env.get("ntasks") or 8)))
        notes = [
            "%d 原子周期板，k_grid=%d×%d×%d，sc_accuracy_forces=%s"
            % (natoms, fc.k_grid[0], fc.k_grid[1], fc.k_grid[2],
               fc.sc_accuracy_forces),
            "弛豫阈值 %s eV/Å；最外层 %d 层由 constrain_relaxation 冻结"
            % (fc.relax_force_threshold, electrode.freeze_outer_layers),
            "control.in 含 compute_forces 与 final_forces_cleaned：aims.out 必须有力块",
            "prepare 只写输入与作业脚本，未提交任何远端命令",
        ]
        ctx.log("thermal_relax: 已写 control.in（%d 原子）" % natoms)
        return self.prepared(
            ctx, stage_dir,
            inputs=[geometry, control, script],
            resources=est_record(estimate),
            notes=notes,
            commands=submit_command(self.script),
            parent_stage="thermal_device",
            parent_fingerprint=sha256_file(source),
        )

    def collect(self, ctx, stage_dir: Path) -> Collected:
        stage_dir = Path(stage_dir)
        output = ctx.require_stage_file(self.name, "aims.out")
        text = output.read_text(encoding="utf-8", errors="replace")
        if "Have a nice day" not in text:
            raise StageError(
                "thermal_relax 的 FHI-aims 未正常结束（aims.out 无 'Have a nice day'）；"
                "拒绝 collect 残缺的弛豫结果")
        status = relax_status(output)
        if status is not None and not status["converged"]:
            raise StageError(
                "thermal_relax 未收敛：%s（力常数必须来自收敛构型，否则虚频与声学支"
                "都会被残余应力污染）" % describe(status))
        try:
            forces = read_forces(text)
        except ForceBlockError as exc:
            raise StageError("thermal_relax 解析力块失败: %s" % exc) from exc
        max_force = float(abs(forces).max())
        next_step = stage_dir / "geometry.in.next_step"
        if not next_step.is_file():
            raise StageError(
                "thermal_relax 缺 geometry.in.next_step：FHI-aims 没有产出弛豫后的几何，"
                "无法进入力常数扫描")
        results = {
            "thermal_relax": {
                "atom_count": self.count_aims_atoms(next_step),
                "normal_termination": True,
                "relax_converged": True,
                "max_force_ev_per_A": max_force,
                "force_threshold_ev_per_A": float(
                    ctx.config.thermal.fc.relax_force_threshold),
                "n_force_blocks": text.count("Total atomic forces"),
            }
        }
        return self.collected(
            ctx, stage_dir,
            artifacts=[name for name in self.outputs
                       if (stage_dir / name).is_file()],
            results=results,
            notes=["弛豫后最大力 %.3e eV/Å（%d 个力块，取最后一块）"
                   % (max_force, results["thermal_relax"]["n_force_blocks"])])
