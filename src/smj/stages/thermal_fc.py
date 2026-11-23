"""thermal_fc（有限位移力常数扫描）：一个 Slurm **数组作业** 算完 6·N_dev 个位移点。

一个器件区原子要 6 个位移点(±x/±y/±z);器件区 N_dev=54 时共 324 个。
中心差分 ``K_ij = -(F_i(+δ) - F_i(-δ)) / 2δ`` 需要成对的正负位移,所以位移数
恒为 6·N_dev,**不是** 3·N_dev(那是只做前向差分的下界,精度差一个量级)。

为什么是数组而不是 324 个独立作业:一次 sbatch、一个 job id、排一次队;
collect 只要数 ``fc/*/forces.json`` 的份数就知道完整不完整。

**每个数组任务的工作**(见 ``sub_thermal_fc.sh``):
``fc/NNNN/`` 下建 geometry.in(施加该位移)、复制 control.in、跑 FHI-aims、
把最后一个力块抽成 ``forces.json``。三个脚本都是 stdlib-only,跟着作业传过去。
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np

from ..errors import StageBlocked, StageError
from ..junction.tcontrol import make_control
from ..phonon.fc import assemble_fc, make_displacements
from ..provenance import sha256_file
from ._helpers import est_record, submit_command
from ._thermal_common import render_thermal_sub, upstream_file, write_script
from .base import Collected, Prepared, Stage, register

_HERE = Path(__file__).resolve().parent.parent / "phonon"
_HELPERS = ("apply_displacement.py", "extract_forces.py")


@register
class ThermalFcStage(Stage):
    name = "thermal_fc"
    depends = ("thermal_relax",)
    script = "sub_thermal_fc.sh"
    input_hint = "消费 thermal_relax/geometry.in.next_step（弛豫后的器件）"
    outputs = ("displacements.json", "force_constants.npy", "fc_summary.json")

    def prepare(self, ctx, stage_dir: Path) -> Prepared:
        stage_dir = Path(stage_dir)
        thermal = ctx.config.thermal
        fc_cfg = thermal.fc
        source = upstream_file(
            ctx, self.name, "geometry.in.next_step", "thermal_relax",
            override=ctx.override_path("geometry"))
        geometry = stage_dir / "geometry.in"
        shutil.copy2(source, geometry)
        natoms = self.count_aims_atoms(geometry)

        device_indices = self._device_indices(ctx)
        if not device_indices:
            raise StageBlocked(
                self.name,
                "缺 thermal_device/device_metadata.json 的 device_indices："
                "无法确定哪些原子要位移（放宽到全部原子会让位移数翻倍且包含冻结 lead）",
                [str(ctx.artifact("thermal_device", "device_metadata.json"))])
        disps = make_displacements(natoms, device_indices,
                                   float(fc_cfg.displacement))
        manifest = stage_dir / "displacements.json"
        manifest.write_text(
            json.dumps(disps, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8")

        species = ctx.species_dir()
        if species is None:
            raise StageBlocked(
                self.name,
                "缺少 FHI-aims species 目录：请设 cluster.species_dir 或 paths.species_defaults")
        control = make_control(
            geometry, stage_dir / "control.in", "thermal_force", species,
            xc="pbe",
            k_grid=tuple(fc_cfg.k_grid),
            sc_accuracy_forces=str(fc_cfg.sc_accuracy_forces),
            # 与 thermal_relax 同一套金属板混参(见该阶段注释):324 个位移作业
            # 每个都是「金属板 + 真空 + Γ」,用 1E-5/0.2/10 会集体卡住。
            sc_accuracy_rho="1E-4",
            charge_mix_param=0.35,
            n_max_pulay=20,
            vdw="hirshfeld" if getattr(ctx.config.electrode, "vdw", True) else None,
        )
        for name in _HELPERS:
            shutil.copy2(_HERE / name, stage_dir / name)
        estimate = ctx.estimate(self.name, natoms)
        env = ctx.stage_env(self.name, estimate)
        script = write_script(
            stage_dir / self.script,
            render_thermal_sub(
                self.script, env, ntasks=int(env.get("ntasks") or 8),
                extra={"last": len(disps) - 1,
                       "parallel": int(fc_cfg.array_parallel),
                       "count": len(disps),
                       "walltime": str(fc_cfg.array_time_limit),
                       "mem": str(fc_cfg.memory)}))
        notes = [
            "%d 原子，器件区 %d 个原子 → %d 个位移点（±%g Å ×3 轴）"
            % (natoms, len(device_indices), len(disps), fc_cfg.displacement),
            "数组作业：--array=0-%d%%%d（并发 %d，每点 %s 内存，墙钟 %s）"
            % (len(disps) - 1, int(fc_cfg.array_parallel),
               int(fc_cfg.array_parallel), fc_cfg.memory, fc_cfg.array_time_limit),
            "每个位移点独立起算，不共享密度（避免把一个位移的密度带给另一个位移）",
            "prepare 只写输入与作业脚本，未提交任何远端命令",
        ]
        ctx.log("thermal_fc: %d 个位移点（器件区 %d 原子）"
                % (len(disps), len(device_indices)))
        return self.prepared(
            ctx, stage_dir,
            inputs=[geometry, control, manifest, script,
                    *[stage_dir / name for name in _HELPERS]],
            resources=est_record(estimate),
            notes=notes,
            commands=submit_command(self.script),
            parent_stage="thermal_relax",
            parent_fingerprint=sha256_file(source),
        )

    def collect(self, ctx, stage_dir: Path) -> Collected:
        stage_dir = Path(stage_dir)
        manifest = stage_dir / "displacements.json"
        if not manifest.is_file():
            raise StageError("thermal_fc 缺 displacements.json")
        disps = json.loads(manifest.read_text(encoding="utf-8"))
        n_atoms = self.count_aims_atoms(stage_dir / "geometry.in")
        forces = self._read_forces(stage_dir, len(disps), n_atoms)
        try:
            k_dev = assemble_fc(disps, forces, n_atoms)
        except ValueError as exc:
            raise StageError("thermal_fc 组装力常数失败: %s" % exc) from exc
        np.save(stage_dir / "force_constants.npy", k_dev)

        elements = self._elements(stage_dir / "geometry.in")
        results = {
            "thermal_fc": {
                "n_atoms": n_atoms,
                "n_displacements": len(disps),
                "displacement_A": float(disps[0]["delta_A"]) if disps else 0.0,
                "n_device_atoms": len({int(d["atom"]) for d in disps}),
                "force_constants_shape": list(k_dev.shape),
                "fc_max_abs": float(np.abs(k_dev).max()) if k_dev.size else 0.0,
                "n_elements": len(set(elements)),
            }
        }
        notes = ["力常数 %d×%d，max|K|=%.4g eV/Å²（ASR 残差在 thermal_analyze 里"
                 "按**完整**矩阵算：这里只有器件行，求和必然不闭合）"
                 % (k_dev.shape[0], k_dev.shape[1],
                    results["thermal_fc"]["fc_max_abs"])]
        return self.collected(
            ctx, stage_dir,
            artifacts=[name for name in ("displacements.json", "force_constants.npy")
                       if (stage_dir / name).is_file()],
            results=results, notes=notes)

    # -- 内部 -------------------------------------------------------------- #
    @staticmethod
    def _device_indices(ctx) -> list[int]:
        path = ctx.artifact("thermal_device", "device_metadata.json")
        try:
            meta = json.loads(Path(path).read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return []
        return [int(i) for i in (meta.get("device_indices") or [])]

    @staticmethod
    def _read_forces(stage_dir: Path, expected: int, n_atoms: int) -> np.ndarray:
        """按 slot 顺序读 ``fc/NNNN/forces.json``;缺任何一个都拒绝。"""
        found: dict[int, np.ndarray] = {}
        missing: list[int] = []
        for slot in range(expected):
            record = stage_dir / "fc" / ("%04d" % slot) / "forces.json"
            if not record.is_file():
                missing.append(slot)
                continue
            try:
                payload = json.loads(record.read_text(encoding="utf-8"))
                array = np.asarray(payload["forces"], dtype=float)
            except (OSError, ValueError, KeyError, TypeError) as exc:
                raise StageError("thermal_fc 读取 %s 失败: %s" % (record, exc)) from exc
            if array.shape != (n_atoms, 3):
                raise StageError(
                    "thermal_fc 的 %s 力数组形状是 %r，应为 (%d, 3)"
                    % (record, array.shape, n_atoms))
            found[slot] = array
        if missing:
            raise StageError(
                "thermal_fc 缺 %d 个位移点的力（前几个 slot: %s）；"
                "数组作业没跑完就 collect 会让力常数矩阵少行"
                % (len(missing), missing[:5]))
        return np.stack([found[slot] for slot in range(expected)], axis=0)

    @staticmethod
    def _elements(geometry: Path) -> list[str]:
        elements = []
        for line in Path(geometry).read_text(encoding="utf-8").splitlines():
            parts = line.split()
            if len(parts) >= 5 and parts[0] == "atom":
                elements.append(parts[4])
        return elements
