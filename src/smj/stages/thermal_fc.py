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
from ._thermal_common import (
    THERMAL_OCCUPATION_WIDTH_EV,
    render_thermal_sub,
    upstream_file,
    write_script,
)
from .base import Collected, Prepared, Stage, register

_HERE = Path(__file__).resolve().parent.parent / "phonon"
_HELPERS = ("apply_displacement.py", "extract_forces.py", "concat_forces.py")


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
            occupation_width=THERMAL_OCCUPATION_WIDTH_EV,
            k_grid=tuple(fc_cfg.k_grid),
            sc_accuracy_forces=str(fc_cfg.sc_accuracy_forces),
            # 与 thermal_relax 同一套金属板混参(见该阶段注释):324 个位移作业
            # 每个都是「金属板 + 真空 + Γ」,用 1E-5/0.2/10 会集体卡住。
            sc_accuracy_rho="1E-4",
            charge_mix_param=0.35,
            n_max_pulay=20,
            vdw="hirshfeld" if getattr(ctx.config.electrode, "vdw", True) else None,
        )
        estimate = ctx.estimate(self.name, natoms)
        script = self.write_job_files(ctx, stage_dir, len(disps), estimate)
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

    @classmethod
    def write_job_files(cls, ctx, stage_dir: Path, count: int, estimate: dict) -> Path:
        """刷新作业脚本及 Python 工具，直接复用已准备的计算输入。"""
        fc_cfg = ctx.config.thermal.fc
        for name in _HELPERS:
            shutil.copy2(_HERE / name, stage_dir / name)
        env = ctx.stage_env(cls.name, estimate)
        return write_script(stage_dir / cls.script, render_thermal_sub(
            cls.script, env, ntasks=int(env.get("ntasks") or 8),
            extra={"last": count - 1, "count": count, "parallel": int(fc_cfg.array_parallel),
                   "walltime": str(fc_cfg.array_time_limit), "mem": str(fc_cfg.memory)}))

    def collect(self, ctx, stage_dir: Path) -> Collected:
        stage_dir = Path(stage_dir)
        manifest = stage_dir / "displacements.json"
        if not manifest.is_file():
            raise StageError("thermal_fc 缺 displacements.json")
        disps = json.loads(manifest.read_text(encoding="utf-8"))
        n_atoms = self.count_aims_atoms(stage_dir / "geometry.in")
        forces, summary = self._read_force_data(stage_dir, len(disps), n_atoms)
        try:
            k_dev = assemble_fc(disps, forces, n_atoms)
        except ValueError as exc:
            raise StageError("thermal_fc 组装力常数失败: %s" % exc) from exc
        np.save(stage_dir / "force_constants.npy", k_dev)

        elements = self._elements(stage_dir / "geometry.in")
        cycles = [s["n_cycles"] for s in (summary or {}).get("slots", [])
                  if type(s.get("n_cycles")) is int and s["n_cycles"] > 0]
        results = {
            "thermal_fc": {
                "n_atoms": n_atoms,
                "n_displacements": len(disps),
                "displacement_A": float(disps[0]["delta_A"]) if disps else 0.0,
                "n_device_atoms": len({int(d["atom"]) for d in disps}),
                "force_constants_shape": list(k_dev.shape),
                "fc_max_abs": float(np.abs(k_dev).max()) if k_dev.size else 0.0,
                "n_elements": len(set(elements)),
                # 收敛证据:没有 fc_summary.json 时是 None(老 run 兼容),
                # 有则是远端逐槽位读 aims.out 的结果。
                "scf_converged": True if summary is not None else None,
                "scf_cycles_min": min(cycles) if cycles else None,
                "scf_cycles_max": max(cycles) if cycles else None,
            }
        }
        notes = ["力常数 %d×%d，max|K|=%.4g eV/Å²（ASR 残差在 thermal_analyze 里"
                 "按**完整**矩阵算：这里只有器件行，求和必然不闭合）"
                 % (k_dev.shape[0], k_dev.shape[1],
                    results["thermal_fc"]["fc_max_abs"])]
        if cycles:
            notes.append("SCF：%d 个位移点收敛，其中 %d 个有圈数记录（%d–%d）"
                         % (len(disps), len(cycles), min(cycles), max(cycles)))
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
    def _read_json(path: Path) -> dict:
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise StageError("thermal_fc 读取 %s 失败: %s" % (path, exc)) from exc
        if not isinstance(data, dict):
            raise StageError("thermal_fc 的 %s 格式错误，应为 JSON 对象" % path)
        return data

    @staticmethod
    def _retry_hint(slots) -> str:
        return "smj thermal --stage thermal_fc --retry " + ",".join(map(str, sorted(slots)))

    @classmethod
    def _check_scf(cls, summary: dict, expected: int) -> None:
        if (not isinstance(summary, dict) or type(summary.get("count")) is not int
                or summary["count"] != expected or not isinstance(summary.get("slots"), list)):
            raise StageError("fc_summary.json 的槽位数或 slots 格式有误，请重新回取汇总")
        retry = set()
        for key in ("missing", "not_converged", "unknown"):
            indices = summary.get(key, [])
            if not isinstance(indices, list) or any(
                    type(s) is not int or not 0 <= s < expected for s in indices):
                raise StageError("fc_summary.json 的 %s 格式有误，请重新回取汇总" % key)
            retry.update(indices)
        seen = set()
        for item in summary["slots"]:
            slot = item.get("slot") if isinstance(item, dict) else None
            if type(slot) is not int or not 0 <= slot < expected or slot in seen:
                raise StageError("fc_summary.json 含重复或无效槽位，请重新回取汇总")
            seen.add(slot)
            if item.get("converged") is not True:
                retry.add(slot)
        retry.update(set(range(expected)) - seen)
        if retry:
            raise StageError(
                "thermal_fc 槽位缺力、SCF 未收敛或缺少收敛证据。用 `%s` 重跑后再 collect"
                % cls._retry_hint(retry))

    @classmethod
    def _read_forces(cls, stage_dir: Path, expected: int, n_atoms: int) -> np.ndarray:
        forces, _ = cls._read_force_data(stage_dir, expected, n_atoms)
        return forces

    @classmethod
    def _read_force_data(cls, stage_dir: Path, expected: int,
                         n_atoms: int) -> tuple[np.ndarray, dict | None]:
        """远端重试只回取汇总；没有汇总的旧目录可从逐槽位文件读取。"""
        summary_path = stage_dir / "fc_summary.json"
        summary = cls._read_json(summary_path) if summary_path.is_file() else None
        bundle = stage_dir / "forces_all.json"
        if bundle.is_file():
            payload = cls._read_json(bundle)
            if "summary" in payload:
                embedded = payload["summary"]
                if not isinstance(embedded, dict) or (summary is not None and summary != embedded):
                    raise StageError("forces_all.json 与 fc_summary.json 的收敛证据不符，请重新回取")
                summary = embedded
            blocks = payload.get("forces")
            if not isinstance(blocks, list) or len(blocks) != expected:
                raise StageError("thermal_fc 的 forces_all.json 位移点数应为 %d" % expected)
        else:
            blocks = []
            for slot in range(expected):
                record = stage_dir / "fc" / ("%04d" % slot) / "forces.json"
                blocks.append(cls._read_json(record).get("forces") if record.is_file() else None)
        if summary is not None:
            cls._check_scf(summary, expected)

        missing = [i for i, block in enumerate(blocks) if block is None]
        if missing:
            raise StageError("thermal_fc 缺 %d 个位移点的力。用 `%s` 重跑后再 collect"
                             % (len(missing), cls._retry_hint(missing)))
        forces = np.empty((expected, n_atoms, 3))
        for slot, block in enumerate(blocks):
            try:
                array = np.asarray(block, dtype=float)
            except (ValueError, TypeError) as exc:
                raise StageError("thermal_fc 槽位 %d 力格式错误: %s" % (slot, exc)) from exc
            if array.shape != (n_atoms, 3) or not np.isfinite(array).all():
                raise StageError("thermal_fc 槽位 %d 力形状或数值无效。用 `%s` 重跑"
                                 % (slot, cls._retry_hint([slot])))
            forces[slot] = array
        return forces, summary

    @staticmethod
    def _elements(geometry: Path) -> list[str]:
        elements = []
        for line in Path(geometry).read_text(encoding="utf-8").splitlines():
            parts = line.split()
            if len(parts) >= 5 and parts[0] == "atom":
                elements.append(parts[4])
        return elements
