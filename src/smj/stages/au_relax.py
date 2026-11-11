"""au_relax：固定两侧背层的几何弛豫（远端 AIMS 阶段，0.1.2 新增）。

手册 §5.3.1 要求：分子经锚定基接到预定义 Au 团簇后，**固定两个背层**做一次
几何弛豫，防止电极在结构优化时塌陷；弛豫后的几何才是 `aims_scf` 的输入。

本阶段：
- prepare 复制 ``au60_build/geometry.in``，给两侧背层整层加
  ``constrain_relaxation .true.``，写 ``control.in``（bfgs 模板）+ ``sub_aims.sh``；
- collect 要求 ``aims.out`` 正常结束且 ``geometry.in.next_step`` 存在，用
  next_step 覆盖 ``geometry.in``；随后按 det-v1 重算表面三元组，必须与
  ``au60_build`` 的 metadata 完全一致（背层冻结后表面原子不该移动）——
  不一致即 StageError，不把"看起来跑完了"的几何交下游；
- ``electrode.relax: false`` 时本阶段退化为本地直通（``local_only``）：
  只把 au60 的几何与 metadata 复制过来，不提交任何作业。
"""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

from .._vendor.junction.electrode import backplane_indices as _backplane_indices
from .._vendor.junction.fingerprint import attach_fingerprint
from .._vendor.junction.geometry import read_aims_geometry
from .._vendor.junction.inputs import make_au_relax_input
from .._vendor.junction.surface_index import verify_surface_indices
from ..errors import StageBlocked, StageError
from ._helpers import est_record, submit_command, vdw_mode
from .base import Collected, Prepared, Stage, register

SURFACE_KEYS = ("left_surface_indices", "right_surface_indices")
BACKPLANE_KEYS = ("left_backplane_indices", "right_backplane_indices")
SKIP_MARKER = "relax_skipped.json"
MOVED_TOL = 1e-4  # 被冻结的原子在弛豫前后应保持原位（Å）


@register
class AuRelaxStage(Stage):
    name = "au_relax"
    depends = ("au60_build",)
    script = "sub_aims.sh"
    input_hint = "消费 au60_build 的 geometry.in 与背层索引"
    outputs = ("aims.out",)

    # -- prepare ----------------------------------------------------------- #
    def prepare(self, ctx, stage_dir: Path) -> Prepared:
        stage_dir = Path(stage_dir)
        upstream = ctx.upstream(
            self.name, "geometry.in", override=ctx.override_path("geometry"))
        geometry = stage_dir / "geometry.in"
        shutil.copy2(upstream, geometry)
        metadata = self._read_metadata(ctx, stage_dir)
        natoms = self.count_aims_atoms(geometry)

        if not bool(getattr(ctx.config.electrode, "relax", True)):
            shutil.copy2(metadata["_path"], stage_dir / "electrode_metadata.json")
            (stage_dir / SKIP_MARKER).write_text(
                json.dumps({"reason": "electrode.relax=false",
                            "note": "跳过背层固定弛豫；下游直接消费 au60_build 几何"},
                           ensure_ascii=False, indent=2),
                encoding="utf-8")
            ctx.log("au_relax: electrode.relax=false，本地直通（%d 原子）" % natoms)
            return self.prepared(
                ctx, stage_dir,
                inputs=[geometry],
                artifacts=["geometry.in", "electrode_metadata.json", SKIP_MARKER],
                resources=est_record({"local": True}),
                notes=["electrode.relax=false：不做 AIMS 弛豫，直接沿用 au60_build 几何",
                       "如需按手册 §5.3.1 固定背层弛豫，把 electrode.relax 改回 true"],
                commands=[],
                local_only=True,
            )

        constrained = self._constrained_atoms(ctx, geometry, metadata)
        species = self._species_dir(ctx)
        estimate = ctx.estimate(self.name, natoms)
        env = ctx.stage_env(self.name, estimate)
        written = make_au_relax_input(
            geometry, stage_dir, constrained=constrained,
            species_dir=species,
            relax_force_threshold=ctx.config.electrode.relax_force_threshold,
            env=env, vdw=vdw_mode(ctx.config.electrode.vdw),
        )
        ctx.log("au_relax: 已写 control.in（%d 原子，冻结 %d 个背层原子）"
                % (natoms, len(constrained)))
        return self.prepared(
            ctx, stage_dir,
            inputs=[geometry, Path(written["control"]), Path(written["submit_script"]),
                    stage_dir / "input_summary.json"],
            resources=est_record(estimate),
            notes=[
                "%d 原子，冻结背层 %d 个（left=%s right=%s）"
                % (natoms, len(constrained),
                   metadata.get("left_backplane_indices"), metadata.get("right_backplane_indices")),
                "force_threshold=%s，ntasks=%s，partition=%s"
                % (ctx.config.electrode.relax_force_threshold,
                   env.get("ntasks"), env.get("partition")),
                "prepare 只写输入与作业脚本，未提交任何远端命令",
            ],
            commands=submit_command(self.script),
        )

    # -- collect ----------------------------------------------------------- #
    def collect(self, ctx, stage_dir: Path) -> Collected:
        stage_dir = Path(stage_dir)
        metadata = self._read_metadata(ctx, stage_dir)

        if (stage_dir / SKIP_MARKER).is_file():
            attach_fingerprint(stage_dir)
            return self.collected(
                ctx, stage_dir,
                artifacts=["geometry.in", "electrode_metadata.json", SKIP_MARKER],
                results={"au_relax": {"enabled": False,
                                      "atom_count": self.count_aims_atoms(stage_dir / "geometry.in")}},
                notes=["electrode.relax=false：本阶段为本地直通，未做弛豫"],
            )

        output = ctx.require_stage_file(self.name, "aims.out")
        try:
            text = output.read_text(encoding="utf-8", errors="replace")
        except OSError:
            text = ""
        if "Have a nice day" not in text:
            raise StageError(
                "au_relax 的 aims.out 没有正常结束（缺 'Have a nice day'）；"
                "几何未收敛，不能交给 aims_scf（先看 %s）" % output)
        next_step = stage_dir / "geometry.in.next_step"
        if not next_step.is_file():
            raise StageError(
                "au_relax 缺 geometry.in.next_step：AIMS 弛豫没有产出新几何（检查 aims.out）")

        before = read_aims_geometry(stage_dir / "geometry.in")
        after = read_aims_geometry(next_step)
        if len(before) != len(after):
            raise StageError("au_relax 原子数变化（%d → %d）；弛豫不应增删原子"
                             % (len(before), len(after)))
        moved = self._moved_frozen(before, after, metadata)
        if moved:
            raise StageError(
                "au_relax 有 %d 个被背层冻结的原子发生了位移（如 1-based %s，最大 %.3g Å）；"
                "AIMS 没按 constrain_relaxation 固定它们，不能继续"
                % (len(moved), sorted(moved)[:5], max(moved.values())))
        shutil.copy2(next_step, stage_dir / "geometry.in")

        indices_ok, mismatch = self._verify_surface(stage_dir / "geometry.in", metadata)
        if not indices_ok:
            raise StageError(
                "au_relax 后 det-v1 表面三元组与 au60_build metadata 不一致: %s；"
                "背层冻结失效或电极塌陷，拒绝交下游" % mismatch)

        meta = {k: v for k, v in metadata.items() if k != "_path"}
        input_geometry = ctx.artifact("au60_build", "geometry.in")
        meta["relax"] = {
            "enabled": True,
            "converged": True,
            "force_threshold": ctx.config.electrode.relax_force_threshold,
            "input_fingerprint": _sha256(input_geometry),
            "relaxed_fingerprint": _sha256(stage_dir / "geometry.in"),
            "frozen_atoms": len(self._constrained_atoms(ctx, stage_dir / "geometry.in", metadata)),
        }
        (stage_dir / "electrode_metadata.json").write_text(
            json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
        attach_fingerprint(stage_dir)

        return self.collected(
            ctx, stage_dir,
            artifacts=[name for name in ("geometry.in", "geometry.in.next_step", "control.in",
                                         "aims.out", "sub_aims.sh", "input_summary.json",
                                         "electrode_metadata.json")
                       if (stage_dir / name).is_file()],
            results={"au_relax": {
                "enabled": True,
                "atom_count": len(after),
                "normal_termination": True,
                "force_threshold": ctx.config.electrode.relax_force_threshold,
                "surface_indices_verified": True,
            }},
            notes=["弛豫后几何与表面三元组已核对：%s" % (metadata.get("left_surface_indices"),)],
        )

    # -- 内部 -------------------------------------------------------------- #
    def _read_metadata(self, ctx, stage_dir: Path) -> dict:
        path = ctx.artifact("au60_build", "electrode_metadata.json")
        try:
            data = json.loads(Path(path).read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise StageBlocked(self.name, "读不到 au60_build 的 electrode_metadata.json: %s" % exc)
        data["_path"] = str(path)
        return data

    def _constrained_atoms(self, ctx, geometry: Path, metadata: dict) -> list[int]:
        """两侧背层整层（1-based）；metadata 缺字段时按几何现算（旧 run 兼容）。"""
        out: list[int] = []
        for key in BACKPLANE_KEYS:
            values = metadata.get(key) or []
            out.extend(int(i) for i in values)
        if out:
            return sorted(set(out))
        atoms = read_aims_geometry(geometry)
        for side in ("left", "right"):
            indices = metadata.get("%s_electrode_indices" % side) or []
            tip = metadata.get("%s_tip_index" % side)
            linker = metadata.get("%s_linker_index" % side)
            if not indices or tip is None or linker is None:
                continue
            out.extend(_backplane_indices(atoms, [int(i) for i in indices],
                                          int(tip), int(linker)))
        if not out:
            raise StageBlocked(
                self.name,
                "au60_build metadata 既没有背层索引，也无法按电极索引现算；"
                "请重跑 au60_build 生成 left/right_backplane_indices")
        return sorted(set(out))

    def _moved_frozen(self, before, after, metadata: dict) -> dict[int, float]:
        moved: dict[int, float] = {}
        frozen: list[int] = []
        for key in BACKPLANE_KEYS:
            frozen.extend(int(i) for i in (metadata.get(key) or []))
        for index in sorted(set(frozen)):
            if not (1 <= index <= len(before)):
                continue
            a = before[index - 1].coord()
            b = after[index - 1].coord()
            dist = ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2) ** 0.5
            if dist > MOVED_TOL:
                moved[index] = dist
        return moved

    def _verify_surface(self, geometry: Path, metadata: dict) -> tuple[bool, str]:
        """用 det-v1 规则重算两侧表面三元组并与 metadata 比对。

        注意 metadata 里的索引是 1-based，``electrode.surface_indices`` 内部按
        0-based 处理原子列表，直接传会整体错一位（实测会把 tip 自己选成"外层"）。
        这里用 ``surface_index.verify_surface_indices``，它负责 1-based 换算与
        atomes 分侧。
        """
        atoms = read_aims_geometry(geometry)
        for side in ("left", "right"):
            info = verify_surface_indices(atoms, metadata, side)
            if not info.get("ok"):
                return False, "%s: %s" % (side, "; ".join(info.get("reasons") or []))
        return True, ""

    def _species_dir(self, ctx) -> Path:
        species = ctx.species_dir()
        if species is None:
            raise StageBlocked(
                self.name,
                "缺少 FHI-aims species 目录：请在 smj.yaml 设 cluster.species_dir 或 paths.species_defaults")
        if not species.is_dir():
            raise StageBlocked(self.name, "FHI-aims species 目录不存在: %s" % species, [str(species)])
        return species


def _sha256(path: str | Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
