"""au_relax：可选的贴装后几何弛豫（远端 AIMS 阶段，0.1.2 新增）。

默认协议不在电极扩展后再弛豫：``aims_tip_relax`` 已做全原子优化，
``au60_build`` 只把每个 tip 扩展成完整金字塔团簇并沿用弛豫后的接触几何，
此后结几何完全冻结。因此 ``electrode.relax`` 默认 **false**，本阶段只做本地
直通（复制 au60 的几何与 metadata，写 ``relax_skipped.json``）。

仅当显式 ``electrode.relax: true``（旧口径/收敛检查）时：
- prepare 复制 ``au60_build/geometry.in``，给两侧背层整层加
  ``constrain_relaxation .true.``，写 ``control.in``（bfgs 模板）+ ``sub_aims.sh``；
- collect 要求 ``aims.out`` 正常结束且 ``geometry.in.next_step`` 存在，用
  next_step 覆盖 ``geometry.in``；随后按 det-v1 重算表面三元组，必须与
  ``au60_build`` 的 metadata 完全一致（背层冻结后表面原子不该移动）——
  不一致即 StageError，不把"看起来跑完了"的几何交下游。
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
from ..provenance import sha256_file, verify_parent_geometry
from ._helpers import est_record, submit_command, vdw_mode
from .base import Collected, Prepared, Stage, register

SURFACE_KEYS = ("left_surface_indices", "right_surface_indices")
BACKPLANE_KEYS = ("left_backplane_indices", "right_backplane_indices")
SKIP_MARKER = "relax_skipped.json"
MOVED_TOL = 1e-4  # 被冻结的原子在弛豫前后应保持原位（Å）


def _gold_indices(geometry: Path) -> list[int]:
    """geometry.in 里所有 Au 原子的 1-based 下标。"""
    out: list[int] = []
    lines = Path(geometry).read_text(encoding="utf-8", errors="replace").splitlines()
    order = 0
    for line in lines:
        parts = line.split()
        if len(parts) >= 5 and parts[0] == "atom":
            order += 1
            if parts[4].strip().lower() == "au":
                out.append(order)
    return out


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
        # P0-3：先证明"要弛豫的几何"就是 au60_build 指纹记录的那一份。
        verify_parent_geometry(self.name, upstream, metadata)
        natoms = self.count_aims_atoms(geometry)

        if not bool(getattr(ctx.config.electrode, "relax", False)):
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
                notes=["electrode.relax=false：电极扩展后几何冻结，直接沿用 au60_build 几何",
                       "如需旧口径的贴装后弛豫，显式设 electrode.relax: true"],
                commands=[],
                local_only=True,
                parent_stage="au60_build",
                parent_fingerprint=sha256_file(upstream),
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
        freeze_mode = str(getattr(ctx.config.electrode, "relax_freeze", "backplane") or "backplane")
        ctx.log("au_relax: 已写 control.in（%d 原子，relax_freeze=%s，冻结 %d 个原子）"
                % (natoms, freeze_mode, len(constrained)))
        return self.prepared(
            ctx, stage_dir,
            inputs=[geometry, Path(written["control"]), Path(written["submit_script"]),
                    stage_dir / "input_summary.json"],
            resources=est_record(estimate),
            notes=[
                "%d 原子，relax_freeze=%s，冻结 %d 个（left 背层=%s right 背层=%s）"
                % (natoms, freeze_mode, len(constrained),
                   metadata.get("left_backplane_indices"), metadata.get("right_backplane_indices")),
                "force_threshold=%s，ntasks=%s，partition=%s"
                % (ctx.config.electrode.relax_force_threshold,
                   env.get("ntasks"), env.get("partition")),
                "prepare 只写输入与作业脚本，未提交任何远端命令",
            ],
            commands=submit_command(self.script),
            parent_stage="au60_build",
            parent_fingerprint=sha256_file(upstream),
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
        frozen = self._constrained_atoms(ctx, stage_dir / "geometry.in", metadata)
        moved = self._moved_frozen(before, after, frozen)
        if moved:
            raise StageError(
                "au_relax 有 %d 个被冻结的原子发生了位移（如 1-based %s，最大 %.3g Å）；"
                "AIMS 没按 constrain_relaxation 固定它们，不能继续"
                % (len(moved), sorted(moved)[:5], max(moved.values())))
        shutil.copy2(next_step, stage_dir / "geometry.in")

        indices_ok, mismatch = self._verify_surface(before, after, metadata)
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
            "frozen_atoms": len(frozen),
            "freeze_mode": str(getattr(ctx.config.electrode, "relax_freeze", "backplane") or "backplane"),
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
        """要冻结的原子（1-based）。

        ``electrode.relax_freeze``：

        - ``backplane``（默认，手册 §5.3.1 口径）：两侧背层整层；
        - ``electrodes``：全部 Au 原子（含 tip），只让分子骨架弛豫 ——
          适合电极已按理想晶格搭好、不想让 60+ 个 Au 参与 BFGS 的情形。

        锚定原子（Au–S 接触的 S）**不冻结**：薁 demo 实测，冻住 S 会让 C–S–Au
        的应变全压到分子骨架上（0.4–0.5 eV/Å 一直不降），而放开 S、只冻 Au 时
        S 会自然扭到 Au 表面上方（Au–S ≈ 2.56 Å），力很快降到 0.02 量级。
        S 移动会让 tip→linker 轴偏转，所以 collect 的表面校验改成"三个表面原子
        坐标没动即通过"（见 :meth:`_verify_surface`），不再用被带偏的轴去重算。
        """
        mode = str(getattr(ctx.config.electrode, "relax_freeze", "backplane") or "backplane")
        if mode == "electrodes":
            out = _gold_indices(geometry)
            if not out:
                raise StageBlocked(self.name, "geometry.in 里没有 Au 原子，relax_freeze=electrodes 无意义")
            return sorted(set(out))
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

    @staticmethod
    def _anchor_indices(metadata: dict) -> list[int]:
        """两侧锚定原子（linker，1-based）: 冻结它们以保持 Au–S 接触距离。"""
        out: list[int] = []
        for key in ("left_linker_index", "right_linker_index"):
            value = metadata.get(key)
            if value is not None:
                out.append(int(value))
        return sorted(set(out))


    def _moved_frozen(self, before, after, frozen_indices) -> dict[int, float]:
        moved: dict[int, float] = {}
        for index in sorted(set(int(i) for i in frozen_indices)):
            if not (1 <= index <= len(before)):
                continue
            a = before[index - 1].coord()
            b = after[index - 1].coord()
            dist = ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2) ** 0.5
            if dist > MOVED_TOL:
                moved[index] = dist
        return moved

    def _verify_surface(self, before, after, metadata: dict) -> tuple[bool, str]:
        """核对 det-v1 表面三元组在弛豫后仍然有效。

        判据分两层：

        1. metadata 给的三个索引必须是**同一元素（电极元素）、非共线**的真实三角；
        2. 三个原子在弛豫前后坐标未动时直接通过（电极被冻结时本来就不该动）——
           这直接说明三元组仍指向同一块几何，且不受 tip→linker 轴被分子扭转带偏的影响。
           只有表面原子确实移动了（backplane 模式只冻背层时）才按 det-v1 重算，
           并要求与 metadata 一致。
        """
        for side in ("left", "right"):
            given = [int(i) for i in (metadata.get("%s_surface_indices" % side) or [])]
            if len(given) != 3:
                return False, "%s: metadata 表面索引不是 3 个: %r" % (side, given)
            if not self._triple_ok(after, metadata, side, given):
                return False, "%s: metadata 表面三元组 %s 不是有效的电极三角" % (side, given)
            moved = self._moved_frozen(before, after, given)
            if not moved:
                continue
            info = verify_surface_indices(after, metadata, side)
            if not info.get("ok"):
                return False, "%s: 表面原子移动 %s 且重算失败: %s" % (
                    side, sorted(moved), "; ".join(info.get("reasons") or []))
        return True, ""

    @staticmethod
    def _triple_ok(atoms, metadata: dict, side: str, given: list[int]) -> bool:
        tip = int(metadata["%s_tip_index" % side]) - 1
        if any(not (1 <= i <= len(atoms)) for i in given):
            return False
        if any(atoms[i - 1].element != atoms[tip].element for i in given):
            return False
        a, b, c = (atoms[i - 1].coord() for i in given)
        u = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
        v = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
        cross = (u[1] * v[2] - u[2] * v[1],
                 u[2] * v[0] - u[0] * v[2],
                 u[0] * v[1] - u[1] * v[0])
        area2 = (cross[0] ** 2 + cross[1] ** 2 + cross[2] ** 2) ** 0.5
        return area2 > 2.0e-6

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
