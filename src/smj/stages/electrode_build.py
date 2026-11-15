"""electrode_build（构建电极）：给分子两端装配电极团簇（本地阶段，不提交）。

阶段参数全部来自 ``electrode.*``：
- ``source``: ``builtin`` = 本地按晶格常数生成 fcc(111) 三角岛；
  ``library`` = 调用 FHI-aims ``electrodes.library`` 预定义团簇；
- ``element`` / ``a_lattice``: 电极金属（Au/Ag/Cu/Pt/Pd/Ni/Al…）与晶格常数。

输入是 ``aims_tip_relax/geometry.in.next_step``（分子 + 两枚电极金属 tip），
输出是补全两侧电极后的 ``geometry.in`` 与 ``electrode_metadata.json``。
prepare 即产出；collect 只做一致性解析/汇总，不再改文件。
"""

from __future__ import annotations

import json
from pathlib import Path

from .._vendor.junction.electrode import append_au60
from .._vendor.junction.fingerprint import attach_fingerprint
from ..electrode_library import resolve_cluster
from ..errors import StageError
from ..provenance import sha256_file
from ._helpers import est_record
from .base import Collected, Prepared, Stage, register


@register
class ElectrodeBuildStage(Stage):
    name = "electrode_build"
    depends = ("aims_tip_relax",)
    needs_submit = False
    input_hint = "消费 aims_tip_relax/geometry.in.next_step"
    outputs = ("geometry.in", "electrode_metadata.json")

    # -- prepare（本地阶段：直接产出） -------------------------------------- #
    def prepare(self, ctx, stage_dir: Path) -> Prepared:
        stage_dir = Path(stage_dir)
        # --geometry 覆盖：阶段独立运行时直接给"分子 + 两枚电极 tip"的几何文件。
        source = ctx.upstream(
            self.name, "geometry.in.next_step", override=ctx.override_path("geometry"))
        electrode = ctx.config.electrode
        element = electrode.element
        geometry = stage_dir / "geometry.in"
        metadata_path = stage_dir / "electrode_metadata.json"
        notes: list[str] = []
        cluster = None
        if getattr(electrode, "source", "builtin") == "library":
            # 0.1.2: 用 FHI-aims electrodes.library 的预定义团簇（手册 §5.3.1）。
            cluster = resolve_cluster(ctx, electrode)
            if cluster.element != element:
                raise StageError(
                    "electrode_build 电极元素不一致：electrode.element=%s，但库团簇 %s 是 %s；"
                    "请改 electrode.library_cluster 或 electrode.element（fail-closed）"
                    % (element, cluster.relpath, cluster.element))
            if (electrode.nlayers and cluster.nlayers_suggested
                    and int(electrode.nlayers) != int(cluster.nlayers_suggested)):
                notes.append(
                    "警告: tcontrol nlayers=%s 与库文件建议值 %s 不一致 —— 手册要求取"
                    "库文件头部建议值，请核对 electrode.nlayers"
                    % (electrode.nlayers, cluster.nlayers_suggested))
        # 0.2.0rc1：library 模式完全忽略 au_per_electrode（团簇大小由库文件决定）。
        au_target = None if cluster is not None else electrode.au_per_electrode
        try:
            metadata = append_au60(
                source, geometry, metadata_path,
                au_per_electrode=au_target,
                element=element,
                a_lattice=electrode.a_lattice,
                cluster=cluster,
                provenance=cluster.provenance() if cluster is not None else None,
            )
        except ValueError as exc:
            raise StageError("electrode_build 无法装配电极: %s" % exc) from exc
        # 指纹写回 metadata：surface_fingerprint = sha256(geometry.in)，
        # 下游 review 用它确认"确认过的表面"和"实际算的表面"是同一份文件。
        metadata = attach_fingerprint(stage_dir)
        natoms = self.count_aims_atoms(geometry)
        estimate = ctx.estimate(self.name, natoms)
        surface = metadata.get("surface_indices") or {}
        left_surface = surface.get("left_surface_indices") or metadata.get("left_surface_indices")
        right_surface = surface.get("right_surface_indices") or metadata.get("right_surface_indices")
        ctx.log("electrode_build: %d → %d 原子（本地阶段，%s电极，%s）"
                % (metadata.get("original_atom_count", "?"),
                   metadata.get("final_atom_count", natoms),
                   element,
                   "库团簇" if cluster is not None else "自建三角岛"))
        if cluster is not None:
            notes.append("库团簇: %s（%d 原子，%s，建议 nlayers=%s，sha256 %s…）"
                         % (cluster.relpath, cluster.natoms, cluster.element,
                            cluster.nlayers_suggested, cluster.sha256[:12]))
            notes.append(
                "au_per_electrode 仅对 source: builtin 生效；library 模式按库团簇整簇贴装"
                "（该键已忽略）")
        notes.extend([
            "本地阶段：在 aims_tip_relax 的几何上补全两侧 %s(111) 电极，无需集群提交" % element,
            "%s 每侧目标 %s 个（%s），最终 %d 原子"
            % (element, metadata.get("au_per_electrode"),
               "库团簇整簇贴装" if cluster is not None else "按完整原子层截断",
               metadata.get("final_atom_count", natoms)),
            "表面原子（1-based）left=%s right=%s，背层 left=%s right=%s，指纹 %s…"
            % (left_surface, right_surface,
               metadata.get("left_backplane_indices"), metadata.get("right_backplane_indices"),
               str(metadata.get("surface_fingerprint"))[:12]),
        ])
        return self.prepared(
            ctx, stage_dir,
            inputs=[source],
            artifacts=list(self.outputs),
            resources=est_record(estimate),
            notes=notes,
            commands=[],
            parent_stage="aims_tip_relax",
            parent_fingerprint=sha256_file(source),
        )

    # -- collect（解析/汇总） ---------------------------------------------- #
    def collect(self, ctx, stage_dir: Path) -> Collected:
        stage_dir = Path(stage_dir)
        geometry = ctx.require_stage_file(self.name, "geometry.in")
        meta_path = ctx.require_stage_file(self.name, "electrode_metadata.json")
        try:
            metadata = json.loads(meta_path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise StageError("electrode_build 无法解析 electrode_metadata.json: %s" % exc) from exc
        if not isinstance(metadata, dict):
            raise StageError("electrode_build 的 electrode_metadata.json 不是 JSON 对象")
        expected = metadata.get("final_atom_count")
        actual = self.count_aims_atoms(geometry)
        if expected is not None and int(expected) != actual:
            raise StageError(
                "electrode_build 原子数不一致：metadata 记 %s，geometry.in 实际 %d"
                % (expected, actual))
        fingerprint = metadata.get("surface_fingerprint")
        surface = metadata.get("surface_indices") or {}
        results = {
            "electrode_build": {
                "final_atom_count": actual,
                "original_atom_count": metadata.get("original_atom_count"),
                "au_per_electrode": metadata.get("au_per_electrode"),
                "left_tip_index": metadata.get("left_tip_index"),
                "right_tip_index": metadata.get("right_tip_index"),
                "surface_fingerprint": fingerprint,
                "electrode_source": metadata.get("electrode_source") or "builtin",
                "electrode_element": metadata.get("electrode_element")
                or ctx.config.electrode.element,
                "electrode_a_lattice": metadata.get("electrode_a_lattice")
                or ctx.config.electrode.a_lattice,
                "library_cluster": metadata.get("library_cluster"),
                "backplane_indices": {
                    "left": metadata.get("left_backplane_indices"),
                    "right": metadata.get("right_backplane_indices"),
                },
                "surface_indices": {
                    "left": surface.get("left_surface_indices")
                    or metadata.get("left_surface_indices"),
                    "right": surface.get("right_surface_indices")
                    or metadata.get("right_surface_indices"),
                },
            }
        }
        notes = ["本地阶段产物已就绪；surface_fingerprint 供审阅确认使用" if fingerprint
                 else "electrode_metadata.json 缺 surface_fingerprint：审阅确认无法绑定几何"]
        return self.collected(
            ctx, stage_dir, artifacts=list(self.outputs), results=results, notes=notes)
