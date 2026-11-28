"""Protocol and artifact checks shared by the SIESTA runner and report."""
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from ..errors import StageError
from ..provenance import canonical_json, sha256_text
from .io import sha


def stage_protocol(config, name):
    """Only settings affecting this calculation; parent digests cover inherited data."""
    settings = dict(config.siesta)
    pseudos = settings.pop("pseudopotentials", {})
    settings.pop("electrode_reuse", None)
    # Routing switches only: they select collected upstream artifacts, they do not
    # change the physical inputs of the stage whose protocol is being compared.
    settings.pop("reuse_relax_from", None)
    settings.pop("reuse_scf_from", None)
    contact = settings.pop("contact", {})
    electrode = settings.pop("electrode", {})
    tbgrid = settings.pop("tbtrans_kgrid", settings.get("device_kgrid"))
    if name != "molecule_opt":
        settings.pop("molecule_vacuum_ang", None)
    if name not in ("molecule_opt", "siesta_junction_relax"):
        settings.pop("relax_steps", None)
        settings.pop("force_ev_ang", None)
    if name != "siesta_electrode":
        settings.pop("electrode_kgrid", None)
    if name in ("molecule_opt", "siesta_electrode"):
        settings.pop("device_kgrid", None)
        settings.pop("device_mixing_weight", None)
    if name == "siesta_electrode":
        settings["basis"] = {"Au": settings.get("basis", {}).get("Au")}
    settings["pseudopotentials"] = {
        el: sha(path) if path and path.is_file() else None
        for el, raw in pseudos.items() if name != "siesta_electrode" or el == "Au"
        for path in [config.resolve(raw)]
    }
    if name == "molecule_opt":
        path = config.molecule_xyz()
        settings.update(molecule=asdict(config.molecule), xyz=sha(path) if path and path.is_file() else None,
                        anchors=[config.electrode.left_anchor, config.electrode.right_anchor])
    elif name == "siesta_electrode":
        settings["geometry"] = {k: v for k, v in electrode.items()
                                if k not in ("slab_layers", "movable_layers")}
    elif name == "siesta_junction_relax":
        settings.update(electrode=electrode, contact=contact, strip_h=config.electrode.strip_thiol_h,
                        anchors=[config.electrode.left_anchor, config.electrode.right_anchor],
                        initial_moments=config.molecule.initial_moments)
    elif name == "tbtrans_te":
        settings = {"kgrid": tbgrid, "window": config.transport.relative_window()}
    elif name == "junction_analyze":
        settings = {"temperature": config.transport.temperature}
    return sha256_text(canonical_json(settings))


def check_outputs(layout, name, record):
    """Include large binary files (HSX and TBT.nc) in the chain of custody."""
    if record.status != "collected" or not record.outputs_sha256:
        raise StageError(f"{name}: 上游尚未完成采集")
    for filename, expected in record.outputs_sha256.items():
        path = layout.artifact(name, filename)
        if not path.is_file() or sha(path) != expected:
            raise StageError(f"{name}: 已采集产物缺失或发生变化: {filename}")


def parents(ctx, stage, verified=None):
    from . import get_stage
    verified = {} if verified is None else verified
    records = {}
    deps = (stage.dependencies(ctx) if hasattr(stage, "dependencies")
            else stage.depends)
    for name in deps:
        if name in verified:
            records[name] = verified[name]
            continue
        record = ctx.state.stage(name)
        check_outputs(ctx.layout, name, record)
        if record.extra.get("siesta_protocol") != stage_protocol(ctx.config, name):
            raise StageError(f"{name}: 协议已变化，请先重新运行该阶段")
        if parents(ctx, get_stage(name), verified) != record.extra.get("siesta_parents"):
            raise StageError(f"{name}: 上游产物已变化，请重新运行该阶段")
        verified[name] = records[name] = record.outputs_sha256
    return records


def check_prepared(ctx, stage):
    record = ctx.state.stage(stage.name)
    if record.extra.get("siesta_protocol") != stage_protocol(ctx.config, stage.name):
        raise StageError(f"{stage.name}: 输入协议已变化，请重新 prepare")
    if parents(ctx, stage) != record.extra.get("siesta_parents"):
        raise StageError(f"{stage.name}: 上游产物已变化，请重新 prepare")
    for name, expected in record.inputs_sha256.items():
        path = Path(name)
        if not path.is_absolute():
            path = ctx.layout.artifact(stage.name, name)
        if not path.is_file() or sha(path) != expected:
            raise StageError(f"{stage.name}: 准备后的输入发生变化: {name}")
