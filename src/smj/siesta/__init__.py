"""Zero-bias SIESTA / TranSIESTA / TBtrans workflow."""

from ..errors import StageError


def get_stage(name):
    if name in ("tbtrans_te", "junction_analyze"):
        from .spectrum import STAGES
    else:
        from .stages import STAGES
    if name not in STAGES:
        raise StageError(f"{name} 不属于 SIESTA 电子输运路线")
    return STAGES[name]


def remote_outputs(ctx, stage):
    definition = get_stage(stage)
    if hasattr(definition, "remote_files_for"):
        return tuple(definition.remote_files_for(ctx))
    return tuple(definition.remote_files)
