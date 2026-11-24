"""声子链三个远端阶段共用的作业脚本渲染。

与电子链的 ``junction.inputs._render_aims_sub`` 保持同一套环境变量口径
(``aims_module`` / ``aims_bin`` / ``aims_prefix`` / ``aims_mem`` / ``aims_walltime``),
但**不重复实现**分支逻辑:这里只拼 ``module_load`` 与 ``aims_bin`` 两段,
模板本身放在 ``junction/templates/`` 下,与 ``sub_aims.sh`` 并列。
"""

from __future__ import annotations

from pathlib import Path

from ..errors import StageBlocked
from ..junction.tcontrol import TEMPLATE_DIR, read_template
from ..slurm import gres_line as slurm_gres

#: 热导链的 Gaussian 占据展宽 [eV]。
#:
#: 电子链算的是有限团簇(有 HOMO-LUMO 隙),用 0.01 eV 没问题;热导器件是
#: Au(111) 二维板,是金属。0.01 eV 的展宽等于没有展宽,费米面附近的占据
#: 在自洽迭代里反复跳变,job 26375/26376 都是 200 多圈才能收敛到 1e-4。
#: 0.1 eV 是金属板的常规值(300 K 物理展宽是 0.026 eV):力常数要的是
#: ±0.01 Å 两个位移之间的**力差**,展宽只要在两者间保持一致就不引入偏差。
THERMAL_OCCUPATION_WIDTH_EV = 0.1


def upstream_file(ctx, stage: str, filename: str, prev: str, *,
                  override=None) -> Path:
    """取上游产物文件。

    热导链不在 :data:`smj.state.STAGE_ORDER` 里,``Context.upstream`` 内部靠
    ``previous_stage()`` 推导上游,对它必然报错;所以这里把 ``prev`` 显式传进来,
    其余行为(``--geometry`` 覆盖、缺失即 ``StageBlocked``)与电子链完全一致。
    """
    if override is not None:
        path = ctx.resolve(override)
        if path is None or not path.exists():
            raise StageBlocked(stage, "显式输入不存在: %s" % (override,),
                               [str(path)] if path else [])
        return path
    return ctx.require_stage_file(prev, filename, target=stage)


def aims_launcher(env: dict, ntasks: int) -> tuple[str, str, str]:
    """返回 ``(module_load, exe_define, aims_bin)``,与电子链分支一致。"""
    if str(env.get("aims_module") or "").strip():
        module_load = "module load %s" % str(env["aims_module"]).strip()
        exe_define = 'AIMS_EXE="${FHI_AIMS_BIN:?FHI-aims module did not define FHI_AIMS_BIN}"'
        command = 'mpirun -np %d "$AIMS_EXE" > aims.out 2>&1' % ntasks
    elif str(env.get("aims_bin") or "").strip():
        ldmod = str(env.get("aims_ld_module") or "").strip()
        module_load = ("module load %s" % ldmod) if ldmod else ""
        exe_define = ""
        command = 'mpirun -np %d "%s" > aims.out 2>&1' % (
            ntasks, str(env["aims_bin"]).strip())
    else:
        prefix = env.get("aims_prefix") or "UNCONFIGURED"
        module_load = ("export PATH={p}/install/bin:$PATH\n"
                       "export LD_LIBRARY_PATH={p}/install/lib:$LD_LIBRARY_PATH\n"
                       "export INCLUDE={p}/install/include:$INCLUDE").format(p=prefix)
        exe_define = ""
        command = "mpirun -np %d aims.x > aims.out 2>&1" % ntasks
    return module_load, exe_define, command


def _gres_line(gres) -> str:
    """生成 ``#SBATCH --gres`` 行 —— 规则本体在 :func:`smj.slurm.gres_line`。

    这里只是把共用的 ``--gres=...`` 片段补上 ``#SBATCH`` 前缀。纯 CPU 时必须
    显式要 ``gpu:0``：gpu1 节点声明了 ``Gres=gpu:4``，不写这一行作业会按整节点
    GRES 分配，把后面的作业全堵住（实测依据见 :data:`smj.slurm.DEFAULT_CPU_GRES`）。
    """
    fragment = slurm_gres(gres)
    return ("#SBATCH %s" % fragment) if fragment else ""


def render_thermal_sub(template: str, env: dict, *, ntasks: int,
                       extra: dict | None = None) -> str:
    """渲染 ``sub_thermal_*.sh``。"""
    module_load, exe_define, command = aims_launcher(env, ntasks)
    fields = {
        "partition": env.get("partition") or "UNCONFIGURED",
        "ntasks": ntasks,
        "mem": env.get("aims_mem") or "12G",
        "walltime": env.get("aims_walltime") or "06:00:00",
        "gres_line": _gres_line(env.get("gres")),
        "module_load": module_load,
        "exe_define": exe_define,
        "aims_bin": command,
    }
    fields.update(extra or {})
    return read_template(template).format(**fields)


def write_script(path: Path, text: str) -> Path:
    """写作业脚本。

    **必须显式 ``newline="\\n"``**：``Path.write_text`` 在 Windows 上会把
    ``\\n`` 翻成 ``\\r\\n``，而 gpu1 的 ``sbatch`` 直接拒收 DOS 行尾
    (``Batch script contains DOS line breaks``)。本仓其余写远端正文本的地方
    (``inputs._write_script`` / ``geometry.write_aims_geometry`` / ``tcontrol``)
    也都显式指定了 ``newline="\\n"``。
    """
    path = Path(path)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)
    return path


def template_path(name: str) -> Path:
    return Path(TEMPLATE_DIR) / name
