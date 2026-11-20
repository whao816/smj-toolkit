# -*- coding: utf-8 -*-
"""SMJ 干跑演示的内置夹具与预期边界。

`agenticdft demo --smj` 用真实引擎跑一条单分子结链，但只到文档记录的干跑边界：
多阶段 smj 链在 ``--dry-run`` 下没有上游 ``gaussian_opt`` 的 ``opt.xyz``，
``aims_tip_relax`` 会按 fail-closed 停在「缺上游 opt.xyz」
(``docs/smj-implementation-status.md`` §2.9)。演示把「到达该边界且
``gaussian_opt`` 输入已生成」作为成功口径，而不是伪装整链跑通。

内置分子取 BDT(1,4-苯二硫酚)中性起始几何，与
``docs/smj-real-cases-config.md`` §1.1 的 benchmark 输入逐原子一致；
锚定原子为 1-based 的 S11/S12，``strip_thiol_h`` 会移除 H13/H14。
"""
from pathlib import Path

BDT_NAME = 'BDT'

#: 14 原子 BDT 起始几何(与 docs/smj-real-cases-config.md §1.1 逐原子一致)。
BDT_XYZ = """14
BDT neutral HS-C6H4-SH start geometry (Gaussian b3lyp/6-31g(d) opt+freq input); thiol H at idx 13/14 to be stripped after opt
C  1.40000000  0.00000000  0.00000000
C  0.70000000  1.21243557  0.00000000
C -0.70000000  1.21243557  0.00000000
C -1.40000000  0.00000000  0.00000000
C -0.70000000 -1.21243557  0.00000000
C  0.70000000 -1.21243557  0.00000000
H  1.24500000  2.15640326  0.00000000
H -1.24500000  2.15640326  0.00000000
H -1.24500000 -2.15640326  0.00000000
H  1.24500000 -2.15640326  0.00000000
S  3.17000000  0.00000000  0.00000000
S -3.17000000  0.00000000  0.00000000
H  3.31006800 -1.33266000  0.00000000
H -3.31006800 -1.33266000  0.00000000
"""

#: 1-based 锚定原子: 两个 S(硫醇盐口径, 去 H 后 S 直连 Au tip)。
BDT_ANCHORS = (11, 12)

#: 干跑必然停下的阶段(见模块 docstring)。
EXPECTED_STOP_STAGE = 'aims_tip_relax'

EXPECTED_STOP_REASON = (
    '多阶段 smj 链的 --dry-run 不产出 gaussian_opt 的 opt.xyz，'
    'aims_tip_relax 按 fail-closed 停在「缺上游 opt.xyz」；'
    '这是演练语义限制，真跑不受影响'
)

#: 工作流计划里应出现的阶段链(仅用于展示; 以 results.workflow_plan 为准)。
DEFAULT_STAGE_CHAIN = (
    'gaussian_opt', 'aims_tip_relax', 'au60_build',
    'aims_scf', 'aitranss_te', 'junction_analyze',
)


def write_bdt_xyz(path):
    """把内置 BDT 几何写成 xyz，返回写入路径。"""
    path = Path(path)
    path.write_text(BDT_XYZ, encoding='utf-8')
    return path


def demo_overrides(structure_file):
    """`demo --smj` 的 TaskCard 覆盖: BDT + 硫醇盐口径 + 与验收台账同口径的窗口。"""
    return {
        'material': BDT_NAME,
        'properties': ['molecular_junction_transport'],
        'material_class': 'molecule',
        'structure_source': 'file',
        'structure_file': str(structure_file),
        'left_anchor': BDT_ANCHORS[0],
        'right_anchor': BDT_ANCHORS[1],
        'au_tip_distance': 2.5,
        'au_per_electrode': 20,
        'nlayers': 2,
        'energy_window': [-0.3, 0.0002, 0.0],
        'junction_temperature': 300.0,
        'strip_thiol_h': True,
    }


def stage_chain(results):
    """从 run 的 workflow_plan 取阶段名; 缺失时退回静态链。"""
    plan = (results or {}).get('workflow_plan') or {}
    names = []
    for stage in plan.get('stages') or []:
        if not isinstance(stage, dict):
            continue
        name = stage.get('name') or stage.get('stage')
        if name:
            names.append(str(name))
    return names or list(DEFAULT_STAGE_CHAIN)
