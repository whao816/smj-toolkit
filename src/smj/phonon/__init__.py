"""单分子结声子热导(``G_ph``):力常数 → 动力学矩阵 → 透射 → Landauer 积分。

与 ``smj.junction``(电子输运)并列的独立内核,不共享代码也不改变电子链的任何数字。
设计文档见 ``docs/thermal-conductance-plan-2026-11-23.md``,
开发计划见 ``docs/thermal-conductance-devplan-2026-11-23.md``。

模块划分(方案 §6.1):
- ``units``       唯一换算源(力常数 → meV、能量 → cm⁻¹、单通道量子 G_Q)
- ``negf``        表面格林函数(decimation)与 Caroli 透射
- ``landauer``    声子 Landauer 线性响应积分(唯一公式源)
- ``dynamical``   质量加权动力学矩阵、虚频、声学求和残差
- ``forces``      FHI-aims 力块解析
- ``fc``          有限位移 → 力常数矩阵
- ``device``      周期板器件几何(含 S 去质子化)
"""
