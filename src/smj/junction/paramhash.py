"""高通量参数指纹(W3-3): 同一套物理输入必须得到同一个哈希, 否则"重复率"没法审计。

设计约束(与框架 §3 数据模型对齐):

* **只哈希物理输入**, 不哈希时间戳 / 路径 / 主机名 / 代码版本里易变的部分;
* 键顺序固定(``sort_keys=True``), 浮点先量化(``_q``), 避免 1e-16 噪声改变哈希;
* 计算方式与"人类可读的参数字典"共用同一份 ``canonical_params()``, 报告里能直接打印;
* 仓库里**已有** ``surface_fingerprint``(几何字节) —— 那是"这一份几何"的指纹,
  本模块是"这次计算的定义"的指纹, 两者互补, 都进 provenance。

用法::

    h = params_hash(molecule_id='azobenzene', isomer='E', anchor_pair='S-Au',
                    contact_motif='atop', electrode='Au', d_elec=2.45,
                    method='au60-pbe-light-consteta', extra={'estep': 1e-4})

``extra`` 里的键会被压实并纳入哈希(例如 ``sigma_near``、``layers``), 但**不**接受
``None`` 之外的容器嵌套 —— 保持扁平才能在报告里逐行列出。
"""

from __future__ import annotations

import hashlib
import json
import math
from typing import Any, Mapping

#: 参与哈希的标量字段(顺序即报告里的展示顺序)。
FIELDS = (
    'molecule_id',
    'isomer',
    'anchor_pair',
    'contact_motif',
    'electrode',
    'd_elec',
    'method',
)

#: 浮点量化位数: 距离写到 1e-3 A, 能量/网格等写进 extra 时由调用方自行量化。
_DIGITS = 6


def _q(value: Any) -> Any:
    """把浮点压到 6 位小数; 其余类型原样返回(便于 JSON 序列化)。"""
    if isinstance(value, bool):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError('参数里有非有限浮点数: %r' % (value,))
        return round(value, _DIGITS)
    if isinstance(value, (int, str)) or value is None:
        return value
    raise TypeError('参数值必须是标量(数字/字符串/None), 实为 %r' % (type(value).__name__,))


def canonical_params(molecule_id: str, isomer: str, anchor_pair: str,
                     contact_motif: str, electrode: str, d_elec: float,
                     method: str,
                     extra: Mapping[str, Any] | None = None) -> dict:
    """返回**规范化后**的参数字典(报告里直接打印这一份)。"""
    for name, value in (('molecule_id', molecule_id), ('isomer', isomer),
                        ('anchor_pair', anchor_pair), ('contact_motif', contact_motif),
                        ('electrode', electrode), ('method', method)):
        if not isinstance(value, str) or not value.strip():
            raise ValueError('%s 必须是非空字符串, 实为 %r' % (name, value))
    if not isinstance(d_elec, (int, float)) or isinstance(d_elec, bool):
        raise ValueError('d_elec 必须是数字, 实为 %r' % (d_elec,))
    if not 1.0 <= float(d_elec) <= 5.0:
        raise ValueError('d_elec 超出合理范围 1.0-5.0 A: %r' % (d_elec,))

    out = {'molecule_id': molecule_id.strip().lower(),
           'isomer': isomer.strip().upper(),
           'anchor_pair': anchor_pair.strip(),
           'contact_motif': contact_motif.strip().lower(),
           'electrode': electrode.strip().capitalize(),
           'd_elec': _q(float(d_elec)),
           'method': method.strip()}
    if out['isomer'] not in ('E', 'Z'):
        raise ValueError("isomer 只能是 E 或 Z, 实为 %r" % (out['isomer'],))
    if out['contact_motif'] not in ('atop', 'bridge', 'hollow'):
        raise ValueError('contact_motif 只能是 atop/bridge/hollow, 实为 %r'
                         % (out['contact_motif'],))
    if extra:
        for key, value in extra.items():
            if key in out:
                raise ValueError('extra 不能覆盖标准字段: %s' % key)
            if isinstance(value, (dict, list, tuple, set)):
                raise TypeError('extra[%s] 必须是标量(扁平化后再传)' % key)
            out[key] = _q(value)
    return out


def params_hash(*args: Any, **kwargs: Any) -> str:
    """规范化参数的 sha256(取前 16 位十六进制, 足够审计且不刷屏)。"""
    canon = canonical_params(*args, **kwargs)
    blob = json.dumps(canon, sort_keys=True, ensure_ascii=True,
                      separators=(',', ':')).encode('utf-8')
    return hashlib.sha256(blob).hexdigest()[:16]


def file_hash(path: str) -> str:
    """内容哈希(用于 ``result_hash``: 同一份 TE.dat/报告必须同哈希)。"""
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()[:16]
