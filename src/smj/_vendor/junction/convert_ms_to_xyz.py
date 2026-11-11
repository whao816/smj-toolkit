"""Materials Studio XSD → XYZ 低耦合转换器(AGD-2026-09-22-SMJ, 可选适配)。

来源: 第三方工具 convert_ms_to_xyz.py(198 行, tmp/自动化/, 只读来源不删除不入库)
按计划 §5.1 收窄迁入:**仅覆盖 XSD → XYZ**, 不含 CIF 转换能力(CIF 输入不在本期
承诺, 需另立课题)。复用 tools/junction/geometry.py 的 Atom/write_xyz, 避免重复承载。
纯本地工具, 不做真实计算。
"""
from __future__ import annotations

import re
from pathlib import Path
import xml.etree.ElementTree as ET

from .geometry import Atom, write_xyz


def clean_element(text: str) -> str:
    match = re.match(r"([A-Za-z]{1,2})", (text or "").strip())
    if not match:
        raise ValueError("XSD 原子元素名无法解析(非字母): %r" % (text or "")[:20])
    raw = match.group(1)
    return raw[0].upper() + raw[1:].lower()


def read_xsd(path: str | Path) -> list[Atom]:
    """解析 Materials Studio XSD(Atom3d 节点, XYZ 属性)为笛卡尔原子列表。"""
    root = ET.parse(path).getroot()
    atoms: list[Atom] = []
    for node in root.iter():
        tag = node.tag.split("}")[-1]
        if tag != "Atom3d":
            continue
        element = (node.attrib.get("Components") or node.attrib.get("ElementSymbol")
                   or node.attrib.get("Element") or node.attrib.get("ForcefieldType"))
        xyz = node.attrib.get("XYZ") or node.attrib.get("Xyz")
        if not xyz:
            continue
        values = [float(v) for v in re.split(r"[,\s]+", xyz.strip()) if v]
        if len(values) >= 3:
            atoms.append(Atom(clean_element(element), values[0], values[1], values[2]))
    if not atoms:
        raise ValueError("XSD 文件中未找到 Atom3d XYZ 坐标")
    return atoms


def convert_xsd_to_xyz(input_path: str | Path, output_path: str | Path | None = None) -> Path:
    """XSD → XYZ。output_path 缺省为输入同名 .xyz。"""
    inp = Path(input_path)
    if inp.suffix.lower() != ".xsd":
        raise ValueError("仅支持 .xsd 输入(CIF 输入不在本期承诺, 见计划 §5.1)")
    atoms = read_xsd(inp)
    out = Path(output_path) if output_path else inp.with_suffix(".xyz")
    write_xyz(out, atoms, comment="converted from %s" % inp.name)
    return out