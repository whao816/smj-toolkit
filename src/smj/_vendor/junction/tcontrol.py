"""单分子结输入生成与解析:tcontrol 填充、FHI-aims control.in、species 查找、efermi/nsaos 解析。

来源: 第三方工具 qmtransport/control.py(约 112 行)+ qmtransport/tcontrol.py(54 行)
(tmp/自动化/, 只读来源不删除不入库)按本项目协议适配迁入(AGD-2026-09-22-SMJ, P2)。
修改(§4.3):efermi 解析顺序**维持 TE.dat 表头优先**(aitranss 实际使用的后验值, 与
透射谱自洽), 缺失回退 aims.out, 均缺 fail-closed 报错; 占位值(nsaos/valence_electrons
反解失败)一律 error 拦截, 不得静默填充。模板目录: 本仓库 agenticdft/agenticdft/tools/
junction/templates/(P4 输入写出器落模板; 缺失即报错, 不静默)。
纯本地工具: 只做输入生成与解析, 不做真实计算。
"""
from __future__ import annotations

from pathlib import Path
import json
import re
from .geometry import Atom, count_atoms, elements_in_geometry, read_aims_geometry


TEMPLATE_DIR = Path(__file__).resolve().parent / "templates"


def read_template(name: str) -> str:
    path = TEMPLATE_DIR / name
    if not path.exists():
        raise FileNotFoundError(f"模板缺失: {path} (P4 输入写出器应已落模板; 缺失即 fail-closed)")
    return path.read_text(encoding="utf-8")


def find_species_file(species_dir: str | Path, element: str) -> Path:
    root = Path(species_dir)
    patterns = [
        f"*_{element}_default",
        f"*_{element}_default.in",
        f"{element}_default",
        f"{element}.in",
        element,
    ]
    for pattern in patterns:
        matches = sorted(root.glob(pattern))
        if matches:
            return matches[0]
    raise FileNotFoundError(f"Cannot find FHI-aims species block for {element} in {species_dir}")


def species_blocks(atoms: list[Atom], species_dir: str | Path) -> str:
    blocks = []
    for element in elements_in_geometry(atoms):
        path = find_species_file(species_dir, element)
        blocks.append(path.read_text(encoding="utf-8").rstrip())
    return "\n\n".join(blocks) + "\n"


def vdw_line(vdw: str | None) -> str:
    """W0-5(AGD-2026-11-03-SMJ-HTP): 色散修正行 —— 缺省空串(输出与旧版逐字节一致)。

    `vdw_correction_hirshfeld` = Tkatchenko-Scheffler 色散修正(FHI-aims 221103 内建
    TS 参数表, 不需要 species 块里写 vdw_params —— 2026-11-05 在 gpu1 上用 Au2 实测
    确认)。**它是后自洽的**: 只改总能量与力, 不改 KS 哈密顿量(实测同几何、开/关两跑
    的 SCF 总能之差恰等于 TS 修正项)。因此它影响的是**弛豫几何**, 不影响固定几何下的
    T(E)。传未知值一律 ValueError(fail-closed, 不静默忽略)。
    """
    raw = str(vdw or '').strip().lower()
    if not raw or raw in ('none', 'off', '0', 'false'):
        return ''
    if raw in ('hirshfeld', 'ts', 'on', '1', 'true'):
        return '\n  vdw_correction_hirshfeld'
    raise ValueError('未知 vdw 取值 %r: 只支持 hirshfeld/ts(缺省关); 见 W0-5 裁决' % (vdw,))


def make_control(
    geometry_path: str | Path,
    output_path: str | Path,
    mode: str,
    species_dir: str | Path,
    xc: str = "pbe",
    spin: str = "none",
    charge: float = 0.0,
    force_threshold: str = "1.e-2",
    occupation_width: float = 0.01,
    fixed_spin_moment: float = 0.0,
    sc_iter_limit: int = 500,
    states: list[int] | None = None,
    vdw: str | None = None,
) -> None:
    atoms = read_aims_geometry(geometry_path)
    if mode == "relax":
        template_name = "control_relax.in"
        extra = {}
    elif mode == "scf":
        template_name = "control_scf.in"
        extra = {}
    elif mode == "orbital":
        template_name = "control_orbital.in"
        if not states:
            raise ValueError("--states is required for orbital mode")
        extra = {"state_outputs": "\n".join(f"  output cube eigenstate {state}" for state in states)}
    else:
        raise ValueError("mode must be relax, scf, or orbital")
    text = read_template(template_name).format(
        xc=xc,
        spin=spin,
        charge=charge,
        force_threshold=force_threshold,
        occupation_width=occupation_width,
        fixed_spin_moment=fixed_spin_moment,
        sc_iter_limit=sc_iter_limit,
        species_blocks=species_blocks(atoms, species_dir),
        vdw_line=vdw_line(vdw),
        **extra,
    )
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    Path(output_path).write_text(text, encoding="utf-8", newline="\n")


def parse_mos_nsaos(path: str | Path) -> int:
    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    patterns = [
        r"nsaos\s*[:=]?\s*(\d+)",
        r"number\s+of\s+(?:basis\s+functions|states|saos)\s*[:=]?\s*(\d+)",
        r"^\s*(\d+)\s+.*(?:nsaos|basis)",
    ]
    for pattern in patterns:
        match = re.search(pattern, text, flags=re.IGNORECASE | re.MULTILINE)
        if match:
            return int(match.group(1))
    raise ValueError(f"Cannot determine nsaos from {path}; pass --nsaos explicitly")


_VALENCE_BY_ELEMENT = {"H": 1, "He": 2, "Li": 1, "Be": 2, "B": 3, "C": 4, "N": 5, "O": 6,
                       "F": 7, "Ne": 8, "Na": 1, "Mg": 2, "Al": 3, "Si": 4, "P": 5, "S": 6,
                       "Cl": 7, "Ar": 8, "K": 1, "Ca": 2, "Sc": 3, "Ti": 4, "V": 5, "Cr": 6,
                       "Mn": 7, "Fe": 8, "Co": 9, "Ni": 10, "Cu": 11, "Zn": 12, "Ga": 3,
                       "Ge": 4, "As": 5, "Se": 6, "Br": 7, "Kr": 8, "Rb": 1, "Sr": 2,
                       "Y": 3, "Zr": 4, "Nb": 5, "Mo": 6, "Tc": 7, "Ru": 8, "Rh": 9,
                       "Pd": 10, "Ag": 11, "Cd": 12, "In": 3, "Sn": 4, "Sb": 5, "Te": 6,
                       "I": 7, "Xe": 8, "Cs": 1, "Ba": 2, "La": 3, "Ce": 4, "Pr": 5,
                       "Nd": 6, "Pm": 7, "Sm": 8, "Eu": 9, "Gd": 10, "Tb": 11, "Dy": 12,
                       "Ho": 13, "Er": 14, "Tm": 15, "Yb": 16, "Lu": 17, "Hf": 4, "Ta": 5,
                       "W": 6, "Re": 7, "Os": 8, "Ir": 9, "Pt": 10, "Au": 31, "Hg": 12,
                       "Tl": 3, "Pb": 4, "Bi": 5, "Po": 6, "At": 7, "Rn": 8}


def _auto_valence_electrons(geometry: str | Path) -> int:
    """按 FHI-aims light 默认价电子表自动计算总价电子。
    Au=31(6s1 5p6 5d10 4f14, FHI-aims light 将 4f 计入价层),
    其余元素按主族价电子。用于 $valence_electrons(aitranss 电荷守恒)。"""
    from .geometry import read_aims_geometry
    atoms = read_aims_geometry(geometry)
    total = sum(_VALENCE_BY_ELEMENT.get(a.element, 0) for a in atoms)
    if total <= 0:
        raise ValueError("Cannot determine valence electrons from geometry; pass explicitly")
    return total


def parse_efermi(paths: list[str | Path]) -> float:
    """efermi 解析(§4.3 纪律): 按传入顺序逐个文件匹配; 调用方须把 TE.dat 放最前。
    TE.dat 表头 efermi 是 aitranss 实际使用的后验值(与透射谱自洽), 优先于 aims.out。
    支持 Ha/hartree(原样)与 eV(除以 27.2114 换算成 Ha, aitranss 期望 Ha)。
    全部缺失/无匹配即 fail-closed 报错, 不得用占位值静默填充。"""
    patterns = [
        (r"efermi\s*=\s*([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)", 1.0),
        (r"Fermi\s+level.*?([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)\s*(?:Ha|hartree)", 1.0),
        (r"Chemical\s+potential.*?([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)\s*(?:Ha|hartree)", 1.0),
        (r"Chemical\s+potential.*?([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)\s*eV", 1.0 / 27.2114),
        (r"Fermi\s+level.*?([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)\s*eV", 1.0 / 27.2114),
    ]
    for path in paths:
        p = Path(path)
        if not p.exists():
            continue
        text = p.read_text(encoding="utf-8", errors="ignore")
        for pattern, scale in patterns:
            matches = re.findall(pattern, text, flags=re.IGNORECASE)
            if matches:
                return float(matches[-1]) * scale
    raise ValueError("Cannot determine efermi; pass --efermi explicitly")


def make_tcontrol(
    scf_dir: str | Path,
    output: str | Path,
    metadata: str | Path | None = None,
    nsaos: int | None = None,
    efermi: float | None = None,
    surfaces: list[int] | None = None,
    nlayers: int = 4,
    ener: float = -0.4,
    estep: float = 0.0001,
    eend: float = 0.0,
    valence_electrons: int = 0,
) -> dict:
    scf = Path(scf_dir)
    geometry = scf / "geometry.in"
    if nsaos is None:
        nsaos = parse_mos_nsaos(scf / "mos.aims")
    if efermi is None:
        # §4.3: TE.dat 表头优先(aitranss 实际使用值), 缺失回退 aims.out, 均缺 fail-closed。
        efermi = parse_efermi([scf / "TE.dat", scf / "aims.out", scf / "scf.out", scf / "out"])
    if surfaces is None:
        if not metadata:
            raise ValueError("Provide --metadata or six --surfaces values")
        data = json.loads(Path(metadata).read_text(encoding="utf-8"))
        surfaces = data["left_surface_indices"] + data["right_surface_indices"]
    if len(surfaces) != 6:
        raise ValueError("surfaces must contain six atom indices: lsurc lsurx lsury rsurc rsurx rsury")
    if valence_electrons <= 0:
        # FHI-aims light 默认价电子(与 control.in species 块 valence 定义一致)
        valence_electrons = _auto_valence_electrons(geometry)
    values = {
        "natoms": count_atoms(geometry),
        "nsaos": nsaos,
        "lsurc": surfaces[0],
        "lsurx": surfaces[1],
        "lsury": surfaces[2],
        "rsurc": surfaces[3],
        "rsurx": surfaces[4],
        "rsury": surfaces[5],
        "nlayers": nlayers,
        "ener": f"{ener:.10f}",
        "estep": f"{estep:.10f}",
        "eend": f"{eend:.10f}",
        "valence_electrons": valence_electrons,
        "efermi": f"{efermi:.12f}",
    }
    text = read_template("tcontrol.template").format(**values)
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    Path(output).write_text(text, encoding="utf-8", newline="\n")
    return values