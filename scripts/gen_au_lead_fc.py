"""生成冻结的体相 Au(111) lead 力常数（M4-T2）。

为什么需要它（方案 §3.5 / D4）：热导链只对**器件区**做有限位移，lead 区的
力常数行拿不到；而 ``thermal_analyze`` 需要 ``D_L``/``V_L`` 才能做表面自能
decimation。lead 区力常数与分子无关，算一次冻结入库，之后每个分子结复用。

为什么是"体相"而不是"板的表面"：decimation 描述的是半无限 lead，器件里那两层
lead 只是它的头两段。用带真空的板会把表面弛豫效应混进本该是体相的量。Pauly 组
的做法相同（"a big cluster of several hundred atoms" 另算一次 + decimation），
器件区与 lead 化学环境不一致是显式接受的近似（方案 §3.5）。

为什么只位移一层：6 层周期胞里，层 0 与层 1 之间的力常数只有一个直接耦合
（另一个方向的近邻在 5 层之外，超出 10 Å 截断），所以不存在 ± 两个像叠加的
双计问题。层 1 的行块由平移对称性给出：``K(1,1) = K(0,0)``、``K(1,0) = K(0,1)^T``。

为什么 k_grid 1 1 1：3×3 超胞 + Γ 等价于原胞 3×3 横向 k 网格；再叠 3×3×1
k 点相当于原胞 9×9，对力常数没有意义，只是把机时乘 4（实测单迭代 53 s → 15 s）。

用法::

    python scripts/gen_au_lead_fc.py --out runs/au-lead-fc --species-dir demo/.../species
    # 跑完后（作业目录里 324/54 个 force 文件齐全）
    python scripts/gen_au_lead_fc.py --assemble runs/au-lead-fc --out src/smj/junction/data/au111_3x3_l4
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / "src"))

from smj.junction.geometry import Atom, write_aims_geometry          # noqa: E402
from smj.junction.tcontrol import make_control                       # noqa: E402
from smj.phonon.fc import assemble_fc, make_displacements, save_lead_fc  # noqa: E402
from smj.stages._thermal_common import THERMAL_OCCUPATION_WIDTH_EV  # noqa: E402
from smj.phonon.forces import read_forces                            # noqa: E402
from smj.phonon.units import ELEMENT_MASS_AMU                        # noqa: E402

A_LATTICE_AU = 4.078230
DELTA_Z_PER_LAYER = A_LATTICE_AU / np.sqrt(3.0)      # 2.35455 Å，与 device.py 同源


def in_plane_vectors(a_lattice: float, in_plane: tuple[int, int]):
    """(a1, a2)：2D 六角胞的两个面内基矢（原子数 = n0 × n1）。"""
    d = a_lattice / np.sqrt(2.0)
    n0, n1 = in_plane
    a1 = np.array([n0 * d, 0.0, 0.0])
    a2 = np.array([n1 * d / 2.0, n1 * d * np.sqrt(3.0) / 2.0, 0.0])
    return a1, a2


def build_bulk_cell(layers: int = 6, in_plane: tuple[int, int] = (3, 3),
                    a_lattice: float = A_LATTICE_AU):
    """Au(111) 沿 z 周期的体相胞：``layers`` 层 × ``in_plane`` 个原子。

    ``layers`` 必须是 3 的倍数，否则 ABC 堆垛在周期边界处错层。
    """
    if layers < 3 or layers % 3:
        raise SystemExit("layers 必须是 3 的倍数（ABC 堆垛周期），得到 %d" % layers)
    n0, n1 = in_plane
    a1, a2 = in_plane_vectors(a_lattice, in_plane)
    shift = (a1 + a2) / (3.0 * n0)          # 面内 ABC 堆垛位移（原胞 (a1+a2)/3）
    atoms = []
    for layer in range(layers):
        z = layer * DELTA_Z_PER_LAYER - (layers - 1) * DELTA_Z_PER_LAYER / 2.0
        offset = shift * (layer % 3)
        for i in range(n0):
            for j in range(n1):
                p = a1 * (i / n0) + a2 * (j / n1) + offset
                atoms.append(Atom("Au", float(p[0]), float(p[1]), float(z)))
    lattice = (tuple(a1), tuple(a2), (0.0, 0.0, layers * DELTA_Z_PER_LAYER))
    return atoms, lattice


def layer_indices(in_plane: tuple[int, int], layers: int, layer: int):
    per = in_plane[0] * in_plane[1]
    return list(range(layer * per, layer * per + per))


def generate(out: Path, *, species_dir: Path, layers: int, in_plane: tuple[int, int],
             ntasks: int, partition: str, mem: str, walltime: str,
             array_parallel: int, module_load: str) -> Path:
    atoms, lattice = build_bulk_cell(layers, in_plane)
    out.mkdir(parents=True, exist_ok=True)
    write_aims_geometry(out / "geometry.in", atoms, lattice=lattice)

    surface = layer_indices(in_plane, layers, 0)
    disps = make_displacements(len(atoms), surface, 0.01)
    # 与 thermal_fc 阶段同格式:裸数组。apply_displacement.load_manifest 只收
    # 数组,写成 {"displacements": [...]} 会让每个数组任务在第一步就退出。
    (out / "displacements.json").write_text(
        json.dumps(disps, indent=2) + "\n", encoding="utf-8")

    make_control(out / "geometry.in", out / "control.in", "thermal_force",
                 str(species_dir), xc="pbe", k_grid=(1, 1, 1),
                 sc_accuracy_forces="5E-3",
                 # 体相 Au 是金属,和器件侧用同一个展宽,免得 lead 的化学环境
                 # 与器件里的那两层因为展宽不同而对不上。
                 occupation_width=THERMAL_OCCUPATION_WIDTH_EV,
                 sc_accuracy_rho="1E-4", charge_mix_param=0.35,
                 n_max_pulay=20)

    # 作业脚本：与热导链数组脚本同构（每任务一个目录 + 汇总哨兵）
    template = (_ROOT / "src" / "smj" / "junction" / "templates"
                / "sub_thermal_fc.sh").read_text(encoding="utf-8")
    script = template.format(
        partition=partition, ntasks=ntasks, mem=mem, walltime=walltime,
        module_load=module_load, exe_define="",
        # gpu1 的节点带 Gres=gpu:4，不显式要 gpu:0 时 Slurm 会把整节点 4 张卡
        # 都算给该作业，于是数组任务之间无法并行（详见 _thermal_common._gres_line）。
        gres_line="--gres=gpu:0",
        aims_bin='mpirun -np %d "%s" > aims.out 2>&1' % (
            ntasks, "/data/apps/fhi-aims.221103-cuda12.8-oneapi2025.3/bin/"
                    "aims.221103.scalapack.mpi.x"),
        last=len(disps) - 1, parallel=array_parallel, count=len(disps))
    with open(out / "sub_lead_fc.sh", "w", encoding="utf-8", newline="\n") as handle:
        handle.write(script)
    for name in ("apply_displacement.py", "extract_forces.py"):
        shutil.copyfile(_ROOT / "src" / "smj" / "phonon" / name, out / name)

    meta = {"layers": layers, "in_plane": list(in_plane), "n_atoms": len(atoms),
            "a_lattice_A": A_LATTICE_AU, "n_displacements": len(disps),
            "surface_layer_indices": surface}
    (out / "lead_case.json").write_text(json.dumps(meta, indent=2) + "\n",
                                        encoding="utf-8")
    return out


def assemble(case: Path, out: Path, *, cutoff: float = 10.0) -> int:
    meta = json.loads((case / "lead_case.json").read_text(encoding="utf-8"))
    in_plane = tuple(meta["in_plane"])
    layers = int(meta["layers"])
    n_atoms = int(meta["n_atoms"])
    per = in_plane[0] * in_plane[1]
    _manifest = json.loads((case / "displacements.json").read_text(encoding="utf-8"))
    disps = _manifest["displacements"] if isinstance(_manifest, dict) else _manifest

    forces = []
    for rec in disps:
        path = case / "fc" / ("%04d" % rec["slot"]) / "forces.json"
        if not path.exists():
            raise SystemExit("数组作业没跑完: %s 不存在" % path)
        payload = json.loads(path.read_text(encoding="utf-8"))
        # extract_forces 写的是 {"n_atoms", "forces", "block_index", "n_blocks"}
        forces.append(np.asarray(payload["forces"], dtype=float))
    full = assemble_fc(disps, forces, n_atoms)      # (3*n_surface, 3*n_atoms)
    rows_l0 = np.arange(3 * per)
    cols_l0 = np.arange(3 * per)
    cols_l1 = np.arange(3 * per, 6 * per)
    k00 = full[np.ix_(rows_l0, cols_l0)]
    k01 = full[np.ix_(rows_l0, cols_l1)]

    # 层 1 的行块由平移对称性给出：K(1,1) = K(0,0)、K(1,0) = K(0,1)^T
    block = np.zeros((6 * per, 6 * per), dtype=float)
    block[:3 * per, :3 * per] = k00
    block[:3 * per, 3 * per:] = k01
    block[3 * per:, :3 * per] = k01.T
    block[3 * per:, 3 * per:] = k00
    block = 0.5 * (block + block.T)                 # 张量对称化

    # 截断：面内按最小像折算，再叠上真实层间距
    positions = np.array([[a.x, a.y, a.z] for a in build_bulk_cell(layers, in_plane)[0]])
    a1, a2 = in_plane_vectors(A_LATTICE_AU, in_plane)
    cell2 = np.array([[a1[0], a1[1]], [a2[0], a2[1]]])
    inverse = np.linalg.inv(cell2)
    dropped = 0
    pair = 2 * per
    for row in range(pair):
        for col in range(pair):
            atom_a = row // 3 + (per if row >= 3 * per else 0)
            atom_b = col // 3 + (per if col >= 3 * per else 0)
            delta = positions[atom_a] - positions[atom_b]
            frac = inverse @ delta[:2]
            wrapped = (frac - np.round(frac)) @ cell2
            distance = float(np.linalg.norm(np.array([wrapped[0], wrapped[1], delta[2]])))
            if distance > cutoff:
                block[3 * row:3 * row + 3, 3 * col:3 * col + 3] = 0.0
                dropped += 1

    save_lead_fc(out, block, element="Au", a_lattice_A=A_LATTICE_AU,
                 n_atoms_per_layer=per, layers=2, in_plane=in_plane,
                 fc_cutoff_A=cutoff,
                 source="gen_au_lead_fc.py: bulk 3x3x%d, Gamma, 0.01 A central differences"
                        % layers)
    residual = float(np.max(np.abs(block.sum(axis=1))))
    print("lead 力常数已写入 %s（%d×%d，截断半径 %.1f Å 丢掉 %d/%d 个原子对，"
          "ASR 行和残差 %.3e eV/Å²）"
          % (out, 3 * pair, 3 * pair, cutoff, dropped, pair ** 2, residual))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True, help="生成目录 / 输出目录")
    ap.add_argument("--assemble", metavar="CASE", help="从已有工况目录组装")
    ap.add_argument("--species-dir", help="FHI-aims species 目录（生成时必需）")
    ap.add_argument("--layers", type=int, default=6)
    ap.add_argument("--in-plane", type=int, nargs=2, default=(3, 3), metavar=("N0", "N1"))
    ap.add_argument("--ntasks", type=int, default=8)
    ap.add_argument("--array-parallel", type=int, default=8)
    ap.add_argument("--partition", default="main")
    ap.add_argument("--mem", default="12G")
    ap.add_argument("--walltime", default="06:00:00")
    ap.add_argument("--module-load", default="module load fhi-aims/221103-cuda12.8-oneapi2025.3")
    ap.add_argument("--cutoff", type=float, default=10.0)
    args = ap.parse_args(argv)

    if args.assemble:
        return assemble(Path(args.assemble), Path(args.out), cutoff=args.cutoff)
    if not args.species_dir:
        ap.error("生成工况时必须给 --species-dir")
    path = generate(Path(args.out), species_dir=Path(args.species_dir),
                    layers=args.layers, in_plane=tuple(args.in_plane),
                    ntasks=args.ntasks, partition=args.partition, mem=args.mem,
                    walltime=args.walltime, array_parallel=args.array_parallel,
                    module_load=args.module_load)
    print("工况已生成: %s" % path)
    print("提交: cd %s && sbatch sub_lead_fc.sh" % path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
