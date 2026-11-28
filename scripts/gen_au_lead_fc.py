"""生成冻结的体相 Au(111) lead 力常数（M4-T2）。

为什么需要它（方案 §3.5 / D4）：热导链只对**器件区**做有限位移，lead 区的
力常数行拿不到；而 ``thermal_analyze`` 需要 ``D_L``/``V_L`` 才能做表面自能
decimation。lead 区力常数与分子无关，算一次冻结入库，之后每个分子结复用。

为什么是"体相"而不是"板的表面"：decimation 描述的是半无限 lead，器件里那两层
lead 只是它的头两段。用带真空的板会把表面弛豫效应混进本该是体相的量。Pauly 组
的做法相同（"a big cluster of several hundred atoms" 另算一次 + decimation），
器件区与 lead 化学环境不一致是显式接受的近似（方案 §3.5）。

分别位移层 0（A）与层 1（B）时，组装先比较实测 AB 与 BA 转置的差异，再对称化。
只有 A 层位移时，沿用 BB=AA、BA=AB 转置的近似，并在资产来源中标明。

为什么 k_grid 1 1 1：3×3 超胞 + Γ 等价于原胞 3×3 横向 k 网格；再叠 3×3×1
k 点相当于原胞 9×9，对力常数没有意义，只是把机时乘 4（实测单迭代 53 s → 15 s）。

用法::

    python scripts/gen_au_lead_fc.py --out lead-a --species-dir demo/.../species
    python scripts/gen_au_lead_fc.py --out lead-b --reference-layer 1 --species-dir demo/.../species
    # 两套位移作业均有完整力文件后
    python scripts/gen_au_lead_fc.py --assemble lead-a lead-b --out lead-asset
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / "src"))

from smj.junction.geometry import Atom, read_aims_geometry, write_aims_geometry  # noqa: E402
from smj.junction.tcontrol import make_control                       # noqa: E402
from smj.phonon.fc import assemble_fc, make_displacements, save_lead_fc  # noqa: E402
from smj.phonon.geometry import atom_distances                        # noqa: E402
from smj.phonon.concat_forces import scan_aims                         # noqa: E402
from smj.stages.thermal_fc import ThermalFcStage                      # noqa: E402
from smj.stages._thermal_common import THERMAL_OCCUPATION_WIDTH_EV  # noqa: E402
from smj.phonon.forces import read_forces                            # noqa: E402

A_LATTICE_AU = 4.078230

# Both the login-node submit command and each array worker use the same check.
_LEAD_JOB_DIR = r"""
if [[ ${1:-} != --submit ]]; then
    cd -- "${SLURM_SUBMIT_DIR:?missing Slurm submit directory}" || exit 1
fi
_home=$(cd -P -- "${HOME:?missing home directory}" && pwd -P) || exit 1
_stage=$(pwd -P) || exit 1
_parent=${PWD%/*}
_run=${_parent##*/}
_rel="smj-runs/$_run/lead_fc"
check_lead_dir() {
    local LC_ALL=C
    [[ $_run =~ ^[A-Za-z0-9]+(-[A-Za-z0-9]+)*$ && ${#_run} -le 40 &&
       $_stage == "$_home/$_rel" && ( $PWD == "$HOME/$_rel" || $PWD == "$_stage" ) ]]
}
if ! check_lead_dir; then
    echo 'Run this job from $HOME/smj-runs/<run-name>/lead_fc.' >&2
    exit 1
fi
# An old linked output or slot could redirect the calculation outside the stage.
_link=$(find . -ignore_readdir_race -type l -print -quit) || exit 1
if [[ -n $_link ]]; then
    echo "Refusing symlink in lead_fc: $_link" >&2
    exit 1
fi
if [[ ${1:-} == --submit ]]; then
    if [[ $# -ne 1 ]]; then
        echo 'Usage: bash sub_lead_fc.sh --submit' >&2
        exit 1
    fi
    exec sbatch --chdir="$_stage" --output='slurm-%A_%a.out' --error='slurm-%A_%a.out' sub_lead_fc.sh
fi
"""


def _geometry(path: Path):
    atoms = read_aims_geometry(path)
    if any(atom.element != "Au" for atom in atoms):
        raise ValueError(f"{path}: expected only Au atoms")
    vectors = [line.split()[1:4] for line in path.read_text(encoding="utf-8").splitlines()
               if line.split()[:1] == ["lattice_vector"]]
    return np.asarray(vectors, dtype=float), np.array([atom.coord() for atom in atoms])


def _read_reference(case: Path):
    meta = json.loads((case / "lead_case.json").read_text(encoding="utf-8"))
    in_plane = meta["in_plane"]
    n_atoms = meta["n_atoms"]
    layers = meta["layers"]
    if (len(in_plane) != 2 or any(type(v) is not int or v < 1 for v in in_plane)
            or type(layers) is not int or layers < 3 or layers % 3
            or n_atoms != layers * math.prod(in_plane)
            or not math.isfinite(meta["a_lattice_A"]) or meta["a_lattice_A"] <= 0):
        raise ValueError(f"{case}: invalid Au supercell metadata")
    lattice, positions = _geometry(case / "geometry.in")
    if (positions.shape != (n_atoms, 3) or lattice.shape != (3, 3)
            or not np.isfinite(positions).all() or not np.isfinite(lattice).all()):
        raise ValueError(f"{case}: invalid geometry.in")

    records = json.loads((case / "displacements.json").read_text(encoding="utf-8"))
    if isinstance(records, dict):
        records = records["displacements"]
    per = math.prod(in_plane)
    if len(records) != 6 * per or meta["n_displacements"] != len(records):
        raise ValueError(f"{case}: incomplete displacement manifest")
    atoms = sorted({record["atom"] for record in records})
    layer = atoms[0] // per
    if (layer not in (0, 1) or atoms != layer_indices(in_plane, layers, layer)
            or sorted(meta["surface_layer_indices"]) != atoms
            or meta.get("reference_layer", layer) != layer):
        raise ValueError(f"{case}: displaced atoms must match reference layer 0 or 1")
    slots = [rec["slot"] for rec in records]
    if sorted(slots) != list(range(len(records))):
        raise ValueError(f"{case}: invalid or duplicated displacement slots")
    records = sorted(records, key=lambda rec: rec["slot"])
    atoms = sorted({record["atom"] for record in records})
    step = abs(records[0]["delta_A"])
    if not math.isfinite(step) or step == 0.0:
        raise ValueError(f"{case}: non-finite or zero displacement step")
    bundle_path = case / "forces_all.json"
    summary_path = case / "fc_summary.json"
    if bundle_path.is_file():
        force_blocks, summary = ThermalFcStage._read_force_data(
            case, len(records), n_atoms)
        source = "forces_all.json"
        evidence = "converged (summary)" if summary is not None else "unknown (bundle without summary)"
    else:
        force_blocks = []
        summary = ThermalFcStage._read_json(summary_path) if summary_path.is_file() else None
        if summary is not None:
            ThermalFcStage._check_scf(summary, len(records))
            evidence = "converged (summary)"
        else:
            evidence = "unknown"
        sources = set()
        log_status = []
        for rec in records:
            slot_dir = case / "fc" / ("%04d" % rec["slot"])
            forces_path = slot_dir / "forces.json"
            log_path = slot_dir / "aims.out"
            if forces_path.is_file():
                payload = json.loads(forces_path.read_text(encoding="utf-8"))
                if payload.get("n_atoms", n_atoms) != n_atoms:
                    raise ValueError(f"{forces_path}: atom count does not match")
                force = np.asarray(payload["forces"], dtype=float)
                sources.add("forces.json")
            elif log_path.is_file():
                force = read_forces(log_path, n_atoms=n_atoms)
                sources.add("aims.out")
            else:
                raise FileNotFoundError(f"{slot_dir}: missing forces.json and aims.out")
            force_blocks.append(force)
            if summary is None and log_path.is_file():
                log_status.append(scan_aims(
                    log_path.read_text(encoding="utf-8", errors="replace"))["converged"])
        if summary is None and log_status:
            if False in log_status:
                raise ValueError(f"{case}: explicit unconverged SCF evidence in aims.out")
            if len(log_status) == len(records) and all(status is True for status in log_status):
                evidence = "converged (aims.out)"
        source = "+".join(sorted(sources))
        if len(sources) > 1:
            source = "mixed " + source
    forces = []
    for rec, force in zip(records, force_blocks):
        if (type(rec["atom"]) is not int or type(rec["axis"]) is not int
                or rec["axis"] not in (0, 1, 2) or type(rec["slot"]) is not int
                or rec["sign"] not in (-1, 1) or rec["delta_A"] != rec["sign"] * step):
            raise ValueError(f"{case}: invalid displacement at slot {rec['slot']}")
        if force.shape != (n_atoms, 3) or not np.isfinite(force).all():
            raise ValueError(f"{case}: invalid force data at slot {rec['slot']}")
        forces.append(force)
    # assemble_fc enumerates the measured atoms in sorted order. The equality
    # above establishes the actual row IDs even if slots are interleaved.
    fc = assemble_fc(records, forces, n_atoms)
    return layer, fc, meta, step, positions, lattice, source, evidence


def _directional_asr(matrix):
    """逐笛卡尔平移方向检查实测矩阵行，可接受未测全的矩形矩阵。"""
    sums = np.asarray(matrix).reshape(matrix.shape[0], -1, 3).sum(axis=1)
    return float(np.max(np.abs(sums)))


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
    if n0 < 1 or n1 < 1 or not math.isfinite(a_lattice) or a_lattice <= 0:
        raise ValueError("面内尺寸和晶格常数必须为正")
    a1, a2 = in_plane_vectors(a_lattice, in_plane)
    shift = a1 / (3.0 * n0) + a2 / (3.0 * n1)
    dz = a_lattice / np.sqrt(3.0)
    atoms = []
    for layer in range(layers):
        z = layer * dz - (layers - 1) * dz / 2.0
        offset = shift * (layer % 3)
        for i in range(n0):
            for j in range(n1):
                p = a1 * (i / n0) + a2 * (j / n1) + offset
                atoms.append(Atom("Au", float(p[0]), float(p[1]), float(z)))
    lattice = (tuple(a1), tuple(a2), (0.0, 0.0, layers * dz))
    return atoms, lattice


def layer_indices(in_plane: tuple[int, int], layers: int, layer: int):
    per = in_plane[0] * in_plane[1]
    return list(range(layer * per, layer * per + per))


def generate(out: Path, *, species_dir: Path, layers: int, in_plane: tuple[int, int],
             ntasks: int, partition: str, mem: str, walltime: str,
             array_parallel: int, module_load: str, reference_layer: int = 0) -> Path:
    atoms, lattice = build_bulk_cell(layers, in_plane)
    out.mkdir(parents=True, exist_ok=True)
    write_aims_geometry(out / "geometry.in", atoms, lattice=lattice)

    if reference_layer not in (0, 1) or reference_layer >= layers:
        raise ValueError("reference_layer 只能为 0 或 1")
    surface = layer_indices(in_plane, layers, reference_layer)
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
        module_load=_LEAD_JOB_DIR + module_load, exe_define="",
        # gpu1 的节点带 Gres=gpu:4，不显式要 gpu:0 时 Slurm 会把整节点 4 张卡
        # 都算给该作业，于是数组任务之间无法并行（详见 _thermal_common._gres_line）。
        gres_line="--gres=gpu:0",
        aims_bin='mpirun -np %d "%s" > aims.out 2>&1' % (
            ntasks, "/data/apps/fhi-aims.221103-cuda12.8-oneapi2025.3/bin/"
                    "aims.221103.scalapack.mpi.x"),
        last=len(disps) - 1, parallel=array_parallel, count=len(disps))
    with open(out / "sub_lead_fc.sh", "w", encoding="utf-8", newline="\n") as handle:
        handle.write(script)
    for name in ("apply_displacement.py", "extract_forces.py", "concat_forces.py"):
        shutil.copyfile(_ROOT / "src" / "smj" / "phonon" / name, out / name)

    meta = {"layers": layers, "in_plane": list(in_plane), "n_atoms": len(atoms),
            "a_lattice_A": A_LATTICE_AU, "n_displacements": len(disps),
            "surface_layer_indices": surface, "reference_layer": reference_layer}
    (out / "lead_case.json").write_text(json.dumps(meta, indent=2) + "\n",
                                        encoding="utf-8")
    return out


def assemble(case: Path | list[Path], out: Path, *, cutoff: float = 10.0) -> int:
    cases = [case] if isinstance(case, (Path, str)) else case
    if len(cases) not in (1, 2) or not math.isfinite(cutoff) or cutoff <= 0.0:
        raise ValueError("Provide one or two reference layers and a positive finite cutoff")
    measured = {}
    sources = {}
    force_inputs = {}
    for path in map(Path, cases):
        layer, rows, meta, step, positions, lattice, force_source, scf_evidence = _read_reference(path)
        if layer in measured:
            raise ValueError("Duplicate reference layer %d" % layer)
        settings = (meta["layers"], meta["in_plane"], meta["a_lattice_A"], step)
        if measured:
            if (settings != first_settings or not np.array_equal(positions, first_positions)
                    or not np.array_equal(lattice, first_lattice)):
                raise ValueError("Reference cases have different supercells or displacement settings")
            control = (path / "control.in").read_text(encoding="utf-8")
            if control != (Path(cases[0]) / "control.in").read_text(encoding="utf-8"):
                raise ValueError("Reference cases have different control.in settings")
        else:
            first_settings, first_positions, first_lattice = settings, positions, lattice
        measured[layer] = rows
        sources[layer] = str(path)
        force_inputs[layer] = {"source": force_source, "scf_evidence": scf_evidence}
    if 0 not in measured:
        raise ValueError("A layer-1 case requires its matching layer-0 case")

    layers, in_plane, a_lattice, step = first_settings
    per = math.prod(in_plane)
    split = 3 * per
    aa, ab = measured[0][:, :split], measured[0][:, split:2 * split]
    intra_reciprocity = None
    inter_reciprocity = None
    if 1 in measured:
        ba, bb = measured[1][:, :split], measured[1][:, split:2 * split]
        intra_a = float(np.max(np.abs(aa - aa.T)))
        intra_b = float(np.max(np.abs(bb - bb.T)))
        intra_reciprocity = max(intra_a, intra_b)
        inter_reciprocity = float(np.max(np.abs(ab - ba.T)))
        method = "measured A+B rows, reciprocal pairs averaged"
    else:
        ba, bb = ab.T, aa
        method = "A-only approximation (BB=AA, BA=AB.T)"
    raw = np.block([[aa, ab], [ba, bb]])
    block = 0.5 * (raw + raw.T)
    pair = 2 * per
    drop = atom_distances(first_positions[:pair], first_lattice) > cutoff
    block[drop.repeat(3, axis=0).repeat(3, axis=1)] = 0.0

    source = (f"gen_au_lead_fc.py: {in_plane[0]}x{in_plane[1]}x{layers} bulk, "
              f"{step:g} A central differences, {method}")
    save_lead_fc(out, block, element="Au", a_lattice_A=a_lattice,
                 n_atoms_per_layer=per, layers=2, in_plane=in_plane,
                 fc_cutoff_A=cutoff, source=source)
    full_rows = np.vstack(list(measured.values()))
    residual = _directional_asr(full_rows)
    relative = residual / (float(np.max(np.abs(full_rows))) or 1.0)
    details = {
        "measured_layers": sorted(measured), "source_cases_by_layer": sources,
        "force_inputs_by_layer": force_inputs,
        "raw_intralayer_reciprocity_max_ev_per_A2": intra_reciprocity,
        "raw_intralayer_AA_reciprocity_max_ev_per_A2": intra_a if 1 in measured else None,
        "raw_intralayer_BB_reciprocity_max_ev_per_A2": intra_b if 1 in measured else None,
        "raw_interlayer_reciprocity_max_ev_per_A2": inter_reciprocity,
        "raw_reciprocity_max_ev_per_A2": max(
            value for value in (intra_reciprocity, inter_reciprocity)
            if value is not None) if 1 in measured else None,
        "asr_residual_full_cell_ev_per_A2": residual,
        "asr_relative_full_cell": relative,
        "asr_residual_two_layer_block_ev_per_A2": _directional_asr(block),
        "n_layers": layers,
        "note": "ASR uses each measured row with all bulk atoms before truncation and symmetrization. "
                "Raw intralayer reciprocity is max abs(AA-AA.T, BB-BB.T); raw interlayer "
                "reciprocity is max abs(AB-BA.T).",
    }
    (out / "lead_asr.json").write_text(
        json.dumps(details, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"lead 力常数已写入 {out}（{method}，截断半径 {cutoff:g} Å "
          f"丢掉 {np.count_nonzero(drop)}/{pair ** 2} 个原子对）")
    print(f"  实测行对完整 {layers} 层周期胞的逐方向 ASR 残差 {residual:.3e} eV/Å²")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True, help="生成目录 / 输出目录")
    ap.add_argument("--assemble", metavar="CASE", nargs="+", help="从一个或两个参考层工况组装")
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
    ap.add_argument("--reference-layer", type=int, choices=(0, 1), default=0,
                    help="生成第 0 层或第 1 层的独立有限位移")
    args = ap.parse_args(argv)

    if args.assemble:
        return assemble([Path(case) for case in args.assemble],
                        Path(args.out), cutoff=args.cutoff)
    if not args.species_dir:
        ap.error("生成工况时必须给 --species-dir")
    path = generate(Path(args.out), species_dir=Path(args.species_dir),
                    layers=args.layers, in_plane=tuple(args.in_plane),
                    ntasks=args.ntasks, partition=args.partition, mem=args.mem,
                    walltime=args.walltime, array_parallel=args.array_parallel,
                    module_load=args.module_load, reference_layer=args.reference_layer)
    print("工况已生成: %s" % path)
    example = "au-lead-a" if args.reference_layer == 0 else "au-lead-b"
    print("将工况文件复制到集群 ~/smj-runs/<run-name>/lead_fc/。示例提交命令：\n"
          f'  cd "$HOME/smj-runs/{example}/lead_fc" && bash sub_lead_fc.sh --submit')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
