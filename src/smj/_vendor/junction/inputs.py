"""单分子结程序族输入写出器(AGD-2026-09-22-SMJ P4b)。

来源: 第三方工具 qmtransport/gaussian.py 的 make_gaussian_input + cli.py 的
add-tip-au/append-au60 语义, 按本项目协议适配(模板自 tools/junction/templates/,
env 显式传入, 缺省从环境变量读取)。四个集群阶段:
- make_gaussian_input:  分子 XYZ -> Gaussian opt+freq gjf + sub_g16.sh
- make_tip_relax_input: 消费 opt.xyz, add_tip_au 加 2 Au tip -> geometry.in +
                        control_relax.in(bfgs) + sub_aims.sh
- make_aims_scf_input:  消费 au60_build geometry.in -> control_scf.in + sub_aims.sh
- make_aitranss_input:  消费 aims_scf 产物(omat/mos/basis/aims.out) -> tcontrol + sub_te.sh

纪律(AGENTS.md 来源规则 / 计划 §4.3):
- species blocks 只来自用户提供目录(find_species_file), 绝不虚构;
- 虚频/表面原子不静默: 虚频交给 validate 报错, 表面 6 索引需人工复核;
- 集群部署路径(g16root/aims_prefix/aitranss_bin)从 env 参数取, 缺省回退环境变量,
  再缺省写 UNCONFIGURED 占位(提交前 validate_junction_run 拦截, 不静默伪造);
- 每阶段记录参数清单进 summary JSON(charge/spin/method/basis/species tier/partition/cores)。
纯本地输入生成, 不做真实计算。
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path

from .geometry import (Atom, add_tip_au, read_aims_geometry, read_xyz,
                       strip_thiol_h as _strip_thiol_h, write_aims_geometry, write_xyz)
from .tcontrol import find_species_file, make_control, make_tcontrol, read_template

#: DEFAULT_ENV 的集群变量**单源清单**: (DEFAULT_ENV 键, 环境变量名, 缺省值)。
#: 2026-10-24(缺口⑦ 教训): 原先这里是手写字典 —— 会话脚本新增 AIMS_BIN / AIMS_LD_MODULE /
#: MKL_LIB 时没人同步映射, 真实 run 的值被静默丢成 UNCONFIGURED/aims.x(slurm 秒退)。
#: 现改为「声明式表 + 派生 DEFAULT_ENV」, 并由 tests/test_smj_aims_binary_env.py::
#: test_session_script_env_vars_are_all_mapped 对账 agenticdft/scripts/smj-cluster-env-session.cmd
#: 的每个非空 `set NAME=`: 脚本新增变量却忘了映射 -> 测试立刻红, 而不是等真跑失败。
ENV_SPECS = (
    # Gaussian
    ('g16root', 'G16ROOT', ''),
    ('g16_scratch', 'G16_SCRATCH', ''),
    ('g16_module', 'G16_MODULE', ''),          # module load 名(如 Gaussian/g16.c02)
    ('g16_mem', 'G16_MEM', '16G'),             # 单作业内存申请(避免贪整节点导致排队)
    ('g16_walltime', 'G16_WALLTIME', '01:00:00'),   # 时长上限(短作业利于 backfill)
    ('g16_gres', 'G16_GRES', ''),              # GPU 请求(如 'gpu:1'; 空=纯 CPU)
    # FHI-aims
    ('aims_prefix', 'AIMS_PREFIX', ''),
    ('aims_module', 'AIMS_MODULE', ''),        # module 名(遮蔽 aims_bin; 见 junction/preflight.py)
    ('aims_bin', 'AIMS_BIN', ''),              # 直接可执行文件(scalapack/mpi 单文件版)
    ('aims_ld_module', 'AIMS_LD_MODULE', ''),  # 只提供运行时库的 module(scalapack 版链 libcublas/libmkl/libimf)
    ('aims_gpus', 'AIMS_GPUS', ''),            # GPU 数(非空=GPU 模式, 1 GPU 对应 1 MPI rank)
    ('aims_mem', 'AIMS_MEM', '16G'),
    ('aims_walltime', 'AIMS_WALLTIME', '02:00:00'),
    # aitranss
    ('aitranss_bin', 'AITRANSS_BIN', ''),
    ('mkl_lib', 'MKL_LIB', ''),                # aitranss 运行时 LD_LIBRARY_PATH(oneAPI mkl 库目录)
    ('te_mem', 'TE_MEM', '8G'),
    ('te_gres', 'TE_GRES', ''),
    ('te_walltime', 'TE_WALLTIME', '02:00:00'),
)

#: 会话脚本里出现但**有意不映射**进 DEFAULT_ENV 的环境变量 -> 理由。
#: 目标态是空表; 每加一条必须写清「值从哪条路径显式传入」, 否则就是缺口⑦ 类静默丢配置。
UNMAPPED_ENV_ALLOWED = {
    'AIMS_SPECIES_DIR': 'species 目录由 CLI --aims-species-tier 显式传入, 不属集群部署路径映射',
}

#: ENV_SPECS 里登记的环境变量名(测试对账用; 唯一来源仍是 ENV_SPECS)
MAPPED_ENV_VARS = frozenset(env for _, env, _ in ENV_SPECS)

#: 默认 env = 环境变量快照(import 时) + 非环境变量缺省(调度资源类)。
DEFAULT_ENV = {key: os.environ.get(env, default) for key, env, default in ENV_SPECS}
DEFAULT_ENV.update({'partition': 'main', 'ntasks': 1, 'nproc': 64})


def _env(env: dict | None) -> dict:
    merged = dict(DEFAULT_ENV)
    if env:
        merged.update({k: v for k, v in env.items() if v is not None})
    return merged


def _require(env: dict, keys: list[str]) -> None:
    missing = [k for k in keys if not str(env.get(k) or '').strip()]
    if missing:
        raise ValueError('smj 集群部署路径未配置(环境变量 %s 未设): %s —— 提交前须配置, 不静默伪造'
                         % (', '.join(k.upper() for k in missing), ', '.join(missing)))


def _write(path: Path, text: str) -> None:
    """写文件, 强制 LF 换行(集群作业脚本/gjf 必须 Unix 换行, CRLF 会被 sbatch/Gaussian 拒绝)。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.replace('\r\n', '\n').replace('\r', '\n'), encoding='utf-8', newline='\n')


def make_gaussian_input(xyz_path: str | Path, outdir: str | Path, name: str,
                        charge: int = 0, multiplicity: int = 1,
                        method: str = 'b3lyp', basis: str = '6-31g(d)',
                        env: dict | None = None) -> dict:
    """分子 XYZ -> Gaussian opt+freq 输入(gjf + sub_g16.sh)。"""
    e = _env(env)
    atoms = read_xyz(xyz_path)
    coordinates = '\n'.join(' %-2s %16.8f %16.8f %16.8f' % (a.element, a.x, a.y, a.z)
                            for a in atoms)
    out = Path(outdir)
    gjf = out / ('%s.gjf' % name)
    _write(gjf,
        read_template('gaussian_opt_freq.gjf').format(
            name=name, nproc=e['nproc'], method=method, basis=basis,
            charge=charge, multiplicity=multiplicity, coordinates=coordinates))
    sub = out / 'sub_g16.sh'
    # 高斯加载方式: env.g16_module(如 Gaussian/g16.c02)非空 -> module load;
    # 否则回退 g16root 路径方式(兼容既有测试/无 module 系统环境)。
    if str(e.get('g16_module') or '').strip():
        g16_setup = 'module load %s' % str(e['g16_module']).strip()
    else:
        g16_setup = (
            'export COM_DIR={g16root}\n'
            'export g16root=$COM_DIR\n'
            'export PATH=$g16root:$PATH\n'
            'source $g16root/bsd/g16.profile\n'
            'export GAUSS_EXEDIR=$g16root'
        ).format(g16root=e['g16root'] or 'UNCONFIGURED')
    _write(sub,
        read_template('sub_g16.sh').format(
            partition=e['partition'], nproc=e['nproc'],
            mem=e.get('g16_mem') or '16G',
            walltime=e.get('g16_walltime') or '01:00:00',
            gres_line=('#SBATCH --gres=%s' % str(e['g16_gres']).strip())
            if str(e.get('g16_gres') or '').strip() else '',
            g16_setup=g16_setup,
            scratch=e['g16_scratch'] or 'UNCONFIGURED'))
    summary = {'name': name, 'method': method, 'basis': basis, 'charge': charge,
               'multiplicity': multiplicity, 'nproc': e['nproc'],
               'partition': e['partition'], 'atom_count': len(atoms),
               'gaussian_load': 'module:%s' % str(e['g16_module']).strip()
               if str(e.get('g16_module') or '').strip() else 'g16root:%s' % (e['g16root'] or 'UNCONFIGURED')}
    _write(out / 'input_summary.json', json.dumps(summary, indent=2, ensure_ascii=False))
    return {'gjf': str(gjf), 'submit_script': str(sub), 'summary': summary}


def _render_aims_sub(e: dict, ntasks: int) -> str:
    """渲染 sub_aims.sh: 双模式(CUDA GPU / 纯 CPU)。
    GPU 模式(env.aims_gpus 非空): 1 GPU 对应 1 MPI rank(1 核心), ntasks==gpus 数;
    CPU 模式(aims_gpus 空): ntasks 核。module load 优先 aims_module, 回退 aims_bin/aims_prefix。
    直接二进制(aims_bin)分支: 若同时给 aims_ld_module 则仍插 module load(只提供运行时库,
    不选二进制; 实测该 scalapack 版链 libcublas/libmkl/libimf, 不 load 会 execvp/缺库失败)。
    aims_gres_sched(仅调度): 加 #SBATCH --gres=gpu:N 但不改 rank 数(占空闲 GPU 快速调度, 计算仍 CPU)。"""
    gpus = str(e.get('aims_gpus') or '').strip()
    gres_sched = str(e.get('aims_gres_sched') or '').strip()
    if gpus:
        nrank = int(gpus)
        gres_line = '#SBATCH --gpus=%s' % gpus
    else:
        nrank = int(ntasks)
        gres_line = ('#SBATCH --gres=%s' % gres_sched) if gres_sched else ''
    if str(e.get('aims_module') or '').strip():
        module_load = 'module load %s' % str(e['aims_module']).strip()
        exe_define = 'AIMS_EXE="${FHI_AIMS_BIN:?FHI-aims module did not define FHI_AIMS_BIN}"'
        aims_bin = 'mpirun -np %d "$AIMS_EXE" > aims.out 2>&1' % nrank
    elif str(e.get('aims_bin') or '').strip():
        # 直接指定 FHI-aims 可执行文件路径(如 scalapack/mpi 单文件版)。
        # 2026-10-24 缺口⑦: 该二进制仍链 CUDA/MKL/Intel 运行时库, 只给路径不 load module
        # 会 execvp/缺库失败 -> aims_ld_module 只负责提供 LD_LIBRARY_PATH, 不改变所选二进制。
        ldmod = str(e.get('aims_ld_module') or '').strip()
        module_load = ('module load %s' % ldmod) if ldmod else ''
        exe_define = ''
        aims_bin = 'mpirun -np %d "%s" > aims.out 2>&1' % (nrank, str(e['aims_bin']).strip())
    else:
        module_load = ('export PATH={p}/install/bin:$PATH\n'
                       'export LD_LIBRARY_PATH={p}/install/lib:$LD_LIBRARY_PATH\n'
                       'export INCLUDE={p}/install/include:$INCLUDE').format(p=e['aims_prefix'] or 'UNCONFIGURED')
        exe_define = ''
        aims_bin = 'mpirun -np %d aims.x' % nrank
    return read_template('sub_aims.sh').format(
        partition=e['partition'], ntasks=ntasks,
        mem=e.get('aims_mem') or '16G',
        gres_line=gres_line,
        walltime=e.get('aims_walltime') or '02:00:00',
        module_load=module_load,
        exe_define=exe_define,
        aims_bin=aims_bin)


def measure_au_s_bonds(geometry_path: str | Path, left_anchor: int, right_anchor: int,
                       au_tip_count: int = 2) -> dict:
    """量 geometry.in 里两条锚定-S ↔ Au tip 键长(A)。tip 约定为**文件末尾** `au_tip_count` 个原子
    (add_tip_au 的写入顺序: 分子原子在前, 左 tip、右 tip 依次追加)。

    用于两处: ① 建结时自校验"落盘的距离就是目标距离"; ② 弛豫后量**实际**Au–S,
    不让"目标 2.45 Å"停留在写输入那一刻(H1 归因要的是弛豫后的真实值)。
    """
    atoms = read_aims_geometry(geometry_path)
    n = len(atoms)
    if au_tip_count < 1 or n <= au_tip_count:
        raise ValueError('geometry 原子数 %d 不足以识别 %d 个 Au tip' % (n, au_tip_count))
    for idx in (left_anchor, right_anchor):
        if not 1 <= int(idx) <= n - au_tip_count:
            raise ValueError('锚定原子序号 %r 落在 tip 区(分子只有 %d 个原子): 拒绝静默量错原子'
                             % (idx, n - au_tip_count))
    tips = [n - au_tip_count + 1, n]  # 1-based: 左 tip、右 tip
    return {'left': math.dist(atoms[int(left_anchor) - 1].coord(), atoms[tips[0] - 1].coord()),
            'right': math.dist(atoms[int(right_anchor) - 1].coord(), atoms[tips[1] - 1].coord())}


def make_tip_relax_input(xyz_path: str | Path, outdir: str | Path,
                         left_anchor: int, right_anchor: int, au_tip_distance: float,
                         left_direction: tuple | None = None,
                         right_direction: tuple | None = None,
                         species_dir: str | Path | None = None,
                         xc: str = 'pbe', spin: str = 'none', charge: float = 0.0,
                         env: dict | None = None,
                         strip_thiol_h: bool = False,
                         anchor_au_s_target: float | None = None,
                         pin_contact: bool = True,
                         vdw: str | None = None) -> dict:
    """消费 opt.xyz, 加 2 Au tip(锚定原子/距离/方向), 写 geometry.in +
    control_relax.in(bfgs) + sub_aims.sh。species blocks 只来自用户提供目录。

    strip_thiol_h=True 时先移除与 S 成键的 H(二硫醇盐, TaskCard §1.4.1 路线), 锚定
    索引随之重映射; 去质子前几何另存 opt_neutral.xyz 供复核。缺省 False = 原样保留。

    W0-1(anchor_au_s_target): 给定目标 Au–S 键长(埃)时, tip 直接放到该距离;
    `pin_contact=True`(缺省)再把两侧锚定原子与两个 tip 一起冻结
    (`constrain_relaxation .true.`), 使 bfgs 无法把接触压短 —— 这正是"距离扫描结"的
    口径(分子刚性 + 电极间距可设)。`pin_contact=False` 则只改初始距离、仍全自由弛豫,
    用来回答"光把 tip 放远, 弛豫会不会自己压回去"。
    `anchor_au_s_target=None` 时完全不进该分支: geometry.in 与旧版**逐字节一致**。
    """
    e = _env(env)
    if not species_dir:
        raise ValueError('aims_tip_relax 需 species_dir(FHI-aims species 目录, 用户提供; 不虚构)')
    target = None
    if anchor_au_s_target is not None:
        target = float(anchor_au_s_target)
        if not 1.5 <= target <= 3.5:
            raise ValueError('anchor_au_s_target 需在 1.5–3.5 A 之间(实得 %.3f): Au–S 键长约 '
                             '2.3–2.5 A, 超出即 fail-closed, 不静默照建' % target)
        distance = target
    elif au_tip_distance is None:
        raise ValueError('aims_tip_relax 需 au_tip_distance 或 anchor_au_s_target(二选一)')
    else:
        distance = float(au_tip_distance)
    atoms = read_xyz(xyz_path)
    raw_atoms = atoms
    strip_info = None
    if strip_thiol_h:
        # 去 H 后原子索引前移 -> 锚定索引必须同步重映射(不静默沿用旧索引)
        atoms, strip_info = _strip_thiol_h(atoms, anchors=(left_anchor, right_anchor))
        left_anchor = strip_info['left_anchor']
        right_anchor = strip_info['right_anchor']
    junction = add_tip_au(atoms, left_anchor, right_anchor, distance,
                          left_direction, right_direction)
    n_mol = len(atoms)
    constrained: tuple[int, ...] = ()
    if target is not None and pin_contact:
        # 分子两个锚定原子 + 两个 tip: 冻结后 bfgs 改不动这条键, 实际键长 = 目标键长
        constrained = (int(left_anchor), int(right_anchor), n_mol + 1, n_mol + 2)
    out = Path(outdir)
    geom = out / 'geometry.in'
    write_aims_geometry(geom, junction, constrained=constrained)
    placed = measure_au_s_bonds(geom, left_anchor, right_anchor)
    if target is not None:
        for side in ('left', 'right'):
            if abs(placed[side] - target) > 1e-6:
                raise RuntimeError('Au–S 落盘距离与目标不符(%s: %.6f vs %.6f A): 建结逻辑有 bug, '
                                   '拒绝把这个几何交给 aims 弛豫' % (side, placed[side], target))
    if strip_info is not None:
        neutral = out / 'opt_neutral.xyz'
        write_xyz(neutral, raw_atoms,
                  comment='Gaussian opt 原样(含 S-H), 去质子前存档')
        strip_info['opt_neutral_xyz'] = str(neutral)
    make_control(geom, out / 'control.in', 'relax', species_dir,
                 xc=xc, spin=spin, charge=charge, vdw=vdw)
    sub = out / 'sub_aims.sh'
    _write(sub, _render_aims_sub(e, e['ntasks']))
    summary = {'left_anchor': left_anchor, 'right_anchor': right_anchor,
               'au_tip_distance_A': au_tip_distance,
               'placement_distance_A': distance,
               'anchor_au_s_target': target,
               'contact_pinned': bool(constrained),
               'placed_au_s_A': placed if target is not None else None,
               'left_direction': list(left_direction) if left_direction else 'auto',
               'right_direction': list(right_direction) if right_direction else 'auto',
               'species_dir': str(species_dir), 'xc': xc, 'spin': spin, 'charge': charge,
               'ntasks': e['ntasks'], 'partition': e['partition'],
               'aims_gpus': str(e.get('aims_gpus') or '').strip() or 'cpu',
               'atom_count': len(junction), 'au_tip_count': 2,
               'vdw': str(vdw or ''),
               'strip_thiol_h': strip_info}

    _write(out / 'input_summary.json', json.dumps(summary, indent=2, ensure_ascii=False))
    return {'geometry': str(geom), 'control': str(out / 'control.in'),
            'submit_script': str(sub), 'summary': summary}


def make_aims_scf_input(geometry_in: str | Path, outdir: str | Path,
                        species_dir: str | Path | None = None,
                        xc: str = 'pbe', spin: str = 'none', charge: float = 0.0,
                        fixed_spin_moment: float = 0.0, occupation_width: float = 0.01,
                        sc_iter_limit: int = 500, env: dict | None = None,
                        vdw: str | None = None) -> dict:
    """消费 au60_build 的 geometry.in -> control_scf.in(aitranss 输出) + sub_aims.sh。"""
    e = _env(env)
    if not species_dir:
        raise ValueError('aims_scf 需 species_dir(FHI-aims species 目录, 用户提供; 不虚构)')
    geom = Path(geometry_in)
    if not geom.exists():
        raise FileNotFoundError('aims_scf 缺上游 geometry.in: %s' % geom)
    out = Path(outdir)
    make_control(geom, out / 'control.in', 'scf', species_dir,
                 xc=xc, spin=spin, charge=charge, occupation_width=occupation_width,
                 fixed_spin_moment=fixed_spin_moment, sc_iter_limit=sc_iter_limit, vdw=vdw)
    sub = out / 'sub_aims.sh'
    _write(sub, _render_aims_sub(e, e['ntasks']))
    summary = {'species_dir': str(species_dir), 'xc': xc, 'spin': spin, 'charge': charge,
               'fixed_spin_moment': fixed_spin_moment, 'occupation_width': occupation_width,
               'sc_iter_limit': sc_iter_limit, 'ntasks': e['ntasks'],
               'partition': e['partition'],
               'aims_gpus': str(e.get('aims_gpus') or '').strip() or 'cpu',
               'vdw': str(vdw or ''),
               'atom_count': len(read_aims_geometry(geom))}
    _write(out / 'input_summary.json', json.dumps(summary, indent=2, ensure_ascii=False))
    return {'control': str(out / 'control.in'), 'submit_script': str(sub), 'summary': summary}


def make_aitranss_input(scf_dir: str | Path, outdir: str | Path,
                        metadata: str | Path | None = None,
                        nsaos: int | None = None, efermi: float | None = None,
                        surfaces: list[int] | None = None,
                        nlayers: int = 4, ener: float = -0.4, estep: float = 0.0001,
                        eend: float = 0.0, valence_electrons: int = 0,
                        env: dict | None = None) -> dict:
    """消费 aims_scf 产物 -> tcontrol + sub_te.sh(aitranss 透射)。
    efermi 解析 TE.dat 表头优先(§4.3), 缺失回退 aims.out, 均缺 fail-closed。
    表面 6 索引从 electrode_metadata.json(left+right)取, 属首猜, 需人工复核。"""
    e = _env(env)
    scf = Path(scf_dir)
    out = Path(outdir)
    values = make_tcontrol(scf, out / 'tcontrol', metadata=metadata, nsaos=nsaos,
                           efermi=efermi, surfaces=surfaces, nlayers=nlayers,
                           ener=ener, estep=estep, eend=eend,
                           valence_electrons=valence_electrons)
    sub = out / 'sub_te.sh'
    gres = str(e.get('te_gres') or '').strip()
    mkl_lib = str(e.get('mkl_lib') or '').strip()
    mkl_setup = ('export LD_LIBRARY_PATH=%s:$LD_LIBRARY_PATH' % mkl_lib) if mkl_lib else ''
    _write(sub,
        read_template('sub_te.sh').format(
            partition=e['partition'],
            mem=e.get('te_mem') or '8G',
            gres_line=('#SBATCH --gres=%s' % gres) if gres else '',
            walltime=e.get('te_walltime') or '02:00:00',
            mkl_setup=mkl_setup,
            cpus=e.get('te_cpus') or 16,
            aitranss_bin=e['aitranss_bin'] or 'UNCONFIGURED'))
    summary = {**values, 'aitranss_bin': e['aitranss_bin'] or 'UNCONFIGURED',
               'partition': e['partition'],
               'surface_indices_note': '首猜, 需人工复核 6 表面原子索引(24 章 §5)'}
    _write(out / 'input_summary.json', json.dumps(summary, indent=2, ensure_ascii=False))
    return {'tcontrol': str(out / 'tcontrol'), 'submit_script': str(sub),
            'summary': summary, 'values': values}