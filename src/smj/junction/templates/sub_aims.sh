#!/bin/bash
#SBATCH -p {partition}
#SBATCH -N 1
#SBATCH -n {ntasks}
#SBATCH --mem={mem}
{gres_line}
#SBATCH --time={walltime}

# FHI-aims loading: module load + FHI_AIMS_BIN (CUDA GPU build 221103),
#   or direct binary path (aims_bin, e.g. scalapack/mpi single-file build).
#   GPU mode (env.aims_gpus non-empty, e.g. 1): 1 GPU maps to 1 MPI rank (1 core).
#   CPU mode (env.aims_gpus empty): use {ntasks} cores.
{module_load}
{exe_define}
export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export MKL_NUM_THREADS=1
ulimit -s unlimited || true
ulimit -c 0
{aims_bin}

# --- completion sentinel (P0-3, AGD-2026-09-25-SMJ): 退出码+时间戳 ---
# 谓词优先读本文件; 失败作业据此被拦, 不再因 "产物非空" 误判完成。
# 注意: 渲染器用 str.format, 此处禁用花括号。
_dsh_rc=$?
_dsh_dir=$SLURM_SUBMIT_DIR
[ -n "$_dsh_dir" ] || _dsh_dir=.
printf '%s %s\n' "$_dsh_rc" "$(date -Is)" > "$_dsh_dir/.dsh_done"
exit "$_dsh_rc"