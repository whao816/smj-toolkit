#!/bin/bash
#SBATCH -p {partition}
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c {nproc}
#SBATCH --mem={mem}
#SBATCH --time={walltime}
{gres_line}

# Gaussian loading mode, rendered by make_gaussian_input from env:
#   module-load mode (env.g16_module non-empty, e.g. Gaussian/g16.c02):
#       module load Gaussian/g16.c02
#   g16root-path mode (env.g16_module empty):
#       export COM_DIR=...
#       source $g16root/bsd/g16.profile
{g16_setup}
# Scratch must exist before g16 opens its files; a stale cluster path used to
# abort the job with FIO/stdio error 2 before the first SCF step.
{scratch_setup}
srun g16 "$1"

# --- completion sentinel (P0-3, AGD-2026-09-25-SMJ): 退出码+时间戳 ---
# 注意: 渲染器用 str.format, 此处禁用花括号。
_dsh_rc=$?
_dsh_dir=$SLURM_SUBMIT_DIR
[ -n "$_dsh_dir" ] || _dsh_dir=.
printf '%s %s\n' "$_dsh_rc" "$(date -Is)" > "$_dsh_dir/.dsh_done"
exit "$_dsh_rc"
