#!/bin/bash
#SBATCH -p {partition}
#SBATCH -N 1
#SBATCH -n {ntasks}
#SBATCH --mem={mem}
#SBATCH --time={walltime}
#SBATCH {gres_line}

# 热导链的弛豫作业:control.in 已由 prepare 按 control_thermal_relax.in 写好
# (compute_forces + final_forces_cleaned),这里只负责跑。
{module_load}
{exe_define}
export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export MKL_NUM_THREADS=1
ulimit -s unlimited || true
ulimit -c 0
{aims_bin}

# --- completion sentinel ---
_dsh_rc=$?
_dsh_dir=$SLURM_SUBMIT_DIR
[ -n "$_dsh_dir" ] || _dsh_dir=.
printf '%s %s\n' "$_dsh_rc" "$(date -Is)" > "$_dsh_dir/.dsh_done"
exit "$_dsh_rc"
