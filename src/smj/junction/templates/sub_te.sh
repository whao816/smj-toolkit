#!/bin/bash
#SBATCH -p {partition}
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c {cpus}
#SBATCH --mem={mem}
{gres_line}
#SBATCH --time={walltime}

ulimit -s unlimited
# P1-4: aitranss 是 OpenMP 构建 —— 只给 1 核会慢约一个数量级(实测: 单线程 ~10.8 s/点,
# 全窗 4001 点 ~12 h, 即本会话 24513 被墙钟截断的成因; 16 线程 ~1.46 点/s 全窗 ~42 min;
# 精确倍率未严格测定, 故只声明"约一个数量级")。
export OMP_NUM_THREADS={cpus}
{mkl_setup}
{aitranss_bin} < tcontrol > te.out

# --- completion sentinel (P0-3, AGD-2026-09-25-SMJ): 退出码+时间戳 ---
# 谓词优先读本文件; 截断/被杀的 aitranss 据此被拦, 不再因 TE.dat 非空误判。
# 注意: 渲染器用 str.format, 此处禁用花括号。
_dsh_rc=$?
_dsh_dir=$SLURM_SUBMIT_DIR
[ -n "$_dsh_dir" ] || _dsh_dir=.
printf '%s %s\n' "$_dsh_rc" "$(date -Is)" > "$_dsh_dir/.dsh_done"
exit "$_dsh_rc"