#!/bin/bash
#SBATCH -p {partition}
#SBATCH -N 1
#SBATCH -n {ntasks}
#SBATCH --mem={mem}
#SBATCH --time={walltime}
#SBATCH --array=0-{last}%{parallel}

# 有限位移力常数扫描:每个数组任务算一个 (原子, 轴, 符号) 位移点。
#
# 为什么用数组而不是 324 个独立作业:同一次 sbatch 提交、同一个 job id,
# 排一次队;而且 collect 只需要数 fc/*/forces.json 的份数就能判完整。
#
# 目录约定:fc/NNNN/{{geometry.in, control.in, aims.out, forces.json}}
#   NNNN = 四位补零的数组下标,与 displacements.json 的 slot 一一对应。
{module_load}
{exe_define}
export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export MKL_NUM_THREADS=1
ulimit -s unlimited || true
ulimit -c 0

cd "$SLURM_SUBMIT_DIR" || exit 1
_i=$SLURM_ARRAY_TASK_ID
_d=$(printf 'fc/%04d' "$_i")
mkdir -p "$_d" || exit 1
cd "$_d" || exit 1

cp ../../control.in . || exit 1
python ../../apply_displacement.py --geometry ../../geometry.in \
    --manifest ../../displacements.json --index "$_i" --out geometry.in || exit 1

{aims_bin}
_rc=$?
if [ "$_rc" -eq 0 ]; then
    python ../../extract_forces.py aims.out > forces.json || _rc=1
fi
printf '%s %s\n' "$_rc" "$(date -Is)" > .done

# --- 阶段级哨兵:最后一个完成的任务负责汇总 ---
# 每个任务写完自己的 .done 后数一次 fc/*/done 的份数;数到 {count} 的那个
# 任务再抢 .finalize.lock(原子性 mkdir),抢到的在锁内**重新数一次**并汇总退出码。
_count=$(ls -1 ../*/done 2>/dev/null | wc -l)
if [ "$_count" -ge {count} ]; then
    if mkdir ../.finalize.lock 2>/dev/null; then
        _total=$(ls -1 ../*/done 2>/dev/null | wc -l)
        _agg=0
        if [ "$_total" -ge {count} ]; then
            for _f in ../*/done; do
                _v=$(cat "$_f" 2>/dev/null)
                _v=${{_v%% *}}
                case "$_v" in
                    0) ;;
                    *) _agg=1 ;;
                esac
            done
        else
            _agg=1
        fi
        printf '%s %s\n' "$_agg" "$(date -Is)" > ../../.dsh_done
        rmdir ../.finalize.lock 2>/dev/null || true
    fi
fi
exit "$_rc"
