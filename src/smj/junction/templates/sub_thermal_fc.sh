#!/bin/bash
#SBATCH -p {partition}
#SBATCH -N 1
#SBATCH -n {ntasks}
#SBATCH --mem={mem}
#SBATCH --time={walltime}
#SBATCH {gres_line}
#SBATCH --array=0-{last}%{parallel}

# 有限位移力常数扫描:每个数组任务算一个 (原子, 轴, 符号) 位移点。
#
# 为什么用数组而不是 324 个独立作业:同一次 sbatch 提交、同一个 job id,
# 排一次队;而且 collect 只需要数 fc/*/forces.json 的份数就能判完整。
#
# 目录约定:fc/NNNN/{{geometry.in, control.in, aims.out, forces.json, .done}}
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

# A direct sbatch rerun (e.g. sub_lead_fc.sh) bypasses the submission helper.
# No worker can finalize until every current worker has removed the old bundle.
_job=${{SLURM_ARRAY_JOB_ID:?missing Slurm array job id}}
_tasks=${{SLURM_ARRAY_TASK_COUNT:?missing Slurm array size}}
rm -f .done forces.json aims.out ../../.dsh_done ../../forces_all.json ../../fc_summary.json || exit 1

run_slot() {{
    cp ../../control.in . || return 1
    python ../../apply_displacement.py --geometry ../../geometry.in \
        --manifest ../../displacements.json --index "$_i" --out geometry.in || return 1
    {aims_bin}
    _aims_rc=$?
    [ "$_aims_rc" -eq 0 ] || return "$_aims_rc"
    python ../../extract_forces.py aims.out forces.json
}}

run_slot
_rc=$?
printf '%s %s %s\n' "$_rc" "$(date -Is)" "$_job" > .done.tmp && mv .done.tmp .done || exit 1

# 逐个检查本批次的槽位，旧批次的额外目录不能顶替未完成的位移点。
_agg=0
_current=0
for ((_slot=0; _slot<{count}; _slot++)); do
    printf -v _marker '../%04d/.done' "$_slot"
    [ -f "$_marker" ] || exit "$_rc"
    read -r _status _stamp _owner < "$_marker" || exit 1
    [ "$_status" = 0 ] || _agg=1
    if [ "$_owner" = "$_job" ]; then ((_current+=1)); fi
done
[ "$_current" -eq "$_tasks" ] || exit "$_rc"

# 每批次只允许一个任务汇总，即使慢任务在发布后才走到这里。
if mkdir "../.finalize.lock.$_job" 2>/dev/null; then
    python ../../concat_forces.py ../ {count} ../../forces_all.json \
        --summary ../../fc_summary.json --allow-missing || _agg=1
    printf '%s %s\n' "$_agg" "$(date -Is)" > ../../.dsh_done.tmp \
        && mv ../../.dsh_done.tmp ../../.dsh_done || exit 1
fi
exit "$_rc"
