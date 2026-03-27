#!/bin/bash
# 链式提交：3 个 solve job + 1 个 merge job
# 用法：bash submit/submit_all.sh
# 在 case3_30sites_Np12/ 目录下运行

set -e

mkdir -p output

J0=$(sbatch --parsable submit/submit_solve.sh 0  4)
J1=$(sbatch --parsable submit/submit_solve.sh 5  9)
J2=$(sbatch --parsable submit/submit_solve.sh 10 14)
JM=$(sbatch --parsable --dependency=afterok:${J0}:${J1}:${J2} submit/submit_merge.sh)

echo "========================================"
echo "Solve jobs : $J0  $J1  $J2"
echo "Merge job  : $JM"
echo "监控命令   : squeue -j ${J0},${J1},${J2},${JM}"
echo "========================================"
