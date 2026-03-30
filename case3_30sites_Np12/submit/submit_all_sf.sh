#!/bin/bash
# 谱流计算：3 个节点各负责 5 个 k-sector，完成后自动合并+绘图
# 用法：bash submit/submit_all_sf.sh
# 在 case3_30sites_Np12/ 目录下运行

set -e

mkdir -p output_sf

J0=$(sbatch --parsable submit/submit_sf.sh 0  4)
J1=$(sbatch --parsable submit/submit_sf.sh 5  9)
J2=$(sbatch --parsable submit/submit_sf.sh 10 14)
JM=$(sbatch --parsable --dependency=afterok:${J0}:${J1}:${J2} submit/submit_merge_sf.sh)

echo "========================================"
echo "谱流计算 job（按 k-sector 分节点）:"
echo "  扇区[ 0- 4]: $J0"
echo "  扇区[ 5- 9]: $J1"
echo "  扇区[10-14]: $J2"
echo "  合并+绘图  : $JM  (依赖以上三个完成)"
echo "监控命令: squeue -j ${J0},${J1},${J2},${JM}"
echo "========================================"
