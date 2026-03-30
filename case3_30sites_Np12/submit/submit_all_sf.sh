#!/bin/bash
# 谱流计算：3 个节点各负责 5 个 k-sector，对全部通量点计算
# 用法：bash submit/submit_all_sf.sh
# 在 case3_30sites_Np12/ 目录下运行

set -e

mkdir -p output_sf

# 3 个节点，各负责 5 个扇区（0-based），跑所有 40 个 φ 点
J0=$(sbatch --parsable submit/submit_sf.sh 0  4)
J1=$(sbatch --parsable submit/submit_sf.sh 5  9)
J2=$(sbatch --parsable submit/submit_sf.sh 10 14)

echo "========================================"
echo "谱流计算 job（按 k-sector 分节点）:"
echo "  扇区[ 0- 4]: $J0"
echo "  扇区[ 5- 9]: $J1"
echo "  扇区[10-14]: $J2"
echo "监控命令: squeue -j ${J0},${J1},${J2}"
echo ""
echo "完成后运行合并和绘图："
echo "  julia merge_sf.jl && julia plot_sf.jl"
echo "========================================"
