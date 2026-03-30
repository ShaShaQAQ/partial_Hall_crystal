#!/bin/bash
# 谱流计算：3 个节点并行，各处理 10 个 φ 点，完成后合并
# 用法：bash submit/submit_all_sf.sh
# 在 case3_30sites_Np12/ 目录下运行

set -e

mkdir -p output_sf

# 3 个节点，各处理 φ 索引 [0-9], [10-19], [20-29]
J0=$(sbatch --parsable submit/submit_sf.sh 0  9)
J1=$(sbatch --parsable submit/submit_sf.sh 10 19)
J2=$(sbatch --parsable submit/submit_sf.sh 20 29)

echo "========================================"
echo "谱流计算 job:"
echo "  φ[ 0- 9]: $J0"
echo "  φ[10-19]: $J1"
echo "  φ[20-29]: $J2"
echo "监控命令: squeue -j ${J0},${J1},${J2}"
echo ""
echo "完成后运行合并和绘图："
echo "  julia merge_sf.jl && julia plot_sf.jl"
echo "========================================"
