# Multi-Node Parallel ED Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** 将 15 个 k 扇区拆分到 3 个 SLURM 节点并行求解，通过 JLD2 中间文件传递基态波函数，merge job 汇总能谱、计算结构因子并绘图。

**Architecture:** `main.jl` 读取环境变量 `SECTOR_START`/`SECTOR_END` 决定本节点求解哪些扇区，结果写入 `output/partial_N.jld2`；`merge.jl` 加载所有 partial 文件，合并排序，计算结构因子，输出 dat 文件和 PDF；三个 submit 脚本 + 一个链式提交入口。

**Tech Stack:** Julia 1.12, KrylovKit.jl, JLD2.jl, SparseArrays, SLURM (amd_256 partition)

---

### Task 1: 降低 kd（feature/csr-sparse worktree，立即修复 OOM）

**Files:**
- Modify: `case3_30sites_Np12/main.jl:36`

**Step 1: 修改 kd**

将第 36 行：
```julia
kd  = 120    # CSR 版：256 GB 内存约束（15×kd×5.77M×16B + 31GB CSR ≈ 203 GB）
             # 原 matrix-free 版 kd=200 需 314 GB（超限）
```
改为：
```julia
kd  = 100    # 降低峰值内存：15×100×5.77M×16B + 31GB CSR ≈ 169 GB（留 ~87 GB 余量）
```

同时更新第 84 行注释中的估算数字：
```julia
println("      内存预估：CSR ~31 GB + Krylov ~$(round(15*kd*5.77e6*16/1e9,digits=0)) GB")
```
（此行无需改动，kd 是变量，会自动反映新值）

**Step 2: 提交**

```bash
git add case3_30sites_Np12/main.jl
git commit -m "fix: reduce kd 120→100 to avoid OOM on 256GB node"
```

---

### Task 2: 创建 feature/multi-node 分支

**Step 1: 从 feature/csr-sparse 创建新分支**

在主 repo 目录（不是 worktree）执行：
```bash
cd /Users/shajianyu/CMP_manybody/partial_Hall_crystal
git worktree add .worktrees/multi-node -b feature/multi-node feature/csr-sparse
```

此后所有 Task 3–6 的文件操作均在 `.worktrees/multi-node/` 下进行。

**Step 2: 确认分支创建成功**

```bash
git worktree list
```
预期输出包含 `.worktrees/multi-node  [feature/multi-node]`

---

### Task 3: 修改 main.jl 支持扇区范围参数

**Files:**
- Modify: `case3_30sites_Np12/main.jl`（在 `.worktrees/multi-node/` 下）

**Step 1: 在文件顶部（include 之后，参数定义之前）添加 JLD2 依赖和扇区范围读取**

在 `using Printf, LinearAlgebra, Dates, SparseArrays` 这行后面加：
```julia
using JLD2
```

**Step 2: 在 BLAS 线程设置块之前添加扇区范围读取**

```julia
# ── 多节点参数：本节点负责的扇区索引范围（0-based，闭区间）──
seg_start = parse(Int, get(ENV, "SECTOR_START", "0"))
seg_end   = parse(Int, get(ENV, "SECTOR_END",   "14"))
seg_id    = seg_start   # 用于命名输出文件
```

**Step 3: 在 [3/5] 扇区构建完成后，添加扇区切片**

在 `total_reps = sum(...)` 那行之后加：
```julia
# ── 本节点只求解分配的扇区 ──
secs_local = secs[seg_start+1 : seg_end+1]   # Julia 1-based
@printf("  本节点负责扇区索引 %d–%d（共 %d 个）\n",
        seg_start, seg_end, length(secs_local))
flush(stdout)
```

**Step 4: 将 [4/5] 中的 `secs` 替换为 `secs_local`**

```julia
# 原：
all_ev, gs_vecs = compute_spectrum_sparse_with_vecs(
    secs, lat, hops0, V1, V2, V3;
    nev=nev, krylovdim=kd, verbose=true)
# 改为：
all_ev, gs_vecs = compute_spectrum_sparse_with_vecs(
    secs_local, lat, hops0, V1, V2, V3;
    nev=nev, krylovdim=kd, verbose=true)
```

**Step 5: 删除 [5/5] 结构因子部分**（merge.jl 负责，main.jl 不再计算）

删除从 `# ── 5. 结构因子 ──` 到文件末尾的所有内容，替换为：

```julia
# ── 5. 保存中间结果（供 merge.jl 使用）──
mkpath("output")
out_file = "output/partial_$(seg_id).jld2"
jldsave(out_file; ev_pairs=all_ev, gs_vecs=gs_vecs)
@printf("  中间结果保存: %s\n", out_file)

@printf("\n总耗时: %.1f s (%.2f min)\n", time()-t_total, (time()-t_total)/60)
println("完成时间: ", now())
println("完成!")
```

**Step 6: 更新步骤编号注释**（可选，保持注释清晰）

将 `[4/5]` 改为 `[4/4]`，将 `[5/5]` 改为 `[4/4 保存]`，或直接更新顶部说明。

**Step 7: 提交**

```bash
git add case3_30sites_Np12/main.jl
git commit -m "feat: main.jl reads SECTOR_START/SECTOR_END, saves partial JLD2"
```

---

### Task 4: 新增 merge.jl

**Files:**
- Create: `case3_30sites_Np12/merge.jl`

**完整内容：**

```julia
# ============================================================
# merge.jl — 汇总所有节点的 partial_*.jld2，计算结构因子，绘图
# 用法：julia merge.jl
# 前置条件：output/partial_0.jld2  partial_5.jld2  partial_10.jld2 均存在
# ============================================================
include("../shared/lattice.jl")
include("../shared/hoppings.jl")
include("../shared/basis.jl")
include("../shared/ksector.jl")
include("../shared/structure_factor.jl")

using Printf, Dates, JLD2

Np  = 12
V1  = 10.0
V2  = 2.0
V3  = 2.0
t1  = 1.0
t3  = 0.2

println("="^60)
println("Merge — 汇总 3 个节点的结果")
println("开始时间: ", now())
println("="^60)

t_total = time()

# ── 1. 重建晶格与完整扇区列表（供结构因子使用）──
print("[1/4] 重建晶格与扇区... "); flush(stdout)
t = @elapsed begin
    lat   = TiltedLat30()
    basis = gen_basis(lat.Ns, Np)
    secs  = KSector[]
    for m in lat.ktab
        push!(secs, build_ksector(basis, lat, m))
    end
end
@printf("done  Ns=%d  Nuc=%d  (%.2f s)\n", lat.Ns, lat.Nuc, t)

# ── 2. 加载所有 partial 文件 ──
println("[2/4] 加载 partial JLD2 文件...")
partial_files = [
    ("output/partial_0.jld2",  "k=0..4"),
    ("output/partial_5.jld2",  "k=5..9"),
    ("output/partial_10.jld2", "k=10..14"),
]

all_ev   = Tuple{Int,Float64}[]
gs_vecs  = Dict{Int,Vector{ComplexF64}}()

for (fpath, label) in partial_files
    isfile(fpath) || error("找不到文件: $fpath")
    data = load(fpath)
    append!(all_ev, data["ev_pairs"])
    merge!(gs_vecs, data["gs_vecs"])
    @printf("  ✓ %s  (%s)\n", fpath, label)
end
sort!(all_ev, by=x->x[2])
E0   = all_ev[1][2]
gs_k = all_ev[1][1]
@printf("  全局基态: k=%d  E0=%.10f\n", gs_k, E0)

# ── 3. 输出能谱文件 ──
println("[3/4] 输出结果文件...")
mkpath("output")

open("output/spectrum_Np$(Np).dat","w") do f
    println(f,"# k  E-E0  [30sites Np=$Np V1=$V1 V2=$V2 V3=$V3 t'=$t3 multi-node]")
    for (m,e) in all_ev
        @printf(f,"%d  %.10f\n", m, e-E0)
    end
end
println("  能谱保存: output/spectrum_Np$(Np).dat")

# ── 结构因子 ──
gs_sec = secs[findfirst(s -> s.m == gs_k, secs)]
print("  全局基态结构因子... "); flush(stdout)
t_sq = @elapsed begin
    gs_vec = gs_vecs[gs_k]
    results_ord, results_srt, norm2, Np_check =
        structure_factor_survey(gs_sec, gs_vec, lat, Np)
end
@printf("done (%.1f s)\n", t_sq)
print_sq_report(results_srt, norm2, Np_check, Np, lat.Ns, gs_k)

open("output/sq_Np$(Np).dat","w") do f
    println(f,"# ik  N(k)  kx  ky  [30sites Np=$Np V1=$V1 t'=$t3]")
    println(f,"# 基态扇区 k=$gs_k  E0=$(round(E0,digits=8))")
    for (m, sq, kx, ky) in results_ord
        @printf(f,"%d  %.6f  %.6f  %.6f\n", m, sq, kx, ky)
    end
end
println("  结构因子保存: output/sq_Np$(Np).dat")

print("  全扇区结构因子矩阵... "); flush(stdout)
t_all = @elapsed open("output/sq_all_Np$(Np).dat","w") do f
    println(f,"# sector_k  q  N(q)  kx  ky")
    for sec in secs
        haskey(gs_vecs, sec.m) || continue
        res_ord, _, _, _ = structure_factor_survey(sec, gs_vecs[sec.m], lat, Np)
        for (m, sq, kx, ky) in res_ord
            @printf(f,"%d  %d  %.6f  %.6f  %.6f\n", sec.m, m, sq, kx, ky)
        end
    end
end
@printf("done (%.1f s)\n", t_all)
println("  全扇区结构因子保存: output/sq_all_Np$(Np).dat")

# ── 4. 绘图 ──
println("[4/4] 绘图...")
# 将输出路径传给 plot.jl 的约定：plot.jl 读取 output/ 下的 dat 文件
include("plot.jl")

@printf("\n总耗时: %.1f s (%.2f min)\n", time()-t_total, (time()-t_total)/60)
println("完成时间: ", now())
println("完成!")
```

> **注意**：`plot.jl` 目前可能硬编码读取当前目录下的 dat 文件。如果 `plot.jl` 需要适配 `output/` 路径，在 Task 4 commit 后检查并做最小改动。

**Step 1: 提交**

```bash
git add case3_30sites_Np12/merge.jl
git commit -m "feat: add merge.jl for multi-node result aggregation"
```

---

### Task 5: 新增 submit/ 目录和三个脚本

**Files:**
- Create: `case3_30sites_Np12/submit/submit_solve.sh`
- Create: `case3_30sites_Np12/submit/submit_merge.sh`
- Create: `case3_30sites_Np12/submit/submit_all.sh`

#### submit_solve.sh

```bash
#!/bin/bash
#SBATCH --job-name=ED_Np12_solve
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=amd_256
#SBATCH --time=240:00:00
#SBATCH --output=output/solve_%j_k${1}-${2}.out
#SBATCH --error=output/solve_%j_k${1}-${2}.err

export JULIA_NUM_THREADS=64
export OMP_NUM_THREADS=1
export SECTOR_START=$1
export SECTOR_END=$2

JULIA=~/sc71394/shajy/software/env_soft/julia_1.12.4/julia-1.12.4/bin/julia

echo "========================================"
echo "Job ID:    $SLURM_JOB_ID"
echo "Node:      $SLURMD_NODENAME"
echo "Sectors:   $SECTOR_START – $SECTOR_END"
echo "Start:     $(date)"
$JULIA --version
echo "========================================"

mkdir -p output
cd $SLURM_SUBMIT_DIR

$JULIA --check-bounds=no --threads=64 main.jl

EXIT_CODE=$?
echo "End: $(date)   Exit: $EXIT_CODE"
exit $EXIT_CODE
```

#### submit_merge.sh

```bash
#!/bin/bash
#SBATCH --job-name=ED_Np12_merge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=amd_256
#SBATCH --time=04:00:00
#SBATCH --output=output/merge_%j.out
#SBATCH --error=output/merge_%j.err

export JULIA_NUM_THREADS=16
export OMP_NUM_THREADS=1

JULIA=~/sc71394/shajy/software/env_soft/julia_1.12.4/julia-1.12.4/bin/julia

echo "========================================"
echo "Job ID:    $SLURM_JOB_ID"
echo "Node:      $SLURMD_NODENAME"
echo "Start:     $(date)"
$JULIA --version
echo "========================================"

cd $SLURM_SUBMIT_DIR

$JULIA --check-bounds=no --threads=16 merge.jl

EXIT_CODE=$?
echo "End: $(date)   Exit: $EXIT_CODE"
exit $EXIT_CODE
```

#### submit_all.sh

```bash
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
```

**Step 1: 提交**

```bash
git add case3_30sites_Np12/submit/
git commit -m "feat: add submit scripts for multi-node parallel ED"
```

---

### Task 6: 添加 output/.gitkeep 和 .gitignore 条目

**Files:**
- Create: `case3_30sites_Np12/output/.gitkeep`
- Modify: `.gitignore`（根目录）

**Step 1: 创建占位文件**

```bash
mkdir -p case3_30sites_Np12/output
touch case3_30sites_Np12/output/.gitkeep
```

**Step 2: 在 .gitignore 中追加**

```
# multi-node 运行产物
case3_30sites_Np12/output/*.jld2
case3_30sites_Np12/output/*.dat
case3_30sites_Np12/output/*.pdf
case3_30sites_Np12/output/*.out
case3_30sites_Np12/output/*.err
```

**Step 3: 确认 .gitkeep 不被忽略**（上述 gitignore 用 `*.ext` 精确匹配，不会忽略 `.gitkeep`）

**Step 4: 提交**

```bash
git add case3_30sites_Np12/output/.gitkeep .gitignore
git commit -m "chore: add output/.gitkeep and gitignore rules for multi-node artifacts"
```

---

### Task 7: 检查 plot.jl 路径兼容性

**Files:**
- Read: `case3_30sites_Np12/plot.jl`
- Modify（如需）: `case3_30sites_Np12/plot.jl` 或 `merge.jl`

**Step 1: 检查 plot.jl 读取的文件路径**

查看 `plot.jl` 是否硬编码读取 `spectrum_Np12.dat`（当前目录）或已有路径参数。

若 `plot.jl` 读取当前目录，有两种选项：
- **最小改动**：在 `merge.jl` 中，于 `include("plot.jl")` 之前临时 `cd("output")`，绘图后 `cd("..")`
- **更清晰**：给 `plot.jl` 添加可选的 `DATA_DIR` 变量支持

选择最小改动方案即可。

**Step 2: 若有改动则提交**

```bash
git add case3_30sites_Np12/plot.jl case3_30sites_Np12/merge.jl
git commit -m "fix: align plot.jl data path with output/ directory"
```

若无需改动，跳过此步。

---

### Task 8: 本地冒烟测试（小规模验证）

在本地用 `test_local.jl` 的参数规模（Np=4~6，少量扇区）验证多节点逻辑：

**Step 1: 模拟 solve job 0**

```bash
cd case3_30sites_Np12
SECTOR_START=0 SECTOR_END=4 julia --threads=4 main.jl
```

预期：生成 `output/partial_0.jld2`，无报错

**Step 2: 模拟 solve job 1 和 2**

```bash
SECTOR_START=5  SECTOR_END=9  julia --threads=4 main.jl
SECTOR_START=10 SECTOR_END=14 julia --threads=4 main.jl
```

**Step 3: 运行 merge**

```bash
julia --threads=4 merge.jl
```

预期：生成 `output/spectrum_Np12.dat`、`output/sq_Np12.dat`、PDF 图

**Step 4: 与单节点结果对比**（可选）

若已有单节点跑出的 `spectrum_Np5.dat`（test_local 结果），对比能量值确认一致。

---

## 执行顺序总结

| Task | 分支 | 说明 |
|------|------|------|
| 1 | feature/csr-sparse | 降 kd，立即可重提服务器 |
| 2 | — | 创建 multi-node worktree |
| 3–7 | feature/multi-node | 代码实现 |
| 8 | feature/multi-node | 本地验证 |
