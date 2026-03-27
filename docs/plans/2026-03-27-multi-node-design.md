# 多节点并行 ED 设计文档

**日期**: 2026-03-27
**分支**: feature/multi-node（从 feature/csr-sparse 分叉）
**背景**: case3 30 sites Np=12，15 个 k 扇区，单节点 256 GB 内存不足以同时跑所有扇区的 Lanczos。

---

## 问题

单节点并行（`Threads.@threads` 跑全部 15 扇区）峰值内存 ~220 GB，加上 OS/Julia 开销超出 amd_256 节点 256 GB 限制，导致 SIGKILL（退出码 9:0）。

---

## 方案

**方案 A（采用）**：固定拆分，3 个 solve job + 1 个 merge job。

每个 solve job 独立重建晶格和动量扇区（耗时短，不值得共享），只解自己负责的 5 个扇区，输出 JLD2 中间文件。merge job 依赖所有 solve job 成功后启动，汇总结果并绘图。

---

## Job 拓扑

```
[job_solve_0]  k=0..4   → output/partial_0.jld2
[job_solve_1]  k=5..9   → output/partial_1.jld2   ──→  [job_merge]
[job_solve_2]  k=10..14 → output/partial_2.jld2
```

- 3 个 solve job 同时提交，互相独立
- merge job：`--dependency=afterok:$J0:$J1:$J2`
- `submit_all.sh` 负责链式提交并打印所有 job ID

---

## 中间文件格式（JLD2）

每个 `output/partial_N.jld2` 包含：

| 键 | 类型 | 内容 |
|----|------|------|
| `ev_pairs` | `Vector{Tuple{Int,Float64}}` | (k, E) 对，本节点所有扇区 |
| `gs_vecs` | `Dict{Int,Vector{ComplexF64}}` | k → 基态波函数 |

---

## 内存估算（每节点 5 扇区）

| 部分 | 大小 |
|------|------|
| CSR 矩阵（5 扇区） | ~10 GB |
| Krylov 向量（5 × kd × dim × 16B） | ~55 GB（kd=120）|
| 合计 | **~65 GB**，远低于 256 GB |

---

## 文件结构

```
case3_30sites_Np12/
├── main.jl              # 改动：读 SECTOR_START/SECTOR_END 环境变量
├── merge.jl             # 新增：加载 partial_*.jld2，合并能谱，结构因子，绘图
├── plot.jl              # 不变（merge.jl include）
├── test_local.jl        # 不变
├── submit/
│   ├── submit_all.sh    # 入口脚本
│   ├── submit_solve.sh  # solve job 模板
│   └── submit_merge.sh  # merge job 模板
└── output/              # .gitignore 忽略，运行时产生
    ├── partial_0.jld2
    ├── partial_1.jld2
    ├── partial_2.jld2
    ├── spectrum_Np12.dat
    ├── sq_Np12.dat
    └── *.pdf
```

---

## 代码改动要点

### main.jl

```julia
seg_start = parse(Int, get(ENV, "SECTOR_START", "0"))
seg_end   = parse(Int, get(ENV, "SECTOR_END",   "4"))
secs_local = secs[seg_start+1 : seg_end+1]   # Julia 1-based

# ... compute_spectrum_sparse_with_vecs(secs_local, ...) ...

using JLD2
mkpath("output")
jldsave("output/partial_$(seg_start).jld2"; ev_pairs=all_ev, gs_vecs=gs_vecs)
```

### merge.jl

1. 加载 `output/partial_0.jld2`、`partial_5.jld2`、`partial_10.jld2`
2. 合并 `ev_pairs`，全局排序得完整能谱
3. 合并 `gs_vecs`，调用 `structure_factor_survey` 计算 N(q)
4. 写 `spectrum_Np12.dat`、`sq_Np12.dat`、`sq_all_Np12.dat`
5. include `plot.jl` 生成 PDF

### submit_all.sh

```bash
J0=$(sbatch --parsable submit/submit_solve.sh 0 4)
J1=$(sbatch --parsable submit/submit_solve.sh 5 9)
J2=$(sbatch --parsable submit/submit_solve.sh 10 14)
JM=$(sbatch --parsable --dependency=afterok:$J0:$J1:$J2 submit/submit_merge.sh)
echo "Solve jobs: $J0 $J1 $J2  |  Merge job: $JM"
```

---

## 分支管理

- `feature/csr-sparse`：单节点 CSR 版本，降 kd=100 修复 OOM
- `feature/multi-node`：从 `feature/csr-sparse` 分叉，多节点拆分方案
- `master`：稳定版本，待两个分支验证后合并
