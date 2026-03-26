# CSR 稀疏矩阵重构实现计划

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** 将 matrix-free `Hv!`（随机哈希查找）替换为预建 CSR 稀疏矩阵的 `mul!(w, H, v)`，在 256 GB 单节点服务器上加速 Lanczos 求解。

**Architecture:** 新增 `build_sparse_H`（COO→CSC，含对角+hopping）；新增 `compute_spectrum_sparse_with_vecs`（并行建 CSR → 并行 Lanczos → GC）；更新 `case3_30sites_Np12/main.jl` 使用新流程。`krylovdim` 从 200 降至 120（Np=12）以适应 256 GB 内存。

**Tech Stack:** Julia, SparseArrays（标准库），KrylovKit.jl，Threads.@threads，BLAS.set_num_threads

---

### Task 1: 在 `hamiltonian.jl` 中添加 `build_sparse_H`

**Files:**
- Modify: `shared/hamiltonian.jl`（在文件末尾追加）

**Step 1: 确认 SparseArrays 可用**

```bash
cd /Users/shajianyu/CMP_manybody/partial_Hall_crystal/.worktrees/csr-sparse
julia -e "using SparseArrays; println(\"ok\")"
```
Expected: `ok`

**Step 2: 在 `shared/hamiltonian.jl` 末尾追加 `build_sparse_H`**

在文件末尾（第 126 行之后）添加：

```julia
# ============================================================
# CSR 预建稀疏哈密顿量
#
# 一次性将 Hv! 的所有操作展开为 SparseMatrixCSC 的非零元，
# 后续 Hv = mul!(w, H, v)，顺序内存访问，BLAS 优化。
# ============================================================
using SparseArrays

"""
    build_sparse_H(sec, lat, hops, V1, V2, V3) -> SparseMatrixCSC{ComplexF64,Int32}

将动量扇区 sec 的哈密顿量预建为稀疏矩阵（CSR 格式）。
包含：对角势能（内联 precompute_diag_H 逻辑）+ off-diagonal 跳跃项。
fock2rep 哈希查找只发生一次（此函数内），之后 Lanczos 全程无 Dict 访问。
"""
function build_sparse_H(sec::KSector, lat::GenLat,
                        hops::Vector{Tuple{Int,Int,ComplexF64}},
                        V1::Float64, V2::Float64, V3::Float64)
    Nrep = length(sec.reps)
    Nrep == 0 && return spzeros(ComplexF64, Int32, 0, 0)
    Ns = lat.Ns

    Is = Int32[]
    Js = Int32[]
    Vs = ComplexF64[]
    sizehint!(Is, Nrep * 50)
    sizehint!(Js, Nrep * 50)
    sizehint!(Vs, Nrep * 50)

    for i in 1:Nrep
        ni  = 1.0 / sec.norms[i]
        ni2 = ni * ni

        # ── 对角势能（内联 precompute_diag_H） ──
        diag_val = 0.0
        for (c_α, F) in sec.orbit_data[i]
            Eint = 0.0
            for s in 1:Ns
                (F>>(s-1))&1 == 0 && continue
                for k in 1:6
                    Eint += V1 * ((F>>(lat.nb1[s,k]-1))&1) * 0.5
                    Eint += V2 * ((F>>(lat.nb2[s,k]-1))&1) * 0.5
                    Eint += V3 * ((F>>(lat.nb3[s,k]-1))&1) * 0.5
                end
            end
            diag_val += abs2(c_α) * ni2 * Eint
        end
        push!(Is, Int32(i)); push!(Js, Int32(i)); push!(Vs, ComplexF64(diag_val))

        # ── 跳跃（off-diagonal） ──
        for (c_α, F) in sec.orbit_data[i]
            coeff_base = c_α * ni
            for (tgt, src, tval) in hops
                (F>>(src-1))&1 == 0 && continue
                (F>>(tgt-1))&1 == 1 && continue
                newF, sgn_hop = hop(F, tgt, src)
                info = get(sec.fock2rep, newF, nothing)
                info === nothing && continue
                (j, c_conj_norm) = info
                push!(Is, Int32(j))
                push!(Js, Int32(i))
                push!(Vs, coeff_base * tval * sgn_hop * c_conj_norm)
            end
        end
    end

    # Julia sparse() 自动合并重复 (i,j) 索引（求和），正确处理多轨道贡献
    return sparse(Is, Js, Vs, Nrep, Nrep)
end
```

**Step 3: Commit**

```bash
cd /Users/shajianyu/CMP_manybody/partial_Hall_crystal/.worktrees/csr-sparse
git add shared/hamiltonian.jl
git commit -m "feat: add build_sparse_H to hamiltonian.jl"
```

---

### Task 2: 写测试脚本验证 CSR 与 matrix-free 结果一致

**Files:**
- Create: `test_csr.jl`（worktree 根目录）

**Step 1: 创建测试脚本**

```julia
# test_csr.jl — 验证 CSR 与 matrix-free 本征值一致（在本地 Mac 运行）
# 使用 RectLat3x6 (18 sites) + Np=4，C(18,4)=3060，秒级完成
include("shared/lattice.jl")
include("shared/hoppings.jl")
include("shared/basis.jl")
include("shared/ksector.jl")
include("shared/hamiltonian.jl")
include("shared/solver.jl")

using SparseArrays, LinearAlgebra, KrylovKit, Printf

t1 = 1.0; t3 = 0.2; V1 = 1.0; V2 = 0.0; V3 = 0.0; Np = 4
lat   = RectLat3x6()
basis = gen_basis(lat.Ns, Np)
hops  = build_hops(lat, t1, t3, 0.0)
secs  = build_all_ksectors(basis, lat)

@printf("RectLat3x6  Ns=%d  Nuc=%d  C(%d,%d)=%d\n",
        lat.Ns, lat.Nuc, lat.Ns, Np, length(basis))

max_err = 0.0
for sec in secs
    length(sec.reps) == 0 && continue

    # matrix-free（原版）
    diag_H = precompute_diag_H(sec, lat, V1, V2, V3)
    v0 = normalize(randn(ComplexF64, length(sec.reps)))
    vals_mf, _, _ = eigsolve(
        v -> (w = similar(v); Hv!(w, v, sec, lat, hops, V1, V2, V3, diag_H); w),
        v0, 4, :SR; ishermitian=true, krylovdim=50, tol=1e-12)

    # CSR（新版）
    H_csr = build_sparse_H(sec, lat, hops, V1, V2, V3)
    # H 应为 Hermitian：检验 ||H - H'|| < 1e-10
    herm_err = maximum(abs.(H_csr - H_csr'))
    herm_err > 1e-10 && error("k=$(sec.m): H 不是 Hermitian，误差 $herm_err")

    vals_csr, _, _ = eigsolve(
        v -> H_csr * v,
        v0, 4, :SR; ishermitian=true, krylovdim=50, tol=1e-12)

    n   = min(length(vals_mf), length(vals_csr), 4)
    err = maximum(abs.(sort(real.(vals_mf[1:n])) .- sort(real.(vals_csr[1:n]))))
    max_err = max(max_err, err)
    @printf("  k=%2d  dim=%4d  nnz=%6d  |ΔE|_max=%.2e  %s\n",
            sec.m, length(sec.reps), nnz(H_csr), err,
            err < 1e-8 ? "✓" : "✗ FAIL")
end

println()
println(max_err < 1e-8 ? "✓ ALL PASS  max_err=$(max_err)" :
                         "✗ FAIL     max_err=$(max_err)")
```

**Step 2: 运行测试**

```bash
cd /Users/shajianyu/CMP_manybody/partial_Hall_crystal/.worktrees/csr-sparse
julia test_csr.jl
```

Expected output（每行 `✓`，最后一行 `✓ ALL PASS`，`max_err < 1e-8`）：
```
RectLat3x6  Ns=18  Nuc=9  C(18,4)=3060
  k= 0  dim= 340  nnz= ...  |ΔE|_max=...e-...  ✓
  ...
✓ ALL PASS  max_err=...e-...
```

**Step 3: Commit**

```bash
git add test_csr.jl
git commit -m "test: add test_csr.jl verifying CSR vs matrix-free eigenvalues"
```

---

### Task 3: 在 `solver.jl` 中添加 `compute_spectrum_sparse_with_vecs`

**Files:**
- Modify: `shared/solver.jl`（在文件末尾追加）

**Step 1: 在文件末尾追加新函数**

```julia
# ============================================================
# CSR 稀疏矩阵求解器
# 流程：
#   A. Threads.@threads 并行 build_sparse_H（一次性，fock2rep 查找在此完成）
#   B. Threads.@threads 并行 Lanczos（mul!(w,H,v)，顺序内存访问）
#   C. 返回 (all_ev, gs_vecs)，供后处理直接使用
# ============================================================

"""
    compute_spectrum_sparse_with_vecs(secs, lat, hops, V1, V2, V3; ...)

用预建 CSR 稀疏矩阵做并行 Lanczos。
步骤 A 并行建矩阵，步骤 B 并行求解，返回格式与 compute_spectrum_with_vecs 相同。

krylovdim 建议值：
  Np=12, 30 sites（dim≈5.77M）：krylovdim=120（203 GB peak）
  Np=13, 30 sites（dim≈7.98M）：krylovdim=80 （204 GB peak）
"""
function compute_spectrum_sparse_with_vecs(
        secs::Vector{KSector}, lat::GenLat,
        hops::Vector{Tuple{Int,Int,ComplexF64}},
        V1::Float64, V2::Float64, V3::Float64;
        nev::Int=4, krylovdim::Int=120, verbose::Bool=false)

    n  = length(secs)
    lk = ReentrantLock()

    # ── A. 并行建 CSR ──
    @printf("  [CSR-A] 并行构建 %d 个稀疏 H 矩阵 (%d 线程)...\n",
            n, Threads.nthreads())
    flush(stdout)
    H_csrs = Vector{SparseMatrixCSC{ComplexF64,Int32}}(undef, n)
    t_build = @elapsed Threads.@threads for i in 1:n
        H_csrs[i] = build_sparse_H(secs[i], lat, hops, V1, V2, V3)
    end
    nnz_total = sum(nnz(H) for H in H_csrs)
    @printf("  [CSR-A] 完成  耗时 %.1f s  总非零元 %d  内存估算 %.2f GB\n",
            t_build, nnz_total, nnz_total * 12 / 1e9)
    flush(stdout)

    # ── B. 并行 Lanczos ──
    @printf("  [CSR-B] 并行 Lanczos  krylovdim=%d  nev=%d...\n", krylovdim, nev)
    flush(stdout)
    all_res  = Vector{Vector{Tuple{Int,Float64}}}(undef, n)
    all_vecs = Vector{Vector{ComplexF64}}(undef, n)

    t_solve = @elapsed Threads.@threads for i in 1:n
        sec  = secs[i]
        H    = H_csrs[i]
        Nrep = length(sec.reps)

        if Nrep == 0
            all_res[i]  = Tuple{Int,Float64}[]
            all_vecs[i] = ComplexF64[]
            continue
        end

        nev_actual = min(nev, Nrep)
        kd  = min(krylovdim, Nrep)
        v0  = normalize(randn(ComplexF64, Nrep))
        n_hv = Ref(0); t_hv = Ref(0.0)

        function Hf(v)
            result = similar(v)
            dt = @elapsed mul!(result, H, v)
            n_hv[] += 1; t_hv[] += dt
            if verbose && n_hv[] % 50 == 0
                msg = @sprintf("    [k=%2d] H·v #%3d  均 %.5f s/次\n",
                               sec.m, n_hv[], t_hv[] / n_hv[])
                lock(lk) do; print(msg); flush(stdout); end
            end
            return result
        end

        t_eig = @elapsed vals, vecs, info = eigsolve(
            Hf, v0, nev_actual, :SR;
            ishermitian=true, krylovdim=kd, maxiter=500, tol=1e-10)

        info.converged < nev_actual &&
            @warn "k=$(sec.m) 只收敛 $(info.converged)/$nev_actual 个本征值"

        nhv = n_hv[]
        msg = @sprintf(
            "  ✓ [k=%2d] dim=%7d  H·v×%3d  均 %.5f s/次  Lanczos %.1f s  E0=%.8f\n",
            sec.m, Nrep, nhv, nhv > 0 ? t_hv[] / nhv : 0.0, t_eig,
            isempty(vals) ? NaN : real(vals[1]))
        lock(lk) do; print(msg); flush(stdout); end

        n_vals = min(nev_actual, length(vals))
        all_res[i]  = [(sec.m, e) for e in sort(real.(vals[1:n_vals]))]
        all_vecs[i] = isempty(vecs) ? ComplexF64[] : vecs[1]
    end
    @printf("  [CSR-B] Lanczos 完成  耗时 %.1f s (%.2f min)\n",
            t_solve, t_solve / 60)
    flush(stdout)

    all_ev = Tuple{Int,Float64}[]
    for res in all_res
        isassigned(all_res, findfirst(==(res), all_res)) && append!(all_ev, res)
    end
    sort!(all_ev, by=x->x[2])

    gs_vecs = Dict{Int,Vector{ComplexF64}}(
        secs[i].m => all_vecs[i]
        for i in 1:n if isassigned(all_vecs, i) && !isempty(all_vecs[i])
    )
    return all_ev, gs_vecs
end
```

**Step 2: 运行原有测试确认未破坏旧接口**

```bash
cd /Users/shajianyu/CMP_manybody/partial_Hall_crystal/.worktrees/csr-sparse
julia test_csr.jl
```

Expected: `✓ ALL PASS`（测试脚本 include 了 solver.jl，确认加载不报错）

**Step 3: Commit**

```bash
git add shared/solver.jl
git commit -m "feat: add compute_spectrum_sparse_with_vecs to solver.jl"
```

---

### Task 4: 更新 `case3_30sites_Np12/main.jl`

**Files:**
- Modify: `case3_30sites_Np12/main.jl`（在 worktree 内）

**Step 1: 修改参数区（第 29-36 行）**

将：
```julia
t1  = 1.0
t3  = 0.2
V1  = 1.0
V2  = 0.0
V3  = 0.0
Np  = 12
nev = 8
kd  = 200    # 大扇区需要更大 Krylov 子空间保证收敛
```

改为：
```julia
t1  = 1.0
t3  = 0.2
V1  = 1.0
V2  = 0.0
V3  = 0.0
Np  = 12
nev = 8
kd  = 120    # CSR 版：256 GB 内存约束（15×kd×5.77M×16B + 31GB CSR ≈ 203 GB）
             # 原 matrix-free 版 kd=200 需 314 GB（超限）
```

**Step 2: 在 `using` 行添加 SparseArrays**

将第 18 行：
```julia
using Printf, LinearAlgebra, Dates
```
改为：
```julia
using Printf, LinearAlgebra, Dates, SparseArrays
```

**Step 3: 替换步骤 4（第 81-111 行）**

将从 `# ── 4. 零通量能谱` 到 `println("  能谱保存: spectrum_Np$(Np).dat")` 的内容替换为：

```julia
# ── 4. CSR 稀疏矩阵能谱（步骤 A 建矩阵 + 步骤 B Lanczos）──
println("\n[4/5] CSR 稀疏矩阵能谱（krylovdim=$(kd)，$(Threads.nthreads()) 线程）...")
println("      内存预估：CSR ~31 GB + Krylov ~$(round(15*kd*5.77e6*16/1e9,digits=0)) GB")
println("      保留各扇区基态波函数，供结构因子直接使用")
flush(stdout)
hops0 = build_hops(lat, t1, t3, 0.0)
t_spec = @elapsed begin
    all_ev, gs_vecs = compute_spectrum_sparse_with_vecs(
        secs, lat, hops0, V1, V2, V3;
        nev=nev, krylovdim=kd, verbose=true)
end
E0 = all_ev[1][2]
@printf("  能谱总耗时: %.1f s (%.2f min)\n", t_spec, t_spec/60)

println("\n  最低 20 个能量（相对基态）:")
for (i,(m,e)) in enumerate(all_ev[1:min(20,end)])
    @printf("    %2d. k=%2d  ΔE=%.8f\n", i, m, e-E0)
end
n15 = min(15, length(all_ev))
gap = length(all_ev) > 15 ? all_ev[16][2]-all_ev[15][2] : NaN
@printf("  前15态散布: %.8f\n", all_ev[n15][2]-E0)
@printf("  15→16 gap:  %.8f\n", gap)

open("spectrum_Np$(Np).dat","w") do f
    println(f,"# k  E-E0  [30sites Np=$Np V1=$V1 V2=$V2 V3=$V3 t'=$t3 CSR]")
    for (m,e) in all_ev
        @printf(f,"%d  %.10f\n", m, e-E0)
    end
end
println("  能谱保存: spectrum_Np$(Np).dat")

# ── 释放 CSR 矩阵（步骤 5 结构因子不需要，释放 ~31 GB）──
# （CSR 变量在 compute_spectrum_sparse_with_vecs 内部，已离开作用域，GC 可回收）
GC.gc()
@printf("  GC 后内存：已释放 CSR 矩阵\n")
flush(stdout)
```

**Step 4: 运行语法检查**

```bash
cd /Users/shajianyu/CMP_manybody/partial_Hall_crystal/.worktrees/csr-sparse
julia --compile=min -e "include(\"case3_30sites_Np12/main.jl\")" 2>&1 | head -20
```

Expected: 仅出现参数/执行输出，无 `ERROR` 或 `LoadError`

**Step 5: Commit**

```bash
git add case3_30sites_Np12/main.jl
git commit -m "feat: switch case3 main.jl to CSR sparse pipeline (kd=120)"
```

---

### Task 5: 端到端本地验证

**Step 1: 创建本地小规模测试入口**

```bash
cd /Users/shajianyu/CMP_manybody/partial_Hall_crystal/.worktrees/csr-sparse
julia test_csr.jl
```

Expected: `✓ ALL PASS`

**Step 2: 运行完整测试**

```bash
julia test_csr.jl 2>&1
```

检查：
- 所有扇区 `|ΔE|_max < 1e-8` ✓
- H Hermitian 检验通过 ✓

**Step 3: Final commit（如有未提交的改动）**

```bash
git status
# 若有遗漏文件：
git add -p
git commit -m "chore: finalize CSR sparse refactor"
```

---

## 内存使用总结

| 阶段 | 峰值内存（Np=12） |
|------|-----------------|
| 建 k-sector | ~7 GB |
| 建 CSR（A 阶段） | ~38 GB（CSR+KSec） |
| 并行 Lanczos（B 阶段） | ~203 GB（CSR+Krylov） |
| 结构因子（GC 后） | ~25 GB（KSec+eigvec） |

## 服务器运行命令

```bash
julia --threads 64 case3_30sites_Np12/main.jl 2>&1 | tee run_Np12.log
```
