# CSR 预存稀疏矩阵重构设计

**分支:** `feature/csr-sparse`
**目标系统:** 30 sites，Np=12（填充 2/5），服务器 256 GB 单节点 64 核

## 问题

当前 `Hv!` 每次 Lanczos 迭代都实时计算哈希查找（`fock2rep` Dict 随机内存访问）。
对 C(30,12)=86.5M 的 Hilbert 空间，这是主要性能瓶颈。

## 解法

一次性将每个动量扇区的哈密顿量构建为 `SparseMatrixCSC`（CSR 格式），
之后每次 `H*v` 变为标准稀疏矩阵-向量乘法（顺序内存访问，BLAS 优化）。

## 内存约束（256 GB）

| 组件 | Np=12 | Np=13 |
|------|-------|-------|
| KSector (fock2rep+orbit) | ~6 GB | ~8 GB |
| 15 个 CSR 矩阵 | ~31 GB | ~43 GB |
| Krylov 向量（15×kd×dim×16B）| 166 GB (kd=120) | 153 GB (kd=80) |
| **合计** | **~203 GB** | **~204 GB** |

`krylovdim` 从 200 降至 120（Np=12）/ 80（Np=13）。
`nev=8` 时理论最低 krylovdim≈26，原来 200 是保守值，CSR 加速后更多重启仍更快。

## 新流程

```
[1] 建晶格             (不变)
[2] 生成基矢           (不变)
[3] 建 k-sector        (不变，fock2rep 仍需构建)
[4] NEW: Threads@threads 并行建 15 个 H_csr (build_sparse_H)
[5] NEW: Threads@threads 并行 Lanczos (mul!(w, H_csr, v), krylovdim=120)
[6] NEW: H_csrs = nothing; GC.gc()  → 释放 ~31 GB
[7] 结构因子           (不变，仍用 fock2rep + eigenvectors)
```

## 改动文件

| 文件 | 改动 |
|------|------|
| `shared/hamiltonian.jl` | 新增 `build_sparse_H()` |
| `shared/solver.jl` | 新增 `compute_spectrum_sparse_with_vecs()` |
| `case3_30sites_Np12/main.jl` | 步骤 4/5/6 替换为 CSR 流程 |

不改动：`ksector.jl`、`basis.jl`、`structure_factor.jl`、`lattice.jl`、`hoppings.jl`

## build_sparse_H 逻辑

```
对每个代表元 i：
  push! (i, i, diag_H[i])                     ← 对角势能（内联 precompute_diag_H）
  对 orbit_data[i] 每个 (c_α, F)：
    对每个 hopping (tgt, src, tval)：
      Pauli 检验 → hop(F,tgt,src) → fock2rep 查 j
      push! (j, i, c_α * ni * tval * sgn * c_conj_norm)
sparse(Is, Js, Vs, Nrep, Nrep)               ← Julia 自动合并重复项
```

## 备选扩展（暂不实现）

若单节点内存不足，可将 15 个扇区分配到多个计算节点，
各节点独立运行部分扇区后合并能谱。
