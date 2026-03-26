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
    global max_err = max(max_err, err)
    @printf("  k=%2d  dim=%4d  nnz=%6d  |ΔE|_max=%.2e  %s\n",
            sec.m, length(sec.reps), nnz(H_csr), err,
            err < 1e-8 ? "✓" : "✗ FAIL")
end

println()
println(max_err < 1e-8 ? "✓ ALL PASS  max_err=$(max_err)" :
                         "✗ FAIL     max_err=$(max_err)")
