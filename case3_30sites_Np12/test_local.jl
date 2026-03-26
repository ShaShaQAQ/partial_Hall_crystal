# 本地测试：同一晶格，Np=5，验证代码正确性
# 运行：julia --threads 4 test_local.jl
include("../shared/lattice.jl")
include("../shared/hoppings.jl")
include("../shared/basis.jl")
include("../shared/ksector.jl")
include("../shared/hamiltonian.jl")
include("../shared/solver.jl")
include("../shared/structure_factor.jl")

using Printf, LinearAlgebra, Dates

BLAS.set_num_threads(1)

t1 = 1.0; t3 = 0.2; V1 = 1.0; V2 = 0.0; V3 = 0.0
Np = 5; nev = 8; kd = 100

println("="^60)
println("本地测试 — 30 sites, Np=$Np, C(30,$Np)=$(binomial(30,Np))")
println("t=$t1  t'=$t3  V1=$V1  V2=$V2  V3=$V3")
println("Julia threads: $(Threads.nthreads())")
println("开始: ", now())
println("="^60)

t_total = time()

print("[1/5] 构建晶格... "); flush(stdout)
lat = TiltedLat30()
println("done  Ns=$(lat.Ns)  Nuc=$(lat.Nuc)")

print("[2/5] 生成基矢... "); flush(stdout)
basis = gen_basis(lat.Ns, Np)
println("done  $(length(basis)) 个 Fock 态")

println("[3/5] 构建动量扇区...")
secs = KSector[]
for (i, m) in enumerate(lat.ktab)
    local t = @elapsed sec = build_ksector(basis, lat, m)
    push!(secs, sec)
    @printf("  k=%2d (%2d/%d): dim=%5d  %.2f s\n", m, i, lat.Nuc, length(sec.reps), t)
end
total_reps = sum(length(s.reps) for s in secs)
@printf("  总代表元: %d / %d\n", total_reps, length(basis))

println("\n[4/5] 能谱（krylovdim=$(kd)，保留各扇区基态波函数）...")
hops0 = build_hops(lat, t1, t3, 0.0)
t_spec = @elapsed begin
    all_ev, gs_vecs = compute_spectrum_with_vecs(secs, lat, hops0, V1, V2, V3;
                                                  nev=nev, krylovdim=kd, verbose=true)
end
E0 = all_ev[1][2]
@printf("  能谱耗时: %.2f s\n", t_spec)

println("\n  最低 20 个能量（相对基态）:")
for (i,(m,e)) in enumerate(all_ev[1:min(20,end)])
    @printf("    %2d. k=%2d  ΔE=%.8f\n", i, m, e-E0)
end

open("spectrum_Np$(Np).dat","w") do f
    println(f,"# k  E-E0  [30sites Np=$Np V1=$V1 V2=$V2 V3=$V3 t'=$t3]")
    for (m,e) in all_ev
        @printf(f,"%d  %.10f\n", m, e-E0)
    end
end
println("  能谱保存: spectrum_Np$(Np).dat")

println("\n[5/5] 结构因子（直接使用已保存波函数）...")
gs_k   = all_ev[1][1]
gs_sec = secs[findfirst(s -> s.m == gs_k, secs)]

# ── 全局基态 ──
gs_vec = gs_vecs[gs_k]
results_ord, results_srt, norm2, Np_check =
    structure_factor_survey(gs_sec, gs_vec, lat, Np)
print_sq_report(results_srt, norm2, Np_check, Np, lat.Ns, gs_k)

open("sq_Np$(Np).dat","w") do f
    println(f,"# ik  N(k)  kx  ky  [30sites Np=$Np V1=$V1 t'=$t3]")
    println(f,"# 基态扇区 k=$gs_k  E0=$(round(E0,digits=8))")
    for (m, sq, kx, ky) in results_ord
        @printf(f,"%d  %.6f  %.6f  %.6f\n", m, sq, kx, ky)
    end
end
println("  结构因子保存: sq_Np$(Np).dat")

# ── 全扇区 ──
open("sq_all_Np$(Np).dat","w") do f
    println(f,"# sector_k  q  N(q)  kx  ky  [30sites Np=$Np V1=$V1 t'=$t3]")
    for sec in secs
        haskey(gs_vecs, sec.m) || continue
        res_ord, _, _, _ = structure_factor_survey(sec, gs_vecs[sec.m], lat, Np)
        for (m, sq, kx, ky) in res_ord
            @printf(f,"%d  %d  %.6f  %.6f  %.6f\n", sec.m, m, sq, kx, ky)
        end
    end
end
println("  全扇区结构因子保存: sq_all_Np$(Np).dat")

@printf("\n总耗时: %.2f s\n", time()-t_total)
println("完成!")
