# ============================================================
# Case 3 能谱图 + 结构因子图
# 用法：julia plot.jl（需先跑完 main.jl 生成 .dat 文件）
# ============================================================
using Plots, Printf

Np  = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 12
Nuc = 15
# 预期简并度：Np=12(ν=2/5)→15，其余默认标最低1个
n_gs = Np == 12 ? 15 : 1

# ── 1. 能谱图 ──
dat = readlines("spectrum_Np$(Np).dat")
ks = Int[]; es = Float64[]
for line in dat
    startswith(line, "#") && continue
    parts = split(strip(line))
    length(parts) < 2 && continue
    push!(ks, parse(Int, parts[1]))
    push!(es, parse(Float64, parts[2]))
end

xs = ks .+ 1    # k=0..14 → x=1..15

# 最低 n_gs 个态标红
gs_set = Set(sortperm(es)[1:min(n_gs, length(es))])
colors = [i in gs_set ? :red : :steelblue for i in 1:length(es)]
sizes  = [i in gs_set ? 8 : 5 for i in 1:length(es)]

p = scatter(xs, es;
    marker=:circle, markersize=sizes, markerstrokewidth=0,
    color=colors,
    xlabel="momentum k",
    ylabel="E - E₀",
    title="Tilted 30-site (4×4-1), Np=$Np,  ν=2/5\n" *
          "t=1, t'=0.2, V₁=1, V₂=V₃=0",
    xticks=(1:15, string.(0:14)),
    legend=false,
    ylims=(-0.02, 0.40),
    xlims=(0.3, 15.7),
    framestyle=:box,
    size=(600, 520),
    left_margin=8Plots.mm, bottom_margin=10Plots.mm
)

for i in sortperm(es)[1:min(n_gs, length(es))]
    annotate!(xs[i], es[i] + 0.008, text(@sprintf("%.4f", es[i]), 6, :red, :center))
end

savefig(p, "spectrum_Np$(Np).pdf")
println("保存：spectrum_Np$(Np).pdf")

# ── 2. 全局基态结构因子折线图 ──
sq_dat = readlines("sq_Np$(Np).dat")
sq_ks = Int[]; sq_vals = Float64[]
for line in sq_dat
    startswith(line, "#") && continue
    parts = split(strip(line))
    length(parts) < 2 && continue
    push!(sq_ks, parse(Int, parts[1]))
    push!(sq_vals, parse(Float64, parts[2]))
end

sq_xs = sq_ks .+ 1
ord   = sortperm(sq_xs)
sx, sy = sq_xs[ord], sq_vals[ord]

p2 = plot(sx, sy;
    marker=:circle, markersize=6, markerstrokewidth=0,
    color=:steelblue, linewidth=1.5,
    xlabel="momentum q",
    ylabel="n(q)",
    title="Tilted 30-site (4×4-1), Np=$Np: Structure Factor n(q)  [global GS]\n" *
          "t=1, t'=0.2, V₁=1, V₂=V₃=0",
    xticks=(1:15, string.(0:14)),
    legend=false,
    xlims=(0.3, 15.7),
    ylims=(0, max(maximum(sq_vals[sq_ks .!= 0]), 0.3) * 1.15),
    framestyle=:box,
    size=(600, 520),
    left_margin=8Plots.mm, bottom_margin=10Plots.mm
)

nonzero_idx = findall(k -> k != 0, sq_ks)
top3 = sort(nonzero_idx, by=i -> -sq_vals[i])[1:min(3, end)]
for i in top3
    annotate!(sq_xs[i], sq_vals[i] + 0.015,
        text(@sprintf("%.4f", sq_vals[i]), 7, :navy, :center))
end

savefig(p2, "sq_Np$(Np).pdf")
println("保存：sq_Np$(Np).pdf")

# ── 3. 全扇区结构因子热图 ──
# sq_all 文件格式：sector_k  q  N(q)  kx  ky
isfile("sq_all_Np$(Np).dat") || begin
    println("sq_all_Np$(Np).dat 不存在，跳过热图")
    exit(0)
end

sq_all_dat = readlines("sq_all_Np$(Np).dat")
sec_ks = Int[]; q_ks = Int[]; nq_vals = Float64[]
for line in sq_all_dat
    startswith(line, "#") && continue
    parts = split(strip(line))
    length(parts) < 3 && continue
    push!(sec_ks, parse(Int, parts[1]))
    push!(q_ks,   parse(Int, parts[2]))
    push!(nq_vals, parse(Float64, parts[3]))
end

# 整理成矩阵 Z[sector_k+1, q+1]
Z = fill(NaN, Nuc, Nuc)
for (sk, qk, nq) in zip(sec_ks, q_ks, nq_vals)
    Z[sk+1, qk+1] = nq
end

p3 = heatmap(0:Nuc-1, 0:Nuc-1, Z;
    xlabel="momentum q",
    ylabel="sector k",
    title="Structure Factor n(q) — all sectors\nNp=$Np, t=1, t'=0.2, V₁=1",
    xticks=0:Nuc-1,
    yticks=0:Nuc-1,
    color=:viridis,
    aspect_ratio=:equal,
    framestyle=:box,
    size=(620, 580),
    left_margin=10Plots.mm, bottom_margin=10Plots.mm
)

savefig(p3, "sq_all_Np$(Np).pdf")
println("保存：sq_all_Np$(Np).pdf")
