# ============================================================
# Case 1 能谱图（t'=0.22, V1=1.6）
# 用法：julia plot.jl（需先跑完 main.jl 生成 spectrum_Np6.dat）
# ============================================================
using Plots, Printf

Np   = 6
Nuc  = 9
n_gs = Nuc ÷ gcd(Np, Nuc)   # = 3

dat = readlines("spectrum_Np$(Np).dat")
ks = Int[]; es = Float64[]
for line in dat
    startswith(line, "#") && continue
    parts = split(strip(line))
    length(parts) < 2 && continue
    push!(ks, parse(Int, parts[1]))
    push!(es, parse(Float64, parts[2]))
end

# x 轴：m1 外层(0:2)，m2 内层(0:2)，共 9 个 k 点
order   = [(m1, m2) for m1 in 0:2 for m2 in 0:2]
klabels = ["($m1,$m2)" for (m1, m2) in order]
new_x   = Dict((m1, m2) => i for (i, (m1, m2)) in enumerate(order))
xs      = [new_x[(k % 3, k ÷ 3)] for k in ks]

# 最低 n_gs 个态标红
gs_set = Set(sortperm(es)[1:n_gs])
colors = [i in gs_set ? :red : :steelblue for i in 1:length(es)]
sizes  = [i in gs_set ? 8 : 5 for i in 1:length(es)]

p = scatter(xs, es;
    marker=:circle, markersize=sizes, markerstrokewidth=0,
    color=colors,
    xlabel="momentum (m₁, m₂)",
    ylabel="E - E₀",
    title="3×6 lattice, Np=$Np,  ν=$(Np÷gcd(Np,Nuc))/$( Nuc÷gcd(Np,Nuc))\n" *
          "t=1, t'=0.22, V₁=1.6, V₂=V₃=0",
    xticks=(1:9, klabels),
    xrotation=45,
    legend=false,
    ylims=(-0.02, 0.40),
    xlims=(0.3, 9.7),
    framestyle=:box,
    size=(520, 520),
    left_margin=8Plots.mm, bottom_margin=12Plots.mm
)

for i in sortperm(es)[1:n_gs]
    annotate!(xs[i], es[i] + 0.006, text(@sprintf("%.4f", es[i]), 7, :red, :center))
end

savefig(p, "spectrum_Np$(Np).pdf")
println("保存：spectrum_Np$(Np).pdf")

# ── 结构因子 N(k) 图 ──
sq_dat = readlines("sq_Np$(Np).dat")
sq_ks = Int[]; sq_vals = Float64[]
for line in sq_dat
    startswith(line, "#") && continue
    parts = split(strip(line))
    length(parts) < 2 && continue
    push!(sq_ks, parse(Int, parts[1]))
    push!(sq_vals, parse(Float64, parts[2]))
end

sq_xs = [new_x[(k % 3, k ÷ 3)] for k in sq_ks]
ord   = sortperm(sq_xs)
sx, sy = sq_xs[ord], sq_vals[ord]

p2 = plot(sx, sy;
    marker=:circle, markersize=6, markerstrokewidth=0,
    color=:steelblue, linewidth=1.5,
    xlabel="momentum (m₁, m₂)",
    ylabel="n(q)",
    title="3×6 lattice, Np=$Np: Structure Factor n(q)\n" *
          "t=1, t'=0.22, V₁=1.6, V₂=V₃=0",
    xticks=(1:9, klabels),
    xrotation=45,
    legend=false,
    xlims=(0.3, 9.7),
    ylims=(0, max(maximum(sq_vals[sq_ks .!= 0]), 0.5) * 1.15),
    framestyle=:box,
    size=(520, 520),
    left_margin=8Plots.mm, bottom_margin=12Plots.mm
)

nonzero_idx = findall(k -> k != 0, sq_ks)
top2 = sort(nonzero_idx, by=i -> -sq_vals[i])[1:min(2, end)]
for i in top2
    annotate!(sq_xs[i], sq_vals[i] + 0.015,
        text(@sprintf("%.4f", sq_vals[i]), 7, :navy, :center))
end

savefig(p2, "sq_Np$(Np).pdf")
println("保存：sq_Np$(Np).pdf")
