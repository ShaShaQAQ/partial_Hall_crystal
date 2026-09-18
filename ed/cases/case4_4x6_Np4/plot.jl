# ============================================================
# Case 4 能谱图
# 用法：julia plot.jl（需先跑完 main.jl 生成 spectrum_Np4.dat）
# ============================================================
using Plots, Printf

Np   = 4
Nuc  = 12
n_gs = Nuc ÷ Np   # 基态简并度 = 3

dat = readlines("spectrum_Np$(Np).dat")
ks = Int[]; es = Float64[]
for line in dat
    startswith(line, "#") && continue
    parts = split(strip(line))
    length(parts) < 2 && continue
    push!(ks, parse(Int, parts[1]))
    push!(es, parse(Float64, parts[2]))
end

# x轴：按本地 ED 的 (m1,m2) 约定排列，m1 外层，m2 内层
order   = [(m1, m2) for m1 in 0:3 for m2 in 0:2]
klabels = ["($m1,$m2)" for (m1, m2) in order]
new_x   = Dict((m1, m2) => i for (i, (m1, m2)) in enumerate(order))
xs      = [new_x[(k % 4, k ÷ 4)] for k in ks]

# 最低 n_gs 个态标红
gs_set = Set(sortperm(es)[1:n_gs])
colors = [i in gs_set ? :red : :steelblue for i in 1:length(es)]
sizes  = [i in gs_set ? 8 : 5 for i in 1:length(es)]

p = scatter(xs, es;
    marker=:circle, markersize=sizes, markerstrokewidth=0,
    color=colors,
    xlabel="momentum (m₁, m₂)",
    ylabel="E - E₀",
    title="4×6 lattice, Np=$Np,  ν=$(Np÷gcd(Np,Nuc))/$( Nuc÷gcd(Np,Nuc))\n" *
          "t=1, t'=0.2, V₁=1, V₂=V₃=0",
    xticks=(1:12, klabels),
    xrotation=45,
    legend=false,
    ylims=(-0.02, 0.20),
    xlims=(0.3, 12.7),
    framestyle=:box,
    size=(520, 520),
    left_margin=8Plots.mm, bottom_margin=12Plots.mm
)

# 标注基态能量
for i in sortperm(es)[1:n_gs]
    annotate!(xs[i], es[i] + 0.006, text(@sprintf("%.4f", es[i]), 7, :red, :center))
end

savefig(p, "spectrum_Np$(Np).pdf")
println("保存：spectrum_Np$(Np).pdf")

# ── 结构因子 n(k) 图 ──
sq_dat = readlines("sq_Np$(Np).dat")
sq_ks = Int[]; sq_vals = Float64[]
for line in sq_dat
    startswith(line, "#") && continue
    parts = split(strip(line))
    length(parts) < 2 && continue
    push!(sq_ks, parse(Int, parts[1]))
    push!(sq_vals, parse(Float64, parts[2]))
end

sq_xs = [new_x[(k % 4, k ÷ 4)] for k in sq_ks]
ord   = sortperm(sq_xs)
sx, sy = sq_xs[ord], sq_vals[ord]

p2 = plot(sx, sy;
    marker=:circle, markersize=6, markerstrokewidth=0,
    color=:steelblue, linewidth=1.5,
    xlabel="momentum (m₁, m₂)",
    ylabel="n(q)",
    title="4×6 lattice, Np=$Np: Structure Factor n(q)\n" *
          "t=1, t'=0.2, V₁=1, V₂=V₃=0",
    xticks=(1:12, klabels),
    xrotation=45,
    legend=false,
    xlims=(0.3, 12.7),
    ylims=(0, max(maximum(sq_vals), 0.3) * 1.15),
    framestyle=:box,
    size=(520, 520),
    left_margin=8Plots.mm, bottom_margin=12Plots.mm
)

nonzero_idx = findall(k -> k != 0, sq_ks)
top2 = sort(nonzero_idx, by=i -> -sq_vals[i])[1:min(2, end)]
for i in top2
    annotate!(sq_xs[i], sq_vals[i] + 0.008,
        text(@sprintf("%.4f", sq_vals[i]), 7, :navy, :center))
end

savefig(p2, "sq_Np$(Np).pdf")
println("保存：sq_Np$(Np).pdf")

# ── 全扇区结构因子热图 ──
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

# 热图矩阵：行 = x轴位置（sector_k），列 = y轴位置（q）
# 用同一 new_x 映射使两轴均按 (m1,m2) 顺序排列
Z = fill(NaN, Nuc, Nuc)
for (sk, qk, nq) in zip(sec_ks, q_ks, nq_vals)
    xi = new_x[(sk % 4, sk ÷ 4)]
    yi = new_x[(qk % 4, qk ÷ 4)]
    Z[xi, yi] = nq
end

p3 = heatmap(1:Nuc, 1:Nuc, Z';
    xlabel="sector k",
    ylabel="momentum q",
    title="Structure Factor n(q) — all sectors\n4×6, Np=$Np, t=1, t'=0.2, V₁=1",
    xticks=(1:Nuc, klabels),
    yticks=(1:Nuc, klabels),
    xrotation=45,
    yrotation=0,
    color=:viridis,
    aspect_ratio=:equal,
    framestyle=:box,
    size=(680, 640),
    left_margin=14Plots.mm, bottom_margin=14Plots.mm
)

savefig(p3, "sq_all_Np$(Np).pdf")
println("保存：sq_all_Np$(Np).pdf")
