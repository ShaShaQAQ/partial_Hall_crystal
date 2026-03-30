# ============================================================
# 谱流图：sf_Np$(Np).dat → sf_Np$(Np).pdf
# 用法：julia plot_sf.jl
# ============================================================
using Plots, Printf

Np = @isdefined(Np) ? Np : 12
V1 = @isdefined(V1) ? V1 : 10.0
V2 = @isdefined(V2) ? V2 : 2.0
V3 = @isdefined(V3) ? V3 : 2.0
t1 = @isdefined(t1) ? t1 : 1.0
t3 = @isdefined(t3) ? t3 : 0.2
n_show = @isdefined(n_show) ? n_show : 2   # 每扇区保留最低 n_show 个态

dat = readlines("sf_Np$(Np).dat")

phi_idx_all = Int[]
phi_all     = Float64[]
k_all       = Int[]
e_all       = Float64[]

for line in dat
    startswith(line, "#") && continue
    parts = split(strip(line))
    length(parts) < 4 && continue
    push!(phi_idx_all, parse(Int,     parts[1]))
    push!(phi_all,     parse(Float64, parts[2]))
    push!(k_all,       parse(Int,     parts[3]))
    push!(e_all,       parse(Float64, parts[4]))
end

N_phi = maximum(phi_idx_all) + 1
Nuc   = maximum(k_all) + 1

# ── 每个 (phi_idx, k) 扇区保留最低 n_show 个态 ──
keep = fill(false, length(e_all))
from_pk = Dict{Tuple{Int,Int}, Vector{Int}}()
for (i, (pi_, k)) in enumerate(zip(phi_idx_all, k_all))
    push!(get!(from_pk, (pi_, k), Int[]), i)
end
for (_, idxs) in from_pk
    sorted = sort(idxs, by=i->e_all[i])
    for i in sorted[1:min(n_show, length(sorted))]
        keep[i] = true
    end
end

phi_plot = phi_all[keep]
k_plot   = k_all[keep]
e_plot   = e_all[keep]

# ── 着色：Nuc 个 k 扇区用不同颜色 ──
# 15 种颜色，循环使用 tab20 / 手动序列
palette_colors = [
    :steelblue, :crimson, :forestgreen, :darkorange, :purple,
    :brown, :deeppink, :teal, :goldenrod, :navy,
    :olive, :coral, :slategray, :mediumorchid, :darkcyan
]
colors = [palette_colors[mod(k, length(palette_colors)) + 1] for k in k_plot]

# y 轴范围：覆盖每个 phi 点的第 16 态（基态流中有 15 态简并群）
# 用最低 20 个态的最高能量 + 5% 留白
e_sorted_unique = sort(unique(round.(e_plot, digits=8)))
ymax = length(e_sorted_unique) >= 20 ? e_sorted_unique[20] * 1.1 + 0.05 :
       maximum(e_plot) * 1.1 + 0.05

p = scatter(phi_plot, e_plot;
    color      = colors,
    markersize = 3,
    markerstrokewidth = 0,
    alpha      = 0.8,
    xlabel     = "φ_y / 2π",
    ylabel     = "E - E₀",
    title      = "Spectral Flow — Tilted 30-site (4×4-1), Np=$Np (ν=2/5)\n" *
                 "t=$t1, t'=$t3, V₁=$V1, V₂=$V2, V₃=$V3",
    legend     = false,
    xlims      = (-0.04, 2.04),
    ylims      = (-0.02, ymax),
    framestyle = :box,
    size       = (900, 600),
    left_margin  = 8Plots.mm,
    bottom_margin = 10Plots.mm
)

savefig(p, "sf_Np$(Np).pdf")
println("保存: sf_Np$(Np).pdf")
