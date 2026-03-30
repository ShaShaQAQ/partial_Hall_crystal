# ============================================================
# 合并谱流结果：output_sf/sf_*.jld2 → sf_Np$(Np).dat
# 用法：julia merge_sf.jl
# ============================================================
using JLD2, Printf

# 找所有分段结果
jld2_files = sort(filter(f -> endswith(f, ".jld2"),
    readdir("output_sf", join=true)))
isempty(jld2_files) && error("output_sf/ 中没有 .jld2 文件")

println("找到 $(length(jld2_files)) 个 JLD2 文件:")
for f in jld2_files; println("  $f"); end

# 从第一个文件读取元数据
d0   = load(jld2_files[1])
N_phi     = d0["N_phi"]
phi_grid  = d0["phi_grid"]
Np        = d0["Np"]
V1        = d0["V1"]
V2        = get(d0, "V2", 0.0)
V3        = get(d0, "V3", 0.0)
t1        = d0["t1"]
t3        = d0["t3"]

# 合并所有 sf_data
all_sf = Dict{Int, Vector{Tuple{Int,Float64}}}()
for f in jld2_files
    d = load(f)
    for (k, v) in d["sf_data"]
        all_sf[k] = v
    end
end

n_found = length(all_sf)
@printf("合并后覆盖 φ 点数: %d / %d\n", n_found, N_phi)
if n_found < N_phi
    missing_idx = setdiff(0:N_phi-1, keys(all_sf))
    println("缺失 φ 索引: $missing_idx")
end

# 全局最低能量（用于相对坐标）
E_min = minimum(
    minimum(e for (_, e) in ev)
    for (_, ev) in all_sf if !isempty(ev)
)
@printf("全局最低能量 E_min = %.10f\n", E_min)

# 写 .dat 文件
outfile = "sf_Np$(Np).dat"
open(outfile, "w") do f
    println(f, "# 谱流数据  Np=$Np  V1=$V1 V2=$V2 V3=$V3 t1=$t1 t3=$t3  N_phi=$N_phi")
    println(f, "# phi_idx  phi/2pi  k  E-Emin")
    for phi_idx in sort(collect(keys(all_sf)))
        phi_val = phi_grid[phi_idx + 1]
        for (k, e) in sort(all_sf[phi_idx], by=x->x[2])
            @printf(f, "%d  %.8f  %d  %.10f\n",
                    phi_idx, phi_val / (2π), k, e - E_min)
        end
    end
end
println("保存: $outfile")
