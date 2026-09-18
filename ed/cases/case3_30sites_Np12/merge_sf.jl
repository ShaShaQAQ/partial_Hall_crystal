# ============================================================
# 合并谱流结果：output_sf/sf_*.jld2 → sf_Np$(Np).dat
# 用法：julia merge_sf.jl
# ============================================================
using JLD2, Printf
cd(@__DIR__)

# 找所有分段结果
jld2_files = sort(filter(f -> endswith(f, ".jld2"),
    readdir("output_sf", join=true)))
isempty(jld2_files) && error("output_sf/ 中没有 .jld2 文件")

println("找到 $(length(jld2_files)) 个 JLD2 文件:")
for f in jld2_files; println("  $f"); end

# 从第一个文件读取元数据
d0       = load(jld2_files[1])
N_phi    = d0["N_phi"]
phi_grid = d0["phi_grid"]
Np       = d0["Np"]
V1       = d0["V1"]
V2       = get(d0, "V2", 0.0)
V3       = get(d0, "V3", 0.0)
t1       = d0["t1"]
t3       = d0["t3"]

# 合并所有 sf_data
# 每个文件的 sf_data[phi_idx] 只含该节点负责的扇区 → append 而非覆盖
all_sf = Dict{Int, Vector{Tuple{Int,Float64}}}()
for f in jld2_files
    d = load(f)
    for (phi_idx, ev) in d["sf_data"]
        append!(get!(all_sf, phi_idx, Tuple{Int,Float64}[]), ev)
    end
end
# 每个 φ 点内按能量排序
for ev in values(all_sf)
    sort!(ev, by=x->x[2])
end

# 检查每个 φ 点是否已被所有文件覆盖（每个文件贡献部分扇区）
# 以 φ 点为键，收集来自不同节点的 (k, E) 对
phi_covered = Set(keys(all_sf))
n_found = length(phi_covered)
@printf("合并后覆盖 φ 点数: %d / %d\n", n_found, N_phi)
if n_found < N_phi
    missing_idx = setdiff(0:N_phi-1, phi_covered)
    println("缺失 φ 索引: $missing_idx")
end

# 全局最低能量（用于相对坐标）
E_min = minimum(
    minimum(e for (_, e) in ev)
    for (_, ev) in all_sf if !isempty(ev)
)
@printf("全局最低能量 E_min = %.10f\n", E_min)

# 写 .dat 文件
outfile = "output_sf/sf_Np$(Np).dat"
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
