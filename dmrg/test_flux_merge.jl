include("merge_flux_pump.jl")
using .FluxPumpMerge
using Test

function write_fake_density(path, values)
    open(path, "w") do io
        println(io, "# site x y rx ry density")
        for (i, n) in enumerate(values)
            println(io, "$i 0 0 0 0 $n")
        end
    end
end

@testset "flux pump merge" begin
    root = mktempdir()
    mkpath(joinpath(root, "phi_000"))
    mkpath(joinpath(root, "phi_001"))
    write_fake_density(joinpath(root, "phi_000", "density.dat"), [0.5, 0.5, 0.5, 0.5])
    write_fake_density(joinpath(root, "phi_001", "density.dat"), [0.6, 0.7, 0.4, 0.3])
    FluxPumpMerge.write_summary(joinpath(root, "phi_000", "summary.dat"), 0.0, -1.0)
    FluxPumpMerge.write_summary(joinpath(root, "phi_001", "summary.dat"), 2π, -0.9)

    rows = merge_flux_data(root; Lx=2, Ly=2, steps=1)
    @test length(rows) == 2
    @test rows[1].mid ≈ 0.0 atol=1e-12
    @test rows[2].mid ≈ 0.3 atol=1e-12
    @test rows[2].right_edge ≈ 0.0 atol=1e-12
    @test isfile(joinpath(root, "pumping.dat"))
    @test isfile(joinpath(root, "cumulative_phi_001.dat"))
end
