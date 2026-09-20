using Test
using SparseArrays

function spectrum_driver_output(driver, expression, arguments=String[])
    project = dirname(Base.active_project())
    command = `$(Base.julia_cmd()) --startup-file=no --threads=1
        --project=$project -e $expression $driver $arguments`
    return read(command, String)
end

@testset "30-site spectrum driver loads" begin
    driver = joinpath(
        @__DIR__, "..", "cases", "case3_30sites_Np12",
        "run_spectrum.jl")
    @test isfile(driver)
    if isfile(driver)
        include(driver)
        @test isdefined(Main, :run_large_spectrum)
        expected_np = parse(
            Int, parse_spectrum_argument("--Np", "12"))
        @test SPECTRUM_NP == expected_np
        @test SPECTRUM_V1 == 10.0
        @test SPECTRUM_V2 == 2.0
        @test SPECTRUM_V3 == 2.0
        @test SPECTRUM_LATTICE == "corrected"
        @test SPECTRUM_SECTOR_END == SPECTRUM_SECTOR_START
        @test endswith(spectrum_partial_path(5), "partial_5.jld2")
        @test isdefined(Main, :spectrum_lanczos_matrix)
        if isdefined(Main, :spectrum_lanczos_matrix)
            matrix = sparse(ComplexF64[1 0; 0 2])
            threaded = spectrum_lanczos_matrix(matrix)
            @test threaded isa ThreadedCSR
            product = zeros(ComplexF64, 2)
            mul!(product, threaded, ComplexF64[1, 2])
            @test product == ComplexF64[1, 4]
            @test applicable(
                lanczos_sparse_sectors, KSector[], ThreadedCSR[])
        end
    end
end

@testset "spectrum driver propagates Np in a fresh process" begin
    driver = joinpath(
        @__DIR__, "..", "cases", "case3_30sites_Np12",
        "run_spectrum.jl")
    expression =
        "driver = popfirst!(ARGS); include(driver); print(SPECTRUM_NP)"
    output = spectrum_driver_output(
        driver, expression, ["--Np", "13"])
    @test strip(output) == "13"
end

@testset "single-sector threaded spectrum Lanczos" begin
    lattice = RectLat4x6()
    sector = build_ksector(gen_basis(lattice.Ns, 2), lattice, 0)
    matrix = build_sparse_H(
        sector, lattice, build_hops(lattice, 1.0, 0.2, 0.0),
        1.0, 0.0, 0.0)
    threaded = spectrum_lanczos_matrix(matrix)
    levels, states = lanczos_sparse_sectors(
        [sector], [threaded]; nev=2, krylovdim=10)
    @test length(levels) == 2
    @test haskey(states, 0)
    @test norm(states[0]) ≈ 1.0 atol=1e-12
    product = similar(states[0])
    mul!(product, threaded, states[0])
    @test norm(product .- levels[1][2] .* states[0]) < 1e-8
end
