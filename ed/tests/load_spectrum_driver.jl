using Test
using SparseArrays

@testset "30-site spectrum driver loads" begin
    driver = joinpath(
        @__DIR__, "..", "cases", "case3_30sites_Np12",
        "run_spectrum.jl")
    @test isfile(driver)
    if isfile(driver)
        include(driver)
        @test isdefined(Main, :run_large_spectrum)
        @test SPECTRUM_NP == 12
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
