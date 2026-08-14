using Test
using LinearAlgebra
using SparseArrays
using InfiniteCylinderDMRG

@testset "fermionic ring translation" begin
    for Ny in 2:4
        dimension = 1 << (2Ny)
        translation = fermionic_ring_translation(Ny)

        @test issparse(translation)
        @test size(translation) == (dimension, dimension)
        @test nnz(translation) == dimension
        @test translation' * translation ≈ I atol=1e-14 rtol=0
        @test translation^Ny ≈ I atol=1e-14 rtol=0

        for bits in 0:(dimension - 1)
            target, sign = translated_fock_basis(bits, Ny)
            @test translation[target + 1, bits + 1] == sign
            @test count(!iszero, translation[:, bits + 1]) == 1
        end
    end

    @test translated_fock_basis(0b0101, 2) == (0b0101, -1)
    @test translated_fock_basis(0b1100, 2) == (0b0011, 1)
    @test translated_fock_basis(0b1111, 2) == (0b1111, 1)
end

@testset "fermionic translation action and geometry gate" begin
    Ny = 2
    bits = 0b0101
    coefficients = zeros(ComplexF64, 1 << (2Ny))
    coefficients[bits + 1] = 1 + 2im
    translated = apply_fermionic_ring_translation(coefficients, Ny)
    target, sign = translated_fock_basis(bits, Ny)

    @test translated[target + 1] == sign * coefficients[bits + 1]
    @test count(!iszero, translated) == 1

    paper = InfiniteCylinderConfig(;
        geometry=:paper_straight,
        Ny,
        x_period=3,
        filling_num=1,
        filling_den=3,
    )
    legacy = InfiniteCylinderConfig(; Ly=4, x_period=3, filling_num=1, filling_den=3)
    @test fermionic_ring_translation(paper) == fermionic_ring_translation(Ny)
    @test apply_fermionic_ring_translation(coefficients, paper) == translated
    @test_throws ArgumentError fermionic_ring_translation(legacy)
    @test_throws ArgumentError apply_fermionic_ring_translation(coefficients, legacy)
end

@testset "fermionic translation input validation" begin
    @test_throws ArgumentError translated_fock_basis(-1, 2)
    @test_throws ArgumentError translated_fock_basis(1 << 4, 2)
    @test_throws ArgumentError translated_fock_basis(0, 0)
    @test_throws ArgumentError translated_fock_basis(0, true)
    @test_throws DimensionMismatch apply_fermionic_ring_translation(zeros(15), 2)
end
