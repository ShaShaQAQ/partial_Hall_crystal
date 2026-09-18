@testset "projected response source" begin
    J = ComplexF64[2 1-im 0; 1+im 0 2; 0 2 0]
    g = ComplexF64[1, 0, 0]
    f = projected_response_source(J, g)
    @test dot(g, f) ≈ 0.0 atol=1e-14
    @test f ≈ ComplexF64[0, 1+im, 0]
    @test_throws ArgumentError projected_response_source(
        Matrix{ComplexF64}(I, 3, 3), g)
end

@testset "response Lanczos continued fraction" begin
    H = ComplexF64[0.0 0.2 0.0 0.0;
                   0.2 1.0 0.3 0.0;
                   0.0 0.3 2.0 0.4;
                   0.0 0.0 0.4 4.0]
    E0 = 0.0
    f = ComplexF64[0, 1, 1im, 0]
    kernel = response_lanczos(H, E0, f; mmax=4, breakdown_tol=1e-14)
    @test length(kernel.alpha) == 4
    @test length(kernel.beta) == 3
    @test kernel.source_norm2 ≈ real(dot(f, f))

    for z in (0.3 + 0.2im, 1.4 + 0.1im, -0.7 - 0.15im)
        direct = dot(f, (z * I - (H - E0 * I)) \ f)
        @test lanczos_resolvent(kernel, z) ≈ direct atol=1e-11 rtol=1e-11
    end

    closed = response_lanczos(Diagonal(ComplexF64[0, 1, 2]), 0.0,
                              ComplexF64[0, 1, 0]; mmax=3)
    @test closed.breakdown
    @test length(closed.alpha) == 1
    @test isempty(closed.beta)
end

@testset "exact spectral reference and Kubo decomposition" begin
    H = ComplexF64[0 0 0; 0 1 0; 0 0 3]
    J = ComplexF64[0 1-im 2; 1+im 0 0; 2 0 0]
    K = ComplexF64[5 0 0; 0 0 0; 0 0 0]
    eig = eigen(Hermitian(H))
    Eg = eig.values[1]
    g = eig.vectors[:, 1]
    f = projected_response_source(J, g)
    data = exact_spectral_data(eig.values, eig.vectors, Eg, f)

    @test data.deltas ≈ [1.0, 3.0]
    @test data.weights ≈ [2.0, 4.0]
    for z in (0.5 + 0.2im, -0.8 - 0.1im)
        direct = dot(f, (z * I - (H - Eg * I)) \ f)
        @test exact_resolvent(data, z) ≈ direct atol=1e-12
    end

    Kexp = real(dot(g, K * g))
    G(z) = exact_resolvent(data, z)
    parts = conductivity_parts(G, Kexp, 7.5, 0.8, 0.065)
    @test parts.total ≈ parts.drude + parts.regular atol=1e-13
    @test parts.drude_weight ≈ Kexp - 2sum(data.weights ./ data.deltas)

    kernel = response_lanczos(H, Eg, f; mmax=3)
    @test lanczos_drude_weight(kernel, Kexp) ≈ parts.drude_weight atol=1e-12
end

@testset "system area" begin
    @test system_area(RectLat4x6()) ≈ 12sqrt(3) atol=1e-12
end
