@testset "uniform Ax hopping derivatives" begin
    lat = RectLat4x6()
    t1, t3 = 1.0, 0.2
    hops0, jx_hops, kxx_hops = build_hops_x_derivatives(lat, t1, t3)

    asdict(hops) = Dict((tgt, src) => amp for (tgt, src, amp) in hops)
    @test asdict(hops0) == asdict(build_hops(lat, t1, t3, 0.0))

    delta = 1e-5
    hp = asdict(build_hops_Ax(lat, t1, t3, delta))
    hm = asdict(build_hops_Ax(lat, t1, t3, -delta))
    h0 = asdict(hops0)
    jx = asdict(jx_hops)
    kxx = asdict(kxx_hops)

    for key in keys(h0)
        @test (hp[key] - hm[key]) / (2delta) ≈ jx[key] atol=1e-9 rtol=1e-9
        @test (hp[key] - 2h0[key] + hm[key]) / delta^2 ≈ kxx[key] atol=2e-6 rtol=2e-6
    end

    basis = gen_basis(lat.Ns, 1)
    sec = build_ksector(basis, lat, 2)
    function materialize(hops)
        n = length(sec.reps)
        operator = zeros(ComplexF64, n, n)
        input = zeros(ComplexF64, n)
        diagonal = zeros(Float64, n)
        for column in 1:n
            input[column] = 1
            Hv!(view(operator, :, column), input, sec, lat, hops,
                0.0, 0.0, 0.0, diagonal)
            input[column] = 0
        end
        return operator
    end
    H = materialize(hops0)
    Jx = materialize(jx_hops)
    Kxx = materialize(kxx_hops)
    @test norm(H - H') / max(norm(H), 1.0) < 1e-12
    @test norm(Jx - Jx') / max(norm(Jx), 1.0) < 1e-12
    @test norm(Kxx - Kxx') / max(norm(Kxx), 1.0) < 1e-12
end
