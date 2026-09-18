@testset "analytic Bloch kx derivatives" begin
    t1, t3 = 1.0, 0.2
    delta1 = 1e-6
    delta2 = 1e-4
    for k in ([0.173, -0.291], [1.137, 0.419], [-0.733, 1.271])
        data = get_Hk_x_derivatives(k, t1, t3)
        hp1 = get_Hk(k .+ [delta1, 0.0], t1, t3)
        hm1 = get_Hk(k .- [delta1, 0.0], t1, t3)
        hp2 = get_Hk(k .+ [delta2, 0.0], t1, t3)
        hm2 = get_Hk(k .- [delta2, 0.0], t1, t3)
        d1_fd = (hp1 - hm1) / (2delta1)
        d2_fd = (hp2 - 2data.Hk + hm2) / delta2^2

        @test data.Hk == get_Hk(k, t1, t3)
        @test norm(data.dHdkx - d1_fd) /
              max(norm(data.dHdkx), 1.0) < 1e-9
        @test norm(data.d2Hdkx2 - d2_fd) /
              max(norm(data.d2Hdkx2), 1.0) < 2e-7
    end
end

function single_particle_matrix(lat, hops)
    operator = zeros(ComplexF64, lat.Ns, lat.Ns)
    for (target, source, amplitude) in hops
        operator[target, source] += amplitude
    end
    return operator
end

function torus_matrix_from_bloch(lat, bloch_matrix)
    operator = zeros(ComplexF64, lat.Ns, lat.Ns)
    a1_uc = collect(lat.a1)
    a2_uc = 2 .* collect(lat.a2)
    sublat = sublat_pos(lat)
    Nk = length(lat.kpoints)
    for (source, (ix, iy)) in enumerate(lat.sites)
        source_orb = mod(iy, 2) + 1
        source_cell_y = div(iy - (source_orb - 1), 2)
        for (n1_uc, n2_uc) in lat.uc_trans
            Rcart = n1_uc .* a1_uc .+ n2_uc .* a2_uc
            for target_orb in 1:2
                delta = sublat[target_orb] .- sublat[source_orb]
                amplitude = 0.0 + 0.0im
                for k in lat.kpoints
                    phase = exp(-1im * dot(k, Rcart .+ delta))
                    amplitude +=
                        phase * bloch_matrix(k)[target_orb, source_orb]
                end
                amplitude /= Nk
                tix = ix + n1_uc
                tiy = 2(source_cell_y + n2_uc) + target_orb - 1
                target_site = lat.site_idx[canon_prim(lat, tix, tiy)]
                operator[target_site, source] += amplitude
            end
        end
    end
    return operator
end

@testset "physical finite-torus Ax operators" begin
    lat = RectLat4x6()
    t1, t3 = 1.0, 0.2
    hops0, jx_hops, kxx_hops = build_hops_x_derivatives(lat, t1, t3)
    H0 = single_particle_matrix(lat, hops0)
    Jx = single_particle_matrix(lat, jx_hops)
    Kxx = single_particle_matrix(lat, kxx_hops)

    data(k) = get_Hk_x_derivatives(k, t1, t3)
    H0_ref = torus_matrix_from_bloch(lat, k -> data(k).Hk)
    Jx_ref = torus_matrix_from_bloch(lat, k -> data(k).dHdkx)
    Kxx_ref = torus_matrix_from_bloch(lat, k -> data(k).d2Hdkx2)
    relerr(actual, reference) =
        norm(actual - reference) / max(norm(reference), 1.0)

    @test relerr(H0, H0_ref) < 1e-12
    @test relerr(Jx, Jx_ref) < 1e-12
    @test relerr(Kxx, Kxx_ref) < 1e-12

    Ax = 0.071
    HAx = single_particle_matrix(lat, build_hops_Ax(lat, t1, t3, Ax))
    HAx_ref = torus_matrix_from_bloch(
        lat, k -> get_Hk(k .+ [Ax, 0.0], t1, t3))
    @test relerr(HAx, HAx_ref) < 1e-12

    delta1 = 1e-6
    delta2 = 1e-4
    Hp1 = single_particle_matrix(
        lat, build_hops_Ax(lat, t1, t3, delta1))
    Hm1 = single_particle_matrix(
        lat, build_hops_Ax(lat, t1, t3, -delta1))
    Hp2 = single_particle_matrix(
        lat, build_hops_Ax(lat, t1, t3, delta2))
    Hm2 = single_particle_matrix(
        lat, build_hops_Ax(lat, t1, t3, -delta2))
    @test relerr((Hp1 - Hm1) / (2delta1), Jx) < 1e-8
    @test relerr((Hp2 - 2H0 + Hm2) / delta2^2, Kxx) < 1e-5

    Hlegacy = single_particle_matrix(lat, build_hops(lat, t1, t3, 0.0))
    @test relerr(H0, Hlegacy) < 1e-12
    for operator in (H0, Jx, Kxx, HAx)
        @test norm(operator - operator') / max(norm(operator), 1.0) < 1e-12
    end
end
