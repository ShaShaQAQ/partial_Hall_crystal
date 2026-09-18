@testset "cached fermion translations" begin
    available = all(isdefined(Main, name) for name in
                    (:translation_site_maps, :translation_inversion_masks,
                     :translate_fock_cached))
    @test available
    if available
        lat = RectLat4x6()
        maps = translation_site_maps(lat)
        masks = translation_inversion_masks(maps, lat.Ns)
        states = Int64[
            (Int64(1) << 0) | (Int64(1) << 7) | (Int64(1) << 13),
            (Int64(1) << 1) | (Int64(1) << 2) | (Int64(1) << 14),
            (Int64(1) << 3) | (Int64(1) << 8) | (Int64(1) << 17),
        ]
        for state in states, ti in eachindex(lat.uc_trans)
            reference = translate_fock(
                state, lat, lat.Tnx[ti], lat.Tny[ti])
            cached = translate_fock_cached(
                state, maps, masks, ti, lat.Ns)
            @test cached == reference
        end

        basis = gen_basis(lat.Ns, 2)
        reference_sector = build_ksector(basis, lat, 2)
        cached_sector = build_ksector(basis, lat, 2; trans_maps=maps)
        @test cached_sector.reps == reference_sector.reps
        @test cached_sector.norms == reference_sector.norms
        @test cached_sector.fock2rep == reference_sector.fock2rep
    end
end

@testset "legacy tilted-30 saved-state convention" begin
    available = isdefined(Main, :LegacyTiltedLat30)
    @test available
    if available
        corrected = TiltedLat30()
        legacy = LegacyTiltedLat30()
        expected_step = (
            4 .* collect(legacy.b1) .- collect(legacy.b2) ./ 2
        ) ./ 15
        @test legacy.kpoints[2] ≈ expected_step atol=1e-14
        @test legacy.kpoints[2] != corrected.kpoints[2]
        @test legacy.phase_table == corrected.phase_table
        @test length(build_hops(legacy, 1.0, 0.2, 0.0)) == 300
        _, current_hops, diamagnetic_hops =
            build_hops_x_derivatives(legacy, 1.0, 0.2)
        @test length(current_hops) == 300
        @test length(diamagnetic_hops) == 300
    end
end

@testset "threaded matrix-free Hamiltonian" begin
    available = all(isdefined(Main, name) for name in
                    (:threaded_Hv_workspace, :Hv_threaded!))
    @test available
    if available
        lat = RectLat4x6()
        sector = build_ksector(gen_basis(lat.Ns, 2), lat, 2)
        hops = build_hops(lat, 1.0, 0.2, 0.0)
        diagonal = precompute_diag_H(sector, lat, 1.0, 0.0, 0.0)
        source = normalize(ComplexF64.(1:length(sector.reps)))
        serial = similar(source)
        threaded = similar(source)
        Hv!(serial, source, sector, lat, hops, 1.0, 0.0, 0.0,
            diagonal; skip_tol=0.0)
        Hv_threaded!(threaded, source, sector, lat, hops, diagonal,
                     threaded_Hv_workspace(length(source)); skip_tol=0.0)
        @test threaded ≈ serial atol=1e-12
    end
end

@testset "threaded CSR response kernel" begin
    available = isdefined(Main, :threaded_csr)
    @test available
    if available
        matrix = sparse(ComplexF64[
            0 0 0 0;
            0 1 0 0;
            0 0 2 0;
            0 0 0 3
        ])
        csr = threaded_csr(matrix)
        source = ComplexF64[1, 2, 3, 4]
        product = similar(source)
        mul!(product, csr, source)
        @test product ≈ matrix * source
        @test size(csr) == size(matrix)
        @test nnz(csr) == nnz(matrix)

        ground = ComplexF64[1, 0, 0, 0]
        response_source = ComplexF64[0, 1, 1, 1]
        seen = Int[]
        kernel = response_lanczos(
            csr, 0.0, response_source; mmax=3, ground_state=ground,
            progress=(step, args...) -> push!(seen, step))
        @test !isempty(seen)
        for z in (0.4 + 0.065im, 2.0 + 0.065im)
            exact = dot(response_source, (z * I - matrix) \ response_source)
            @test lanczos_resolvent(kernel, z) ≈ exact atol=1e-11
        end
    end
end

@testset "representative-orbit structure factor" begin
    available = isdefined(Main, :structure_factor_survey_representatives)
    @test available
    if available
        lat = RectLat4x6()
        sector = build_ksector(gen_basis(lat.Ns, 2), lat, 2)
        state = normalize(ComplexF64.(1:length(sector.reps)))
        coefficients = expand_to_fock(sector, state)
        filling = 2 / lat.Ns
        reference = [compute_sq(
                         coefficients, lat, momentum[1], momentum[2], filling)
                     for momentum in lat.kpoints]
        ordered, _, norm2, particle_number =
            structure_factor_survey_representatives(
                sector, state, lat, 2)
        @test [row[2] for row in ordered] ≈ reference atol=1e-12
        @test norm2 ≈ 1.0 atol=1e-12
        @test particle_number ≈ 2.0 atol=1e-12
    end
end

@testset "30-site response driver is versioned" begin
    driver = joinpath(
        @__DIR__, "..", "cases", "case3_30sites_Np12",
        "run_optical_response.jl")
    @test isfile(driver)
    syntax_tree = Meta.parseall(read(driver, String))
    function has_parse_error(node)
        node isa Expr || return false
        node.head in (:error, :incomplete) && return true
        return any(has_parse_error, node.args)
    end
    @test !has_parse_error(syntax_tree)
end
