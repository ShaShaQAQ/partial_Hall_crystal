using Test
using ITensors
using ITensorMPS
using ITensorInfiniteMPS
using InfiniteCylinderDMRG
using Random

@testset "VUMPS runner helpers" begin
    cfg = InfiniteCylinderConfig(; Ly=2, x_period=3)
    sites, initstate, psi = initial_infinite_mps(cfg)
    @test length(sites) == sites_per_cell(cfg)
    @test flux(psi.AL) == QN()
    @test maximum(link_dimensions(psi)) == 1

    @test_throws ArgumentError initial_infinite_mps(cfg; occupied_sites=Int[])
    @test_throws ArgumentError initial_infinite_mps(cfg; occupied_sites=[1, 1])
    @test_throws ArgumentError initial_infinite_mps(cfg; occupied_sites=[0])
    @test_throws ArgumentError initial_infinite_mps(cfg; occupied_sites=[7])
    @test_throws ArgumentError initial_infinite_mps(cfg; occupied_sites=[true])

    @test centered_site_charges(cfg) == (empty=-1, occupied=5)
    @test isnothing(validate_centered_site_charges(cfg, sites))
    uncentered = siteinds("Fermion", sites_per_cell(cfg); conserve_qns=true)
    @test_throws ArgumentError validate_centered_site_charges(cfg, uncentered)
    dense_sites = siteinds("Fermion", sites_per_cell(cfg))
    @test_throws ArgumentError validate_centered_site_charges(cfg, dense_sites)
    @test_throws ArgumentError validate_centered_site_charges(cfg, sites[1:(end - 1)])

    @test vumps_converged(
        9e-6,
        8e-6,
        2;
        vumps_tol=1e-5,
        energy_tol=1e-4,
        energy_mismatch_tol=1e-4,
        stable_iterations=2,
    )
    @test !vumps_converged(
        2e-5,
        8e-6,
        2;
        vumps_tol=1e-5,
        energy_tol=1e-4,
        energy_mismatch_tol=1e-4,
        stable_iterations=2,
    )
    @test !vumps_converged(
        9e-6,
        2e-4,
        8e-6,
        2;
        vumps_tol=1e-5,
        energy_tol=1e-4,
        energy_mismatch_tol=1e-4,
        stable_iterations=2,
    )
    @test !vumps_converged(
        9e-6,
        8e-6,
        2e-4,
        2;
        vumps_tol=1e-5,
        energy_tol=1e-4,
        energy_mismatch_tol=1e-4,
        stable_iterations=2,
    )
    @test !vumps_converged(
        9e-6,
        8e-6,
        8e-6,
        1;
        vumps_tol=1e-5,
        energy_tol=1e-4,
        energy_mismatch_tol=1e-4,
        stable_iterations=2,
    )
    @test !vumps_converged(
        NaN,
        8e-6,
        8e-6,
        2;
        vumps_tol=1e-5,
        energy_tol=1e-4,
        energy_mismatch_tol=1e-4,
        stable_iterations=2,
    )
    @test !vumps_converged(
        -1e-6,
        8e-6,
        8e-6,
        2;
        vumps_tol=1e-5,
        energy_tol=1e-4,
        energy_mismatch_tol=1e-4,
        stable_iterations=2,
    )
    @test_throws ArgumentError vumps_converged(
        1e-6,
        1e-6,
        1e-6,
        2;
        vumps_tol=0.0,
        energy_tol=1e-4,
        energy_mismatch_tol=1e-4,
        stable_iterations=2,
    )

    for energy_normalization_sites in (6, 18)
        below_tol_per_site = 9e-5 * energy_normalization_sites
        above_tol_per_site = 1.1e-4 * energy_normalization_sites
        normalized = try
            vumps_converged(
                9e-6,
                below_tol_per_site,
                below_tol_per_site,
                2;
                vumps_tol=1e-5,
                energy_tol=1e-4,
                energy_mismatch_tol=1e-4,
                stable_iterations=2,
                energy_normalization_sites,
            )
        catch error
            error
        end
        @test normalized === true
        rejected = try
            vumps_converged(
                9e-6,
                above_tol_per_site,
                below_tol_per_site,
                2;
                vumps_tol=1e-5,
                energy_tol=1e-4,
                energy_mismatch_tol=1e-4,
                stable_iterations=2,
                energy_normalization_sites,
            )
        catch error
            error
        end
        @test rejected === false
    end
    for invalid_sites in (true, 0, -1)
        @test_throws ArgumentError vumps_converged(
            9e-6,
            8e-6,
            8e-6,
            2;
            vumps_tol=1e-5,
            energy_tol=1e-4,
            energy_mismatch_tol=1e-4,
            stable_iterations=2,
            energy_normalization_sites=invalid_sites,
        )
    end
end

@testset "VUMPS state canonicalization boundary" begin
    cfg = InfiniteCylinderConfig(; Ly=2, x_period=3)
    sites, _, psi = initial_infinite_mps(cfg)
    inconsistent = copy(psi)
    inconsistent.C[1] *= -1

    repaired = InfiniteCylinderDMRG._canonicalize_vumps_state(inconsistent)
    @test isnothing(InfiniteCylinderDMRG._validate_checkpoint_state(repaired, cfg))
    @test siteinds(only, repaired.AL) == sites
    @test link_dimensions(repaired) == link_dimensions(psi)
    @test flux(repaired.AL) == QN()
    @test all(
        isapprox(before.density, after.density; atol=1e-12, rtol=0) for
        (before, after) in zip(density_data(psi, cfg), density_data(repaired, cfg))
    )
end

@testset "paper D=2 VUMPS canonicalization regression" begin
    manifest = joinpath(@__DIR__, "..", "benchmarks", "fqahc_fig2.toml")
    spec = load_fig2_benchmark(manifest)
    dimension = 2
    point = 1
    candidate_id = first(
        InfiniteCylinderDMRG._default_fig2_candidate_ids(
            spec, dimension, point, nothing
        ),
    )

    saved_default_rng = copy(Random.default_rng())
    saved_index_rng = copy(ITensors.index_id_rng())
    try
        Random.seed!(Random.default_rng(), 0)
        Random.seed!(ITensors.index_id_rng(), 0)
        prepared = InfiniteCylinderDMRG._prepare_fig2_candidate_state(
            spec, spec.config, candidate_id, nothing
        )
        hamiltonian = build_infinite_mpo(
            spec.config, spec.model, prepared.sites
        )
        current = expand_subspace(
            prepared.psi, hamiltonian, dimension; cutoff=1.0e-8
        )
        @test isnothing(
            InfiniteCylinderDMRG._validate_checkpoint_state(
                current, spec.config
            ),
        )
        expansion_canonical =
            InfiniteCylinderDMRG._canonicalize_vumps_state(current)
        @test isnothing(
            InfiniteCylinderDMRG._validate_checkpoint_state(
                expansion_canonical, spec.config
            ),
        )
        @test siteinds(only, expansion_canonical.AL) ==
            siteinds(only, current.AL)
        @test link_dimensions(expansion_canonical) ==
            link_dimensions(current)
        @test flux(expansion_canonical.AL) == flux(current.AL)
        for _ in 1:3
            current = vumps_iteration(
                hamiltonian, current; vumps_tol=1.0e2, imaginary_tol=1.0e-12
            ).psi
        end

        expected_sites = siteinds(only, current.AL)
        expected_dimensions = link_dimensions(current)
        expected_flux = flux(current.AL)
        Random.seed!(Random.default_rng(), 0)
        Random.seed!(ITensors.index_id_rng(), 0)
        canonical = InfiniteCylinderDMRG._canonicalize_vumps_state(current)

        @test isnothing(
            InfiniteCylinderDMRG._validate_checkpoint_state(
                canonical, spec.config
            ),
        )
        @test siteinds(only, canonical.AL) == expected_sites
        @test link_dimensions(canonical) == expected_dimensions
        @test flux(canonical.AL) == expected_flux
        left_ids = Set(id.(linkinds(only, canonical.AL)))
        right_ids = Set(id.(linkinds(only, canonical.AR)))
        @test isempty(intersect(left_ids, right_ids))
        for site in 1:nsites(canonical)
            left_index = only(commoninds(canonical.C[site], canonical.AL[site]))
            right_index = only(commoninds(canonical.C[site], canonical.AR[site]))
            left_space = InfiniteCylinderDMRG._index_qn_dimensions(left_index)
            right_space = InfiniteCylinderDMRG._index_qn_dimensions(right_index)
            @test right_space == Dict(
                -charge => block_dimension for
                (charge, block_dimension) in left_space
            )
        end
    finally
        copy!(Random.default_rng(), saved_default_rng)
        copy!(ITensors.index_id_rng(), saved_index_rng)
    end
end

@testset "paper QN subspace expansion reaches requested target" begin
    manifest = joinpath(@__DIR__, "..", "benchmarks", "fqahc_fig2.toml")
    spec = load_fig2_benchmark(manifest)
    candidate_id = first(
        InfiniteCylinderDMRG._default_fig2_candidate_ids(
            spec, 4, 1, nothing
        ),
    )
    saved_default_rng = copy(Random.default_rng())
    saved_index_rng = copy(ITensors.index_id_rng())
    try
        Random.seed!(Random.default_rng(), 0)
        Random.seed!(ITensors.index_id_rng(), 0)
        prepared = InfiniteCylinderDMRG._prepare_fig2_candidate_state(
            spec, spec.config, candidate_id, nothing
        )
        hamiltonian = build_infinite_mpo(
            spec.config, spec.model, prepared.sites
        )
        expected_sites = siteinds(only, prepared.psi.AL)
        expected_flux = flux(prepared.psi.AL)
        result = run_vumps(
            hamiltonian,
            prepared.psi;
            maxdim_schedule=[4],
            cutoff=1.0e-9,
            max_iterations=2,
            vumps_tol=1.0e-6,
            energy_tol=1.0e-6,
            energy_mismatch_tol=1.0e-6,
            stable_iterations=2,
            canonical_seed=0,
        )
        expanded = result.psi

        @test maximum(link_dimensions(expanded)) == 4
        @test nsites(expanded) == sites_per_cell(spec.config)
        @test siteinds(only, expanded.AL) == expected_sites
        @test flux(expanded.AL) == expected_flux
    finally
        copy!(Random.default_rng(), saved_default_rng)
        copy!(ITensors.index_id_rng(), saved_index_rng)
    end
end

@testset "pinned VUMPS iteration boundary" begin
    initspin(_) = "↑"
    sites = infsiteinds("S=1/2", 1; initstate=initspin)
    psi = InfMPS(sites, initspin)
    H = InfiniteSum{MPO}(Model("ising"), sites; J=1.0, h=1.2)

    step = vumps_iteration(H, psi; vumps_tol=1e-5)
    @test step.psi isa InfiniteCanonicalMPS
    @test length(step.eps_left) == 1
    @test length(step.eps_right) == 1
    @test all(isfinite, step.eps_left)
    @test all(isfinite, step.eps_right)
    @test step.energy_left isa Float64
    @test step.energy_right isa Float64
    @test isfinite(step.elapsed_seconds)
    @test step.elapsed_seconds >= 0
    parallel_step = vumps_iteration(
        H,
        psi;
        vumps_tol=1e-5,
        multisite_update_alg=:parallel,
    )
    @test parallel_step.psi isa InfiniteCanonicalMPS
    @test all(isfinite, parallel_step.eps_left)
    @test all(isfinite, parallel_step.eps_right)
    @test_throws ArgumentError vumps_iteration(
        H,
        psi;
        vumps_tol=1e-5,
        multisite_update_alg=:unsupported,
    )
    @test_throws ArgumentError vumps_iteration(
        H, psi; vumps_tol=1e-5, imaginary_tol=-1e-12
    )
end

@testset "actual cylinder complex-energy boundary" begin
    cfg = InfiniteCylinderConfig(; Ly=2, x_period=3)
    params = CylinderModelParams(; t1=1.0, t3=0.2, V1=1.0, V2=0.0, V3=0.0)
    sites, _, psi = initial_infinite_mps(cfg)
    H = build_infinite_mpo(cfg, params, sites)

    result = run_vumps(
        H,
        psi;
        maxdim_schedule=[1],
        cutoff=1e-8,
        max_iterations=1,
        vumps_tol=1e-5,
        energy_tol=1e-4,
        energy_mismatch_tol=1e-4,
        stable_iterations=2,
        imaginary_tol=1e-12,
    )
    @test !result.converged
    @test length(result.records) == 1
    @test result.records[1].stage == 1
    @test result.records[1].iteration == 1
    @test isfinite(result.records[1].energy_left)
    @test isfinite(result.records[1].energy_right)
    @test result.records[1].energy_left ≈ 0.4 atol=1e-12
    @test result.records[1].energy_right ≈ 0.4 atol=1e-12
    @test result.records[1].energy_mismatch <= 1e-12
    @test isnothing(InfiniteCylinderDMRG._validate_checkpoint_state(result.psi, cfg))
    @test_throws ArgumentError run_vumps(
        H,
        psi;
        maxdim_schedule=[1],
        cutoff=1e-8,
        max_iterations=1,
        vumps_tol=1e-5,
        imaginary_tol=Inf,
    )
end

@testset "checked subspace expansion" begin
    @test InfiniteCylinderDMRG._subspace_expansion_progressed(
        [1, 1, 2, 2, 2, 1], [1, 2, 2, 2, 2, 2]
    )
    @test !InfiniteCylinderDMRG._subspace_expansion_progressed(
        [1, 1, 2, 2, 2, 1], [1, 1, 2, 2, 2, 1]
    )
    @test_throws DimensionMismatch InfiniteCylinderDMRG._subspace_expansion_progressed(
        [1, 2], [1, 2, 3]
    )

    initspin(_) = "↑"
    sites = infsiteinds("S=1/2", 1; initstate=initspin)
    psi = InfMPS(sites, initspin)
    H = InfiniteSum{MPO}(Model("ising"), sites; J=1.0, h=1.2)
    expanded = expand_subspace(psi, H, 2; cutoff=1e-8)
    @test maximum(link_dimensions(expanded)) > maximum(link_dimensions(psi))
    grown = run_vumps(
        H,
        psi;
        maxdim_schedule=[4],
        cutoff=1e-8,
        max_iterations=3,
        vumps_tol=1e-5,
        energy_tol=1e-4,
        energy_mismatch_tol=1e-4,
        stable_iterations=2,
    )
    @test maximum(link_dimensions(grown.psi)) == 4
    @test length(grown.expansions) >= 2
    @test all(record -> record.stage == 1 && record.target == 4, grown.expansions)
    @test !InfiniteCylinderDMRG._subspace_expansion_progressed(
        [2, 2], [1, 3]
    )

    stalled_H = InfiniteSum{MPO}(Model("ising"), sites; J=0.0, h=0.0)
    error = try
        expand_subspace(psi, stalled_H, 2; cutoff=1e-8)
        nothing
    catch caught
        caught
    end
    @test error isa ErrorException
    @test occursin("stalled at achieved maxdim=1", sprint(showerror, error))
    @test occursin("target maxdim=2", sprint(showerror, error))
    @test expand_subspace(expanded, H, 2; cutoff=1e-8) === expanded
    @test_throws ArgumentError expand_subspace(psi, H, 0; cutoff=1e-8)
    @test_throws ArgumentError expand_subspace(psi, H, 2; cutoff=0.0)
end

@testset "subspace expansion records" begin
    before = [1, 2, 3]
    after = [2, 3, 4]
    record = try
        SubspaceExpansionRecord(2, 8, before, after, true, 0.25)
    catch error
        error
    end
    @test record isa SubspaceExpansionRecord
    before[1] = 99
    after[end] = 99
    @test record.stage == 2
    @test record.target == 8
    @test record.before == [1, 2, 3]
    @test record.after == [2, 3, 4]
    @test record.progressed
    @test record.elapsed_seconds == 0.25
end

@testset "dense canonical states receive fresh right link IDs" begin
    initspin(_) = "↑"
    sites = infsiteinds("S=1/2", 1; initstate=initspin)
    psi = InfMPS(sites, initspin)
    canonical = ITensorMPS.orthogonalize(psi.AL, :)

    @test all(index -> !hasqns(index), inds(canonical.C[1]))
    normalized = InfiniteCylinderDMRG._normalize_right_link_convention(canonical)
    @test normalized !== canonical
    @test id.(linkinds(only, normalized.AL)) == id.(linkinds(only, canonical.AL))
    @test id.(linkinds(only, normalized.AR)) != id.(linkinds(only, canonical.AR))
    @test isempty(
        intersect(
            Set(id.(linkinds(only, normalized.AL))),
            Set(id.(linkinds(only, normalized.AR))),
        ),
    )
end

@testset "staged VUMPS runner" begin
    initspin(_) = "↑"
    sites = infsiteinds("S=1/2", 1; initstate=initspin)
    psi = InfMPS(sites, initspin)
    H = InfiniteSum{MPO}(Model("ising"), sites; J=1.0, h=1.2)
    options = (
        maxdim_schedule=[1],
        cutoff=1e-8,
        vumps_tol=1e-5,
        energy_tol=1e-4,
        energy_mismatch_tol=1e-4,
        stable_iterations=2,
    )

    stopped = run_vumps(H, psi; options..., max_iterations=1)
    @test !stopped.converged
    @test stopped.reason == "stage 1 reached maximum iterations (1) without convergence"
    @test length(stopped.records) == 1
    @test stopped.records[1].stage == 1
    @test stopped.records[1].iteration == 1
    @test !stopped.records[1].converged
    @test isempty(stopped.expansions)

    parallel_stopped = run_vumps(
        H,
        psi;
        options...,
        max_iterations=1,
        multisite_update_alg=:parallel,
    )
    @test !parallel_stopped.converged
    @test length(parallel_stopped.records) == 1
    @test_throws ArgumentError run_vumps(
        H,
        psi;
        options...,
        max_iterations=1,
        multisite_update_alg=:unsupported,
    )

    expanded_stop = run_vumps(
        H, psi; options..., maxdim_schedule=[2], max_iterations=1
    )
    @test !expanded_stop.converged
    @test expanded_stop.records[1].stage == 1
    @test expanded_stop.records[1].maxlinkdim == 2
    @test length(expanded_stop.expansions) == 1
    expansion = only(expanded_stop.expansions)
    @test expansion.stage == 1
    @test expansion.target == 2
    @test expansion.before == [1]
    @test expansion.after == [2]
    @test expansion.progressed
    @test isfinite(expansion.elapsed_seconds)
    @test expansion.elapsed_seconds >= 0

    @test maximum(link_dimensions(expanded_stop.psi)) == 2
    oversized_error = try
        run_vumps(
            H,
            expanded_stop.psi;
            options...,
            maxdim_schedule=[1],
            max_iterations=1,
        )
        nothing
    catch caught
        caught
    end
    @test oversized_error isa ArgumentError
    if oversized_error isa ArgumentError
        message = sprint(showerror, oversized_error)
        @test occursin("initial/achieved maxdim=2", message)
        @test occursin("first target=1", message)
    end

    converged = run_vumps(H, psi; options..., max_iterations=3)
    @test converged.converged
    @test converged.reason == "converged after 1 stage"
    @test length(converged.records) == 3
    @test converged.records[end].converged
    @test converged.records[end].precision_error < 1e-5
    @test converged.records[end].delta_energy < 1e-4
    @test converged.records[end].energy_mismatch < 1e-4

    @test_throws ArgumentError run_vumps(
        H, psi; options..., maxdim_schedule=Int[], max_iterations=1
    )
    @test_throws ArgumentError run_vumps(
        H, psi; options..., maxdim_schedule=[2, 1], max_iterations=1
    )
    @test_throws ArgumentError run_vumps(
        H, psi; options..., maxdim_schedule=[true], max_iterations=1
    )
    @test_throws ArgumentError run_vumps(H, psi; options..., max_iterations=0)
    @test_throws ArgumentError run_vumps(H, psi; options..., max_iterations=true)
    @test_throws ArgumentError run_vumps(
        H, psi; options..., stable_iterations=0, max_iterations=1
    )
    @test_throws ArgumentError run_vumps(
        H, psi; options..., stable_iterations=true, max_iterations=1
    )
end

@testset "VUMPS records and energy normalization" begin
    cfg = InfiniteCylinderConfig(; Ly=2, x_period=3)
    _, _, psi = initial_infinite_mps(cfg)
    record = VUMPSRecord(
        2, 3, 4, -1.0, -0.9, 0.1, 0.01, 1e-6, 2e-6, 2e-6, 0.25, false
    )
    @test record.stage == 2
    @test record.iteration == 3
    @test record.maxlinkdim == 4
    @test record.energy_mismatch == 0.1
    result = VUMPSResult(psi, [record], false, "maximum iterations reached")
    @test result.psi === psi
    @test result.records == [record]
    @test !result.converged
    @test result.reason == "maximum iterations reached"
    @test isempty(result.expansions)

    @test unit_cell_energy(2.5, 6) == 15.0
    @test unit_cell_energy([1.0, 2.0, 3.0], 3) == 6.0
    @test unit_cell_energy(2.5 + 1e-13im, 6; imaginary_tol=1e-12) == 15.0
    @test unit_cell_energy(
        ComplexF64[1.0 + 1e-13im, 2.0 - 2e-13im, 3.0],
        3;
        imaginary_tol=1e-12,
    ) == 6.0
    @test_throws ArgumentError unit_cell_energy([1.0, 2.0], 3)
    @test_throws ArgumentError unit_cell_energy(Float64[], 3)
    @test_throws ArgumentError unit_cell_energy(Inf, 3)
    @test_throws ArgumentError unit_cell_energy(floatmax(Float64), 2)
    @test_throws ArgumentError unit_cell_energy(
        [floatmax(Float64), floatmax(Float64)], 2
    )
    @test_throws ArgumentError unit_cell_energy(
        BigFloat(floatmax(Float64)) * 2, 1
    )
    @test_throws ArgumentError unit_cell_energy(
        2.5 + 2e-12im, 6; imaginary_tol=1e-12
    )
    @test_throws ArgumentError unit_cell_energy(
        ComplexF64[1.0, 2.0 + 2e-12im, 3.0], 3; imaginary_tol=1e-12
    )
    for invalid_energy in (ComplexF64(Inf, 0.0), ComplexF64(0.0, Inf), ComplexF64(NaN, 0.0))
        @test_throws ArgumentError unit_cell_energy(invalid_energy, 3)
        @test_throws ArgumentError unit_cell_energy([invalid_energy], 1)
    end
    oversized_tol = BigFloat(floatmax(Float64)) * 2
    for invalid_tol in (-1e-12, NaN, Inf, oversized_tol)
        @test_throws ArgumentError unit_cell_energy(
            1.0 + 0.0im, 1; imaginary_tol=invalid_tol
        )
        @test_throws ArgumentError unit_cell_energy(
            ComplexF64[1.0], 1; imaginary_tol=invalid_tol
        )
    end
end
