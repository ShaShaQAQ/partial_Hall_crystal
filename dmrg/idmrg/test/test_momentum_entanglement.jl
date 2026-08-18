using Test
using LinearAlgebra
using Random
using ITensors
using InfiniteCylinderDMRG

@testset "deterministic momentum RNG" begin
    @test isdefined(InfiniteCylinderDMRG, :DEFAULT_MOMENTUM_RNG_SEED)
    @test isdefined(InfiniteCylinderDMRG, :_resolve_momentum_rng)

    if isdefined(InfiniteCylinderDMRG, :DEFAULT_MOMENTUM_RNG_SEED) &&
            isdefined(InfiniteCylinderDMRG, :_resolve_momentum_rng)
        default_seed = InfiniteCylinderDMRG.DEFAULT_MOMENTUM_RNG_SEED
        resolve_rng = InfiniteCylinderDMRG._resolve_momentum_rng
        rng_bytes(rng) = collect(reinterpret(UInt8, rand(rng, UInt64, 16)))

        @test default_seed isa UInt64
        @test default_seed == 0x4d4f4d454e54554d

        first_rng = resolve_rng(nothing, default_seed)
        second_rng = resolve_rng(nothing, default_seed)
        @test first_rng isa Random.Xoshiro
        @test second_rng isa Random.Xoshiro
        @test first_rng !== second_rng
        @test rng_bytes(first_rng) == rng_bytes(second_rng)

        parent_task = current_task()
        before_global_change = fetch(@async begin
            Random.seed!(0x1111)
            rand(UInt64, 37)
            (task=current_task(), bytes=rng_bytes(resolve_rng(nothing, nothing)))
        end)
        after_global_change = fetch(@async begin
            Random.seed!(0x2222)
            rand(UInt64, 91)
            (task=current_task(), bytes=rng_bytes(resolve_rng(nothing, nothing)))
        end)
        @test before_global_change.bytes == after_global_change.bytes
        @test before_global_change.task !== parent_task
        @test after_global_change.task !== parent_task
        @test before_global_change.task !== after_global_change.task

        @test rng_bytes(resolve_rng(nothing, UInt64(1))) !=
            rng_bytes(resolve_rng(nothing, UInt64(2)))
        @test resolve_rng(nothing, UInt64(0)) isa Random.Xoshiro
        @test resolve_rng(nothing, typemax(UInt64)) isa Random.Xoshiro

        explicit_rng = Random.Xoshiro(0x7a57)
        @test resolve_rng(explicit_rng, nothing) === explicit_rng
        @test_throws ArgumentError resolve_rng(explicit_rng, default_seed)
        @test_throws ArgumentError resolve_rng(nothing, true)
        @test_throws ArgumentError resolve_rng(nothing, -1)
        @test_throws ArgumentError resolve_rng(
            nothing,
            UInt128(typemax(UInt64)) + one(UInt128),
        )
    end

    momentum_binding = Base.Docs.Binding(
        InfiniteCylinderDMRG,
        :momentum_entanglement_data,
    )
    momentum_doc = Base.Docs._doc(momentum_binding)
    @test momentum_doc isa Base.Docs.DocStr
    if momentum_doc isa Base.Docs.DocStr
        doc_text = lowercase(join(string.(momentum_doc.text), "\n"))
        @test occursin("0x4d4f4d454e54554d", doc_text)
        @test occursin("consum", doc_text)
        @test occursin("mutually exclusive", doc_text)
        @test occursin("uint64", doc_text)
    end

    config = InfiniteCylinderConfig(;
        geometry=:paper_straight,
        Ny=2,
        x_period=1,
        filling_num=1,
        filling_den=1,
    )
    _, _, psi = initial_infinite_mps(config; occupied_sites=[1, 3])
    seed = UInt64(0x7a57)
    seeded_result = try
        momentum_entanglement_data(
            psi,
            config;
            cut_x=1,
            residual_tol=1e-8,
            transfer_tol=1e-8,
            rng_seed=seed,
        )
    catch error
        error
    end
    explicit_result = try
        momentum_entanglement_data(
            psi,
            config;
            cut_x=1,
            residual_tol=1e-8,
            transfer_tol=1e-8,
            rng=Random.Xoshiro(seed),
        )
    catch error
        error
    end
    @test seeded_result isa MomentumEntanglementData
    @test explicit_result isa MomentumEntanglementData
    if seeded_result isa MomentumEntanglementData &&
            explicit_result isa MomentumEntanglementData
        level_snapshot(data) = [
            (
                level.probability,
                level.entanglement_energy,
                level.raw_charge,
                level.ky_index,
                level.ky,
                level.eigenphase,
                level.phase_residual,
            ) for level in data.levels
        ]
        validation_snapshot(data) = (
            data.validation.unitary_residual,
            data.validation.power_residual,
            data.validation.commutator_residual,
            data.validation.transfer_residual,
            data.validation.valid,
            data.validation.reason,
        )
        @test level_snapshot(seeded_result) == level_snapshot(explicit_result)
        @test validation_snapshot(seeded_result) ==
            validation_snapshot(explicit_result)
        @test seeded_result.momentum_convention ===
            explicit_result.momentum_convention
        @test seeded_result.sector_ky_offsets ==
            explicit_result.sector_ky_offsets
    end

    conflicting_keywords = try
        momentum_entanglement_data(
            psi,
            config;
            cut_x=1,
            rng=Random.Xoshiro(seed),
            rng_seed=seed,
        )
        nothing
    catch error
        error
    end
    @test conflicting_keywords isa ArgumentError
end

@testset "sector-relative momentum convention" begin

@testset "momentum label validation" begin
    rho = Diagonal([0.4, 0.3, 0.2, 0.1])
    Ty = Diagonal(ComplexF64[1, cis(pi), 1, cis(pi)])
    result = validate_momentum_blocks(rho, Ty; Ny=2, residual_tol=1e-12)

    @test result.valid
    @test result.validation.valid
    @test hasproperty(result, :momentum_convention)
    @test hasproperty(result, :sector_ky_offsets)
    if hasproperty(result, :momentum_convention) &&
            hasproperty(result, :sector_ky_offsets)
        @test result.momentum_convention ===
            :sector_relative_canonical_cyclic_orbit
        @test !isempty(result.sector_ky_offsets)
    end
    @test occursin("sector-relative Δky", result.reason)
    @test all(level -> level.ky_index in 0:1, result.levels)
    @test maximum(level.phase_residual for level in result.levels) < 1e-12
    @test result.validation.unitary_residual < 1e-12
    @test result.validation.power_residual < 1e-12
    @test result.validation.commutator_residual < 1e-12

    bad_unitary = validate_momentum_blocks(rho, 1.01Ty; Ny=2, residual_tol=1e-12)
    @test !bad_unitary.valid
    @test all(level -> ismissing(level.ky_index), bad_unitary.levels)
    if hasproperty(bad_unitary, :momentum_convention) &&
            hasproperty(bad_unitary, :sector_ky_offsets)
        @test bad_unitary.momentum_convention === :invalid
        @test isempty(bad_unitary.sector_ky_offsets)
    end

    broken = ComplexF64[1 0.1; 0 1]
    broken_result = validate_momentum_blocks(
        Diagonal([0.6, 0.4]),
        broken;
        Ny=2,
        residual_tol=1e-12,
    )
    @test !broken_result.valid
    @test all(level -> ismissing(level.ky_index), broken_result.levels)
end

@testset "sector-relative momentum is invariant under charge-sector gauge" begin
    Ny = 3
    probabilities = [0.35, 0.15, 0.30, 0.20]
    charges = [0, 0, 1, 1]
    raw_momenta = [1, 2, 2, 0]
    translation = Diagonal(ComplexF64[cis(2pi * momentum / Ny) for momentum in raw_momenta])
    gauge_phase = cis(2pi / Ny)
    gauge_rotated = Matrix(translation)
    gauge_rotated[3:4, 3:4] .*= gauge_phase

    reference = validate_momentum_blocks(
        Diagonal(probabilities),
        translation;
        Ny,
        raw_charges=charges,
        residual_tol=1e-12,
    )
    rotated = validate_momentum_blocks(
        Diagonal(probabilities),
        gauge_rotated;
        Ny,
        raw_charges=charges,
        residual_tol=1e-12,
    )

    @test reference.valid
    @test rotated.valid
    @test [(level.raw_charge, level.ky_index) for level in reference.levels] ==
        [(level.raw_charge, level.ky_index) for level in rotated.levels]
    @test [level.ky_index for level in reference.levels] == [0, 1, 0, 1]
    if hasproperty(reference, :sector_ky_offsets) &&
            hasproperty(rotated, :sector_ky_offsets)
        @test reference.sector_ky_offsets[0] == 2
        @test reference.sector_ky_offsets[1] == 1
        @test rotated.sector_ky_offsets[0] == 2
        @test rotated.sector_ky_offsets[1] == 0
    end
    @test !isapprox(
        [level.eigenphase for level in reference.levels if level.raw_charge == 1],
        [level.eigenphase for level in rotated.levels if level.raw_charge == 1];
        atol=1e-12,
        rtol=0,
    )
end

@testset "canonical cyclic orbit handles a degenerate wrap" begin
    Ny = 3
    probabilities = [0.4, 0.4, 0.2]
    raw_momenta = [0, 1, 2]
    shifted_momenta = mod.(raw_momenta .+ 2, Ny)
    translation = Diagonal(ComplexF64[cis(2pi * momentum / Ny) for momentum in raw_momenta])
    shifted_translation = Diagonal(
        ComplexF64[cis(2pi * momentum / Ny) for momentum in shifted_momenta]
    )

    reference = validate_momentum_blocks(
        Diagonal(probabilities),
        translation;
        Ny,
        raw_charges=zeros(Int, 3),
        residual_tol=1e-12,
    )
    shifted = validate_momentum_blocks(
        Diagonal(probabilities),
        shifted_translation;
        Ny,
        raw_charges=zeros(Int, 3),
        residual_tol=1e-12,
    )

    @test reference.valid
    @test shifted.valid
    @test reference.levels[1].probability == reference.levels[2].probability
    @test [level.ky_index for level in reference.levels] == [0, 1, 2]
    @test [level.ky_index for level in shifted.levels] == [0, 1, 2]
    @test [
        (level.probability, level.raw_charge, level.ky_index) for
        level in reference.levels
    ] == [
        (level.probability, level.raw_charge, level.ky_index) for
        level in shifted.levels
    ]
    if hasproperty(reference, :sector_ky_offsets) &&
            hasproperty(shifted, :sector_ky_offsets)
        @test reference.sector_ky_offsets[0] == 0
        @test shifted.sector_ky_offsets[0] == 1
    end
end

@testset "zero-weight momentum sectors are excluded without renormalizing" begin
    Ny = 3
    zero_sector = validate_momentum_blocks(
        Diagonal([0.55, 0.45, 0.0, 0.0]),
        Diagonal(ComplexF64[cis(2pi / 3), cis(4pi / 3), 1, cis(2pi / 3)]);
        Ny,
        raw_charges=[0, 0, 1, 1],
        residual_tol=1e-12,
    )

    @test zero_sector.valid
    @test length(zero_sector.levels) == 2
    @test all(level -> level.raw_charge == 0, zero_sector.levels)
    @test all(level -> !ismissing(level.ky_index), zero_sector.levels)
    if hasproperty(zero_sector, :sector_ky_offsets)
        @test sort!(collect(keys(zero_sector.sector_ky_offsets))) == [0]
    end

    thresholded = try
        validate_momentum_blocks(
            Diagonal([0.50, 0.30, 0.15, 0.05]),
            Diagonal(ComplexF64[1, cis(2pi / 3), cis(4pi / 3), 1]);
            Ny,
            raw_charges=[0, 0, 1, 1],
            residual_tol=1e-12,
            probability_tol=0.10,
        )
    catch error
        error
    end
    @test thresholded isa MomentumEntanglementData
    if thresholded isa MomentumEntanglementData
        @test thresholded.valid
        @test [level.probability for level in thresholded.levels] == [0.50, 0.30, 0.15]
        @test sum(level.probability for level in thresholded.levels) ≈ 0.95 atol=1e-15
        @test sort!(collect(keys(thresholded.sector_ky_offsets))) == [0, 1]
        @test all(level -> !ismissing(level.ky_index), thresholded.levels)
    end

    empty_support = try
        validate_momentum_blocks(
            Diagonal([0.6, 0.4]),
            Diagonal(ComplexF64[1, cis(2pi / 3)]);
            Ny,
            probability_tol=1.0,
        )
    catch error
        error
    end
    @test empty_support isa ArgumentError
    if empty_support isa ArgumentError
        @test occursin(
            "probability_tol removed all positive Schmidt support",
            sprint(showerror, empty_support),
        )
    end
end

end

@testset "degenerate charge blocks and strict invalidation" begin
    rho = Diagonal([0.4, 0.4, 0.1, 0.1])
    Ty = ComplexF64[
        0 1 0 0
        1 0 0 0
        0 0 1 0
        0 0 0 -1
    ]
    result = validate_momentum_blocks(
        rho,
        Ty;
        Ny=2,
        raw_charges=[0, 0, 1, 1],
        residual_tol=1e-12,
    )
    @test result.valid
    @test count(level -> level.ky_index == 0, result.levels) == 2
    @test count(level -> level.ky_index == 1, result.levels) == 2
    @test [level.raw_charge for level in result.levels] == [0, 0, 1, 1]

    charge_mixing = validate_momentum_blocks(
        Diagonal([0.5, 0.5]),
        ComplexF64[0 1; 1 0];
        Ny=2,
        raw_charges=[0, 1],
        residual_tol=1e-12,
    )
    @test !charge_mixing.valid
    @test occursin("charge", lowercase(charge_mixing.reason))
    @test all(level -> ismissing(level.ky_index), charge_mixing.levels)

    cross_charge_density = validate_momentum_blocks(
        ComplexF64[0.5 0.1; 0.1 0.5],
        Matrix{ComplexF64}(I, 2, 2);
        Ny=2,
        raw_charges=[0, 1],
        residual_tol=1e-12,
    )
    @test !cross_charge_density.valid
    @test occursin("charge-density", lowercase(cross_charge_density.reason))
    @test all(level -> ismissing(level.ky_index), cross_charge_density.levels)

    transfer_failure = validate_momentum_blocks(
        Diagonal([0.6, 0.4]),
        Diagonal(ComplexF64[1, -1]);
        Ny=2,
        transfer_residual=1e-4,
        residual_tol=1e-12,
    )
    @test !transfer_failure.valid
    @test occursin("transfer", lowercase(transfer_failure.reason))
    @test all(level -> ismissing(level.ky_index), transfer_failure.levels)
end

@testset "momentum validation input checks" begin
    rho = Diagonal([0.6, 0.4])
    Ty = Diagonal(ComplexF64[1, -1])
    @test_throws ArgumentError validate_momentum_blocks(rho, Ty; Ny=0)
    @test_throws ArgumentError validate_momentum_blocks(rho, Ty; Ny=true)
    @test_throws ArgumentError validate_momentum_blocks(rho, Ty; Ny=2, residual_tol=0)
    @test_throws DimensionMismatch validate_momentum_blocks(rho, Matrix{ComplexF64}(I, 3, 3); Ny=2)
    @test_throws DimensionMismatch validate_momentum_blocks(
        rho,
        Ty;
        Ny=2,
        raw_charges=[0],
    )
end

@testset "Task 7D3 blockwise production contract" begin

@testset "required blockwise helpers" begin
    @test isdefined(InfiniteCylinderDMRG, :_validate_momentum_block_vectors)
    @test isdefined(InfiniteCylinderDMRG, :_schmidt_density_blocks)
    @test isdefined(InfiniteCylinderDMRG, :_twisted_block_matrices)
end

@testset "blockwise momentum validation matches the dense compatibility API" begin
    @test isdefined(InfiniteCylinderDMRG, :_validate_momentum_block_vectors)
    if isdefined(InfiniteCylinderDMRG, :_validate_momentum_block_vectors)
        validate_blocks =
            InfiniteCylinderDMRG._validate_momentum_block_vectors
        density_blocks = Matrix{ComplexF64}[
            ComplexF64[0.28 0; 0 0.22],
            reshape(ComplexF64[0.20], 1, 1),
            ComplexF64[0.18 0; 0 0.12],
        ]
        translation_blocks = Matrix{ComplexF64}[
            ComplexF64[1 0; 0 -1],
            reshape(ComplexF64[1], 1, 1),
            ComplexF64[-1 0; 0 1],
        ]
        block_charges = [-2, 0, 3]
        dense_density = cat(density_blocks...; dims=(1, 2))
        dense_translation = cat(translation_blocks...; dims=(1, 2))
        level_charges = reduce(
            vcat,
            [fill(charge, size(block, 1)) for
             (charge, block) in zip(block_charges, density_blocks)],
        )

        dense = validate_momentum_blocks(
            dense_density,
            dense_translation;
            Ny=2,
            raw_charges=level_charges,
            residual_tol=1e-12,
        )
        blocked = validate_blocks(
            density_blocks,
            translation_blocks;
            Ny=2,
            raw_charges=block_charges,
            residual_tol=1e-12,
        )

        @test blocked.valid == dense.valid
        @test blocked.reason == dense.reason
        @test blocked.momentum_convention == dense.momentum_convention
        @test blocked.sector_ky_offsets == dense.sector_ky_offsets
        @test blocked.validation.valid == dense.validation.valid
        @test blocked.validation.reason == dense.validation.reason
        @test blocked.validation.unitary_residual ≈
            dense.validation.unitary_residual atol=1e-15 rtol=0
        @test blocked.validation.power_residual ≈
            dense.validation.power_residual atol=1e-15 rtol=0
        @test blocked.validation.commutator_residual ≈
            dense.validation.commutator_residual atol=1e-15 rtol=0
        @test blocked.validation.transfer_residual ==
            dense.validation.transfer_residual
        @test length(blocked.levels) == length(dense.levels)
        for (blocked_level, dense_level) in zip(blocked.levels, dense.levels)
            @test blocked_level.level == dense_level.level
            @test blocked_level.probability ≈ dense_level.probability atol=1e-15 rtol=0
            @test blocked_level.entanglement_energy ≈
                dense_level.entanglement_energy atol=1e-15 rtol=0
            @test blocked_level.raw_charge == dense_level.raw_charge
            @test blocked_level.ky_index == dense_level.ky_index
            @test blocked_level.ky == dense_level.ky
            @test blocked_level.eigenphase ≈ dense_level.eigenphase atol=1e-15 rtol=0
            @test blocked_level.phase_residual ≈
                dense_level.phase_residual atol=1e-15 rtol=0
        end

        bad_translation_blocks = copy(translation_blocks)
        bad_translation_blocks[2] = 1.01 .* bad_translation_blocks[2]
        bad_dense_translation = cat(bad_translation_blocks...; dims=(1, 2))
        dense_bad = validate_momentum_blocks(
            dense_density,
            bad_dense_translation;
            Ny=2,
            raw_charges=level_charges,
            residual_tol=1e-12,
        )
        blocked_bad = validate_blocks(
            density_blocks,
            bad_translation_blocks;
            Ny=2,
            raw_charges=block_charges,
            residual_tol=1e-12,
        )
        @test !blocked_bad.valid
        @test blocked_bad.reason == dense_bad.reason
        @test blocked_bad.validation.unitary_residual ≈
            dense_bad.validation.unitary_residual atol=1e-15 rtol=0
        @test blocked_bad.validation.power_residual ≈
            dense_bad.validation.power_residual atol=1e-15 rtol=0
        @test all(level -> ismissing(level.ky_index), blocked_bad.levels)
    end
end

@testset "blockwise validation scales with stored block elements" begin
    @test isdefined(InfiniteCylinderDMRG, :_validate_momentum_block_vectors)
    if isdefined(InfiniteCylinderDMRG, :_validate_momentum_block_vectors)
        block_count = 600
        block_dimension = 4
        total_dimension = block_count * block_dimension
        probability = inv(Float64(total_dimension))
        density_blocks = [
            Matrix{ComplexF64}(I, block_dimension, block_dimension) .* probability for
            _ in 1:block_count
        ]
        translation_blocks = [
            Matrix{ComplexF64}(I, block_dimension, block_dimension) for
            _ in 1:block_count
        ]
        stored_density_elements = sum(length, density_blocks)

        result = InfiniteCylinderDMRG._validate_momentum_block_vectors(
            density_blocks,
            translation_blocks;
            Ny=2,
            raw_charges=collect(0:(block_count - 1)),
            residual_tol=1e-12,
        )

        @test result.valid
        @test length(result.levels) == total_dimension
        @test stored_density_elements == block_count * block_dimension^2
        @test stored_density_elements < total_dimension^2 ÷ 100
        @test sum(level.probability for level in result.levels) ≈ 1.0 atol=1e-12
    end
end

@testset "QN Schmidt density is extracted blockwise" begin
    @test isdefined(InfiniteCylinderDMRG, :_schmidt_density_blocks)
    if isdefined(InfiniteCylinderDMRG, :_schmidt_density_blocks)
        left_index = Index(
            [QN("Nf", 0) => 2, QN("Nf", 1) => 1],
            "Link,Left",
        )
        right_index = dag(Index(
            [QN("Nf", 0) => 2, QN("Nf", 1) => 1],
            "Link,Right",
        ))
        center = ITensor(ComplexF64, QN(), left_index, right_index)
        center[
            left_index => 1:2,
            right_index => 1:2,
        ] = ComplexF64[1 2im; 3 4]
        center[left_index => 3, right_index => 3] = 5im
        dense_center = Matrix{ComplexF64}(
            Array(center, left_index, right_index)
        )
        dense_reference = dense_center * dense_center'

        density_blocks, block_charges =
            InfiniteCylinderDMRG._schmidt_density_blocks(
                center, left_index, right_index
            )
        assembled_for_test = cat(density_blocks...; dims=(1, 2))

        @test density_blocks isa Vector{Matrix{ComplexF64}}
        @test block_charges == [0, 1]
        @test size.(density_blocks) == [(2, 2), (1, 1)]
        @test assembled_for_test ≈ dense_reference atol=1e-14 rtol=0
        @test sum(length, density_blocks) == 5
        @test all(block -> all(isfinite, block), density_blocks)
    end
end

end

@testset "positive transverse translation phase convention" begin
    Ny = 3
    phase = cis(2pi / Ny)
    state_vector = zeros(ComplexF64, 1 << (2Ny))
    for ring in 0:(Ny - 1)
        one_particle_bits = 1 << (2ring)
        state_vector[one_particle_bits + 1] = cis(-2pi * ring / Ny) / sqrt(Ny)
    end
    translated_vector = apply_fermionic_ring_translation(state_vector, Ny)

    @test translated_vector ≈ phase * state_vector atol=1e-14 rtol=0
    @test dot(state_vector, translated_vector) ≈ phase atol=1e-14 rtol=0
    @test dot(translated_vector, state_vector) ≈ conj(phase) atol=1e-14 rtol=0
    @test !isapprox(phase, conj(phase); atol=1e-14, rtol=0)
end

@testset "finite-cell FSWAP helpers mutate only their input" begin
    @test isdefined(InfiniteCylinderDMRG, :_validate_finite_cell_topology)
    @test isdefined(InfiniteCylinderDMRG, :_validate_fermionic_swap)
    @test isdefined(InfiniteCylinderDMRG, :_validate_fermionic_move)
    @test isdefined(InfiniteCylinderDMRG, :_swap_fermionic_mps_sites!)
    @test isdefined(InfiniteCylinderDMRG, :_move_fermionic_mps_site!)
    @test isdefined(InfiniteCylinderDMRG, :_translate_left_cell!)

    if isdefined(InfiniteCylinderDMRG, :_swap_fermionic_mps_sites!) &&
            isdefined(InfiniteCylinderDMRG, :_move_fermionic_mps_site!)
        config = InfiniteCylinderConfig(;
            geometry=:paper_straight,
            Ny=3,
            x_period=1,
            filling_num=2,
            filling_den=3,
        )
        _, _, psi = initial_infinite_mps(config; occupied_sites=[1, 5])
        bond = config.Ly
        configuration_tag = InfiniteCylinderDMRG._configuration_signature(config)
        finite_cell = InfiniteCylinderDMRG.MPS(
            [psi.AL[bond + offset] for offset in 1:InfiniteCylinderDMRG.nsites(psi)]
        )

        function mps_snapshot(state)
            return (
                tensors=[copy(tensor) for tensor in state],
                leftlim=InfiniteCylinderDMRG.ITensorMPS.leftlim(state),
                rightlim=InfiniteCylinderDMRG.ITensorMPS.rightlim(state),
            )
        end

        function test_snapshot_unchanged(state, snapshot)
            @test InfiniteCylinderDMRG.ITensorMPS.leftlim(state) ==
                snapshot.leftlim
            @test InfiniteCylinderDMRG.ITensorMPS.rightlim(state) ==
                snapshot.rightlim
            @test all(
                inds(actual) == inds(expected) &&
                    isapprox(actual, expected; atol=0, rtol=0) for
                (actual, expected) in zip(state, snapshot.tensors)
            )
        end

        function test_rejected_without_mutation(operation)
            state = copy(finite_cell)
            snapshot = mps_snapshot(state)
            error = try
                operation(state)
                nothing
            catch caught
                caught
            end
            @test error isa ArgumentError
            test_snapshot_unchanged(state, snapshot)
        end

        function invalidate_physical_tag!(state, site, tag)
            physical = InfiniteCylinderDMRG._physical_site_index(
                state, site, tag
            )
            invalid = removetags(physical, tag)
            state[site] = replaceinds(state[site], physical => invalid)
            return state
        end

        function test_translation_rejected_without_mutation(target, original)
            target_snapshot = mps_snapshot(target)
            original_snapshot = mps_snapshot(original)
            error = try
                InfiniteCylinderDMRG._translate_left_cell!(
                    target,
                    original,
                    config;
                    configuration_tag,
                    cutoff=0.0,
                )
                nothing
            catch caught
                caught
            end
            @test error isa ArgumentError
            test_snapshot_unchanged(target, target_snapshot)
            test_snapshot_unchanged(original, original_snapshot)
        end

        unchanged = copy(finite_cell)
        unchanged_snapshot = mps_snapshot(unchanged)
        move_result = InfiniteCylinderDMRG._move_fermionic_mps_site!(
            unchanged,
            2,
            2;
            configuration_tag,
            cutoff=0.0,
        )
        @test move_result === unchanged
        test_snapshot_unchanged(unchanged, unchanged_snapshot)

        swapped = copy(finite_cell)
        swap_result = InfiniteCylinderDMRG._swap_fermionic_mps_sites!(
            swapped,
            1;
            configuration_tag,
            ortho="left",
            cutoff=0.0,
        )
        @test swap_result === swapped

        moved = copy(finite_cell)
        moved_result = InfiniteCylinderDMRG._move_fermionic_mps_site!(
            moved,
            config.Ly - 1,
            1;
            configuration_tag,
            cutoff=0.0,
        )
        @test moved_result === moved

        for invalid_bond in (0, length(finite_cell))
            test_rejected_without_mutation() do state
                InfiniteCylinderDMRG._swap_fermionic_mps_sites!(
                    state,
                    invalid_bond;
                    configuration_tag,
                    ortho="left",
                    cutoff=0.0,
                )
            end
        end
        for invalid_ortho in (:left, "center")
            test_rejected_without_mutation() do state
                InfiniteCylinderDMRG._swap_fermionic_mps_sites!(
                    state,
                    1;
                    configuration_tag,
                    ortho=invalid_ortho,
                    cutoff=0.0,
                )
            end
        end
        for invalid_cutoff in (true, -1.0, Inf, NaN, "0")
            test_rejected_without_mutation() do state
                InfiniteCylinderDMRG._swap_fermionic_mps_sites!(
                    state,
                    1;
                    configuration_tag,
                    ortho="left",
                    cutoff=invalid_cutoff,
                )
            end
        end
        test_rejected_without_mutation() do state
            InfiniteCylinderDMRG._swap_fermionic_mps_sites!(
                state,
                1;
                configuration_tag="not-the-configuration-tag",
                ortho="left",
                cutoff=0.0,
            )
        end

        count = length(finite_cell)
        for (invalid_source, invalid_target) in (
            (0, 0),
            (count + 1, count + 1),
            (0, 1),
            (1, 0),
            (count + 1, 1),
            (1, count + 1),
        )
            test_rejected_without_mutation() do state
                InfiniteCylinderDMRG._move_fermionic_mps_site!(
                    state,
                    invalid_source,
                    invalid_target;
                    configuration_tag,
                    cutoff=0.0,
                )
            end
        end


        late_invalid = copy(finite_cell)
        invalidate_physical_tag!(late_invalid, 5, configuration_tag)
        late_snapshot = mps_snapshot(late_invalid)
        late_error = try
            InfiniteCylinderDMRG._move_fermionic_mps_site!(
                late_invalid,
                1,
                5;
                configuration_tag,
                cutoff=0.0,
            )
            nothing
        catch caught
            caught
        end
        @test late_error isa ArgumentError
        test_snapshot_unchanged(late_invalid, late_snapshot)

        truncated_original = InfiniteCylinderDMRG.MPS(
            [finite_cell[site] for site in 1:(length(finite_cell) - 1)]
        )
        test_translation_rejected_without_mutation(
            copy(truncated_original), truncated_original
        )

        aliased_cell = copy(finite_cell)
        test_translation_rejected_without_mutation(aliased_cell, aliased_cell)

        late_invalid_original = copy(finite_cell)
        invalidate_physical_tag!(
            late_invalid_original,
            length(late_invalid_original),
            configuration_tag,
        )
        test_translation_rejected_without_mutation(
            copy(finite_cell), late_invalid_original
        )

        missing_target_boundary = copy(finite_cell)
        target_right_boundary = only(filter(
            index -> !hastags(index, configuration_tag),
            uniqueinds(
                missing_target_boundary[end],
                missing_target_boundary[end - 1],
            ),
        ))
        missing_target_boundary[end] *= onehot(dag(target_right_boundary) => 1)
        test_translation_rejected_without_mutation(
            missing_target_boundary, copy(finite_cell)
        )

        extra_original_boundary = copy(finite_cell)
        extra_original_boundary[end] *= onehot(
            Index(QN() => 1; tags="invalid-extra-right-boundary") => 1
        )
        test_translation_rejected_without_mutation(
            copy(finite_cell), extra_original_boundary
        )

        extra_interior = copy(finite_cell)
        extra_interior[3] *= onehot(
            Index(QN() => 1; tags="invalid-extra-interior-index") => 1
        )
        test_translation_rejected_without_mutation(
            extra_interior, copy(finite_cell)
        )

        partial_route_invalid = copy(finite_cell)
        partial_route_invalid[5] *= onehot(
            Index(QN() => 1; tags="invalid-off-route-interior-index") => 1
        )
        partial_route_snapshot = mps_snapshot(partial_route_invalid)
        partial_route_error = try
            InfiniteCylinderDMRG._move_fermionic_mps_site!(
                partial_route_invalid,
                1,
                2;
                configuration_tag,
                cutoff=0.0,
            )
            nothing
        catch caught
            caught
        end
        @test partial_route_error isa ArgumentError
        test_snapshot_unchanged(
            partial_route_invalid, partial_route_snapshot
        )
    end
end

@testset "finite-cell FSWAP translation matches exact Fock translation" begin
    function finite_cell_coefficients(cell, configuration_tag)
        physical_indices = [
            InfiniteCylinderDMRG._physical_site_index(
                cell, site, configuration_tag
            ) for site in eachindex(cell)
        ]
        wavefunction = foldl(*, cell)
        boundary_indices = filter(
            index -> !hastags(index, configuration_tag),
            collect(inds(wavefunction)),
        )
        return vec(Array(wavefunction, physical_indices..., boundary_indices...))
    end

    function finite_cell_snapshot(cell)
        return (
            tensors=[copy(tensor) for tensor in cell],
            leftlim=InfiniteCylinderDMRG.ITensorMPS.leftlim(cell),
            rightlim=InfiniteCylinderDMRG.ITensorMPS.rightlim(cell),
        )
    end

    function finite_cell_matches_snapshot(cell, snapshot)
        return InfiniteCylinderDMRG.ITensorMPS.leftlim(cell) ==
                snapshot.leftlim &&
            InfiniteCylinderDMRG.ITensorMPS.rightlim(cell) ==
            snapshot.rightlim &&
            all(
                inds(actual) == inds(expected) &&
                    isapprox(actual, expected; atol=0, rtol=0) for
                (actual, expected) in zip(cell, snapshot.tensors)
            )
    end

    function relabel_physical_indices!(target, source, configuration_tag)
        for site in eachindex(target)
            actual = InfiniteCylinderDMRG._physical_site_index(
                target, site, configuration_tag
            )
            desired = InfiniteCylinderDMRG._physical_site_index(
                source, site, configuration_tag
            )
            actual == desired && continue
            target[site] = replaceinds(target[site], actual => desired)
        end
        return target
    end

    function apply_fermionic_site_move(coefficients, site_count, source, target)
        dimension = 1 << site_count
        length(coefficients) == dimension || throw(
            DimensionMismatch("coefficient vector has the wrong length")
        )
        moved = similar(coefficients)
        for bits in 0:(dimension - 1)
            target_bits = 0
            occupied_targets = Int[]
            for mode in 1:site_count
                iszero(bits & (1 << (mode - 1))) && continue
                target_mode = if mode == source
                    target
                elseif source < mode <= target
                    mode - 1
                elseif target <= mode < source
                    mode + 1
                else
                    mode
                end
                target_bits |= 1 << (target_mode - 1)
                push!(occupied_targets, target_mode)
            end
            inversions = 0
            for first in 1:(length(occupied_targets) - 1)
                for second in (first + 1):length(occupied_targets)
                    inversions +=
                        occupied_targets[first] > occupied_targets[second]
                end
            end
            moved[target_bits + 1] =
                (isodd(inversions) ? -1 : 1) * coefficients[bits + 1]
        end
        return moved
    end

    cases = (
        (
            config=InfiniteCylinderConfig(;
                geometry=:paper_straight,
                Ny=2,
                x_period=1,
                filling_num=1,
                filling_den=1,
            ),
            occupied_sites=[1, 3],
        ),
        (
            config=InfiniteCylinderConfig(;
                geometry=:paper_straight,
                Ny=3,
                x_period=1,
                filling_num=2,
                filling_den=3,
            ),
            occupied_sites=[1, 5],
        ),
    )

    for (; config, occupied_sites) in cases
        _, _, psi = initial_infinite_mps(config; occupied_sites)
        bond = config.Ly
        count = InfiniteCylinderDMRG.nsites(psi)
        before = [copy(psi.AL[bond + offset]) for offset in 1:count]
        source = InfiniteCylinderDMRG.MPS(
            [psi.AL[bond + offset] for offset in 1:count]
        )
        source_snapshot = finite_cell_snapshot(source)
        translated_target = copy(source)
        configuration_tag = InfiniteCylinderDMRG._configuration_signature(config)
        if isdefined(InfiniteCylinderDMRG, :_translate_left_cell!)
            translated_result = InfiniteCylinderDMRG._translate_left_cell!(
                translated_target,
                source,
                config;
                configuration_tag,
                cutoff=0.0,
            )
            @test translated_result === translated_target
            @test finite_cell_coefficients(
                translated_target, configuration_tag
            ) ≈ apply_fermionic_ring_translation(
                finite_cell_coefficients(source, configuration_tag), config
            ) atol=1e-12 rtol=0
            @test finite_cell_matches_snapshot(source, source_snapshot)
        end

        original, translated = InfiniteCylinderDMRG._translated_left_cell(
            psi, config, bond; cutoff=0.0
        )
        original_coefficients = finite_cell_coefficients(
            original, configuration_tag
        )
        translated_coefficients = finite_cell_coefficients(
            translated, configuration_tag
        )

        @test translated_coefficients ≈ apply_fermionic_ring_translation(
            original_coefficients, config
        ) atol=1e-12 rtol=0
        @test all(
            isapprox(actual, expected; atol=0, rtol=0) for
            (actual, expected) in zip(original, before)
        )
        @test all(
            isapprox(psi.AL[bond + offset], before[offset]; atol=0, rtol=0) for
            offset in 1:count
        )
    end
    forward_config = InfiniteCylinderConfig(;
        geometry=:paper_straight,
        Ny=3,
        x_period=1,
        filling_num=2,
        filling_den=3,
    )
    _, _, forward_psi = initial_infinite_mps(
        forward_config; occupied_sites=[1, 3]
    )
    forward_bond = forward_config.Ly
    forward_source = InfiniteCylinderDMRG.MPS(
        [
            forward_psi.AL[forward_bond + offset] for
            offset in 1:InfiniteCylinderDMRG.nsites(forward_psi)
        ]
    )
    forward_tag = InfiniteCylinderDMRG._configuration_signature(forward_config)
    forward_coefficients = finite_cell_coefficients(forward_source, forward_tag)
    forward_snapshot = finite_cell_snapshot(forward_source)
    forward_target = copy(forward_source)
    forward_result = InfiniteCylinderDMRG._move_fermionic_mps_site!(
        forward_target,
        1,
        4;
        configuration_tag=forward_tag,
        cutoff=0.0,
    )
    @test forward_result === forward_target
    relabel_physical_indices!(forward_target, forward_source, forward_tag)
    @test finite_cell_coefficients(forward_target, forward_tag) ≈
        apply_fermionic_site_move(forward_coefficients, 6, 1, 4) atol=1e-12 rtol=0
    @test finite_cell_matches_snapshot(forward_source, forward_snapshot)
end

@testset "AL twisted transfer uses the left fixed point" begin
    config = InfiniteCylinderConfig(;
        geometry=:paper_straight,
        Ny=3,
        x_period=1,
        filling_num=1,
        filling_den=3,
    )
    _, _, psi = initial_infinite_mps(config; occupied_sites=[1])
    bond = config.Ly
    original, translated = InfiniteCylinderDMRG._translated_left_cell(
        psi, config, bond; cutoff=0.0
    )
    configuration_tag = InfiniteCylinderDMRG._configuration_signature(config)
    original_bra = InfiniteCylinderDMRG._dag_prime_finite_virtual_inds(
        original, configuration_tag
    )
    positive_right_transfer = InfiniteCylinderDMRG.ITensorMap(
        translated, original_bra
    )
    expected_left_transfer = transpose(positive_right_transfer)
    _, _, schmidt_index = InfiniteCylinderDMRG._schmidt_density_on_left_link(
        psi, config, bond
    )

    @test isdefined(InfiniteCylinderDMRG, :_left_twisted_transfer_map)
    if isdefined(InfiniteCylinderDMRG, :_left_twisted_transfer_map)
        left_map = InfiniteCylinderDMRG._left_twisted_transfer_map
        signature = Tuple{
            typeof(original),
            typeof(translated),
            typeof(configuration_tag),
        }
        @test hasmethod(left_map, signature)
        if hasmethod(left_map, signature)
            left_transfer = left_map(
                original, translated, configuration_tag
            )
            @test length(left_transfer.itensors) ==
                length(expected_left_transfer.itensors)
            @test all(
                isapprox(actual, expected; atol=0, rtol=0) for
                (actual, expected) in zip(
                    left_transfer.itensors, expected_left_transfer.itensors
                )
            )
            @test InfiniteCylinderDMRG.input_inds(left_transfer) ==
                InfiniteCylinderDMRG.output_inds(positive_right_transfer)
            @test InfiniteCylinderDMRG.output_inds(left_transfer) ==
                InfiniteCylinderDMRG.input_inds(positive_right_transfer)

            vector_indices = collect(
                dag(InfiniteCylinderDMRG.input_inds(left_transfer))
            )
            row_index, column_index = sort(vector_indices; by=plev, rev=true)
            @test dir(row_index) == dir(schmidt_index)
            @test dir(column_index) == dir(dag(schmidt_index))

            @test isdefined(InfiniteCylinderDMRG, :_primal_twisted_fixed_point)
            if isdefined(InfiniteCylinderDMRG, :_primal_twisted_fixed_point)
                dual_fixed_point = ITensor(
                    ComplexF64, QN(), vector_indices...
                )
                primal_fixed_point = InfiniteCylinderDMRG._primal_twisted_fixed_point(
                    dual_fixed_point
                )
                primal_row, primal_column =
                    InfiniteCylinderDMRG._twisted_matrix_indices(
                        primal_fixed_point, schmidt_index
                    )
                @test dir(primal_row) == dir(dag(schmidt_index))
                @test dir(primal_column) == dir(schmidt_index)
            end
        end
    end
end

@testset "twisted transfer charge blocks normalize independently" begin
    @test isdefined(InfiniteCylinderDMRG, :_twisted_block_matrix)
    @test isdefined(InfiniteCylinderDMRG, :_twisted_block_matrices)
    if isdefined(InfiniteCylinderDMRG, :_twisted_block_matrix)
        schmidt_index = Index(
            [QN("Nf", 0) => 2, QN("Nf", 1) => 1],
            "Link,Schmidt",
        )
        row_index = prime(dag(schmidt_index))
        column_index = schmidt_index
        fixed_point = ITensor(
            ComplexF64,
            QN(),
            row_index,
            column_index,
        )
        fixed_point[
            row_index => 1:2,
            column_index => 1:2,
        ] = 2im * ComplexF64[0 1; -1 0]
        fixed_point[row_index => 3, column_index => 3] = -3.0

        if isdefined(InfiniteCylinderDMRG, :_primal_twisted_fixed_point)
            dual_fixed_point = ITensors.setinds(
                fixed_point,
                dag.([row_index, column_index]),
            )
            primal_fixed_point = InfiniteCylinderDMRG._primal_twisted_fixed_point(
                dual_fixed_point
            )
            @test Array(primal_fixed_point, row_index, column_index) ==
                Array(fixed_point, row_index, column_index)
        end

        matrix = InfiniteCylinderDMRG._twisted_block_matrix(
            fixed_point,
            schmidt_index,
            2,
        )
        identity_matrix = Matrix{ComplexF64}(I, 3, 3)

        if isdefined(InfiniteCylinderDMRG, :_twisted_block_matrices)
            matrices = InfiniteCylinderDMRG._twisted_block_matrices(
                fixed_point,
                schmidt_index,
                2,
            )
            @test matrices isa Vector{Matrix{ComplexF64}}
            @test size.(matrices) == [(2, 2), (1, 1)]
            @test cat(matrices...; dims=(1, 2)) ≈ matrix atol=1e-14 rtol=0
            @test sum(length, matrices) == 5
        end

        @test norm(matrix[1:2, 1:2]) / sqrt(2) ≈ 1.0 atol=1e-12
        @test norm(matrix[3:3, 3:3]) ≈ 1.0 atol=1e-12
        @test norm(matrix^2 - identity_matrix) ≤ 1e-12
        @test iszero(matrix[1:2, 3:3])
        @test iszero(matrix[3:3, 1:2])
    end
end

@testset "repeated identical QN blocks form one logical sector" begin
    @test isdefined(InfiniteCylinderDMRG, :_neutral_qn_block_groups)
    @test isdefined(InfiniteCylinderDMRG, :_random_twisted_group)
    @test isdefined(InfiniteCylinderDMRG, :_project_twisted_group)
    @test isdefined(InfiniteCylinderDMRG, :_twisted_group_dimension)
    schmidt_index = Index(
        [
            QN("Nf", 0) => 1,
            QN("Nf", 0) => 2,
            QN("Nf", 1) => 1,
        ],
        "Link,RepeatedSchmidt",
    )
    right_index = dag(Index(
        [
            QN("Nf", 0) => 1,
            QN("Nf", 0) => 2,
            QN("Nf", 1) => 1,
        ],
        "Link,RepeatedRight",
    ))
    center = ITensor(ComplexF64, QN(), schmidt_index, right_index)
    zero_charge_center = ComplexF64[
        1 2 3
        4 5 6
        7 8 10
    ]
    center[schmidt_index => 1, right_index => 1] = zero_charge_center[1, 1]
    center[schmidt_index => 1, right_index => 2:3] = zero_charge_center[1, 2:3]
    center[schmidt_index => 2:3, right_index => 1] = zero_charge_center[2:3, 1]
    center[schmidt_index => 2:3, right_index => 2:3] = zero_charge_center[2:3, 2:3]
    center[schmidt_index => 4, right_index => 4] = 2im
    dense_center = Matrix{ComplexF64}(
        Array(center, schmidt_index, right_index)
    )
    dense_density_reference = dense_center * dense_center'

    density_result = try
        InfiniteCylinderDMRG._schmidt_density_blocks(
            center, schmidt_index, right_index
        )
    catch error
        error
    end
    @test density_result isa Tuple
    if density_result isa Tuple
        density_blocks, block_charges = density_result
        @test block_charges == [0, 1]
        @test size.(density_blocks) == [(3, 3), (1, 1)]
        @test cat(density_blocks...; dims=(1, 2)) ≈
            dense_density_reference atol=1e-13 rtol=0
        @test sum(length, density_blocks) == 10
    end

    row_index = prime(dag(schmidt_index))
    fixed_point = ITensor(
        ComplexF64,
        QN(),
        row_index,
        schmidt_index,
    )
    zero_charge_translation = ComplexF64[
        0 1 0
        1 0 0
        0 0 -1
    ]
    fixed_point[row_index => 1, schmidt_index => 1] = zero_charge_translation[1, 1]
    fixed_point[row_index => 1, schmidt_index => 2:3] = zero_charge_translation[1, 2:3]
    fixed_point[row_index => 2:3, schmidt_index => 1] = zero_charge_translation[2:3, 1]
    fixed_point[row_index => 2:3, schmidt_index => 2:3] = zero_charge_translation[2:3, 2:3]
    fixed_point[row_index => 4, schmidt_index => 4] = -2.0

    if all(
        helper -> isdefined(InfiniteCylinderDMRG, helper),
        (
            :_neutral_qn_block_groups,
            :_random_twisted_group,
            :_project_twisted_group,
            :_twisted_group_dimension,
        ),
    )
        qn_groups, structural_groups =
            InfiniteCylinderDMRG._neutral_qn_block_groups(
                row_index, schmidt_index
            )
        @test qn_groups == [[1, 2], [3]]
        @test [
            [Tuple(Int.(Tuple(block))) for block in group] for
            group in structural_groups
        ] == [
            [(1, 1), (2, 1), (1, 2), (2, 2)],
            [(3, 3)],
        ]
        @test InfiniteCylinderDMRG._twisted_group_dimension(
            [row_index, schmidt_index], structural_groups[1]
        ) == 9

        initial = InfiniteCylinderDMRG._random_twisted_group(
            Random.Xoshiro(0x7d3),
            [row_index, schmidt_index],
            structural_groups[1],
        )
        ordered_initial = ITensors.permute(initial, row_index, schmidt_index)
        @test all(structural_groups[1]) do block
            value = ordered_initial[block]
            !isnothing(value) && norm(value) > 0
        end
        @test isnothing(ordered_initial[structural_groups[2][1]])

        projected = InfiniteCylinderDMRG._project_twisted_group(
            fixed_point,
            [row_index, schmidt_index],
            structural_groups[1],
        )
        ordered_projected = ITensors.permute(
            projected, row_index, schmidt_index
        )
        @test all(structural_groups[1]) do block
            projected_value = ordered_projected[block]
            source_value = fixed_point[block]
            !isnothing(projected_value) &&
                !isnothing(source_value) &&
                Array(projected_value) == Array(source_value)
        end
        @test isnothing(ordered_projected[structural_groups[2][1]])
    end

    translation_result = try
        InfiniteCylinderDMRG._twisted_block_matrices(
            fixed_point, schmidt_index, 2
        )
    catch error
        error
    end
    @test translation_result isa Vector{Matrix{ComplexF64}}
    if translation_result isa Vector{Matrix{ComplexF64}}
        @test size.(translation_result) == [(3, 3), (1, 1)]
        @test sum(length, translation_result) == 10
        @test all(
            block -> norm(block' * block - I) <= 1e-12,
            translation_result,
        )
        @test all(
            block -> norm(block^2 - I) <= 1e-12,
            translation_result,
        )
    end
end

@testset "Ny=2 tensor momentum smoke" begin
    config = InfiniteCylinderConfig(;
        geometry=:paper_straight,
        Ny=2,
        x_period=1,
        filling_num=1,
        filling_den=1,
    )
    @test isdefined(InfiniteCylinderDMRG, :_mixed_symmetry_eigenvalue_residual)
    if isdefined(InfiniteCylinderDMRG, :_mixed_symmetry_eigenvalue_residual)
        residual = InfiniteCylinderDMRG._mixed_symmetry_eigenvalue_residual
        common_phase = cis(0.37)

        @test residual(ComplexF64[common_phase, common_phase]) ≈ 0.0 atol=1e-14
        @test residual(ComplexF64[0.9, 0.9]) ≈ 0.1 atol=1e-14
        @test residual(ComplexF64[1, cis(0.1)]) > 0
        @test_throws ArgumentError residual(ComplexF64[])
        @test_throws ArgumentError residual(ComplexF64[NaN + 0im])
        @test_throws ArgumentError residual(ComplexF64[Inf + 0im])
    end
    @test isdefined(InfiniteCylinderDMRG, :_dominant_mixed_eigenpair)
    if isdefined(InfiniteCylinderDMRG, :_dominant_mixed_eigenpair)
        select_dominant = InfiniteCylinderDMRG._dominant_mixed_eigenpair
        unique_pair = select_dominant(
            ComplexF64[1, 0.8], [:first, :second], 1e-8
        )
        @test unique_pair.value == 1
        @test unique_pair.vector == :first
        @test unique_pair.gap ≈ 0.2 atol=1e-14

        unsorted_pair = select_dominant(
            ComplexF64[0.8, 1], [:subleading, :dominant], 1e-8
        )
        @test unsorted_pair.value == 1
        @test unsorted_pair.vector == :dominant

        one_pair = select_dominant(ComplexF64[cis(0.2)], [:only], 1e-8)
        @test one_pair.vector == :only
        @test isinf(one_pair.gap)
        @test isdefined(InfiniteCylinderDMRG, :MomentumConstructionError)
        if isdefined(InfiniteCylinderDMRG, :MomentumConstructionError)
            construction_error = InfiniteCylinderDMRG.MomentumConstructionError(
                "expected momentum construction failure"
            )
            @test sprint(showerror, construction_error) ==
                "expected momentum construction failure"
            @test_throws InfiniteCylinderDMRG.MomentumConstructionError select_dominant(
                ComplexF64[1, 1 - 1e-12], [:first, :second], 1e-8
            )
        end
        @test_throws ArgumentError select_dominant(
            ComplexF64[], Symbol[], 1e-8
        )
        @test_throws ArgumentError select_dominant(
            ComplexF64[Inf + 0im], [:only], 1e-8
        )
    end
    _, _, psi = initial_infinite_mps(config; occupied_sites=[1, 3])
    _, _, schmidt_index =
        InfiniteCylinderDMRG._schmidt_density_on_left_link(
            psi, config, config.Ly
        )
    virtual = InfiniteCylinderDMRG._twisted_virtual_translation(
        psi,
        config,
        config.Ly;
        schmidt_index,
        tol=1e-8,
        rng=Random.Xoshiro(0x7a57),
        cutoff=0.0,
    )
    result = momentum_entanglement_data(
        psi,
        config;
        cut_x=1,
        residual_tol=1e-8,
        transfer_tol=1e-8,
        rng=Random.Xoshiro(0x7a57),
    )

    @test result isa MomentumEntanglementData
    @test !isempty(result.levels)
    @test !hasproperty(virtual, :matrix)
    @test hasproperty(virtual, :blocks)
    if hasproperty(virtual, :blocks)
        @test virtual.blocks isa Vector{Matrix{ComplexF64}}
        @test !isempty(virtual.blocks)
        @test sum(length, virtual.blocks) ==
            sum(block -> size(block, 1)^2, virtual.blocks)
    end
    @test hasproperty(virtual, :eigenvalues)
    @test hasproperty(virtual, :eigenpair_residual)
    @test hasproperty(virtual, :symmetry_residual)
    @test hasproperty(virtual, :requested_eigenpairs)
    @test hasproperty(virtual, :dominant_gaps)
    if all(
        property -> hasproperty(virtual, property),
        (:eigenvalues, :eigenpair_residual, :symmetry_residual),
    )
        @test virtual.eigenvalues isa Vector{ComplexF64}
        @test virtual.residual ≈ max(
            virtual.eigenpair_residual, virtual.symmetry_residual
        ) atol=1e-14 rtol=0
        @test result.validation.transfer_residual ≈ virtual.residual atol=1e-14 rtol=0
    end
    if hasproperty(virtual, :requested_eigenpairs) &&
            hasproperty(virtual, :dominant_gaps)
        @test length(virtual.requested_eigenpairs) == length(virtual.eigenvalues)
        @test length(virtual.dominant_gaps) == length(virtual.eigenvalues)
        @test all(==(1), virtual.requested_eigenpairs)
        @test all(isinf, virtual.dominant_gaps)
    end
    @test result.reason ==
        "valid charge-resolved Schmidt spectrum with sector-relative Δky labels"
    @test result.valid
    @test result.validation.valid
    @test result.momentum_convention ===
        :sector_relative_canonical_cyclic_orbit
    @test !isempty(result.sector_ky_offsets)
    @test all(level -> !ismissing(level.ky_index), result.levels)
    @test maximum(level.phase_residual for level in result.levels) <= 1e-8

    schmidt_lowered = join(
        sprint(show, lowered) for lowered in code_lowered(
            InfiniteCylinderDMRG._schmidt_density_on_left_link,
            Tuple{typeof(psi),typeof(config),Int},
        )
    )
    translation_lowered = join(
        sprint(show, lowered) for lowered in code_lowered(
            InfiniteCylinderDMRG._twisted_virtual_translation,
            Tuple{typeof(psi),typeof(config),Int},
        )
    )
    momentum_lowered = join(
        sprint(show, lowered) for lowered in code_lowered(
            momentum_entanglement_data,
            Tuple{typeof(psi),typeof(config)},
        )
    )
    @test !occursin("Array", schmidt_lowered)
    @test !occursin("_twisted_block_matrix", translation_lowered)
    @test !occursin("validate_momentum_blocks", momentum_lowered)

    _, _, broken_psi = initial_infinite_mps(config; occupied_sites=[1, 2])
    broken_result = momentum_entanglement_data(
        broken_psi,
        config;
        cut_x=1,
        residual_tol=1e-8,
        transfer_tol=1e-8,
        rng=Random.Xoshiro(0xbad5),
    )
    @test !broken_result.valid
    @test !broken_result.validation.valid
    @test broken_result.momentum_convention === :invalid
    @test isempty(broken_result.sector_ky_offsets)
    @test all(level -> ismissing(level.ky_index), broken_result.levels)
    @test occursin(
        "transverse-translation symmetry",
        lowercase(broken_result.reason),
    ) || occursin("mixed-symmetry", lowercase(broken_result.reason))
    @test isinf(broken_result.validation.transfer_residual)

    legacy = InfiniteCylinderConfig(; Ly=4, x_period=1, filling_num=1, filling_den=1)
    _, _, legacy_psi = initial_infinite_mps(legacy; occupied_sites=[1, 3])
    @test_throws ArgumentError momentum_entanglement_data(
        legacy_psi,
        legacy;
        cut_x=1,
    )
end

@testset "Ny=3 tensor translation fixes the positive phase direction" begin
    config = InfiniteCylinderConfig(;
        geometry=:paper_straight,
        Ny=3,
        x_period=1,
        filling_num=1,
        filling_den=3,
    )
    _, _, product = initial_infinite_mps(config; occupied_sites=[1])
    cell = InfiniteCylinderDMRG.MPS(
        [product.AL[site] for site in 1:InfiniteCylinderDMRG.nsites(product)]
    )
    configuration_tag = InfiniteCylinderDMRG._configuration_signature(config)
    physical_indices = [
        InfiniteCylinderDMRG._physical_site_index(
            cell, site, configuration_tag
        ) for site in eachindex(cell)
    ]
    left_boundary = only(filter(
        index -> !hastags(index, configuration_tag),
        uniqueinds(cell[1], cell[2]),
    ))
    right_boundary = only(filter(
        index -> !hastags(index, configuration_tag),
        uniqueinds(cell[end], cell[end - 1]),
    ))

    expected_phase = cis(2pi / config.Ny)
    wavefunction = ITensor(ComplexF64, QN(), physical_indices...)
    for ring in 0:(config.Ny - 1)
        selectors = [
            physical_indices[site] => (site == 2ring + 1 ? 2 : 1) for
            site in eachindex(physical_indices)
        ]
        wavefunction[selectors...] =
            cis(-2pi * ring / config.Ny) / sqrt(config.Ny)
    end
    combined = wavefunction *
        onehot(left_boundary => 1) *
        onehot(right_boundary => 1)
    link_tags = TagSet[
        tags(only(commoninds(cell[site], cell[site + 1]))) for
        site in 1:(length(cell) - 1)
    ]
    original = InfiniteCylinderDMRG.MPS(
        combined,
        physical_indices;
        leftinds=left_boundary,
        orthocenter=length(physical_indices),
        tags=link_tags,
        cutoff=0.0,
    )
    translated = copy(original)
    InfiniteCylinderDMRG._translate_left_cell!(
        translated,
        original,
        config;
        configuration_tag,
        cutoff=0.0,
    )

    function finite_cell_coefficients(state)
        physical = [
            InfiniteCylinderDMRG._physical_site_index(
                state, site, configuration_tag
            ) for site in eachindex(state)
        ]
        state_tensor = foldl(*, state)
        boundary = filter(
            index -> !hastags(index, configuration_tag),
            collect(inds(state_tensor)),
        )
        return vec(Array(state_tensor, physical..., boundary...))
    end

    source = finite_cell_coefficients(original)
    image = finite_cell_coefficients(translated)
    measured_phase = dot(source, image) / dot(source, source)
    @test image ≈ expected_phase * source atol=1e-13 rtol=0
    @test measured_phase ≈ expected_phase atol=1e-13 rtol=0
    @test !isapprox(
        measured_phase, conj(expected_phase); atol=1e-13, rtol=0
    )
end

@testset "entangled multi-QN tensor smoke fails momentum labels explicitly" begin
    config = InfiniteCylinderConfig(;
        geometry=:paper_straight,
        Ny=2,
        x_period=1,
        filling_num=1,
        filling_den=1,
    )
    canonical = fetch(@async begin
        Random.seed!(0x7d4)
        sites, _, product = initial_infinite_mps(
            config; occupied_sites=[1, 3]
        )
        model = CylinderModelParams(;
            t1=1.0,
            t3=0.2,
            V1=1.0,
            V2=0.0,
            V3=0.0,
        )
        H = build_infinite_mpo(config, model, sites)
        expanded = expand_subspace(product, H, 4; cutoff=1e-10)
        step = vumps_iteration(H, expanded; vumps_tol=1e-5)
        InfiniteCylinderDMRG._canonicalize_vumps_state(step.psi)
    end)

    density_blocks, block_charges, _ =
        InfiniteCylinderDMRG._schmidt_density_on_left_link(
            canonical, config, config.Ly
        )
    sector_weights = Float64.(real.(tr.(density_blocks)))
    significant_charges = [
        charge for (charge, weight) in zip(block_charges, sector_weights) if
        weight > 0.05
    ]
    @test maximum(link_dimensions(canonical)) > 1
    @test length(density_blocks) >= 2
    @test length(unique(block_charges)) >= 2
    @test length(unique(significant_charges)) >= 2
    @test sum(sector_weights) ≈ 1.0 atol=1e-12 rtol=0

    result = momentum_entanglement_data(
        canonical,
        config;
        cut_x=1,
        residual_tol=1e-7,
        transfer_tol=1e-7,
        rng_seed=0x7d4d4,
    )
    expected_level_charges = reduce(
        vcat,
        [fill(charge, size(block, 1)) for
         (charge, block) in zip(block_charges, density_blocks)],
    )
    @test length(result.levels) == sum(size(block, 1) for block in density_blocks)
    @test [level.raw_charge for level in result.levels] == expected_level_charges
    @test !result.valid
    @test !result.validation.valid
    @test all(level -> ismissing(level.ky_index), result.levels)
    @test occursin("twisted-transfer", lowercase(result.reason))
    @test occursin("residual", lowercase(result.reason))
    @test isinf(result.validation.transfer_residual)
end
