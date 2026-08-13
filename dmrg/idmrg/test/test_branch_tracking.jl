using Test
using ITensorInfiniteMPS: InfMPS
using InfiniteCylinderDMRG

@testset "branch metrics" begin
    @test bhattacharyya_distance(Dict(0 => 1.0), Dict(0 => 1.0)) ≈ 0.0
    @test bhattacharyya_distance(Dict(0 => 1.0), Dict(1 => 1.0)) ≈ 1.0
    @test bhattacharyya_distance(
        Dict(0 => 2.0, 1 => 2.0),
        Dict(0 => 1.0),
    ) ≈ 1 - sqrt(0.5)
    @test bhattacharyya_distance(
        Dict(0 => 0.0, 1 => 3.0),
        Dict(1 => 7.0, 2 => 0.0),
    ) ≈ 0.0

    left = Dict(0 => 0.75, 1 => 0.25)
    right = Dict(2 => 3.0, 3 => 1.0)
    left_before = copy(left)
    right_before = copy(right)
    aligned = best_sector_shift(left, right)
    @test aligned.shift == -2
    @test aligned.distance ≈ 0.0
    @test left == left_before
    @test right == right_before

    tied = best_sector_shift(Dict(0 => 0.5, 2 => 0.5), Dict(1 => 1.0))
    @test tied.shift == -1
    @test tied.distance ≈ 1 - inv(sqrt(2))

    levels1 = [
        (raw_charge=0, entanglement_energy=0.1, probability=0.5),
        (raw_charge=1, entanglement_energy=0.3, probability=0.3),
        (raw_charge=0, entanglement_energy=0.7, probability=0.2),
    ]
    levels2 = [
        (raw_charge=0, entanglement_energy=0.2, probability=0.45),
        (raw_charge=1, entanglement_energy=0.4, probability=0.35),
        (raw_charge=0, entanglement_energy=1.1, probability=0.2),
    ]
    spectrum = charge_matched_spectrum_distance(levels1, levels2)
    @test spectrum.valid
    @test spectrum.matched == 3
    @test spectrum.requested == 20
    @test spectrum.distance ≈ 0.4
    @test occursin("3 of 20", spectrum.reason)

    limited = charge_matched_spectrum_distance(levels1, levels2; nlevels=2)
    @test limited.valid
    @test limited.matched == 2
    @test limited.requested == 2
    @test limited.distance ≈ 0.1

    unequal_multiplicities = charge_matched_spectrum_distance(
        [
            (raw_charge=0, entanglement_energy=0.1, probability=0.5),
            (raw_charge=0, entanglement_energy=0.7, probability=0.5),
        ],
        [
            (raw_charge=0, entanglement_energy=0.2, probability=0.4),
            (raw_charge=0, entanglement_energy=1.1, probability=0.3),
            (raw_charge=0, entanglement_energy=1.5, probability=0.3),
        ],
    )
    @test unequal_multiplicities.valid
    @test unequal_multiplicities.matched == 2
    @test unequal_multiplicities.requested == 20
    @test unequal_multiplicities.distance ≈ 0.4

    unmatched = charge_matched_spectrum_distance(
        [(raw_charge=0, entanglement_energy=0.1, probability=1.0)],
        [(raw_charge=1, entanglement_energy=0.1, probability=1.0)],
    )
    @test !unmatched.valid
    @test unmatched.matched == 0
    @test unmatched.requested == 20
    @test isnan(unmatched.distance)
    @test occursin("no charge-matched", unmatched.reason)

    flags = branch_flags(
        ;
        fidelity=0.9,
        charge_step=1 / 3,
        sector_distance=0.2,
        spectrum_distance=0.6,
        energy_crossing=false,
    )
    @test flags.fidelity
    @test flags.charge
    @test flags.sector
    @test flags.spectrum
    @test !flags.energy_crossing

    boundary_flags = branch_flags(
        ;
        fidelity=1 - 0.5e-3,
        charge_step=-0.1,
        sector_distance=0.1,
        spectrum_distance=0.5,
        energy_crossing=true,
    )
    @test !boundary_flags.fidelity
    @test !boundary_flags.charge
    @test !boundary_flags.sector
    @test !boundary_flags.spectrum
    @test boundary_flags.energy_crossing

    for weights in (
        Dict{Int,Float64}(),
        Dict(0 => 0.0),
        Dict(0 => -1.0),
        Dict(0 => NaN),
        Dict(0 => Inf),
    )
        @test_throws ArgumentError bhattacharyya_distance(weights, Dict(0 => 1.0))
        @test_throws ArgumentError best_sector_shift(weights, Dict(0 => 1.0))
    end
    @test_throws ArgumentError bhattacharyya_distance(
        Dict(0.5 => 1.0), Dict(0 => 1.0)
    )
    @test_throws ArgumentError best_sector_shift(Dict(true => 1.0), Dict(0 => 1.0))

    invalid_levels = (
        [(raw_charge=0, entanglement_energy=NaN, probability=1.0)],
        [(raw_charge=0, entanglement_energy=Inf, probability=1.0)],
        [(raw_charge=0, entanglement_energy=0.0, probability=NaN)],
        [(raw_charge=0, entanglement_energy=0.0, probability=Inf)],
        [(raw_charge=0, entanglement_energy=0.0, probability=-0.1)],
        [(raw_charge=true, entanglement_energy=0.0, probability=1.0)],
    )
    for levels in invalid_levels
        @test_throws ArgumentError charge_matched_spectrum_distance(levels, levels2)
    end
    @test_throws ArgumentError charge_matched_spectrum_distance(
        levels1, levels2; nlevels=0
    )
    @test_throws ArgumentError charge_matched_spectrum_distance(
        levels1, levels2; nlevels=true
    )

    invalid_flag_arguments = (
        (; fidelity=NaN),
        (; fidelity=-0.1),
        (; fidelity=1.1),
        (; charge_step=Inf),
        (; sector_distance=-0.1),
        (; spectrum_distance=NaN),
        (; fidelity_drop_tol=-1.0),
        (; charge_jump_tol=Inf),
        (; sector_tol=-0.1),
        (; spectrum_tol=NaN),
        (; energy_crossing=1),
    )
    defaults = (
        fidelity=1.0,
        charge_step=0.0,
        sector_distance=0.0,
        spectrum_distance=0.0,
        energy_crossing=false,
    )
    for overrides in invalid_flag_arguments
        arguments = merge(defaults, overrides)
        @test_throws ArgumentError branch_flags(; arguments...)
    end
end

@testset "candidate selection" begin
    cfg = InfiniteCylinderConfig(; Ly=6, x_period=1)
    sites, _, previous = initial_infinite_mps(cfg; occupied_sites=[1])
    previous_entanglement = [entanglement_data(previous, cfg; cut_x=1)]
    second_init(n::Integer) = mod1(n, sites_per_cell(cfg)) == 2 ? "Occ" : "Emp"
    different = InfMPS(sites, second_init)
    different_entanglement = [entanglement_data(different, cfg; cut_x=1)]

    candidates = [
        BranchCandidate(
            "invalid-low",
            different,
            normalize_energy(cfg, -12.0),
            different_entanglement,
            true,
            false,
            "observable validation failed",
        ),
        BranchCandidate(
            "adiabatic",
            previous,
            normalize_energy(cfg, -6.0),
            previous_entanglement,
            true,
            true,
            "valid",
        ),
        BranchCandidate(
            "ground",
            different,
            normalize_energy(cfg, -9.0),
            different_entanglement,
            true,
            true,
            "valid",
        ),
        BranchCandidate(
            "unconverged-low",
            previous,
            normalize_energy(cfg, -18.0),
            previous_entanglement,
            false,
            true,
            "not converged",
        ),
    ]
    raw_before = [candidate.entanglement for candidate in candidates]

    ground = select_candidate(candidates; mode=:ground)
    @test ground isa CandidateSelection
    @test ground.mode == :ground
    @test ground.selected_index == 3
    @test ground.candidates == candidates
    @test all(isnothing, ground.fidelities)
    @test occursin("minimum", ground.reason)

    adiabatic = select_candidate(
        candidates;
        mode=:adiabatic,
        previous_state=previous,
        config=cfg,
        tol=1e-10,
    )
    @test adiabatic.mode == :adiabatic
    @test adiabatic.selected_index == 2
    @test adiabatic.candidates == candidates
    @test adiabatic.fidelities[1] === nothing
    @test adiabatic.fidelities[2] isa MixedTransferFidelity
    @test adiabatic.fidelities[2].valid
    @test adiabatic.fidelities[2].fidelity_cell ≈ 1.0 atol=1e-10
    @test adiabatic.fidelities[3] isa MixedTransferFidelity
    @test adiabatic.fidelities[4] === nothing
    @test occursin("maximum", adiabatic.reason)

    @test all(
        candidate.entanglement === before for
        (candidate, before) in zip(candidates, raw_before)
    )
    @test ground.candidates !== candidates

    tied_ground = select_candidate(
        [candidates[2], candidates[2]];
        mode=:ground,
    )
    @test tied_ground.selected_index == 1
    tied_adiabatic = select_candidate(
        [candidates[2], candidates[2]];
        mode=:adiabatic,
        previous_state=previous,
        config=cfg,
        tol=1e-10,
    )
    @test tied_adiabatic.selected_index == 1

    @test_throws ArgumentError select_candidate(BranchCandidate[]; mode=:ground)
    @test_throws ArgumentError select_candidate(
        [candidates[1], candidates[4]]; mode=:ground
    )
    @test_throws ArgumentError select_candidate(candidates; mode=:adiabatic)
    @test_throws ArgumentError select_candidate(
        candidates; mode=:adiabatic, previous_state=previous
    )
    @test_throws ArgumentError select_candidate(candidates; mode=:unsupported)
    @test_throws ArgumentError BranchCandidate(
        "",
        previous,
        normalize_energy(cfg, -6.0),
        previous_entanglement,
        true,
        true,
        "valid",
    )
end

@testset "mixed transfer fidelity" begin
    cfg = InfiniteCylinderConfig(; Ly=6, x_period=1)
    sites, _, psi1 = initial_infinite_mps(cfg; occupied_sites=[1])
    same = mixed_transfer_fidelity(psi1, psi1, cfg; tol=1e-10)
    @test same isa MixedTransferFidelity
    @test same.valid
    @test same.converged
    @test same.fidelity_cell ≈ 1.0 atol=1e-10
    @test abs(same.mixed_value) ≈ 1.0 atol=1e-10
    @test abs(same.self_value1) ≈ 1.0 atol=1e-10
    @test abs(same.self_value2) ≈ 1.0 atol=1e-10
    @test same.mixed_residual <= 1e-10
    @test same.self_residual1 <= 1e-10
    @test same.self_residual2 <= 1e-10
    @test occursin("valid", same.reason)

    second_init(n::Integer) = mod1(n, sites_per_cell(cfg)) == 2 ? "Occ" : "Emp"
    psi2 = InfMPS(sites, second_init)
    different = mixed_transfer_fidelity(psi1, psi2, cfg; tol=1e-10)
    @test different isa MixedTransferFidelity
    if different.valid
        @test 0 <= different.fidelity_cell < 1 - 1e-10
    else
        @test !different.converged || occursin("invalid", lowercase(different.reason))
        @test isnan(different.fidelity_cell)
    end

    _, _, fresh_indices = initial_infinite_mps(cfg; occupied_sites=[1])
    @test_throws ArgumentError mixed_transfer_fidelity(psi1, fresh_indices, cfg)
    wrong_geometry = InfiniteCylinderConfig(; Ly=2, x_period=3)
    @test_throws ArgumentError mixed_transfer_fidelity(psi1, psi1, wrong_geometry)
    @test_throws ArgumentError mixed_transfer_fidelity(psi1, psi1, cfg; tol=0.0)

    invalid = InfiniteCylinderDMRG._mixed_transfer_result(
        ComplexF64(NaN),
        1.0 + 0.0im,
        1.0 + 0.0im,
        0.0,
        0.0,
        0.0,
        true;
        residual_tolerance=1e-10,
    )
    @test !invalid.valid
    @test invalid.converged
    @test isnan(invalid.fidelity_cell)
    @test occursin("finite", invalid.reason)

    over_normalized = InfiniteCylinderDMRG._mixed_transfer_result(
        1.01 + 0.0im,
        1.0 + 0.0im,
        1.0 + 0.0im,
        0.0,
        0.0,
        0.0,
        true;
        residual_tolerance=1e-10,
    )
    @test !over_normalized.valid
    @test over_normalized.converged
    @test isnan(over_normalized.fidelity_cell)
    @test occursin("outside [0, 1]", over_normalized.reason)

    roundoff = InfiniteCylinderDMRG._mixed_transfer_result(
        1.0 + 10eps(Float64) + 0.0im,
        1.0 + 0.0im,
        1.0 + 0.0im,
        0.0,
        0.0,
        0.0,
        true;
        residual_tolerance=1e-10,
    )
    @test roundoff.valid
    @test roundoff.fidelity_cell == 1.0
end
