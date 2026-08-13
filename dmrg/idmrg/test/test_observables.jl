using Test
using ITensors
using ITensorMPS
using InfiniteCylinderDMRG
using Random

@testset "observable helpers" begin
    cfg = InfiniteCylinderConfig(; Ly=6, x_period=1)
    energy = normalize_energy(cfg, -6.0)
    @test energy.per_cell == -6.0
    @test energy.per_x == -6.0
    @test energy.per_unit_cell == -2.0
    @test energy.per_site == -1.0
    for nonfinite_energy in (
        NaN,
        Inf,
        -Inf,
        BigFloat(floatmax(Float64)) * 2,
    )
        @test_throws ArgumentError normalize_energy(cfg, nonfinite_energy)
    end

    sites, _, psi = initial_infinite_mps(cfg)
    signature = InfiniteCylinderDMRG._configuration_signature(cfg)
    @test all(hastags(site, signature) for site in sites)
    @test all(!hastags(site, "Site") for site in sites)
    flux_cfg = InfiniteCylinderConfig(; Ly=6, x_period=1, phi_y=0.37π)
    @test InfiniteCylinderDMRG._configuration_signature(flux_cfg) == signature
    oversized_signature_cfg = InfiniteCylinderConfig(
        2,
        1,
        typemax(Int),
        typemax(Int),
        0.0,
    )
    @test_throws ArgumentError initial_infinite_mps(oversized_signature_cfg)

    params = CylinderModelParams(; t1=1.0, t3=0.2, V1=1.0, V2=0.4, V3=0.3)
    H = build_infinite_mpo(cfg, params, sites)
    raw_energy = sum(real, expect(psi, H))
    @test energy_data(psi, H, cfg) == normalize_energy(cfg, raw_energy)

    densities = density_data(psi, cfg)
    @test length(densities) == sites_per_cell(cfg)
    @test [row.site for row in densities] == collect(1:sites_per_cell(cfg))
    @test [(row.x, row.y) for row in densities] ==
        [site_coordinates(cfg, site) for site in 1:sites_per_cell(cfg)]
    @test [row.density for row in densities] ==
        [site in default_occupied_sites(cfg) ? 1.0 : 0.0 for site in 1:sites_per_cell(cfg)]
    size_mismatch = InfiniteCylinderConfig(; Ly=6, x_period=3)
    @test_throws ArgumentError energy_data(psi, H, size_mismatch)
    @test_throws ArgumentError density_data(psi, size_mismatch)

    same_size_mismatches = (
        InfiniteCylinderConfig(; Ly=2, x_period=3),
        InfiniteCylinderConfig(6, 1, 2, 3, 0.0),
    )
    for mismatch in same_size_mismatches
        @test_throws ArgumentError energy_data(psi, H, mismatch)
        @test_throws ArgumentError density_data(psi, mismatch)
        @test_throws ArgumentError entanglement_data(psi, mismatch; cut_x=1)
    end

    entanglement = entanglement_data(psi, cfg; cut_x=1)
    @test isapprox(
        sum(level.probability for level in entanglement.levels),
        1.0;
        atol=1e-12,
    )
    @test isapprox(entanglement.entropy, 0.0; atol=1e-12)
    @test !isempty(entanglement.sectors)
    @test all(
        level.physical_charge == level.raw_charge / charge_scale(cfg) for
        level in entanglement.levels
    )
    @test entanglement.raw_schmidt_polarization == sum(
        level.probability * level.physical_charge for level in entanglement.levels
    )
    @test sum(sector.weight for sector in entanglement.sectors) ≈ 1.0 atol=1e-12

    cfg3 = InfiniteCylinderConfig(; Ly=6, x_period=3)
    _, _, charged_product =
        initial_infinite_mps(cfg3; occupied_sites=[1, 2, 3])
    for (cut_x, raw_charge) in zip(1:3, (12, 6, 0))
        charged_entanglement =
            entanglement_data(charged_product, cfg3; cut_x)
        @test only(charged_entanglement.levels).raw_charge == raw_charge
        @test only(charged_entanglement.levels).physical_charge ==
            raw_charge / charge_scale(cfg3)
        @test only(charged_entanglement.sectors).raw_charge == raw_charge
    end

    @test_throws ArgumentError entanglement_data(psi, cfg; cut_x=0)
    @test_throws ArgumentError entanglement_data(psi, cfg; cut_x=2)
    @test_throws ArgumentError entanglement_data(psi, cfg; cut_x=true)
    @test_throws ArgumentError InfiniteCylinderDMRG._normalized_schmidt_probabilities(
        Float64[]
    )
    @test_throws ArgumentError InfiniteCylinderDMRG._normalized_schmidt_probabilities(
        [0.0]
    )
    @test_throws ArgumentError InfiniteCylinderDMRG._normalized_schmidt_probabilities(
        [NaN]
    )
    @test_throws ArgumentError InfiniteCylinderDMRG._normalized_schmidt_probabilities(
        [-1.0]
    )
    schmidt_index = Index(
        [
            QN("Nf", -6, -1) => 2,
            QN("Nf", 0, -1) => 1,
            QN("Nf", 6, -1) => 3,
        ],
        "Schmidt";
        dir=ITensors.Out,
    )
    level_qns = InfiniteCylinderDMRG._schmidt_level_qns(schmidt_index)
    @test val.(level_qns, "Nf") == [-6, -6, 0, 6, 6, 6]
    mapped_weights = InfiniteCylinderDMRG._schmidt_sector_weights(
        level_qns,
        [0.1, 0.2, 0.3, 0.05, 0.15, 0.2],
    )
    @test mapped_weights[QN("Nf", -6, -1)] ≈ 0.3
    @test mapped_weights[QN("Nf", 0, -1)] ≈ 0.3
    @test mapped_weights[QN("Nf", 6, -1)] ≈ 0.4

    transfer = InfiniteCylinderDMRG._neutral_transfer_result(
        cfg,
        ComplexF64[1.0, 0.5],
        [1e-12, 2e-12],
        2,
        residual_tolerance=1e-10,
    )
    @test transfer.valid
    @test transfer.converged == 2
    @test transfer.eigenvalues == ComplexF64[1.0, 0.5]
    @test transfer.phases == [0.0, 0.0]
    @test transfer.residual_norms == [1e-12, 2e-12]
    @test transfer.ratio == 0.5
    @test transfer.xi_cell ≈ -1 / log(0.5)
    @test transfer.xi_x ≈ cfg.x_period * transfer.xi_cell

    invalid_transfer_cases = (
        (ComplexF64[1.0], [1e-12], 1, "fewer than two"),
        (ComplexF64[1.0, 0.5], [1e-4, 1e-12], 2, "residual"),
        (ComplexF64[1.0, 0.0], [1e-12, 1e-12], 2, "positive"),
        (ComplexF64[1.0, 1.0], [1e-12, 1e-12], 2, "strictly below one"),
        (ComplexF64[0.5, 1.0], [1e-12, 1e-12], 2, "strictly below one"),
        (ComplexF64[NaN, 0.5], [1e-12, 1e-12], 2, "finite"),
        (ComplexF64[1.0, 0.5], [NaN, 1e-12], 2, "finite"),
    )
    for (values, residuals, converged, message) in invalid_transfer_cases
        invalid = InfiniteCylinderDMRG._neutral_transfer_result(
            cfg,
            values,
            residuals,
            converged;
            residual_tolerance=1e-10,
        )
        @test !invalid.valid
        @test occursin(message, invalid.reason)
        @test isnan(invalid.xi_cell)
        @test isnan(invalid.xi_x)
    end
    @test_throws ArgumentError InfiniteCylinderDMRG._neutral_transfer_result(
        cfg,
        ComplexF64[1.0, 0.5],
        [1e-12, 1e-12],
        3;
        residual_tolerance=1e-10,
    )

    product_transfer = neutral_transfer_data(psi, cfg; neigs=2, tol=1e-10)
    @test product_transfer isa NeutralTransferData
    if product_transfer.valid
        @test product_transfer.converged >= 2
        @test all(<=(1e-10), product_transfer.residual_norms[1:2])
        @test 0 < product_transfer.ratio < 1
        @test isfinite(product_transfer.xi_cell)
    else
        @test !isempty(product_transfer.reason)
        @test isnan(product_transfer.xi_cell)
    end
    @test_throws ArgumentError neutral_transfer_data(psi, cfg; neigs=1)
    @test_throws ArgumentError neutral_transfer_data(psi, cfg; tol=0.0)

    seeded_first = neutral_transfer_data(
        psi,
        cfg;
        neigs=2,
        tol=1e-10,
        rng=Random.Xoshiro(2026),
    )
    seeded_second = neutral_transfer_data(
        psi,
        cfg;
        neigs=2,
        tol=1e-10,
        rng=Random.Xoshiro(2026),
    )
    @test seeded_first.eigenvalues == seeded_second.eigenvalues
    @test seeded_first.residual_norms == seeded_second.residual_norms
end
