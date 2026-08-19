using Test
using ITensors
using ITensorMPS
using ITensorInfiniteMPS
using InfiniteCylinderDMRG
using LinearAlgebra: Diagonal
using Random

saved_canonical_test_default_rng = copy(Random.default_rng())
saved_canonical_test_index_rng = copy(ITensors.index_id_rng())
try
    @testset "VUMPS canonicalization owns and restores tensor RNG state" begin
        config = InfiniteCylinderConfig(; Ly=2, x_period=3)
        _, _, psi = initial_infinite_mps(config)
        inconsistent = copy(psi)
        inconsistent.C[1] *= -1

        Random.seed!(Random.default_rng(), 0x51a7e)
        Random.seed!(ITensors.index_id_rng(), 0x1d51)
        default_before = repr(copy(Random.default_rng()))
        index_before = repr(copy(ITensors.index_id_rng()))
        first_canonical = InfiniteCylinderDMRG._canonicalize_vumps_state(
            inconsistent; rng_seed=0x1234
        )
        @test repr(copy(Random.default_rng())) == default_before
        @test repr(copy(ITensors.index_id_rng())) == index_before

        Random.seed!(Random.default_rng(), 0x9f31)
        Random.seed!(ITensors.index_id_rng(), 0xa471)
        second_canonical = InfiniteCylinderDMRG._canonicalize_vumps_state(
            inconsistent; rng_seed=0x1234
        )
        @test id.(linkinds(only, second_canonical.AR)) ==
            id.(linkinds(only, first_canonical.AR))
        @test isnothing(
            InfiniteCylinderDMRG._validate_checkpoint_state(
                second_canonical, config
            ),
        )

        for invalid_seed in (true, -1, big(typemax(UInt64)) + 1)
            @test_throws ArgumentError InfiniteCylinderDMRG._canonicalize_vumps_state(
                inconsistent; rng_seed=invalid_seed
            )
        end
    end
finally
    copy!(Random.default_rng(), saved_canonical_test_default_rng)
    copy!(ITensors.index_id_rng(), saved_canonical_test_index_rng)
end

default_rng_restored = repr(copy(Random.default_rng())) ==
    repr(saved_canonical_test_default_rng)
index_rng_restored = repr(copy(ITensors.index_id_rng())) ==
    repr(saved_canonical_test_index_rng)
@testset "canonicalization test preserves its caller RNG state" begin
    @test default_rng_restored
    @test index_rng_restored
end

@testset "deterministic Arnoldi fixed point has an audited eigen residual" begin
    config = InfiniteCylinderConfig(; Ly=2, x_period=3)
    _, _, psi = initial_infinite_mps(config)
    fixed_point, eigenvalue =
        InfiniteCylinderDMRG._arnoldi_transfer_fixed_point(psi.AL)
    transfer = ITensorInfiniteMPS.TransferMatrix(psi.AL)
    applied = transfer(fixed_point)
    residual = norm(applied - eigenvalue * fixed_point) /
        max(norm(applied), abs(eigenvalue) * norm(fixed_point), 1.0)
    @test residual <= 1e-12
    adjoint_fixed_point = swapinds(
        dag(fixed_point), reverse(Pair(inds(fixed_point)...))
    )
    @test norm(fixed_point - adjoint_fixed_point) / norm(fixed_point) <= 1e-10
end

@testset "positive transfer square root audits the fixed-point spectrum" begin
    @test isdefined(InfiniteCylinderDMRG, :_positive_transfer_sqrt)
    if isdefined(InfiniteCylinderDMRG, :_positive_transfer_sqrt)
        positive_transfer_sqrt =
            InfiniteCylinderDMRG._positive_transfer_sqrt

        @testset "dense tolerance-scale negative mode is clamped" begin
            link = Index(3, "Transfer,Dense")
            row_index = link'
            column_index = dag(link)
            fixed_point = ITensor(
                Diagonal([-5.0e-13, 4.0, 9.0]),
                row_index,
                column_index,
            )

            positive_sqrt = positive_transfer_sqrt(
                fixed_point; tol=1.0e-12
            )
            sqrt_matrix = Array(
                positive_sqrt, row_index, column_index
            )

            @test inds(positive_sqrt) == inds(fixed_point)
            @test all(isfinite, sqrt_matrix)
            @test sqrt_matrix[1, 1] == 0.0
            @test sqrt_matrix ≈ Diagonal([0.0, 2.0, 3.0]) atol=1.0e-13 rtol=0
        end

        @testset "dense material negative mode fails closed with context" begin
            link = Index(2, "Transfer,Indefinite")
            fixed_point = ITensor(
                Diagonal([-1.0e-4, 1.0]), link', dag(link)
            )
            caught = try
                positive_transfer_sqrt(fixed_point; tol=1.0e-12)
                nothing
            catch error
                error
            end

            @test caught isa Exception
            if caught isa Exception
                message = lowercase(sprint(showerror, caught))
                @test occursin("canonical transfer", message)
                @test occursin("negative", message)
                @test occursin("eigenvalue", message)
            end
        end

        @testset "QN blocks preserve their index metadata and neutral flux" begin
            link = Index(
                [QN("Transfer", -1) => 2, QN("Transfer", 1) => 2];
                tags="Transfer,QNFixedPoint",
                dir=ITensors.Out,
            )
            row_index = link'
            column_index = dag(link)
            fixed_point_matrix = [
                -5.0e-13 0.0 0.0 0.0
                0.0 4.0 0.0 0.0
                0.0 0.0 5.0 4.0
                0.0 0.0 4.0 5.0
            ]
            clamped_matrix = [
                0.0 0.0 0.0 0.0
                0.0 4.0 0.0 0.0
                0.0 0.0 5.0 4.0
                0.0 0.0 4.0 5.0
            ]
            fixed_point = ITensor(
                fixed_point_matrix, row_index, column_index
            )
            input_indices = inds(fixed_point)
            input_spaces = [
                [
                    qn(index, block) => blockdim(index, block) for
                    block in 1:nblocks(index)
                ] for index in input_indices
            ]

            positive_sqrt = positive_transfer_sqrt(
                fixed_point; tol=1.0e-12
            )
            output_indices = inds(positive_sqrt)
            sqrt_matrix = Array(
                positive_sqrt, row_index, column_index
            )

            @test output_indices == input_indices
            @test id.(output_indices) == id.(input_indices)
            @test tags.(output_indices) == tags.(input_indices)
            @test plev.(output_indices) == plev.(input_indices)
            @test dir.(output_indices) == dir.(input_indices)
            @test [
                [
                    qn(index, block) => blockdim(index, block) for
                    block in 1:nblocks(index)
                ] for index in output_indices
            ] == input_spaces
            @test all(index -> nblocks(index) == 2, output_indices)
            @test flux(positive_sqrt) == QN()
            @test all(isfinite, sqrt_matrix)
            @test sqrt_matrix[1, 1] == 0.0
            @test sqrt_matrix * sqrt_matrix ≈ clamped_matrix atol=1.0e-12 rtol=0
        end
    end
end

@testset "canonical right links always receive fresh IDs" begin
    initspin(_) = "↑"
    dense_sites = infsiteinds("S=1/2", 1; initstate=initspin)
    dense_state = InfMPS(dense_sites, initspin)
    dense_canonical = ITensorMPS.orthogonalize(dense_state.AL, :)
    dense_before = id.(linkinds(only, dense_canonical.AR))
    dense_normalized =
        InfiniteCylinderDMRG._normalize_right_link_convention(dense_canonical)
    @test id.(linkinds(only, dense_normalized.AL)) ==
        id.(linkinds(only, dense_canonical.AL))
    @test id.(linkinds(only, dense_normalized.AR)) != dense_before
    @test isempty(
        intersect(
            Set(id.(linkinds(only, dense_normalized.AL))),
            Set(id.(linkinds(only, dense_normalized.AR))),
        ),
    )

    config = InfiniteCylinderConfig(; Ly=2, x_period=3)
    _, _, qn_canonical = initial_infinite_mps(config)
    qn_before = id.(linkinds(only, qn_canonical.AR))
    qn_normalized =
        InfiniteCylinderDMRG._normalize_right_link_convention(qn_canonical)
    @test id.(linkinds(only, qn_normalized.AL)) ==
        id.(linkinds(only, qn_canonical.AL))
    @test id.(linkinds(only, qn_normalized.AR)) != qn_before
    @test isempty(
        intersect(
            Set(id.(linkinds(only, qn_normalized.AL))),
            Set(id.(linkinds(only, qn_normalized.AR))),
        ),
    )
    @test isnothing(
        InfiniteCylinderDMRG._validate_checkpoint_state(
            qn_normalized, config
        ),
    )
end

@testset "workflow canonical RNG streams are stable and candidate-specific" begin
    @test InfiniteCylinderDMRG._single_point_canonical_seed(17) ==
        InfiniteCylinderDMRG._single_point_canonical_seed(17)
    @test InfiniteCylinderDMRG._single_point_canonical_seed(17) !=
        InfiniteCylinderDMRG._single_point_canonical_seed(18)

    @test isdefined(InfiniteCylinderDMRG, :_flux_candidate_seed)
    if isdefined(InfiniteCylinderDMRG, :_flux_candidate_seed)
        flux_candidate_seed = InfiniteCylinderDMRG._flux_candidate_seed
        flux_reference = flux_candidate_seed(17, 2, "cold_01")
        @test flux_reference isa Int
        @test 0 <= flux_reference <= typemax(Int)
        @test flux_reference == flux_candidate_seed(17, 2, "cold_01")
        @test flux_reference != flux_candidate_seed(18, 2, "cold_01")
        @test flux_reference != flux_candidate_seed(17, 3, "cold_01")
        @test flux_reference != flux_candidate_seed(17, 2, "cold_02")
    end

    reference = InfiniteCylinderDMRG._fig2_candidate_seed(32, 1, "cold_a")
    @test reference == InfiniteCylinderDMRG._fig2_candidate_seed(
        32, 1, "cold_a"
    )
    @test reference != InfiniteCylinderDMRG._fig2_candidate_seed(
        64, 1, "cold_a"
    )
    @test reference != InfiniteCylinderDMRG._fig2_candidate_seed(
        32, 2, "cold_a"
    )
    @test reference != InfiniteCylinderDMRG._fig2_candidate_seed(
        32, 1, "cold_b"
    )
    @test 0 <= reference <= typemax(Int)
end
