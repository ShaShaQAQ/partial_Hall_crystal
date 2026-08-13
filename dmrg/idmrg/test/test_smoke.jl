using Test
using TOML
using HDF5
using InfiniteCylinderDMRG

const SMOKE_PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SMOKE_OUTPUT_ROOT = joinpath(@__DIR__, "..", "smoke_output")
const SMOKE_OUTPUT_FILENAMES = Set([
    "summary.toml",
    "convergence.tsv",
    "expansion.tsv",
    "density.tsv",
    "entanglement_spectrum.tsv",
    "schmidt_sectors.tsv",
    "transfer_spectrum.tsv",
    "state.h5",
])

function smoke_arguments(Ly::Int, x_period::Int, tolerance::Float64, label::String)
    return [
        "--Ly=$Ly",
        "--x_period=$x_period",
        "--filling_num=1",
        "--filling_den=3",
        "--t1=1",
        "--t3=.2",
        "--V1=1",
        "--V2=0",
        "--V3=0",
        "--phi_y=0",
        "--maxdim=1,4,8",
        "--cutoff=1e-8",
        "--vumps_tol=$tolerance",
        "--maxiter=50",
        "--transfer_tol=1e-8",
        "--output=$(joinpath(SMOKE_OUTPUT_ROOT, label))",
    ]
end

function parse_tsv(path::AbstractString)
    lines = readlines(path)
    @test length(lines) >= 2
    header = split(first(lines), '\t')
    return [Dict(zip(header, split(line, '\t'))) for line in lines[2:end]]
end

function assert_smoke_contract(settings::SinglePointSettings, result)
    output = settings.output
    config = settings.config
    @test !settings.allow_nonconverged
    @test result.valid
    @test result.optimization.converged
    @test !occursin("maximum iterations", result.optimization.reason)
    @test !isempty(result.optimization.records)
    @test last(result.optimization.records).converged
    @test maximum(record.maxlinkdim for record in result.optimization.records) > 1
    @test maximum(link_dimensions(result.optimization.psi)) > 1

    energies = (
        result.energy.per_cell,
        result.energy.per_x,
        result.energy.per_unit_cell,
        result.energy.per_site,
    )
    @test all(isfinite, energies)
    @test result.transfer.valid
    @test result.transfer.converged >= 2
    @test all(
        residual -> isfinite(residual) && residual <= settings.transfer_tol,
        result.transfer.residual_norms[1:2],
    )
    @test isfinite(result.transfer.xi_cell) && result.transfer.xi_cell > 0
    @test isfinite(result.transfer.xi_x) && result.transfer.xi_x > 0

    @test length(result.entanglements) == config.x_period
    for data in result.entanglements
        @test !isempty(data.levels)
        @test !isempty(data.sectors)
        @test isapprox(
            sum(level.probability for level in data.levels),
            1.0;
            atol=1e-10,
            rtol=0,
        )
        @test isapprox(
            sum(sector.weight for sector in data.sectors),
            1.0;
            atol=1e-10,
            rtol=0,
        )
    end

    @test isdir(output)
    @test SMOKE_OUTPUT_FILENAMES ⊆ Set(readdir(output))
    convergence = parse_tsv(joinpath(output, "convergence.tsv"))
    @test parse(Int, convergence[end]["maxlinkdim"]) > 1
    @test convergence[end]["converged"] == "true"
    expansions = parse_tsv(joinpath(output, "expansion.tsv"))
    @test !isempty(expansions)
    @test all(row -> row["progressed"] == "true", expansions)
    @test all(row -> !isempty(row["before"]) && !isempty(row["after"]), expansions)
    transfer_rows = parse_tsv(joinpath(output, "transfer_spectrum.tsv"))
    @test count(row -> row["converged"] == "true", transfer_rows) >= 2
    @test count(row -> row["valid"] == "true", transfer_rows) >= 2
    @test all(
        isfile(joinpath(output, filename)) && filesize(joinpath(output, filename)) > 0
        for filename in SMOKE_OUTPUT_FILENAMES
    )

    summary = TOML.parsefile(joinpath(output, "summary.toml"))
    @test summary["format"] == "infinite_cylinder_vumps_text_v1"
    @test summary["algorithm"] == "VUMPS"
    @test summary["valid"] === true
    @test summary["converged"] === true
    @test summary["configuration"] == Dict(
        "Ly" => config.Ly,
        "x_period" => config.x_period,
        "filling_num" => config.filling_num,
        "filling_den" => config.filling_den,
        "phi_y" => config.phi_y,
        "sites_per_cell" => sites_per_cell(config),
        "particles_per_cell" => particles_per_cell(config),
        "charge_scale" => charge_scale(config),
        "signature" => InfiniteCylinderDMRG._configuration_signature(config),
    )
    @test summary["dependencies"]["itensor_infinite_mps_commit"] ==
        InfiniteCylinderDMRG.ITENSOR_INFINITE_MPS_COMMIT
    @test summary["optimization"]["maxdim_schedule"] == [1, 4, 8]
    @test summary["optimization"]["vumps_tol"] == settings.vumps_tol
    @test summary["observables"]["transfer_valid"] === true
    @test summary["observables"]["correlation_length_neutral_cell"] > 0
    @test summary["observables"]["correlation_length_neutral_x"] > 0

    checkpoint = settings.checkpoint
    h5open(checkpoint, "r") do file
        metadata = attributes(file)
        @test read(metadata["format"]) == "infinite_cylinder_vumps_v1"
        @test read(metadata["backend_commit"]) ==
            InfiniteCylinderDMRG.ITENSOR_INFINITE_MPS_COMMIT
        @test read(metadata["configuration_signature"]) ==
            InfiniteCylinderDMRG._configuration_signature(config)
        @test read(metadata["charge_scale"]) == charge_scale(config)
    end
    reloaded = load_checkpoint(checkpoint, config)
    reloaded_density = density_data(reloaded, config)
    @test all(
        isapprox(before.density, after.density; rtol=1e-10, atol=1e-12)
        for (before, after) in zip(result.densities, reloaded_density)
    )
    for before in result.entanglements
        after = entanglement_data(reloaded, config; cut_x=before.cut_x)
        @test isapprox(before.entropy, after.entropy; rtol=1e-10, atol=1e-12)
        @test isapprox(
            before.raw_schmidt_polarization,
            after.raw_schmidt_polarization;
            rtol=1e-10,
            atol=1e-12,
        )
    end
    return nothing
end

@testset "PBS smoke workflow contract" begin
    pbs_path = joinpath(@__DIR__, "..", "jobs", "run_smoke.pbs")
    @test isfile(pbs_path)
    if isfile(pbs_path)
        script = read(pbs_path, String)
        @test occursin(r"#PBS -l nodes=1:ppn=24", script)
        @test occursin(r"#PBS -o .*dmrg/idmrg/logs/", script)
        @test occursin("--project=dmrg/idmrg", script)
        @test occursin("OPENBLAS_NUM_THREADS=1", script)
        @test occursin("MKL_NUM_THREADS=1", script)
        @test occursin("OMP_NUM_THREADS=1", script)
        @test occursin("--Ly=2", script)
        @test occursin("--x_period=3", script)
        @test occursin("--vumps_tol=1e-5", script)
        @test occursin("smoke_output/Ly2_xp3", script)
        @test occursin("--Ly=6", script)
        @test occursin("--x_period=1", script)
        @test occursin("--vumps_tol=1e-4", script)
        @test occursin("smoke_output/Ly6_xp1", script)
        @test occursin(r"Ly2_xp3[\s\S]*&&[\s\S]*--Ly=6", script)
        @test !occursin("--allow_nonconverged", script)
    end

    ignore_path = joinpath(SMOKE_PROJECT_ROOT, ".gitignore")
    @test "dmrg/idmrg/smoke_output/" in readlines(ignore_path)
end

@testset "opt-in VUMPS smoke contract" begin
    if get(ENV, "IDMRG_RUN_SMOKE", "false") == "true"
        for (Ly, x_period, tolerance, label) in (
            (2, 3, 1e-5, "Ly2_xp3"),
            (6, 1, 1e-4, "Ly6_xp1"),
        )
            @testset "$label" begin
                settings = parse_single_point_args(
                    smoke_arguments(Ly, x_period, tolerance, label)
                )
                result = run_single_point(settings)
                assert_smoke_contract(settings, result)
            end
        end
    else
        @test true # Expensive smoke runs are intentionally opt-in.
    end
end
