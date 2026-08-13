using Test
using InfiniteCylinderDMRG
using ITensors
using ITensorMPS
using ITensorInfiniteMPS
using KrylovKit
using HDF5
using Random
using TOML

function argument_error_message(f)
    error = try
        f()
        nothing
    catch exception
        exception
    end
    @test error isa ArgumentError
    return sprint(showerror, error)
end

@testset "single-point CLI parsing" begin
    args = [
        "--Ly=6",
        "--x_period=1",
        "--filling_num=1",
        "--filling_den=3",
        "--phi_y=0",
        "--maxdim=1,4,8",
        "--t1=1.25",
        "--t3=0.3",
        "--V1=2.0",
        "--V2=0.1",
        "--V3=-0.2",
        "--cutoff=1e-9",
        "--vumps_tol=2e-5",
        "--energy_tol=3e-5",
        "--mismatch_tol=4e-5",
        "--imaginary_tol=5e-12",
        "--transfer_tol=6e-9",
        "--transfer_neigs=3",
        "--maxiter=17",
        "--stability=4",
        "--threads=2",
        "--output=out/single",
        "--checkpoint=out/checkpoint.h5",
        "--load=in/seed.h5",
        "--occupied_sites=2",
        "--seed=1234",
        "--allow_nonconverged=true",
    ]
    settings = parse_single_point_args(args)

    @test settings isa SinglePointSettings
    @test settings.config == InfiniteCylinderConfig(6, 1, 1, 3, 0.0)
    @test settings.model == CylinderModelParams(1.25, 0.3, 2.0, 0.1, -0.2)
    @test settings.maxdim_schedule == [1, 4, 8]
    @test settings.cutoff == 1e-9
    @test settings.vumps_tol == 2e-5
    @test settings.energy_tol == 3e-5
    @test settings.energy_mismatch_tol == 4e-5
    @test settings.imaginary_tol == 5e-12
    @test settings.transfer_tol == 6e-9
    @test settings.transfer_neigs == 3
    @test settings.max_iterations == 17
    @test settings.stable_iterations == 4
    @test settings.threads == 2
    @test settings.output == "out/single"
    @test settings.checkpoint == "out/checkpoint.h5"
    @test settings.load == "in/seed.h5"
    @test settings.occupied_sites == [2]
    @test settings.seed == 1234
    @test settings.allow_nonconverged
    @test InfiniteCylinderDMRG._optimization_metadata(settings).model == (
        t1=1.25,
        t3=0.3,
        V1=2.0,
        V2=0.1,
        V3=-0.2,
    )

    defaults = parse_single_point_args([
        "--Ly=2",
        "--x_period=3",
        "--filling_num=1",
        "--filling_den=3",
        "--phi_y=0",
        "--maxdim=1,4,8",
        "--output=out/defaults",
    ])
    @test defaults.model == CylinderModelParams(1.0, 0.2, 1.0, 0.0, 0.0)
    @test defaults.cutoff == 1e-8
    @test defaults.vumps_tol == 1e-5
    @test defaults.energy_tol == 1e-4
    @test defaults.energy_mismatch_tol == 1e-4
    @test defaults.imaginary_tol == 1e-12
    @test defaults.transfer_tol == 1e-8
    @test defaults.transfer_neigs == 4
    @test defaults.max_iterations == 50
    @test defaults.stable_iterations == 2
    @test defaults.threads == 1
    @test defaults.checkpoint == joinpath("out/defaults", "state.h5")
    @test isnothing(defaults.load)
    @test defaults.occupied_sites == default_occupied_sites(defaults.config)
    @test defaults.seed == 0
    @test !defaults.allow_nonconverged
end

@testset "single-point CLI rejects malformed input" begin
    valid = [
        "--Ly=6",
        "--x_period=1",
        "--filling_num=1",
        "--filling_den=3",
        "--phi_y=0",
        "--maxdim=1,4,8",
        "--output=out",
    ]
    bad_cases = [
        (["Ly=6"; valid[2:end]], "--key=value"),
        ([valid; "--mystery=1"], "unknown option --mystery"),
        ([valid; "--Ly=8"], "duplicate option --Ly"),
        ([valid[1:5]; "--maxdim=1,bad"; valid[7:end]], "--maxdim"),
        ([valid; "--allow_nonconverged=yes"], "--allow_nonconverged"),
        ([valid[1:4]; "--phi_y=NaN"; valid[6:end]], "--phi_y must be finite"),
        ([valid; "--cutoff=Inf"], "--cutoff must be finite"),
        ([valid[1:5]; "--maxdim=1,1"; valid[7:end]], "strictly increasing"),
        ([valid[1:5]; "--maxdim="; valid[7:end]], "must not be empty"),
        ([valid; "--seed=-1"], "--seed must be nonnegative"),
        ([valid; "--seed=bad"], "--seed must be an integer"),
        ([valid; "--occupied_sites="], "--occupied_sites must not be empty"),
        ([valid; "--occupied_sites=1,2"], "particle count"),
        ([valid; "--occupied_sites=7"], "outside"),
        ([valid; "--occupied_sites=1,1"], "particle count"),
    ]
    for (arguments, message) in bad_cases
        @test occursin(message, argument_error_message(() -> parse_single_point_args(arguments)))
    end

    for required in ("Ly", "x_period", "filling_num", "filling_den", "phi_y", "maxdim", "output")
        arguments = filter(argument -> !startswith(argument, "--$required="), valid)
        @test occursin(
            "missing required option --$required",
            argument_error_message(() -> parse_single_point_args(arguments)),
        )
    end

    invalid_required = [
        ("Ly", "odd", "--Ly must be an integer"),
        ("x_period", "0", "x_period must be positive"),
        ("filling_num", "0", "filling_num must be positive"),
        ("filling_den", "0", "filling_den must be positive"),
        ("output", "", "--output must not be empty"),
    ]
    for (key, value, message) in invalid_required
        arguments = replace.(valid, Regex("^--$key=.*") => "--$key=$value")
        @test occursin(message, argument_error_message(() -> parse_single_point_args(arguments)))
    end

    @test occursin(
        "stability",
        argument_error_message(() -> parse_single_point_args([
            valid;
            "--maxiter=2";
            "--stability=3";
        ])),
    )
    @test occursin(
        "energy_tol",
        argument_error_message(() -> parse_single_point_args([
            valid;
            "--vumps_tol=1e308";
        ])),
    )
end

@testset "checkpoint optimization restart is rejected before callbacks" begin
    settings = parse_single_point_args([
        "--Ly=6",
        "--x_period=1",
        "--filling_num=1",
        "--filling_den=3",
        "--phi_y=0",
        "--maxdim=1",
        "--output=unused",
        "--load=checkpoint.h5",
    ])
    calls = Ref(0)
    callback(args...) = (calls[] += 1; error("callback must not run"))
    operations = SinglePointOperations(
        configure_threads=callback,
        initialize=callback,
        load_state=callback,
        state_sites=callback,
        build_hamiltonian=callback,
        optimize=callback,
        energy=callback,
        density=callback,
        entanglement=callback,
        transfer=callback,
        fidelity=callback,
        write_outputs=callback,
        save_state=callback,
        cold_state=callback,
    )
    error = try
        run_single_point(settings; operations)
        nothing
    catch exception
        exception
    end
    @test error isa WorkflowValidationError
    message = sprint(showerror, error)
    @test occursin("765f2777703bc1138b009adbed1b97bde1973402", message)
    @test occursin("unsupported", lowercase(message))
    @test occursin("audit", lowercase(message))
    @test occursin("observable-only", lowercase(message))
    @test occursin("fresh deterministic", lowercase(message))
    @test occursin("cold", lowercase(message))
    @test calls[] == 0
end

@testset "flux-scan checkpoint restart is rejected before callbacks" begin
    settings = parse_flux_scan_args([
        "--Ly=6",
        "--x_period=1",
        "--filling_num=1",
        "--filling_den=3",
        "--phi_start=0",
        "--phi_stop=1",
        "--phi_steps=2",
        "--maxdim=1",
        "--output=unused",
        "--load=checkpoint.h5",
    ])
    calls = Ref(0)
    callback(args...) = (calls[] += 1; error("callback must not run"))
    operations = SinglePointOperations(
        configure_threads=callback,
        initialize=callback,
        load_state=callback,
        state_sites=callback,
        build_hamiltonian=callback,
        optimize=callback,
        energy=callback,
        density=callback,
        entanglement=callback,
        transfer=callback,
        fidelity=callback,
        write_outputs=callback,
        save_state=callback,
        cold_state=callback,
    )
    error = try
        run_flux_scan(settings; operations)
        nothing
    catch exception
        exception
    end
    @test error isa WorkflowValidationError
    message = sprint(showerror, error)
    @test occursin("765f2777703bc1138b009adbed1b97bde1973402", message)
    @test occursin("unsupported", lowercase(message))
    @test occursin("audit", lowercase(message))
    @test occursin("observable-only", lowercase(message))
    @test occursin("fresh deterministic", lowercase(message))
    @test occursin("cold", lowercase(message))
    @test calls[] == 0
end

@testset "safe CLI threading" begin
    configured = configure_cli_threads(1)
    @test configured.requested == 1
    @test configured.blas == 1
    @test configured.strided == 1
    @test !configured.blocksparse
    @test_throws ArgumentError configure_cli_threads(0)
    @test_throws ArgumentError configure_cli_threads(true)
end

@testset "deterministic runtime dependency metadata" begin
    dependencies = InfiniteCylinderDMRG._runtime_dependencies()
    @test dependencies == (
        julia_version=string(VERSION),
        itensor_infinite_mps_version=string(Base.pkgversion(ITensorInfiniteMPS)),
        itensor_mps_version=string(Base.pkgversion(ITensorMPS)),
        itensors_version=string(Base.pkgversion(ITensors)),
        krylovkit_version=string(Base.pkgversion(KrylovKit)),
        hdf5_version=string(Base.pkgversion(HDF5)),
    )

    mktempdir() do directory
        settings = parse_single_point_args([
            "--Ly=6",
            "--x_period=1",
            "--filling_num=1",
            "--filling_den=3",
            "--phi_y=0",
            "--maxdim=1",
            "--output=$directory",
        ])
        cfg = settings.config
        _, _, psi = initial_infinite_mps(cfg)
        result = VUMPSResult(psi, VUMPSRecord[], true, "converged")
        InfiniteCylinderDMRG._default_write_outputs(
            directory,
            cfg,
            result,
            normalize_energy(cfg, -1.0),
            density_data(psi, cfg),
            [entanglement_data(psi, cfg; cut_x=1)],
            InfiniteCylinderDMRG._neutral_transfer_result(
                cfg,
                ComplexF64[1.0, 0.5],
                [1e-12, 1e-12],
                2;
                residual_tolerance=1e-10,
            ),
            settings,
        )
        summary = TOML.parsefile(joinpath(directory, "summary.toml"))
        actual = summary["dependencies"]
        @test actual["julia_version"] == string(VERSION)
        @test actual["itensor_infinite_mps_version"] ==
            string(Base.pkgversion(ITensorInfiniteMPS))
        @test actual["itensor_mps_version"] == string(Base.pkgversion(ITensorMPS))
        @test actual["itensors_version"] == string(Base.pkgversion(ITensors))
        @test actual["krylovkit_version"] == string(Base.pkgversion(KrylovKit))
        @test actual["hdf5_version"] == string(Base.pkgversion(HDF5))
        @test actual["itensor_infinite_mps_commit"] ==
            InfiniteCylinderDMRG.ITENSOR_INFINITE_MPS_COMMIT
    end
end

@testset "guarded bin entry points" begin
    vumps_bin = joinpath(@__DIR__, "..", "bin", "run_vumps.jl")
    flux_bin = joinpath(@__DIR__, "..", "bin", "run_flux_scan.jl")
    @test include(vumps_bin) === nothing
    @test include(flux_bin) === nothing
end

@testset "raw flux-scan CLI and inclusive grid" begin
    arguments = [
        "--Ly=6",
        "--x_period=1",
        "--filling_num=1",
        "--filling_den=3",
        "--phi_start=-1",
        "--phi_stop=1",
        "--phi_steps=5",
        "--maxdim=1,4,8",
        "--output=out/scan",
        "--branch_mode=adiabatic",
        "--cold_patterns=1;2",
        "--fidelity_drop_tol=0.002",
        "--charge_jump_tol=0.2",
        "--sector_tol=0.3",
        "--spectrum_tol=0.4",
        "--spectrum_levels=12",
    ]
    settings = parse_flux_scan_args(arguments)
    @test settings isa FluxScanSettings
    @test settings.point.config == InfiniteCylinderConfig(6, 1, 1, 3, -1.0)
    @test settings.point.output == "out/scan"
    @test settings.phi_stop == 1.0
    @test settings.phi_steps == 5
    @test settings.branch_mode === :adiabatic
    @test settings.cold_patterns == [[1], [2]]
    @test settings.fidelity_drop_tol == 0.002
    @test settings.charge_jump_tol == 0.2
    @test settings.sector_tol == 0.3
    @test settings.spectrum_tol == 0.4
    @test settings.spectrum_levels == 12
    @test flux_grid(settings) == [-1.0, -0.5, 0.0, 0.5, 1.0]
    @test first(flux_grid(1.0, -1.0, 5)) == 1.0
    @test last(flux_grid(1.0, -1.0, 5)) == -1.0
    @test flux_grid(0.25, 0.25, 1) == [0.25]

    defaults = parse_flux_scan_args(filter(
        argument -> !startswith(argument, "--branch_mode=") &&
                    !startswith(argument, "--cold_patterns=") &&
                    !startswith(argument, "--fidelity_drop_tol=") &&
                    !startswith(argument, "--charge_jump_tol=") &&
                    !startswith(argument, "--sector_tol=") &&
                    !startswith(argument, "--spectrum_tol=") &&
                    !startswith(argument, "--spectrum_levels="),
        arguments,
    ))
    @test defaults.branch_mode === :ground
    @test isempty(defaults.cold_patterns)
    @test defaults.fidelity_drop_tol == 1e-3
    @test defaults.charge_jump_tol == 0.1
    @test defaults.sector_tol == 0.1
    @test defaults.spectrum_tol == 0.5
    @test defaults.spectrum_levels == 20
    metadata = InfiniteCylinderDMRG._scan_metadata(defaults)
    @test metadata["branch"]["mode"] == "ground"
    @test metadata["branch"]["cut_x"] == 1
    @test metadata["branch"]["fidelity_drop_tol"] == 1e-3
    @test metadata["branch"]["charge_jump_tol"] == 0.1
    @test metadata["branch"]["sector_tol"] == 0.1
    @test metadata["branch"]["spectrum_tol"] == 0.5
    @test metadata["branch"]["spectrum_levels"] == 20

    bad_cases = [
        ([arguments; "--phi_start=0"], "duplicate option --phi_start"),
        (replace.(arguments, "--branch_mode=adiabatic" => "--branch_mode=follow"), "--branch_mode"),
        (replace.(arguments, "--phi_steps=5" => "--phi_steps=1"), "equal start and stop"),
        (replace.(arguments, "--phi_steps=5" => "--phi_steps=0"), "--phi_steps must be positive"),
        (replace.(arguments, "--cold_patterns=1;2" => "--cold_patterns=1,2"), "particle count"),
        (replace.(arguments, "--cold_patterns=1;2" => "--cold_patterns=1;7"), "outside"),
        (replace.(arguments, "--cold_patterns=1;2" => "--cold_patterns=1;1,1"), "particle count"),
    ]
    for (bad, message) in bad_cases
        @test occursin(message, argument_error_message(() -> parse_flux_scan_args(bad)))
    end
end

@testset "pure raw flux rows and branch events" begin
    levels_a = [
        (raw_charge=0, entanglement_energy=0.1, probability=0.6),
        (raw_charge=1, entanglement_energy=0.3, probability=0.4),
    ]
    levels_b = [
        (raw_charge=0, entanglement_energy=0.15, probability=0.5),
        (raw_charge=1, entanglement_energy=0.9, probability=0.5),
    ]
    previous = [
        FluxCandidateData("warm", -1.0, 0.20, Dict(0 => 0.6, 1 => 0.4), levels_a),
        FluxCandidateData("cold_01", -0.9, 0.10, Dict(0 => 0.7, 1 => 0.3), levels_a),
    ]
    current = [
        FluxCandidateData(
            "warm",
            -0.8,
            -0.25,
            Dict(0 => 0.2, 1 => 0.8),
            levels_b;
            fidelity_to_previous=0.7,
            fidelity_valid=true,
        ),
        FluxCandidateData(
            "cold_01",
            -1.1,
            0.15,
            Dict(0 => 0.65, 1 => 0.35),
            levels_a;
            fidelity_to_previous=0.95,
            fidelity_valid=true,
        ),
    ]

    @test select_flux_candidate_data(current; mode=:ground) == 2
    @test select_flux_candidate_data(current; mode=:adiabatic) == 2
    tied = [current[2], current[2]]
    @test select_flux_candidate_data(tied; mode=:ground) == 1
    @test select_flux_candidate_data(tied; mode=:adiabatic) == 1

    first_row = flux_scan_row(
        1,
        0.0,
        previous[1],
        1;
        reference_polarization=0.20,
        x_period=3,
        selection_mode=:ground,
        ground_energy=-1.0,
        ground_energy_valid=true,
    )
    next_row = flux_scan_row(
        2,
        0.5,
        current[1],
        1;
        reference_polarization=0.20,
        x_period=3,
        selection_mode=:adiabatic,
        ground_energy=-1.1,
        ground_energy_valid=true,
    )
    @test first_row.delta_raw_schmidt_polarization == 0.0
    @test next_row.raw_schmidt_polarization == -0.25
    @test next_row.delta_raw_schmidt_polarization == -0.45
    @test next_row.polarization_valid
    @test first_row.cut_x == 1
    @test first_row.selection_mode == "ground"
    @test next_row.selection_mode == "adiabatic"
    @test next_row.energy_delta_to_ground_raw ≈ 0.3
    @test next_row.energy_delta_to_ground_valid
    @test next_row.fidelity_density_x_raw ≈ -log(0.7) / 3
    @test next_row.fidelity_density_x_valid

    gauge_rows = InfiniteCylinderDMRG._same_flux_sector_gauge_rows(
        2,
        0.5,
        [
            current[1],
            FluxCandidateData(
                "cold_01",
                -1.1,
                0.1,
                Dict(2 => 0.2, 3 => 0.8),
                levels_a,
            ),
            FluxCandidateData(
                "cold_02",
                -1.0,
                0.2,
                Dict(0 => 0.2, 1 => 0.8),
                levels_a,
            ),
        ],
    )
    @test length(gauge_rows) == 1
    @test only(gauge_rows).first_label == "cold_01"
    @test only(gauge_rows).second_label == "cold_02"
    @test only(gauge_rows).sector_shift_raw == 2
    @test only(gauge_rows).sector_distance_raw ≈ 0.0
    @test only(gauge_rows).valid
    cut_mismatch = InfiniteCylinderDMRG._same_flux_sector_gauge_rows(
        2,
        0.5,
        [
            FluxCandidateData(
                "cold_01", -1.1, 0.1, Dict(0 => 1.0), levels_a; cut_x=1
            ),
            FluxCandidateData(
                "cold_02", -1.0, 0.2, Dict(0 => 1.0), levels_a; cut_x=2
            ),
        ],
    )
    @test !only(cut_mismatch).valid
    @test occursin("cuts do not match", only(cut_mismatch).reason)

    event = flux_branch_event(
        first_row,
        next_row,
        previous,
        current;
        fidelity_drop_tol=0.1,
        charge_jump_tol=0.4,
        sector_tol=0.08,
        spectrum_tol=0.5,
        spectrum_levels=20,
    )
    @test event.fidelity_raw == 0.7
    @test event.fidelity_valid
    @test event.charge_step_raw == -0.45
    @test event.charge_step_valid
    @test event.sector_distance_raw > 0
    @test event.sector_distance_valid
    @test event.sector_shift_raw == 0
    @test event.sector_shift_valid
    @test event.spectrum_distance_raw ≈ 0.6
    @test event.spectrum_distance_valid
    @test event.energy_crossing_raw
    @test event.energy_crossing_valid
    @test event.fidelity_drop_tol == 0.1
    @test event.charge_jump_tol == 0.4
    @test event.sector_tol == 0.08
    @test event.spectrum_tol == 0.5
    @test event.spectrum_levels == 20
    @test event.flag_fidelity
    @test event.flag_charge
    @test event.flag_sector
    @test event.flag_spectrum
    @test event.flag_energy_crossing
    @test event.flags_valid

    invalid = FluxCandidateData(
        "invalid",
        NaN,
        NaN,
        Dict{Int,Float64}(),
        NamedTuple[];
        converged=false,
        valid=false,
        energy_valid=false,
        polarization_valid=false,
        sector_valid=false,
        spectrum_valid=false,
    )
    @test_throws ArgumentError select_flux_candidate_data([invalid]; mode=:ground)
    @test_throws ArgumentError select_flux_candidate_data(current; mode=:other)

    raw_nonconverged = FluxCandidateData(
        "raw_nonconverged",
        -2.0,
        0.0,
        Dict(0 => 1.0),
        levels_a;
        converged=false,
        valid=false,
        fidelity_to_previous=0.99,
        fidelity_valid=true,
    )
    @test select_flux_candidate_data([current[1], raw_nonconverged]; mode=:ground) == 1
    @test select_flux_candidate_data([current[1], raw_nonconverged]; mode=:adiabatic) == 1

    ancillary_invalid = FluxCandidateData(
        "ancillary_invalid",
        -2.0,
        NaN,
        Dict{Int,Float64}(),
        NamedTuple[];
        converged=true,
        valid=false,
        energy_valid=true,
        polarization_valid=false,
        sector_valid=false,
        spectrum_valid=false,
        fidelity_to_previous=0.99,
        fidelity_valid=true,
        reason="neutral transfer failed",
    )
    @test select_flux_candidate_data([current[1], ancillary_invalid]; mode=:ground) == 2
    @test select_flux_candidate_data([current[1], ancillary_invalid]; mode=:adiabatic) == 2
    ancillary_row = flux_scan_row(
        2,
        0.5,
        FluxCandidateData(
            "ancillary_invalid",
            -2.0,
            0.25,
            Dict(0 => 1.0),
            levels_a;
            converged=true,
            valid=false,
            energy_valid=true,
            polarization_valid=true,
            sector_valid=true,
            spectrum_valid=true,
            fidelity_to_previous=0.9,
            fidelity_valid=true,
            reason="neutral transfer failed",
        ),
        1;
        reference_polarization=0.2,
        x_period=1,
        selection_mode=:adiabatic,
        ground_energy=-2.0,
        ground_energy_valid=true,
    )
    @test ancillary_row.polarization_valid
    @test ancillary_row.delta_raw_schmidt_polarization ≈ 0.05

    invalid_crossing = flux_branch_event(
        first_row,
        next_row,
        [previous[1], invalid],
        [current[1], invalid];
        fidelity_drop_tol=0.1,
        charge_jump_tol=0.4,
        sector_tol=0.08,
        spectrum_tol=0.5,
        spectrum_levels=20,
    )
    @test !invalid_crossing.energy_crossing_valid
    @test !invalid_crossing.flag_energy_crossing
end

function synthetic_transfer(cfg; valid=true)
    return InfiniteCylinderDMRG._neutral_transfer_result(
        cfg,
        valid ? ComplexF64[1.0, 0.5] : ComplexF64[1.0],
        valid ? [1e-12, 1e-12] : [1e-12],
        valid ? 2 : 1;
        residual_tolerance=1e-10,
    )
end

function fake_single_point_operations(events; converged=true, transfer_valid=true)
    return SinglePointOperations(
        configure_threads=threads -> push!(events, (:threads, threads)),
        initialize=(config, occupied_sites) -> begin
            push!(events, (:initialize, config.phi_y, copy(occupied_sites)))
            sites, _, psi = initial_infinite_mps(config; occupied_sites)
            (; sites, psi)
        end,
        load_state=(path, config) -> error("unexpected load"),
        state_sites=psi -> siteinds(only, psi.AL),
        build_hamiltonian=(config, model, sites) -> begin
            push!(events, (:build, config.phi_y, [ITensors.id(site) for site in sites]))
            (:hamiltonian, config.phi_y)
        end,
        optimize=(H, psi, settings) -> begin
            push!(events, (:optimize, H[2]))
            VUMPSResult(psi, VUMPSRecord[], converged, converged ? "converged" : "stopped")
        end,
        energy=(psi, H, config, settings) -> begin
            push!(events, (:energy, H[2], settings.imaginary_tol))
            normalize_energy(config, -1.0 - H[2])
        end,
        density=(psi, config) -> begin
            push!(events, (:density, config.phi_y))
            density_data(psi, config)
        end,
        entanglement=(psi, config, cut_x) -> begin
            push!(events, (:entanglement, config.phi_y, cut_x))
            entanglement_data(psi, config; cut_x)
        end,
        transfer=(psi, config, settings, rng) -> begin
            push!(events, (:transfer, config.phi_y, rand(rng, UInt64)))
            synthetic_transfer(config; valid=transfer_valid)
        end,
        fidelity=(previous, current, config, settings, rng) -> begin
            push!(events, (:fidelity, config.phi_y, rand(rng, UInt64)))
            InfiniteCylinderDMRG._mixed_transfer_result(
                0.9 + 0im,
                1.0 + 0im,
                1.0 + 0im,
                0.0,
                0.0,
                0.0,
                true;
                residual_tolerance=1e-10,
            )
        end,
        write_outputs=(directory, config, result, energy, density, entanglement, transfer, settings) -> begin
            push!(events, (:write, directory, result.converged, transfer.valid))
            mkpath(directory)
            write(joinpath(directory, "raw.txt"), "raw")
            directory
        end,
        save_state=(path, psi, config) -> begin
            push!(events, (:save, path, config.phi_y))
            mkpath(dirname(path))
            write(path, "checkpoint")
            path
        end,
    )
end

@testset "callable single-point workflow preserves raw invalid output" begin
    mktempdir() do directory
        arguments = [
            "--Ly=6",
            "--x_period=1",
            "--filling_num=1",
            "--filling_den=3",
            "--phi_y=0",
            "--maxdim=1",
            "--output=$(joinpath(directory, "point"))",
        ]
        settings = parse_single_point_args(arguments)
        events = Any[]
        result = run_single_point(settings; operations=fake_single_point_operations(events))
        @test result isa SinglePointRunResult
        @test result.valid
        @test result.optimization.converged
        @test result.transfer.valid
        @test (:initialize, 0.0, [1]) in events
        @test length(result.entanglements) == settings.config.x_period
        @test findfirst(event -> event[1] === :optimize, events) <
            findfirst(event -> event[1] === :energy, events) <
            findfirst(event -> event[1] === :transfer, events)
        @test isfile(joinpath(settings.output, "raw.txt"))
        @test isfile(settings.checkpoint)

        stopped_settings = parse_single_point_args([arguments; "--allow_nonconverged=false"])
        stopped_events = Any[]
        error = try
            run_single_point(
                stopped_settings;
                operations=fake_single_point_operations(
                    stopped_events;
                    converged=false,
                    transfer_valid=false,
                ),
            )
            nothing
        catch exception
            exception
        end
        @test error isa WorkflowValidationError
        @test occursin("nonconverged", sprint(showerror, error))
        @test isfile(joinpath(stopped_settings.output, "raw.txt"))
        @test isfile(stopped_settings.checkpoint)
        @test any(event -> event[1] === :write && !event[3] && !event[4], stopped_events)

        allowed = parse_single_point_args(replace.(
            [arguments; "--allow_nonconverged=false"],
            "--allow_nonconverged=false" => "--allow_nonconverged=true",
        ))
        raw = run_single_point(
            allowed;
            operations=fake_single_point_operations(
                Any[];
                converged=false,
                transfer_valid=false,
            ),
        )
        @test !raw.valid
        @test !raw.optimization.converged
        @test !raw.transfer.valid
    end
end

@testset "callable flux scan reuses site IDs and writes deterministic tables" begin
    mktempdir() do directory
        settings = parse_flux_scan_args([
            "--Ly=6",
            "--x_period=1",
            "--filling_num=1",
            "--filling_den=3",
            "--phi_start=0",
            "--phi_stop=1",
            "--phi_steps=3",
            "--maxdim=1",
            "--output=$(joinpath(directory, "scan"))",
            "--branch_mode=adiabatic",
            "--cold_patterns=2;3",
            "--allow_nonconverged=true",
        ])
        events = Any[]
        scan = run_flux_scan(settings; operations=fake_single_point_operations(events))
        @test scan isa FluxScanRunResult
        @test length(scan.rows) == 3
        @test length(scan.events) == 2
        @test scan.rows[1].phi_y == 0.0
        @test scan.rows[end].phi_y == 1.0
        @test all(row -> row.delta_raw_schmidt_polarization == 0.0, scan.rows)
        @test all(row -> row.selected_label == "warm", scan.rows)
        @test [row.selection_mode for row in scan.rows] ==
            ["ground", "adiabatic", "adiabatic"]
        @test all(row -> row.cut_x == 1, scan.rows)
        @test scan.rows[1].fidelity_density_x_valid === false
        @test all(row -> row.fidelity_density_x_valid, scan.rows[2:end])
        @test all(row -> row.energy_delta_to_ground_valid, scan.rows)
        @test all(row -> row.energy_delta_to_ground_raw == 0.0, scan.rows)
        @test length(scan.sector_gauge) == 3
        @test all(row -> row.valid, scan.sector_gauge)

        builds = filter(event -> event[1] === :build, events)
        @test length(builds) == 3
        @test all(event -> event[3] == builds[1][3], builds)
        @test count(event -> event[1] === :initialize, events) == 1
        @test only(filter(event -> event[1] === :initialize, events))[3] == [1]
        @test count(event -> event[1] === :fidelity, events) == 6

        for point in 0:2, candidate in
            ("candidate_warm", "candidate_cold_01", "candidate_cold_02")
            candidate_directory = joinpath(
                settings.point.output,
                "phi_$(lpad(point, 3, '0'))",
                candidate,
            )
            @test isfile(joinpath(candidate_directory, "raw.txt"))
            @test isfile(joinpath(candidate_directory, "state.h5"))
        end
        summary_lines = readlines(joinpath(settings.point.output, "scan_summary.tsv"))
        event_lines = readlines(joinpath(settings.point.output, "branch_events.tsv"))
        gauge_lines = readlines(joinpath(settings.point.output, "sector_gauge.tsv"))
        metadata = TOML.parsefile(joinpath(settings.point.output, "scan_metadata.toml"))
        @test first(summary_lines) == SCAN_SUMMARY_HEADER
        @test first(event_lines) == BRANCH_EVENTS_HEADER
        @test first(gauge_lines) == SECTOR_GAUGE_HEADER
        @test length(summary_lines) == 4
        @test length(event_lines) == 3
        @test length(gauge_lines) == 4
        @test metadata["branch"]["mode"] == "adiabatic"
        @test metadata["branch"]["cut_x"] == 1
        @test metadata["scan"]["phi_steps"] == 3
        @test occursin("selected_candidate", first(summary_lines))
        @test occursin("selection_mode", first(summary_lines))
        @test occursin("fidelity_density_x_raw", first(summary_lines))
        @test occursin("energy_delta_to_ground_raw", first(summary_lines))
        @test occursin("raw_schmidt_polarization", first(summary_lines))
        @test occursin("fidelity_raw", first(event_lines))
        @test occursin("fidelity_valid", first(event_lines))
        @test occursin("spectrum_levels", first(event_lines))
        @test !occursin("unwrap", read(joinpath(settings.point.output, "scan_summary.tsv"), String))
    end
end

@testset "one-point scan persists metadata without branch events" begin
    mktempdir() do directory
        settings = parse_flux_scan_args([
            "--Ly=6",
            "--x_period=1",
            "--filling_num=1",
            "--filling_den=3",
            "--phi_start=0",
            "--phi_stop=0",
            "--phi_steps=1",
            "--maxdim=1",
            "--output=$(joinpath(directory, "single_flux"))",
            "--branch_mode=adiabatic",
        ])
        scan = run_flux_scan(
            settings;
            operations=fake_single_point_operations(Any[]),
        )
        @test length(scan.rows) == 1
        @test isempty(scan.events)
        metadata = TOML.parsefile(
            joinpath(settings.point.output, "scan_metadata.toml")
        )
        @test metadata["branch"]["mode"] == "adiabatic"
        @test metadata["branch"]["cut_x"] == 1
        @test metadata["scan"] == Dict(
            "phi_start" => 0.0,
            "phi_stop" => 0.0,
            "phi_steps" => 1,
        )
        @test readlines(joinpath(settings.point.output, "branch_events.tsv")) ==
            [BRANCH_EVENTS_HEADER]
        @test readlines(joinpath(settings.point.output, "sector_gauge.tsv")) ==
            [SECTOR_GAUGE_HEADER]
    end
end

@testset "flux scan preserves invalid raw failure paths" begin
    mktempdir() do directory
        base_arguments = [
            "--Ly=6",
            "--x_period=1",
            "--filling_num=1",
            "--filling_den=3",
            "--phi_start=0",
            "--phi_stop=1",
            "--phi_steps=2",
            "--maxdim=1",
            "--output=$(joinpath(directory, "allowed"))",
            "--branch_mode=adiabatic",
            "--allow_nonconverged=true",
        ]
        settings = parse_flux_scan_args(base_arguments)
        operations = fake_single_point_operations(
            Any[];
            converged=false,
            transfer_valid=false,
        )
        allowed = run_flux_scan(settings; operations)
        @test length(allowed.rows) == 2
        @test all(row -> !row.valid && !row.converged, allowed.rows)
        @test all(row -> !row.selection_valid, allowed.rows)
        @test isfile(joinpath(settings.point.output, "scan_summary.tsv"))
        @test isfile(joinpath(settings.point.output, "branch_events.tsv"))

        rejected_settings = parse_flux_scan_args(replace.(
            base_arguments,
            "--output=$(joinpath(directory, "allowed"))" =>
                "--output=$(joinpath(directory, "rejected"))",
            "--allow_nonconverged=true" => "--allow_nonconverged=false",
        ))
        rejected_error = try
            run_flux_scan(rejected_settings; operations)
            nothing
        catch exception
            exception
        end
        @test rejected_error isa WorkflowValidationError
        @test isfile(joinpath(rejected_settings.point.output, "scan_summary.tsv"))
        @test isfile(joinpath(rejected_settings.point.output, "branch_events.tsv"))
        @test occursin(
            "selection_valid",
            first(readlines(joinpath(rejected_settings.point.output, "scan_summary.tsv"))),
        )

        events = Any[]
        valid_ops = fake_single_point_operations(events)
        throwing_fidelity_ops = SinglePointOperations(
            configure_threads=valid_ops.configure_threads,
            initialize=valid_ops.initialize,
            load_state=valid_ops.load_state,
            state_sites=valid_ops.state_sites,
            build_hamiltonian=valid_ops.build_hamiltonian,
            optimize=valid_ops.optimize,
            energy=valid_ops.energy,
            density=valid_ops.density,
            entanglement=valid_ops.entanglement,
            transfer=valid_ops.transfer,
            fidelity=(args...) -> error("synthetic fidelity failure"),
            write_outputs=valid_ops.write_outputs,
            save_state=valid_ops.save_state,
            cold_state=valid_ops.cold_state,
        )
        fidelity_settings = parse_flux_scan_args(replace.(
            base_arguments,
            "--output=$(joinpath(directory, "allowed"))" =>
                "--output=$(joinpath(directory, "fidelity"))",
        ))
        fidelity_scan = run_flux_scan(fidelity_settings; operations=throwing_fidelity_ops)
        @test isnan(fidelity_scan.rows[2].fidelity_raw)
        @test !fidelity_scan.rows[2].fidelity_valid
        @test occursin("synthetic fidelity failure", fidelity_scan.rows[2].reason)
        @test !fidelity_scan.rows[2].selection_valid

        cfg = InfiniteCylinderConfig(; Ly=6, x_period=1)
        _, _, psi = initial_infinite_mps(cfg)
        stopped = VUMPSResult(psi, VUMPSRecord[], true, "converged")
        invalid_result = SinglePointRunResult(
            parse_single_point_args([
                "--Ly=6",
                "--x_period=1",
                "--filling_num=1",
                "--filling_den=3",
                "--phi_y=0",
                "--maxdim=1",
                "--output=unused",
            ]),
            stopped,
            normalize_energy(cfg, -1.0),
            density_data(psi, cfg),
            EntanglementData[],
            synthetic_transfer(cfg; valid=false),
            false,
            ["entanglement failed: synthetic observable failure"],
        )
        reasoned = InfiniteCylinderDMRG._candidate_data(invalid_result, "invalid")
        @test occursin("synthetic observable failure", reasoned.reason)
        @test !reasoned.valid
    end
end
