using Test
using TOML
using HDF5
using ITensors
using ITensorMPS
using InfiniteCylinderDMRG

function captured_exception(f)
    try
        f()
        return nothing
    catch error
        return error
    end
end

@testset "validated atomic HDF5 checkpoints" begin
    mktempdir() do directory
        cfg = InfiniteCylinderConfig(; Ly=6, x_period=1, phi_y=0.25)
        _, _, psi = initial_infinite_mps(cfg)
        path = joinpath(directory, "nested", "state.h5")

        save_checkpoint(path, psi, cfg)
        @test isfile(path)
        @test !ispath(path * ".tmp")

        loaded = load_checkpoint(path, cfg)
        @test loaded isa InfiniteCanonicalMPS
        before_entanglement = entanglement_data(psi, cfg; cut_x=1)
        after_entanglement = entanglement_data(loaded, cfg; cut_x=1)
        @test isapprox(
            before_entanglement.entropy,
            after_entanglement.entropy;
            atol=1e-10,
            rtol=0,
        )
        before_density = density_data(psi, cfg)
        after_density = density_data(loaded, cfg)
        @test all(
            isapprox(before.density, after.density; atol=1e-10, rtol=0) for
            (before, after) in zip(before_density, after_density)
        )

        h5open(path, "r") do file
            metadata = attributes(file)
            @test read(metadata["format"]) == "infinite_cylinder_vumps_v1"
            @test read(metadata["backend_commit"]) ==
                "765f2777703bc1138b009adbed1b97bde1973402"
            @test read(metadata["Ly"]) == cfg.Ly
            @test read(metadata["x_period"]) == cfg.x_period
            @test read(metadata["filling_num"]) == cfg.filling_num
            @test read(metadata["filling_den"]) == cfg.filling_den
            @test read(metadata["phi_y"]) === cfg.phi_y
            @test read(metadata["sites_per_cell"]) == sites_per_cell(cfg)
            @test read(metadata["particles_per_cell"]) == particles_per_cell(cfg)
            @test read(metadata["charge_scale"]) == charge_scale(cfg)
        end

        incompatible = (
            InfiniteCylinderConfig(; Ly=2, x_period=3, phi_y=0.25),
            InfiniteCylinderConfig(6, 1, 2, 3, 0.25),
            InfiniteCylinderConfig(; Ly=6, x_period=1, phi_y=nextfloat(0.25)),
        )
        for mismatch in incompatible
            error = captured_exception(() -> load_checkpoint(path, mismatch))
            @test error isa CheckpointCompatibilityError
            @test occursin("checkpoint configuration mismatch", sprint(showerror, error))
        end

        missing_metadata = joinpath(directory, "missing.h5")
        h5open(missing_metadata, "w") do file
            attributes(file)["format"] = "infinite_cylinder_vumps_v1"
        end
        missing_error =
            captured_exception(() -> load_checkpoint(missing_metadata, cfg))
        @test missing_error isa CheckpointFormatError
        @test occursin("backend_commit", sprint(showerror, missing_error))

        corrupt_metadata = joinpath(directory, "corrupt.h5")
        cp(path, corrupt_metadata)
        h5open(corrupt_metadata, "r+") do file
            delete_attribute(file, "sites_per_cell")
            attributes(file)["sites_per_cell"] = "six"
        end
        corrupt_error =
            captured_exception(() -> load_checkpoint(corrupt_metadata, cfg))
        @test corrupt_error isa CheckpointFormatError
        @test occursin("sites_per_cell", sprint(showerror, corrupt_error))

        boolean_metadata = joinpath(directory, "boolean.h5")
        cp(path, boolean_metadata)
        h5open(boolean_metadata, "r+") do file
            delete_attribute(file, "x_period")
            attributes(file)["x_period"] = true
        end
        boolean_error =
            captured_exception(() -> load_checkpoint(boolean_metadata, cfg))
        @test boolean_error isa CheckpointFormatError
        @test occursin("x_period", sprint(showerror, boolean_error))

        corrupt_state = joinpath(directory, "noncanonical.h5")
        cp(path, corrupt_state)
        noncanonical = copy(psi)
        noncanonical.C[1] *= 2
        h5open(corrupt_state, "r+") do file
            delete_object(file, "state")
            write(file, "state", noncanonical)
        end
        state_error = captured_exception(() -> load_checkpoint(corrupt_state, cfg))
        @test state_error isa CheckpointFormatError
        @test occursin("canonical", sprint(showerror, state_error))

        zero_state = joinpath(directory, "zero.h5")
        cp(path, zero_state)
        zero_center = copy(psi)
        for site in 1:sites_per_cell(cfg)
            zero_center.C[site] *= 0
        end
        h5open(zero_state, "r+") do file
            delete_object(file, "state")
            write(file, "state", zero_center)
        end
        zero_error = captured_exception(() -> load_checkpoint(zero_state, cfg))
        @test zero_error isa CheckpointFormatError
        @test occursin("center", sprint(showerror, zero_error))
    end
end

@testset "concurrent atomic replacements use private temporary files" begin
    mktempdir() do directory
        path = joinpath(directory, "shared.txt")
        entered = Channel{Nothing}(2)
        release = Channel{Nothing}(2)

        function coordinated_write(contents)
            try
                return InfiniteCylinderDMRG._atomic_replace(path) do temporary
                    open(temporary, "w") do io
                        write(io, contents)
                        put!(entered, nothing)
                        take!(release)
                        write(io, "-complete")
                    end
                end
            catch error
                return error
            end
        end

        first = @async coordinated_write("first")
        take!(entered)
        second = @async coordinated_write("second")
        take!(entered)
        put!(release, nothing)
        yield()
        put!(release, nothing)

        @test fetch(first) == path
        @test fetch(second) == path
        @test read(path, String) in ("first-complete", "second-complete")
        @test isempty(filter(name -> occursin(".tmp", name), readdir(directory)))
    end
end

@testset "deterministic text output contract" begin
    mktempdir() do directory
        cfg = InfiniteCylinderConfig(; Ly=6, x_period=1, phi_y=0.125)
        _, _, psi = initial_infinite_mps(cfg)
        records = [
            VUMPSRecord(
                1,
                1,
                4,
                -1.25,
                -1.24,
                0.01,
                Inf,
                2e-4,
                3e-4,
                3e-4,
                0.5,
                false,
            ),
            VUMPSRecord(
                1,
                2,
                4,
                -1.26,
                -1.26,
                0.0,
                0.01,
                1e-6,
                1e-6,
                1e-6,
                0.75,
                true,
            ),
        ]
        result = VUMPSResult(psi, records, true, "converged after 1 stage")
        energy = normalize_energy(cfg, -1.26)
        densities = density_data(psi, cfg)
        entanglements = [entanglement_data(psi, cfg; cut_x=1)]
        transfer = InfiniteCylinderDMRG._neutral_transfer_result(
            cfg,
            ComplexF64[1.0, 0.5im],
            [1e-12, 2e-12],
            2;
            residual_tolerance=1e-10,
        )

        write_output_files(
            directory,
            cfg,
            result,
            energy,
            densities,
            entanglements,
            transfer;
            optimization=(
                maxdim_schedule=[1, 4],
                cutoff=1e-8,
                vumps_tol=1e-5,
            ),
            dependencies=(julia_version=string(VERSION),),
        )

        expected_headers = Dict(
            "convergence.tsv" =>
                "stage\titeration\tmaxlinkdim\tenergy_left\tenergy_right\tenergy_mismatch\tdelta_energy\teps_left\teps_right\tprecision_error\telapsed_seconds\tconverged",
            "density.tsv" => "site\tx\ty\tdensity",
            "entanglement_spectrum.tsv" =>
                "cut_x\tbond\tlevel\tsingular_value\tprobability\tentanglement_energy\tqn\traw_charge\tphysical_charge",
            "schmidt_sectors.tsv" =>
                "cut_x\tbond\tqn\traw_charge\tphysical_charge\tweight",
            "transfer_spectrum.tsv" =>
                "level\teigenvalue_real\teigenvalue_imag\tmagnitude\tphase\tresidual_norm\tconverged\tvalid",
        )
        for (filename, header) in expected_headers
            lines = readlines(joinpath(directory, filename))
            @test first(lines) == header
            @test !isempty(lines[2:end])
            @test !ispath(joinpath(directory, filename * ".tmp"))
        end
        @test occursin("\tfalse", read(joinpath(directory, "convergence.tsv"), String))
        @test occursin("\ttrue", read(joinpath(directory, "convergence.tsv"), String))

        summary_path = joinpath(directory, "summary.toml")
        summary_text = read(summary_path, String)
        summary = TOML.parse(summary_text)
        @test summary["format"] == "infinite_cylinder_vumps_text_v1"
        @test summary["algorithm"] == "VUMPS"
        @test summary["valid"] === true
        @test summary["converged"] === true
        @test summary["configuration"]["phi_y"] === cfg.phi_y
        @test summary["energy"]["per_cell"] === energy.per_cell
        @test summary["observables"]["transfer_valid"] === true
        @test summary["entanglement"][1]["cut_x"] == entanglements[1].cut_x
        @test summary["entanglement"][1]["bond"] == entanglements[1].bond
        @test summary["entanglement"][1]["valid"] === true
        @test summary["entanglement"][1]["entropy"] ===
            entanglements[1].entropy
        @test summary["entanglement"][1]["raw_schmidt_polarization"] ===
            entanglements[1].raw_schmidt_polarization
        @test summary["optimization"]["maxdim_schedule"] == [1, 4]
        @test !ispath(summary_path * ".tmp")

        convergence_fields =
            split(readlines(joinpath(directory, "convergence.tsv"))[2], '\t')
        @test length(convergence_fields) == 12
        @test convergence_fields == [
            "1",
            "1",
            "4",
            "-1.25",
            "-1.24",
            "0.01",
            "Inf",
            "0.0002",
            "0.0003",
            "0.0003",
            "0.5",
            "false",
        ]
        density_fields = split(readlines(joinpath(directory, "density.tsv"))[2], '\t')
        @test length(density_fields) == 4
        @test parse(Int, density_fields[1]) == densities[1].site
        @test parse(Int, density_fields[2]) == densities[1].x
        @test parse(Int, density_fields[3]) == densities[1].y
        @test parse(Float64, density_fields[4]) == densities[1].density
        entanglement_fields = split(
            readlines(joinpath(directory, "entanglement_spectrum.tsv"))[2], '\t'
        )
        @test length(entanglement_fields) == 9
        @test parse(Int, entanglement_fields[1]) == entanglements[1].levels[1].cut_x
        @test parse(Int, entanglement_fields[2]) == entanglements[1].levels[1].bond
        @test parse(Int, entanglement_fields[3]) == entanglements[1].levels[1].level
        @test parse(Float64, entanglement_fields[4]) ==
            entanglements[1].levels[1].singular_value
        sector_fields =
            split(readlines(joinpath(directory, "schmidt_sectors.tsv"))[2], '\t')
        @test length(sector_fields) == 6
        @test parse(Int, sector_fields[1]) == entanglements[1].sectors[1].cut_x
        @test parse(Int, sector_fields[2]) == entanglements[1].sectors[1].bond
        @test parse(Float64, sector_fields[6]) == entanglements[1].sectors[1].weight
        transfer_fields =
            split(readlines(joinpath(directory, "transfer_spectrum.tsv"))[2], '\t')
        @test length(transfer_fields) == 8
        @test parse(Float64, transfer_fields[2]) == real(transfer.eigenvalues[1])
        @test parse(Float64, transfer_fields[3]) == imag(transfer.eigenvalues[1])
        @test transfer_fields[7:8] == ["true", "true"]

        first_render = Dict(
            filename => read(joinpath(directory, filename), String) for
            filename in ("summary.toml", keys(expected_headers)...)
        )
        write_output_files(
            directory,
            cfg,
            result,
            energy,
            densities,
            entanglements,
            transfer;
            optimization=(
                maxdim_schedule=[1, 4],
                cutoff=1e-8,
                vumps_tol=1e-5,
            ),
            dependencies=(julia_version=string(VERSION),),
        )
        for (filename, contents) in first_render
            @test read(joinpath(directory, filename), String) == contents
        end

        invalid_transfer = InfiniteCylinderDMRG._neutral_transfer_result(
            cfg,
            ComplexF64[1.0],
            [1e-12],
            1;
            residual_tolerance=1e-10,
        )
        stopped = VUMPSResult(psi, records[1:1], false, "iteration limit reached")
        invalid_directory = joinpath(directory, "invalid")
        write_output_files(
            invalid_directory,
            cfg,
            stopped,
            energy,
            densities,
            entanglements,
            invalid_transfer,
        )
        invalid_summary = TOML.parsefile(joinpath(invalid_directory, "summary.toml"))
        @test invalid_summary["valid"] === false
        @test invalid_summary["converged"] === false
        @test invalid_summary["energy"]["valid"] === false
        @test invalid_summary["observables"]["density_valid"] === false
        @test invalid_summary["observables"]["entanglement_valid"] === false
        @test invalid_summary["observables"]["transfer_valid"] === false
        @test invalid_summary["entanglement"][1]["valid"] === false
        @test isnan(invalid_summary["entanglement"][1]["entropy"])
        @test isnan(
            invalid_summary["entanglement"][1]["raw_schmidt_polarization"]
        )
        @test isnan(
            invalid_summary["observables"]["correlation_length_neutral_cell"]
        )
        invalid_transfer_text =
            read(joinpath(invalid_directory, "transfer_spectrum.tsv"), String)
        @test occursin("\tfalse\tfalse", invalid_transfer_text)
        @test endswith(
            readlines(joinpath(invalid_directory, "density.tsv"))[2], "\tNaN"
        )
        invalid_entanglement_fields = split(
            readlines(joinpath(invalid_directory, "entanglement_spectrum.tsv"))[2],
            '\t',
        )
        @test all(field -> field == "NaN", invalid_entanglement_fields[4:6])
        @test invalid_entanglement_fields[9] == "NaN"
        invalid_sector_fields = split(
            readlines(joinpath(invalid_directory, "schmidt_sectors.tsv"))[2], '\t'
        )
        @test invalid_sector_fields[5:6] == ["NaN", "NaN"]

        original_level = entanglements[1].levels[1]
        invalid_level = EntanglementLevel(
            original_level.cut_x,
            original_level.bond,
            original_level.level,
            original_level.singular_value,
            original_level.probability,
            NaN,
            original_level.qn,
            original_level.raw_charge,
            original_level.physical_charge,
        )
        invalid_entanglement = EntanglementData(
            entanglements[1].cut_x,
            entanglements[1].bond,
            entanglements[1].entropy,
            entanglements[1].raw_schmidt_polarization,
            [invalid_level; entanglements[1].levels[2:end]],
            entanglements[1].sectors,
        )
        invalid_entanglement_directory = joinpath(directory, "invalid_entanglement")
        write_output_files(
            invalid_entanglement_directory,
            cfg,
            result,
            energy,
            densities,
            [invalid_entanglement],
            transfer,
        )
        invalid_entanglement_summary =
            TOML.parsefile(joinpath(invalid_entanglement_directory, "summary.toml"))
        @test invalid_entanglement_summary["valid"] === false
        @test invalid_entanglement_summary["observables"]["entanglement_valid"] ===
            false
        @test invalid_entanglement_summary["entanglement"][1]["valid"] === false

        invalid_phase_transfer = NeutralTransferData(
            transfer.eigenvalues,
            [NaN; transfer.phases[2:end]],
            transfer.residual_norms,
            transfer.converged,
            transfer.valid,
            transfer.reason,
            transfer.ratio,
            transfer.xi_cell,
            transfer.xi_x,
        )
        invalid_phase_directory = joinpath(directory, "invalid_phase")
        write_output_files(
            invalid_phase_directory,
            cfg,
            result,
            energy,
            densities,
            entanglements,
            invalid_phase_transfer,
        )
        invalid_phase_summary =
            TOML.parsefile(joinpath(invalid_phase_directory, "summary.toml"))
        @test invalid_phase_summary["valid"] === false
        @test invalid_phase_summary["observables"]["transfer_valid"] === false
        @test occursin(
            "\tfalse\tfalse",
            read(joinpath(invalid_phase_directory, "transfer_spectrum.tsv"), String),
        )
    end
end
