using Test

include(joinpath(@__DIR__, "plot_optical_components.jl"))


function write_sample_curves(path; columns=9)
    rows = [
        0.0  -4.0  0.0  -4.0  0.0  0.06  0.0  0.06  0.0;
        0.5  -0.2 -1.0  -0.2 -1.0  0.30 -0.4  0.30 -0.4;
        1.0   0.5  0.2   0.5  0.2  0.10  0.1  0.10  0.1;
    ]
    open(path, "w") do io
        println(io, "# synthetic optical-response table")
        for row in eachrow(rows)
            println(io, join(row[1:columns], ' '))
        end
    end
end


@testset "optical component plot data" begin
    mktempdir() do directory
        data_path = joinpath(directory, "curves.dat")
        write_sample_curves(data_path)

        data = load_optical_curves(data_path)
        @test size(data) == (3, 9)

        specs = component_specs()
        @test getproperty.(specs, :key) == [
            "regular_real", "regular_imag", "total_real", "total_imag"
        ]
        @test [(spec.exact_column, spec.lanczos_column) for spec in specs] == [
            (6, 8), (7, 9), (2, 4), (3, 5)
        ]

        malformed_path = joinpath(directory, "malformed.dat")
        write_sample_curves(malformed_path; columns=8)
        @test_throws ArgumentError load_optical_curves(malformed_path)

        output_dir = joinpath(directory, "figures")
        outputs = generate_optical_component_plots(
            data_path, output_dir; lanczos_m=892
        )
        @test length(outputs) == 4
        @test all(isfile, outputs)
        @test all(path -> filesize(path) > 1000, outputs)
    end
end
