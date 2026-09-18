using DelimitedFiles
using Plots

const FIXED_PDF_CREATION_DATE = "20000101000000"


function normalize_pdf_metadata!(path)
    bytes = read(path)
    marker = collect(codeunits("/CreationDate (D:"))
    marker_range = findfirst(marker, bytes)
    marker_range === nothing &&
        error("generated PDF has no CreationDate metadata: $path")
    date_start = last(marker_range) + 1
    date_end = date_start + length(FIXED_PDF_CREATION_DATE) - 1
    bytes[date_start:date_end] .= codeunits(FIXED_PDF_CREATION_DATE)
    open(path, "w") do io
        write(io, bytes)
    end
    return path
end


function load_optical_curves(path)::Matrix{Float64}
    isfile(path) || throw(ArgumentError("optical curve data not found: $path"))
    data = readdlm(path, Float64; comments=true, comment_char='#')
    ndims(data) == 2 ||
        throw(ArgumentError("optical curve data must be a matrix"))
    size(data, 2) == 9 || throw(ArgumentError(
        "optical curve data must have exactly 9 numeric columns; " *
        "found $(size(data, 2))"
    ))
    size(data, 1) >= 2 || throw(ArgumentError(
        "optical curve data must contain at least two frequency samples"
    ))
    all(isfinite, data) ||
        throw(ArgumentError("optical curve data contain non-finite values"))
    all(diff(view(data, :, 1)) .> 0) ||
        throw(ArgumentError("frequency samples must be strictly increasing"))
    return Matrix{Float64}(data)
end


function component_specs()
    return [
        (
            key="regular_real",
            exact_column=6,
            lanczos_column=8,
            ylabel="Re sigma_xx^reg",
        ),
        (
            key="regular_imag",
            exact_column=7,
            lanczos_column=9,
            ylabel="Im sigma_xx^reg",
        ),
        (
            key="total_real",
            exact_column=2,
            lanczos_column=4,
            ylabel="Re sigma_xx^total",
        ),
        (
            key="total_imag",
            exact_column=3,
            lanczos_column=5,
            ylabel="Im sigma_xx^total",
        ),
    ]
end


function generate_optical_component_plots(
        data_path,
        output_dir;
        lanczos_m=892)
    data = load_optical_curves(data_path)
    mkpath(output_dir)
    omega = view(data, :, 1)
    outputs = String[]

    for spec in component_specs()
        exact = view(data, :, spec.exact_column)
        lanczos = view(data, :, spec.lanczos_column)
        absolute_error = max.(abs.(lanczos .- exact), 1e-18)

        panel = plot(
            omega,
            exact;
            color=:blue,
            linewidth=2.8,
            linestyle=:solid,
            label="exact (full sum)",
            xlabel="omega",
            ylabel=spec.ylabel,
            legend=:topleft,
            legendfontsize=7,
            guidefontsize=10,
            tickfontsize=8,
            gridalpha=0.25,
            size=(600, 400),
            inset=(1, bbox(0.54, 0.13, 0.38, 0.29)),
        )
        plot!(
            panel,
            omega,
            lanczos;
            color=:red,
            linewidth=1.5,
            linestyle=:dash,
            label="Lanczos M=$lanczos_m",
            subplot=1,
        )
        plot!(
            panel,
            omega,
            absolute_error;
            color=:black,
            linewidth=0.9,
            label=nothing,
            yscale=:log10,
            xlabel="",
            ylabel="abs diff",
            title="pointwise error",
            titlefontsize=7,
            guidefontsize=6,
            tickfontsize=5,
            gridalpha=0.20,
            framestyle=:box,
            subplot=2,
        )

        output_path = joinpath(
            output_dir, "corrected_sigma_xx_$(spec.key).pdf"
        )
        savefig(panel, output_path)
        normalize_pdf_metadata!(output_path)
        push!(outputs, output_path)
    end
    return outputs
end


if abspath(PROGRAM_FILE) == @__FILE__
    report_dir = @__DIR__
    data_path = joinpath(
        report_dir,
        "data",
        "optical_response",
        "optical_response_curves.dat",
    )
    output_dir = joinpath(report_dir, "figures", "optical_response")
    for output in generate_optical_component_plots(
            data_path, output_dir; lanczos_m=892)
        println(output)
    end
end
