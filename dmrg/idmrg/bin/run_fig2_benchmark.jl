using InfiniteCylinderDMRG

function run_fig2_benchmark_main(
    args=ARGS;
    operations=nothing,
    operations_for_backend=fig2_operations_for_backend,
    configure_threads=configure_cli_threads,
)
    settings = parse_fig2_benchmark_args(args)
    configure_threads(settings.threads)
    spec = load_fig2_benchmark(settings.manifest)
    resolved_operations = isnothing(operations) ?
        operations_for_backend(spec) : operations
    run = run_fig2_benchmark(
        spec,
        settings.output;
        stage=settings.stage,
        dimensions=settings.dimensions,
        fluxes=settings.fluxes,
        operations=resolved_operations,
    )
    report = write_fig2_acceptance_report!(
        spec,
        settings.output;
        checkpoint_audit=resolved_operations.checkpoint_audit,
        progress_audit=resolved_operations.progress_audit,
        candidate_ids_provider=resolved_operations.candidate_ids,
    )
    return (; run, report)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_fig2_benchmark_main()
end
