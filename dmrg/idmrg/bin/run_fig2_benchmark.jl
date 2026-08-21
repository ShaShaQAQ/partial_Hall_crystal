using InfiniteCylinderDMRG

function run_fig2_benchmark_main(
    args=ARGS;
    operations=Fig2BenchmarkOperations(),
    configure_threads=configure_cli_threads,
)
    settings = parse_fig2_benchmark_args(args)
    configure_threads(settings.threads)
    spec = load_fig2_benchmark(settings.manifest)
    run = run_fig2_benchmark(
        spec,
        settings.output;
        stage=settings.stage,
        dimensions=settings.dimensions,
        fluxes=settings.fluxes,
        operations,
    )
    report = write_fig2_acceptance_report!(
        spec,
        settings.output;
        checkpoint_audit=operations.checkpoint_audit,
        progress_audit=operations.progress_audit,
        candidate_ids_provider=operations.candidate_ids,
    )
    return (; run, report)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_fig2_benchmark_main()
end
