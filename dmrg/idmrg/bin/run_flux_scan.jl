using InfiniteCylinderDMRG

function run_flux_scan_main(args=ARGS)
    settings = parse_flux_scan_args(args)
    run_flux_scan(settings)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_flux_scan_main()
end
