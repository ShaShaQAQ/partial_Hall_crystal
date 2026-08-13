using InfiniteCylinderDMRG

function run_vumps_main(args=ARGS)
    settings = parse_single_point_args(args)
    run_single_point(settings)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_vumps_main()
end
