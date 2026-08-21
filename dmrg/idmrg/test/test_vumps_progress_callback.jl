using Test
using ITensors
using ITensorMPS
using ITensorInfiniteMPS
using InfiniteCylinderDMRG

struct VUMPSProgressStop <: Exception end

@testset "VUMPS event callback exposes restartable boundaries" begin
    initspin(_) = "↑"
    sites = infsiteinds("S=1/2", 1; initstate=initspin)
    psi = InfMPS(sites, initspin)
    H = InfiniteSum{MPO}(Model("ising"), sites; J=1.0, h=1.2)
    events = Any[]

    result = run_vumps(
        H,
        psi;
        maxdim_schedule=[2],
        cutoff=1.0e-8,
        max_iterations=1,
        vumps_tol=1.0e-5,
        energy_tol=1.0e-4,
        energy_mismatch_tol=1.0e-4,
        stable_iterations=2,
        progress_callback=event -> push!(events, event),
    )

    @test !result.converged
    @test length(events) == 2
    expansion, iteration = events
    @test expansion isa VUMPSProgressEvent
    @test expansion.sequence == 1
    @test expansion.kind == :expansion
    @test expansion.stage == 1
    @test expansion.iteration == 1
    @test expansion.target == 2
    @test expansion.psi isa InfiniteCanonicalMPS
    @test expansion.expansion isa SubspaceExpansionRecord
    @test isnothing(expansion.record)
    @test iteration isa VUMPSProgressEvent
    @test iteration.sequence == 2
    @test iteration.kind == :iteration
    @test iteration.stage == 1
    @test iteration.iteration == 1
    @test iteration.target == 2
    @test iteration.psi isa InfiniteCanonicalMPS
    @test isnothing(iteration.expansion)
    @test iteration.record isa VUMPSRecord

    calls = Ref(0)
    @test_throws VUMPSProgressStop run_vumps(
        H,
        psi;
        maxdim_schedule=[1],
        cutoff=1.0e-8,
        max_iterations=1,
        vumps_tol=1.0e-5,
        energy_tol=1.0e-4,
        energy_mismatch_tol=1.0e-4,
        stable_iterations=2,
        progress_callback=event -> begin
            calls[] += 1
            throw(VUMPSProgressStop())
        end,
    )
    @test calls[] == 1
end
