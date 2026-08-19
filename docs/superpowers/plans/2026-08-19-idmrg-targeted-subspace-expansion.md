# Targeted iDMRG Subspace Expansion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make every VUMPS subspace-expansion stage reach its requested bond dimension or fail before optimization, then prove the behavior on the paper QN Hamiltonian and a real D=32 candidate.

**Architecture:** Strengthen the existing `expand_subspace` boundary instead of changing Fig. 2 scheduling or injecting artificial sectors. The function repeatedly applies the upstream Hamiltonian-generated expansion until the target maximum link dimension is reached, with monotonic-progress and overshoot checks on every pass. Existing candidate, checkpoint, observable, and selection code remains unchanged.

**Tech Stack:** Julia 1.12.5, ITensors.jl, ITensorMPS.jl, ITensorInfiniteMPS.jl, HDF5, PBS/Torque on W003 `cmt`.

---

### Task 1: Add target-reach regressions

**Files:**
- Modify: `dmrg/idmrg/test/test_vumps_runner.jl`

- [ ] **Step 1: Extend the small Ising expansion test**

After the existing target-2 assertion in `checked subspace expansion`, add:

```julia
expanded_to_four = expand_subspace(psi, H, 4; cutoff=1e-8)
@test maximum(link_dimensions(expanded_to_four)) == 4
@test !InfiniteCylinderDMRG._subspace_expansion_progressed(
    [2, 2], [1, 3]
)
```

- [ ] **Step 2: Add the real paper-QN regression**

Add this testset after `paper D=2 VUMPS canonicalization regression`:

```julia
@testset "paper QN subspace expansion reaches requested target" begin
    manifest = joinpath(@__DIR__, "..", "benchmarks", "fqahc_fig2.toml")
    spec = load_fig2_benchmark(manifest)
    candidate_id = first(
        InfiniteCylinderDMRG._default_fig2_candidate_ids(
            spec, 4, 1, nothing
        ),
    )
    saved_default_rng = copy(Random.default_rng())
    saved_index_rng = copy(ITensors.index_id_rng())
    try
        Random.seed!(Random.default_rng(), 0)
        Random.seed!(ITensors.index_id_rng(), 0)
        prepared = InfiniteCylinderDMRG._prepare_fig2_candidate_state(
            spec, spec.config, candidate_id, nothing
        )
        hamiltonian = build_infinite_mpo(
            spec.config, spec.model, prepared.sites
        )
        expected_sites = siteinds(only, prepared.psi.AL)
        expected_flux = flux(prepared.psi.AL)
        expanded = expand_subspace(
            prepared.psi, hamiltonian, 4; cutoff=1.0e-9
        )

        @test maximum(link_dimensions(expanded)) == 4
        @test nsites(expanded) == sites_per_cell(spec.config)
        @test siteinds(only, expanded.AL) == expected_sites
        @test flux(expanded.AL) == expected_flux
    finally
        copy!(Random.default_rng(), saved_default_rng)
        copy!(ITensors.index_id_rng(), saved_index_rng)
    end
end
```

- [ ] **Step 3: Commit and synchronize the RED test**

```bash
git add dmrg/idmrg/test/test_vumps_runner.jl
git commit -m "test: require targeted iDMRG subspace growth"
```

Transfer the exact commit to the dedicated W003 checkout with a Git bundle,
fast-forward it there, and push `DMRG` using W003's existing SSH key. Do not
include the user's finite-DMRG worktree changes.

- [ ] **Step 4: Run the targeted test through PBS and verify RED**

Submit `test/test_vumps_runner.jl` using a read-only launcher derived from
`dmrg/idmrg/jobs/run_tests.pbs` with a 45-minute walltime.

Expected: the new Ising and/or paper-QN assertion fails because the current
one-pass implementation returns maximum link dimension 2 for target 4. No
Julia command runs on the Mac.

### Task 2: Enforce the target in `expand_subspace`

**Files:**
- Modify: `dmrg/idmrg/src/VUMPSRunner.jl:276-293`
- Test: `dmrg/idmrg/test/test_vumps_runner.jl`

- [ ] **Step 1: Make the progress predicate monotone**

After its existing length check, replace the predicate body with:

```julia
return all(after[index] >= before[index] for index in eachindex(before, after)) &&
    any(after[index] > before[index] for index in eachindex(before, after))
```

This makes each accepted pass increase the bounded sum of link dimensions;
growth of one link cannot hide shrinkage of another.

- [ ] **Step 2: Replace one-pass expansion with the bounded progress loop**

Keep the existing argument validation and early return, then implement:

```julia
current = psi
while maximum(link_dimensions(current)) < target
    before = link_dimensions(current)
    expanded = subspace_expansion(current, H; maxdim=target, cutoff)
    after = link_dimensions(expanded)
    _subspace_expansion_progressed(before, after) || error(
        "subspace expansion stalled at achieved maxdim=$(maximum(after)) " *
        "before target maxdim=$target; before=$before after=$after"
    )
    maximum(after) <= target || error(
        "subspace expansion exceeded target maxdim=$target: after=$after"
    )
    current = expanded
end
maximum(link_dimensions(current)) == target || error(
    "subspace expansion did not reach target maxdim=$target"
)
return current
```

- [ ] **Step 3: Tighten the stalled-Hamiltonian assertion**

Replace the existing broad message check with:

```julia
@test occursin("stalled at achieved maxdim=1", sprint(showerror, error))
@test occursin("target maxdim=2", sprint(showerror, error))
```

- [ ] **Step 4: Run the same targeted PBS test and verify GREEN**

Expected: all `test_vumps_runner.jl` testsets pass; the paper-QN target-4 case
preserves sites and flux; the zero Hamiltonian still fails before VUMPS.

- [ ] **Step 5: Commit and synchronize the implementation**

```bash
git add dmrg/idmrg/src/VUMPSRunner.jl dmrg/idmrg/test/test_vumps_runner.jl
git commit -m "fix: reach requested iDMRG expansion dimension"
```

Synchronize the exact commit through W003 and push `DMRG` with the cluster SSH
key.

### Task 3: Verify one real D=32 Fig. 2 candidate

**Files:**
- Read: `dmrg/idmrg/src/Fig2Benchmark.jl`
- Output: `/home/public/shajy/codex/results/fqahc-fig2/debug/targeted-d32-candidate`

- [ ] **Step 1: Submit a single-candidate PBS probe**

Use a read-only PBS launcher on `cmt`, one node and 24 cores, with a 12-hour
walltime. The Julia entrypoint loads the immutable Fig. 2 manifest and calls
`InfiniteCylinderDMRG._default_fig2_run_candidate` for dimension 32, point 1,
flux 0, and candidate `cdw_t0_dopant_y0` into the dedicated debug output.
Export the same Julia depot and thread limits as the production runner.

- [ ] **Step 2: Audit the candidate gate**

Require all of the following before resuming the pilot:

```text
requested_maxdim = 32
achieved_maxlinkdim = 32
checkpoint_maxlinkdim = 32
converged = true
valid = true
restart_valid = true
```

Also require nonempty `state.h5`, `convergence.tsv`, `transfer_spectrum.tsv`,
`entanglement_spectrum.tsv`, `schmidt_sectors.tsv`, and candidate diagnostics.
Record energy density, leading correlation length, raw Schmidt polarization,
sector weights, momentum validity/reason, walltime, CPU utilization, and peak
RSS. A momentum diagnostic may be invalid, but it must report a concrete raw
reason rather than be silently omitted.

- [ ] **Step 3: If the probe stalls, preserve evidence and return to design**

Do not add random sectors or relax the target gate. Preserve the output and
error, then compare powers-of-two VUMPS staging against a QN-aware deterministic
enrichment design before further implementation.

### Task 4: Run full verification and restore pilot production

**Files:**
- Read: `dmrg/idmrg/test/runtests.jl`
- Output: `/home/public/shajy/codex/results/fqahc-fig2/pilot`
- Archive: `/home/public/shajy/codex/results/fqahc-fig2/pilot-d2-red-1869087`

- [ ] **Step 1: Run the complete iDMRG suite on W003**

Submit `dmrg/idmrg/jobs/run_tests.pbs` to `cmt` with 24 cores. Expected: every
test passes, including checkpoint restart, canonical residual `1e-10`, Fig. 2
job contracts, and the new target-growth regression.

- [ ] **Step 2: Archive the invalid pilot without deleting it**

After the single-candidate and full-suite gates pass, move the stopped pilot
directory atomically to `pilot-d2-red-1869087`. Verify that it still contains
the three original `state.h5` files and its ledger. Create a fresh empty
`pilot` path through the production runner; do not edit or reuse the invalid
ledger.

- [ ] **Step 3: Submit and monitor the fresh pilot**

Run `dmrg/idmrg/jobs/submit_fig2_stage.sh`. Verify checkout audit, exact Git
commit, exported depot, immutable config/launcher, per-candidate target bond
dimension, checkpoint reload, energy density, transfer spectrum/correlation
length, entanglement spectrum, Schmidt sectors, raw/tracked pump, branch
fidelity, and final acceptance report.

- [ ] **Step 4: Review and commit any documentation updates**

Update the iDMRG README/report with the expansion-target failure mode, the
strong target contract, W003 job IDs, and measured benchmark resources. Run
`git diff --check`, request spec and code-quality review, and synchronize the
final commit without touching finite-DMRG changes.
