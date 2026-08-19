# Targeted subspace expansion for the infinite-cylinder runner

## Problem and evidence

The Fig. 2 pilot job `1869087.w003` requested bond dimension 32, but its
first three completed candidates all reported `achieved_maxlinkdim = 2`.
Their `expansion.tsv` files show one expansion from product-state links of
dimension 1 to links of dimension at most 2, followed by a converged VUMPS
stage. The selection contract correctly rejects these states because the
achieved and requested dimensions differ.

`expand_subspace` currently treats growth of any link as success. For the
long-range, QN-conserving paper Hamiltonian, one application of
`subspace_expansion` can add only part of the Hamiltonian-generated Krylov
space. Progress therefore does not imply that the requested target was
reached.

## Chosen design

`expand_subspace` will own a strong target contract. If the current maximum
link dimension is below the target, it will repeatedly call the upstream
`subspace_expansion` operation with the same target and cutoff. After every
call it will:

1. validate that the link count is unchanged;
2. require at least one link dimension to grow; and
3. stop only when the maximum link dimension equals the target.

Each successful pass increases an integer link dimension and the upstream
operation is capped by `maxdim=target`, so the loop is finite. If a pass
stalls before the target, the runner fails immediately with the before,
after, achieved, and target dimensions. It must not run VUMPS on a state
that merely made partial expansion progress.

The public `SubspaceExpansionRecord` remains one record per VUMPS schedule
stage. Its `before` and `after` fields describe the complete targeted
expansion, and its elapsed time includes all internal passes. No Hamiltonian,
QN convention, optimizer tolerance, candidate selection rule, flux rule, or
raw pump observable changes.

## Rejected alternatives

- A powers-of-two Fig. 2 schedule would interleave full VUMPS convergence at
  every intermediate dimension. It is much more expensive at production
  dimensions and leaves the lower-level partial-target bug available to
  other callers.
- Random or manually constructed QN-sector enrichment could open unsupported
  sectors and complicate reproducibility. It is unnecessary unless repeated
  Hamiltonian-generated expansion proves unable to reach the target.

## Verification

All Julia tests run through PBS on W003's `cmt` queue.

1. Add a regression using the real paper QN Hamiltonian and a Fig. 2 product
   state with target 4. Before the fix it must fail because the result reaches
   only dimension 2; after the fix it must reach dimension 4 while preserving
   the unit cell, site indices, and conserved flux.
2. Retain a stalled-Hamiltonian test. A target that cannot be reached must
   throw before VUMPS starts, with the achieved and requested dimensions in
   the error.
3. Run the scoped runner tests and then the complete iDMRG suite on W003.
4. Run one real D=32 Fig. 2 candidate. Its candidate metadata and reloaded
   checkpoint must both report maxlinkdim 32, and it must emit energy density,
   correlation length/transfer spectrum, entanglement spectrum, Schmidt
   sectors, and momentum diagnostics.
5. Resume the pilot from its immutable ledger only after the single-candidate
   gate passes. The three D=2 candidates from job `1869087.w003` are retained
   as failure evidence but cannot be reused because they violate the target
   contract.

## Failure handling

A targeted expansion stall is a hard numerical error, not a converged low-rank
state. Candidate completion and ledger updates occur only after the target is
reached and the existing checkpoint/artifact audits pass. Existing finite-DMRG
jobs and files remain out of scope.
