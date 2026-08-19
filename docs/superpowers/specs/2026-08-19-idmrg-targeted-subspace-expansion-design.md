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

The first attempted implementation made `expand_subspace` repeatedly call the
upstream operation without an intervening optimizer step. W003 job
`1869130.w003` proved that this stalls at maximum dimension 2 for the real
paper Hamiltonian. The activation probe `1869132.w003` then established:

- a second expansion without VUMPS stays at maximum dimension 2;
- after one complete VUMPS unit-cell iteration, the same target-4 expansion
  reaches maximum dimension 4; and
- additional VUMPS iterations do not change that result.

This agrees with the pinned upstream example, which alternates subspace
expansion and TDVP/VUMPS rather than applying expansion repeatedly to an
unoptimized enlarged state.

`expand_subspace` will therefore represent exactly one Hamiltonian-generated,
auditable expansion pass. It will:

1. validate that the link count is unchanged;
2. require every link dimension to be nondecreasing and at least one to grow;
   and
3. require that no link exceeds the requested target.

`run_vumps` owns the strong stage target contract. Within one requested stage
it alternates one expansion pass and one VUMPS iteration until the maximum
link dimension equals the target. Partial growth is legal, but a zero-progress
pass is a hard error. The VUMPS iteration is a complete sequential update of
the reference cell and activates the newly generated Krylov directions before
the next pass.

`max_iterations` counts every VUMPS iteration within a stage, including these
activation iterations. Energy history and the stable-iteration counter reset
after every expansion because energies across different variational spaces
are not consecutive convergence samples. A low-dimensional state must never
be declared converged even if its residual and energy tests pass. Only records
at the exact target dimension participate in the final convergence decision.
If the iteration budget expires before the target or before target-dimension
convergence, the result is nonconverged with a specific reason.

Every expansion pass gets its own `SubspaceExpansionRecord` with the same
stage and target plus its exact before/after dimensions and elapsed time. No
Hamiltonian, QN convention, optimizer tolerance, candidate selection rule,
flux rule, or raw pump observable changes.

## Rejected alternatives

- Repeated expansion without VUMPS activation was rejected by job
  `1869130.w003`: the second pass made no progress beyond dimension 2.
- A powers-of-two Fig. 2 schedule with full convergence at every intermediate
  dimension is more expensive than necessary. The activation probe shows one
  complete unit-cell VUMPS iteration is sufficient to generate the next
  expansion direction; final convergence is still enforced at the requested
  target.
- Random or manually constructed QN-sector enrichment could open unsupported
  sectors and complicate reproducibility. It is unnecessary unless repeated
  Hamiltonian-generated expansion proves unable to reach the target.

## Verification

All Julia tests run through PBS on W003's `cmt` queue.

1. Add a regression using the real paper QN Hamiltonian and a Fig. 2 product
   state with target 4. `run_vumps` must reach dimension 4 within an activation
   iteration plus a target-dimension iteration while preserving the unit cell,
   site indices, and conserved flux. The regression need not require final
   convergence with that deliberately short budget.
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

A zero-progress expansion is a hard numerical error. Exhausting the VUMPS
iteration budget before the exact target produces a nonconverged result, not a
converged low-rank state. Candidate completion and ledger updates occur only
after target-dimension convergence and the existing checkpoint/artifact audits
pass. Existing finite-DMRG jobs and files remain out of scope.
