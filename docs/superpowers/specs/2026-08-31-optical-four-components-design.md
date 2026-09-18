# Four-Component Optical-Conductivity Figure Design

Date: 2026-08-31

## Goal

Extend the corrected `4x6, Np=4, m=2` optical benchmark in the DMRG summary so
that all four longitudinal components are plotted and displayed:

1. `Re sigma_xx^reg`
2. `Im sigma_xx^reg`
3. `Re sigma_xx^total`
4. `Im sigma_xx^total`

The figures must also make the exact reference visibly distinguishable from the
final response-Lanczos result, even though their numerical difference is near
floating-point accuracy.

## Data Source and Scope

Use the existing corrected benchmark data in
`optical_response_curves.dat`. Its columns already contain the exact and final
Lanczos real and imaginary parts of both total and regular conductivity. The
frequency grid, broadening, Hamiltonian, current operator, diamagnetic operator,
area, and normalization therefore remain exactly those documented in the
current summary.

Do not rerun the expensive ED benchmark and do not modify any response formula.
The change is a reproducible replot of stored data plus a LaTeX report update.
Shift current remains outside the scope.

## Reproducible Plot Generation

Add a small report-side plotting script and a report-owned copy of the corrected
curve data. The script will validate that the data have nine numeric columns and
then generate four one-page vector PDFs with stable names under
`dmrg/report/figures/optical_response/`.

Each main panel will use:

- exact: thick blue solid line;
- final response Lanczos, `M=892`: thinner red dashed line;
- common horizontal axis `omega` over the stored range;
- a component-specific vertical-axis label;
- an inset showing the pointwise absolute exact-Lanczos difference on a
  logarithmic scale.

The dashed gaps expose the exact blue curve in the main panel. The inset gives a
direct quantitative check that the two calculations agree, instead of relying
only on visual overlap.

## Summary Layout

Replace the current two selected-component subplots in the corrected benchmark
figure with a `2x2` group containing all four final component PDFs. Keep the
existing Lanczos-convergence plot as a separate figure immediately after the
four-component group.

The summary text and captions will state how the exact curve is computed from
the complete finite-sector eigenstate sum, that Lanczos evaluates the same
resolvent, and that the red dashed curve nearly covers the blue exact curve.
The existing benchmark table and conductivity definitions remain unchanged.

## Verification

- A plot-data test must fail when any required column or output is missing, then
  pass after the plotting implementation is added.
- Confirm that four stable component PDFs are generated and each has one page.
- Check that the numeric curves read from the report-owned data agree with the
  source benchmark data exactly.
- Compile `dmrg_summary.tex` with `latexmk` and reject fatal LaTeX errors,
  undefined controls, overfull boxes, and PDF-string warnings.
- Render and visually inspect every page affected by the new figure layout for
  legibility, clipping, overlap, and correct captions.
- Preserve unrelated untracked files in the DMRG worktree.
