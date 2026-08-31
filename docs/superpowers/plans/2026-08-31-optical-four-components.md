# Four-Component Optical-Conductivity Figures Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate and document all four corrected longitudinal regular/total optical-conductivity components with visibly distinguishable exact and final response-Lanczos curves.

**Architecture:** Keep the expensive ED result fixed, import its nine-column curve table into the report, and add one focused Julia plotting script plus one focused test. The script validates the stored table and creates four vector PDFs; LaTeX then displays them as a `2x2` component group while retaining the existing convergence plot separately.

**Tech Stack:** Julia (`DelimitedFiles`, `Plots`, `Test`), LaTeX (`graphicx`, `subcaption`), `latexmk`, Poppler, Git.

---

## File map

- Create `dmrg/report/data/optical_response/optical_response_curves.dat`: report-owned corrected exact/Lanczos curves.
- Create `dmrg/report/plot_optical_components.jl`: data validation and four-component PDF generation.
- Create `dmrg/report/test_plot_optical_components.jl`: output and column-mapping regression tests.
- Create four PDFs under `dmrg/report/figures/optical_response/`: the regular/total real/imaginary component figures.
- Modify `dmrg/report/dmrg_summary.tex`: explain exact visibility and display the four figures plus convergence.

### Task 1: Test and implement four-component plotting

**Files:**
- Create: `dmrg/report/test_plot_optical_components.jl`
- Create: `dmrg/report/plot_optical_components.jl`
- Create: `dmrg/report/data/optical_response/optical_response_curves.dat`

- [ ] **Step 1: Write the failing test**

Create a Julia test that includes `plot_optical_components.jl`, constructs a
small valid nine-column table, verifies the four column mappings, generates the
four expected PDF names in a temporary directory, and verifies that malformed
eight-column input throws `ArgumentError`.

- [ ] **Step 2: Run the test and verify RED**

Run:

```bash
julia --startup-file=no dmrg/report/test_plot_optical_components.jl
```

Expected: failure because `plot_optical_components.jl` does not exist.

- [ ] **Step 3: Implement the minimal plotting module**

Implement:

```julia
load_optical_curves(path)::Matrix{Float64}
component_specs()::Vector
generate_optical_component_plots(data_path, output_dir; lanczos_m=892)
```

Require exactly nine numeric columns and at least two frequency samples. Map
columns `(2,4)`, `(3,5)`, `(6,8)`, and `(7,9)` to total real, total imaginary,
regular real, and regular imaginary respectively. Draw exact as a thick blue
solid line and Lanczos as a thinner red dashed line. Add a log-scale inset of
`max(abs(lanczos-exact), 1e-18)` to each figure.

- [ ] **Step 4: Import the fixed source curve data**

Copy the existing benchmark table from:

```text
/Users/shajianyu/CMP_manybody/partial_Hall_crystal/.worktrees/optical-response-lanczos/case4_4x6_Np4/optical_response_output/optical_response_curves.dat
```

to:

```text
dmrg/report/data/optical_response/optical_response_curves.dat
```

- [ ] **Step 5: Run the test and verify GREEN**

Run:

```bash
julia --startup-file=no dmrg/report/test_plot_optical_components.jl
```

Expected: all tests pass and four temporary PDFs are nonempty.

- [ ] **Step 6: Generate the report figures**

Run:

```bash
julia --startup-file=no dmrg/report/plot_optical_components.jl
```

Expected outputs:

```text
corrected_sigma_xx_regular_real.pdf
corrected_sigma_xx_regular_imag.pdf
corrected_sigma_xx_total_real.pdf
corrected_sigma_xx_total_imag.pdf
```

- [ ] **Step 7: Verify data identity and PDF structure**

Use `cmp` to require byte-identical report/source DAT files. Use `pdfinfo` to
require one page for each of the four new PDFs, then render them with
`pdftoppm` and inspect exact/Lanczos styles, insets, axes, legends, and clipping.

- [ ] **Step 8: Commit plotting code, test, data, and figures**

```bash
git add dmrg/report/data/optical_response \
  dmrg/report/plot_optical_components.jl \
  dmrg/report/test_plot_optical_components.jl \
  dmrg/report/figures/optical_response/corrected_sigma_xx_*.pdf
git diff --cached --check
git commit -m "docs: plot four optical conductivity components"
```

### Task 2: Display four components in the summary

**Files:**
- Modify: `dmrg/report/dmrg_summary.tex`

- [ ] **Step 1: Write a failing document-content check**

Run `rg` for all four new PDF names in `dmrg_summary.tex` and verify at least one
name is absent before editing.

- [ ] **Step 2: Update corrected-benchmark explanation and figures**

State that the blue exact curve is the complete finite-sector eigenstate sum
and is covered by the nearly identical later Lanczos curve. Replace the current
two selected component subfigures with all four new PDFs in a `2x2` group. Use
captions for `Re/Im sigma_xx^reg` and `Re/Im sigma_xx^total`, including the
pointwise absolute-error inset. Move the existing convergence PDF to a separate
figure immediately after the four-component group.

- [ ] **Step 3: Verify the four LaTeX references**

Run:

```bash
for name in regular_real regular_imag total_real total_imag; do
  rg -n "corrected_sigma_xx_${name}.pdf" dmrg/report/dmrg_summary.tex
done
```

Expected: all four commands find one figure reference.

- [ ] **Step 4: Compile and inspect the updated report**

Run `latexmk -g -pdf -interaction=nonstopmode -halt-on-error dmrg_summary.tex`
from `dmrg/report`. Reject LaTeX errors, undefined controls, overfull boxes, and
PDF-string warnings. Render every report page and inspect the corrected
benchmark pages for readable panels, captions, insets, and page flow.

- [ ] **Step 5: Commit the summary update**

```bash
git add dmrg/report/dmrg_summary.tex
git diff --cached --check
git commit -m "docs: show all corrected optical components"
```

### Task 3: Final verification

**Files:**
- Verify all files from Tasks 1-2.

- [ ] **Step 1: Run fresh automated checks**

Run the Julia test, regenerate the four PDFs, compile LaTeX from scratch, confirm
the four PDFs are one page each, confirm all four are referenced once, and run
`git diff --check`.

- [ ] **Step 2: Run final PDF visual QA**

Render the final summary at 140 DPI and inspect every page, with particular
attention to the four-component and convergence pages. Confirm no clipping,
overlap, missing glyphs, empty panels, or illegible inset labels.

- [ ] **Step 3: Record repository state**

Run `git status --short --branch` and `git log -5 --oneline`. Preserve all
unrelated untracked files and report the two feature commits plus the design and
plan commits.
