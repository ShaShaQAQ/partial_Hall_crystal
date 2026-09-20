# 30-site Np=13 Excitation Workflow Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compute and document the full-ED spectrum, saved sector ground states, static structure factor, and response-Lanczos optical conductivity for the two 30-site models at `Np=13`.

**Architecture:** Parameterize the already-validated 30-site `Np=12` spectrum, manifold-analysis, optical, packaging, and plotting pipeline without changing the shared Hamiltonian or response kernels. Run one W003 PBS job per momentum sector, determine each actual low-energy manifold from the completed spectrum, then analyze and average only those states while retaining all fifteen saved sector ground states remotely.

**Tech Stack:** Julia/JLD2/KrylovKit, Python/NumPy/unittest, PBS/Torque on W003, gnuplot, XeLaTeX/latexmk, Git/GitHub.

---

### Task 1: Parameterize particle number in the spectrum and optical drivers

**Files:**
- Modify: `ed/tests/load_spectrum_driver.jl`
- Modify: `ed/tests/load_optical_response_driver.jl`
- Modify: `ed/cases/case3_30sites_Np12/run_spectrum.jl`
- Modify: `ed/cases/case3_30sites_Np12/run_optical_response.jl`

- [ ] **Step 1: Write failing driver tests**

Add tests that run the driver in a fresh Julia process with `--Np 13`, assert that
`SPECTRUM_NP == 13` and `RESPONSE_NP == 13`, and save a synthetic response whose
model tuple has `Np=13`:

```julia
model = (Np=13, t1=1.0, t3=0.2, V1=10.0, V2=2.0, V3=2.0)
save_response_curve(...; model_parameters=model, requested_mmax=600)
@test load(output_path)["Np"] == 13
```

- [ ] **Step 2: Run the focused tests on W003 and verify RED**

Run:

```bash
/home/shajy/bin/julia --threads=1 --project=. ed/tests/load_spectrum_driver.jl
/home/shajy/bin/julia --threads=1 --project=. ed/tests/load_optical_response_driver.jl
```

Expected: failure because the spectrum driver rejects `Np=13` and the optical
driver has no `RESPONSE_NP` parameter.

- [ ] **Step 3: Implement minimal particle-number propagation**

In `run_spectrum.jl`, replace the fixed `Np == 12` assertion with
`1 <= SPECTRUM_NP < 30`. In `run_optical_response.jl`, parse `--Np`, require the
saved JLD2 `Np` to equal it, and replace both hard-coded basis/model values:

```julia
const RESPONSE_NP = parse(Int, parse_response_argument("--Np", "12"))
Int(data["Np"]) == RESPONSE_NP || error("saved particle number does not match --Np")
basis = gen_basis(lattice.Ns, RESPONSE_NP)
model_parameters=(Np=RESPONSE_NP, ...)
```

- [ ] **Step 4: Run the focused tests on W003 and verify GREEN**

Run both commands from Step 2. Expected: all test sets pass.

- [ ] **Step 5: Commit**

```bash
git add ed/tests/load_spectrum_driver.jl ed/tests/load_optical_response_driver.jl \
  ed/cases/case3_30sites_Np12/run_spectrum.jl \
  ed/cases/case3_30sites_Np12/run_optical_response.jl
git commit -m "feat: parameterize 30-site particle number"
```

### Task 2: Generalize ground-manifold and structure-factor analysis

**Files:**
- Modify: `ed/tests/moderate_fqahc_analysis_test.jl`
- Modify: `ed/tests/load_ground_manifold_driver.jl`
- Modify: `ed/cases/case3_30sites_Np12/analyze_ground_manifold.jl`

- [ ] **Step 1: Write failing metadata and manifold tests**

Create a reusable expected-model tuple and test validation for `Np=13`, both
interaction sets, an arbitrary result ID/phase label, and non-fifteen manifold
sizes. Include rejection of a synthetic partial carrying `Np=12` when `Np=13`
was requested.

```julia
expected = (Np=13, t1=1.0, t3=0.2, V1=100.0, V2=0.0, V3=0.0,
            lattice_convention="corrected")
validated = validate_manifold_partial(metadata, "synthetic.jld2", expected)
@test validated.Np == 13
```

- [ ] **Step 2: Run the focused tests on W003 and verify RED**

Run:

```bash
/home/shajy/bin/julia --threads=1 --project=. ed/tests/moderate_fqahc_analysis_test.jl
/home/shajy/bin/julia --threads=1 --project=. ed/tests/load_ground_manifold_driver.jl
```

Expected: method/constant failures because the analyzer is fixed to the old
moderate `Np=12` result.

- [ ] **Step 3: Parameterize the analyzer**

Add CLI arguments `--Np`, `--t1`, `--t3`, `--V1`, `--V2`, `--V3`,
`--result-id`, `--phase-label`, and `--manifold-size`. Propagate them through
partial validation, basis construction, structure-factor particle number,
headers, diagnostics, `parameters.toml`, and `manifest.toml`. Keep
`classify_candidate_manifold(levels, size)` explicit so the size is chosen only
after inspecting the completed spectrum.

- [ ] **Step 4: Run the focused tests on W003 and verify GREEN**

Run both commands from Step 2. Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add ed/tests/moderate_fqahc_analysis_test.jl \
  ed/tests/load_ground_manifold_driver.jl \
  ed/cases/case3_30sites_Np12/analyze_ground_manifold.jl
git commit -m "feat: generalize 30-site manifold analysis"
```

### Task 3: Generalize optical packaging and plotting

**Files:**
- Modify: `reports/fqahc_results/test_package_phc30_moderate_optics.py`
- Modify: `reports/fqahc_results/test_plot_phc30_moderate_fqahc.py`
- Modify: `reports/fqahc_results/package_phc30_moderate_optics.py`
- Modify: `reports/fqahc_results/plot_phc30_moderate_fqahc.py`

- [ ] **Step 1: Write failing arbitrary-manifold tests**

Use synthetic `Np=13` files for sectors `(2, 7, 12)` and assert that packaging
copies exactly three text curves, records `Np=13` and the supplied result ID,
and creates a three-state average. Test plotting with configurable
`manifold_size`, title, and filename prefix.

- [ ] **Step 2: Run Python tests on W003 and verify RED**

Run:

```bash
python3 -m unittest \
  reports.fqahc_results.test_package_phc30_moderate_optics \
  reports.fqahc_results.test_plot_phc30_moderate_fqahc
```

Expected: failures from fixed sectors `0:14`, `Np=12`, result ID, labels, and
filenames.

- [ ] **Step 3: Implement configurable packaging and figures**

Make sector lists, expected metadata, result ID, manifold size, figure prefix,
and title explicit keyword/CLI inputs while retaining the existing `Np=12`
defaults. Average over `len(sectors)` and write that number into the data header.
Require spectrum files to retain eight levels in each of all 15 momentum
sectors, but allow structure-factor and optical files to contain the selected
manifold sectors only.

- [ ] **Step 4: Run the focused tests on W003 and verify GREEN**

Run the command from Step 2. Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add reports/fqahc_results/test_package_phc30_moderate_optics.py \
  reports/fqahc_results/test_plot_phc30_moderate_fqahc.py \
  reports/fqahc_results/package_phc30_moderate_optics.py \
  reports/fqahc_results/plot_phc30_moderate_fqahc.py
git commit -m "feat: package arbitrary 30-site manifolds"
```

### Task 4: Add W003 Np=13 PBS workflows

**Files:**
- Modify: `ed/tests/test_repository_layout.py`
- Create: `ed/cases/case3_30sites_Np12/submit/run_np13_spectrum_w003.pbs`
- Create: `ed/cases/case3_30sites_Np12/submit/run_np13_analyze_w003.pbs`
- Create: `ed/cases/case3_30sites_Np12/submit/run_np13_optical_w003.pbs`

- [ ] **Step 1: Write failing repository-layout tests**

Require all three jobs to use `short`, `select=1:ncpus=24:mem=90gb`, corrected
lattice, `--Np 13`, isolated `phc30_np13_*` run directories, and environment
variables for the interaction parameters/result ID. Require the spectrum and
optical jobs to map one PBS array element to one sector.

- [ ] **Step 2: Run the test on W003 and verify RED**

Run:

```bash
python3 -m unittest ed.tests.test_repository_layout
```

Expected: failure because the three PBS files do not exist.

- [ ] **Step 3: Add the parameterized jobs**

The spectrum job invokes `run_spectrum.jl --Np 13 --nev 8 --krylovdim 60`.
The analysis job requires an explicit `PHC_MANIFOLD_SIZE`. The optical job maps
`PBS_ARRAY_INDEX` through a comma-separated `PHC_SECTORS` list and invokes
`run_optical_response.jl --Np 13 --mmax 600 --eta 0.065`.

- [ ] **Step 4: Run the test on W003 and verify GREEN**

Run the command from Step 2. Expected: all repository-layout tests pass.

- [ ] **Step 5: Commit**

```bash
git add ed/tests/test_repository_layout.py \
  ed/cases/case3_30sites_Np12/submit/run_np13_*_w003.pbs
git commit -m "ops: add W003 Np13 excitation jobs"
```

### Task 5: Run both Np=13 spectra and identify the physical manifolds

**Files:**
- Create remotely: `/home/public/shajy/codex_runs/phc30_np13_v1_100_v2_0_v3_0/spectrum/partial_*.jld2`
- Create remotely: `/home/public/shajy/codex_runs/phc30_np13_v1_10_v2_2_v3_2/spectrum/partial_*.jld2`

- [ ] **Step 1: Push the feature branch and synchronize W003**

Push via the authenticated GitHub remote, fetch it on W003, and verify local and
W003 `git rev-parse HEAD` agree without deleting the existing PBS logs.

- [ ] **Step 2: Submit two 15-sector arrays**

Use `qsub -J 0-14` with the appropriate `PHC_RESULT_ID`, `PHC_V1`, `PHC_V2`, and
`PHC_V3` for the strong and moderate points.

- [ ] **Step 3: Monitor through completion**

Use `qstat` and the PBS logs. A sector is complete only if its log reaches
`saved:` and `done:` and the JLD2 file exists.

- [ ] **Step 4: Verify all raw states**

Check exactly 15 partial files per result, all sector labels `0:14`, `Np=13`,
matching parameters and corrected-lattice fingerprint, normalized saved states,
and residuals below `1e-7`.

- [ ] **Step 5: Determine each manifold size from the spectrum**

Sort all 120 levels, tabulate adjacent gaps and cumulative widths, and select a
low-energy group only where the gap above it exceeds its internal width. Record
the selected size and sector list before submitting analysis or optics.

### Task 6: Compute structure factors and optical responses

**Files:**
- Create: `results/phc30_np13_v1_100_v2_0_v3_0/**`
- Create: `results/phc30_np13_v1_10_v2_2_v3_2/**`
- Create remotely: `/home/public/shajy/codex_runs/<result_id>/optical/**`

- [ ] **Step 1: Run one analysis job per parameter point**

Pass the measured manifold size and a conservative phase label. Verify spectrum,
structure-factor table, diagnostics, parameters, and SHA256 manifest are written.

- [ ] **Step 2: Submit optical arrays over the measured sector lists**

Use one 24-core/90-GB W003 node per sector and preserve the same frequency mesh,
broadening, and response-Lanczos depth as the `Np=12` benchmark.

- [ ] **Step 3: Verify optical outputs**

For every selected sector require both JLD2 and DAT, `Np=13`, matching model
parameters, state residual below `1e-7`, finite seven-column curves, and a
strictly increasing frequency grid.

- [ ] **Step 4: Package small results and generate figures**

Create the optical manifest/diagnostics and six figures per parameter point:
spectrum, structure factor, regular real/imaginary, and total real/imaginary.

- [ ] **Step 5: Commit the traceable small result packages and figures**

Do not add large JLD2 files. Commit TOML manifests, text data, diagnostics, and
PNG figures only.

### Task 7: Update and compile the Chinese result note

**Files:**
- Modify: `reports/fqahc_results/test_report_sources.py`
- Modify: `reports/fqahc_results/dmrg_summary.tex`
- Modify: `reports/fqahc_results/dmrg_summary.pdf`

- [ ] **Step 1: Write failing report tests**

Require two new `Np=13` subsections, both complete parameter sets, measured
manifold sizes/gaps/sectors, explicit ED/response-Lanczos method wording, and all
twelve new figure references.

- [ ] **Step 2: Run report tests on W003 and verify RED**

Run:

```bash
python3 -m unittest reports.fqahc_results.test_report_sources
```

Expected: failure because the `Np=13` result sections are absent.

- [ ] **Step 3: Add result-centered Chinese sections**

For each parameter point state exactly how each plot was computed and what it
shows. Report the numerical low-energy width/gap, dominant structure-factor
momenta, optical peak positions, Drude statistics, and avoid topology claims not
supported by this calculation.

- [ ] **Step 4: Compile and visually verify the PDF on W003**

Run:

```bash
latexmk -g -xelatex -interaction=nonstopmode -halt-on-error dmrg_summary.tex
pdfinfo dmrg_summary.pdf
```

Render all pages with Poppler, inspect for missing/cropped figures, and require
zero overfull boxes.

- [ ] **Step 5: Run report tests and commit**

Expected: all report tests pass and the tracked PDF is rebuilt.

### Task 8: Full verification and synchronization

**Files:**
- Verify all files changed by Tasks 1-7

- [ ] **Step 1: Run the complete W003 test suite**

```bash
/home/shajy/bin/julia --threads=4 --project=. ed/tests/runtests.jl
python3 -m unittest discover -s ed/tests -p 'test_*.py'
python3 -m unittest discover -s reports/fqahc_results -p 'test_*.py'
git diff --check
```

- [ ] **Step 2: Audit result traceability**

Recompute remote SHA256 values, compare every manifest entry, confirm no JLD2 is
tracked, and confirm the two result packages and twelve new figures are present.

- [ ] **Step 3: Merge the feature branch into `integration/ed-excitations`**

Use a non-destructive merge, rerun the focused tests after merge, and push the
integration branch to GitHub.

- [ ] **Step 4: Synchronize W003 and the canonical local checkout**

Verify local, W003, and GitHub all resolve to the same final SHA while preserving
W003's pre-existing untracked PBS logs.

