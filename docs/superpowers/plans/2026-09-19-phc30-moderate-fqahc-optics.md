# 30-site Moderate-Coupling FQAHC Optics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Recompute the translation-consistent 30-site spectrum and saved states at `V1=10,V2=V3=2`, diagnose the 15-state FQAHC candidate manifold, and calculate representative-sector optical conductivity on W003.

**Architecture:** Add a parameterized spectrum driver and a separate result analyzer while reusing the shared ED kernels. Saved JLD2 states remain in a W003 run directory; the analyzer emits small traceable files. Both CLI lattice labels use the translation-consistent physical mesh, while retaining distinct old grouped and new single-sector file layouts.

**Tech Stack:** Julia 1.10, KrylovKit, JLD2, PBS on W003, Python `unittest`, gnuplot, Git.

---

### Task 1: Parameterized translation-consistent spectrum driver

**Files:**
- Create: `ed/cases/case3_30sites_Np12/run_spectrum.jl`
- Create: `ed/tests/load_spectrum_driver.jl`
- Modify: `ed/tests/runtests.jl`

- [ ] **Step 1: Write a failing load test**

The test must include the driver without running the large calculation and assert defaults `Np=12`, `V1=10`, `V2=V3=2`, `lattice=corrected`, plus the `partial_<start>.jld2` naming contract.

- [ ] **Step 2: Run the test on W003 and verify RED**

Run:

```bash
/home/shajy/bin/julia --threads=4 --project=. ed/tests/load_spectrum_driver.jl
```

Expected: failure because `run_spectrum.jl` does not exist.

- [ ] **Step 3: Implement the driver**

The driver must parse `--sector-start`, `--sector-end`, `--output-dir`, `--lattice`, `--t1`, `--t3`, `--V1`, `--V2`, `--V3`, `--nev`, and `--krylovdim`. It must select

```julia
lattice = SPECTRUM_LATTICE == "corrected" ? TiltedLat30() : LegacyTiltedLat30()
```

and save `ev_pairs`, `gs_vecs`, all model parameters, `lattice_convention`, `nev`, and `krylovdim` to `partial_<sector-start>.jld2`.

- [ ] **Step 4: Verify GREEN on W003**

Run the load test and the complete Julia test suite.

- [ ] **Step 5: Commit**

```bash
git commit -m "feat: add parameterized 30-site spectrum driver"
```

### Task 2: Generic 15-state manifold analyzer

**Files:**
- Create: `ed/cases/case3_30sites_Np12/analyze_ground_manifold.jl`
- Create: `ed/tests/load_ground_manifold_driver.jl`
- Create: `ed/tests/moderate_fqahc_analysis_test.jl`
- Modify: `ed/tests/runtests.jl`

- [ ] **Step 1: Write failing tests**

Test parameter consistency checks, sorting of synthetic sector energies, the 15-state width/gap calculation, and driver loading. The classifier API must return `width = E[15]-E[1]`, `next_energy = E[16]-E[1]`, and `separation = E[16]-E[15]` without calling the manifold topological.

- [ ] **Step 2: Verify RED on W003**

Expected: missing analyzer functions.

- [ ] **Step 3: Implement the analyzer**

It must load the three partial files, validate corrected lattice and model parameters, write `data/spectrum.dat`, reconstruct all 15 lowest states sequentially, compute representative-orbit structure factors, and write:

```text
data/structure_factor_ground_manifold.dat
data/manifold_diagnostics.txt
parameters.toml
manifest.toml
```

- [ ] **Step 4: Verify GREEN on W003**

Run focused tests and the full suite.

- [ ] **Step 5: Commit**

```bash
git commit -m "feat: analyze 30-site FQAHC candidate manifold"
```

### Task 3: Preserve old optics while enabling corrected-lattice response

**Files:**
- Modify: `ed/cases/case3_30sites_Np12/run_optical_response.jl`
- Modify: `ed/tests/load_optical_response_driver.jl`
- Modify: `ed/tests/large_sector_optimizations_test.jl`

- [ ] **Step 1: Write a failing convention test**

Assert that `--lattice corrected` selects `TiltedLat30()` and that the default remains `legacy`, preserving the existing `V1=100` result.

- [ ] **Step 2: Verify RED on W003**

Expected: unknown lattice argument or missing selector.

- [ ] **Step 3: Implement explicit lattice selection**

Add `--lattice legacy|corrected`, validate the saved JLD2 `lattice_convention` when present, and record the selected convention in response JLD2/text headers.

- [ ] **Step 4: Verify GREEN and rerun the 4x6 benchmark**

The exact and Lanczos Drude weights must remain equal within the existing tolerance.

- [ ] **Step 5: Commit**

```bash
git commit -m "feat: support corrected lattice in 30-site optics"
```

### Task 4: W003 PBS workflow

**Files:**
- Create: `ed/cases/case3_30sites_Np12/submit/run_spectrum_w003.pbs`
- Create: `ed/cases/case3_30sites_Np12/submit/analyze_moderate_fqahc_w003.pbs`
- Create: `ed/cases/case3_30sites_Np12/submit/run_moderate_optical_w003.pbs`
- Modify: `ed/tests/test_repository_layout.py`

- [ ] **Step 1: Extend the layout test and verify RED**

Require `short`, `select=1:ncpus=24:mem=90gb`, `--threads=24`, corrected lattice, and the dedicated W003 run/result directories.

- [ ] **Step 2: Add the PBS scripts**

Spectrum jobs accept `PHC_SECTOR_START/END`; analyzer consumes all three outputs; optical jobs accept `PHC_SECTOR` and pass `V1=10,V2=V3=2,lattice=corrected`.

- [ ] **Step 3: Verify GREEN**

Run Python tests and `bash -n` for every script.

- [ ] **Step 4: Commit**

```bash
git commit -m "ops: add W003 moderate FQAHC workflow"
```

### Task 5: Run and classify the spectrum on W003

- [ ] **Step 1: Sync the committed source to W003**
- [ ] **Step 2: Run a one-sector preflight and verify residual/output schema**
- [ ] **Step 3: Submit single-sector array jobs for `k=0:14`**
- [ ] **Step 4: Wait for all jobs and inspect convergence logs**
- [ ] **Step 5: Run the manifold analyzer**
- [ ] **Step 6: Check all saved-state residuals, 15-state width/gap, structure-factor peaks, file hashes, and large-file exclusion**

No optical job is submitted unless this step yields an isolated 15-state FQAHC candidate manifold.

### Task 6: Representative optical response and result report

**Files:**
- Create after data exist: `results/phc30_np12_v1_10_v2_2_v3_2_fqahc/*`
- Create or modify: `reports/fqahc_results/plot_phc30_moderate_fqahc.gnuplot`
- Modify: `reports/fqahc_results/dmrg_summary.tex`

- [ ] **Step 1: Submit corrected-lattice optical jobs for `k=5,10,0`**
- [ ] **Step 2: Wait for completion and validate residuals and Lanczos convergence**
- [ ] **Step 3: Compare regular/total real and imaginary components across sectors**
- [ ] **Step 4: If sector dependence is material, submit the remaining 12 sectors and compute the 15-state average**
- [ ] **Step 5: Generate spectrum, structure-factor, and optical figures on W003**
- [ ] **Step 6: Add a Chinese result subsection that states how each figure was calculated and what it shows**
- [ ] **Step 7: Run all source, result-package, Julia, and plotting tests on W003**
- [ ] **Step 8: Commit small results, sync W003/GitHub, and verify three-way SHA equality**
