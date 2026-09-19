# Restore Moderate-Coupling FQAHC Classification Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Restore the `V1=10, V2=V3=2` result to an FQAHC candidate with an isolated `3 x 5` quasi-degenerate manifold, while keeping the distinction between established CDW order and unconfirmed topological order explicit.

**Architecture:** Treat the existing design spec and packaged numerical data as the source of truth. Update report regression tests first, then align the Chinese note and spectrum plot labels with those assertions; regenerate only the affected spectrum figure and compile the PDF on W003.

**Tech Stack:** Python `unittest`, gnuplot through the existing Python plotting script, LaTeX/Tectonic, Poppler, Git.

---

### Task 1: Define the corrected report language

**Files:**
- Modify: `reports/fqahc_results/test_report_sources.py`
- Modify: `reports/fqahc_results/test_plot_phc30_moderate_fqahc.py`

- [ ] **Step 1: Write failing assertions**

Require the report to contain `温和参数的 30-site FQAHC 候选态`, `3\\times5`, `隔离的 15 态准简并低能流形`, and the statement that the fivefold topological part still requires flux insertion or a many-body Chern number. Require the spectrum plot source to use `3x5 候选低能流形放大`. Reject the plain-CDW-only conclusions introduced by commit `60a4b72`.

- [ ] **Step 2: Verify RED on W003**

Run:

```bash
python3 -m unittest \
  reports.fqahc_results.test_report_sources \
  reports.fqahc_results.test_plot_phc30_moderate_fqahc
```

Expected: assertion failures because the note and plot still classify the moderate point as plain CDW.

### Task 2: Restore the physical classification

**Files:**
- Modify: `reports/fqahc_results/dmrg_summary.tex`
- Modify: `reports/fqahc_results/plot_phc30_moderate_fqahc.py`
- Regenerate: `reports/fqahc_results/figures/optical_response/phc30_moderate_manybody_spectrum.png`

- [ ] **Step 1: Update the note**

Describe the strong-coupling `100,0,0` spectrum as an isolated three-state CDW ground manifold. Describe the moderate `10,2,2` spectrum as an isolated 15-state quasi-degenerate manifold spanning all momentum sectors, organized into the five CDW momentum triplets

```text
(0,5,10), (1,6,11), (2,7,12), (3,8,13), (4,9,14).
```

State that this is compatible with `3 x 5`, that the structure factor establishes the period-three factor, and that the fivefold topological factor remains a candidate until flux insertion, charge pumping, or a many-body Chern number is computed.

- [ ] **Step 2: Update and regenerate the spectrum figure on W003**

Change the zoom-panel title to `3x5 候选低能流形放大`, run the existing plotting script against the packaged result data, and verify that only the intended spectrum image changes.

- [ ] **Step 3: Verify GREEN on W003**

Run the two targeted `unittest` modules again. Expected: all targeted tests pass.

### Task 3: Verify and publish the corrected note

**Files:**
- Build output: `reports/fqahc_results/dmrg_summary.pdf` (ignored generated artifact)

- [ ] **Step 1: Run complete tests on W003**

```bash
python3 -m unittest discover -s reports/fqahc_results -p 'test_*.py'
python3 -m unittest discover -s ed/tests -p 'test_*.py'
git diff --check
```

Expected: 15 report tests and 15 ED tests pass with no whitespace errors.

- [ ] **Step 2: Compile and visually inspect the PDF on W003**

Compile with Tectonic into the run result directory, confirm 18 pages and zero overfull boxes, render pages 15-18 with Poppler, and check the spectrum, structure factor, four optical panels, captions, and comparison text for clipping or overlap.

- [ ] **Step 3: Commit and synchronize**

Commit the regression tests, note, plotting source, regenerated spectrum image, and this plan. Transfer the commit to W003 without touching untracked PBS logs, push `integration/ed-excitations` from W003, and verify that local, W003, and GitHub resolve to the same SHA.
