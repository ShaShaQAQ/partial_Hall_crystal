# Optical Conductivity Note Chapter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a self-contained linear optical-conductivity chapter to the DMRG living note, including the legacy explicit Kubo figures, the corrected exact/response-Lanczos benchmark, and unambiguous definitions of every real and imaginary conductivity component.

**Architecture:** Copy nine lightweight one-page PDF assets into a report-owned figures directory, then add one top-level LaTeX section before `Open Checks`. The section derives the old sum-over-states expression and the current resolvent/Lanczos expression from the same spectral data, separates regular/Drude/total conductivity, and labels legacy plots as historical rather than validated.

**Tech Stack:** LaTeX (`article`, `amsmath`, `booktabs`, `graphicx`, `subcaption`), latexmk/pdfTeX, Poppler (`pdfinfo`, `pdftoppm`), Git.

---

## Workspace and file map

Work only in:

```text
/Users/shajianyu/CMP_manybody/partial_Hall_crystal/.worktrees/dmrg-module
```

The branch is `DMRG`. Preserve these unrelated untracked paths:

```text
.codex_diagnostics/
dmrg/plot_ground_density.py
dmrg/test_plot_ground_density.py
tmp/
```

Files for this task:

- Create: `dmrg/report/figures/optical_response/legacy_optical_conductivity.pdf`
- Create: `dmrg/report/figures/optical_response/legacy_optical_conductivity_eta_0p001_tag_m0p359.pdf`
- Create: `dmrg/report/figures/optical_response/legacy_optical_conductivity_eta_0p065.pdf`
- Create: `dmrg/report/figures/optical_response/legacy_optical_conductivity_eta_0p065_tag_m0p373.pdf`
- Create: `dmrg/report/figures/optical_response/legacy_longitudinal_x.pdf`
- Create: `dmrg/report/figures/optical_response/legacy_longitudinal_y.pdf`
- Create: `dmrg/report/figures/optical_response/corrected_optical_conductivity_regular.pdf`
- Create: `dmrg/report/figures/optical_response/corrected_optical_conductivity_total.pdf`
- Create: `dmrg/report/figures/optical_response/corrected_lanczos_convergence.pdf`
- Modify: `dmrg/report/dmrg_summary.tex`
- Generate but do not track: `dmrg/report/dmrg_summary.pdf` and normal LaTeX auxiliary files, already ignored by `.gitignore`.

### Task 1: Import and verify report-owned figure assets

**Files:**
- Create: the nine files under `dmrg/report/figures/optical_response/` listed above.

- [ ] **Step 1: Verify the report-owned asset check is initially RED**

Run:

```bash
test -f dmrg/report/figures/optical_response/legacy_longitudinal_x.pdf
```

Expected: exit status `1`, because the report-owned asset directory does not yet exist.

- [ ] **Step 2: Copy the six historical linear-response PDFs with stable names**

Run:

```bash
mkdir -p dmrg/report/figures/optical_response
cp '/Users/shajianyu/CMP_manybody/ED/组会/optical_conductivity.pdf' \
  dmrg/report/figures/optical_response/legacy_optical_conductivity.pdf
cp '/Users/shajianyu/CMP_manybody/ED/组会/optical_conductivity_0.001_-0.359.pdf' \
  dmrg/report/figures/optical_response/legacy_optical_conductivity_eta_0p001_tag_m0p359.pdf
cp '/Users/shajianyu/CMP_manybody/ED/组会/optical_conductivity_0.065.pdf' \
  dmrg/report/figures/optical_response/legacy_optical_conductivity_eta_0p065.pdf
cp '/Users/shajianyu/CMP_manybody/ED/组会/optical_conductivity_0.065_-0.373.pdf' \
  dmrg/report/figures/optical_response/legacy_optical_conductivity_eta_0p065_tag_m0p373.pdf
cp '/Users/shajianyu/CMP_manybody/ED/组会/long_x.pdf' \
  dmrg/report/figures/optical_response/legacy_longitudinal_x.pdf
cp '/Users/shajianyu/CMP_manybody/ED/组会/longitudinal_y.pdf' \
  dmrg/report/figures/optical_response/legacy_longitudinal_y.pdf
```

- [ ] **Step 3: Copy the three corrected benchmark PDFs with stable names**

Run:

```bash
cp /Users/shajianyu/CMP_manybody/partial_Hall_crystal/.worktrees/optical-response-lanczos/case4_4x6_Np4/optical_response_output/optical_conductivity_regular.pdf \
  dmrg/report/figures/optical_response/corrected_optical_conductivity_regular.pdf
cp /Users/shajianyu/CMP_manybody/partial_Hall_crystal/.worktrees/optical-response-lanczos/case4_4x6_Np4/optical_response_output/optical_conductivity_total.pdf \
  dmrg/report/figures/optical_response/corrected_optical_conductivity_total.pdf
cp /Users/shajianyu/CMP_manybody/partial_Hall_crystal/.worktrees/optical-response-lanczos/case4_4x6_Np4/optical_response_output/lanczos_convergence.pdf \
  dmrg/report/figures/optical_response/corrected_lanczos_convergence.pdf
```

- [ ] **Step 4: Verify all assets exist and are one-page PDFs**

Run:

```bash
find dmrg/report/figures/optical_response -maxdepth 1 -type f -name '*.pdf' | sort
for figure in dmrg/report/figures/optical_response/*.pdf; do
  test "$(pdfinfo "$figure" | awk '/^Pages:/ {print $2}')" = 1 || exit 1
done
test "$(find dmrg/report/figures/optical_response -maxdepth 1 -type f -name '*.pdf' | wc -l | tr -d ' ')" = 9
```

Expected: nine paths are printed and the command exits `0`.

- [ ] **Step 5: Commit the report-owned assets**

```bash
git add dmrg/report/figures/optical_response
git diff --cached --check
git commit -m "docs: import optical response comparison figures"
```

### Task 2: Add the self-contained optical-conductivity chapter

**Files:**
- Modify: `dmrg/report/dmrg_summary.tex`

- [ ] **Step 1: Verify the chapter check is RED**

Run:

```bash
rg -n '\\section\{Optical Conductivity: Kubo Sum and Response Lanczos\}' \
  dmrg/report/dmrg_summary.tex
```

Expected: exit status `1`; the section is absent.

- [ ] **Step 2: Add the figure-layout package**

Add after `\usepackage{graphicx}`:

```tex
\usepackage{subcaption}
```

- [ ] **Step 3: Insert the complete chapter before `\section{Open Checks}`**

Insert this exact LaTeX block:

```tex
\clearpage
\section{Optical Conductivity: Kubo Sum and Response Lanczos}
\label{sec:optical-conductivity}

This section records the linear optical-response conventions used by the older
explicit sum-over-states code and by the corrected response-Lanczos benchmark.
Only linear conductivity is discussed.  The two numerical methods evaluate the
same many-body resolvent when supplied with the same Hamiltonian and current
operators; their distinction is computational rather than physical.

\subsection{Electromagnetic convention and operators}

For a uniform Cartesian vector potential along (x), define
\begin{equation}
  H(A_x)=H(0)+A_xJ_x+\frac{A_x^2}{2}K_{xx}+O(A_x^3),
  \qquad
  J_x=\left.\frac{\partial H}{\partial A_x}\right|_{0},
  \qquad
  K_{xx}=\left.\frac{\partial^2H}{\partial A_x^2}\right|_{0}.
  \label{eq:optical-operators}
\end{equation}
The corrected finite-torus implementation uses
\begin{equation}
  H_A(k)=H(k+A_x\hat{x})
  \label{eq:bloch-ax}
\end{equation}
as the common source of (H(A_x)), (J_x), and (K_{xx}).  Here (A_x) is a
uniform vector potential per unit Cartesian length, not the total twist.  For
the (4\times6) torus, (\Phi_x=4A_x).  We use (e=\hbar=1) and the operator
sign convention (J_x=\partial H/\partial A_x).

For the ground state (\lvert g\rangle), introduce
\begin{equation}
  Q=1-\lvert g\rangle\langle g\rvert,
  \qquad
  \lvert f\rangle=QJ_x\lvert g\rangle,
  \qquad
  z=\omega+i\eta,
  \label{eq:response-source}
\end{equation}
and, for every excited state that couples to the source,
\begin{equation}
  \Delta_n=E_n-E_g,
  \qquad
  w_n=\left|\langle n\rvert J_x\lvert g\rangle\right|^2.
\end{equation}

\subsection{Explicit sum-over-states Kubo formula}

Full diagonalization gives the finite-system resolvent
\begin{align}
  G(z)
  &=\langle f\rvert\left[z-(H-E_g)\right]^{-1}\lvert f\rangle
    =\sum_{n\ne g}\frac{w_n}{z-\Delta_n}, \\
  \chi(z)
  &=G(z)+G(-z)
    =\sum_{n\ne g}w_n
      \left(\frac{1}{z-\Delta_n}-\frac{1}{z+\Delta_n}\right), \\
  \chi_0
  &\equiv\chi(0)=2\operatorname{Re}G(0)
    =-2\sum_{n\ne g}\frac{w_n}{\Delta_n}.
  \label{eq:exact-resolvent}
\end{align}
The active longitudinal expression in the older ED script can be rearranged
state by state as
\begin{align}
 &\frac{w_n}{iz}
  \left(\frac{1}{z+\Delta_n}-\frac{1}{z-\Delta_n}\right)
  +\frac{2iw_n}{z\Delta_n} \\
 &\hspace{1.4cm}=\frac{iw_n}{z}
  \left[
    \frac{1}{z-\Delta_n}-\frac{1}{z+\Delta_n}
    +\frac{2}{\Delta_n}
  \right].
\end{align}
After summing over (n\), its response therefore has the structure
\begin{equation}
  \sigma_{xx}^{\mathrm{old}}(z)
  =\frac{2\pi i}{A_{\mathrm{old}}}
   \frac{\chi(z)-\chi_0}{z},
  \label{eq:old-kubo}
\end{equation}
which is the regular, finite-frequency conductivity in the definitions below.
It is not the total response including an independently evaluated diamagnetic
term.

\subsection{Regular, Drude, and total conductivity}

With the physical finite-torus area (A), define
\begin{align}
  D_{xx}&=\langle g\rvert K_{xx}\lvert g\rangle+\chi_0, \\
  \sigma_{xx}^{\mathrm{Drude}}(z)
    &=\frac{2\pi i}{A}\frac{D_{xx}}{z}, \\
  \sigma_{xx}^{\mathrm{reg}}(z)
    &=\frac{2\pi i}{A}\frac{\chi(z)-\chi_0}{z}, \\
  \sigma_{xx}^{\mathrm{total}}(z)
    &=\sigma_{xx}^{\mathrm{Drude}}(z)
      +\sigma_{xx}^{\mathrm{reg}}(z)
      =\frac{2\pi i}{A}
       \frac{\langle K_{xx}\rangle+\chi(z)}{z}.
  \label{eq:conductivity-parts}
\end{align}
These are complex conductivities.  To make the origin of every plotted real or
imaginary part explicit, write (X=X'+iX'') and
\begin{equation}
  \sigma_X(z)=\frac{2\pi i}{A}\frac{X}{\omega+i\eta}.
\end{equation}
Then
\begin{align}
  \operatorname{Re}\sigma_X
   &=\frac{2\pi}{A}
     \frac{\eta X'-\omega X''}{\omega^2+\eta^2}, \\
  \operatorname{Im}\sigma_X
   &=\frac{2\pi}{A}
     \frac{\omega X'+\eta X''}{\omega^2+\eta^2}.
  \label{eq:real-imag-source}
\end{align}
For the three conductivities, the corresponding numerators are
\begin{equation}
  X_{\mathrm{Drude}}=D_{xx},\qquad
  X_{\mathrm{reg}}=\chi(z)-\chi_0,\qquad
  X_{\mathrm{total}}=\langle K_{xx}\rangle+\chi(z).
\end{equation}

\begin{table}[h]
\centering
\small
\caption{Meaning and source of the complex longitudinal conductivities.}
\begin{tabular}{p{0.20\linewidth}p{0.34\linewidth}p{0.34\linewidth}}
\toprule
Quantity & Real part & Imaginary part \\
\midrule
(\sigma^{\mathrm{reg}}) & Finite-frequency absorption spectrum, broadened by
(\eta) & Reactive/dispersive regular response; the Kramers--Kronig partner of
the absorption spectrum \\
(\sigma^{\mathrm{Drude}}) & Broadened zero-frequency Drude contribution; it
approaches a delta contribution as (\eta\to0^+) & Reactive Drude pole,
approaching (1/\omega) away from zero frequency \\
(\sigma^{\mathrm{total}}) & Sum of the regular and Drude real parts & Sum of
the regular and Drude imaginary parts \\
\bottomrule
\end{tabular}
\end{table}

The finite-torus value of (D_{xx}) is the curvature combination for the
selected finite-size twist and momentum branch.  A negative value at one fixed
twist must not be read directly as a negative thermodynamic dissipative Drude
weight; that interpretation requires following the global ground state versus
twist and performing a size extrapolation.

The benchmark data file stores all four complex curves, even though the current
PDFs plot only selected components:
\begin{table}[h]
\centering
\small
\caption{Columns of \path{optical_response_curves.dat}.}
\begin{tabular}{cl}
\toprule
Columns & Stored response \\
\midrule
2--3 & Exact (\operatorname{Re/Im}\sigma^{\mathrm{total}}) \\
4--5 & Lanczos (\operatorname{Re/Im}\sigma^{\mathrm{total}}) \\
6--7 & Exact (\operatorname{Re/Im}\sigma^{\mathrm{reg}}) \\
8--9 & Lanczos (\operatorname{Re/Im}\sigma^{\mathrm{reg}}) \\
\bottomrule
\end{tabular}
\end{table}

\subsection{Exact diagonalization versus response Lanczos}

The exact reference materializes and fully diagonalizes the Hamiltonian in the
finite (m=2) sector, then evaluates Eq.~\eqref{eq:exact-resolvent} using every
stored (\Delta_n) and (w_n).  It is exact for that finite sector and chosen
broadening, not for the thermodynamic limit.

Response Lanczos starts from
(q_1=\lvert f\rangle/\lVert f\rVert) and generates
\begin{equation}
  (H-E_g)q_j
  =\beta_{j-1}q_{j-1}+\alpha_jq_j+\beta_jq_{j+1}.
  \label{eq:response-lanczos-recurrence}
\end{equation}
If (T_M) is the resulting tridiagonal matrix, the same resolvent is
approximated by
\begin{align}
  G_M(z)
  &=\lVert f\rVert^2 e_1^{\mathsf T}(zI-T_M)^{-1}e_1 \\
  &=\frac{\lVert f\rVert^2}
  {z-\alpha_1-\dfrac{\beta_1^2}
  {z-\alpha_2-\dfrac{\beta_2^2}{\ddots}}}.
  \label{eq:lanczos-continued-fraction}
\end{align}
Thus exact diagonalization and response Lanczos are two evaluations of the same
Kubo resolvent.  Lanczos avoids storing all excited eigenstates and is the route
intended for larger ED sectors.

\subsection{Historical linear-response plots}

\paragraph{Historical-status warning.}
The following plots are preserved as records of the earlier linear-response
calculation, not as validation data for the corrected operators.  In the old
script the active (x) velocity was (v_1), while the Cartesian combination
was commented out; the velocity construction reused a Hamiltonian apply path
that could add interactions; and the normalization used (L_1L_2) rather than
the physical torus area.  In addition, folding the (d_x=\pm L_x/2) hopping
paths before differentiating produced the wrong (J_x) and (K_{xx}), even
though (H(0)) and its spectrum were unchanged.

\begin{figure}[p]
\centering
\begin{subfigure}{0.48\linewidth}
  \includegraphics[width=\linewidth]{figures/optical_response/legacy_optical_conductivity.pdf}
  \caption{Historical transverse plot with no broadening recorded in its file name.}
\end{subfigure}\hfill
\begin{subfigure}{0.48\linewidth}
  \includegraphics[width=\linewidth]{figures/optical_response/legacy_optical_conductivity_eta_0p001_tag_m0p359.pdf}
  \caption{Historical transverse plot; the original file records (\eta=0.001)
  and the unresolved tag (-0.359).}
\end{subfigure}
\begin{subfigure}{0.48\linewidth}
  \includegraphics[width=\linewidth]{figures/optical_response/legacy_optical_conductivity_eta_0p065.pdf}
  \caption{Historical transverse plot with recorded (\eta=0.065).}
\end{subfigure}\hfill
\begin{subfigure}{0.48\linewidth}
  \includegraphics[width=\linewidth]{figures/optical_response/legacy_optical_conductivity_eta_0p065_tag_m0p373.pdf}
  \caption{Historical general plot; the original file records (\eta=0.065)
  and the unresolved tag (-0.373).}
\end{subfigure}
\caption{Earlier linear optical-conductivity outputs.  Labels and broadening
values are transcribed from the plots and file names; unresolved numeric tags
are intentionally not assigned a physical interpretation.}
\label{fig:legacy-optical-general}
\end{figure}

\begin{figure}[p]
\centering
\begin{subfigure}{0.48\linewidth}
  \includegraphics[width=\linewidth]{figures/optical_response/legacy_longitudinal_x.pdf}
  \caption{Historical longitudinal (x)-direction result.}
\end{subfigure}\hfill
\begin{subfigure}{0.48\linewidth}
  \includegraphics[width=\linewidth]{figures/optical_response/legacy_longitudinal_y.pdf}
  \caption{Historical longitudinal (y)-direction result.}
\end{subfigure}
\caption{Earlier longitudinal Kubo plots.  Their real and imaginary curves are
the two components of the old regular-response expression,
Eq.~\eqref{eq:old-kubo}, under the old operator and area conventions.}
\label{fig:legacy-optical-longitudinal}
\end{figure}

\clearpage
\subsection{Corrected (4\times6) FCI benchmark}

The corrected benchmark uses the (4\times6) torus with (N_p=4), momentum
sector (m=2), (t_1=1), (t_3=0.2), (V_1=1), and (V_2=V_3=0).  The
frequency grid is (0:0.001:10) and (\eta=0.065).

\begin{table}[h]
\centering
\caption{Corrected finite-torus optical benchmark.}
\begin{tabular}{lr}
\toprule
Quantity & Value \\
\midrule
(E_g) & (-9.313575887962) \\
(\langle K_{xx}\rangle) & (4.47923524291) \\
(\lVert QJ_x\lvert g\rangle\rVert^2) & (10.8443081433) \\
(D_{xx}) & (-0.885036768331) \\
Regular peak ((\omega,\operatorname{Re}\sigma_{xx}^{\mathrm{reg}})) &
((4.107,9.34982602717)) \\
Final response dimension & (892) \\
Exact--Lanczos scaled maximum error & (5.5927404437300286\times10^{-14}) \\
\bottomrule
\end{tabular}
\end{table}

\begin{figure}[p]
\centering
\begin{subfigure}{0.48\linewidth}
  \includegraphics[width=\linewidth]{figures/optical_response/corrected_optical_conductivity_regular.pdf}
  \caption{(\operatorname{Re}\sigma_{xx}^{\mathrm{reg}}): exact sum and
  successive response-Lanczos dimensions.}
\end{subfigure}\hfill
\begin{subfigure}{0.48\linewidth}
  \includegraphics[width=\linewidth]{figures/optical_response/corrected_optical_conductivity_total.pdf}
  \caption{(\operatorname{Im}\sigma_{xx}^{\mathrm{total}}): exact and final
  response-Lanczos curves.}
\end{subfigure}
\par\medskip
\begin{subfigure}{0.58\linewidth}
  \includegraphics[width=\linewidth]{figures/optical_response/corrected_lanczos_convergence.pdf}
  \caption{Scaled maximum error of the regular curve versus response dimension.}
\end{subfigure}
\caption{Corrected exact/response-Lanczos benchmark.  The exact curve is drawn
first; at the largest response dimensions it is hidden by the later Lanczos
curve because the two agree to floating-point accuracy.}
\label{fig:corrected-optical-benchmark}
\end{figure}

\clearpage
\subsection{What can and cannot be compared}

\begin{table}[h]
\centering
\small
\caption{Status of the historical and corrected optical calculations.}
\begin{tabular}{p{0.23\linewidth}p{0.33\linewidth}p{0.33\linewidth}}
\toprule
Item & Historical ED plots & Corrected benchmark \\
\midrule
Electromagnetic operator & Derivative of the old folded-bond construction;
active (x) component was (v_1) & Analytic derivatives of
(H(k+A_x\hat{x})) before finite-torus folding \\
Area & (A_{\mathrm{old}}=L_1L_2) & Physical supercell area
(A=12\sqrt{3}) \\
Response shown & Explicit regular Kubo expression; plot-dependent tensor
component & Regular, Drude, and total complex conductivity \\
Numerical method & Stored excited states and explicit state sum & Both complete
ED sum and response Lanczos \\
Validation status & Historical only; peak positions and heights are not a
corrected reference & Operator finite differences, exact sum, and Lanczos
continued fraction are mutually checked \\
\bottomrule
\end{tabular}
\end{table}

The qualitative location of spectral structures may be inspected across the
historical and corrected figures, but absolute peak heights, low-frequency
behavior, and tensor labels are not directly comparable until the historical
calculation is rerun with the corrected operators and physical area.
```

- [ ] **Step 4: Verify the section and required definitions are present**

Run:

```bash
rg -n -F \
  -e '\section{Optical Conductivity: Kubo Sum and Response Lanczos}' \
  -e '\sigma_{xx}^{\mathrm{reg}}' \
  -e 'Response Lanczos starts from' \
  -e '\paragraph{Historical-status warning.}' \
  dmrg/report/dmrg_summary.tex
! rg -ni 'shift[ -]current' dmrg/report/dmrg_summary.tex
```

Expected: the section, regular-conductivity formula, Lanczos text, and warning
are found; the phrase `shift current` is absent because it is out of scope.

- [ ] **Step 5: Compile the report with fatal errors enabled**

Run from `dmrg/report`:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error dmrg_summary.tex
```

Expected: exit `0` and an updated `dmrg_summary.pdf`.

- [ ] **Step 6: Commit the chapter source after compilation succeeds**

```bash
git add dmrg/report/dmrg_summary.tex
git diff --cached --check
git commit -m "docs: add optical conductivity methods and benchmark"
```

### Task 3: Render and visually verify the final report

**Files:**
- Verify: `dmrg/report/dmrg_summary.pdf`
- Modify if visual defects are found: `dmrg/report/dmrg_summary.tex`

- [ ] **Step 1: Inspect LaTeX diagnostics**

Run:

```bash
rg -n 'LaTeX Error|Undefined control sequence|Overfull \\hbox|Overfull \\vbox' \
  dmrg/report/dmrg_summary.log || true
pdfinfo dmrg/report/dmrg_summary.pdf | rg 'Pages|Page size|PDF version'
```

Expected: no LaTeX errors or undefined controls. Any overfull box introduced by
the new chapter must be fixed before completion.

- [ ] **Step 2: Render every final page to PNG**

Run:

```bash
REPORT_RENDER_DIR=$(mktemp -d /tmp/dmrg-optical-report-render.XXXXXX)
pdftoppm -png -r 140 dmrg/report/dmrg_summary.pdf \
  "$REPORT_RENDER_DIR/page"
find "$REPORT_RENDER_DIR" -maxdepth 1 -type f -name 'page-*.png' | sort
```

Expected: one PNG per PDF page.

- [ ] **Step 3: Visually inspect all rendered pages**

Use `view_image` on every PNG and verify:

- equations and continued fraction fit within the text block;
- response-component and comparison tables do not cross margins;
- all nine plots are sharp and their captions are readable;
- no float-only blank page, clipped figure, overlapping caption, or orphaned
  heading is present;
- the pre-existing DMRG pages remain intact.

If a defect is found, edit only `dmrg/report/dmrg_summary.tex`, rerun `latexmk`,
rerender every page, and repeat this check.

- [ ] **Step 4: Run fresh final structural and Git verification**

Run:

```bash
test "$(find dmrg/report/figures/optical_response -maxdepth 1 -type f -name '*.pdf' | wc -l | tr -d ' ')" = 9
rg -n '\\section\{Optical Conductivity: Kubo Sum and Response Lanczos\}' \
  dmrg/report/dmrg_summary.tex
git diff --check
git status --short --branch
```

Expected: nine tracked input figures, one optical section, clean diff checks,
and only the pre-existing unrelated untracked paths outside the committed task.

- [ ] **Step 5: Commit visual-QA fixes if any were required**

If Step 3 required source changes:

```bash
git add dmrg/report/dmrg_summary.tex
git diff --cached --check
git commit -m "docs: refine optical response report layout"
```

If no source changes were required after the Task 2 commit, do not create an
empty commit.
