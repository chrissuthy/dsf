# dsf – development notes

Running record of decisions and status for the `dsf` package and the
simulation study behind *Destination Selection Function: a competing hazards
model for redistribution data* (Sutherland et al., in prep., for *Ecology*
Statistical Innovations). Commits reference this file.

## Model

- Hazard λ_ij = exp(−d(s_i, v_j)² / 2σ_i²) θ_j; π_ij = λ_ij / Σ_k λ_ik over
  destinations operational in the individual's observation year;
  e_i ~ Categorical(π_i); likelihood is the product over individuals.
- log σ_i = γ₀ + γ′x_i (individual covariates, with intercept);
  log θ_j = β′z_j (destination covariates, **no intercept** – cancels in the
  normalisation); optional log(effort_j) offset.
- Non-Euclidean: least-cost distance with cell cost exp(α·c(x)); α = 0
  recovers Euclidean (up to lattice-path error).
- Identifiability: only terms that vary across destinations within an
  individual are identifiable; individual main effects and the θ intercept
  cancel. σ is a *redistribution scale*, not a pure dispersal distance.
- Relationship: equals the SCR multi-catch likelihood conditioned on capture
  with a known activity centre; conditional-logit/SSF form with the choice set
  defined by the sampling design.
- Dispersal range D = 2σ (~86% of Gaussian dispersal events fall within D).

## Package

- `simulate_dsf()` returns only what a field study records (`ind_df`,
  `trap_df`, `op_mat`), plus the `landscape` raster and a reusable `design`.
- `dsf()` fits by maximum likelihood; returns class `dsf` with `coef`, `se`,
  `logLik`, `AIC`.
- `make_cost_fun()` (internal) builds least-cost distances with igraph on a
  raster aggregated to `cost_res`.
- Key options: `trap_layout = "random" | "grid"`; `jitter` (fraction of the
  lattice spacing, 0 = regular lattice, default); `births = "traps" |
  "uniform"`; `surface = TRUE` only when a Gaussian random field covariate is
  needed; `landscape_seed` fixes the design across replicates while the
  dispersal draw varies; `exclude_natal`; `cost = ~ cov` with
  `beta_cost` / `cost_res`.
- Tests: `tests/testthat/test-roundtrip.R`.
- Workflow after an edit: Document → Test → Install → commit → push.

## Simulation study (scripts live outside the package)

| Script | Purpose |
|---|---|
| `sim00_extent_density.R` | Rayleigh-estimator bias across landscape extent (2×2, 4×4, 6×6 D) × sampling intensity (1, 4, 9 per D²; spacing D, D/2, D/3); 100 individuals, uniform origins, regular lattice, 100 reps; DSF recorded for contrast |
| `fig01_naive_bias.R` | Figure 1 from sim00: relative bias vs intensity, one curve per extent |
| `sim01_null.R` | Null model recovery |
| `sim02_age.R` | Age covariate on σ |
| `sim03_theta.R` | Landscape covariate on θ |
| `sim04_noneuc.R` | Non-Euclidean, `cost = ~ cov`, α = 1 |
| `sim05_naive.R` | Rayleigh vs DSF across trap counts (random traps) |

## Design decisions

- **Densities 1, 4, 9 per D² rather than 1, 3, 5:** a square lattice at
  spacing D/√k on a square of side m·D needs m√k traps per side, an integer
  only for perfect-square k.
- **Grid, not random traps, in sim00:** removes between-replicate variation
  from trap configuration; uncertainty then reflects the dispersal process
  alone. Same lattice per scenario via `landscape_seed`.
- **No jitter:** uniform random origins already break lattice symmetry.
- **Fixed abundance across scenarios (100 individuals):** cells differ in
  bias, not sample size.
- **Natal trap:** sim00 uses uniform origins with no exclusion; sims 1–5 use
  births at traps with `exclude_natal = TRUE`.
- **Choice set = traps operational in the capture year** (`op_mat`).
- **Sibling origins** inflate σ by ~√2 (simulation: 1.43); record pair type.

## Technical notes

- Log-sum-exp (subtract row max) in the likelihood; −Inf mask for
  unavailable traps; natal trap = `which.min` of Euclidean d².
- GRF via `fields::matern.image.cov` / `sim.rf`; `aRange` in km converted to
  cells; embedding size 2^ceiling(log2(max(2·ncell, 4·range_cells))) avoids
  the negative-FFT error.
- `identical(f, ~1)` fails across environments → `is_null_formula()` checks
  `term.labels`.

## Log

- **2026-10-17** – Model appropriateness, SCR/SSF relationship,
  identifiability; naming; paper structure.
- **2026-10-18** – Package skeleton, `simulate_dsf()` / `dsf()`, GitHub repo,
  tests; sims 1–5 scripts; non-Euclidean cost distance; grid layout.
- **2026-10-19** – Methods to open with the extent × intensity bias
  demonstration; sim06 renamed `sim00`; densities 1/4/9; `jitter` argument
  (default 0); `fig01` script; these notes added.
