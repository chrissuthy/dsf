# dsf – development notes

Running record of decisions and status for the `dsf` package and the
simulation study behind *Destination Selection Function: a competing hazards
model for redistribution data* (Sutherland et al., in prep., for *Ecology*
Statistical Innovations). Commits reference this file.

## Model

- Hazard λ_ij = exp(−d(s_i, v_j)² / 2σ_i²) θ_{j,t_i}; π_ij = λ_ij / Σ_k λ_ik
  over destinations operational in the individual's observation year t_i;
  e_i ~ Categorical(π_i); likelihood is the product over individuals.
- log σ_i = γ₀ + γ′x_i (individual covariates, with intercept);
  log θ_jt = β′z_jt + log(effort_jt) (destination covariates, time-invariant
  or varying by year, **no intercept**).
- Non-Euclidean: least-cost distance with cell cost exp(α·c(x)); α = 0
  recovers Euclidean (up to lattice-path error).
- Identifiability: the normalisation is within individual and each
  individual belongs to one year, so anything constant across destinations
  within a year cancels: the θ intercept, individual main effects, and
  **year main effects on θ**. Identifiable θ terms are those that vary
  across destinations within a year (effort_t, bait_t, habitat). Year may
  enter σ as an individual covariate. σ is a *redistribution scale*.
- Relationship: equals the SCR multi-catch likelihood conditioned on capture
  with a known activity centre; conditional-logit/SSF form with the choice set
  defined by the sampling design.
- Dispersal range D = 2σ (~86% of Gaussian dispersal events fall within D).

## Data structure (settled 2026-10-21)

Four field objects, assembled by `dsf_frame()`:

| object | shape | content |
|---|---|---|
| `ind` | n rows | `src_x`, `src_y`, `dest` (row of `traps`), `obs_year`, individual covariates |
| `traps` | J rows | `trap_id`, `x`, `y`, time-invariant destination covariates (optional `effort`) |
| `op` | J × T, colnames = years | 1 where operational |
| `trap_tv` | J × T × p array (dimnames: years, covariate names) or long data frame (`trap`, `year`, …) | time-varying destination covariates (optional `effort`) |

`dsf_frame()` returns: year index per individual, n × J squared distances,
n × J availability (op in year t_i, optional natal exclusion), observed
index, a long destination-year table `Z` (J·T rows, destination within year;
time-invariant columns recycled, time-varying columns filled in), and a J × T
log-effort offset. NAs at non-operational destination-years are set to 0
(never used). Checks: dest indexes traps; observed destination operational in
its year; no empty choice sets; effort > 0 where operational.

`dsf()` accepts a frame or the raw objects. θ design matrix is built on `Z`
((J·T) × q); the n × J log θ matrix is `t(matrix(Z β, J, T)[, t_i])`.
`theta_design()` refuses an intercept and any column constant across
destinations within every year. `simulate_dsf()` draws from the same frame.

## Package

- `simulate_dsf()` returns only what a field study records (`ind_df`,
  `trap_df`, `op_mat`, `trap_tv`), plus `landscape`, `design` and `frame`.
- `dsf()` fits by maximum likelihood; returns class `dsf` with `coef`, `se`,
  `logLik`, `AIC`, `frame`.
- `dsf_frame()` assembles and checks the data (above).
- `make_cost_fun()` (internal) builds least-cost distances with igraph on a
  raster aggregated to `cost_res`.
- Simulator options: `trap_layout = "random" | "grid"`; `jitter` (fraction of
  lattice spacing, 0 = regular, default); `births = "traps" | "uniform"`;
  `surface = TRUE` only when a Gaussian random field is needed; `trap_cov`
  (vector, data frame, or function of ntraps) for user-supplied destination
  covariates; `trap_tv` (array or function of (ntraps, years)) for
  time-varying ones; `landscape_seed` fixes the design across replicates
  (functions for `trap_cov`/`trap_tv` are evaluated inside its scope);
  `exclude_natal`; `cost = ~ cov` with `beta_cost` / `cost_res`.
- Tests: `tests/testthat/test-roundtrip.R` (13).
- Citation: `inst/CITATION`.
- Workflow after an edit: Document → Test → Install → commit → push.

## Simulation study (scripts live outside the package)

| Script | Purpose |
|---|---|
| `sim00_extent_density.R` | Paper Sim 1: Rayleigh-estimator bias across extent (2×2, 4×4, 6×6 D) × intensity (1, 4, 9 per D²); 100 individuals, uniform origins, regular lattice, 100 reps |
| `fig01_naive_bias.R` | Figure 1 from sim00 |
| `sim02_age.R` | Paper Sim 2.1: age covariate on σ (γ₀ = log 5, γ₁ = 0.3); null vs age DSF; n = 50 and 500 |
| `sim03_theta.R` | Paper Sim 2.2: z_j ~ U(−2, 2) on θ, β = 0.5; null vs hazard DSF; n = 50 and 500 |
| `sim04_noneuc.R` | Paper Sim 2.3: resistance surface exp(α c(x)), α = 1, 200 × 200 raster (500 m), exponential GRF range 10 km; Euclidean vs non-Euclidean DSF; 50 reps (parallel over replicates) |
| `sim01_null.R`, `sim05_naive.R` | earlier null and trap-count studies; superseded by sim00/sim02 |

Shared Sim 2 design: 100 × 100 km, 250 random traps fixed across years and
replicates, 5 cohorts of 10 or 100 born at uniform locations, ages 1/2/3 with
0.5/0.3/0.2, no natal exclusion.

## Results so far

- Sim 1: Rayleigh σ̂ biased −17/−13/−13% (2×2 D), −8/−7/−6% (4×4), −5/−4/−4%
  (6×6) across intensities 1/4/9 — edge truncation, governed by extent not
  intensity; DSF unbiased throughout.
- Sim 2.1: null DSF σ̂ ≈ 6.5 km (+35% at age 1, −27% at age 3), narrow
  intervals; age DSF unbiased at both n.
- Sim 2.2: β̂ unbiased at n = 500 (+6% at n = 50, wide); σ̂ unbiased under
  both models — θ heterogeneity unrelated to distance does not bias σ.
- Sim 2.3: Euclidean σ̂ ≈ +90% (cost-scale vs km; model attributes ease of
  passage to distance); non-Euclidean DSF recovers σ and α.

## Design decisions

- Densities 1, 4, 9 per D² (perfect squares) so the lattice is square on
  every landscape; grid not random traps; no jitter; fixed abundance.
- Grid lattice spacing D/√k; 2×2 D landscape is "end to end".
- Natal trap: sims with uniform origins use no exclusion.
- Choice set = traps operational in the capture year; capture-timing vs age
  confounding is a Discussion point.
- Sibling origins inflate σ by ~√2; record pair type.
- Euclidean and non-Euclidean σ are on different scales; report the
  discrepancy as misattribution, not simple bias.

## Technical notes

- Log-sum-exp in the likelihood; −Inf mask for unavailable traps; natal trap
  = `which.min` of Euclidean d².
- GRF via `fields::matern.image.cov` / `sim.rf`; `aRange` in km converted to
  cells; embedding size 2^ceiling(log2(max(2·ncell, 4·range_cells))).
- `identical(f, ~1)` fails across environments → `is_null_formula()`.
- Non-Euclidean fits: ~1.5 min (n = 50) to ~10 min (n = 500) on a 40 000-node
  graph; parallelise over replicates with `parallel::parLapply`, seeding per
  replicate.

## Log

- **2026-10-17** – Model, SCR/SSF relationship, identifiability; naming; paper
  structure.
- **2026-10-18** – Package skeleton, `simulate_dsf()` / `dsf()`, GitHub repo,
  tests; sims 1–5; non-Euclidean cost distance; grid layout.
- **2026-10-19** – Methods open with the extent × intensity demonstration;
  sim06 → `sim00`; densities 1/4/9; `jitter`; `trap_cov`; sims 2–4 reworked
  to the shared design with two cohort sizes; `inst/CITATION`; Figures 1–3
  and Table 1 drafted.
- **2026-10-21** – Data structure settled: `dsf_frame()` with four field
  objects; θ_jt (time-varying destination covariates) threaded through
  `dsf()` and `simulate_dsf()` (`trap_tv`); year main effects on θ refused;
  13 tests.
