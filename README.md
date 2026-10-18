# dsf: Destination Selection Functions

A competing-hazards model for redistribution data in which an individual's
origin, its destination, and the full set of candidate destinations are known
(kill traps in an invasive-species removal, nest boxes, translocation release
sites). The likelihood has the form of a conditional logit over the candidate
destinations, derived from the sampling process rather than assumed.

- `simulate_dsf()` generates a landscape, trap array, individuals, and
  redistribution events, returning only what a field programme would record.
- `dsf()` fits the model: covariates on the redistribution scale
  (individual level, `sigma ~`) and on destination attractiveness
  (destination level, `theta ~ 0 + ...`), with a trap-by-year operational
  matrix defining each individual's choice set.

```r
# install.packages("remotes")
remotes::install_github("chrissuthy/dsf")

sim <- simulate_dsf(~ age_c, ~ 0 + cov,
                    beta_sigma = c(log(5), 0.3), beta_theta = 0.8,
                    exclude_natal = TRUE, seed = 1)
fit <- dsf(~ age_c, ~ 0 + cov, sim$ind_df, sim$trap_df, sim$op_mat,
           exclude_natal = TRUE)
fit
```

Companion to Sutherland et al., *Destination Selection Function: a competing
hazards model for redistribution data* (in prep.).
