# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

`adjustr` is an R package for Bayesian model sensitivity analysis using Pareto-smoothed importance sampling (PSIS). It works with Stan models (rstan, brms, cmdstanr) to estimate how posterior quantities change under alternative prior/likelihood specifications without refitting.

## Development Commands

```bash
# Check package (full R CMD check)
R CMD check .

# Run all tests
Rscript -e 'testthat::test_dir("tests/testthat")'

# Run a single test file
Rscript -e 'testthat::test_file("tests/testthat/test_weights.R")'

# Build documentation (roxygen2)
Rscript -e 'devtools::document()'

# Install locally for testing
R CMD INSTALL .

# Build pkgdown site
Rscript -e 'pkgdown::build_site()'
```

## Architecture

The package follows a three-step pipeline: **specify → weight → summarize**.

### Core Pipeline (in execution order)

1. **`R/make_spec.R`** — `make_spec()` creates `adjustr_spec` objects from formula-based sampling statements and parameter grids. Also implements dplyr verb S3 methods for `adjustr_spec`.

2. **`R/adjust_weights.R`** — `adjust_weights()` is the main workhorse. It:
   - Parses Stan model code via `parse_model()`
   - Computes log-probabilities under original and alternative specs
   - Uses `loo::psis()` for Pareto-smoothed importance sampling
   - Returns an `adjustr_weighted` tibble with `.weights` and `.pareto_k` columns
   - Also contains `extract_samp_stmts()` and `get_fit_obj()` (normalizes stanfit/brmsfit/cmdstanr inputs)

3. **`R/use_weights.R`** — Post-weighting analysis:
   - `summarize.adjustr_weighted()` — S3 method using tidy eval to compute weighted posterior summaries (mean, var, sd, quantile, wasserstein)
   - `spec_plot()` — ggplot2-based visualization of posterior quantities vs specification parameters
   - `get_resampling_idxs()` — importance resampling indices
   - Internal weighted ECDF and Wasserstein distance implementations

### Supporting Modules

- **`R/parsing.R`** — Stan model code parser using regex. `parse_model()` extracts variable declarations (with their block types) and sampling statements from Stan code. `match_sampling_stmts()` matches user-provided specs to parsed statements.

- **`R/logprob.R`** — Log-probability computation engine. Maps Stan distribution names to R density functions (with curried evaluation via `make_dens()`). `calc_lp()` evaluates log-densities using tidy eval. Additional distributions loaded from `extraDistr` in `.onLoad()`.

### Key Design Patterns

- Stan distributions are mapped to R density functions in the `distrs` list (`R/logprob.R`), turned into a curried `distr_env` environment at load time
- Heavy use of rlang tidy evaluation — sampling specs are R formulas, summary expressions are quosures evaluated against MCMC draws
- Two custom S3 classes: `adjustr_spec` (specification) and `adjustr_weighted` (results tibble)
- Tests use a pre-fitted Stan model stored in `tests/test_model.rda`, loaded in `tests/testthat/setup.R`
