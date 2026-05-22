# Compute Pareto-smoothed Importance Weights for Alternative Model Specifications

Given a set of new sampling statements, which can be parametrized by a
data frame or list, compute Pareto-smoothed importance weights and
attach them to the specification object, for further calculation and
plotting.

## Usage

``` r
adjust_weights(spec, object, data = NULL, keep_bad = FALSE, incl_orig = TRUE)
```

## Arguments

- spec:

  An object of class `adjustr_spec`, probably produced by
  [`make_spec`](https://corymccartan.github.io/adjustr/reference/make_spec.md),
  containing the new sampling sampling statements to replace their
  counterparts in the original Stan model, and the data, if any, by
  which these sampling statements are parametrized.

- object:

  A model object, either of type
  [`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.html),
  `stanreg` (from rstanarm), `brmsfit` (from brms), or a list with two
  elements: `model` containing a `CmdStanModel`, and `fit` containing a
  `CmdStanMCMC` object (both from the cmdstanr package).

- data:

  The data that was used to fit the model in `object`. Required only if
  one of the new sampling specifications involves Stan data variables.

- keep_bad:

  When `FALSE` (the strongly recommended default), alternate
  specifications which deviate too much from the original posterior, and
  which as a result cannot be reliably estimated using importance
  sampling (i.e., if the Pareto shape parameter is larger than 0.7),
  have their weights discarded—weights are set to `NA_real_`.

- incl_orig:

  When `TRUE`, include a row for the original model specification, with
  all weights equal. Can facilitate comparison and plotting later.

## Value

A tibble, produced by converting the provided `specs` to a tibble (see
[`as.data.frame.adjustr_spec`](https://corymccartan.github.io/adjustr/reference/as.data.frame.adjustr_spec.md)),
and adding columns `.weights`, containing vectors of weights for each
draw, and `.pareto_k`, containing the diagnostic Pareto shape
parameters. Values greater than 0.7 indicate that importance sampling is
not reliable. If `incl_orig` is `TRUE`, a row is added for the original
model specification. Weights can be extracted with the
[`pull.adjustr_weighted`](https://corymccartan.github.io/adjustr/reference/pull.adjustr_weighted.md)
method. The returned object also includes the model sample draws, in the
`draws` attribute.

## Details

This function does the bulk of the sensitivity analysis work. It
operates by parsing the model code from the provided Stan object,
extracting the parameters and their sampling statements. It then uses R
metaprogramming/tidy evaluation tools to flexibly evaluate the log
density for each draw and each sampling statement, under the original
and alternative specifications. From these, the function computes the
overall importance weight for each draw and performs Pareto-smoothed
importance sampling. All of the work is performed in R, without
recompiling or refitting the Stan model.

## References

Vehtari, A., Simpson, D., Gelman, A., Yao, Y., & Gabry, J. (2024).
Pareto smoothed importance sampling. *Journal of Machine Learning
Research*, 25(72), 1–58.
[doi:10.48550/arXiv.1507.02646](https://doi.org/10.48550/arXiv.1507.02646)

## See also

[`make_spec`](https://corymccartan.github.io/adjustr/reference/make_spec.md),
[`summarize.adjustr_weighted`](https://corymccartan.github.io/adjustr/reference/summarize.adjustr_weighted.md),
[`spec_plot`](https://corymccartan.github.io/adjustr/reference/spec_plot.md)

## Examples

``` r
# \donttest{
spec = make_spec(eta ~ student_t(df, 0, 1), df=4:10)
adjust_weights(spec, eightschools_m, keep_bad=TRUE)
#> Loading required namespace: rstan
#> # A tibble: 8 × 4
#>      df .samp                     .weights   .pareto_k
#>   <int> <chr>                     <list>         <dbl>
#> 1     4 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf
#> 2     5 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf
#> 3     6 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf
#> 4     7 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf
#> 5     8 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf
#> 6     9 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf
#> 7    10 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf
#> 8    NA <original model>          <dbl [20]>      -Inf
# }
```
