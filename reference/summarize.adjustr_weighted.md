# Summarize Posterior Distributions Under Alternative Model Specifications

Uses weights computed in
[`adjust_weights`](https://corymccartan.github.io/adjustr/reference/adjust_weights.md)
to compute posterior summary statistics. These statistics can be
compared against their reference values to quantify the sensitivity of
the model to aspects of its specification.

## Usage

``` r
# S3 method for class 'adjustr_weighted'
summarise(.data, ..., .resampling = FALSE, .model_data = NULL)

# S3 method for class 'adjustr_weighted'
summarize(.data, ..., .resampling = FALSE, .model_data = NULL)
```

## Arguments

- .data:

  An `adjustr_weighted` object.

- ...:

  Name-value pairs of expressions. The name of each argument will be the
  name of a new variable, and the value will be computed for the
  posterior distribution of eight alternative specification. For
  example, a value of `mean(theta)` will compute the posterior mean of
  `theta` for each alternative specification.

  Also supported is the custom function `wasserstein`, which computes
  the Wasserstein-p distance between the posterior distribution of the
  provided expression under the new model and under the original model,
  with `p=1` the default. Lower the `spacing` parameter from the default
  of 0.005 to compute a finer (but slower) approximation.

  The arguments in `...` are automatically quoted and evaluated in the
  context of `.data`. They support unquoting and splicing.

- .resampling:

  Whether to compute summary statistics by first resampling the data
  according to the weights. Defaults to `FALSE`, but will be used for
  any summary statistic that is not `mean`, `var` or `sd`.

- .model_data:

  Stan model data, if not provided in the earlier call to
  [`adjust_weights`](https://corymccartan.github.io/adjustr/reference/adjust_weights.md).

## Value

An `adjustr_weighted` object, with the new columns specified in `...`
added.

## See also

[`adjust_weights`](https://corymccartan.github.io/adjustr/reference/adjust_weights.md),
[`spec_plot`](https://corymccartan.github.io/adjustr/reference/spec_plot.md)

## Examples

``` r
# \donttest{
spec = make_spec(eta ~ student_t(df, 0, 1), df=4:10)
adjusted = adjust_weights(spec, eightschools_m, keep_bad=TRUE)

summarize(adjusted, mean(mu), var(mu))
#> # A tibble: 8 × 6
#>      df .samp                     .weights   .pareto_k `mean(mu)` `var(mu)`
#>   <int> <chr>                     <list>         <dbl>      <dbl>     <dbl>
#> 1     4 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf       7.46      15.3
#> 2     5 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf       7.44      15.5
#> 3     6 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf       7.42      15.7
#> 4     7 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf       7.41      15.8
#> 5     8 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf       7.39      15.9
#> 6     9 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf       7.38      16.0
#> 7    10 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf       7.38      16.1
#> 8    NA <original model>          <dbl [20]>      -Inf       7.27      17.1
# }
```
