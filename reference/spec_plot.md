# Plot Posterior Quantities of Interest Under Alternative Model Specifications

Uses weights computed in
[`adjust_weights`](https://corymccartan.github.io/adjustr/reference/adjust_weights.md)
to plot posterior quantities of interest versus specification parameters

## Usage

``` r
spec_plot(
  x,
  by,
  post,
  only_mean = FALSE,
  ci_level = 0.8,
  outer_level = 0.95,
  ...
)
```

## Arguments

- x:

  An `adjustr_weighted` object.

- by:

  The x-axis variable, which is usually one of the specification
  parameters. Can be set to `1` if there is only one specification.
  Automatically quoted and evaluated in the context of `x`.

- post:

  The posterior quantity of interest, to be computed for each resampled
  draw of each specification. Should evaluate to a single number for
  each draw. Automatically quoted and evaluated in the context of `x`.

- only_mean:

  Whether to only plot the posterior mean. May be more stable.

- ci_level:

  The inner credible interval to plot. Central 100\*ci_level posterior
  draws.

- outer_level:

  The outer credible interval to plot.

- ...:

  Ignored.

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html) object
which can be further customized with the **ggplot2** package.

## See also

[`adjust_weights`](https://corymccartan.github.io/adjustr/reference/adjust_weights.md),
[`summarize.adjustr_weighted`](https://corymccartan.github.io/adjustr/reference/summarize.adjustr_weighted.md)

## Examples

``` r
spec = make_spec(eta ~ student_t(df, 0, 1), df=4:10)
adjusted = adjust_weights(spec, eightschools_m, keep_bad=TRUE)

spec_plot(adjusted, df, mu, only_mean=TRUE)

```
