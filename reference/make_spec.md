# Set Up Model Adjustment Specifications

Takes a set of new sampling statements, which can be parametrized by
other arguments, data frames, or lists, and creates an `adjustr_spec`
object suitable for use in
[`adjust_weights`](https://corymccartan.github.io/adjustr/reference/adjust_weights.md).

## Usage

``` r
make_spec(...)

# S3 method for class 'adjustr_spec'
print(x, ...)

# S3 method for class 'adjustr_spec'
length(x)
```

## Arguments

- ...:

  Model specification. Each argument can either be a formula, a named
  vector, data frames, or lists.

  Formula arguments provide new sampling statements to replace their
  counterparts in the original Stan model. All such formulas must be of
  the form `variable ~ distribution(parameters)`, where `variable` and
  `parameters` are Stan data variables or parameters, or are provided by
  other arguments to this function (see below), and where `distribution`
  matches one of the univariate [Stan
  distributions](https://mc-stan.org/docs/2_22/functions-reference/conventions-for-probability-functions.html).
  Arithmetic expressions of parameters are also allowed, but care must
  be taken with multivariate parameter arguments. Since specifications
  are passed as formulas, R's arithmetic operators are used, not Stan's.
  As a result, matrix and elementwise multipilcation in Stan sampling
  statements may not be interpreted correctly. Moving these computations
  out of sampling statements and into local variables will ensure
  correct results.

  For named vector arguments, each entry of the vector will be
  substituted into the corresponding parameter in the sampling
  statements. For data frames, each entry in each column will be
  substituted into the corresponding parameter in the sampling
  statements.

  List arguments are coerced to data frames. They can either be lists of
  named vectors, or lists of lists of single-element named vectors.

  The lengths of all parameter arguments must be consistent. Named
  vectors can have length 1 or must have length equal to the number of
  rows in all data frame arguments and the length of list arguments.

- x:

  An `adjustr_spec` object (for `print` and `length` methods).

## Value

An object of class `adjustr_spec`, which is essentially a list with two
elements: `samp`, which is a list of sampling formulas, and `params`,
which is a list of lists of parameters. Core [dplyr
verbs](https://corymccartan.github.io/adjustr/reference/dplyr.adjustr_spec.md)
which don't involve grouping
([`filter`](https://dplyr.tidyverse.org/reference/filter.html),
[`arrange`](https://dplyr.tidyverse.org/reference/arrange.html),
[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html),
[`select`](https://dplyr.tidyverse.org/reference/select.html),
[`rename`](https://dplyr.tidyverse.org/reference/rename.html), and
[`slice`](https://dplyr.tidyverse.org/reference/slice.html)) are
supported and operate on the underlying table of specification
parameters.

## See also

[`adjust_weights`](https://corymccartan.github.io/adjustr/reference/adjust_weights.md),
[`summarize.adjustr_weighted`](https://corymccartan.github.io/adjustr/reference/summarize.adjustr_weighted.md),
[`spec_plot`](https://corymccartan.github.io/adjustr/reference/spec_plot.md)

## Examples

``` r
make_spec(eta ~ cauchy(0, 1))
#> Sampling specifications:
#> eta ~ cauchy(0, 1)
#> <environment: 0x55873de21b60>

make_spec(eta ~ student_t(df, 0, 1), df=1:10)
#> Sampling specifications:
#> eta ~ student_t(df, 0, 1)
#> <environment: 0x55873de21b60>
#> 
#> Specification parameters:
#>  df
#>   1
#>   2
#>   3
#>   4
#>   5
#>   6
#>   7
#>   8
#>   9
#>  10

params = tidyr::crossing(df=1:10, infl=c(1, 1.5, 2))
make_spec(eta ~ student_t(df, 0, 1),
          y ~ normal(theta, infl*sigma),
          params)
#> Sampling specifications:
#> eta ~ student_t(df, 0, 1)
#> <environment: 0x55873de21b60>
#> y ~ normal(theta, infl * sigma)
#> <environment: 0x55873de21b60>
#> 
#> Specification parameters:
#>  df infl
#>   1    1
#>   1  1.5
#>   1    2
#>   2    1
#>   2  1.5
#>   2    2
#>   3    1
#>   3  1.5
#>   3    2
#>   4    1
#>   4  1.5
#>   4    2
#>   5    1
#>   5  1.5
#>   5    2
#>  [ reached 'max' / getOption("max.print") -- omitted 15 rows ]
```
