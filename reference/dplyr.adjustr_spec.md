# `dplyr` Methods for `adjustr_spec` Objects

Core [`dplyr`](https://dplyr.tidyverse.org/reference/dplyr-package.html)
verbs which don't involve grouping
([`filter`](https://dplyr.tidyverse.org/reference/filter.html),
[`arrange`](https://dplyr.tidyverse.org/reference/arrange.html),
[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html),
[`select`](https://dplyr.tidyverse.org/reference/select.html),
[`rename`](https://dplyr.tidyverse.org/reference/rename.html), and
[`slice`](https://dplyr.tidyverse.org/reference/slice.html)) are
implemented and operate on the underlying table of specification
parameters.

## Usage

``` r
# S3 method for class 'adjustr_spec'
filter(.data, ..., .preserve = FALSE)

# S3 method for class 'adjustr_spec'
arrange(.data, ...)

# S3 method for class 'adjustr_spec'
rename(.data, ...)

# S3 method for class 'adjustr_spec'
select(.data, ...)

# S3 method for class 'adjustr_spec'
slice(.data, ..., .preserve = FALSE)
```

## Arguments

- .data:

  the `adjustr_spec` object

- ...:

  additional arguments to underlying method

- .preserve:

  as in `filter` and `slice`

## Value

A modified `adjustr_spec` object.

## Examples

``` r
spec = make_spec(eta ~ student_t(df, 0, 1), df=1:10)

arrange(spec, desc(df))
#> Sampling specifications:
#> eta ~ student_t(df, 0, 1)
#> <environment: 0x55eba5f15680>
#> 
#> Specification parameters:
#>  df
#>  10
#>   9
#>   8
#>   7
#>   6
#>   5
#>   4
#>   3
#>   2
#>   1
slice(spec, 4:7)
#> Sampling specifications:
#> eta ~ student_t(df, 0, 1)
#> <environment: 0x55eba5f15680>
#> 
#> Specification parameters:
#>  df
#>   4
#>   5
#>   6
#>   7
```
