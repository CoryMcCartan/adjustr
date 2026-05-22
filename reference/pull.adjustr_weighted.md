# Extract Weights From an `adjustr_weighted` Object

This function modifies the default behavior of
[`dplyr::pull`](https://dplyr.tidyverse.org/reference/pull.html) to
extract the `.weights` column.

## Usage

``` r
# S3 method for class 'adjustr_weighted'
pull(.data, var = ".weights", name = NULL, ...)
```

## Arguments

- .data:

  A table of data

- var:

  A variable, as in
  [`pull`](https://dplyr.tidyverse.org/reference/pull.html). The default
  returns the `.weights` column, and if there is only one row, it
  returns the first element of that column

- name:

  Ignored

- ...:

  Ignored

## Value

A numeric vector of weights (if a single row) or a list of numeric
vectors.
