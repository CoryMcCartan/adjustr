# Convert an `adjustr_spec` Object Into a Data Frame

Returns the data frame of specification parameters, with added columns
of the form `.samp_1`, `.samp_2`, ... for each sampling statement (or
just `.samp` if there is only one sampling statement).

## Usage

``` r
# S3 method for class 'adjustr_spec'
as.data.frame(x, ...)
```

## Arguments

- x:

  the `adjustr_spec` object

- ...:

  additional arguments to underlying method

## Value

A data frame with one row per specification and columns for each
parameter, plus `.samp` (or `.samp_1`, `.samp_2`, ...) columns
containing the sampling statement formulas.
