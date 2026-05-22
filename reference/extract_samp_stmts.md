# Extract Model Sampling Statements From a Stan Model.

Prints a list of sampling statements extracted from the `model` block of
a Stan program, with each labelled "parameter" or "data" depending on
the type of variable being sampled.

## Usage

``` r
extract_samp_stmts(object)
```

## Arguments

- object:

  A [`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.html)
  model object.

## Value

Invisibly returns a list of sampling formulas.

## Examples

``` r
# \donttest{
extract_samp_stmts(eightschools_m)
#> Sampling statements for model 2c8d1d8a30137533422c438f23b83428:
#>   parameter   eta ~ std_normal()
#>   parameter   mu ~ uniform(-1e+100, 1e+100)
#>   parameter   tau ~ uniform(-1e+100, 1e+100)
#>   data        y ~ normal(theta, sigma)
#> Sampling statements for model 2c8d1d8a30137533422c438f23b83428:
#>   parameter   eta ~ std_normal()
#>   data        y ~ normal(theta, sigma)
# }
```
