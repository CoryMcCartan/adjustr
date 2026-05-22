# Pre-fitted Eight Schools Model

A small
[`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.html)
object from the classic "eight schools" hierarchical model (Rubin,
1981), included for use in examples and tests. The model was fit with 2
chains of 1000 iterations (500 warmup) for a compact object size.

## Usage

``` r
eightschools_m
```

## Format

A [`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.html)
object.

## Examples

``` r
extract_samp_stmts(eightschools_m)
#> Sampling statements for model 2c8d1d8a30137533422c438f23b83428:
#>   parameter   eta ~ std_normal()
#>   parameter   mu ~ uniform(-1e+100, 1e+100)
#>   parameter   tau ~ uniform(-1e+100, 1e+100)
#>   data        y ~ normal(theta, sigma)
```
