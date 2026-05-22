# Get Importance Resampling Indices From Weights

Takes a vector of weights, or data frame or list containing sets of
weights, and resamples indices for use in later computation.

## Usage

``` r
get_resampling_idxs(x, frac = 1, replace = TRUE)
```

## Arguments

- x:

  A vector of weights, a list of weight vectors, or a data frame of type
  `adjustr_weighted` containing a `.weights` list-column of weights.

- frac:

  A real number giving the fraction of draws to resample; the default,
  1, resamples all draws. Smaller values should be used when
  `replace=FALSE`.

- replace:

  Whether sampling should be with replacement. When weights are extreme
  it may make sense to use `replace=FALSE`, but accuracy is not
  guaranteed in these cases.

## Value

A vector, list, or data frame, depending of the type of `x`, containing
the sampled indices. If any weights are `NA`, the indices will also be
`NA`.

## Examples

``` r
spec = make_spec(eta ~ student_t(df, 0, 1), df=4:10)
adjusted = adjust_weights(spec, eightschools_m, keep_bad=TRUE)

get_resampling_idxs(adjusted)
#> # A tibble: 8 × 5
#>      df .samp                     .weights   .pareto_k .idxs     
#>   <int> <chr>                     <list>         <dbl> <list>    
#> 1     4 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf <int [20]>
#> 2     5 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf <int [20]>
#> 3     6 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf <int [20]>
#> 4     7 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf <int [20]>
#> 5     8 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf <int [20]>
#> 6     9 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf <int [20]>
#> 7    10 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf <int [20]>
#> 8    NA <original model>          <dbl [20]>      -Inf <int [20]>
get_resampling_idxs(adjusted, frac=0.5, replace=FALSE)
#> # A tibble: 8 × 5
#>      df .samp                     .weights   .pareto_k .idxs     
#>   <int> <chr>                     <list>         <dbl> <list>    
#> 1     4 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf <int [10]>
#> 2     5 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf <int [10]>
#> 3     6 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf <int [10]>
#> 4     7 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf <int [10]>
#> 5     8 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf <int [10]>
#> 6     9 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf <int [10]>
#> 7    10 eta ~ student_t(df, 0, 1) <dbl [20]>       Inf <int [10]>
#> 8    NA <original model>          <dbl [20]>      -Inf <int [10]>
```
