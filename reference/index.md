# Package index

## Model Adjustments

Core functions for sensitivity anlysis workflow.

- [`make_spec()`](https://corymccartan.github.io/adjustr/reference/make_spec.md)
  [`print(`*`<adjustr_spec>`*`)`](https://corymccartan.github.io/adjustr/reference/make_spec.md)
  [`length(`*`<adjustr_spec>`*`)`](https://corymccartan.github.io/adjustr/reference/make_spec.md)
  : Set Up Model Adjustment Specifications
- [`adjust_weights()`](https://corymccartan.github.io/adjustr/reference/adjust_weights.md)
  : Compute Pareto-smoothed Importance Weights for Alternative Model
  Specifications
- [`summarise(`*`<adjustr_weighted>`*`)`](https://corymccartan.github.io/adjustr/reference/summarize.adjustr_weighted.md)
  [`summarize(`*`<adjustr_weighted>`*`)`](https://corymccartan.github.io/adjustr/reference/summarize.adjustr_weighted.md)
  : Summarize Posterior Distributions Under Alternative Model
  Specifications
- [`spec_plot()`](https://corymccartan.github.io/adjustr/reference/spec_plot.md)
  : Plot Posterior Quantities of Interest Under Alternative Model
  Specifications

## Helper Functions

Various helper functions for examining a model or building sampling
specifications.

- [`extract_samp_stmts()`](https://corymccartan.github.io/adjustr/reference/extract_samp_stmts.md)
  : Extract Model Sampling Statements From a Stan Model.

- [`as.data.frame(`*`<adjustr_spec>`*`)`](https://corymccartan.github.io/adjustr/reference/as.data.frame.adjustr_spec.md)
  :

  Convert an `adjustr_spec` Object Into a Data Frame

- [`filter(`*`<adjustr_spec>`*`)`](https://corymccartan.github.io/adjustr/reference/dplyr.adjustr_spec.md)
  [`arrange(`*`<adjustr_spec>`*`)`](https://corymccartan.github.io/adjustr/reference/dplyr.adjustr_spec.md)
  [`rename(`*`<adjustr_spec>`*`)`](https://corymccartan.github.io/adjustr/reference/dplyr.adjustr_spec.md)
  [`select(`*`<adjustr_spec>`*`)`](https://corymccartan.github.io/adjustr/reference/dplyr.adjustr_spec.md)
  [`slice(`*`<adjustr_spec>`*`)`](https://corymccartan.github.io/adjustr/reference/dplyr.adjustr_spec.md)
  :

  `dplyr` Methods for `adjustr_spec` Objects

- [`get_resampling_idxs()`](https://corymccartan.github.io/adjustr/reference/get_resampling_idxs.md)
  : Get Importance Resampling Indices From Weights

- [`pull(`*`<adjustr_weighted>`*`)`](https://corymccartan.github.io/adjustr/reference/pull.adjustr_weighted.md)
  :

  Extract Weights From an `adjustr_weighted` Object

- [`adjustr`](https://corymccartan.github.io/adjustr/reference/adjustr-package.md)
  [`adjustr-package`](https://corymccartan.github.io/adjustr/reference/adjustr-package.md)
  : adjustr: 'Stan' Model Adjustments and Sensitivity Analyses using
  Importance Sampling
