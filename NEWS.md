# adjustr 0.2.0

* First CRAN submission.

* Replaced fragile regex-based Stan model parser with a robust token-based
  parser. Supports modern Stan syntax including the `array[]` declaration
  style, new constrained types (`sum_to_zero_vector`, etc.), and `target +=`
  statements with `_lupdf`/`_lupmf` suffixes.

* Expanded test suite from 108 to 197 tests, covering parsing, log probability
  calculations, specification creation, weight computation, and summary
  functions.

# adjustr 0.1.2

* Add support for `cmdstanr` objects by passing a list containing the fit and the model object.

* Fix bug in parsing code that caused an error with some `target +=` model statements.


# adjustr 0.1.1

* Improved documentation and additional references.

* Fix bug in `extract_samp_stmts()` which prevented `brmsfit` objects from being used directly.


# adjustr 0.1.0

* Initial release.

* Basic workflow implemented: `make_spec()`, `adjust_weights()`, and `summarize()`/`spec_plot()`.
