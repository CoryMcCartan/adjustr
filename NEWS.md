# adjustr 0.1.2

* Add support for `cmdstanr` objects by passing a list containing the fit and the model object.

* Fix bug in parsing code that caused an error with some `target +=` model statements.


# adjustr 0.1.1

* Improved documentation and additional references.

* Fix bug in `extract_samp_stmts()` which prevented `brmsfit` objects from being used directly.


# adjustr 0.1.0

* Initial release.

* Basic workflow implemented: `make_spec()`, `adjust_weights()`, and `summarize()`/`spec_plot()`.