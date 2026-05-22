# adjustr: Stan Model Adjustments and Sensitivity Analyses using Importance Sampling

Assess the sensitivity of a Bayesian model (fitted using 'Stan' via
'rstan', 'brms', or 'cmdstanr') to the specification of its likelihood
and priors. Users provide a series of alternate sampling specifications,
and the package uses Pareto-smoothed importance sampling (PSIS) to
estimate posterior quantities of interest under each specification,
without needing to refit the model. Methods are based on Vehtari,
Simpson, Gelman, Yao, and Gabry (2024)
[doi:10.48550/arXiv.1507.02646](https://doi.org/10.48550/arXiv.1507.02646)
.

## See also

Useful links:

- <https://corymccartan.com/adjustr/>

- Report bugs at <https://github.com/CoryMcCartan/adjustr/issues>

## Author

**Maintainer**: Cory McCartan <mccartan@psu.edu> \[copyright holder\]

Authors:

- Cory McCartan <mccartan@psu.edu> \[copyright holder\]
