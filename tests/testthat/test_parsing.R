context("Stan model parsing")

test_that("Empty model handled correctly", {
    parsed = parse_model("")
    expect_equal(length(parsed$vars), 0)
    expect_equal(length(parsed$samps), 0)
})

test_that("Correct parsed variables", {
    correct_vars = c(J= "data", y= "data", sigma="data", mu="parameters",
                     tau="parameters", eta="parameters",
                     theta="transformed parameters")
    code = eightschools_m@stanmodel@model_code
    parsed = parse_model(code)
    expect_equal(parsed$vars, correct_vars)
})

test_that("Correct parsed sampling statements", {
    correct_samp = list(eta ~ std_normal(), y ~ normal(theta, sigma),
                        mu ~ uniform(-1e+100, 1e+100), tau ~ uniform(-1e+100, 1e+100))
    code = eightschools_m@stanmodel@model_code
    parsed = parse_model(code)
    expect_equal(parsed$samp, correct_samp)
})

test_that("Provided sampling statements can be matched to model", {
    model_samp = list(eta ~ std_normal(), y ~ normal(theta, sigma))
    prov_samp = list(eta ~ exponential(5))
    matched = match_sampling_stmts(prov_samp, model_samp)
    expect_false(identical(matched, prov_samp))
    expect_length(matched, 1)
    expect_equal(matched[[1]], model_samp[[1]])
})

test_that("Extra sampling statements not in model throw an error", {
    model_samp = list(eta ~ std_normal(), y ~ normal(theta, sigma))
    prov_samp = list(eta ~ exponential(5), x ~ normal(theta, sigma))
    expect_error(match_sampling_stmts(prov_samp, model_samp),
                 "No matching sampling statement found for x ~ normal\\(theta, sigma\\)")
})

test_that("Variables are correctly extracted from sampling statements", {
    expect_equal(get_stmt_vars(y ~ normal(theta %*% eta/2, 2%%4 + sigma)),
                 c("y", "theta", "eta", "sigma"))
    expect_error(get_stmt_vars(y ~ .), "y ~ \\. does not")
})

# --- Tests for new token-based parser ---

test_that("New array syntax is parsed correctly", {
    code = "
    data {
        array[N] real y;
        array[N, K] int x;
    }
    parameters {
        array[3] vector[N] beta;
    }
    model {
        y ~ normal(0, 1);
    }"
    parsed = parse_model(code)
    expect_equal(parsed$vars[c("y", "x", "beta")],
                 c(y = "data", x = "data", beta = "parameters"))
})

test_that("Old-style array syntax is still parsed correctly", {
    code = "
    data {
        int<lower=0> N;
        real y[N];
        real<lower=0> sigma[N];
    }
    parameters {
        real mu;
    }
    model {
        y ~ normal(mu, sigma);
    }"
    parsed = parse_model(code)
    expect_equal(parsed$vars,
                 c(N = "data", y = "data", sigma = "data", mu = "parameters"))
})

test_that("target += with lpdf/lpmf is parsed correctly", {
    code = "
    data {
        int<lower=1> N;
        array[N] int Y;
    }
    parameters {
        real mu;
    }
    model {
        target += normal_lpdf(mu | 0, 10);
        target += poisson_lpmf(Y | mu);
    }"
    parsed = parse_model(code)
    expect_length(parsed$samp, 2)
    expect_equal(parsed$samp[[1]], mu ~ normal(0, 10))
    expect_equal(parsed$samp[[2]], Y ~ poisson(mu))
})

test_that("Statements inside control flow are skipped", {
    code = "
    data {
        int<lower=1> N;
        array[N] real y;
    }
    parameters {
        real mu;
        real<lower=0> sigma;
    }
    model {
        mu ~ normal(0, 10);
        if (N > 0) {
            y ~ normal(mu, sigma);
        }
    }"
    parsed = parse_model(code)
    # Only the top-level statement should be found
    sampled = vapply(parsed$samp, function(s) deparse(rlang::f_lhs(s)), "")
    expect_true("mu" %in% sampled)
    expect_false("y" %in% sampled)
})

test_that("Profile blocks are handled (statements inside skipped)", {
    code = '
    data {
        int N;
        array[N] real y;
    }
    parameters {
        real mu;
        real<lower=0> sigma;
    }
    model {
        mu ~ normal(0, 10);
        profile("likelihood") {
            y ~ normal(mu, sigma);
        }
    }'
    parsed = parse_model(code)
    sampled = vapply(parsed$samp, function(s) deparse(rlang::f_lhs(s)), "")
    expect_true("mu" %in% sampled)
    # profile block increases brace depth, so y ~ ... is skipped
    expect_false("y" %in% sampled)
})

test_that("Block comments and line comments are stripped", {
    code = "
    data {
        /* this is a block comment
           spanning multiple lines */
        int N;
        // real ghost_var;
    }
    model {
    }"
    parsed = parse_model(code)
    expect_equal(names(parsed$vars), "N")
    expect_true(! "ghost_var" %in% names(parsed$vars))
})

test_that("New constrained types are parsed", {
    code = "
    parameters {
        sum_to_zero_vector[5] beta;
        simplex[K] theta;
        ordered[3] cutpoints;
        cholesky_factor_corr[K] L;
    }
    model {
        beta ~ normal(0, 1);
        theta ~ dirichlet(alpha);
        cutpoints ~ normal(0, 5);
        L ~ lkj_corr_cholesky(2);
    }"
    parsed = parse_model(code)
    expect_equal(parsed$vars[c("beta", "theta", "cutpoints", "L")],
                 c(beta = "parameters", theta = "parameters",
                   cutpoints = "parameters", L = "parameters"))
})

test_that("Tokenizer handles Stan operators correctly", {
    tokens = tokenize_stan("target += normal_lpdf(y | mu, sigma);")
    expect_true("target" %in% tokens$value)
    expect_true("+=" %in% tokens$value)
    expect_true("|" %in% tokens$value)
})

test_that("Functions block is parsed without errors", {
    code = '
    functions {
        real my_fun(real x) {
            return x * 2;
        }
    }
    data {
        real y;
    }
    parameters {
        real mu;
    }
    model {
        mu ~ normal(0, 10);
        y ~ normal(mu, 1);
    }'
    parsed = parse_model(code)
    expect_equal(parsed$vars[c("y", "mu")], c(y = "data", mu = "parameters"))
    expect_length(parsed$samp, 2)
    expect_equal(parsed$samp[[2]], y ~ normal(mu, 1))
})

test_that("Truncated sampling statements are parsed (truncation stripped)", {
    code = "
    data {
        real<lower=0> y;
    }
    parameters {
        real mu;
        real<lower=0> sigma;
    }
    model {
        mu ~ normal(0, 10);
        sigma ~ exponential(1);
        y ~ normal(mu, sigma) T[0, ];
    }"
    parsed = parse_model(code)
    sampled = vapply(parsed$samp, function(s) deparse(rlang::f_lhs(s)), "")
    expect_true("y" %in% sampled)
    y_samp = parsed$samp[[which(sampled == "y")]]
    expect_equal(y_samp, y ~ normal(mu, sigma))
})

test_that("Multiple variables declared on separate lines in same block", {
    code = "
    parameters {
        real alpha;
        real beta;
        real<lower=0> sigma;
        vector[10] gamma;
    }
    model {
        alpha ~ normal(0, 1);
        beta ~ normal(0, 1);
        sigma ~ exponential(1);
        gamma ~ normal(0, 1);
    }"
    parsed = parse_model(code)
    expect_equal(sort(names(parsed$vars)),
                 sort(c("alpha", "beta", "sigma", "gamma")))
    expect_true(all(parsed$vars == "parameters"))
})

test_that("Transformed data and generated quantities variables are found", {
    code = "
    data {
        int N;
    }
    transformed data {
        real log_N = log(N);
    }
    parameters {
        real mu;
    }
    generated quantities {
        real y_rep = normal_rng(mu, 1);
    }
    model {
        mu ~ normal(0, 1);
    }"
    parsed = parse_model(code)
    expect_equal(parsed$vars["log_N"], c(log_N = "transformed data"))
    expect_equal(parsed$vars["y_rep"], c(y_rep = "generated quantities"))
})

test_that("Variable with assignment in declaration is parsed", {
    code = "
    transformed parameters {
        real<lower=0> disc = 1;
        vector[3] theta = mu + tau * eta;
    }
    parameters {
        real mu;
        real<lower=0> tau;
        vector[3] eta;
    }
    model {
        mu ~ normal(0, 1);
        tau ~ exponential(1);
        eta ~ std_normal();
    }"
    parsed = parse_model(code)
    expect_equal(parsed$vars["disc"], c(disc = "transformed parameters"))
    expect_equal(parsed$vars["theta"], c(theta = "transformed parameters"))
})

test_that("target += with _lupdf and _lupmf suffixes are parsed", {
    code = "
    data {
        real y;
    }
    parameters {
        real mu;
    }
    model {
        target += normal_lupdf(y | mu, 1);
    }"
    parsed = parse_model(code)
    expect_length(parsed$samp, 2)  # 1 explicit + 1 implicit uniform for mu
    expect_equal(parsed$samp[[1]], y ~ normal(mu, 1))
})

test_that("Mixed tilde and target += in same model", {
    code = "
    data {
        int N;
        array[N] real y;
    }
    parameters {
        real mu;
        real<lower=0> sigma;
    }
    model {
        mu ~ normal(0, 10);
        sigma ~ exponential(1);
        target += normal_lpdf(y | mu, sigma);
    }"
    parsed = parse_model(code)
    # 3 explicit statements, no implicit uniforms (all params have priors)
    expect_length(parsed$samp, 3)
    expect_equal(parsed$samp[[1]], mu ~ normal(0, 10))
    expect_equal(parsed$samp[[2]], sigma ~ exponential(1))
    expect_equal(parsed$samp[[3]], y ~ normal(mu, sigma))
})
