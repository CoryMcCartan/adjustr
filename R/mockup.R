if (F) {
library(rlang)
library(dplyr)
library(purrr)
library(stringr)
library(dparser)
library(rstan)
rstan_options(auto_write = TRUE)

model_code = "data {
    int<lower=0> J;         // number of schools
    real y[J];              // estimated treatment effects
    real<lower=0> sigma[J]; // standard error of effect estimates
}
parameters {
    real mu;                // population treatment effect
    real<lower=0> tau;      // standard deviation in treatment effects
    vector[J] eta;          // unscaled deviation from mu by school
}
transformed parameters {
    vector[J] theta = mu + tau * eta;        // school treatment effects
}
model {
    eta ~ std_normal();
    y ~ normal(theta, sigma);
}"

model_d = list(J = 8,
               y = c(28,  8, -3,  7, -1,  1, 18, 12),
               sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
eightschools_m = stan(model_code=model_code, chains=2, data=model_d, warmup=500,
                      iter=510, save_dso=F, save_warmup=F)
eightschools_m@stanmodel@dso = new("cxxdso")
save(eightschools_m, file="tests/test_model.rda")

#slot(eightschools_m@stanmodel, "dso", F) = NULL
draws = extract(eightschools_m)



grammar = paste(readLines("R/stan.dpg"), collapse="\n")
parse_func = dparse(grammar, set_op_priority_from_rule=T, longest_match=T)


model_code = readr::read_file("~/Documents/Analyses/elections/president/stan/polls.stan")



prog_sections = filter(d, name=="program", value != "")
sec_names = str_trim(str_extract(prog_sections$value, "^.+(?=\\{)"))
id_section = Vectorize(function(i)
    sec_names[which.max(i <= c(prog_sections$i, Inf)) - 1])
prog_vars_d = filter(d, name=="var_decl", pos==1)
prog_vars = id_section(prog_vars_d$i)
names(prog_vars) = prog_vars_d$value

samp_stmts = d %>%
    filter(name == "sampling_statement", pos == -2) %>%
    pmap(function(value, ...) as.formula(value, env=global_env()))

samp_lhs = map_chr(samp_stmts, ~ prog_vars[as.character(f_lhs(.))])
get_rhs_deps = function(func) {
    map_chr(call_args(f_rhs(func)), ~ prog_vars[as.character(.)])
}
samp_rhs = map(samp_stmts, get_rhs_deps)



make_dens = function(f) {
    function(x) {
        function(...) {
            f(x, ..., log=T)
        }
    }
}
distrs = list(
    normal = dnorm,
    std_normal = dnorm,
    student_t = dt,
    exponential = dexp,
    gamma = function(x, alpha, beta, ...) dgamma(x, shape=alpha, rate=beta, ...)
)
distr_env = new_environment(map(distrs, make_dens), parent=global_env())


form = sigma_natl ~ gamma(alpha, 3/mean)
ref_form = samp_stmts[[5]]
xr = f_rhs(form)
xl = f_lhs(form)
draws = list(sigma_natl = rgamma(500, 2, 3/0.06))

combos = list(alpha = c(1, 2, 3),
                  mean = seq(0.01, 0.1, 0.01)) %>%
    cross_df

xr

eval_tidy(xr, combos, distr_env)
eval_tidy(xl, combos, distr_env)
call_fn(xr, distr_env)(eval_tidy(xl, draws))(4, 3)

call_fn(f_rhs(ref_form), distr_env)(eval_tidy(f_lhs(ref_form), draws))



specs = cross(list(df=1:4, mean=0))
make_rs_weights(sm, eta ~ student_t(df, mean, 1), specs)

}

