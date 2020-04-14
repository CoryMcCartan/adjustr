# build the parser and store it
get_parser = function() {
    if (env_has(pkg_env, "parse_func")) {
        return(pkg_env$parse_func)
    } else {
        grammar_file = system.file("stan.g", package="adjustr")
        grammar = paste(readLines(grammar_file), collapse="\n")
        pkg_env$parse_func = dparser::dparse(grammar)
        return(pkg_env$parse_func)
    }
}

# Parse Satan `model_code` into a data frame which represents the parsing tree
parse_model = function(model_code) {
    capture.output(
    get_parser()(model_code, function(name, value, pos, depth) {
        cat(paste0('"', name, '","', value, '","', pos, '","', depth, '"\n'))
    })
    ) %>%
        paste(collapse="\n") %>%
        paste0("name,value,pos,depth\n", .) %>%
        read.csv(text=., as.is=T) %>%
        mutate(i = 1:n())
}

# Take a parsing tree and return a named vector, with the names matching the
# model's variable names and the values representing the program blocks they
# are defined in
get_variables = function(parsed_model) {
    prog_sections = filter(parsed_model, stringr::str_starts(name, "program__"), pos==-2)
    sec_names = stringr::str_extract(prog_sections$value, "^.+(?=\\{)") %>%
        stringr::str_trim()
    id_section = Vectorize(function(i)
        sec_names[which.max(i <= c(prog_sections$i, Inf)) - 1])

    prog_vars_d = filter(parsed_model, name=="var_decl", pos==1)
    prog_vars = id_section(prog_vars_d$i)
    names(prog_vars) = prog_vars_d$value
    prog_vars
}

# Take a parsing tree and return a list of formulas, one for each samping statement
# in the model
get_sampling_stmts_list = function(parsed_model) {
    parsed_model %>%
        filter(name == "sampling_statement", pos == -2) %>%
        purrr::pmap(function(value, ...) as.formula(value, env=global_env()))
}

# Take a list of prior formulas and return a matching list of sampling statements
match_priors_to_samp = function(priors, parsed_samp) {
    ref_vars = map_chr(parsed_samp, ~ as.character(f_lhs(.)))
    pri_vars = map_chr(priors, ~ as.character(f_lhs(.)))
    indices = match(pri_vars, ref_vars)
    # check that every prior was matched
    if (any(is.na(indices))) {
        stop(paste("No matching sampling statement found for prior",
                   priors[which.max(is.na(indices))] ))
    }
    parsed_samp[indices]
}
