# build the parser and store it
get_parser = function() { # nocov start
    if (env_has(pkg_env, "parse_func")) {
        return(pkg_env$parse_func)
    } else {
        grammar_file = system.file("stan.g", package="adjustr")
        grammar = paste(readLines(grammar_file), collapse="\n")
        pkg_env$parse_func = dparser::dparse(grammar)
        return(pkg_env$parse_func)
    }
} # nocov end

# Parse Satan `model_code` into a data frame which represents the parsing tree
parse_model = function(model_code) {
    parser_output = utils::capture.output(
    get_parser()(model_code, function(name, value, pos, depth) {
        cat(stringr::str_glue('"{name}","{value}",{pos},{depth}\n\n'))
    })
    )
    parser_csv = paste0("name,value,pos,depth\n",
                        paste(parser_output, collapse="\n"))
    parsed = utils::read.csv(text=parser_csv, as.is=T)
    parsed$i = 1:nrow(parsed)
    parsed
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
get_sampling_stmts = function(parsed_model) {
    parsed_model %>%
        filter(name == "sampling_statement", pos == -2) %>%
        mutate(value = stringr::str_replace(value, " \\.\\*", "*")) %>%
        purrr::pmap(function(value, ...) stats::as.formula(value, env=empty_env()))
}
# Take a list of provided sampling formulas and return a matching list of
# sampling statements from a reference list
match_sampling_stmts = function(new_samp, ref_samp) {
    ref_vars = map_chr(ref_samp, ~ as.character(f_lhs(.)))
    new_vars = map_chr(new_samp, ~ as.character(f_lhs(.)))
    indices = match(new_vars, ref_vars)
    # check that every prior was matched
    if (any(is.na(indices))) {
        stop("No matching sampling statement found for prior ",
             new_samp[which.max(is.na(indices))],
             "\n  Check sampling statements and ensure that model data ",
             "has been provided.")
    }
    ref_samp[indices]
}

# Extract a list of variables from a sampling statement
# R versions of mathematical operators must be used
get_stmt_vars = function(stmt) {
    get_ast = function(x) purrr::map_if(as.list(x), is_call, get_ast)
    if (!is_call(f_rhs(stmt)))
        stop("Sampling statment ", format(stmt),
             " does not contain a distribution on the right-hand side.")
    # pull out variables from RHS
    rhs_vars = call_args(f_rhs(stmt)) %>%
        get_ast %>%
        unlist %>%
        purrr::discard(~ is(., "numeric")) %>%
        as.character %>%
        purrr::discard(~ . %in% c("`+`", "`-`", "`*`", "`/`", "`^`", "`%*%`", "`%%`"))
    c(as.character(f_lhs(stmt)), rhs_vars)
}

