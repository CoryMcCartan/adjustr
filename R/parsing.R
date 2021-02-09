# regexes
identifier = "[a-zA-Z][a-zA-Z0-9_]*"
re_stmt = paste0("(int|real|(?:unit_|row_)?vector|(?:positive_)?ordered|simplex",
    "|(?:cov_|corr_)?matrix|cholesky_factor(?:_corr|_cov)?)(?:<.+>)?(?:\\[.+\\])?",
    " (", identifier, ")(?:\\[.+\\])? ?=?")
#re_block = paste0(block_names, " ?\\{ ?(.+) ?\\} ?", block_names, "")
re_block = "((?:transformed )?data|(?:transformed )?parameters|model|generated quantities)"
re_samp = paste0("(", identifier, " ?~[^~{}]+)")
re_samp2 = paste0("target ?\\+= ?(", identifier, ")_lp[md]f\\((",
                  identifier, ")(?:| ?[|]? ?(.+))\\)")

# Extract variable name from variable declaration, or return NA if no declaration
get_variables = function(statement) {
    matches = stringr::str_match(statement, re_stmt)[,3]
    matches[!is.na(matches)]
}

get_sampling = function(statement) {
    samps = stringr::str_match(statement, re_samp)[,2]
    samps2 = stringr::str_match(statement, re_samp2)#[,,3]
    samps2_rearr = paste0(samps2[,3], " ~ ", samps2[,2], "(", coalesce(samps2[,4], ""), ")")
    stmts = c(samps[!is.na(samps)], samps2_rearr[!is.na(samps2[,1])])
    map(stmts, ~ stats::as.formula(., env=empty_env()))
}

# Parse Stan `model_code` into a list with two elements: `vars` named
# vector, with the names matching the model's variable names and the values
# representing the program blocks they are defined in; `samp` is a list of
# sampling statements (as formulas)
parse_model = function(model_code) {
    clean_code = stringr::str_replace_all(model_code, "//.*", "") %>%
        stringr::str_replace_all("/\\*[^*]*\\*+(?:[^/*][^*]*\\*+)*/", "") %>%
        stringr::str_replace_all("\\n", " ") %>%
        stringr::str_replace_all("\\s\\s+", " ")

    block_names = stringr::str_extract_all(clean_code, re_block)[[1]]
    if (length(block_names)==0) return(list(vars=character(0), samps=list()))

    block_locs = rbind(stringr::str_locate_all(clean_code, re_block)[[1]],
                       c(nchar(clean_code), NA))
    blocks = map(1:length(block_names), function(i) {
        block = stringr::str_sub(clean_code, block_locs[i,2]+1, block_locs[i+1,1])
        start = stringr::str_locate_all(block, stringr::fixed("{"))[[1]][1,1] + 1
        end = utils::tail(stringr::str_locate_all(block, stringr::fixed("}"))[[1]][,1], 1) - 1
        stringr::str_trim(stringr::str_sub(block, start+1, end-1))
    })
    names(blocks) = block_names

    statements = map(blocks, ~ stringr::str_split(., "; ?", simplify=T)[1,])

    vars = map(statements, get_variables)
    vars = purrr::flatten_chr(purrr::imap(vars, function(name, block) {
        block = rep(block, length(name))
        names(block) = name
        block
    }))


    samps = map(statements, get_sampling)
    names(samps) = NULL
    samps = flatten(samps)

    parameters = names(vars)[vars == "parameters"]
    sampled_pars = map(samps, ~ as.character(f_lhs(.))) %>%
        purrr::as_vector()
    uniform_pars = setdiff(parameters, sampled_pars)
    uniform_samp = paste0(uniform_pars, " ~ uniform(-1e100, 1e100)")
    uniform_samp = map(uniform_samp, ~ stats::as.formula(., env=empty_env()))

    list(vars=vars, samp=c(samps, uniform_samp))
}


# Take a list of provided sampling formulas and return a matching list of
# sampling statements from a reference list
match_sampling_stmts = function(new_samp, ref_samp) {
    ref_vars = map(ref_samp, ~ as.character(f_lhs(.))) %>%
        purrr::as_vector()
    new_vars = map(new_samp, ~ as.character(f_lhs(.))) %>%
        purrr::as_vector()
    indices = match(new_vars, ref_vars)
    # check that every prior was matched
    if (any(is.na(indices))) {
        stop("No matching sampling statement found for ",
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

