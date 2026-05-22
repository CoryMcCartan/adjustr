# Token types used by the Stan tokenizer
# Order matters: longer/more-specific patterns first; word boundaries prevent
# partial matches (e.g. "in" inside "int").
TOKEN_TYPES = c(
    BLOCK_KW   = "(?:transformed\\s+data|transformed\\s+parameters|generated\\s+quantities|functions|parameters|data|model)\\b",
    TYPE_KW    = paste0("(?:complex_row_vector|complex_vector|complex_matrix|",
                        "cholesky_factor_corr|cholesky_factor_cov|",
                        "column_stochastic_matrix|row_stochastic_matrix|",
                        "sum_to_zero_vector|sum_to_zero_matrix|",
                        "positive_ordered|row_vector|cov_matrix|corr_matrix|",
                        "unit_vector|",
                        "ordered|simplex|complex|vector|matrix|array|tuple|void|",
                        "real|int)\\b"),
    KEYWORD    = "(?:continue|profile|jacobian|return|reject|target|break|while|print|else|for|if|in)\\b",
    NUMBER     = "(?:[0-9]+\\.?[0-9]*(?:[eE][+-]?[0-9]+)?|\\.?[0-9]+(?:[eE][+-]?[0-9]+)?)",
    STRING     = '"[^"]*"',
    IDENTIFIER = "[a-zA-Z_][a-zA-Z0-9_]*",
    PLUSASSIGN  = "\\+=",
    TILDE      = "~",
    LBRACE     = "\\{",
    RBRACE     = "\\}",
    LBRACK     = "\\[",
    RBRACK     = "\\]",
    LPAREN     = "\\(",
    RPAREN     = "\\)",
    LANGLE     = "<",
    RANGLE     = ">",
    SEMICOLON  = ";",
    COMMA      = ",",
    BAR        = "\\|",
    COLON      = ":",
    OP         = "(?:\\+|-|\\*|/|%|\\^|\\.\\*|\\./|\\.\\^|==|!=|<=|>=|&&|\\|\\||!|'|\\.)"
)

# Combined pattern: try each token type in order
TOKEN_RE = paste0("(", TOKEN_TYPES, ")", collapse="|")

# Tokenize Stan code into a data frame with columns: type, value
tokenize_stan = function(code) {
    # Strip comments
    code = stringr::str_replace_all(code, "//[^\n]*", "")
    code = stringr::str_replace_all(code, "/\\*[^*]*\\*+(?:[^/*][^*]*\\*+)*/", "")

    # Find all tokens
    matches = stringr::str_match_all(code, TOKEN_RE)[[1]]
    if (nrow(matches) == 0) {
        return(data.frame(type = character(0), value = character(0),
                          stringsAsFactors = FALSE))
    }

    values = matches[, 1]
    # Determine type by which capture group matched
    type_cols = matches[, -1, drop = FALSE]
    type_idx = apply(type_cols, 1, function(row) which(!is.na(row))[1])
    types = names(TOKEN_TYPES)[type_idx]

    # Reclassify: block keywords vs other keywords vs type keywords vs identifiers
    # Block keywords that were captured as KEYWORD or IDENTIFIER need reclassification
    # The regex priority handles this: BLOCK_KW is first, so "data" matches BLOCK_KW
    # But single words like "model" might match KEYWORD if preceded by another token
    # We rely on regex ordering: BLOCK_KW > KEYWORD > TYPE_KW > IDENTIFIER

    # Normalize whitespace in block keywords
    values = stringr::str_replace_all(values, "\\s+", " ")

    data.frame(type = types, value = values, stringsAsFactors = FALSE)
}


# Extract top-level blocks from token stream.
# Returns a named list of token data frames, one per block.
extract_blocks = function(tokens) {
    block_kw_idx = which(tokens$type == "BLOCK_KW")
    if (length(block_kw_idx) == 0) return(list())

    blocks = list()
    for (i in seq_along(block_kw_idx)) {
        bk_idx = block_kw_idx[i]
        block_name = tokens$value[bk_idx]

        # Find the opening brace
        open_idx = bk_idx + 1
        if (open_idx > nrow(tokens) || tokens$value[open_idx] != "{") next

        # Walk forward counting braces to find matching close
        depth = 0
        close_idx = NA
        for (j in open_idx:nrow(tokens)) {
            if (tokens$value[j] == "{") depth = depth + 1
            else if (tokens$value[j] == "}") depth = depth - 1
            if (depth == 0) { close_idx = j; break }
        }
        if (is.na(close_idx)) next

        # Extract tokens inside the block (excluding outer braces)
        if (close_idx > open_idx + 1) {
            block_tokens = tokens[(open_idx + 1):(close_idx - 1), , drop = FALSE]
            rownames(block_tokens) = NULL
        } else {
            block_tokens = data.frame(type = character(0), value = character(0),
                                      stringsAsFactors = FALSE)
        }
        blocks[[block_name]] = block_tokens
    }

    blocks
}


# Parse variable declarations from a block's token stream.
# Returns a character vector of variable names.
# Handles both old-style `real y` and new-style `array[N] real y`.
parse_declarations = function(tokens) {
    if (nrow(tokens) == 0) return(character(0))

    vars = character(0)
    i = 1
    n = nrow(tokens)
    depth = 0  # brace depth -- only parse at top level

    while (i <= n) {
        # Track brace depth
        if (tokens$value[i] == "{") { depth = depth + 1; i = i + 1; next }
        if (tokens$value[i] == "}") { depth = depth - 1; i = i + 1; next }

        # Only parse declarations at top level (depth 0)
        if (depth != 0) { i = i + 1; next }

        # Check for type keyword (possibly preceded by array dims)
        start_i = i
        is_decl = FALSE

        # Optional: array[...] prefix
        if (tokens$value[i] == "array" && tokens$type[i] == "TYPE_KW") {
            # Skip array[dims]
            i = i + 1
            i = skip_brackets(tokens, i, "[", "]")
        }

        # Optional: tuple(...) type
        if (i <= n && tokens$value[i] == "tuple" && tokens$type[i] == "TYPE_KW") {
            i = i + 1
            i = skip_brackets(tokens, i, "(", ")")
            is_decl = TRUE
        } else if (i <= n && tokens$type[i] == "TYPE_KW") {
            # Standard type keyword
            i = i + 1
            is_decl = TRUE
        }

        if (!is_decl) {
            # Not a declaration; skip to next semicolon at this depth
            i = start_i
            i = skip_to_semicolon(tokens, i, depth)
            next
        }

        # Optional: type constraints <lower=..., upper=...> or <...>
        if (i <= n && tokens$value[i] == "<") {
            i = skip_brackets(tokens, i, "<", ">")
        }

        # Optional: dimensions [...]
        if (i <= n && tokens$value[i] == "[") {
            i = skip_brackets(tokens, i, "[", "]")
        }

        # Now we should be at the variable name (an identifier)
        if (i <= n && tokens$type[i] == "IDENTIFIER") {
            vars = c(vars, tokens$value[i])
            i = i + 1
            # Old-style trailing dimensions: real y[J];
            if (i <= n && tokens$value[i] == "[") {
                i = skip_brackets(tokens, i, "[", "]")
            }
        }

        # Skip rest of declaration (possible assignment, more declarators) to semicolon
        i = skip_to_semicolon(tokens, i, depth)
    }

    vars
}


# Parse sampling statements from a block's token stream.
# Returns a list of R formulas.
# Only extracts statements at brace depth 0 (top-level within the block).
parse_sampling = function(tokens) {
    if (nrow(tokens) == 0) return(list())

    samps = list()
    i = 1
    n = nrow(tokens)
    depth = 0

    while (i <= n) {
        if (tokens$value[i] == "{") { depth = depth + 1; i = i + 1; next }
        if (tokens$value[i] == "}") { depth = depth - 1; i = i + 1; next }

        if (depth != 0) { i = i + 1; next }

        # Pattern 1: expr ~ dist(args);
        # Look for TILDE at top level
        tilde_pos = find_token_in_stmt(tokens, i, "TILDE")
        semi_pos = find_token_in_stmt(tokens, i, "SEMICOLON")

        if (!is.na(tilde_pos) && !is.na(semi_pos) && tilde_pos < semi_pos) {
            samp = parse_tilde_stmt(tokens, i, tilde_pos, semi_pos)
            if (!is.null(samp)) samps = c(samps, list(samp))
            i = semi_pos + 1
            next
        }

        # Pattern 2: target += dist_lpdf(expr | args);  or  dist_lpmf(...)
        if (!is.na(semi_pos) && tokens$value[i] == "target" &&
            i + 1 <= n && tokens$value[i + 1] == "+=") {
            samp = parse_target_stmt(tokens, i + 2, semi_pos)
            if (!is.null(samp)) samps = c(samps, list(samp))
            i = semi_pos + 1
            next
        }

        # Not a sampling statement; skip to next semicolon
        if (!is.na(semi_pos)) {
            i = semi_pos + 1
        } else {
            i = i + 1
        }
    }

    samps
}


# Parse `lhs ~ dist(args)` into an R formula
parse_tilde_stmt = function(tokens, start, tilde_pos, semi_pos) {
    # LHS: tokens from start to tilde_pos - 1
    lhs_str = paste(tokens$value[start:(tilde_pos - 1)], collapse = "")
    # Strip any trailing indexing from LHS for the formula variable name
    lhs_var = extract_base_var(lhs_str)

    # RHS: tokens from tilde_pos + 1 to semi_pos - 1
    rhs_tokens = tokens[(tilde_pos + 1):(semi_pos - 1), , drop = FALSE]

    # Handle optional truncation T[...] at the end -- just strip it
    rhs_tokens = strip_truncation(rhs_tokens)

    rhs_str = stan_expr_to_r(rhs_tokens)

    make_formula(lhs_var, rhs_str)
}


# Parse `target += dist_lpdf(expr | args)` into an R formula
# Returns NULL if not a recognized distribution call
parse_target_stmt = function(tokens, start, semi_pos) {
    if (start >= semi_pos) return(NULL)

    expr_tokens = tokens[start:(semi_pos - 1), , drop = FALSE]
    # Look for a function call with _lpdf, _lpmf, _lupdf, _lupmf suffix
    # at the top level of the expression
    if (nrow(expr_tokens) < 1) return(NULL)

    # The expression might be just `dist_lpdf(x | args)` or could have
    # additional terms. We only handle the simple case.
    func_name = expr_tokens$value[1]
    if (expr_tokens$type[1] != "IDENTIFIER") return(NULL)

    suffixes = c("_lpdf", "_lpmf", "_lupdf", "_lupmf", "_lpmf", "_log")
    matched_suffix = NULL
    for (sfx in suffixes) {
        if (endsWith(func_name, sfx)) {
            matched_suffix = sfx
            break
        }
    }
    if (is.null(matched_suffix)) return(NULL)

    dist_name = sub(paste0(matched_suffix, "$"), "", func_name)

    # Expect ( after function name
    if (nrow(expr_tokens) < 3 || expr_tokens$value[2] != "(") return(NULL)

    # Find the matching )
    close_paren = find_matching(expr_tokens, 2, "(", ")")
    if (is.na(close_paren)) return(NULL)

    # Tokens inside parens
    inner = expr_tokens[3:(close_paren - 1), , drop = FALSE]
    if (nrow(inner) == 0) return(NULL)

    # Find the | separator (at top level within these parens)
    bar_pos = find_bar_in_tokens(inner)

    if (!is.na(bar_pos)) {
        # New syntax: dist_lpdf(x | args)
        lhs_tokens = inner[1:(bar_pos - 1), , drop = FALSE]
        if (bar_pos < nrow(inner)) {
            args_tokens = inner[(bar_pos + 1):nrow(inner), , drop = FALSE]
        } else {
            args_tokens = data.frame(type = character(0), value = character(0),
                                     stringsAsFactors = FALSE)
        }
    } else if (matched_suffix == "_log") {
        # Old syntax: dist_log(x, args) -- first arg is the variate
        # Find first top-level comma
        comma_pos = find_comma_in_tokens(inner)
        if (is.na(comma_pos)) {
            lhs_tokens = inner
            args_tokens = data.frame(type = character(0), value = character(0),
                                     stringsAsFactors = FALSE)
        } else {
            lhs_tokens = inner[1:(comma_pos - 1), , drop = FALSE]
            args_tokens = inner[(comma_pos + 1):nrow(inner), , drop = FALSE]
        }
    } else {
        # No | found in lpdf/lpmf -- might be single-arg distribution
        lhs_tokens = inner
        args_tokens = data.frame(type = character(0), value = character(0),
                                 stringsAsFactors = FALSE)
    }

    lhs_str = paste(lhs_tokens$value, collapse = "")
    lhs_var = extract_base_var(lhs_str)

    if (nrow(args_tokens) > 0) {
        rhs_str = paste0(dist_name, "(", stan_expr_to_r(args_tokens), ")")
    } else {
        rhs_str = paste0(dist_name, "()")
    }

    make_formula(lhs_var, rhs_str)
}


# --- Helper functions ---

# Skip balanced brackets starting at position i.
# Returns position after the closing bracket.
skip_brackets = function(tokens, i, open, close) {
    n = nrow(tokens)
    if (i > n || tokens$value[i] != open) return(i)
    depth = 0
    while (i <= n) {
        if (tokens$value[i] == open) depth = depth + 1
        else if (tokens$value[i] == close) depth = depth - 1
        i = i + 1
        if (depth == 0) break
    }
    i
}

# Skip to semicolon at the given brace depth, return position after it
skip_to_semicolon = function(tokens, i, target_depth) {
    n = nrow(tokens)
    depth = target_depth
    while (i <= n) {
        if (tokens$value[i] == "{") depth = depth + 1
        else if (tokens$value[i] == "}") depth = depth - 1
        else if (tokens$value[i] == ";" && depth == target_depth) return(i + 1)
        i = i + 1
    }
    i
}

# Find a token of given type in the current statement (up to next ; at depth 0)
# Only matches at paren/bracket depth 0 within the statement
find_token_in_stmt = function(tokens, start, token_type) {
    n = nrow(tokens)
    depth = 0  # paren/bracket depth
    for (i in start:n) {
        v = tokens$value[i]
        if (v %in% c("(", "[")) depth = depth + 1
        else if (v %in% c(")", "]")) depth = depth - 1
        if (depth == 0 && tokens$type[i] == token_type) return(i)
        if (depth == 0 && v == ";" ) return(if (token_type == "SEMICOLON") i else NA)
        if (v == "{" || v == "}") return(NA)  # hit nested block
    }
    NA
}

# Find matching close bracket
find_matching = function(tokens, open_pos, open, close) {
    n = nrow(tokens)
    depth = 0
    for (i in open_pos:n) {
        if (tokens$value[i] == open) depth = depth + 1
        else if (tokens$value[i] == close) depth = depth - 1
        if (depth == 0) return(i)
    }
    NA
}

# Find | at top level in a token data frame
find_bar_in_tokens = function(tokens) {
    depth = 0
    for (i in seq_len(nrow(tokens))) {
        v = tokens$value[i]
        if (v %in% c("(", "[")) depth = depth + 1
        else if (v %in% c(")", "]")) depth = depth - 1
        else if (v == "|" && depth == 0) return(i)
    }
    NA
}

# Find first comma at top level in a token data frame
find_comma_in_tokens = function(tokens) {
    depth = 0
    for (i in seq_len(nrow(tokens))) {
        v = tokens$value[i]
        if (v %in% c("(", "[")) depth = depth + 1
        else if (v %in% c(")", "]")) depth = depth - 1
        else if (v == "," && depth == 0) return(i)
    }
    NA
}

# Extract base variable name from an expression like "y[n]" or "Y[n, m]"
extract_base_var = function(expr_str) {
    stringr::str_extract(expr_str, "^[a-zA-Z_][a-zA-Z0-9_]*")
}

# Convert Stan expression tokens to an R-compatible expression string.
# Handles Stan-to-R operator differences.
stan_expr_to_r = function(tokens) {
    values = tokens$value
    # Stan integer division is %/% in R (Stan uses %/% too, but sometimes \)
    paste(values, collapse = "")
}

# Strip trailing T[...] truncation from RHS tokens
strip_truncation = function(tokens) {
    n = nrow(tokens)
    if (n < 1) return(tokens)

    # Check if last token is ] and there's a T before a [
    if (tokens$value[n] != "]") return(tokens)

    # Walk back to find the matching [
    depth = 0
    for (i in n:1) {
        if (tokens$value[i] == "]") depth = depth + 1
        else if (tokens$value[i] == "[") depth = depth - 1
        if (depth == 0) {
            # Check if preceded by T
            if (i > 1 && tokens$value[i - 1] == "T" && tokens$type[i - 1] == "IDENTIFIER") {
                return(tokens[1:(i - 2), , drop = FALSE])
            }
            return(tokens)
        }
    }
    tokens
}

# Create a formula from LHS and RHS strings, or return NULL on failure
make_formula = function(lhs, rhs) {
    stmt = paste(lhs, "~", rhs)
    tryCatch(
        stats::as.formula(stmt, env = rlang::empty_env()),
        error = function(e) NULL
    )
}


# Parse Stan `model_code` into a list with two elements:
#   `vars` — named character vector mapping variable names to their block
#   `samp` — list of sampling statements as R formulas
parse_model = function(model_code) {
    tokens = tokenize_stan(model_code)
    blocks = extract_blocks(tokens)

    if (length(blocks) == 0) return(list(vars = character(0), samps = list()))

    # Extract variable declarations from each block
    vars = character(0)
    for (block_name in names(blocks)) {
        block_tokens = blocks[[block_name]]
        block_vars = parse_declarations(block_tokens)
        if (length(block_vars) > 0) {
            new_vars = rep(block_name, length(block_vars))
            names(new_vars) = block_vars
            vars = c(vars, new_vars)
        }
    }

    # Extract sampling statements from model block
    samps = list()
    if ("model" %in% names(blocks)) {
        samps = parse_sampling(blocks[["model"]])
    }

    # Add implicit uniform priors for parameters without sampling statements
    parameters = names(vars)[vars == "parameters"]
    sampled_pars = purrr::map_chr(samps, ~ deparse(rlang::f_lhs(.)))
    uniform_pars = setdiff(parameters, sampled_pars)
    if (length(uniform_pars) > 0) {
        uniform_samp = purrr::map(uniform_pars, function(p) {
            stats::as.formula(paste0(p, " ~ uniform(-1e100, 1e100)"),
                              env = rlang::empty_env())
        })
    } else {
        uniform_samp = NULL
    }

    list(vars = vars, samp = c(samps, uniform_samp))
}


# Take a list of provided sampling formulas and return a matching list of
# sampling statements from a reference list
match_sampling_stmts = function(new_samp, ref_samp) {
    ref_vars = purrr::map_chr(ref_samp, ~ deparse(rlang::f_lhs(.)))
    new_vars = purrr::map_chr(new_samp, ~ deparse(rlang::f_lhs(.)))
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
    get_ast = function(x) purrr::map_if(as.list(x), rlang::is_call, get_ast)
    if (!rlang::is_call(rlang::f_rhs(stmt)))
        stop("Sampling statment ", format(stmt),
             " does not contain a distribution on the right-hand side.")
    # pull out variables from RHS
    rhs_vars = rlang::call_args(rlang::f_rhs(stmt)) %>%
        get_ast %>%
        unlist %>%
        purrr::discard(is.numeric) %>%
        as.character %>%
        purrr::discard(~ . %in% c("`+`", "`-`", "`*`", "`/`", "`^`", "`%*%`", "`%%`"))
    c(deparse(rlang::f_lhs(stmt)), rhs_vars)
}
