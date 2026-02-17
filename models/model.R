# alll right
# so, the formula along with the vehicle is like the spec
# it should generate a function

# Utility parts

new_generator = function(vehicle, data, formula_orig)
{
    # A good sanity check and to avoid any problems with new_var().
    if (anyDuplicated(names(data))) {
        stop("error: duplicated names in data.")
    }
    if (any(names(data) == "gen")) {
        stop("data contains a column `gen`, which is reserved for internal use. Please rename this column.")
    }
  
    rlang::new_environment(list(
        vehicle = vehicle,
        data = data.table::as.data.table(data),
        make = function(gen) {},
        formula_orig = formula_orig,
        terms = list()
    ))
}

add_code = function(gen, code)
{
    code = rlang::quo_get_expr(rlang::enquo(code))
    blen = length(body(gen$make))
    
    # Add statement(s)
    if (is.call(code) && identical(code[[1]], quote(`{`))) {
        elen = length(code) - 1
        body(gen$make)[(blen + 1):(blen + elen)] = code[-1]
    } else {
        body(gen$make)[[blen + 1]] = code
    }
}

# Create a variable name with no clashes in gen$make's body or among 
# gen$data's column names.
new_var = function(gen, ..., sep = "_")
{
    # Assemble variable name
    varname = do.call(paste, c(lapply(list(...), as.character), list(sep = sep)))
    
    # Assemble variable names already used...
    existing = c(
        all.vars(body(gen$make)), # ... in body of generator
        names(formals(gen$make)), # ... as argument to generator
        names(gen$data)           # ... as a column in data
    )
    # Apply make_unique to all; should only modify (at most) "varname" entry.
    uniqued = make.unique(c(existing, varname), sep)
    
    # Return unique-ified "varname" entry.
    return (as.name(uniqued[length(uniqued)]))
}

new_var_by_idx = function(gen, by, x)
{
    if (!is.null(by)) {
        var = new_var(gen, "by", by)
        add_code(gen, {
            !!var = gen$data[gen$valid, match(!!by, unique(!!by))]
            gen$group_lookup = rbind(gen$group_lookup,
                unique(gen$data[gen$valid,
                  .(term = !!as.character(x),
                    group = !!as.character(by),
                    name = !!by,
                    index = !!as.name(paste0("..", var))
                )])
            )
        })
        
        return (var)
    } else {
        return (NULL)
    }
}

new_var_by_fct = function(gen, by, x)
{
    if (!is.null(by)) {
        var = new_var(gen, "by", by)
        add_code(gen, {
            !!var = gen$data[gen$valid, as.factor(!!by)]
            gen$group_lookup = rbind(gen$group_lookup,
                unique(gen$data[gen$valid,
                  .(term = !!as.character(x),
                    group = !!as.character(by),
                    name = !!by,
                    index = NA_integer_
                )])
            )
        })
        return (var)
    } else {
        return (NULL)
    }
}

finalize = function(gen, formula)
{
    add_code(gen, { 
        # remove all names starting with dot
        rm(list = setdiff(ls(all.names = TRUE), ls(all.names = FALSE))); 
        
        # return formula with function environment explicitly attached.
        formula = !!formula
        environment(formula) = environment()
        return (formula)
    })
}

term = function(term, name, ..., .env = parent.frame())
{
    list(
        term = do.call(rlang::expr, list(substitute(term)), envir = .env),
        extras = c(list(name = as.character(name)), list(...))
    )
}

# Model parts
make_model = function(formula, vehicle, train, test = NULL)
{
    # Create generator
    g = new_generator(vehicle, rbind(train, test), formula)
    
    # Keep track of train versus test data
    g$set = factor(
        c(rep("train", nrow(train)), rep("test", if (is.null(test)) 0 else nrow(test))), 
        levels = c("train", "test", "unused")
    )
    setDT(g$data) ## TODO remove?

    # Get outcome variable
    # The rlang::duplicate is needed here due to a perversity of data.table.
    y = rlang::duplicate(eval(formula[[2]], g$data, environment(formula)))
    g$Y = y[!is.na(y)]
    var_Y = new_var(g, "Y")
    add_code(g, {
        !!var_Y = rlang::duplicate(eval(quote(!!formula[[2]]), gen$data, environment(gen$formula)))[gen$valid]
    })
    formula[[2]] = var_Y
    
    # Get valid data indices (non-NA)
    g$valid = !is.na(y)
    stopifnot(length(g$valid) == nrow(g$data))
    g$set[!g$valid] = "unused"
    g$weights = case_match(g$set, "train" ~ 1, "test" ~ 0, "unused" ~ NA_real_)

    # Data and term substitutions
    # Operate on right side of formula y ~ x1 + x2 + ...
    formula[[3]] = elixir::expr_apply(formula[[3]], function(x) {
        # Look for functions starting with a capital letter
        m = elixir::expr_match(x, ~{ `.F:name/^[A-Z]`(...A) });
        if (!is.null(m[[1]])) {
            # Substitute with Funcname.vehicle(args)
            func = match.fun(paste0(m[[1]]$F, ".", vehicle))
            processed = do.call(func, args = c(list(g), m[[1]]$A), quote = TRUE, envir = environment(formula))
            
            # Add terms & extras
            if (!is.null(processed$extras$name)) {
                nm = tail(make.unique(c(names(g$terms), as.character(processed$extras$name)), "_"), 1)
                g$terms[[nm]] = c(list(term = processed$term), processed$extras)
            }

            return (processed$term)
        }
        return (x)
    }, into = TRUE)
    
    # Add final code to make func
    finalize(g, formula)
    
    # Generate formula and environment
    g$formula = g$make(g)
    
    # Save reference to environment
    g$env = environment(g$formula)
    
    # Add class
    class(g) = "generic_model"
    
    # g is now the model object.
    return (g)
}

pretty_formula = function(formula)
{
    cat(paste0(gsub("\\+", "\\+\n   ", deparse1(formula)), "\n"))
}

print.generic_model = function(x)
{
    cat("Generic model:\n")
    pretty_formula(x$formula_orig)
    
    cat('\n"', x$vehicle, '" formula:\n', sep = "")
    pretty_formula(x$formula)

    cat("\nComponents", paste0("$", ls(x), collapse = ", "), "\n")

    cat("\n$data columns:", paste0(names(x$data), collapse = ", "), "\n")
    cat("$env contents:", paste0(ls(x$env), collapse = ", "), "\n")
}

fit_model = function(model, nsamp = 1000)
{
    fit = match.fun(paste0("fit_model.", model$vehicle))
    fit(model, nsamp)
}

# Spatial data helpers
# Read in a csv file with two columns, name and adjoining.
# name is the name of each area, repeated for as many rows as it has neighbours,
# and adjoining is the corresponding set of neighbours.
# e.g. for countries of the UK and their land borders, we would have:
# name, adjoining
# England, Scotland
# England, Wales
# Northern Ireland, NA   <- note there is an entry here with NA
# Scotland, England
# Wales, England
read_graph_mgcv = function(file, keep = NULL)
{
    # Load CSV
    info = fread(file)
    nb_list = split(info[, adjoining], info[, name])
    
    # Put in alphabetical order & omit NAs
    nb_list = lapply(nb_list, \(x) sort(x[!is.na(x)]))
    
    # Sanity check
    invalid = setdiff(unique(unlist(nb_list)), names(nb_list))
    unattached = setdiff(names(nb_list), unique(unlist(nb_list)))
    if (length(invalid)) {
        stop(file, " -- Adjoining areas not in name column: ", paste(invalid, collapse = ", "))
    }
    if (length(unattached)) {
        warning(file, " -- Areas with no neighbours: ", paste(unattached, collapse = ", "))
    }
    
    if (!is.null(keep)) {
        nb_list = nb_list[keep]
    }

    return (nb_list)
}

read_graph_inla = function(file)
{
    nb_list = read_graph_mgcv(file)
    
    # Convert to matrix
    mat = matrix(0L, nrow = length(nb_list), ncol = length(nb_list))
    for (i in seq_along(nb_list)) {
        mat[i, match(nb_list[[i]], names(nb_list))] <- 1
    }

    return (list(names = names(nb_list), mat = mat))
}

get_fixed_lag = function(gen, x_nm, t_nm, lag, window, group_nm)
{
    data = gen$data[, .(t = get(t_nm), x = get(x_nm), 
        group = if (length(group_nm)) get(group_nm) else 0)]
    
    # Ensure data are all sequential
    iota_check = data[, .(sequential = all(diff(t) == 1)), by = group]
    if (any(!iota_check$sequential)) {
        print(iota_check[])
        stop("Non-time-sequential entries found in data.")
    }
    
    # Get lagged data
    lagged = data[, frollmean(shift(x, n = lag), n = window, 
        align = "right", na.rm = TRUE), by = group][gen$valid, V1]
    
    return (lagged)
}

get_lagged = function(gen, x_nm, group_nm, lags)
{
    # Create lag data matrix
    data = gen$data[, .(x = get(x_nm), 
        group = if (length(group_nm)) get(group_nm) else 0)]

    lagged = NULL
    for (lag in lags) {
        lagged = cbind(lagged, data[, .(shift(x, n = lag)), by = group]$V1)
    }

    return (lagged[gen$valid, ])
}

get_dlnm_crossbasis = function(gen, x_nm, group_nm,
    lags, xknots_q = c(0, 1/3, 2/3, 1), tknots_n = 3)
{
    lagged = get_lagged(gen, x_nm, group_nm, lags)

    # Inner knots and boundary points along predictor dimension
    xpoints = quantile(lagged, xknots_q, na.rm = TRUE)
    xboundary = c(xpoints[1], tail(xpoints, 1))
    xknots = head(tail(xpoints, -1), -1)

    # Knots along time dimension: this is in units of lags, not of time per se.
    # TODO old version here or new version below?
    # tk = ncol(lagged)^(1/tknots_n)
    # tknots = tk^(0:(tknots_n - 1))
    tknots = seq(0, ncol(lagged) - 1, length.out = tknots_n + 2)
    tknots = tknots[2:(tknots_n + 1)]

    # Create crossbasis
    cb = crossbasis(lagged,
       argvar = list(fun = "ns", knots = xknots, Boundary.knots = xboundary),
       arglag = list(fun = "ns", knots = tknots, intercept = TRUE)
    )
    colnames(cb) = stringr::str_replace(colnames(cb), "^v", x_nm)
    
    # For plotting
    attr(cb, "plot_center") = quantile(lagged, 0.5, na.rm = TRUE)
    attr(cb, "plot_min_lag") = min(lags)
    attr(cb, "plot_max_lag") = max(lags)

    return (cb)
}
