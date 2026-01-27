library(elixir)
library(data.table)

# TODO get rid of extras (if not needed)

make_model = function(formula, vehicle, data, .env = parent.frame())
{
    model = formula
    
    environment(model) = rlang::new_environment(parent = .env)
    environment(model)$data_full = as.data.table(data)
    environment(model)$vehicle = vehicle
    environment(model)$terms = list()
    
    # Get valid data indices (non-NA)
    environment(model)$indices = environment(model)$data_full[, !is.na(get(model[[2]]))]
    environment(model)$data_fit = environment(model)$data_full[environment(model)$indices]
    
    # Data and term substitutions
    # Operate on right side of formula y ~ x1 + ...
    j = 0
    model[[3]] = expr_apply(model[[3]], function(x) {
        # Look for functions starting with a capital letter
        m = expr_match(x, ~{ `.F:name/^[A-Z]`(...A) });
        if (!is.null(m[[1]])) {
            j <<- j + 1
            # Substitute with Funcname.vehicle(args)
            func = match.fun(paste0(m[[1]]$F, ".", vehicle))
            processed = do.call(func, args = c(list(environment(model)), m[[1]]$A), quote = TRUE, envir = .env)
            
            nm = as.character(processed$extras$name)
            while (is.null(nm) || nm %in% names(environment(model)$terms)) {
                nm = paste0(nm, j)
            }
            environment(model)$terms[[nm]] = c(list(term = processed$term), processed$extras)
            
            return (processed$term)
        }
        return (x)
    }, into = TRUE)
    
    return (model)
}

term = function(term, ..., .env = parent.frame())
{
    list(
        term = do.call(rlang::expr, list(substitute(term)), envir = .env),
        extras = list(...)
    )
}

fit_model = function(model)
{
    fit = match.fun(paste0("fit_model.", environment(model)$vehicle))
    fit(model)
}
