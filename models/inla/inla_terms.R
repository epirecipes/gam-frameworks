inla_by_var = function(env, by, x)
{
    if (!is.null(by)) {
        by_var = paste0(by, "_IDX")
        if (!by_var %in% names(env$data_full)) {
            env$data_fit[, (by_var) := match(get(by), unique(get(by)))]
            env$data_full[[by_var]][env$indices] = env$data_fit[[by_var]]
            env$group_lookup = rbind(env$group_lookup,
                unique(env$data_fit[, .(term = as.character(x), 
                    group = as.character(by), name = get(by), index = get(by_var))])
            )
        }
        return (as.name(by_var))
    } else {
        return (NULL)
    }
}

inla_plot_random_effect = function(term_name) function(mod, fit)
{
    index = which(names(fit$summary.random) == term_name)
    
    env = environment(mod)
    size = fit$size.random[[index]]
    summ = fit$summary.random[[index]]
    group_category = env$group_lookup[term == term_name, unique(group)]
    group_names = env$group_lookup[term == term_name, name]
    
    summ$by = rep(group_names, each = size$N)
    
    ggplot(summ) +
        geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`), fill = "#9999ff", alpha = 0.5) +
        geom_line(aes(x = ID, y = mean), colour = "#9999ff", linewidth = 1) +
        facet_wrap(~by) +
        labs(x = term_name, y = "Effect", title = paste0("Effect of ", term_name, " by ", group_category))
}

inla_plot_distributed_lag = function(cb, pattern, name, cen, minlag, maxlag) function(mod, fit)
{
    # extract full coef and vcov and create indicators for each term
    coef <- fit$summary.fixed$mean
    vcov <- fit$misc$lincomb.derived.covariance.matrix
    
    # find position of the terms associated with crossbasis
    indt <- grep(pattern, fit$names.fixed)

    # extract predictions from the DLNM centred on overall mean 
    predt <- dlnm::crosspred(cb, coef = coef[indt], vcov = vcov[indt,indt],
                       model.link = "log", bylag = 1, cen = cen) 
    
    # contour plot of exposure-lag-response associations
    x <- predt$predvar
    z <- predt$matRRfit
    y <- seq(minlag, maxlag, length.out = dim(z)[2])

    pal <- rev(RColorBrewer::brewer.pal(11, "PRGn"))
    levels <- pretty(z, 20)
    col1 <- colorRampPalette(pal[1:6])
    col2 <- colorRampPalette(pal[6:11])
    cols <- c(col1(sum(levels < 1)), col2(sum(levels > 1)))

    filled.contour(x, y, z,
        xlab = name, ylab = "Lag (days)", main = "",
        col = cols, levels = levels)
}


Offset.inla = function(env, x)
{
    term(offset(!!x),
        name = x)
}

# TODO replicate versus group!!
Trend.inla = function(env, x, by = NULL)
{
    hyper3.rw = quote(list(
        prec = list(prior = "pc.prec", param = c(1, 0.01))
    ))
    
    by = inla_by_var(env, by, x)

    term(f(!!x, model = "rw1", replicate = !!by, 
        scale.model = TRUE, constr = TRUE, hyper = !!hyper3.rw),
        name = x, plot = inla_plot_random_effect(x))
}

Cycle.inla = function(env, x, by = NULL)
{
    hyper3.rw = quote(list(
        prec = list(prior = "pc.prec", param = c(1, 0.01))
    ))
    
    by = inla_by_var(env, by, x)

    term(f(!!x, model = "rw2", cyclic = TRUE, replicate = !!by, 
        scale.model = TRUE, constr = TRUE, hyper = !!hyper3.rw),
        name = x, plot = inla_plot_random_effect(x))
}

Spatial.inla = function(env, x, graph_file)
{
    hyper.bym2 = quote(list(
        theta1 = list(prior = "pc.prec", param = c(1, 0.01)),
        theta2 = list(prior = "pc", param = c(0.5, 0.5))
    ))

    term(f(!!x, model = "bym2", graph = !!graph_file, scale.model = TRUE, 
        constr = TRUE, adjust.for.con.comp = TRUE, hyper = !!hyper.bym2),
        name = x)
}

FixedLag.inla = function(env, x, t, lag, window, group = NULL)
{
    # TODO put into shared code
    new_var = paste0(x, "_lag", lag, "_w", window)
    rlang::inject({
        if (is.null(group)) {
            iota_check = env$data_full[, .(sequential = all(diff(!!t) == 1))]
        } else {
            iota_check = env$data_full[, .(sequential = all(diff(!!t) == 1)), by = !!group]
        }
        if (any(!iota_check$sequential)) {
            print(iota_check[])
            stop("Non-time-sequential entries found in data.")
        }
        env$data_full[, (new_var) := frollmean(shift(!!x, n = !!lag), n = !!window, 
            align = "right", na.rm = TRUE), by = !!group]
        env$data_fit[[new_var]] = env$data_full[[new_var]][env$indices]
    })
    term(!!as.name(new_var),
        name = x)
}

DistributedLag.inla = function(env, x, group = NULL,
  lags = seq(0, 180, by = 30),
  xknots_q = c(0, 1/3, 2/3, 1), tknots_n = 3)
{
    # Create lag data matrix
    data = NULL
    
    rlang::inject({
        for (lag in !!lags) {
            new_var = paste0(x, "_lag", lag)
            data = cbind(data, env$data_full[, .(`@x` = shift(!!x, n = lag)), by = !!group]$`@x`)
        }
    })
    
    data = data[env$indices, ]

    # Inner knots and boundary points along predictor dimension
    xpoints = quantile(data, xknots_q, na.rm = TRUE)
    xboundary = c(xpoints[1], tail(xpoints, 1))
    xknots = head(tail(xpoints, -1), -1)
    center = quantile(data, 0.5, na.rm = TRUE)

    # Knots along time dimension: this is in units of lags, not of time per se.
    tk = ncol(data)^(1/tknots_n)
    tknots = tk^(0:(tknots_n - 1))
    
    # Create crossbasis
    cb = crossbasis(data,
       argvar = list(fun = "ns", knots = xknots, Boundary.knots = xboundary),
       arglag = list(fun = "ns", knots = tknots, intercept = TRUE)
    )
    colnames(cb) = str_replace(colnames(cb), "^v", as.character(x))
    
    # Save crossbasis in env
    env$cb = c(env$cb, list(cb))
    cb_i = length(env$cb)
    
    min_lag = rlang::inject(min(!!lags))
    max_lag = rlang::inject(max(!!lags))
    
    term(cb[[!!cb_i]],
        name = x,
        plot = inla_plot_distributed_lag(cb, paste0("^", x), x, center, min_lag, max_lag))
}
