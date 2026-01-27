library(INLA)
library(dlnm)

# Plot functions for INLA vehicle.

inla_plot_random_effect = function(term_name) function(mod, fit)
{
    index = which(names(fit$summary.random) == term_name)
    
    size = fit$size.random[[index]]
    summ = fit$summary.random[[index]]
    group_category = mod$group_lookup[term == term_name, unique(group)]
    group_names = mod$group_lookup[term == term_name, name]
    
    if (length(group_names)) {
        summ$by = rep(group_names, each = size$N)
        ggplot(summ) +
            geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`), fill = "#9999ff", alpha = 0.5) +
            geom_line(aes(x = ID, y = mean), colour = "#9999ff", linewidth = 1) +
            facet_wrap(~by) +
            labs(x = term_name, y = "Effect", title = paste0("Effect of ", term_name, " by ", group_category))

    } else {
        ggplot(summ) +
            geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`), fill = "#9999ff", alpha = 0.5) +
            geom_line(aes(x = ID, y = mean), colour = "#9999ff", linewidth = 1) +
            labs(x = term_name, y = "Effect", title = paste0("Effect of ", term_name))
    }
}

inla_plot_distributed_lag = function(var_cb, name) function(mod, fit)
{
    # Get crossbasis object
    cb = mod$env[[as.character(var_cb)]]
    
    # extract full coef and vcov and create indicators for each term
    coef = fit$summary.fixed$mean
    vcov = fit$misc$lincomb.derived.covariance.matrix
    
    # find position of the terms associated with crossbasis
    indt = grep(paste0("^", name), fit$names.fixed)

    # extract predictions from the DLNM centered on overall mean 
    predt = dlnm::crosspred(cb, coef = coef[indt], vcov = vcov[indt,indt],
                       model.link = "log", bylag = 1, cen = attr(cb, "plot_center")) 

    # contour plot of exposure-lag-response associations
    x = predt$predvar
    z = predt$matRRfit
    y = seq(attr(cb, "plot_min_lag"), attr(cb, "plot_max_lag"), length.out = dim(z)[2])

    pal = rev(RColorBrewer::brewer.pal(11, "PRGn"))
    levels = pretty(z, 20)
    col1 = colorRampPalette(pal[1:6])
    col2 = colorRampPalette(pal[6:11])
    cols = c(col1(sum(levels < 1)), col2(sum(levels > 1)))

    filled.contour(x, y, z,
        xlab = name, ylab = "Lag (days)", main = "",
        col = cols, levels = levels)
}

inla_plot_spatial = function(term_name, graph_name) function(mod, fit)
{
    summ = fit$summary.random[[as.character(term_name)]]
    effects = summ[, c("ID", "mean", "0.025quant", "0.975quant")]
    area_names = mod$env[[as.character(graph_name)]]$names
    
    # Need to sample effects from posterior...
    samp = inla.posterior.sample(1000, fit)
    
    pattern = paste0("^", term_name, ":[0-9]+$")
    phiname = paste0("Phi for ", term_name)
    tauname = paste0("Precision for ", term_name)
    
    idx = grep(pattern, rownames(samp[[1]]$latent))
    idx_u = head(idx, length(idx)/2)
    idx_v = tail(idx, length(idx)/2)

    eta = sapply(samp, function(s) {
        (sqrt(s$hyperpar[phiname]) * s$latent[idx_u] +
        sqrt(1 - s$hyperpar[phiname]) * s$latent[idx_v]) /
        sqrt(s$hyperpar[tauname])
    })
    
    eta = as.data.table(reshape2::melt(eta, varnames = c("id", "samp")))
    summ = eta[, .(mean = mean(value), q025 = quantile(value, 0.025), q975 = quantile(value, 0.975)), by = id]
    summ[, name := area_names[id]]
    
    ids = unique(mod$env[[as.character(term_name)]])
    summ = summ[id %in% ids]
    
    ggplot(summ) +
        geom_pointrange(aes(x = name, y = mean, ymin = q025, ymax = q975), colour = "#9999ff") +
        labs(x = term_name, y = "Effect", title = paste0("Effect of ", term_name)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}

# Definitions for INLA vehicle.

Offset.inla = function(g, x)
{
    term(
        offset(!!x), 
        name = x,
        plot = 0
    )
}

# TODO replicate versus group!!
Trend.inla = function(g, x, by = NULL)
{
    var_hyper = new_var(g, "hyper", x)
    var_by = new_var_by_idx(g, by, x)
    
    add_code(g, {
        !!var_hyper = list(
            prec = list(prior = "pc.prec", param = c(1, 0.01))
        )
    })
    
    term(
        f(!!x, model = "rw1", replicate = !!var_by, 
            scale.model = TRUE, constr = TRUE, hyper = !!var_hyper),
        name = x,
        plot = inla_plot_random_effect(x)
    )
}

Cycle.inla = function(g, x, by = NULL)
{
    var_hyper = new_var(g, "hyper", x)
    var_by = new_var_by_idx(g, by, x)
    
    add_code(g, {
        !!var_hyper = list(
            prec = list(prior = "pc.prec", param = c(1, 0.01))
        )
    })
    
    term(
        f(!!x, model = "rw2", cyclic = TRUE, replicate = !!var_by, 
            scale.model = TRUE, constr = TRUE, hyper = !!var_hyper),
        name = x,
        plot = inla_plot_random_effect(x)
    )
}

Spatial.inla = function(g, x, graph_file)
{
    var_hyper = new_var(g, "hyper", x)
    var_graph = new_var(g, "graph")
    var_area = new_var(g, "area")
    
    add_code(g, {
        !!var_hyper = list(
            theta1 = list(prior = "pc.prec", param = c(1, 0.01)),
            theta2 = list(prior = "pc", param = c(0.5, 0.5))
        )
        !!var_graph = read_graph_inla(!!graph_file)
        !!var_area = match(gen$data[gen$valid, !!x], (!!var_graph)$names)
    })
    
    term(
        f(!!var_area, model = "bym2", graph = (!!var_graph)$mat, scale.model = TRUE, 
            constr = TRUE, adjust.for.con.comp = TRUE, hyper = !!var_hyper),
        name = x,
        plot = inla_plot_spatial(var_area, var_graph)
    )
}

FixedLag.inla = function(g, x, t, lag, window, group = NULL)
{
    var_lagged = new_var(g, x, "lag", lag)
    
    add_code(g, {
        !!var_lagged = get_fixed_lag(gen, !!as.character(x), !!as.character(t), 
            !!lag, !!window, !!as.character(group))
    })
    
    term(
        !!var_lagged, 
        name = x,
        plot = 0
    )
}

DistributedLag.inla = function(g, x, group = NULL, lags = seq(0, 180, by = 30), 
    xknots_q = c(0, 1/3, 2/3, 1), tknots_n = 3)
{
    var_cb = new_var(g, x, "distlag")
    
    add_code(g, {
        !!var_cb = get_dlnm_crossbasis(gen, !!as.character(x), !!as.character(group),
            !!lags, !!xknots_q, !!tknots_n)
    })
    
    term(
        !!var_cb, 
        name = x,
        plot = inla_plot_distributed_lag(var_cb, x)
    )
}



fit_model.inla = function(model)
{
    # fixed effects priors
    control.fixed1 = list(
      mean.intercept = 0, # prior mean for intercept
      prec.intercept = 0.01, # prior precision for intercept
      mean = 0, # prior mean for fixed effects
      prec = 1,  # prior precision for (scaled) fixed effects
      correlation.matrix = TRUE # compute posterior correlation matrix
    )
    
    inla(model$formula,
        verbose = FALSE, data = model$data[model$valid], family = "nbinomial",
        control.fixed = control.fixed1, 
        control.predictor = list(compute = TRUE, link = 1),
        control.compute = list(cpo = TRUE, waic = TRUE, dic = TRUE, config = TRUE, return.marginals = TRUE),
        control.inla = list(
            strategy = "adaptive", # adaptive gaussian
            cmin = 0 # fixing Q factorisation issue https://groups.google.com/g/r-inla-discussion-group/c/hDboQsJ1Mls)
        ),
        inla.mode = "experimental"
    )
}
