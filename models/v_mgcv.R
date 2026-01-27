library(mgcv)
library(dlnm)
library(gratia)

# Plot functions for mgcv vehicle.

mgcv_plot_smooth = function(term_name, by_name) function(mod, fit)
{
    smooth_name = paste0("s(", term_name, ")")
    x = as.data.table(smooth_estimates(fit, smooth_name, partial_match = TRUE))

    group_category = mod$group_lookup[term == term_name, unique(group)]
    group_names = mod$group_lookup[term == term_name, name]
    
    setnames(x, as.character(term_name), "xx")
    
    if (is.symbol(by_name)) {
        setnames(x, as.character(by_name), "group")
        ggplot(x) + 
            geom_ribbon(aes(xx, ymin = .estimate - 1.96 * .se, ymax = .estimate + 1.96 * .se),
                fill = "#9999ff", alpha = 0.5) +
            geom_line(aes(xx, .estimate), 
                colour = "#9999ff", linewidth = 1) +
            facet_wrap(~group) +
            labs(x = term_name, y = "Effect", title = paste0("Effect of ", term_name, " by ", group_category))
    } else {
        ggplot(x) + 
            geom_ribbon(aes(xx, ymin = .estimate - 1.96 * .se, ymax = .estimate + 1.96 * .se),
                fill = "#9999ff", alpha = 0.5) +
            geom_line(aes(xx, .estimate), 
                colour = "#9999ff", linewidth = 1) +
            labs(x = term_name, y = "Effect", title = paste0("Effect of ", term_name))
    }
}

mgcv_plot_distributed_lag = function(x, var_xlag, var_tlag) function(mod, fit)
{
    smooth_name = paste0("te(", var_xlag, ",", var_tlag, ")")
    data = as.data.table(smooth_estimates(fit, smooth_name))
    
    setnames(data, c(as.character(var_xlag), as.character(var_tlag)), c("x", "t"))
    
    ggplot(data) +
        geom_raster(aes(x = x, y = t, fill = .estimate)) +
        scale_fill_gradient2() +
        labs(x = x, y = "Lag", title = paste0("Effect of ", x, " by lag"))
}

# mgcv_plot_distributed_lag = function(var_cb, name) function(mod, fit)
# {
#     # Get crossbasis object
#     cb = mod$env[[as.character(var_cb)]]
#     
#     # extract full coef and vcov
#     coef = coef(fit)
#     vcov = vcov(fit)
# 
#     # find position of the terms associated with crossbasis
#     indt = grep(paste0("^", var_cb), names(coef))
# 
#     # extract predictions from the DLNM centered on overall mean 
#     predt = dlnm::crosspred(cb, coef = coef[indt], vcov = vcov[indt,indt],
#                        model.link = "log", bylag = 1, cen = attr(cb, "plot_center")) 
# 
#     # contour plot of exposure-lag-response associations
#     x = predt$predvar
#     z = predt$matRRfit
#     y = seq(attr(cb, "plot_min_lag"), attr(cb, "plot_max_lag"), length.out = dim(z)[2])
# 
#     pal = rev(RColorBrewer::brewer.pal(11, "PRGn"))
#     levels = pretty(z, 20)
#     col1 = colorRampPalette(pal[1:6])
#     col2 = colorRampPalette(pal[6:11])
#     cols = c(col1(sum(levels < 1)), col2(sum(levels > 1)))
# 
#     filled.contour(x, y, z,
#         xlab = name, ylab = "Lag (days)", main = "",
#         col = cols, levels = levels)
# }

mgcv_plot_spatial = function(term_name, graph_name) function(mod, fit)
{
    smooth_name = paste0("s(", term_name, ")")
    extract = as.data.table(predict(fit, type = "terms", terms = smooth_name, se.fit = TRUE))
    extract = cbind(extract, name = mod$env[[as.character(term_name)]])
    extract = unique(extract)
  
    ggplot(extract) +
        geom_pointrange(aes(x = name, y = fit, ymin = fit - 1.96 * se.fit, ymax = fit + 1.96 * se.fit), colour = "#9999ff") +
        labs(x = term_name, y = "Effect", title = paste0("Effect of ", term_name)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}



# Definitions for mgcv vehicle.

Offset.mgcv = function(g, x)
{
    term(offset(!!x), name = x)
}

# TODO replicate versus group!!
Trend.mgcv = function(g, x, by = NULL)
{
    var_by = new_var_by_fct(g, by, x)
    if (is.null(var_by)) var_by = NA # for mgcv

    term(
        s(!!x, bs = "ps", m = 1, by = !!var_by), 
        name = x,
        plot = mgcv_plot_smooth(x, var_by)
    )
}

Cycle.mgcv = function(g, x, by = NULL)
{
    var_by = new_var_by_fct(g, by, x)
    if (is.null(var_by)) var_by = NA # for mgcv
    
    term(
        s(!!x, bs = "cc", m = 2, by = !!var_by), 
        name = x,
        plot = mgcv_plot_smooth(x, var_by)
    )
}

Spatial.mgcv = function(g, x, graph_file)
{
    var_graph = new_var(g, "graph")
    var_area = new_var(g, "area")

    add_code(g, {
        !!var_graph = read_graph_mgcv(!!graph_file, unique(gen$data[gen$valid, !!x]))
        !!var_area = factor(gen$data[gen$valid, !!x], names(!!var_graph))
    })
    
    term(
        s(!!var_area, bs = "mrf", xt = list(nb = !!var_graph)),
        name = x,
        plot = mgcv_plot_spatial(var_area, var_graph)
    )
}

FixedLag.mgcv = function(g, x, t, lag, window, group = NULL)
{
    var_lagged = new_var(g, x, "lag", lag)
    
    add_code(g, {
        !!var_lagged = get_fixed_lag(gen, !!as.character(x), !!as.character(t), 
            !!lag, !!window, !!as.character(group))
    })
    
    term(!!var_lagged, name = x) # TODO these bits!
}

DistributedLag.mgcv = function(g, x, group = NULL, lags = seq(0, 180, by = 30), 
    xknots_q = c(0, 1/3, 2/3, 1), tknots_n = 3)
{
    var_xlag = new_var(g, x, "distlag")
    var_tlag = new_var(g, "t_distlag")
    
    add_code(g, {
        !!var_xlag = get_lagged(gen, !!as.character(x), !!as.character(group), !!lags)
        !!var_tlag = matrix(!!lags, nrow(!!var_xlag), length(!!lags), byrow = TRUE)
    })
    
    term(
        te(!!var_xlag, !!var_tlag, bs = c("ps", "ps")), 
        name = x,
        plot = mgcv_plot_distributed_lag(x, var_xlag, var_tlag)
    )

    # var_cb = new_var(g, x, "distlag")
    # 
    # add_code(g, {
    #     !!var_cb = get_dlnm_crossbasis(gen, !!as.character(x), !!as.character(group),
    #         !!lags, !!xknots_q, !!tknots_n)
    # })
    # 
    # term(
    #     !!var_cb, 
    #     name = x,
    #     plot = mgcv_plot_distributed_lag(var_cb, x)
    # )
}






fit_model.mgcv = function(model)
{
    # This seems to be necessary as mgcv doesn't seem to look in the formula's
    # environment. Consider:
    # formula <- y ~ x
    # environment(formula) <- new.env()
    # environment(formula)$y <- (1:10 + rnorm(10))
    # environment(formula)$x <- 1:10
    # glm(formula)
    # gam(formula)
    #
    # The glm call works, but the gam call doesn't.
    with(model$env,
        gam(formula = formula, family = "nb", data = gen$data[gen$valid])
    )
}
