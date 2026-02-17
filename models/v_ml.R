library(xgboost)
library(MASS)

# Plot functions for ML vehicle (TO DO)

ml_plot_effect = function(term_name) function(mod, fit)
{
    stop("not implemented")
}

# Definitions for ML vehicle.

Offset.ml = function(g, x)
{
    warning("Offset isn't implemented")
    term(
        !!x, 
        name = x,
        plot = 0
    )
}

# TODO replicate versus group!!
Trend.ml = function(g, x, by = NULL)
{
    warning("No specific implementation for Trend")

    term(
        !!x,
        name = x,
        plot = ml_plot_effect(x)
    )
}

Cycle.ml = function(g, x, by = NULL)
{
    warning("No specific implementation for Cycle")

    term(
        !!x,
        name = x,
        plot = ml_plot_effect(x)
    )
}

Spatial.ml = function(g, x, graph_file)
{
    warning("No specific implementation for Spatial")

    term(
        !!x,
        name = x,
        plot = ml_plot_effect(x)
    )
}

FixedLag.ml = function(g, x, t, lag, window, group = NULL)
{
    var_lagged = new_var(g, x, "lag", lag)
    
    add_code(g, {
        !!var_lagged = get_fixed_lag(gen, !!as.character(x), !!as.character(t), 
            !!lag, !!window, !!as.character(group))
    })
    
    term(
        !!var_lagged, 
        name = x,
        plot = ml_plot_effect(x)
    )
}

DistributedLag.ml = function(g, x, group = NULL, lags = seq(0, 180, by = 30), 
    xknots_q = c(0, 1/3, 2/3, 1), tknots_n = 3)
{
    warning("PROBLEM - No specific implementation for DistributedLag")

    term(
        !!x,
        name = x,
        plot = ml_plot_effect(x)
    )
}


# negative binomial objective
nb_objective = function(preds, dtrain, theta = 5)
{
    y = getinfo(dtrain, "label")
    mu = exp(preds)
    grad = (mu - y) / (mu + theta)
    hess = (mu * (y + theta)) / ((mu + theta)^2)
    return (list(grad = grad, hess = hess))
}


fit_model.ml = function(model, nsamp)
{
    # Get response and covariates
    covariates = all.vars(model$formula[[3]])
    Y = model$Y
    X = matrix(numeric(0), nrow = length(Y), ncol = 0)
    for (var in covariates) {
        if (var %in% names(model$data)) {
            x = model$data[model$valid, get(var)]
        } else if (var %in% ls(model$env)) {
            x = model$env[[var]]
        } else {
            stop("Variable ", var, " not found")
        }
      
        if (is.numeric(x) || is.logical(x)) {
            X = cbind(X, as.numeric(x))
            colnames(X)[ncol(X)] = var
        } else if (is.factor(x)) {
            levs = levels(x)
            x2 = matrix(0, nrow = nrow(X), ncol = length(levs))
            x2[cbind(1:nrow(X), as.integer(x))] = 1
            colnames(x2) = paste0(var, ":", levs)
            X = cbind(X, x2)
        } else if (is.character(x)) {
            levs = unique(x)
            x2 = matrix(0, nrow = nrow(X), ncol = length(levs))
            x2[cbind(1:nrow(X), match(x, levs))] = 1
            colnames(x2) = paste0(var, ":", levs)
            X = cbind(X, x2)
        }
    }
    
    # Split data
    row_train = which(model$set[model$set != "unused"] == "train")
    row_test = which(model$set[model$set != "unused"] == "test")
    
    method = "xgb"
    
    if (method == "xgb") {
        # Create DMatrices
        ddata = xgb.DMatrix(data = X, label = Y)
        dtrain = xgb.DMatrix(data = X[row_train, ], label = Y[row_train])
        dtest = xgb.DMatrix(data = X[row_test, ], label = Y[row_test])

        # Train xgboost model for point estimates
        cat("Fitting...\n")
        fit = xgb.train(
            params = list(max_depth = 4, eta = 0.01), # lower eta seems to improve fit?
            data = dtrain,
            nrounds = 1000,
            objective = function(preds, dtrain) nb_objective(preds, dtrain, theta = 5)
        )

        # # Do bootstrapping for uncertainty on eta
        # nboot = 500
        # eta_samples = matrix(0, nrow = nrow(ddata), ncol = nboot)
        # for (i in 1:nboot) {
        #     # Create DMatrix
        #     row_train1 = sample(row_train, size = length(row_train), replace = TRUE)
        #     dtrain1 = xgb.DMatrix(data = X[row_train1, ], label = Y[row_train1])
        #     
        #     fit1 = xgb.train(
        #         params = list(max_depth = 4, eta = 1, objective = "count:poisson"), # lower eta seems to improve fit?
        #         data = dtrain1,
        #         nrounds = 100
        #     )
        #     
        #     # eta predictions
        #     eta_samples[, i] = predict(fit1, ddata)
        # }
        # mean(apply(eta_samples, MARGIN = 1, sd))


        # Train xgboost model
        cat("Fitting...\n")
        fit = xgb.train(
            params = list(max_depth = 4, eta = 0.01), # lower eta seems to improve fit?
            data = dtrain,
            nrounds = 100,
            objective = function(preds, dtrain) nb_objective(preds, dtrain, theta = 5)
        )
        
        # Build data
        data = model$data
        data = cbind(data, Set = model$set)
        
        # Get prediction
        cat("Response...\n");
        eta = predict(fit, ddata)
        data[model$valid, "Mean"]   = exp(eta)
        data[model$valid, "Median"] = qnbinom(0.025, mu = exp(eta), size = 5)
        data[model$valid, "Lo95"]   = qnbinom(0.500, mu = exp(eta), size = 5)
        data[model$valid, "Hi95"]   = qnbinom(0.975, mu = exp(eta), size = 5)
        
        # Sample trajectories
        cat("Sampling trajectories...\n")
        traj = rbindlist(lapply(1:nsamp, function(s) {
            data.table(eta = eta, theta = 5, id = seq_len(nrow(X)))
        }), idcol = "sample")
        traj[, row := which(model$valid)[id]]
        traj[, set := model$set[row]]
        traj[, y_obs := model$Y[id]]
        traj[, dist := "nb"]
        setcolorder(traj, c("sample", "id", "row", "set", "y_obs", "dist"))
    
        # Sample predictions
        cat("Sampling predictions...\n")
        traj[, y_pred := rnbinom(1:.N, size = theta, mu = exp(eta))]
    
        return (list(fit = fit, data = data[], traj = traj[]))
    }
    return (NULL)
    
    # Sample trajectories
    cat("Sampling trajectories...\n")
    psamp = inla.posterior.sample(nsamp, fit)
    isize = which(names(psamp[[1]]$hyperpar) %like% "^size for (the )?nbinomial")
    irows = which(rownames(psamp[[1]]$latent) %like% "^Predictor:[0-9]+$")
    traj = rbindlist(lapply(psamp, function(s) {
        data.table(eta = unname(s$latent[irows,]), theta = s$hyperpar[[isize]], id = seq_along(irows))
    }), idcol = "sample")
    traj[, row := which(model$valid)[id]]
    traj[, set := model$set[row]]
    traj[, y_obs := model$Y[id]]
    traj[, dist := "nb"]
    setcolorder(traj, c("sample", "id", "row", "set", "y_obs", "dist"))

    # Sample predictions
    cat("Sampling predictions...\n")
    traj[, y_pred := rnbinom(1:.N, size = theta, mu = exp(eta))]

    return (list(fit = fit, data = data[], traj = traj[]))
}
