library(INLA)
library(dlnm)
library(stringr)

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
    
    inla(model,
        verbose = FALSE, data = environment(model)$data_fit, family = "nbinomial",
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
