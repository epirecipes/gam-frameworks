# inlabru baseline model
Simon Frost

## Libraries

``` r
library(INLA)
```

    Loading required package: Matrix

    This is INLA_25.09.19 built 2025-09-19 08:26:27 UTC.
     - See www.r-inla.org/contact-us for how to get help.
     - List available models/likelihoods/etc with inla.list.models()
     - Use inla.doc(<NAME>) to access documentation
     - Consider upgrading R-INLA to testing[25.10.19] or stable[25.06.07].

``` r
library(inlabru)
```

    Loading required package: fmesher

``` r
library(ggplot2)
```

``` r
set.seed(123)  # for reproducibility
```

## Load data and process

``` r
df <- read.table("../data/lassa_states_endemic.tsv", header=TRUE, row.names=NULL)
# Exclude missing (forecast) datapoints
df <- df[!is.na(df$y),]
nb_loc <- "../data/states_nbmatrix"
```

``` r
vars <- c("y","logpop","polyid","polyid2","region", "region2","year","epiWeek","spi","precip")
df2 <- df[,vars]
```

## Fit model

Priors.

``` r
# fixed effects priors
control.fixed1 = list(mean.intercept=0, # prior mean for intercept
                      prec.intercept=0.01, # prior precision for intercept
                      mean=0, # prior mean for fixed effects
                      prec=1)  # prior precision for (scaled) fixed effects

# hyperpriors for iid, ar1, bym2
hyper.iid = list(theta = list(prior="pc.prec", param=c(1, 0.01)))

# for bym or bym2
hyper.bym = list(theta1 = list(prior="pc.prec", param=c(1, 0.01)),
                 theta2 = list(prior="pc.prec", param=c(1, 0.01)))
hyper.bym2 = list(theta1 = list(prior="pc.prec", param=c(1, 0.01)),
                  theta2 = list(prior="pc", param=c(0.5, 0.5)))

# ar1
hyper.ar1 = list(theta1 = list(prior='pc.prec', param=c(1, 0.01)), 
                 rho = list(prior='pc.cor0', param = c(0.75, 0.5)))

# random walk hyperparameters on smoothing
hyper1.rw = list(prec = list(prior='pc.prec', param=c(0.1, 0.01))) # strictest; sd constrained to be low
hyper2.rw = list(prec = list(prior='pc.prec', param=c(0.5, 0.01))) # medium
hyper3.rw = list(prec = list(prior='pc.prec', param=c(1, 0.01))) # weakest; sd can be quite wide 
```

Components mirror the INLA f() terms.

- `Intercept(1)` gives an explicit intercept component; its prior is set
  via control.fixed
- For BYM2, we pass the same neighborhood graph file and options as in
  INLA

``` r
cmp <- ~
  Intercept(1) +
  spatial(polyid, model = "bym2",
          graph = nb_loc,
          scale.model = TRUE, constr = TRUE, adjust.for.con.comp = TRUE,
          hyper = hyper.bym2) +
  year(year, model = "rw1",
       replicate = polyid2,
       scale.model = TRUE, constr = TRUE,
       hyper = hyper3.rw) +
  epiWeek(epiWeek, model = "rw2",
          cyclic = TRUE,
          replicate = region2,
          scale.model = TRUE, constr = TRUE,
          hyper = hyper3.rw) +
  spi(spi, model = "rw2",
      scale.model = TRUE, constr = TRUE,
      hyper = hyper2.rw) +
  precip(precip, model = "rw2",
         scale.model = TRUE, constr = TRUE,
         hyper = hyper2.rw)
```

This is our observation model, that refers to the above components by
the names we have assigned.

``` r
obs <- bru_obs(
  formula = y ~ Intercept + offset(logpop) + spatial + year + epiWeek + spi + precip,
  family = "nbinomial",
  data = df2
)
```

``` r
fit <- bru(
  cmp,
  obs,
  options = list(
    control.fixed = control.fixed1,
    control.predictor = list(compute = TRUE, link = 1),
    control.compute = list(cpo = TRUE, waic = TRUE, dic = TRUE,
                           config = TRUE, return.marginals = TRUE),
    control.inla = list(strategy = "adaptive", cmin = 0)  # same workaround you used
  )
)
```

## Model summary

``` r
summary(fit)
```

    inlabru version: 2.13.0
    INLA version: 25.09.19
    Components:
    Intercept: main = linear(1), group = exchangeable(1L), replicate = iid(1L), NULL
    spatial: main = bym2(polyid), group = exchangeable(1L), replicate = iid(1L), NULL
    year: main = rw1(year), group = exchangeable(1L), replicate = iid(polyid2), NULL
    epiWeek: main = rw2(epiWeek), group = exchangeable(1L), replicate = iid(region2), NULL
    spi: main = rw2(spi), group = exchangeable(1L), replicate = iid(1L), NULL
    precip: main = rw2(precip), group = exchangeable(1L), replicate = iid(1L), NULL
    Observation models:
      Family: 'nbinomial'
        Tag: <No tag>
        Data class: 'data.frame'
        Response class: 'integer'
        Predictor: 
            y ~ Intercept + offset(logpop) + spatial + year + epiWeek + spi + 
                precip
        Additive/Linear: FALSE/FALSE
        Used components: effects[Intercept, spatial, year, epiWeek, spi, precip], latent[]
    Time used:
        Pre = 14.3, Running = 1.55, Post = 0.94, Total = 16.8 
    Fixed effects:
                mean   sd 0.025quant 0.5quant 0.975quant   mode kld
    Intercept -5.126 0.42     -5.962   -5.124     -4.297 -5.124   0

    Random effects:
      Name    Model
        spatial BYM2 model
       year RW1 model
       epiWeek RW2 model
       spi RW2 model
       precip RW2 model

    Model hyperparameters:
                                                              mean      sd
    size for the nbinomial observations (1/overdispersion)   2.523   0.190
    Precision for spatial                                    0.825   0.332
    Phi for spatial                                          0.267   0.216
    Precision for year                                       1.256   0.283
    Precision for epiWeek                                    0.869   0.265
    Precision for spi                                       24.483  19.803
    Precision for precip                                   146.868 278.607
                                                           0.025quant 0.5quant
    size for the nbinomial observations (1/overdispersion)      2.172    2.515
    Precision for spatial                                       0.347    0.769
    Phi for spatial                                             0.015    0.206
    Precision for year                                          0.786    1.227
    Precision for epiWeek                                       0.456    0.832
    Precision for spi                                           4.435   19.062
    Precision for precip                                        9.129   70.470
                                                           0.975quant   mode
    size for the nbinomial observations (1/overdispersion)       2.92  2.496
    Precision for spatial                                        1.63  0.666
    Phi for spatial                                              0.79  0.038
    Precision for year                                           1.89  1.172
    Precision for epiWeek                                        1.49  0.766
    Precision for spi                                           76.82 11.235
    Precision for precip                                       771.37 22.544

    Deviance Information Criterion (DIC) ...............: 9588.55
    Deviance Information Criterion (DIC, saturated) ....: -9852522.29
    Effective number of parameters .....................: 92.13

    Watanabe-Akaike information criterion (WAIC) ...: 9598.70
    Effective number of parameters .................: 97.46

    Marginal log-Likelihood:  -5449.49 
    CPO, PIT is computed 
    Posterior summaries for the linear predictor and the fitted values are computed
    (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
