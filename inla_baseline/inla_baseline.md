# INLA baseline model
Simon Frost

## Libraries

``` r
library(INLA)
```

    Loading required package: Matrix

    This is INLA_25.06.22-1 built 2025-06-22 16:26:02 UTC.
     - See www.r-inla.org/contact-us for how to get help.
     - List available models/likelihoods/etc with inla.list.models()
     - Use inla.doc(<NAME>) to access documentation
     - Consider upgrading R-INLA to testing[25.09.04] or stable[25.06.07].

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

Define formula, using single `spi` and `precip` variables.

``` r
form = formula(y ~ 1 +
                 offset(logpop) + 
                 f(polyid, model="bym2", graph=nb_loc, scale.model=TRUE, constr=TRUE, adjust.for.con.comp=TRUE, hyper=hyper.bym2) +
                 f(year, model="rw1", replicate=polyid2, scale.model=TRUE, constr=TRUE, hyper=hyper3.rw) +
                 f(epiWeek, model="rw2", cyclic=TRUE, replicate=region2, scale.model=TRUE, constr=TRUE, hyper=hyper3.rw) +
                 f(spi, model="rw2", scale.model=TRUE, constr=TRUE, hyper=hyper2.rw) +
                 f(precip, model="rw2", scale.model=TRUE, constr=TRUE, hyper=hyper2.rw))
```

``` r
fit <- inla(form,
            verbose = FALSE,
            data = df2,
            family = "nbinomial",
            control.fixed = control.fixed1, # prior precision for (scaled) fixed effects
            control.predictor=list(compute=TRUE,
                                   link=1),
            control.compute=list(cpo=TRUE,
                                 waic=TRUE,
                                 dic=TRUE,
                                 config=TRUE,
                                 return.marginals=TRUE),
            control.inla = list(strategy='adaptive', # adaptive gaussian
                                cmin=0), # fixing Q factorisation issue https://groups.google.com/g/r-inla-discussion-group/c/hDboQsJ1Mls)
            inla.mode = "experimental")
```


     *** inla.core.safe:  The inla program failed, but will rerun in case better initial values may help. try=1/1 

     *** inla.core.safe:  rerun with improved initial values 

## Model summary

``` r
summary(fit)
```

    Time used:
        Pre = 11.7, Running = 2.94, Post = 0.702, Total = 15.3 
    Fixed effects:
                  mean    sd 0.025quant 0.5quant 0.975quant   mode    kld
    (Intercept) -5.118 0.418      -5.95   -5.117      -4.29 -5.117 44.917

    Random effects:
      Name    Model
        polyid BYM2 model
       year RW1 model
       epiWeek RW2 model
       spi RW2 model
       precip RW2 model

    Model hyperparameters:
                                                              mean      sd
    size for the nbinomial observations (1/overdispersion)   2.523   0.190
    Precision for polyid                                     0.824   0.330
    Phi for polyid                                           0.267   0.214
    Precision for year                                       1.257   0.283
    Precision for epiWeek                                    0.868   0.265
    Precision for spi                                       24.537  19.965
    Precision for precip                                   146.589 278.204
                                                           0.025quant 0.5quant
    size for the nbinomial observations (1/overdispersion)      2.171    2.515
    Precision for polyid                                        0.346    0.769
    Phi for polyid                                              0.016    0.208
    Precision for year                                          0.788    1.228
    Precision for epiWeek                                       0.455    0.833
    Precision for spi                                           4.482   19.051
    Precision for precip                                        9.157   70.340
                                                           0.975quant   mode
    size for the nbinomial observations (1/overdispersion)      2.921  2.497
    Precision for polyid                                        1.621  0.667
    Phi for polyid                                              0.783  0.040
    Precision for year                                          1.896  1.173
    Precision for epiWeek                                       1.487  0.766
    Precision for spi                                          77.349 11.258
    Precision for precip                                      769.913 22.587

    Deviance Information Criterion (DIC) ...............: 9589.80
    Deviance Information Criterion (DIC, saturated) ....: -9930239.95
    Effective number of parameters .....................: 92.45

    Watanabe-Akaike information criterion (WAIC) ...: 9599.77
    Effective number of parameters .................: 97.59

    Marginal log-Likelihood:  -5510.24 
    CPO, PIT is computed 
    Posterior summaries for the linear predictor and the fitted values are computed
    (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
