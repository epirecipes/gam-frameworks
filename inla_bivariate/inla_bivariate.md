# INLA bivariate model
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
     - Consider upgrading R-INLA to testing[25.09.19] or stable[25.06.07].

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
vars <- c("y","ys","logpop","polyid","polyid2","region", "region2","year","epiWeek","spi","precip")
df2 <- df[,vars]
```

``` r
df_long <- rbind(
  transform(df2, y_bi = y,  k = 1),
  transform(df2, y_bi = ys, k = 2)
)
df_long$k <- factor(df_long$k)
df_long$k2 <- as.integer(df_long$k)
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

# prior on correlations
control.group1 = list(model="exchangeable",hyper=list(rho=list(prior="pc.cor0", param=c(0.5,0.5))))
```

Define formula, using single `spi` and `precip` variables.

``` r
# Outcome-specific intercepts:
# 0 + k yields two intercepts (one per outcome)
form <- formula(
  y_bi ~ 0 + k +
    offset(logpop) +
    f(polyid,  model="bym2", graph=nb_loc, scale.model=TRUE, constr=TRUE,
      adjust.for.con.comp=TRUE, hyper=hyper.bym2,
      group = k2, control.group = control.group1) +
    f(year,    model="rw1", replicate=polyid2, scale.model=TRUE, constr=TRUE, hyper=hyper3.rw,
      group = k2, control.group = control.group1) +
    f(epiWeek, model="rw2", cyclic=TRUE, replicate=region2, scale.model=TRUE, constr=TRUE, hyper=hyper3.rw,
      group = k2, control.group = control.group1) +
    f(spi,     model="rw2", scale.model=TRUE, constr=TRUE, hyper=hyper2.rw,
      group = k2, control.group = control.group1) +
    f(precip,  model="rw2", scale.model=TRUE, constr=TRUE, hyper=hyper2.rw,
      group = k2, control.group = control.group1)
)
```

``` r
fit <- inla(form,
            verbose = FALSE,
            data = df_long,
            family = "nbinomial",
            control.family = list(hyper=list(theta=list(prior="pc.prec", param=c(1,0.01)))),  # shared dispersion
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
        Pre = 12, Running = 29, Post = 2.68, Total = 43.7 
    Fixed effects:
         mean    sd 0.025quant 0.5quant 0.975quant   mode    kld
    k1 -3.688 0.541     -4.664   -3.718     -2.528 -3.797 15.426
    k2 -1.709 0.541     -2.679   -1.739     -0.546 -1.818  3.651

    Random effects:
      Name    Model
        polyid BYM2 model
       year RW1 model
       epiWeek RW2 model
       spi RW2 model
       precip RW2 model

    Model hyperparameters:
                                                             mean     sd 0.025quant
    size for the nbinomial observations (1/overdispersion)  1.549  0.052      1.448
    Precision for polyid                                    0.483  0.204      0.200
    Phi for polyid                                          0.220  0.053      0.127
    GroupRho for polyid                                     0.950  0.039      0.846
    Precision for year                                      0.947  0.158      0.671
    GroupRho for year                                       0.950  0.024      0.888
    Precision for epiWeek                                   0.960  0.218      0.598
    GroupRho for epiWeek                                    0.798  0.075      0.620
    Precision for spi                                      14.136  6.208      5.694
    GroupRho for spi                                        0.040  0.069     -0.074
    Precision for precip                                   90.014 55.439     24.853
    GroupRho for precip                                     0.190  0.231     -0.284
                                                           0.5quant 0.975quant
    size for the nbinomial observations (1/overdispersion)    1.548      1.655
    Precision for polyid                                      0.445      0.987
    Phi for polyid                                            0.217      0.333
    GroupRho for polyid                                       0.960      0.990
    Precision for year                                        0.935      1.291
    GroupRho for year                                         0.955      0.982
    Precision for epiWeek                                     0.938      1.453
    GroupRho for epiWeek                                      0.809      0.911
    Precision for spi                                        12.933     29.618
    GroupRho for spi                                          0.034      0.193
    Precision for precip                                     76.696    233.878
    GroupRho for precip                                       0.200      0.606
                                                             mode
    size for the nbinomial observations (1/overdispersion)  1.546
    Precision for polyid                                    0.378
    Phi for polyid                                          0.211
    GroupRho for polyid                                     0.976
    Precision for year                                      0.914
    GroupRho for year                                       0.963
    Precision for epiWeek                                   0.896
    GroupRho for epiWeek                                    0.831
    Precision for spi                                      10.834
    GroupRho for spi                                        0.003
    Precision for precip                                   55.521
    GroupRho for precip                                     0.225

    Deviance Information Criterion (DIC) ...............: 29371.61
    Deviance Information Criterion (DIC, saturated) ....: -3767912.64
    Effective number of parameters .....................: 127.38

    Watanabe-Akaike information criterion (WAIC) ...: 30137.36
    Effective number of parameters .................: 526.51

    Marginal log-Likelihood:  -16130.69 
    CPO, PIT is computed 
    Posterior summaries for the linear predictor and the fitted values are computed
    (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
