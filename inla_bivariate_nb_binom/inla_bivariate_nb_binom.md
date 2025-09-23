# INLA bivariate model, unequal variances
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
vars <- c("y","ys","yt", "logpop","polyid","polyid2","region", "region2","year","epiWeek","spi","precip")
df2 <- df[,vars]
```

``` r
n <- nrow(df2)

## 1) Build two stacks, one per likelihood
eff1 <- list(
  k1 = rep(1, n), k2 = rep(0, n),
  polyid = df2$polyid, year = df2$year, epiWeek = df2$epiWeek,
  polyid2 = df2$polyid2, region2 = df2$region2,
  spi = df2$spi, precip = df2$precip,
  outcome_id = rep(1L, n),
  logpop = df2$logpop
)
stk1 <- inla.stack(data = list(y = df2$yt), A = 1, effects = eff1, tag = "lik1")

eff2 <- list(
  k1 = rep(0, n), k2 = rep(1, n),
  polyid = df2$polyid, year = df2$year, epiWeek = df2$epiWeek,
  polyid2 = df2$polyid2, region2 = df2$region2,
  spi = df2$spi, precip = df2$precip,
  outcome_id = rep(2L, n),
  logpop = rep(0, n)         # offset = 0 for PROPORTION (Binomial)
)
stk2 <- inla.stack(data = list(y = df2$y), A = 1, effects = eff2, tag = "lik2")

# Join stacks
stk <- inla.stack(stk1, stk2)

## 2) Make a TWO-COLUMN response aligned with the joint stack
idx1 <- inla.stack.index(stk, tag = "lik1")$data
idx2 <- inla.stack.index(stk, tag = "lik2")$data

dat <- inla.stack.data(stk)       # list with all covariates from stacks
A   <- inla.stack.A(stk)
m   <- nrow(A)                    # total length of linear predictor (should be 2n)

Y <- matrix(NA_real_, nrow = m, ncol = 2)
Y[idx1, 1] <- df2$yt               # first likelihood's rows in col 1
Y[idx2, 2] <- df2$y              # second likelihood's rows in col 2
dat$Y <- Y                        # put 2-col response into data

## --- Ntrials only for the Binomial column (second likelihood) -------------
Ntrials <- matrix(0, nrow = m, ncol = 2)
Ntrials[idx2, 2] <- df2$yt
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

# prior on NB parameters
pc.prec1 <- list(hyper = list(theta = list(prior = "pc.prec", param = c(1, 0.01))))
```

Define formula, using single `spi` and `precip` variables.

``` r
form <- Y ~ 0 + k1 + k2 + offset(logpop) +
  f(polyid,  model="bym2", graph=nb_loc, scale.model=TRUE, constr=TRUE,
    adjust.for.con.comp=TRUE, hyper=hyper.bym2,
    group = outcome_id, control.group = control.group1) +
  f(year,    model="rw1", replicate=polyid2, scale.model=TRUE, constr=TRUE, hyper=hyper3.rw,
    group = outcome_id, control.group = control.group1) +
  f(epiWeek, model="rw2", cyclic=TRUE, replicate=region2, scale.model=TRUE, constr=TRUE, hyper=hyper3.rw,
    group = outcome_id, control.group = control.group1) +
  f(spi,     model="rw2", scale.model=TRUE, constr=TRUE, hyper=hyper2.rw,
    group = outcome_id, control.group = control.group1) +
  f(precip,  model="rw2", scale.model=TRUE, constr=TRUE, hyper=hyper2.rw,
    group = outcome_id, control.group = control.group1)
```

``` r
fit <- inla(form,
            verbose = FALSE,
            data = dat,
            family = c("nbinomial","binomial"),
            control.family = list(pc.prec1, list()),
            #Ntrials = Ntrials, # <- only second likelihood uses it
            Ntrials = df2$yt,
            control.fixed = control.fixed1,
            control.predictor = list(compute = TRUE, A = A),
            control.compute = list(cpo=TRUE, waic=TRUE, dic=TRUE,
                                   config=TRUE, return.marginals=TRUE,
                                   return.marginals.predictor=TRUE),
            control.inla = list(strategy="adaptive", cmin=0),
            inla.mode = "experimental")
```

## Model summary

``` r
summary(fit)
```

    Time used:
        Pre = 15.5, Running = 18.3, Post = 1.53, Total = 35.4 
    Fixed effects:
         mean    sd 0.025quant 0.5quant 0.975quant   mode    kld
    k1 -2.576 0.318     -3.198   -2.577     -1.948 -2.577 20.007
    k2 -1.776 0.322     -2.405   -1.777     -1.142 -1.777 10.332

    Random effects:
      Name    Model
        polyid BYM2 model
       year RW1 model
       epiWeek RW2 model
       spi RW2 model
       precip RW2 model

    Model hyperparameters:
                                                              mean       sd
    size for the nbinomial observations (1/overdispersion)   1.351    0.052
    Precision for polyid                                     1.181    0.421
    Phi for polyid                                           0.190    0.105
    GroupRho for polyid                                     -0.004    0.113
    Precision for year                                       1.318    0.206
    GroupRho for year                                       -0.251    0.149
    Precision for epiWeek                                    0.905    0.232
    GroupRho for epiWeek                                     0.343    0.292
    Precision for spi                                      630.139 4380.641
    GroupRho for spi                                         0.032    0.133
    Precision for precip                                    79.560  132.169
    GroupRho for precip                                     -0.068    0.167
                                                           0.025quant 0.5quant
    size for the nbinomial observations (1/overdispersion)      1.252    1.350
    Precision for polyid                                        0.551    1.116
    Phi for polyid                                              0.052    0.167
    GroupRho for polyid                                        -0.230   -0.003
    Precision for year                                          0.971    1.298
    GroupRho for year                                          -0.528   -0.256
    Precision for epiWeek                                       0.532    0.877
    GroupRho for epiWeek                                       -0.287    0.373
    Precision for spi                                           2.280   84.539
    GroupRho for spi                                           -0.220    0.029
    Precision for precip                                        3.752   40.607
    GroupRho for precip                                        -0.411   -0.062
                                                           0.975quant   mode
    size for the nbinomial observations (1/overdispersion)      1.455  1.349
    Precision for polyid                                        2.186  0.997
    Phi for polyid                                              0.456  0.121
    GroupRho for polyid                                         0.213  0.004
    Precision for year                                          1.782  1.251
    GroupRho for year                                           0.055 -0.264
    Precision for epiWeek                                       1.438  0.824
    GroupRho for epiWeek                                        0.818  0.457
    Precision for spi                                        4350.443  3.093
    GroupRho for spi                                            0.302  0.013
    Precision for precip                                      397.272  9.522
    GroupRho for precip                                         0.236 -0.021

    Deviance Information Criterion (DIC) ...............: 29127.62
    Deviance Information Criterion (DIC, saturated) ....: 5809.52
    Effective number of parameters .....................: 180.02

    Watanabe-Akaike information criterion (WAIC) ...: 29491.02
    Effective number of parameters .................: 406.26

    Marginal log-Likelihood:  -16014.93 
    CPO, PIT is computed 
    Posterior summaries for the linear predictor and the fitted values are computed
    (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
