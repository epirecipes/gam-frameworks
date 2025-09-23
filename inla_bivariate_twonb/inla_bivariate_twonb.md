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
vars <- c("y","ys","logpop","polyid","polyid2","region", "region2","year","epiWeek","spi","precip")
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
stk1 <- inla.stack(data = list(y = df2$y), A = 1, effects = eff1, tag = "lik1")

eff2 <- list(
  k1 = rep(0, n), k2 = rep(1, n),
  polyid = df2$polyid, year = df2$year, epiWeek = df2$epiWeek,
  polyid2 = df2$polyid2, region2 = df2$region2,
  spi = df2$spi, precip = df2$precip,
  outcome_id = rep(2L, n),
  logpop = df2$logpop
)
stk2 <- inla.stack(data = list(y = df2$ys), A = 1, effects = eff2, tag = "lik2")

# Join stacks
stk <- inla.stack(stk1, stk2)

## 2) Make a TWO-COLUMN response aligned with the joint stack
idx1 <- inla.stack.index(stk, tag = "lik1")$data
idx2 <- inla.stack.index(stk, tag = "lik2")$data

dat <- inla.stack.data(stk)       # list with all covariates from stacks
A   <- inla.stack.A(stk)
m   <- nrow(A)                    # total length of linear predictor (should be 2n)

Y <- matrix(NA_real_, nrow = m, ncol = 2)
Y[idx1, 1] <- df2$y               # first likelihood's rows in col 1
Y[idx2, 2] <- df2$ys              # second likelihood's rows in col 2
dat$Y <- Y                        # put 2-col response into data
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
            family = c("nbinomial","nbinomial"),
            control.family = list(pc.prec1, pc.prec1),
            control.fixed = control.fixed1,
            control.predictor = list(compute = TRUE, A = A),
            control.compute = list(cpo=TRUE, waic=TRUE, dic=TRUE,
                                   config=TRUE, return.marginals=TRUE,
                                   return.marginals.predictor=TRUE),
            control.inla = list(strategy="adaptive", cmin=0),
            inla.mode = "experimental")
```


     *** inla.core.safe:  The inla program failed, but will rerun in case better initial values may help. try=1/1 

     *** inla.core.safe:  rerun with improved initial values 

## Model summary

``` r
summary(fit)
```

    Time used:
        Pre = 13.2, Running = 23.3, Post = 0.863, Total = 37.4 
    Fixed effects:
         mean    sd 0.025quant 0.5quant 0.975quant   mode    kld
    k1 -3.721 0.454     -4.597   -3.726     -2.813 -3.725 21.014
    k2 -1.709 0.456     -2.591   -1.714     -0.797 -1.712  5.001

    Random effects:
      Name    Model
        polyid BYM2 model
       year RW1 model
       epiWeek RW2 model
       spi RW2 model
       precip RW2 model

    Model hyperparameters:
                                                                 mean      sd
    size for the nbinomial observations (1/overdispersion)      2.514   0.189
    size for the nbinomial observations (1/overdispersion)[2]   1.324   0.052
    Precision for polyid                                        0.604   0.584
    Phi for polyid                                              0.217   0.121
    GroupRho for polyid                                         0.918   0.130
    Precision for year                                          0.950   0.184
    GroupRho for year                                           0.941   0.034
    Precision for epiWeek                                       0.944   0.284
    GroupRho for epiWeek                                        0.805   0.102
    Precision for spi                                          11.563   6.006
    GroupRho for spi                                            0.045   0.076
    Precision for precip                                      109.818 110.526
    GroupRho for precip                                         0.362   0.486
                                                              0.025quant 0.5quant
    size for the nbinomial observations (1/overdispersion)         2.166    2.506
    size for the nbinomial observations (1/overdispersion)[2]      1.223    1.323
    Precision for polyid                                           0.081    0.434
    Phi for polyid                                                 0.064    0.189
    GroupRho for polyid                                            0.541    0.964
    Precision for year                                             0.631    0.935
    GroupRho for year                                              0.854    0.949
    Precision for epiWeek                                          0.492    0.909
    GroupRho for epiWeek                                           0.549    0.826
    Precision for spi                                              3.508   10.358
    GroupRho for spi                                              -0.080    0.038
    Precision for precip                                          16.357   77.424
    GroupRho for precip                                           -0.667    0.473
                                                              0.975quant   mode
    size for the nbinomial observations (1/overdispersion)         2.911  2.487
    size for the nbinomial observations (1/overdispersion)[2]      1.428  1.322
    Precision for polyid                                           2.150  0.214
    Phi for polyid                                                 0.528  0.131
    GroupRho for polyid                                            0.998  0.996
    Precision for year                                             1.351  0.912
    GroupRho for year                                              0.983  0.962
    Precision for epiWeek                                          1.598  0.845
    GroupRho for epiWeek                                           0.940  0.865
    Precision for spi                                             26.488  8.105
    GroupRho for spi                                               0.212  0.005
    Precision for precip                                         402.014 40.443
    GroupRho for precip                                            0.979  0.972

    Deviance Information Criterion (DIC) ...............: 29306.40
    Deviance Information Criterion (DIC, saturated) ....: -10511128.76
    Effective number of parameters .....................: 131.41

    Watanabe-Akaike information criterion (WAIC) ...: 29949.32
    Effective number of parameters .................: 468.28

    Marginal log-Likelihood:  -16100.77 
    CPO, PIT is computed 
    Posterior summaries for the linear predictor and the fitted values are computed
    (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
