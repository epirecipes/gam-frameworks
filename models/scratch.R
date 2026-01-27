# TODO in lassa_sentinel/output/covariates/spatial_model, what is this data?
# what is the "mastomys" indicator here? very interesting
# TODO spi1 contains -Inf
# TODO fine lag seems to outperform coarse lag - makes sense (smoothing)
# TODO knot placement along lag -- 0, 1 too narrow I think

# What do I want to be able to do?
# Plot / summarize effects on linear predictor
# Evaluate predictive accuracy
# Validate predictive accuracy, with differing model terms and hyperparameters

library(data.table)
library(ggplot2)

source("./models/model.R")
source("./models/v_inla.R")
source("./models/v_mgcv.R")

x = source("./models/data.R", local = new.env())$value
data = x$df
nb_loc = "./data/states_nb.csv" # x$nb_loc

formula = y ~
    Offset(logpop) + 
    Trend(year, by = adminName) +
    Cycle(epiWeek, by = region) + 
    Spatial(adminName, graph_file = nb_loc) +
    FixedLag(precip, t = date, lag = 90, window = 60, group = adminName) +
    DistributedLag(spi3, group = adminName, lags = seq(0, 150, by = 1),
      tknots_n = 5, xknots_q = seq(0, 1, by = 0.2))

# INLA
set.seed(123)
mod = make_model(formula, "inla", data)
mod
fit = fit_model(mod)
mod$terms
mod$terms$year$plot(mod, fit)
mod$terms$epiWeek$plot(mod, fit)
mod$terms$spi3$plot(mod, fit)
mod$terms$adminName$plot(mod, fit)

summary(fit) # WAIC 7988
plot(fit)

# MGCV
set.seed(123)
mod = make_model(formula, "mgcv", data)
mod
fit = fit_model(mod)
mod$terms
mod$terms$year$plot(mod, fit)
mod$terms$epiWeek$plot(mod, fit)
mod$terms$spi3$plot(mod, fit)
mod$terms$adminName$plot(mod, fit)

bloop = smooth_estimates(fit, "te(", partial_match = TRUE)
ggplot(bloop) +
    geom_raster(aes(x = spi3_distlag, y = t_distlag, fill = .estimate)) +
    scale_fill_gradient2()


xxx = read_graph_mgcv(file)
mod$data$adminName2 = factor(mod$data$adminName, names(xxx))
mod$data$adminName3 = as.integer(mod$data$adminName2)
table(mod$data$adminName3)
yyy = read_graph_inla(file)
xxx = lapply(xxx, function(x) match(x, names(xxx)))
names(xxx) = NULL
mod$env$nb = xxx

gam(y ~ s(adminName3, bs = "mrf", xt = list(nb = mod$env$nb)), family = "nb", data = mod$data)


ls(x)
x$formula
x$make
# term = function(term, ..., .env = parent.frame())
# {
#     list(
#         term = do.call(rlang::expr, list(substitute(term)), envir = .env),
#         extras = list(...)
#     )
# }
# 
# fit_model = function(model)
# {
#     fit = match.fun(paste0("fit_model.", environment(model)$vehicle))
#     fit(model)
# }



# example data
library(data.table)
library(mgcv)

dt = data.table(t = seq(0, 2 * pi, length.out = 20))
dt[, y := sin(t) + rnorm(20, 0, 0.4)]
plot(dt)

y = dt$y
t = dt$t
yy = dt$y
tt = dt$t

?s

g = new_generator("miscellaneous", dt)
bs = "tp" # "tp"
var_bs = new_var(g, "my_bs")
add_code(g, { !!var_bs = !!bs })
set_lhs(g, quote(y))
add_rhs(g, quote(s(t, bs = my_bs)))
finalize(g)

test = g$make()
test
ls(environment(test))
environment(test)$my_bs
class(environment(test)$my_bs)

plot(gam(formula = test, data = dt))
plot(gam(formula = test))

library(INLA)
summary(inla(yy ~ tt, data = dt))



# Old old old 

source("./models/inla/inla_model.R")
source("./models/inla/inla_terms.R")
source("./models/model.R")
library(ggplot2)

x = source("./models/data.R", local = new.env())$value
data = x$df
nb_loc = x$nb_loc

formula = y ~ 1 + 
    Offset(logpop) + 
    Spatial(polyid, graph_file = nb_loc) +
    Trend(year, by = adminName) +
    Cycle(epiWeek, by = region) #+ 
    # FixedLag(precip, t = date, lag = 90, window = 60, group = adminName)
    DistributedLag(precip, group = adminName, lags = seq(30, 150, by = 1)) +
    DistributedLag(spi3, group = adminName, lags = seq(30, 150, by = 1))
    # FixedLag(spi3, t = date, lag = 90, window = 60, group = adminName)
    # spi3
    # DistributedLag(spi12, group = adminName)
    # DistributedLag(spi1, group = adminName)

set.seed(123)
mod = make_model(formula, "inla", data)

# TODO big problem here, prior is not working with !! approach
# Note the difference between these two: is this why?
environment(mod)$terms$year$term[[7]]
environment(mod)$terms$epiWeek$term[[8]]



fit = fit_model(mod)
environment(mod)$terms$epiWeek$plot(mod, fit)

environment(mod)$terms$year$plot(mod, fit)
environment(mod)$terms$epiWeek$plot(mod, fit)
environment(mod)$terms$precip$plot(mod, fit)
environment(mod)$terms$spi3$plot(mod, fit)

summary(fit) # WAIC = 7967.23
plot(fit)

x = inla.posterior.sample(1000, fit)
xx = inla.posterior.sample.eval(function() return (Predictor), x)



# prior test
max = 144
set.seed(123)
dat = data.table(time = 1:max, temp = 20 + 10 * sin(2*pi*(1:max)/12) + rnorm(max))

proportion = 0.5
num = round(max * proportion)
dat$time[sample.int(max, num)] = dat$time[sample.int(max, num)]
dat[, month := (time - 1) %% 12]

fit = inla(temp ~ f(month, model = "rw2", cyclic = TRUE, scale.model = TRUE, constr = TRUE, 
  hyper = list(theta = list(prior = "pc.prec", param = c(0.8, 0.01)))), data = dat)
plot(fit, plot.random.effects = TRUE, plot.hyperparameters = FALSE, plot.predictor = FALSE, plot.fixed.effects = FALSE,
  plot.lincomb = FALSE, plot.q = FALSE)
plot(dat$time, dat$temp)

fit = inla(y ~ f(epiWeek, model = "rw2", replicate = region2, cyclic = TRUE, scale.model = TRUE, constr = TRUE, hyper = list("prec" = list(prior = "pc.prec", param = c(0.01, 0.01)))), 
  data = data)
plot(fit, plot.random.effects = TRUE, plot.hyperparameters = FALSE, plot.predictor = FALSE, plot.fixed.effects = FALSE,
  plot.lincomb = FALSE, plot.q = FALSE)



model = fit
cb = environment(mod)$cb[[1]]
pattern = "^precip"
name = "precip"
cen = 0



pdf("contour2.pdf", width = 6.5, height = 5, onefile = TRUE)
plot_lag_surface(fit, cb_spi3, "^spi3", "SPI-3", cen = 0)
plot_lag_surface(fit, cb_precip, "^precip", "Precipitation (mm)", cen = 7)


# prediction loop

cutoffs = lubridate::make_date(year = rep(2019:2024, each = 2), month = c(1, 7))

results = list()
for (cut in as.character(cutoffs)) {
  print(cut)
  data2 = copy(data)
  data2[, y_orig := y]
  data2[date >= cut, y := NA]
  mod = make_model(formula, "inla", data2)
  fit = fit_model(mod)
  data2[, y_mean := fit$summary.fitted.values$mean]
  
  results[[cut]] = data2[date >= cut & !is.na(y_orig), .(adminCode, date, y_mean, y_orig, forecast_date = as.Date(cut))]
}
res = rbindlist(results)

ggplot(res) +
  geom_line(aes(x = date, y = y_mean, colour = as.factor(forecast_date))) +
  geom_point(aes(x = date, y = y_orig)) + 
  facet_wrap(~adminCode) + scale_y_log10()

ggplot(res[forecast_date >= "2022-01-01"]) +
  geom_line(aes(x = date, y = y_mean, colour = as.factor(forecast_date))) +
  geom_point(aes(x = date, y = y_orig)) + 
  facet_wrap(~adminCode)

# ok let's compare precip_0 from data to climate.
library(data.table)
library(ggplot2)
library(zoo)

d2 = as.data.table(data)
data0 = d2[, .(adminCode, date, precip_0)]

climate0 = climate[, .(adminCode, date, precip)]
climate0[, precip_c0 := frollmean(precip, n = 60, align = "right", na.rm = TRUE), by = adminCode]

ddd = merge(data0, climate0, by = c("adminCode", "date"))
ggplot(ddd) + geom_point(aes(x = precip_0, y = precip_c0))

lm(precip_0 ~ precip_c0, data = ddd)



# experiment with random forest
data_rf = data[!is.na(y)]
data_rf[, y_rate := (y + 0.5) / exp(logpop)]

data_rf2 = data_rf[date < "2023-10-01"]

# train a default random forest model
rf1 = ranger::ranger(
  y_rate ~ adminCode + year + epiWeek + precip + spi3, 
  data = data_rf2,
  mtry = 2,
  respect.unordered.factors = "order",
  seed = 123
)

x = predict(rf1, data = data_rf)
data_rf[, pred := x$predictions]

library(ggplot2)

ggplot(data_rf) +
  geom_line(aes(x = date, y = y_rate)) +
  geom_line(aes(x = date, y = pred, colour = "prediction")) +
  facet_wrap(~adminCode)

# get OOB RMSE
(default_rmse <- sqrt(ames_rf1$prediction.error))
## [1] 24859.27
