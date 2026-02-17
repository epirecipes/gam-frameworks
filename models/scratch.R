# TODO in lassa_sentinel/output/covariates/spatial_model, what is this data?
# what is the "mastomys" indicator here? very interesting
# TODO spi1 contains -Inf
# TODO fine lag seems to outperform coarse lag - makes sense (smoothing)
# TODO knot placement along lag -- 0, 1 too narrow I think

# What do I want to be able to do?
# Plot / summarize effects on linear predictor
# Evaluate predictive accuracy
# Validate predictive accuracy, with differing model terms and hyperparameters

# For INLA, "zeroinflatednbinomial0" seems to fit better than nbinomial.

library(data.table)
library(ggplot2)

source("./models/model.R")
source("./models/score.R")
source("./models/v_inla.R")
source("./models/v_mgcv.R")
source("./models/v_ml.R")

x = source("./models/data.R", local = new.env())$value
data = x$df
nb_loc = "./data/states_nb.csv" # x$nb_loc

# TEMP: remove data with irregular spacing (because no climate data...)
data = data[date <= "2025-10-17"]

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


# Looking at prediction
schema = fread(
"train_begin,  train_end, test_begin,   test_end
2018-01-01,   2019-10-01, 2019-10-01, 2020-10-01
2018-01-01,   2020-10-01, 2020-10-01, 2021-10-01
2018-01-01,   2021-10-01, 2021-10-01, 2022-10-01
2018-01-01,   2022-10-01, 2022-10-01, 2023-10-01
2018-01-01,   2023-10-01, 2023-10-01, 2024-10-01
2018-01-01,   2024-10-01, 2024-10-01, 2025-10-01")

# to do...
# prediction - yes
# scores - working on it
# random forest
# make consistent the treatment of incomplete cases

# this makes little difference
data[!is.na(year), year := lubridate::year(date + 92)]

# testing this
data[!is.na(year) & year == 2021, covid := TRUE]
data[!is.na(year) & year != 2021, covid := FALSE]

# TODO
# make offset proper for ML model
# test bart for count data (?)

formula = y ~
    Offset(logpop) + 
    covid +
    Trend(year, by = adminName) +
    Cycle(epiWeek, by = region) + 
    Spatial(adminName, graph_file = nb_loc) +
    FixedLag(precip, t = date, lag = 90, window = 60, group = adminName) +
    DistributedLag(spi3, group = adminName, lags = seq(0, 150, by = 1),
      tknots_n = 5, xknots_q = seq(0, 1, by = 0.2))

set.seed(123)

results = list()
presults = list()

for (vehicle in c("mgcv", "inla", "ml")) {
    for (r in seq_len(nrow(schema))) {
        cat("Prediction", r, "...\n")
        
        data_train = data[date >= schema[r, train_begin] & date < schema[r, train_end]]
        data_test  = data[date >= schema[r, test_begin] & date < schema[r, test_end]]
        model = make_model(formula, vehicle, data_train, data_test)
        fit = fit_model(model)
        
        LS = fit$traj[set == "test" & !is.na(eta), log_score_nb(y_obs[1], exp(eta), theta), by = id][, mean(V1)]
        # QS = fit$traj[set == "test" & !is.na(eta), quad_score_nb(y_obs[1], exp(eta), theta), by = id][, mean(V1)]
        # RS = fit$traj[set == "test" & !is.na(eta), rp_score_nb(y_obs[1], exp(eta), theta), by = id][, mean(V1)]
        results[[paste0(vehicle, r)]] = data.table(LS = LS)
        
        # Plotting data
        pdata = fit$traj[, .(
          lo = quantile(y_pred, 0.025, na.rm = TRUE), 
          av = mean(y_pred, na.rm = TRUE),
          hi = quantile(y_pred, 0.975, na.rm = TRUE)), 
          by = .(id, row, set, y_obs)]
        
        pdata = merge(pdata, fit$data[, .(adminCode, date, adminName, row = 1:.N)], by = "row")
        presults[[paste0(vehicle, r)]] = pdata
    }
}

res1 = rbindlist(results, idcol = "part")
res2 = rbindlist(presults, idcol = "part")

res2[, vehicle := stringr::str_remove(part, "[0-9]+$")]
res2[, stage := as.integer(stringr::str_extract(part, "[0-9]+$"))]

ggplot(res2[set == "test"]) + 
    geom_ribbon(aes(x = date, ymin = lo, ymax = hi, fill = stage, group = stage), alpha = 0.5) +
    geom_point(aes(x = date, y = y_obs), size = 0.2) +
    geom_line(aes(x = date, y = av, colour = set)) +
    coord_cartesian(ylim = c(0, 150)) +
    facet_grid(vehicle ~ adminName)

ggplot(res2[vehicle == "ml"]) +
    geom_ribbon(aes(x = date, ymin = lo, ymax = hi, fill = set, group = set), alpha = 0.5) +
    geom_point(aes(x = date, y = y_obs), size = 0.2) +
    geom_line(aes(x = date, y = av, colour = set)) +
    coord_cartesian(ylim = c(0, 150)) +
    facet_grid(part ~ adminName)
  


res2[, vehicle := stringr::str_remove(variable, "[0-9]+$")]
res2[, stage := as.integer(stringr::str_extract(variable, "[0-9]+$"))]

res = melt(as.data.table(results), id.vars = character(0))
res[, vehicle := stringr::str_remove(variable, "[0-9]+$")]
res[, stage := as.integer(stringr::str_extract(variable, "[0-9]+$"))]

pdata = fit$traj[, .(
  lo = quantile(y_pred, 0.025, na.rm = TRUE), 
  av = mean(y_pred, na.rm = TRUE),
  hi = quantile(y_pred, 0.975, na.rm = TRUE)), 
  by = .(id, row, set, y_obs)]

pdata = merge(pdata, fit$data[, .(adminCode, date, adminName, row = 1:.N)], by = "row")

ggplot(pdata) + 
  geom_ribbon(aes(x = date, ymin = lo, ymax = hi, fill = set), alpha = 0.5) +
  geom_point(aes(x = date, y = y_obs), size = 0.2) +
  geom_line(aes(x = date, y = av, colour = set)) +
  facet_wrap(~adminName)


#









## Quantile regression forest
rf2 = ranger::ranger(
    y_rate ~ adminCode + year + epiWeek + precip + spi3, 
    data = data_rf2,
    mtry = 2,
    respect.unordered.factors = "order",
    seed = 123,
    quantreg = TRUE,
    keep.inbag = TRUE # this changes quantile regression - how?
)
pred = predict(rf2, data = data_rf, type = "quantiles", quantiles = c(0.025, 0.5, 0.975, 0.99999))
data_rf$pred_lo = pred$predictions[, 1]
data_rf$pred_md = pred$predictions[, 2]
data_rf$pred_hi = pred$predictions[, 3]

ggplot(data_rf[adminCode == "NG029"]) +
    geom_line(aes(x = date, y = y_rate)) +
    geom_line(aes(x = date, y = pred, colour = "prediction")) +
    geom_line(aes(x = date, y = pred_lo, colour = "quantiles")) +
    geom_line(aes(x = date, y = pred_md, colour = "quantiles")) +
    geom_line(aes(x = date, y = pred_hi, colour = "quantiles")) +
    facet_wrap(~adminCode)
