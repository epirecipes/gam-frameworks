# prepare data table in expected format.
x = source("./models/data.R", local = new.env())$value
df = x$df
df = df[, .(y, date, epiWeek, year, precip, spi3, logpop, region_id = factor(adminName, unique(adminName)))]
df = df[!is.na(y)]

# Cut off last bit.
future_df = df
df = df[date <= "2023-11-15"]
future_df = future_df[date > "2023-11-15" & date < "2024-05-12"]

# Below is from AI.

library(dplyr)
library(lubridate)
library(mgcv)
library(purrr)
library(tidyr)

make_lags <- function(dat, lags = 1:4) {
  dat <- dat %>%
    arrange(region_id, date) %>%
    group_by(region_id)
  for (L in lags) {
    dat <- dat %>% mutate(!!paste0("lag_y_", L) := dplyr::lag(y, L))
  }
  dat %>% ungroup()
}

# Option A (iterated / recursive): 1-step model rolled forward
## 1) Fit an autoregressive NB-GAM

lags <- 1:4

df_feat <- df %>% make_lags(lags)

train <- df_feat %>%
  filter(if_all(all_of(paste0("lag_y_", lags)), ~ !is.na(.x)))

# NB-GAM with cyclic seasonality + covariates + AR lags
fit_A <- gam(
  y ~
    lag_y_1 + lag_y_2 + lag_y_3 + lag_y_4 + 
    s(epiWeek, bs = "cc", k = 12) +
    s(precip, k = 10) + s(spi3, k = 10) +
    region_id +
    offset(logpop),
  family = nb(),
  data = train,
  method = "REML"
)

summary(fit_A)
# plot(fit_A)


## 2) Recursive multi-step forecasting function

iter_forecast_gam_nb <- function(fit, history_df, future_df, lags = 1:4,
                                 nsim = 2000, seed = 1,
                                 recursion = c("mean", "path")) {
  recursion <- match.arg(recursion)
  set.seed(seed)

  history_df <- history_df %>% arrange(region_id, date)
  future_df  <- future_df  %>% arrange(region_id, date)

  out <- future_df %>%
    mutate(pred_mean = NA_real_, pred_lo = NA_real_, pred_hi = NA_real_)

  regions <- unique(out$region_id)

  # Helper: run one region with either deterministic recursion ("mean")
  # or proper path simulation ("path") for intervals.
  forecast_region <- function(hist_r, fut_r) {
    # ensure we have enough history for lags
    hist_r <- hist_r %>% arrange(date)

    if (recursion == "mean") {
      y_buf <- hist_r$y

      for (i in seq_len(nrow(fut_r))) {
        row_i <- fut_r[i, ]

        for (L in lags) row_i[[paste0("lag_y_", L)]] <- y_buf[length(y_buf) - (L - 1)]

        mu <- as.numeric(predict(fit, newdata = row_i, type = "response"))
        # simulate predictive distribution for intervals
        ysim <- rnbinom(nsim, mu = mu, size = fit$family$getTheta(TRUE))

        fut_r$pred_mean[i] <- mean(ysim)
        fut_r$pred_lo[i]   <- unname(quantile(ysim, 0.025, na.rm = TRUE))
        fut_r$pred_hi[i]   <- unname(quantile(ysim, 0.975, na.rm = TRUE))

        # deterministic recursion uses mean
        y_buf <- c(y_buf, fut_r$pred_mean[i])
      }
      fut_r

    } else {
      # PATH simulation: simulate many trajectories through time
      # then summarize across trajectories at each horizon.
      theta <- fit$family$getTheta(TRUE)

      # initial lags from observed history
      y0 <- hist_r$y
      nsteps <- nrow(fut_r)

      # store nsim x nsteps matrix of simulated future values
      simmat <- matrix(NA_real_, nrow = nsim, ncol = nsteps)

      for (s in 1:nsim) {
        if (s %% 10 == 1) { cat(".") }
        y_buf <- y0
        for (i in 1:nsteps) {
          row_i <- fut_r[i, ]
          for (L in lags) row_i[[paste0("lag_y_", L)]] <- y_buf[length(y_buf) - (L - 1)]
          mu <- as.numeric(predict(fit, newdata = row_i, type = "response"))
          y_new <- rnbinom(1, mu = mu, size = theta)
          simmat[s, i] <- y_new
          y_buf <- c(y_buf, y_new)
        }
      }

      fut_r$pred_mean <- colMeans(simmat)
      fut_r$pred_lo   <- apply(simmat, 2, quantile, probs = 0.025, na.rm = TRUE)
      fut_r$pred_hi   <- apply(simmat, 2, quantile, probs = 0.975, na.rm = TRUE)
      fut_r
    }
  }

  # Run region by region
  res <- map_dfr(regions, function(r) {
    hist_r <- history_df %>% filter(region_id == r)
    fut_r  <- out %>% filter(region_id == r) %>% arrange(date)
    forecast_region(hist_r, fut_r)
  })

  res %>% arrange(region_id, date)
}


## Example usage

# forecast_origin <- max(train$date)
# H <- 12

history_df = df

# future_df <- df %>%  # replace with your real future covariates frame
#   filter(date > forecast_origin) %>%
#   arrange(region_id, date) %>%
#   group_by(region_id) %>%
#   slice_head(n = H) %>%
#   ungroup() %>%
#   select(date, epiWeek, year, precip, spi3, logpop, region_id)

# undebug(iter_forecast_gam_nb)
fc_A <- iter_forecast_gam_nb(
  fit = fit_A,
  history_df = history_df,
  future_df = future_df,
  lags = 1:4,
  nsim = 200,
  seed = 42,
  recursion = "path"   # better uncertainty
)

fc_A

library(ggplot2)

ggplot() +
    geom_line(data = history_df[date > "2022-01-01"], aes(date, y), colour = "red") +
    geom_ribbon(data = fc_A, aes(date, ymin = pred_lo, ymax = pmin(100, pred_hi)), fill = "green", alpha = 0.4) + 
    geom_line(data = fc_A, aes(date, pmin(100, pred_mean)), colour = "green") +
    geom_line(data = fc_A, aes(date, y), colour = "blue") +
    facet_wrap(~region_id, scales = "free_y")





# Option B: Model for each lag

## 1) Horizon-specific training sets
make_direct_training <- function(df, h, lags = 1:4) {
  dat <- df %>%
    arrange(region_id, date) %>%
    group_by(region_id) %>%
    mutate(
      y_target    = lead(y, h),
      date_tgt    = lead(date, h),
      epiWeek_tgt = lead(epiWeek, h),
      year_tgt    = lead(year, h),
      precip_tgt  = lead(precip, h),
      spi3_tgt    = lead(spi3, h),
      logpop_tgt  = lead(logpop, h)
    ) %>%
    ungroup()

  dat <- make_lags(dat, lags)

  dat %>%
    filter(
      !is.na(y_target),
      !is.na(epiWeek_tgt),
      !is.na(precip_tgt),
      !is.na(spi3_tgt),
      !is.na(date_tgt),
      !is.na(logpop_tgt),
      if_all(all_of(paste0("lag_y_", lags)), ~ !is.na(.x))
    )
}

## 2) Fit 12 direct NB-GAMS
H <- 20
lags <- 1:4

fits_B <- map(setNames(1:H, paste0("h", 1:H)), function(h) {
  train_h <- make_direct_training(df, h = h, lags = lags)

  gam(
    y_target ~
      lag_y_1 + lag_y_2 + lag_y_3 + lag_y_4 +
      s(epiWeek_tgt, bs = "cc", k = 16) +
      s(precip_tgt, k = 10) + 
      s(spi3_tgt, k = 10) +
      s(as.integer(date_tgt), k = 10, by = region_id) +
      offset(logpop_tgt),
    family = nb(),
    data = train_h,
    method = "REML"
  )
})

summary(fits_B[["h12"]])

## 3) Forecast with the 12 models

direct_forecast_gam_nb <- function(fits, history_df, future_cov, forecast_origin,
                                   lags = 1:4, nsim = 2000, seed = 1) {
  set.seed(seed)

  hist <- history_df %>%
    filter(date <= forecast_origin) %>%
    arrange(region_id, date) %>%
    make_lags(lags)

  # one row per region at the forecast origin
  origin_rows <- hist %>%
    group_by(region_id) %>%
    filter(date == max(date)) %>%
    ungroup() %>%
    select(region_id, date, all_of(paste0("lag_y_", lags)))

  # future_cov: must contain covariates for each target date
  # columns: region_id, date, epiWeek, precip, spi3, logpop (year optional)
  future_cov <- future_cov %>%
    arrange(region_id, date) %>%
    filter(date > forecast_origin)

  out <- map_dfr(seq_along(fits), function(h) {
    fit_h <- fits[[h]]

    targets <- origin_rows %>%
      transmute(
        region_id,
        origin_date = date,
        horizon = h,
        target_date = date + weeks(h),
        lag_y_1, lag_y_2, lag_y_3, lag_y_4
      )

    newdata <- targets %>%
      left_join(
        future_cov %>% select(region_id, date, epiWeek, precip, spi3, logpop),
        by = c("region_id" = "region_id", "target_date" = "date")
      ) %>%
      rename(
        epiWeek_tgt = epiWeek,
        precip_tgt  = precip,
        spi3_tgt    = spi3,
        logpop_tgt  = logpop
      ) %>%
      filter(!is.na(epiWeek_tgt), !is.na(precip_tgt), !is.na(spi3_tgt), !is.na(logpop_tgt))
    
    newdata$date_tgt = newdata$target_date

    mu <- as.numeric(predict(fit_h, newdata = newdata, type = "response"))
    theta <- fit_h$family$getTheta(TRUE)

    # nsim draws per row
    ysim <- replicate(length(mu), rnbinom(nsim, mu = mu, size = theta))
    pred_mean <- colMeans(ysim)
    pred_lo   <- apply(ysim, 2, quantile, probs = 0.025)
    pred_hi   <- apply(ysim, 2, quantile, probs = 0.975)

    tibble(
      region_id   = newdata$region_id,
      origin_date = newdata$origin_date,
      horizon     = h,
      target_date = newdata$target_date,
      pred_mean   = pred_mean,
      pred_lo     = pred_lo,
      pred_hi     = pred_hi
    )
  })

  out %>% arrange(region_id, origin_date, horizon)
}

## 4) Example usage

forecast_origin <- max(df$date)

# supply future covariates for the next 12 weeks per region
fc_B <- direct_forecast_gam_nb(
  fits = fits_B,
  history_df = history_df,
  future_cov = future_df,
  forecast_origin = forecast_origin,
  lags = 1:4,
  nsim = 2000,
  seed = 42
)

fc_B


ggplot() +
    geom_line(data = df, aes(date, y), colour = "red") +
    geom_ribbon(data = fc_B, aes(target_date, ymin = pred_lo, ymax = pred_hi), fill = "black", alpha = 0.4) + 
    geom_line(data = fc_B, aes(target_date, pred_mean), colour = "black") +
    geom_line(data = future_df, aes(date, y), colour = "blue") +
    facet_wrap(~region_id, scales = "free_y")
