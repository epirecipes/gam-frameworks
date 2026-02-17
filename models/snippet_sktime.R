library(reticulate)
library(ggplot2)
library(dplyr)


# Setup (run once)
if (0) {
    # python 3.12
    # reticulate::install_python(version = '3.12')
  
    # create r-sentinel environment
    virtualenv_create("r-sentinel", python = "/opt/homebrew/bin/python3.12")

    # install packages
    virtualenv_install("r-sentinel", c("sktime", "numpy", "pandas", "statsmodels"))
    virtualenv_install("r-sentinel", c("pmdarima", "neuralforecast", "u8darts", "xgboost"))
}

# import packages
use_virtualenv("r-sentinel")
np = import("numpy", convert = FALSE)
pd = import("pandas", convert = FALSE)
sk = import("sktime", convert = FALSE)

# Create synthetic time series data with seasonality and negative binomial distribution
set.seed(123)
n_points = 730
time_index = seq(as.Date("2020-01-01"), by = "day", length.out = n_points)

# Create seasonal pattern (biweekly seasonality)
seasonal_component = 15 + 10 * sin(2 * pi * seq_len(n_points) / 14)

# Add trend
trend_component = seq(0, 5, length.out = n_points)

# Create covariates
temperature = 20 + 10 * sin(2 * pi * seq_len(n_points) / 365) + rnorm(n_points, 0, 2)
promotion = rbinom(n_points, 1, 0.15)  # 15% of days have promotions

# Create count data using poisson distribution
# Mean is determined by trend + seasonality + covariate effects
lambda = exp(log(seasonal_component + trend_component) + 
              0.3 * scale(temperature)[,1] + 
              0.5 * promotion)
counts = rpois(n_points, lambda)

# Create data frame
df = data.frame(
    date = time_index,
    counts = counts,
    temperature = temperature,
    promotion = promotion,
    season = rep(0:13, length.out = n_points),
    sea1 = sin(2 * pi * seq_len(n_points) / 14),
    sea2 = cos(2 * pi * seq_len(n_points) / 14)
)

ggplot(df) +
  geom_line(aes(x = date, y = temperature, colour = "temperature")) +
  geom_line(aes(x = date, y = counts, colour = "count")) +
  geom_step(aes(x = date, y = promotion, colour = "promotion")) +
  geom_line(aes(x = date, y = sea1)) +
  geom_line(aes(x = date, y = sea2)) 


# Split into train and test
train_size = 365
train_df = df[1:train_size, ]
test_df = df[(train_size + 1):n_points, ]

cat("Training data:", nrow(train_df), "observations\n")
cat("Test data:", nrow(test_df), "observations\n")
cat("Count range:", range(df$count), "\n")

# Convert to pandas DataFrames for sktime
# Create proper time index
train_pd = pd$DataFrame(train_df)
train_pd$index = pd$to_datetime(train_df$date)
train_pd = train_pd$drop("date", axis=1L)

test_pd = pd$DataFrame(test_df)
test_pd$index = pd$to_datetime(test_pd$date)
test_pd = test_pd$drop("date", axis=1L)

# Prepare target (y) and features (X) for training
y_train = train_pd$counts
X_train = train_pd$drop("counts", axis = 1L)

y_test = test_pd$counts
X_test = test_pd$drop("counts", axis = 1L)

# ============================================================
# Define and fit all models
# ============================================================

# List to store all model configurations and results
models = list(
    naive = list(
        name = "Naive (Last Value)",
        forecaster = sk$forecasting$naive$NaiveForecaster(strategy = "last"),
        use_X = FALSE
    ),
    arimax = list(
        name = "AutoARIMAX",
        forecaster = sk$forecasting$arima$AutoARIMA(
            seasonal = TRUE,
            suppress_warnings = TRUE,
            error_action = "ignore"
        ),
        use_X = TRUE
    ),
    dartsxgb = list(
        name = "DartsXGB",
        forecaster = sk$forecasting$darts$DartsXGBModel(
            lags = 14L
        ),
        use_X = TRUE
    )
    # nftcn = list(
    #     name = "NFTCN",
    #     forecaster = sk$forecasting$neuralforecast$NeuralForecastTCN(
    #         input_size = 300L
    #     ),
    #     use_X = TRUE
    # ),
    # nfrnn = list(
    #     name = "NFRNN",
    #     forecaster = sk$forecasting$neuralforecast$NeuralForecastRNN(
    #         input_size = 300L
    #     ),
    #     use_X = TRUE
    # ),
    # nfdrnn = list(
    #     name = "NFDRNN",
    #     forecaster = sk$forecasting$neuralforecast$NeuralForecastDilatedRNN(
    #         input_size = 100L
    #     ),
    #     use_X = TRUE
    # ),
    # nfgru = list(
    #     name = "NFGRU",
    #     forecaster = sk$forecasting$neuralforecast$NeuralForecastGRU(
    #         input_size = 300L
    #     ),
    #     use_X = TRUE
    # ),
    # nflstm = list(
    #     name = "NFLSTM",
    #     forecaster = sk$forecasting$neuralforecast$NeuralForecastLSTM(
    #         input_size = 300L
    #     ),
    #     use_X = TRUE
    # )
)

models = mapply(function(m, i) { 
    m$color = hcl.colors(length(models), palette = "Zissou 1")[[i]]; 
    return (m) 
}, models, seq_along(models), SIMPLIFY = FALSE)


# Create forecast horizon (used by all models)
fh = np$arange(1L, as.integer(nrow(test_df) + 1))

# Fit all models and generate predictions
for (model_key in names(models)) {
    cat("\n", model_key, "\n")
    model = models[[model_key]]
    
    # Fit
    if (model$use_X) {
        model$forecaster$fit(fh = fh, y = y_train, X = X_train)
        model$prediction = model$forecaster$predict(X = X_test)
    } else {
        model$forecaster$fit(fh = fh, y = y_train)
        model$prediction = model$forecaster$predict()
    }
    
    # Store back
    models[[model_key]] = model
}

# ============================================================
# Convert predictions to data frame
# ============================================================
predictions_df = data.frame(
    date = test_df$date,
    actual = test_df$count
)

for (model_key in names(models)) {
    predictions_df[[model_key]] = as.numeric(models[[model_key]]$prediction$values)
}

# ============================================================
# Calculate scoring metrics
# ============================================================
mean_absolute_error = sk$performance_metrics$forecasting$mean_absolute_error
mean_squared_error = sk$performance_metrics$forecasting$mean_squared_error

metrics_list = list()
for (model_key in names(models)) {
    pred = models[[model_key]]$prediction
    metrics_list[[model_key]] = list(
        Model = models[[model_key]]$name,
        MAE = py_to_r(mean_absolute_error(y_test, pred)),
        MSE = py_to_r(mean_squared_error(y_test, pred))
    )
}

metrics_summary = do.call(rbind, lapply(metrics_list, as.data.frame))
metrics_summary$RMSE = sqrt(metrics_summary$MSE)

print("Model Performance Metrics:")
print(metrics_summary)

# ============================================================
# Visualization with ggplot2
# ============================================================

# Reshape data for plotting
model_cols = names(models)
plot_data = predictions_df %>%
    tidyr::pivot_longer(
        cols = c(actual, all_of(model_cols)),
        names_to = "series",
        values_to = "value"
    ) %>%
    mutate(
        series = factor(series, 
                       levels = c("actual", model_cols),
                       labels = c("Actual", sapply(models, function(m) m$name)))
    )

# Create color scale
model_colors = c("Actual" = "black", sapply(models, function(m) m$color))
names(model_colors) = c("Actual", sapply(models, function(m) m$name))

# Create the plot
ggplot(plot_data, aes(x = date, y = value, color = series)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = model_colors) +
    labs(
        title = "Time Series Forecasting Comparison",
        subtitle = paste0("Test period: ", nrow(test_df), " days | ",
                         "Best model: ", metrics_summary$Model[which.min(metrics_summary$MAE)]),
        x = "Date",
        y = "Count",
        color = "Series",
        linetype = "Series"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
    )

# Additional plot: Residuals
residuals_df = data.frame(date = test_df$date)
for (model_key in names(models)) {
    residuals_df[[model_key]] = test_df$count - predictions_df[[model_key]]
}

residuals_df = residuals_df %>%
    tidyr::pivot_longer(
        cols = all_of(model_cols),
        names_to = "model",
        values_to = "residual"
    ) %>%
    mutate(
        model = factor(model,
                      levels = model_cols,
                      labels = sapply(models, function(m) m$name))
    )

# Create color scale for residuals
residual_colors = sapply(models, function(m) m$color)
names(residual_colors) = sapply(models, function(m) m$name)

p_residuals = ggplot(residuals_df, aes(x = date, y = residual, color = model)) +
    geom_line(alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    facet_wrap(~model, ncol = 1) +
    scale_color_manual(values = residual_colors) +
    labs(
        title = "Forecast Residuals by Model",
        x = "Date",
        y = "Residual (Actual - Predicted)"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "none",
        strip.text = element_text(face = "bold")
    )

print(p_residuals)






py_run_string("
from sktime.datasets import load_airline
from sktime.forecasting.neuralforecast import NeuralForecastTCN
from sktime.forecasting.base import ForecastingHorizon

# data
y = load_airline()

# define 365-step horizon
fh = ForecastingHorizon(range(1, 366), is_relative=True)

# TCN forecaster
forecaster = NeuralForecastTCN(
    input_size=300,    # history window
    max_steps=500
)

# fit
forecaster.fit(fh = fh, y = y)

# predict full horizon directly
y_pred = forecaster.predict()"
)
