# Snippet for extracting trajectories from gam
library(mgcv)
library(mvtnorm)

dt = data.table(x = 0:100)

dt[, y := rnbinom(101, mu = x - 0.005*x^2, size = 5)]
plot(dt)

# Fit model
m <- gam(y ~ s(x), family = nb(), data = dt)

# Get predictions with standard errors
pred <- predict(m, newdata = dt, se.fit = TRUE, type = "link")
plot(dt$x, log(dt$y), type = "p")
lines(0:100, pred$fit, type = "l")
lines(0:100, pred$fit + 1.96 * pred$se, type = "l")
lines(0:100, pred$fit - 1.96 * pred$se, type = "l")




# Get posterior samples of coefficients
coef_samples <- rmvnorm(1000, mean = coef(m), sigma = vcov(m))

# Get design matrix for new data
Xp <- predict(m, newdata = dt, type = "lpmatrix")

# Compute linear predictor samples
eta_samples <- Xp %*% t(coef_samples)

# Transform to response scale
mu_samples <- exp(eta_samples)  # if using log link

# Now each column is a sample from posterior of μ
# Each defines a NB distribution with parameters (μ_i, θ)
theta <- m$family$getTheta(TRUE)  # extract estimated θ








# Get posterior samples of coefficients
coef_samples <- rmvnorm(1000, mean = coef(m), sigma = vcov(m))

# Get design matrix for new data
Xp <- predict(m, type = "lpmatrix")

# Compute linear predictor samples
eta_samples <- Xp %*% t(coef_samples)

# Transform to response scale
mu_samples <- exp(eta_samples)  # if using log link

# Now each column is a sample from posterior of μ
# Each defines a NB distribution with parameters (μ_i, θ)
theta <- m$family$getTheta(TRUE)  # extract estimated θ
