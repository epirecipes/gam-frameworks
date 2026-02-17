library(mgcv)
library(numDeriv)

# Generate synthetic data
set.seed(123)
n <- 200
x <- sort(runif(n, 0, 10))

# True model: smooth trend + local bumps
eta_true <- sin(x) + ifelse(x > 4 & x < 5, 2, 0) + ifelse(x > 7 & x < 8, -1.5, 0)
mu_true <- exp(eta_true)
y <- rpois(n, lambda = mu_true)

# Plot data
plot(x, y, main = "Synthetic Data", pch = 16, col = rgb(0,0,0,0.3))
lines(x, mu_true, col = "green", lwd = 2)

# ============================================
# JOINT OPTIMIZATION
# ============================================

# Represent GAM with basis functions (simplified)
# Use B-splines with few knots (smooth component)
library(splines)
B_gam <- bs(x, df = 6)  # 6 basis functions for smooth part

# Represent XGBoost-like component with step functions (local component)
# Use indicator basis for "bumps"
n_knots <- 20
knots <- seq(min(x), max(x), length.out = n_knots)
B_xgb <- outer(x, knots[-n_knots], function(x, k) as.numeric(x >= k & x < knots[match(k, knots) + 1]))
colnames(B_xgb) <- paste0("step_", 1:(n_knots-1))

# Joint parameter vector: [GAM coefficients, XGBoost coefficients]
n_gam <- ncol(B_gam)
n_xgb <- ncol(B_xgb)

# Negative log-likelihood + penalties
joint_objective <- function(params, lambda_gam = 0.1, lambda_xgb = 1.0) {
  
  # Split parameters
  beta_gam <- params[1:n_gam]
  beta_xgb <- params[(n_gam + 1):(n_gam + n_xgb)]
  
  # Predictions
  eta_gam <- B_gam %*% beta_gam
  eta_xgb <- B_xgb %*% beta_xgb
  eta_total <- as.vector(eta_gam + eta_xgb)
  mu <- exp(eta_total)
  
  # Negative log-likelihood (Poisson)
  nll <- -sum(dpois(y, lambda = mu, log = TRUE))
  
  # GAM penalty: roughness (second differences as proxy for smoothness)
  # Penalize wigglyness in GAM coefficients
  penalty_gam <- sum(diff(beta_gam, differences = 2)^2)
  
  # XGBoost penalty: L1 on coefficients (sparsity)
  penalty_xgb <- sum(abs(beta_xgb))
  
  # Total objective
  objective <- nll + lambda_gam * penalty_gam + lambda_xgb * penalty_xgb
  
  return(objective)
}

# Gradient (numerical for simplicity)
joint_gradient <- function(params, lambda_gam = 0.1, lambda_xgb = 1.0) {
  grad(joint_objective, params, lambda_gam = lambda_gam, lambda_xgb = lambda_xgb)
}

# Initialize parameters
init_params <- c(rep(0, n_gam), rep(0, n_xgb))

# Optimize jointly
cat("Running joint optimization...\n")
result_joint <- optim(
  par = init_params,
  fn = joint_objective,
  gr = joint_gradient,
  method = "L-BFGS-B",
  lambda_gam = 0.1,   # Penalty on GAM roughness
  lambda_xgb = 1.0,   # Penalty on XGBoost complexity
  control = list(trace = 1, maxit = 200)
)

# Extract fitted parameters
beta_gam_joint <- result_joint$par[1:n_gam]
beta_xgb_joint <- result_joint$par[(n_gam + 1):(n_gam + n_xgb)]

# Predictions
eta_gam_joint <- B_gam %*% beta_gam_joint
eta_xgb_joint <- B_xgb %*% beta_xgb_joint
eta_total_joint <- as.vector(eta_gam_joint + eta_xgb_joint)
mu_joint <- exp(eta_total_joint)

# ============================================
# COMPARE WITH BACKFITTING
# ============================================

backfitting <- function(max_iter = 20, lambda_gam = 0.1, lambda_xgb = 1.0) {
  
  # Initialize
  beta_gam <- rep(0, n_gam)
  beta_xgb <- rep(0, n_xgb)
  
  for (iter in 1:max_iter) {
    
    # Step 1: Fix XGBoost, optimize GAM
    eta_xgb <- B_xgb %*% beta_xgb
    
    obj_gam <- function(beta) {
      eta_gam <- B_gam %*% beta
      eta_total <- as.vector(eta_gam + eta_xgb)
      mu <- exp(eta_total)
      nll <- -sum(dpois(y, lambda = mu, log = TRUE))
      penalty <- lambda_gam * sum(diff(beta, differences = 2)^2)
      nll + penalty
    }
    
    result_gam <- optim(beta_gam, obj_gam, method = "L-BFGS-B", 
                       control = list(maxit = 50))
    beta_gam <- result_gam$par
    
    # Step 2: Fix GAM, optimize XGBoost
    eta_gam <- B_gam %*% beta_gam
    
    obj_xgb <- function(beta) {
      eta_xgb <- B_xgb %*% beta
      eta_total <- as.vector(eta_gam + eta_xgb)
      mu <- exp(eta_total)
      nll <- -sum(dpois(y, lambda = mu, log = TRUE))
      penalty <- lambda_xgb * sum(abs(beta))
      nll + penalty
    }
    
    result_xgb <- optim(beta_xgb, obj_xgb, method = "L-BFGS-B",
                       control = list(maxit = 50))
    beta_xgb <- result_xgb$par
    
    # Check convergence
    eta_total <- as.vector(B_gam %*% beta_gam + B_xgb %*% beta_xgb)
    mu <- exp(eta_total)
    obj_val <- -sum(dpois(y, lambda = mu, log = TRUE)) + 
               lambda_gam * sum(diff(beta_gam, differences = 2)^2) +
               lambda_xgb * sum(abs(beta_xgb))
    
    cat("Backfitting iteration", iter, "- Objective:", obj_val, "\n")
  }
  
  list(beta_gam = beta_gam, beta_xgb = beta_xgb)
}

cat("\nRunning backfitting...\n")
result_backfit <- backfitting(max_iter = 20, lambda_gam = 0.1, lambda_xgb = 1.0)

# Predictions from backfitting
eta_gam_backfit <- B_gam %*% result_backfit$beta_gam
eta_xgb_backfit <- B_xgb %*% result_backfit$beta_xgb
eta_total_backfit <- as.vector(eta_gam_backfit + eta_xgb_backfit)
mu_backfit <- exp(eta_total_backfit)

# ============================================
# VISUALIZE RESULTS
# ============================================

par(mfrow = c(2, 2))

# 1. Joint optimization: total fit
plot(x, y, main = "Joint Optimization: Total Fit", pch = 16, col = rgb(0,0,0,0.3))
lines(x, mu_true, col = "green", lwd = 2, lty = 2)
lines(x, mu_joint, col = "blue", lwd = 2)
legend("topleft", c("True", "Fitted"), col = c("green", "blue"), lwd = 2, lty = c(2, 1))

# 2. Joint optimization: decomposition
plot(x, eta_total_joint, type = "l", lwd = 2, main = "Joint: Decomposition", 
     ylab = "eta", ylim = range(c(eta_gam_joint, eta_xgb_joint, eta_total_joint)))
lines(x, eta_gam_joint, col = "red", lwd = 2)
lines(x, eta_xgb_joint, col = "purple", lwd = 2)
legend("topleft", c("Total", "GAM", "XGB"), col = c("black", "red", "purple"), lwd = 2)

# 3. Backfitting: total fit
plot(x, y, main = "Backfitting: Total Fit", pch = 16, col = rgb(0,0,0,0.3))
lines(x, mu_true, col = "green", lwd = 2, lty = 2)
lines(x, mu_backfit, col = "orange", lwd = 2)
legend("topleft", c("True", "Fitted"), col = c("green", "orange"), lwd = 2, lty = c(2, 1))

# 4. Backfitting: decomposition
plot(x, eta_total_backfit, type = "l", lwd = 2, main = "Backfitting: Decomposition",
     ylab = "eta", ylim = range(c(eta_gam_backfit, eta_xgb_backfit, eta_total_backfit)))
lines(x, eta_gam_backfit, col = "red", lwd = 2)
lines(x, eta_xgb_backfit, col = "purple", lwd = 2)
legend("topleft", c("Total", "GAM", "XGB"), col = c("black", "red", "purple"), lwd = 2)

# ============================================
# COMPARE OBJECTIVES
# ============================================

cat("\n=== COMPARISON ===\n")
cat("Joint optimization final objective:", result_joint$value, "\n")

eta_backfit_total <- as.vector(B_gam %*% result_backfit$beta_gam + 
                               B_xgb %*% result_backfit$beta_xgb)
mu_backfit_check <- exp(eta_backfit_total)
obj_backfit <- -sum(dpois(y, lambda = mu_backfit_check, log = TRUE)) +
               0.1 * sum(diff(result_backfit$beta_gam, differences = 2)^2) +
               1.0 * sum(abs(result_backfit$beta_xgb))
cat("Backfitting final objective:", obj_backfit, "\n")

if (result_joint$value < obj_backfit) {
  cat("Joint optimization achieved BETTER objective (lower is better)\n")
} else {
  cat("Backfitting achieved similar or better objective\n")
}
