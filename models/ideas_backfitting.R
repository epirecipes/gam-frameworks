library(mgcv)
library(xgboost)

# Generate synthetic data
set.seed(123)
n <- 500
x <- sort(runif(n, 0, 10))

# True model: smooth trend + local complex patterns
eta_true <- 0.5 * sin(x) + 
            ifelse(x > 3 & x < 4, 1.5, 0) + 
            ifelse(x > 6 & x < 7, -1.0, 0) +
            0.3 * x
mu_true <- exp(eta_true)
y <- rpois(n, lambda = mu_true)

# Plot data
plot(x, y, main = "Synthetic Data", pch = 16, col = rgb(0,0,0,0.3))
lines(x, mu_true, col = "green", lwd = 2)
legend("topleft", "True mean", col = "green", lwd = 2)

# ============================================
# ROBUST BACKFITTING WITH GAM + XGBOOST
# ============================================

backfitting_gam_xgb <- function(x, y, max_iter = 10, 
                                 xgb_nrounds = 30, 
                                 xgb_eta = 0.3,
                                 xgb_max_depth = 3) {
  
  n <- length(y)
  data_df <- data.frame(x = x, y = y)
  
  # Initialize with more stable starting values
  eta_gam <- rep(log(mean(y) + 1), n)  # Avoid log(0)
  eta_xgb <- rep(0, n)
  
  # Storage for convergence tracking
  objectives <- numeric(max_iter)
  
  for (iter in 1:max_iter) {
    
    cat("\n=== Iteration", iter, "===\n")
    
    # ------------------------------------------------
    # Step 1: Fix XGBoost component, update GAM
    # ------------------------------------------------
    
    # Use XGBoost prediction as offset
    data_df$offset_xgb <- eta_xgb
    
    # Fit GAM with offset
    tryCatch({
      gam_model <- gam(
        y ~ s(x, k = 10) + offset(offset_xgb),
        family = poisson(link = "log"),
        data = data_df,
        method = "REML"
      )
      
      # Get GAM predictions (on eta scale, without offset)
      eta_gam <- predict(gam_model, type = "link") - eta_xgb
      
      cat("GAM fit complete. Deviance:", deviance(gam_model), "\n")
      
    }, error = function(e) {
      cat("GAM fitting failed:", e$message, "\n")
      cat("Keeping previous GAM fit\n")
    })
    
    # ------------------------------------------------
    # Step 2: Fix GAM component, update XGBoost
    # ------------------------------------------------
    
    # Clip eta_gam to prevent extreme values
    eta_gam_clipped <- pmax(pmin(eta_gam, 10), -10)
    
    # Custom XGBoost objective with numerical stability
    xgb_objective_with_offset <- function(preds, dtrain) {
      labels <- getinfo(dtrain, "label")
      offset <- getinfo(dtrain, "base_margin")
      
      # Clip predictions to prevent overflow
      preds_clipped <- pmax(pmin(preds, 10), -10)
      
      # Combined prediction
      eta_total <- preds_clipped + offset
      eta_total_clipped <- pmax(pmin(eta_total, 10), -10)
      
      mu <- exp(eta_total_clipped)
      
      # Add small epsilon for numerical stability
      epsilon <- 1e-10
      mu <- pmax(mu, epsilon)
      
      # Poisson gradients with clipping
      grad <- mu - labels
      hess <- mu
      
      # Ensure hessian is positive and bounded
      hess <- pmax(hess, epsilon)
      hess <- pmin(hess, 1e6)
      
      # Check for any NaN or Inf
      if (any(is.nan(grad)) || any(is.infinite(grad)) ||
          any(is.nan(hess)) || any(is.infinite(hess))) {
        cat("Warning: NaN or Inf detected in gradients\n")
        grad[is.nan(grad) | is.infinite(grad)] <- 0
        hess[is.nan(hess) | is.infinite(hess)] <- 1
      }
      
      return(list(grad = grad, hess = hess))
    }
    
    # Prepare XGBoost data with GAM as offset
    dtrain <- xgb.DMatrix(data = as.matrix(x), label = y)
    setinfo(dtrain, "base_margin", eta_gam_clipped)
    
    # Train XGBoost with error handling
    tryCatch({
      xgb_model <- xgb.train(
        data = dtrain,
        objective = xgb_objective_with_offset,
        nrounds = xgb_nrounds,
        max_depth = xgb_max_depth,
        eta = xgb_eta,
        verbose = 0,
        # Add regularization for stability
        lambda = 1.0,         # L2 regularization
        alpha = 0.5,          # L1 regularization
        min_child_weight = 5  # Prevent overfitting to small groups
      )
      
      # Get XGBoost predictions with clipping
      eta_xgb_raw <- predict(xgb_model, dtrain)
      eta_xgb <- pmax(pmin(eta_xgb_raw, 10), -10)
      
      cat("XGBoost fit complete. Trees:", xgb_nrounds, "\n")
      
    }, error = function(e) {
      cat("XGBoost fitting failed:", e$message, "\n")
      cat("Keeping previous XGBoost fit\n")
    })
    
    # ------------------------------------------------
    # Check convergence
    # ------------------------------------------------
    
    # Combined prediction with clipping
    eta_total <- eta_gam + eta_xgb
    eta_total_clipped <- pmax(pmin(eta_total, 10), -10)
    mu_combined <- exp(eta_total_clipped)
    
    # Add epsilon to prevent log(0)
    mu_combined <- pmax(mu_combined, 1e-10)
    
    # Compute Poisson log-likelihood
    loglik <- sum(dpois(y, lambda = mu_combined, log = TRUE))
    
    if (is.nan(loglik) || is.infinite(loglik)) {
      cat("Warning: Invalid log-likelihood, stopping\n")
      break
    }
    
    objectives[iter] <- -loglik  # Store negative log-likelihood
    
    cat("Combined log-likelihood:", loglik, "\n")
    
    # Check for convergence
    if (iter > 1) {
      change <- abs(objectives[iter] - objectives[iter - 1])
      cat("Change in objective:", change, "\n")
      
      if (change < 1e-3) {
        cat("Converged!\n")
        objectives <- objectives[1:iter]
        break
      }
      
      # Check for divergence
      if (objectives[iter] > objectives[iter - 1] * 1.1) {
        cat("Warning: Objective increasing significantly, may be diverging\n")
      }
    }
  }
  
  list(
    gam_model = gam_model,
    xgb_model = xgb_model,
    eta_gam = eta_gam,
    eta_xgb = eta_xgb,
    eta_total = eta_total,
    mu_combined = mu_combined,
    objectives = objectives[1:iter]
  )
}

# Run backfitting with error handling
cat("Starting backfitting algorithm...\n")

result <- tryCatch({
  backfitting_gam_xgb(
    x = x, 
    y = y, 
    max_iter = 15,
    xgb_nrounds = 20,  # Reduced from 30
    xgb_eta = 0.3,
    xgb_max_depth = 3
  )
}, error = function(e) {
  cat("Backfitting failed completely:", e$message, "\n")
  NULL
})

if (is.null(result)) {
  stop("Could not complete backfitting")
}

# ============================================
# VISUALIZE RESULTS
# ============================================

par(mfrow = c(2, 2))

# 1. Total fit
plot(x, y, main = "Backfitting: Total Fit", pch = 16, col = rgb(0,0,0,0.3))
lines(x, mu_true, col = "green", lwd = 2, lty = 2)
lines(x, result$mu_combined, col = "blue", lwd = 2)
legend("topleft", c("True", "Fitted"), col = c("green", "blue"), 
       lwd = 2, lty = c(2, 1))

# 2. Decomposition on eta scale
plot(x, result$eta_total, type = "l", lwd = 2, 
     main = "Decomposition (eta scale)",
     ylab = "eta = log(mu)", 
     ylim = range(c(result$eta_gam, result$eta_xgb, result$eta_total)))
lines(x, result$eta_gam, col = "red", lwd = 2)
lines(x, result$eta_xgb, col = "purple", lwd = 2)
abline(h = 0, lty = 3, col = "gray")
legend("topleft", c("Total", "GAM component", "XGBoost component"), 
       col = c("black", "red", "purple"), lwd = 2)

# 3. GAM component only
plot(x, exp(result$eta_gam), type = "l", col = "red", lwd = 2,
     main = "GAM Component (smooth trend)",
     ylab = "mu_GAM = exp(eta_GAM)")

# 4. XGBoost component only
plot(x, exp(result$eta_xgb), type = "l", col = "purple", lwd = 2,
     main = "XGBoost Component (local patterns)",
     ylab = "Multiplicative factor = exp(eta_XGB)")
abline(h = 1, lty = 3, col = "gray")

# 5. Convergence plot
par(mfrow = c(1, 1))
plot(1:length(result$objectives), result$objectives, type = "b",
     main = "Convergence: Negative Log-Likelihood",
     xlab = "Iteration", ylab = "Negative Log-Likelihood",
     pch = 16, col = "blue", lwd = 2)
grid()

cat("\n=== FINAL STATISTICS ===\n")
cat("Final log-likelihood:", -tail(result$objectives, 1), "\n")
cat("GAM eta range:", range(result$eta_gam), "\n")
cat("XGBoost eta range:", range(result$eta_xgb), "\n")
cat("Combined mu range:", range(result$mu_combined), "\n")
