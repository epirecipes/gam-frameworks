# Compute discrete predictive coverage for mgcv::gam (Poisson / NB)
# - method = "parametric" (recommended) or "plug-in"
# - For NB, theta is taken from the fit; optionally sample log-theta if theta_log_sd is provided
# - Handles offsets via predict(); include any offsets in `newdata`

# Predictive coverage for mgcv::gam (Poisson / Negative Binomial)
gam_predictive_coverage <- function(
    fit,
    newdata = NULL,
    y = NULL,
    level = 0.95,
    family_override = c("auto","poisson","nb"),
    method = c("parametric","plug-in"),
    nsim = 5000,
    use_shortest = FALSE,
    # parametric options:
    parametric_engine = c("simulate","lpmatrix"),
    unconditional = TRUE,   # include smoothing-parameter uncertainty (simulate + vcov)
    theta_log_sd = NULL     # optional: sample NB theta (only used in lpmatrix mode)
) {
  stopifnot(inherits(fit, "gam"))
  family_override <- match.arg(family_override)
  method <- match.arg(method)
  parametric_engine <- match.arg(parametric_engine)
  
  # Resolve family
  fam_str <- fit$family$family
  fam <- if (family_override == "auto") {
    if (grepl("Negative Binomial", fam_str, ignore.case = TRUE)) "nb" else "poisson"
  } else family_override
  if (!fam %in% c("poisson","nb")) stop("Only Poisson and Negative Binomial are supported.")
  
  # Data and observed y
  if (is.null(newdata)) newdata <- fit$model
  if (is.null(y)) {
    resp <- as.character(formula(fit)[[2]])
    y <- if (!is.null(newdata[[resp]])) newdata[[resp]] else fit$y
    if (is.null(y)) stop("Please supply 'y' or include the response in 'newdata'.")
  }
  y <- as.integer(y)
  n <- NROW(newdata); if (length(y) != n) stop("Length(y) must equal nrow(newdata).")
  
  # (optional) shortest discrete predictive set
  .shortest_discrete_set <- function(pmf, level = 0.95) {
    ord <- order(pmf, decreasing = TRUE)
    S <- sort(as.integer(names(pmf))[seq_len(which(cumsum(pmf[ord]) >= level)[1])])
    c(min(S), max(S))
  }
  
  if (method == "plug-in") {
    mu_hat <- as.numeric(predict(fit, data = newdata, type = "response"))
    if (fam == "poisson") {
      L <- qpois((1 - level)/2, lambda = mu_hat)
      U <- qpois(1 - (1 - level)/2, lambda = mu_hat)
    } else {
      theta_hat <- tryCatch(fit$family$getTheta(TRUE), error = function(e) NULL)
      if (is.null(theta_hat)) {
        theta_hat <- as.numeric(sub(".*\\(([^)]+)\\).*", "\\1", fam_str))
        if (!is.finite(theta_hat)) stop("Could not extract NB theta from the fit.")
      }
      L <- qnbinom((1 - level)/2, size = theta_hat, mu = mu_hat)
      U <- qnbinom(1 - (1 - level)/2, size = theta_hat, mu = mu_hat)
    }
    covered <- (y >= L) & (y <= U)
    return(list(method = "plug-in", family = fam, level = level,
                coverage = mean(covered, na.rm = TRUE),
                covered = covered, L = as.integer(L), U = as.integer(U),
                width = as.integer(U - L)))
  }
  
  # ---------------- PARAMETRIC ----------------
  # Preferred: stats::simulate() dispatches to mgcv's simulate.gam (S3), handles offsets.
  ysim <- NULL
  if (parametric_engine == "simulate") {
    sims <- try(simulate(fit, nsim = nsim, data = newdata,
                         unconditional = unconditional, seed = NULL),
                silent = TRUE)
    if (!inherits(sims, "try-error")) {
      # sims: n x nsim (data.frame or matrix); convert to nsim x n
      ysim <- t(as.matrix(sims))
    }
    # if simulate() failed, fall through to lpmatrix fallback
  }
  
  if (is.null(ysim)) {
    # Fallback: manual parametric bootstrap (must add offset explicitly)
    X <- predict(fit, newdata = newdata, type = "lpmatrix")
    b_hat <- coef(fit)
    Vb <- try(vcov(fit, unconditional = unconditional), silent = TRUE)
    if (inherits(Vb, "try-error") || any(!is.finite(Vb))) Vb <- fit$Vp
    if (NCOL(X) != length(b_hat)) {
      stop("lpmatrix columns (", NCOL(X), ") != length(coef) (", length(b_hat),
           "). Check factor levels and spline bases in 'newdata'.")
    }
    # Rebuild offset for newdata: eta(with offset) - X %*% b_hat
    eta_hat_with_offset <- as.numeric(predict(fit, newdata = newdata, type = "link"))
    off_vec <- eta_hat_with_offset - as.numeric(X %*% b_hat)
    
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Please install 'MASS' for mvrnorm().")
    beta_sim <- MASS::mvrnorm(nsim, mu = b_hat, Sigma = Vb)  # nsim x p
    eta_sim  <- beta_sim %*% t(X)
    eta_sim  <- sweep(eta_sim, 2, off_vec, `+`)              # add offset
    mu_sim   <- exp(eta_sim)
    
    if (fam == "poisson") {
      ysim <- matrix(rpois(nsim * n, lambda = as.vector(mu_sim)),
                     nrow = nsim, ncol = n, byrow = TRUE)
    } else {
      theta_hat <- tryCatch(fit$family$getTheta(TRUE), error = function(e) NULL)
      if (is.null(theta_hat)) {
        theta_hat <- as.numeric(sub(".*\\(([^)]+)\\).*", "\\1", fam_str))
        if (!is.finite(theta_hat)) stop("Could not extract NB theta from the fit.")
      }
      if (!is.null(theta_log_sd) && theta_log_sd > 0) {
        theta_b <- exp(rnorm(nsim, mean = log(theta_hat), sd = theta_log_sd))
      } else {
        theta_b <- rep(theta_hat, nsim)
      }
      ysim <- matrix(NA_integer_, nsim, n)
      for (b in 1:nsim) ysim[b, ] <- rnbinom(n, size = theta_b[b], mu = mu_sim[b, ])
    }
  }
  
  # Discrete PIs + coverage
  alpha <- 1 - level
  if (!use_shortest) {
    L <- apply(ysim, 2, quantile, probs = alpha/2, type = 1, names = FALSE)
    U <- apply(ysim, 2, quantile, probs = 1 - alpha/2, type = 1, names = FALSE)
  } else {
    L <- U <- integer(n)
    for (j in seq_len(n)) {
      tab <- table(ysim[, j]); pmf <- as.numeric(tab) / nsim; names(pmf) <- names(tab)
      lims <- .shortest_discrete_set(pmf, level = level)
      L[j] <- lims[1]; U[j] <- lims[2]
    }
  }
  covered <- (y >= L) & (y <= U)
  
  list(
    method   = paste0("parametric/", if (!is.null(attr(ysim, "from_simulate"))) "simulate" else "lpmatrix",
                      if (unconditional) "+unconditional" else ""),
    family   = fam,
    level    = level,
    coverage = mean(covered, na.rm = TRUE),
    covered  = covered,
    L = as.integer(L), U = as.integer(U),
    width = as.integer(U - L),
    nsim = nsim
  )
}

