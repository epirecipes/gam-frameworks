# ------------------------------------------------------------------------------
# INLA posterior-predictive coverage utilities for count models
#
# What this file does
# -------------------
# - Provides helpers to:
#   * Build shortest discrete predictive sets (.shortest_discrete_set)
#   * Map requested Predictor rows to positions in the INLA latent vector,
#     robustly across INLA versions and configs (.map_predictor_positions)
#   * Extract linear predictor draws (η) at those positions from posterior samples
#     (.get_eta_samples_by_pos)
#   * Extract zero-inflation probabilities and NB size (overdispersion) from
#     hyperparameter draws (.get_zip_pi_draws, .get_nb_size_draws)
# - Main user function:
#   * inla_predictive_coverage(): computes discrete prediction intervals and
#     coverage for Poisson, Negative Binomial, Binomial, and zero-inflated variants
#
# Key ideas
# ---------
# - We work on the *predictive* scale by simulating integer outcomes given
#   posterior uncertainty over parameters (via inla.posterior.sample()).
# - For Poisson/NB families we invert the log link: μ = exp(η) (and optionally
#   multiply μ by exposure E only when a Poisson model was fit with E=).
# - For binomial we invert the logit: p = plogis(η) and need the trials (Ntrials).
# - For zero-inflated models we simulate structural zeros with probability π.
#
# Inputs of inla_predictive_coverage()
# ------------------------------------
# - res      : inla result fitted with control.predictor=list(compute=TRUE)
# - y        : observed integer outcomes of the rows you evaluate (same length as 'indices')
# - family   : "poisson","nbinomial","binomial","zeroinflatedpoisson0/1",
#              "zeroinflatednbinomial0/1" (autodetected if NULL)
# - indices  : which rows in summary.linear.predictor to evaluate (default: all)
# - level    : nominal coverage level (default 0.95)
# - nsim     : number of posterior predictive draws
# - use_shortest : if TRUE, computes shortest discrete predictive set; else equal-tailed quantiles
# - E        : exposure vector (ONLY for Poisson models fitted with E=)
# - Ntrials  : binomial trials per row (required if y not in {0,1})
#
# Outputs
# -------
# A list with coverage, indicator per row (covered), lower/upper bounds, etc.
#
# Notes
# -----
# - The code supports two ways of mapping predictor rows to latent positions:
#   via res$misc$configs$contents (when config=TRUE), else via latent names.
# - Quantiles use type=1 (discrete); this preserves integer-valued intervals.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages(library(INLA))

# ------------------------------------------------------------------------------
# Build the shortest (highest-mass) discrete predictive set at 'level'
# Input: pmf named vector with integer names; Output: [L, U] bounds
# ------------------------------------------------------------------------------
.shortest_discrete_set <- function(pmf, level = 0.95) {
  ord <- order(pmf, decreasing = TRUE)
  k <- which(cumsum(pmf[ord]) >= level)[1]
  S <- sort(as.integer(names(pmf))[ord[seq_len(k)]])
  c(min(S), max(S))
}

# ------------------------------------------------------------------------------
# Given a named pmf vector (names are integers), return the shortest
# contiguous [L,U]
# ------------------------------------------------------------------------------
.shortest_contiguous_interval <- function(pmf, level = 0.95) {
  k <- as.integer(names(pmf))
  o <- order(k); k <- k[o]; p <- pmf[o]
  
  # fill any gaps in support with 0-probability
  full_k <- seq.int(min(k), max(k))
  p_full <- rep(0, length(full_k)); p_full[match(k, full_k)] <- p
  
  cs <- c(0, cumsum(p_full))  # cumulative mass on the filled grid
  n  <- length(full_k)
  best_L <- full_k[1]; best_U <- full_k[n]; best_len <- n + 1
  
  i <- 1
  for (j in 1:n) {
    # advance i while mass in [i, j] is >= level (keep it as tight as possible)
    while (i <= j && (cs[j + 1] - cs[i]) >= level) {
      len <- j - i + 1
      if (len < best_len) { best_len <- len; best_L <- full_k[i]; best_U <- full_k[j] }
      i <- i + 1
    }
  }
  c(best_L, best_U)
}


# ------------------------------------------------------------------------------
# Map global predictor row indices (from summary.linear.predictor) to positions
# in the INLA latent vector. Prefer contents (requires control.compute$config=TRUE),
# else fall back to parsing latent names in a posterior sample.
# ------------------------------------------------------------------------------
.map_predictor_positions <- function(res, indices, samples_for_fallback) {
  contents <- tryCatch(res$misc$configs$contents, error = function(e) NULL)
  if (!is.null(contents) && all(c("tag","start","n") %in% names(contents))) {
    pred_rows <- which(grepl("(?i)^predictor", contents$tag, perl = TRUE))
    if (length(pred_rows) > 0L) {
      starts <- contents$start[pred_rows]
      lens   <- contents$n[pred_rows]
      cumlen <- cumsum(lens)
      block_of <- function(i) which(i <= cumlen)[1]
      pos <- integer(length(indices))
      for (k in seq_along(indices)) {
        i <- indices[k]
        b <- block_of(i); if (is.na(b)) stop("Index ", i, " is out of range.")
        prev <- if (b == 1L) 0L else cumlen[b - 1L]
        local_i <- i - prev
        pos[k] <- starts[b] + local_i - 1L
      }
      return(pos)
    }
  }
  # Fallback: derive positions by matching "Predictor" entries in latent names
  lat <- samples_for_fallback[[1]]$latent
  nm <- names(lat); if (is.null(nm)) nm <- rownames(lat)
  if (is.null(nm)) stop("Latent has no names; cannot locate predictors.")
  pred_idx <- grep("(?i)\\bpredictor\\b", nm, perl = TRUE)
  if (!length(pred_idx)) stop("No predictor entries in latent names.")
  suff <- sub("^.*?([0-9]+)\\s*$", "\\1", nm[pred_idx], perl = TRUE)
  if (all(grepl("^[0-9]+$", suff))) {
    pred_positions_all <- pred_idx[order(as.integer(suff))]
  } else {
    warning("Could not parse predictor indices; assuming current order.")
    pred_positions_all <- pred_idx
  }
  if (max(indices) > length(pred_positions_all)) {
    stop("Requested predictor index exceeds available predictor entries.")
  }
  pred_positions_all[indices]
}

# ------------------------------------------------------------------------------
# Extract matrix of η (linear predictor) samples at given latent positions.
# Returns nsim x n matrix: each row is one posterior draw, each column an index.
# ------------------------------------------------------------------------------
.get_eta_samples_by_pos <- function(samples, latent_positions) {
  nsim <- length(samples); n <- length(latent_positions)
  eta <- matrix(NA_real_, nsim, n)
  for (b in seq_len(nsim)) {
    lat <- samples[[b]]$latent
    if (is.matrix(lat) && ncol(lat) == 1L) lat <- as.numeric(lat[,1L])
    eta[b, ] <- lat[latent_positions]
  }
  if (any(!is.finite(eta))) stop("Failed to extract Predictor entries from samples.")
  eta
}

# ------------------------------------------------------------------------------
# Extract zero-inflation probabilities π from hyperparameter draws.
# Heuristic: if value is not in [0,1], assume logit scale and apply plogis().
# ------------------------------------------------------------------------------
.get_zip_pi_draws <- function(samples) {
  v <- vapply(samples, function(s) {
    hp <- s$hyperpar; if (is.null(hp) || !length(hp)) return(NA_real_)
    nm <- names(hp); idx <- which(grepl("zero|infl", nm, ignore.case = TRUE))
    if (!length(idx)) return(NA_real_)
    val <- unname(hp[idx[1]])
    if (is.na(val)) NA_real_ else if (val < 0 || val > 1) plogis(val) else val
  }, numeric(1))
  if (all(is.na(v))) stop("Could not find zero-inflation probability in hyperpar samples.")
  v
}

# ------------------------------------------------------------------------------
# Extract Negative Binomial 'size' (θ) from hyperparameter draws.
# Heuristic: if name suggests log-scale or value <= 0, exponentiate.
# ------------------------------------------------------------------------------
.get_nb_size_draws <- function(samples) {
  v <- vapply(samples, function(s) {
    hp <- s$hyperpar; if (is.null(hp) || !length(hp)) return(NA_real_)
    nm <- names(hp); idx <- which(grepl("nbinomial|overdisp|size", nm, ignore.case = TRUE))
    if (!length(idx)) return(NA_real_)
    val <- unname(hp[idx[1]])
    if (!is.na(val) && (val <= 0 || grepl("log", nm[idx[1]], ignore.case = TRUE))) exp(val) else val
  }, numeric(1))
  if (all(is.na(v))) stop("Could not find NB 'size' (overdispersion) in hyperpar samples.")
  v
}

# ------------------------------------------------------------------------------
# Main function: posterior-predictive coverage for INLA count families.
# Simulates integer predictive draws y~ (Poisson/NB/Binomial/Z-INF variants)
# and computes discrete prediction intervals + coverage at 'level'.
# ------------------------------------------------------------------------------
inla_predictive_coverage <- function(res, y, family = NULL, indices = NULL,
                                     level = 0.95, nsim = 2000, use_shortest = FALSE,
                                     E = NULL, Ntrials = NULL) {
  stopifnot(inherits(res, "inla"))
  if (is.null(res$summary.linear.predictor))
    stop("Fit must use control.predictor = list(compute = TRUE).")
  
  # Detect family (first if multiple). Normalize to lower-case string.
  if (is.null(family)) {
    fam <- res$.args$family
    family <- if (length(fam) == 1) as.character(fam) else as.character(fam)[1]
  }
  family <- tolower(family)
  
  # Select rows to evaluate. Must supply y of matching length.
  n_all <- nrow(res$summary.linear.predictor)
  if (is.null(indices)) indices <- seq_len(n_all)
  n <- length(indices)
  if (length(y) != n) stop("Length(y) must equal length(indices).")
  
  # Draw from joint posterior (latent + hyperparameters)
  samps <- INLA::inla.posterior.sample(n = nsim, result = res)
  
  # Map summary row indices -> latent vector positions (robust to config)
  latent_pos <- .map_predictor_positions(res, indices, samps)
  
  # Extract η samples at those positions (nsim x n)
  eta_mat <- .get_eta_samples_by_pos(samps, latent_pos)  # nsim x n
  
  # Convenience predicate for family membership
  fam_is <- function(x) family %in% x
  # Container for integer predictive draws
  ysim <- matrix(NA_integer_, nsim, n)
  
  # For log-link count families, invert link to μ = exp(η).
  # If Poisson with E= exposure was used, multiply μ by E (and only then).
  if (fam_is(c("poisson", "zeroinflatedpoisson0", "zeroinflatedpoisson1",
               "zeroinflatednbinomial0", "zeroinflatednbinomial1", "nbinomial"))) {
    mu_mat <- exp(eta_mat)
    
    # IMPORTANT: E is only valid for Poisson models fit via the E= argument.
    if (!is.null(E)) {
      if (!identical(family, "poisson")) {
        stop("Argument 'E' should only be supplied for Poisson models fit via the E= argument.")
      }
      if (length(E) != n) stop("Length(E) must equal length(indices).")
      mu_mat <- sweep(mu_mat, 2, E, `*`)
    }
  }
  
  # Simulate predictive draws by family
  if (family == "poisson") {
    for (j in seq_len(n)) ysim[, j] <- rpois(nsim, lambda = mu_mat[, j])
    
  } else if (fam_is(c("zeroinflatedpoisson0", "zeroinflatedpoisson1"))) {
    # Structural zero with probability π; else Poisson(μ)
    pi_draws <- .get_zip_pi_draws(samps)  # structural-zero probability per draw
    for (j in seq_len(n)) {
      z <- rbinom(nsim, 1L, pi_draws)    # 1 => structural zero
      nz <- which(z == 0L)
      yj <- integer(nsim)
      if (length(nz)) yj[nz] <- rpois(length(nz), lambda = mu_mat[nz, j])
      ysim[, j] <- yj
    }
    
  } else if (fam_is(c("zeroinflatednbinomial0", "zeroinflatednbinomial1"))) {
    # Structural zero with probability π; else NB(size=θ, mean=μ)
    pi_draws   <- .get_zip_pi_draws(samps)
    size_draws <- .get_nb_size_draws(samps)  # NB size (>0) per draw
    for (j in seq_len(n)) {
      z <- rbinom(nsim, 1L, pi_draws)
      nz <- which(z == 0L)
      yj <- integer(nsim)
      if (length(nz)) {
        mu_j <- mu_mat[nz, j]
        yj[nz] <- rnbinom(length(nz), size = size_draws[nz], mu = mu_j)
      }
      ysim[, j] <- yj
    }
    
  } else if (family == "nbinomial") {
    # Plain Negative Binomial (no zero inflation): NB(size=θ, mean=μ)
    size_draws <- .get_nb_size_draws(samps)  # vector length nsim
    for (j in seq_len(n)) {
      ysim[, j] <- rnbinom(nsim, size = size_draws, mu = mu_mat[, j])
    }
    
  } else if (family == "binomial") {
    # Binomial: invert logit to p; require Ntrials per observation
    p_mat <- plogis(eta_mat)      # default link = logit
    if (is.null(Ntrials)) {
      if (all(y %in% c(0L, 1L))) Ntrials <- rep(1L, n)
      else stop("Binomial with counts > 1 requires Ntrials (vector, length = length(indices)).")
    }
    if (length(Ntrials) != n) stop("Length(Ntrials) must equal length(indices).")
    for (j in seq_len(n)) ysim[, j] <- rbinom(nsim, size = Ntrials[j], prob = p_mat[, j])
    
  } else {
    # Guardrail: unsupported family
    stop("Unsupported family: ", family,
         ". Supported: poisson, nbinomial, binomial, ",
         "zeroinflatedpoisson0/1, zeroinflatednbinomial0/1")
  }
  
  # Build discrete prediction intervals at 'level' and compute coverage
  # Equal-tailed (type=1) by default; optional shortest high-mass set.
  alpha <- 1 - level
  if (!use_shortest) {
    L <- apply(ysim, 2, quantile, probs = alpha/2, type = 1, names = FALSE)
    U <- apply(ysim, 2, quantile, probs = 1 - alpha/2, type = 1, names = FALSE)
  } else {
    L <- U <- integer(n)
    for (j in seq_len(n)) {
      tab <- table(ysim[, j])
      pmf <- as.numeric(tab) / nsim; names(pmf) <- names(tab)
      lims <- .shortest_contiguous_interval(pmf, level = level)
      L[j] <- lims[1]; U[j] <- lims[2]
    }
  }
  
  # Coverage indicator and summary proportion
  covered <- (y >= L) & (y <= U)
  list(
    family   = family,
    level    = level,
    coverage = mean(covered),
    covered  = covered,
    L = as.integer(L), U = as.integer(U),
    indices  = indices
  )
}
