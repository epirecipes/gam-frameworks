# Predictive coverage for gamlss models (counts)
# ----------------------------------------------
# - Supports parametric (with parameter uncertainty) and plug-in predictions
# - Works in- or out-of-sample (via newdata)
# - Handles multi-parameter count families (mu, sigma, nu, tau)
# - Produces discrete equal-tailed or shortest HDR intervals and coverage

suppressPackageStartupMessages({
  library(gamlss)
  library(gamlss.dist)
})

.shortest_discrete_set <- function(pmf, level = 0.95) {
  ord <- order(pmf, decreasing = TRUE)
  k <- which(cumsum(pmf[ord]) >= level)[1]
  S <- sort(as.integer(names(pmf))[ord[seq_len(k)]])
  c(min(S), max(S))
}

.get_link_name_gamlss <- function(fit, p) {
  val <- tryCatch(fit[[paste0(p, ".link")]], error = function(e) NULL)
  if (is.character(val) && length(val) == 1) return(val)
  fam <- fit$family
  if (!is.null(fam) && !is.function(fam)) {
    val <- tryCatch(fam[[paste0(p, ".link")]], error = function(e) NULL)
    if (is.character(val) && length(val) == 1) return(val)
  }
  stop("Could not determine the '", p, "' link for this GAMLSS fit.")
}

.extract_link_fit_se <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) {
    cn <- colnames(x)
    fit <- if (!is.null(cn) && "fit"    %in% cn) x[, "fit",    drop = TRUE] else x[, 1, drop = TRUE]
    se  <- if (!is.null(cn) && "se.fit" %in% cn) x[, "se.fit", drop = TRUE] else rep(0, nrow(x))
    return(list(fit = as.numeric(fit), se = as.numeric(se)))
  }
  if (is.list(x) && ("fit" %in% names(x) || "se.fit" %in% names(x))) {
    fit <- x[["fit"]]
    se  <- x[["se.fit"]]; if (is.null(se)) se <- rep(0, length(fit))
    return(list(fit = as.numeric(fit), se = as.numeric(se)))
  }
  if (is.numeric(x)) return(list(fit = as.numeric(x), se = rep(0, length(x))))
  stop("Unexpected structure returned by predictAll(): ", paste(class(x), collapse = "/"))
}

gamlss_predictive_coverage <- function(
    fit,
    y,
    newdata = NULL,
    level = 0.95,
    nsim = 2000,
    method = c("parametric","plug-in"),
    use_shortest = FALSE,
    data_fit = NULL
) {
  stopifnot(inherits(fit, "gamlss"))
  method <- match.arg(method)
  
  fam_obj  <- fit$family
  fam_name <- tryCatch(fam_obj$family, error = function(e) NULL)
  if (is.null(fam_name)) fam_name <- as.character(fam_obj)[1]
  rng_name <- paste0("r", fam_name)
  if (!exists(rng_name, where = asNamespace("gamlss.dist"), inherits = FALSE)) {
    stop("No RNG found for family '", fam_name, "' (looking for gamlss.dist::", rng_name, "()).")
  }
  rng_fun <- get(rng_name, envir = asNamespace("gamlss.dist"))
  
  params <- intersect(c("mu","sigma","nu","tau"), fit$parameters)
  if (!length(params)) stop("Could not detect mu/sigma/nu/tau parameters for this family.")
  
  n <- if (!is.null(newdata)) nrow(newdata) else length(y)
  if (length(y) != n) stop("Length(y) must match the number of rows being evaluated.")
  
  link_names <- vapply(params, function(p) .get_link_name_gamlss(fit, p), character(1))
  invlink <- lapply(link_names, function(ln) make.link.gamlss(ln)$linkinv)
  names(invlink) <- params
  
  get_response_params <- function() {
    pa <- predictAll(fit, newdata = newdata, type = "response", output = "list")
    out <- list(); for (p in params) out[[p]] <- as.numeric(pa[[p]])
    out
  }
  
  get_link_params_with_se <- function() {
    # main attempt: with se.fit
    pa <- tryCatch(
      predictAll(fit, newdata = newdata, type = "link",
                 se.fit = TRUE, use.weights = TRUE, data = data_fit, output = "list"),
      error = function(e) {
        predictAll(fit, newdata = newdata, type = "link",
                   se.fit = TRUE, use.weights = FALSE, data = data_fit, output = "list")
      }
    )
    # fallback means (no se) for NA repair
    pa_fit_only <- predictAll(fit, newdata = newdata, type = "link",
                              se.fit = FALSE, output = "list")
    
    link_mean <- vector("list", length(params)); names(link_mean) <- params
    link_se   <- vector("list", length(params)); names(link_se)   <- params
    
    for (p in params) {
      ext <- .extract_link_fit_se(pa[[p]])
      m <- as.numeric(ext$fit)
      s <- as.numeric(ext$se)
      
      # ensure lengths
      if (length(m) != n) m <- rep(m, length.out = n)
      if (length(s) != n) s <- rep(s, length.out = n)
      
      # backfill non-finite means with no-SE predictions
      m_bad <- !is.finite(m)
      if (any(m_bad)) {
        m0 <- as.numeric(pa_fit_only[[p]])
        if (length(m0) != n) m0 <- rep(m0, length.out = n)
        m[m_bad] <- m0[m_bad]
      }
      # any still bad -> set to median (or 0 if all bad)
      if (any(!is.finite(m))) {
        repl <- if (any(is.finite(m))) stats::median(m[is.finite(m)]) else 0
        m[!is.finite(m)] <- repl
      }
      
      # sanitize SEs: NA/NaN/Inf/negative -> 0 (fixed at mean)
      s[!is.finite(s) | s < 0] <- 0
      
      link_mean[[p]] <- m
      link_se[[p]]   <- s
    }
    list(mean = link_mean, se = link_se)
  }
  
  # parameter draws on response scale (nsim x n)
  par_draws <- switch(method,
                      "plug-in" = {
                        rp <- get_response_params()
                        lapply(rp, function(v) matrix(rep(v, each = nsim), nrow = nsim, ncol = n))
                      },
                      "parametric" = {
                        lp <- get_link_params_with_se()
                        draws <- list()
                        for (p in params) {
                          m <- lp$mean[[p]]; s <- lp$se[[p]]
                          # guard: ensure finite, non-negative sd; finite means
                          m[!is.finite(m)] <- 0
                          s[!is.finite(s) | s < 0] <- 0
                          eta <- matrix(
                            rnorm(nsim * n, mean = rep(m, each = nsim), sd = rep(s, each = nsim)),
                            nrow = nsim, ncol = n
                          )
                          # inverse link to response scale
                          val <- invlink[[p]](eta)
                          
                          # enforce parameter support by link type
                          eps <- 1e-12
                          ln <- link_names[[p]]
                          if (grepl("^log", ln, ignore.case = TRUE)) {
                            # strictly positive
                            val[!is.finite(val) | val <= 0] <- eps
                          } else if (grepl("logit|cloglog|probit|cauchit|loglog", ln, ignore.case = TRUE)) {
                            # probabilities
                            val[!is.finite(val)] <- 0.5
                            val <- pmin(pmax(val, eps), 1 - eps)
                          } else {
                            # generic: replace non-finite with column medians (or 0)
                            bad <- !is.finite(val)
                            if (any(bad)) {
                              col_med <- apply(val, 2, function(col) {
                                ok <- is.finite(col); if (any(ok)) stats::median(col[ok]) else 0
                              })
                              val[bad] <- rep(col_med, each = nsim)[bad]
                            }
                          }
                          draws[[p]] <- val
                        }
                        draws
                      }
  )
  
  # simulate integer outcomes
  ysim <- matrix(NA_integer_, nsim, n)
  for (j in seq_len(n)) {
    args <- list(n = nsim)
    if ("mu"    %in% params) args$mu    <- par_draws$mu[, j]
    if ("sigma" %in% params) args$sigma <- par_draws$sigma[, j]
    if ("nu"    %in% params) args$nu    <- par_draws$nu[, j]
    if ("tau"   %in% params) args$tau   <- par_draws$tau[, j]
    ysim[, j] <- do.call(rng_fun, args)
  }
  
  # prediction intervals + coverage
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
    family   = fam_name,
    method   = method,
    level    = level,
    coverage = mean(covered),
    covered  = covered,
    L = as.integer(L), U = as.integer(U)
  )
}
