dnbmix = function(y, mu, size, log = FALSE)
{
    stopifnot(length(y) == 1)
    stopifnot(length(mu) == length(size))

    # Individual log probabilities
    lpp = dnbinom(y, size = size, mu = mu, log = TRUE)
    # log(sum(exp(individual log probabilities)))
    lse_lp = matrixStats::logSumExp(lpp)
    # get log((1/N)*sum(exp(...))) by subtracting log(N)
    lp = lse_lp - log(length(mu))
    
    if (log) {
        return (lp)
    } else {
        return (exp(lp))
    }
}

pnbmix = function(y, mu, size, log = FALSE)
{
    stopifnot(length(y) == 1)
    stopifnot(length(mu) == length(size))

    # Individual log probabilities
    lPP = pnbinom(y, size = size, mu = mu, log = TRUE)
    # log(sum(exp(individual log probabilities)))
    lse_lP = matrixStats::logSumExp(lPP)
    # get log((1/N)*sum(exp(...))) by subtracting log(N)
    lP = lse_lP - log(length(mu))
    
    if (log) {
        return (lP)
    } else {
        return (exp(lP))
    }
}

# Calculate quadratic score penalty for negative binomial mix
qspnbmix = function(mu, size)
{
    lim = min(50000, ceiling(max(mu) * 5 + 1))
    ps = matrixStats::logSumExp(2 * sapply(0:lim, dnbmix, mu = mu, size = size, log = TRUE))
    ps1 = -Inf
    i = lim
    while (ps1 != ps) {
        i = i + 1
        ps1 = ps
        ps = matrixStats::logSumExp(ps, 2 * dnbmix(i, mu = mu, size = size, log = TRUE))
    }
    return (exp(ps))
}


# Some scoring metrics for count data. 
# y is scalar, has to be added/summed up across multiple points.
# Log score
log_score_nb = function(y, mu, size)
{
    return (-dnbmix(y, mu, size, log = TRUE))
}

# Quadratic score
quad_score_nb = function(y, mu, size)
{
    2 * dnbmix(y, mu, size) - qspnbmix(mu, size)
}

is_wholenumber = function(x, tol = .Machine$double.eps^0.5)
{
    abs(x - round(x)) < tol
}

# Ranked probability score
rp_score_nb = function(y, mu, size)
{
    stopifnot(length(y) == 1)
    stopifnot(is_wholenumber(y))
    stopifnot(length(mu) == length(size))
    y = round(y)

    lim = max(ceiling(max(mu) * 5 + 1), y + 50)
    k0 = seq_len(y) - 1
    k1 = y:lim
    
    ps = sum(vapply(k0, pnbmix, numeric(1), mu, size)^2) + 
        sum((vapply(k1, pnbmix, numeric(1), mu, size) - 1)^2)
    
    ps1 = 0
    i = lim
    while (ps1 != ps) {
        i = i + 1
        ps1 = ps
        ps = ps + (pnbmix(i, mu, size) - 1)^2
    }
    
    return (ps)
}
