library(data.table)
library(mgcv)

dt = data.table(t = 0:48, y = sin((0:48) * pi/24) + rnorm(49, 0, 0.1))
plot(dt)

dt[, ym1 := c(0, head(y, -1))]

fit = gam(y ~ offset(ym1) + s(t), data = dt)

plot(fit)
summary(fit)




dt
##

dt = data.table(t = 0:60, y = sin((0:60) * pi/24) + rnorm(61, 0, 0.1))
plot(dt)


r = 0.5

system.time(
for (r in c(0, 0.2, 0.4, 0.6, 0.8, 1)) {
    dt2 = data.table::copy(dt)
    nsamp = round(r * nrow(dt2))
    ind = sample.int(nrow(dt2), size = nsamp)
    dt2$t[ind] = sample(dt2$t[ind], nsamp)
    dt2 = dt2[order(t)]
    
    fit = gam(y ~ s(t), data = dt2)
    print(plot(fit))
}
)



for (r in c(0, 0.2, 0.4, 0.6, 0.8, 1)) {
    fit = gam(y ~ s(t), data = dt, sp = 1000000*r^30)
    print(plot(fit))
}



# Extreme n^2 edition
system.time(
for (r in c(0, 0.2, 0.4, 0.6, 0.8, 1)) {
    dt2 = data.table::copy(dt)
    nr0 = nrow(dt2)
    dt2 = rep(list(dt2), nr0)
    for (j in seq_along(dt2)) {
        dt2[[j]]$y = shift(dt2[[j]]$y, j - 1, type = "cyclic")
    }
    dt2 = rbindlist(dt2)
    dt2[, w := c(rep(1 - r + r/nr0, nr0), rep(r/nr0, nr0^2 - nr0))]

    fit = gam(y ~ s(t), data = dt2, weights = dt2$w)
    print(plot(fit))
}
)



# Random Kn edition
K = 1
system.time(
for (r in c(0, 0.2, 0.4, 0.6, 0.8, 1)) {
    dt2 = rbindlist(rep(list(dt), K))
    
    nsamp = round(r * nrow(dt2))
    ind = sample.int(nrow(dt2), size = nsamp)
    dt2$y[ind] = sample(dt2$y[ind], nsamp)

    fit = gam(y ~ s(t), data = dt2, weights = rep(1/K, nrow(dt2)))
    print(plot(fit))
}
)


# k edition
system.time(
for (r in c(0, 0.1, 0.2, 0.3, 0.4, 1)) {
    fit = gam(y ~ s(t, k = 2 + 20 * r), data = dt)
    print(plot(fit))
}
)
