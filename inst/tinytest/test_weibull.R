library(tinytest)

## rate parametrization
set.seed(1)
rate <- 0.2
shape <- 3
d <- data.frame(time=rweibullcox(1e4, rate, shape), status=TRUE)
p <- phreg_weibull(Surv(time, status) ~ 1, data=d)
est <- coef(p) |> exp()
expect_true(mean((c(rate, shape) - est)^2) < 1e-2)

## shape parameter regression
set.seed(1)
n <- 1e4
a <- rbinom(n, 1, 0.5)
rate <- 0.2
shape <- 2+a
d <- data.frame(time=rweibullcox(n, rate, shape), status=TRUE)
p <- phreg_weibull(Surv(time, status) ~ 1, ~ factor(a) - 1, data = d)
est <- coef(p) |> exp()
expect_true(mean((c(rate, 2, 3) - est)^2) < 1e-2)

## delayed entry, right-censoring
set.seed(1)
n <- 2e4
rate <- 0.2
shape <- 2
entry <- runif(n)
t0 <- rweibullcox(n, rate, shape)
cens <- rweibullcox(n, 0.5, 1)
time <- pmin(t0, cens)
status <- t0<=cens
d0 <- data.frame(time = time, status = status, entry=entry)
d <- subset(d0, time>entry)

p0 <- phreg_weibull(Surv(time, status) ~ 1, ~ 1, data = d0)
est <- coef(p0) |> exp()
expect_true(mean((c(rate, shape) - est)^2) < 1e-2)

p.wrong <- phreg_weibull(Surv(time, status) ~ 1, ~ 1, data = d)
est <- coef(p.wrong) |> exp()
expect_false(mean((c(rate, shape) - est)^2) < 1e-2)

p <- phreg_weibull(Surv(entry, time, status) ~ 1, ~ 1, data = d)
est <- coef(p) |> exp()
expect_true(mean((c(rate, shape) - est)^2) < 1e-2)
