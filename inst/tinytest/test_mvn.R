
## Censored regression (loglikMVN)
library("tinytest")

cens1 <- function(threshold, type='right') {
  function(x) {
    x <- unlist(x)
    if (type=='left')
      return( survival::Surv(pmax(x,threshold), x>=threshold, type='left') )
      return ( survival::Surv(pmin(x,threshold), x<=threshold) )
  }
}

m0 <- lava::lvm() |>
  lava::covariance(y1 ~ y2, value='r') |>
  lava::regression(y1 + y2 ~ x)
m0 <-  transform(m0, s1 ~ y1, cens1(-2, 'left')) |>
  transform(s2 ~ y2, cens1(2,  'right'))

test_estimatenormal <- function() {
  d <- lava::sim(m0, 500, p=c(r=0.9), seed=1)

  m3 <- lava::lvm() |>
    lava::regression(y1 + s2 ~ x) |>
    lava::covariance(y1 ~ s2, constrain=TRUE, rname='z')

  e3 <- lava::estimate(m3, d)
  expect_true(mean(lava::score(e3)**2)<1e-2)
}
test_estimatenormal()


## pmvn

test_mvn <- function() {
  S <- diag(2)+0.5
  p2 <- mvtnorm::pmvnorm(lower = c(1, 1), upper = c(2, 2), sigma = S)
  p1 <- mvtnorm::pmvnorm(lower = 1, upper = 2, sigma = S[1])
  z2 <- pmvn(lower = c(1, 1), upper = c(2, 2), sigma = S)
  z1 <- pmvn(lower = 1, upper = 2, sigma = S[1])
  # check bivariate
  expect_equivalent(p2, z2)
  # check univariate
  expect_equivalent(p1, z1)
  # check vectorized versions
  z <- pmvn(upper=cbind(c(-3:3)), sigma = 1)
  expect_equivalent(z, pnorm(-3:3))
}
