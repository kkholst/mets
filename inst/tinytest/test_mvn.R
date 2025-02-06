
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
