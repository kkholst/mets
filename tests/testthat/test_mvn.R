
context("Censored regression (loglikMVN)")

cens1 <- function(threshold,type='right') {
  function(x) {
    x <- unlist(x)
    if (type=='left')
      return( survival::Surv(pmax(x,threshold), x>=threshold, type='left') )
      return ( survival::Surv(pmin(x,threshold), x<=threshold) )
  }
}

m0 <- lvm() |>
  covariance(y1 ~ y2, value='r') |>
  regression(y1 + y2 ~ x)
m0 <-  transform(m0, s1 ~ y1, cens1(-2, 'left')) |>
  transform(s2 ~ y2, cens1(2,  'right'))

test_that("phreg.iid", {
  d <- sim(m0, 500, p=c(r=0.9), seed=1)
  head(d)
  d
  m3 <- lvm() |>
    regression(y1 + s2 ~ x) |>
    covariance(y1 ~ s2, constrain=TRUE, rname='z')

  e3 <- estimate(m3, d)
  expect_true(mean(score(e3)**2)<1e-2)
})
