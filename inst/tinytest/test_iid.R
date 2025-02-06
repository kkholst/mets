
## Influence functions, robust standard errors"
library("tinytest")

test_phreg.iid <- function() {
  data(TRACE)
  out1 <- phreg(Surv(time,status==9)~vf+chf, data=TRACE)
  expect_true(norm(estimate(out1)$vcov - solve(-out1$hessian))<1)
}
test_phreg.iid()

test_twinlm.iid <- function() {
  set.seed(1)
  d <- twinsim(100,b1=c(1,-1),b2=c(),acde=c(1,1,0,1))
  l <- twinlm(y ~ 1, data=d, DZ="DZ", zyg="zyg", id="id")
  expect_true(norm(estimate(l)$vcov - vcov(l))<1)
}

test_biprobit <- function() {
  set.seed(1)
  d <- twinsim(100,b1=c(1,-1),b2=c(),acde=c(1,1,0,1))
  l <- biprobit(y0 ~ x1, id="id", rho=~x1, data=d)
  expect_true(norm(estimate(l)$vcov - vcov(l))<1)
}
test_biprobit()
