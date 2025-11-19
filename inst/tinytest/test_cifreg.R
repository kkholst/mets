# cifreg versus crr-cmprsk
library("tinytest")

test_cifreg <- function() {
  library(mets)
  set.seed(100)
  data(bmt)
  bmt$time <-  bmt$time + runif(nrow(bmt)) * 0.001

  fg <- cifregFG(Event(time,cause)~tcell+age+platelet,bmt)

  if (requireNamespace("cmprsk")) {
    library("cmprsk")
    mm <- model.matrix(~ tcell + age + platelet, bmt)[, -1]
    cr <- with(bmt, crr(time,cause,mm))
    mm <-  cbind(cr$coef,
                 diag(cr$var)^.5,
                 fg$coef,
                 fg$se.coef,
                 cr$coef - fg$coef,
                 diag(cr$var)^.5 - fg$se.coef)
    mm
    ## estimate (same) and standard errors (close close)
    expect_true( ((sum(abs(mm[, 5]))) < 0.0001) & ((sum(abs(mm[,6])))<0.0001)   )
  }
}
test_cifreg()



