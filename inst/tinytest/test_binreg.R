# binreg simple-icpw versus prodlim
library("tinytest")

test_binreg_cif <- function() {

  library(mets)
  data(bmt)
  bmt$time <-  bmt$time + runif(nrow(bmt)) * 0.001
  out <- binreg(Event(time, cause) ~ -1+factor(platelet),
                bmt, time=30, cens.model=~strata(platelet))

  summary(out)
  pp <- predict(out, data.frame(platelet = 0:1, tcell = 0:1))
  pp

  cif1 <- cif(Event(time,cause) ~ strata(platelet), bmt)
  scif <- summary(cif1,time=30)
  cif50 <- rbind(scif$pbaseci[[1]],scif$pbaseci[[2]])
  cif50
  dif <- cif50[,2:3] - pp[,1:2]
  dif

  ## estimate (same) and standard errors (close)
  expect_true( ((sum(abs(dif[,1])) < 0.0001) & (sum(abs(dif[,2])) < 0.0002)) )
}
test_binreg_cif()




