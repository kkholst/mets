# cifreg versus crr-cmprsk
library("tinytest")

test_cluster <- function() {
  library(mets)
  set.seed(100)
  data(bmt)
  bmt$id <- 1:nrow(bmt)
  bmt2 <- rbind(bmt,bmt)
  bmt2$time <-  bmt2$time + runif(nrow(bmt2)) * 0.001

  fg <- cifregFG(Event(time,cause)~tcell+age+platelet,bmt)
  estimate(fg)
  ###
  fg2 <- cifregFG(Event(time,cause)~tcell+age+platelet+cluster(id),bmt2)
  estimate(fg2)
  d1 <- max(abs(fg2$se.coef-fg$se.coef))

  fg <- cifreg(Event(time,cause)~tcell+age+platelet,bmt)
  estimate(fg)
  ###
  fg2 <- cifreg(Event(time,cause)~tcell+age+platelet+cluster(id),bmt2)
  estimate(fg2)
  d2 <- max(abs(fg2$se.coef-fg$se.coef))

  fg <- phreg(Event(time,cause!=0)~tcell+age+platelet,bmt)
  estimate(fg)
  ###
  fg2 <- phreg(Event(time,cause!=0)~tcell+age+platelet+cluster(id),bmt2)
  estimate(fg2)
  d3 <- max(abs(fg2$var-fg$var))

  fg <- binreg(Event(time,cause!=0)~tcell+age+platelet,bmt,time=30)
  estimate(fg)
  ###
  fg2 <- binreg(Event(time,cause!=0)~tcell+age+platelet+cluster(id),bmt2,time=30)
  estimate(fg2)
  d4 <- max(abs(fg2$se.coef-fg$se.coef))

  dfactor(bmt) <- tcell.f~tcell
  dfactor(bmt2) <- tcell.f~tcell
###
  fg <- binregATE(Event(time,cause!=0)~tcell.f+age+platelet,bmt,time=30,
		  treat.model=tcell.f~age+platelet,cens.model=~strata(tcell,platelet))
  sfg <- summary(fg)
  ###
  fg2 <- binregATE(Event(time,cause!=0)~tcell.f+age+platelet+cluster(id),bmt2,time=30,
		  treat.model=tcell.f~age+platelet,cens.model=~strata(tcell,platelet))
  sfg2 <- summary(fg2)
  ateg <- max(abs(sfg$ateG[,2]-sfg2$ateG[,2] ))
  ateDR <- max(abs(sfg$ateDR[,2]-sfg2$ateDR[,2] ))

  dd <- max(c(d1,d2,d3,d4,ateg,ateDR))

  expect_true( ((d3 < 0.0011) & (max(c(d1,d2,d4,ateg,ateDR)) < 0.001)) )
} 
test_cluster()


