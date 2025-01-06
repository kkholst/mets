
context("RMST binreg versus resmean.phreg, cif.yearslost")

test_that("resmeanIPCW", {
     set.seed(101)
     data(bmt); bmt$time <- bmt$time+runif(nrow(bmt))*0.001

     ### same as Kaplan-Meier for full censoring model 
     bmt$int <- with(bmt,strata(tcell,platelet))
     out <- resmeanIPCW(Event(time,cause!=0)~-1+int,bmt,time=30,
                                  cens.model=~strata(platelet,tcell),model="lin")
     estimate(out)
     out1 <- phreg(Surv(time,cause!=0)~strata(tcell,platelet),data=bmt)
     rm1 <- resmean.phreg(out1,times=30)
     expect_true( (sum(abs(rm1$intkmtimes[,3] - coef(out)))<0.0001) & 
		 (sum(abs(rm1$intkmtimes[,4] - out$se.coef))<0.05) )
})

     
test_that("cif.yearslost", {
     set.seed(101)
     data(bmt); bmt$time <- bmt$time+runif(nrow(bmt))*0.001

     ### same as cif integral for full censoring model 
     bmt$int <- with(bmt,strata(tcell,platelet))

     ## competing risks years-lost for cause 1  
     outc <- resmeanIPCW(Event(time,cause)~-1+int,bmt,time=30,cause=1,
                                 cens.model=~strata(platelet,tcell),model="lin")
     ## same as integrated cumulative incidence 
     rmc1 <- cif.yearslost(Event(time,cause)~strata(tcell,platelet),data=bmt,times=30,cause=1)
     expect_true( (sum(abs(rmc1$intkmtimes[,3] - coef(outc)))<0.0001) & 
  		 (sum(abs(rmc1$intkmtimes[,5] - outc$se.coef))<0.05) )

})



