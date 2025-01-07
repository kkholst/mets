
context("binreg simple-icpw versus prodlim")

test_that("binreg-cif", {
data(bmt)
bmt$time <-  bmt$time+runif(nrow(bmt))*0.001
out <- binreg(Event(time,cause)~-1+factor(platelet),bmt,time=50,cens.model=~strata(platelet))
summary(out)
pp <- predict(out,data.frame(platelet=0:1,tcell=0:1))

cif1 <- cif(Event(time,cause)~strata(platelet),bmt)
cifs <- cbind(cif1$cumhaz[,1],vecAllStrata(cif1$cumhaz[,2],cif1$strata,cif1$nstrata))
pcif50 <- cpred(cifs,50)[,2]

se.cifs <- cbind(cif1$cumhaz[,1],vecAllStrata(cif1$se.cumhaz[,2],cif1$strata,cif1$nstrata))
se.pcif50 <- cpred(se.cifs,50)[,2]

## estimate (same) and standard errors (close)
expect_true( ((sum(abs(pcif50-pp[,1])))<0.0001) & (sum(abs(pp[,2]-cpred(se.cifs,50)[,2])) < 0.0002))
})

