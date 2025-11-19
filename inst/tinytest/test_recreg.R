
# cifreg versus crr-cmprsk
library("tinytest")

# recreg versus recurrentMarginal for non-parametric 
test_recregMarginal <- function() {
 library(mets)
 data(hfactioncpx12)
 hf <- hfactioncpx12
 hf <- tie.breaker(hf,cause=c(0,1,2),cens.code=4)

 ngl <- recreg(Event(entry,time,status)~strata(treatment)+cluster(id),data=hf,
 cause=1,death.code=2,cens.model=~strata(treatment))
 ### 
 dd <- data.frame(treatment=levels(hf$treatment),id=1,status=0,time=0,entry=0)
 pngl <- predict(ngl,dd,times=3,se=1)
 pgl <- cbind(pngl$cumhaz,pngl$se.cumhaz)
 pgl

 ### 
 meann <- recurrentMarginal(Event(entry,time,status)~strata(treatment)+cluster(id),data=hf,cause=1,death.code=2)
 smeann <- summary(meann,time=3)
 smeann <- rbind(smeann$pbaseci[[1]],smeann$pbaseci[[2]])
 smeann
 pgl
 dif <- pgl - smeann[,2:3]
 dif

 ## estimate (same) and standard errors (not to far)
 expect_true( ((sum(abs(dif[,1])) < 0.000001) & (sum(abs(dif[,2])) < 0.001)) )
}
test_recregMarginal()


# recregIPCW versus recurrentMarginal for non-parametric 
test_recregMarginalIPCW <- function() {
 library(mets)
 data(hfactioncpx12)
 hf <- hfactioncpx12
 ## break all-ties for cause=0,1,2
 hf <- tie.breaker(hf,cause=c(0,1,2),cens.code=4)

 ngl <- recregIPCW(Event(entry,time,status)~treatment+cluster(id),data=hf,cause=1,death.code=2,time=3,
		   cens.model=~strata(treatment))
 dd <- data.frame(treatment=levels(hf$treatment),id=1)
 pngl <- predict(ngl,dd,times=3,se=1)
 ### 
 meann <- recurrentMarginal(Event(entry,time,status)~strata(treatment)+cluster(id),data=hf,cause=1,death.code=2)
 smeann <- summary(meann,time=3)
 smeann <- rbind(smeann$pbaseci[[1]],smeann$pbaseci[[2]])
 smeann
 pngl
 dif <- pngl[,1:2] - smeann[,2:3]
 dif

 ## estimate (same) and standard errors (not to far)
 expect_true( ((sum(abs(dif[,1])) < 0.00001) & (sum(abs(dif[,2])) < 0.001)) )
}
test_recregMarginalIPCW()



