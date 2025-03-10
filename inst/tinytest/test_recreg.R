
# cifreg versus crr-cmprsk
library("tinytest")

# recreg versus recurrentMarginal for non-parametric 
test_recregMarginal <- function() {
 library(mets)
 data(hfaction_cpx12)
 hf <- hfaction_cpx12
 dd <- data.frame(treatment=levels(hf$treatment),id=1)

 ngl <- recreg(Event(entry,time,status)~strata(treatment)+cluster(id),data=hf,
 cause=1,death.code=2,cens.model=~strata(treatment))
 ### 
 pngl <- predict(ngl,dd,times=3,se=1)
 pgl <- cbind(pngl$cumhaz,pngl$se.cumhaz)
 ### 
 meann <- recurrentMarginal(Event(entry,time,status)~strata(treatment)+cluster(id),data=hf,cause=1,death.code=2)
 smeann <- summary(meann,time=3)
 smeann <- rbind(smeann$pbaseci[[1]],smeann$pbaseci[[2]])
 smeann
 pgl
 dif <- pgl - smeann[,2:3]

 ## estimate (same) and standard errors (not to far)
 expect_true( ((sum(abs(dif[,1])) < 0.0001) & (sum(abs(dif[,2])) < 0.1)) )
}
test_recregMarginal()

