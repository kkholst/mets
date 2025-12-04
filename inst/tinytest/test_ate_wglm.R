library("tinytest")
library(mets)
library(riskRegression)
## break ties
set.seed(100)
data(bmt)
bmt$time <- bmt$time+runif(408)*.001
dfactor(bmt)  <-  ~tcell+platelet
bmt$id <- 1:408

test_ate_wglm <- function() { ## {{{
brs <- logitIPCW(Event(time,cause)~tcell+platelet+age,bmt,time=50,cause=1)
estimate(brs)
###
lbrs <- logitIPCWATE(Event(time,cause)~tcell+platelet+age,bmt,time=50,cause=1,treat.model=tcell~platelet+age)
###     
eCR.wglm <- wglm(Surv(time,cause)~ tcell.f+platelet+age,
	      formula.censor =~+1, times = 50, data = bmt, cause = 1, product.limit = TRUE)
summary(eCR.wglm)
estimate(lbrs)$coefmat
estimate(brs)$coefmat
####
awglm <- ate(eCR.wglm, data = bmt, times = 50, treatment = tcell.f~platelet+age, 
	     verbose = FALSE,estimator=c("G-formula","AIPTW","IPTW"))
aa <- confint(awglm)
aa$diffRisk
aa$meanRisk
summary(lbrs)$ateG
summary(lbrs)$ateDR
gg <- aa$meanRisk[2:1,c(4,5)]
dr <- aa$meanRisk[6:5,c(4,5)]
ddr  <-  dr- summary(lbrs)$ateDR[1:2,1:2]
estdr <- ddr[,1]
dgg  <-  gg-summary(lbrs)$ateG[1:2,1:2]
###
Giid <- max(abs(cbind(awglm$iid$GFORMULA[[1]], awglm$iid$GFORMULA[[2]])- lbrs$riskG.iid[,c(2,1)]))
DRiid <- max(abs(cbind(awglm$iid$AIPTW[[1]], awglm$iid$AIPTW[[2]])- lbrs$riskDR.iid[,c(2,1)]))
dgdr <- max( abs(c(unlist(ddr),unlist(dgg))))
c(dgdr,Giid,DRiid)
###
dgg <- max( abs(unlist(dgg)))
ddr <- max( abs(c(unlist(ddr))))
c(dgg,max(abs(estdr)),ddr,Giid,DRiid)
expect_true((dgg < 0.00001) &  (max(abs(estdr)) < 0.00001) & Giid < 0.00001) 
## standard error for DR is missing a censoring term from ate function

brs <- logitIPCW(Event(time,cause)~tcell+platelet+age,bmt,time=50,cause=1,
       cens.model=~strata(tcell,platelet))
lbrs <- logitIPCWATE(Event(time,cause)~tcell+platelet+age,bmt,time=50,cause=1,
       treat.model=tcell~platelet+age,cens.model=~strata(tcell,platelet),se=TRUE)
###     
eCR.wglm <- wglm(Surv(time,cause)~ tcell.f+platelet+age,
	      formula.censor =~strata(tcell,platelet), times = 50, data = bmt, cause = 1, product.limit = TRUE)
d <- summary(eCR.wglm)
estimate(brs)$coefmat
estimate(lbrs)$coefmat
dcoef <- max(abs(d[[1]][,c(1,2)] - estimate(brs)$coefmat[,c(1,2)]))
###
awglm <- ate(eCR.wglm, data = bmt, times = 50, treatment = tcell.f~platelet+age, 
	     verbose = FALSE,estimator=c("G-formula","AIPTW","IPTW"))
summary(awglm)
aa <- confint(awglm)
aa$diffRisk
aa$meanRisk
summary(lbrs)$ateG
summary(lbrs)$ateDR
gg <- aa$meanRisk[2:1,c(4,5)]
dr <- aa$meanRisk[6:5,c(4,5)]
ddr  <-  dr- summary(lbrs)$ateDR[1:2,1:2]
ddr
estdr <- ddr[,1]
dgg  <-  gg-summary(lbrs)$ateG[1:2,1:2]
dgg
###
Giid <- max(abs(cbind(awglm$iid$GFORMULA[[1]], awglm$iid$GFORMULA[[2]])- lbrs$riskG.iid[,c(2,1)]))
DRiid <- max(abs(cbind(awglm$iid$AIPTW[[1]],awglm$iid$AIPTW[[2]])- lbrs$riskDR.iid[,c(2,1)]))
###
dgg <- max( abs(unlist(dgg)))
ddr <- max( abs(c(unlist(ddr))))
c(dgg,max(abs(estdr)),ddr,Giid,DRiid)
expect_true((dgg < 0.00001) &  (max(abs(estdr)) < 0.00001) & Giid < 0.00001) 
## standard error for DR is missing a censoring term from ate function

################################################################
## no censoring case ###########################################
################################################################
bmt$bin <- (bmt$cause %in% c(1,2))
## binregATE based 
lbrs <- logitATE(bin~tcell.f+platelet+age,bmt,treat.model=tcell.f~platelet+age)
## logitIPCWATE
llbrs <- logitATE(bin~tcell.f+platelet+age,bmt,treat.model=tcell.f~platelet+age,binreg=FALSE)
###     
eCR.wglm <- glm(bin~tcell.f+platelet+age,data=bmt,family=binomial())
## robust se's
estimate(eCR.wglm)
estimate(lbrs)$coefmat
estimate(llbrs)$coefmat
###
awglm <- ate(eCR.wglm, data = bmt, times = 50, treatment = tcell.f~platelet+age, verbose = FALSE,estimator=c("G-formula","AIPTW","IPTW"))
aa <- confint(awglm)
aa$diffRisk
aa$meanRisk
summary(lbrs)$ateG
summary(llbrs)$ateG
summary(lbrs)$ateDR
summary(llbrs)$ateDR
gg <- aa$meanRisk[1:2,c(4,5)]
dr <- aa$meanRisk[5:6,c(4,5)]
ddr  <-  dr- summary(lbrs)$ateDR[1:2,1:2]
dgg  <-  gg- summary(lbrs)$ateG[1:2,1:2]
#
###
Giid <- max(abs(cbind(awglm$iid$GFORMULA[[1]],awglm$iid$GFORMULA[[2]])-lbrs$riskG.iid[,c(1,2)]))
DRiid <- max(abs(cbind(awglm$iid$AIPTW[[1]],awglm$iid$AIPTW[[2]])-lbrs$riskDR.iid[,c(1,2)]))
c(Giid,DRiid)
Giid2 <- max(abs(cbind(awglm$iid$GFORMULA[[1]],awglm$iid$GFORMULA[[2]])-llbrs$riskG.iid[,c(2,1)]))
DRiid2 <- max(abs(cbind(awglm$iid$AIPTW[[1]],awglm$iid$AIPTW[[2]])-llbrs$riskDR.iid[,c(2,1)]))
c(Giid2,DRiid2)
dgdr <- max( abs(c(unlist(ddr),unlist(dgg))))
c(dgdr,Giid,DRiid)
dgg <- max( abs(unlist(dgg)))
ddr <- max( abs(c(unlist(ddr))))
c(dgg,ddr,Giid,DRiid)
expect_true((dgg < 0.00001) &  (ddr < 0.00001) & (max(c(Giid,DRiid)) < 0.00001)) 
} ## }}}
test_ate_wglm()



