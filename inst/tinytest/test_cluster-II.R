library("tinytest")
library(mets) 
###
set.seed(100)
###
data(hfactioncpx12)
hf <- hfactioncpx12
## break ties after cause and censoring
hf <- tie.breaker(hf,id="id",cause=c(1,0),cens.code=5)
hf$idn <- hf$id
hf$x <- as.numeric(hf$treatment) 
hf$z <- rnorm(741)[hf$id]
hf$nn <- 1:nrow(hf)
hf$idc <- 5+sample(400,741,TRUE)[hf$id]
###
hf$id2 <- sample(1:741,741)[hf$id]
hf$ido <- hf$id
hf$id <- NULL
###
hf2 <- dsort(hf,~id2+entry+time,return.order=TRUE)

test_phregClust <- function() {
## {{{  works  phreg 

####################################################################################

gl <- phreg(Event(entry,time,status==1)~treatment+z+cluster(idn),data=hf)
gl2 <-phreg(Event(entry,time,status==1)~treatment+z+cluster(idn),data=hf2)
d1 <- max(abs(estimate(gl)$coefmat - estimate(gl2)$coefmat ))
d2 <- max(abs(t(head(iid(gl))- head(iid(gl2)))))
###
gg <- survivalG(gl,hf,time=1,First=TRUE)
gg2 <- survivalG(gl2,hf2,time=1,First=TRUE)
d3 <- max(abs(gg$survivalG$coefmat- gg2$survivalG$coefmat))
d4 <- max(abs(gg$risk$coefmat- gg2$risk$coefmat))
d5 <- max(abs(t(head(gg$risk.iid)- head(gg2$risk.iid))))
###
b1 <- iidBaseline(gl,time=1)
b2 <- iidBaseline(gl2,time=1)
d6 <- max(abs(t(head(cbind(b1$base.iid,b1$beta.iid))- head(cbind(b2$base.iid,b2$beta.iid)))))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
expect_true( (d3< 0.000000001) & ( d5 < 0.000000001 ))
expect_true( (d4< 0.000000001) & ( d5 < 0.000000001 ))
expect_true( (d6< 0.000000001) )
c(d1,d2,d3,d4,d5,d6)


gl <- phreg(Event(entry,time,status==1)~treatment+z,data=hf)
gl2 <-phreg(Event(entry,time,status==1)~treatment+z,data=hf2)
d1 <- max(abs(estimate(gl)$coefmat - estimate(gl2)$coefmat ))
expect_true( (d1<0.0000000001) )
gg <- survivalG(gl,hf,time=1)
gg2 <- survivalG(gl2,hf2,time=1)
d2 <- max(abs(gg$survivalG$coefmat- gg2$survivalG$coefmat))
d3 <- max(abs(gg$risk$coefmat- gg2$risk$coefmat))
d4 <- max(abs(gg$ratio$coefmat- gg2$ratio$coefmat))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
expect_true( (d3< 0.000000001) & ( d4 < 0.000000001 ))
c(d1,d2,d3,d4,d5)
b1 <- iidBaseline(gl,time=1)


gl <- phreg(Event(entry,time,status==1)~treatment+z+cluster(nn),data=hf)
gl2 <-phreg(Event(entry,time,status==1)~treatment+z+cluster(nn),data=hf2)
d1 <- max(abs(estimate(gl)$coefmat - estimate(gl2)$coefmat ))
expect_true( (d1<0.0000000001) )
gg <- survivalG(gl,hf,time=1)
gg2 <- survivalG(gl2,hf2,time=1)
d2 <- max(abs(gg$survivalG$coefmat- gg2$survivalG$coefmat))
d3 <- max(abs(gg$risk$coefmat- gg2$risk$coefmat))
d4 <- max(abs(gg$ratio$coefmat- gg2$ratio$coefmat))
d5 <- max(abs(t(head(gg$risk.iid)- head(gg2$risk.iid))))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
expect_true( (d3< 0.000000001) & ( d4 < 0.000000001 ))
expect_true( (d5< 0.000000001) )
c(d1,d2,d3,d4,d5)

gl <- phreg(Event(entry,time,status==1)~treatment+z+cluster(nn),data=hf)
gl2 <-phreg(Event(entry,time,status==1)~treatment+z,data=hf)
d1 <- max(abs(estimate(gl)$coefmat - estimate(gl2)$coefmat ))
d1
expect_true( (d1<0.0000000001) )
gg <- survivalG(gl,hf,time=1)
## data her skal svare til fit 
gg2 <- survivalG(gl2,hf,time=1)
d2 <- max(abs(gg$survivalG$coefmat- gg2$survivalG$coefmat))
d2
d3 <- max(abs(gg$risk$coefmat- gg2$risk$coefmat))
d4 <- max(abs(gg$ratio$coefmat- gg2$ratio$coefmat))
d5 <- max(abs(gg$risk.iid-gg2$risk.iid))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
expect_true( (d3< 0.000000001) & ( d4 < 0.000000001 ))
expect_true( (d5< 0.000000001) )
c(d1,d2,d3,d4,d5)

###gg <- survivalG(gl,hf2,time=1); gg2 <- survivalG(gl2,hf2,time=1)
###d2 <- max(abs(gg$survivalG$coefmat- gg2$survivalG$coefmat))
###d2
######
###gg <- survivalG(gl,hf2,time=1); gg2 <- survivalG(gl2,hf,time=1)
###d2 <- max(abs(gg$survivalG$coefmat- gg2$survivalG$coefmat))
###cat("survivalG på data i anden order end fit")
###d2
######
###gl2 <-phreg(Event(entry,time,status==1)~treatment+z,data=hf2)
###gg <- survivalG(gl,hf,time=1); gg2 <- survivalG(gl2,hf2,time=1)
###d2 <- max(abs(gg$survivalG$coefmat- gg2$survivalG$coefmat))
###d2

## }}}
}
test_phregClust()

test_recregClust <- function() {
## {{{  works recreg, cifreg, 

gl <- recreg(Event(entry,time,status)~treatment+z+cluster(idn),data=hf,cause=1,death.code=2)
###
gl2 <- recreg(Event(entry,time,status)~treatment+z+cluster(idn),data=hf2,cause=1,death.code=2)
d1 <- max(abs(estimate(gl)$coefmat - estimate(gl2)$coefmat))
d2 <- max(abs(iid(gl)-iid(gl2)))
###
gg <- survivalG(gl,hf,time=1)
gg2 <- survivalG(gl2,hf2,time=1)
d3 <- max(abs(gg$risk$coefmat- gg2$risk$coefmat))
d5 <- max(abs(gg$risk.iid-gg2$risk.iid))
###
b1 <- iidBaseline(gl,time=1)
b2 <- iidBaseline(gl2,time=1)
diid <- c(max(abs(IC(gl)-IC(gl2))), max(abs(iid(gl)-iid(gl2))), max(abs(iid(gl,time=1)-iid(gl2,time=1))),
max(abs(IC(gl,time=1)-IC(gl2,time=1))),
max(abs(b1$base.iid- iid(gl,time=1))),
max(abs(b2$base.iid- iid(gl2,time=1))),
max(abs(b2$beta.iid- iid(gl2))),
max(abs(b1$beta.iid- iid(gl))))
expect_true( (max(diid)< 0.000000001))
###
d4 <- max(abs(t(head(cbind(b1$base.iid,b1$beta.iid))- head(cbind(b2$base.iid,b2$beta.iid)))))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
expect_true( (d3< 0.000000001) & ( d4 < 0.000000001 ))
expect_true( d5 < 0.000000001)
c(d1,d2,d3,d4,max(diid),d5)


gl <- recreg(Event(entry,time,status)~treatment+z+cluster(id2),data=hf,cause=1,death.code=2)
###
gl2 <- recreg(Event(entry,time,status)~treatment+z+cluster(id2),data=hf2,cause=1,death.code=2)
d1 <- max(abs(estimate(gl)$coefmat - estimate(gl2)$coefmat ))
d2 <- max(abs(t(head(iid(gl))- head(iid(gl2)))))
###
gg <- survivalG(gl,hf,time=1)
gg2 <- survivalG(gl2,hf2,time=1)
d3 <- max(abs(gg$risk$coefmat- gg2$risk$coefmat))
###
b1 <- iidBaseline(gl,time=1)
b2 <- iidBaseline(gl2,time=1)
diid <- c(max(abs(IC(gl)-IC(gl2))), max(abs(iid(gl)-iid(gl2))), max(abs(iid(gl,time=1)-iid(gl2,time=1))),
max(abs(IC(gl,time=1)-IC(gl2,time=1))),
max(abs(b1$base.iid- iid(gl,time=1))),
max(abs(b2$base.iid- iid(gl2,time=1))),
max(abs(b2$beta.iid- iid(gl2))),
max(abs(b1$beta.iid- iid(gl))))
expect_true( (max(diid)< 0.000000001))
#
d4 <- max(abs(t(head(cbind(b1$base.iid,b1$beta.iid))- head(cbind(b2$base.iid,b2$beta.iid)))))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
expect_true( (d3< 0.000000001) & ( d4 < 0.000000001 ))
c(d1,d2,d3,d4,max(diid))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
expect_true( (d3< 0.000000001) & ( d4 < 0.000000001 ))
expect_true( max(diid)< 0.000000001)



   ## censoring stratified after treatment 
     gl <- recreg(Event(entry,time,status)~treatment+cluster(idn),data=hf,
     cause=1,death.code=2,cens.model=~strata(treatment))
      gl2 <- recreg(Event(entry,time,status)~treatment+cluster(idn),data=hf2,
     cause=1,death.code=2,cens.model=~strata(treatment))
     d1 <- max(abs(estimate(gl)$coefmat - estimate(gl2)$coefmat ))
     d2 <- max(abs(head(iid(gl))- head(iid(gl2))))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
c(d1,d2)
 

####################################################################################


gl <- cifregFG(Event(time,status)~treatment+z+cluster(id2),data=hf,cause=1)
gl2 <- cifregFG(Event(time,status)~treatment+z+cluster(id2),data=hf2,cause=1)
d1 <- max(abs(estimate(gl)$coefmat - estimate(gl2)$coefmat )) 
class(gl)
d2 <- max(abs(iid(gl)- iid(gl2)))
head(iid(gl))
###
gg <- survivalG(gl,hf,time=1)
gg2 <- survivalG(gl2,hf2,time=1)
d3 <- max(abs(gg$risk$coefmat- gg2$risk$coefmat ))
###
b1 <- iidBaseline(gl,time=1)
b2 <- iidBaseline(gl2,time=1)
diid <- c(max(abs(IC(gl)-IC(gl2))), max(abs(iid(gl)-iid(gl2))), max(abs(iid(gl,time=1)-iid(gl2,time=1))),
max(abs(IC(gl,time=1)-IC(gl2,time=1))),
max(abs(b1$base.iid- iid(gl,time=1))),
max(abs(b2$base.iid- iid(gl2,time=1))),
max(abs(b2$beta.iid- iid(gl2))),
max(abs(b1$beta.iid- iid(gl))))
expect_true( (max(diid)< 0.000000001))
###
d4 <- max(abs(t(head(cbind(b1$base.iid,b1$beta.iid))- head(cbind(b2$base.iid,b2$beta.iid)))))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
expect_true( (d3< 0.000000001) & ( d4 < 0.000000001 ))
expect_true( max(diid)< 0.000000001)
c(d1,d2,d3,d4,max(diid))

gl <- cifregFG(Event(time,status)~treatment+z+cluster(nn),data=hf,cause=1,cens.code=0)
gl2 <- cifregFG(Event(time,status)~treatment+z+cluster(nn),data=hf2,cause=1,cens.code=0)
d1 <- max(abs(estimate(gl)$coefmat - estimate(gl2)$coefmat )) 
###
gg <- survivalG(gl,hf,time=1)
gg2 <- survivalG(gl2,hf2,time=1)
d2 <- max(abs(gg$risk$coefmat- gg2$risk$coefmat ))
b1 <- iidBaseline(gl,time=1)
b2 <- iidBaseline(gl2,time=1)
diid <- c(max(abs(IC(gl)-IC(gl2))), max(abs(iid(gl)-iid(gl2))), max(abs(iid(gl,time=1)-iid(gl2,time=1))),
max(abs(IC(gl,time=1)-IC(gl2,time=1))),
max(abs(b1$base.iid- iid(gl,time=1))),
max(abs(b2$base.iid- iid(gl2,time=1))),
max(abs(b2$beta.iid- iid(gl2))),
max(abs(b1$beta.iid- iid(gl))))
expect_true( (max(diid)< 0.000000001))
d4 <- max(abs(cbind(b1$base.iid,b1$beta.iid)- cbind(b2$base.iid,b2$beta.iid)))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
expect_true( (d3< 0.000000001) & ( d4 < 0.000000001 ))
expect_true( max(diid)< 0.000000001)
c(d1,d2,d3,d4,max(diid))
max(abs(gl2$MGciid-gl$MGciid))



gl <- cifregFG(Event(time,status)~treatment+z+cluster(nn),data=hf,cause=1,cens.code=0)
gl2 <- cifregFG(Event(time,status)~treatment+z,data=hf,cause=1,cens.code=0)
d1 <- max(abs(estimate(gl)$coefmat - estimate(gl2)$coefmat )) 
###
gg <- survivalG(gl,hf,time=1)
gg2 <- survivalG(gl2,hf,time=1)
d2 <- max(abs(gg$risk$coefmat- gg2$risk$coefmat ))
b1 <- iidBaseline(gl,time=1)
b2 <- iidBaseline(gl2,time=1)
diid <- c(max(abs(IC(gl)-IC(gl2))), max(abs(iid(gl)-iid(gl2))), max(abs(iid(gl,time=1)-iid(gl2,time=1))),
max(abs(IC(gl,time=1)-IC(gl2,time=1))),
max(abs(b1$base.iid- iid(gl,time=1))),
max(abs(b2$base.iid- iid(gl2,time=1))),
max(abs(b2$beta.iid- iid(gl2))),
max(abs(b1$beta.iid- iid(gl))))
expect_true( (max(diid)< 0.000000001))
d4 <- max(abs(cbind(b1$base.iid,b1$beta.iid)- cbind(b2$base.iid,b2$beta.iid)))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
expect_true( (d3< 0.000000001) & ( d4 < 0.000000001 ))
expect_true( max(diid)< 0.000000001)
c(d1,d2,d3,d4,max(diid))
max(abs(gl2$MGciid-gl$MGciid))

############################ ok også her 

gl <- cifregFG(Event(time,status)~treatment+z+cluster(nn),data=hf2,cause=1,cens.code=0)
gl2 <- cifregFG(Event(time,status)~treatment+z,data=hf2,cause=1,cens.code=0)
d1 <- max(abs(estimate(gl)$coefmat - estimate(gl2)$coefmat )) 
gg <- survivalG(gl,hf2,time=1)
gg2 <- survivalG(gl2,hf2,time=1)
d2 <- max(abs(gg$risk$coefmat- gg2$risk$coefmat ))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))

## }}}
}
test_recregClust() 

test_binregClust <- function() {
## {{{  works binreg, binregG, binregATE 

gl <- binreg(Event(time,status)~treatment+z+cluster(id2),data=hf,cause=1,time=3)
gl2 <- binreg(Event(time,status)~treatment+z+cluster(id2),data=hf2,cause=1,time=3)
d1 <- max(abs(estimate(gl)$coefmat - estimate(gl2)$coefmat ))
d2 <- max(abs(t(head(iid(gl))- head(iid(gl2)))))
bbg<- binregG(gl,hf)
bbg2<- binregG(gl2,hf2)
head(iid(gl))
head(iid(gl))
head(bbg$risk.iid)
d3 <- max(abs(bbg$risk$coefmat- bbg2$risk$coefmat))
d4 <- max(abs(t((bbg$risk.iid)-(bbg2$risk.iid))))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
expect_true( (d3< 0.000000001) & ( d4 < 0.000000001 ))
c(d1,d2,d3,d4)


gl <- binreg(Event(time,status)~treatment+z+cluster(idn),data=hf,cause=1,time=3)
gl2 <- binreg(Event(time,status)~treatment+z+cluster(idn),data=hf2,cause=1,time=3)
d1 <- max(abs(estimate(gl)$coefmat - estimate(gl2)$coefmat ))
d2 <- max(abs(t((iid(gl))- (iid(gl2)))))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
bbg<- binregG(gl,hf)
bbg2<- binregG(gl2,hf2)
d3 <- max(abs(bbg$risk$coefmat- bbg2$risk$coefmat))
d4 <- max(abs((bbg$risk.iid)-(bbg2$risk.iid)))
expect_true( (d3< 0.000000001) & ( d4 < 0.000000001 ))
c(d1,d2,d3,d4)
###
###

gla <- binregATE(Event(time,status)~treatment+z+cluster(idn),data=hf,cause=1,time=3,
		treat.model=treatment~z)
gla2 <- binregATE(Event(time,status)~treatment+z+cluster(idn),data=hf2,cause=1,time=3,
		treat.model=treatment~z)
sgla2 <- summary(gla2)
sgla <- summary(gla)
gdif <- sgla$ateG- sgla2$ateG
d2 <- max(abs(gla$iid-gla2$iid))
d2
d3 <- max(abs(gdif))
d4 <- max(abs(sgla$ateDR-sgla2$ateDR))
d4
dd2 <- max(abs(gla2$riskG.iid-gla$riskG.iid))
dd3 <- max(abs(gla2$riskDR.iid-gla$riskDR.iid))
expect_true( (max(d1,d2,d3,d4)< 0.000000001) )
expect_true( (dd2< 0.000000001) & ( dd3 < 0.000000001 ))

gl <- binreg(Event(time,status)~treatment+z+cluster(nn),data=hf,cause=1,time=3)
gl2 <- binreg(Event(time,status)~treatment+z+cluster(nn),data=hf2,cause=1,time=3)
d1 <- max(abs(estimate(gl)$coefmat - estimate(gl2)$coefmat ))
d2 <- max(abs(iid(gl)- iid(gl2)))
bbg<- binregG(gl,hf)
bbg2<- binregG(gl2,hf2)
d3 <- max(abs(bbg$risk$coefmat- bbg2$risk$coefmat))
d4 <- max(abs(head(bbg$risk.iid)-head(bbg2$risk.iid)))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
expect_true( (d3< 0.000000001) & ( d4 < 0.000000001 ))
c(d1,d2,d3,d4)
###

###
gla <- binregATE(Event(time,status)~treatment+z+cluster(nn),data=hf,cause=1,time=3,
		treat.model=treatment~z)
gla2 <- binregATE(Event(time,status)~treatment+z+cluster(nn),data=hf2,cause=1,time=3,
		treat.model=treatment~z)
sgla2 <- summary(gla2)
sgla <- summary(gla)
g3 <- max(abs(sgla$ateG- sgla2$ateG))
g3
d2 <- max(abs(gla$iid-gla2$iid))
d4 <- max(abs(sgla$ateDR-sgla2$ateDR))
c(d2,d3,d4)
dd2 <- max(abs(gla2$riskG.iid-gla$riskG.iid))
dd3 <- max(abs(gla2$riskDR.iid-gla$riskDR.iid))
c(dd2,dd3)
expect_true( (max(d1,d2,d3,d4)< 0.000000001) )
expect_true( (dd2< 0.000000001) & ( dd3 < 0.000000001 ))

##################################################################################

gl <- binregATE(Event(time,status)~treatment+z+cluster(idn),data=hf,cause=1,time=3,
		treat.model=treatment~z)
sgl <- summary(gl)
###
gl2 <- binregATE(Event(time,status)~treatment+z+cluster(idn),data=hf2,cause=1,time=3,
		treat.model=treatment~z)
sgl2 <- summary(gl2)
d0 <- max(abs(sgl$coef - sgl2$coef ))
d1 <- max(abs(sgl$ateG-sgl2$ateG))
d2 <- max(abs(sgl$ateDR-sgl2$ateDR))
d3 <- max(abs(cbind(gl$riskG.iid,gl$riskDR.iid)- cbind(gl2$riskG.iid,gl2$riskDR.iid)))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
expect_true( (d3< 0.000000001) )
expect_true( (d0< 0.000000001) )
c(d0,d1,d2,d3)

## }}}
}
test_binregClust()

test_recgetIPCWClust <- function() {
## {{{ recregIPCW 

 ## IPCW at 2 years 
gl <- recregIPCW(Event(entry,time,status)~treatment+z+cluster(idn),data=hf,
     cause=1,death.code=2,time=2,cens.model=~strata(treatment))
gl2 <- recregIPCW(Event(entry,time,status)~treatment+z+cluster(idn),data=hf2,
     cause=1,death.code=2,time=2,cens.model=~strata(treatment))
d1 <- max(abs(estimate(gl)$coefmat - estimate(gl2)$coefmat ))
d2 <- max(abs(head(iid(gl))- head(iid(gl2))))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
c(d1,d2)

## }}} 
}
test_recgetIPCWClust()

test_binregRatioClust <- function() {
## {{{ also works  binregRatio, 

gl <- binregRatio(Event(time,status)~treatment+z+cluster(nn),data=hf,cause=1,time=3,type="II")
gl2 <- binregRatio(Event(time,status)~treatment+z+cluster(nn),data=hf2,cause=1,time=3,type="II")
d1 <- max(abs(estimate(gl)$coefmat - estimate(gl2)$coefmat ))
d2 <- max(abs(head(gl$iid)- head(gl2$iid)))
d3 <- max(abs(head(gl$iidI)- head(gl2$iidI)))
d4 <- max(abs(gl$coefI-gl2$coefI))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
expect_true( (d3< 0.000000001) & ( d4 < 0.000000001 ))
c(d1,d2,d3,d4)


gl <- binregRatio(Event(time,status)~treatment+z+cluster(id2),data=hf,cause=1,time=3)
gl2 <- binregRatio(Event(time,status)~treatment+z+cluster(id2),data=hf2,cause=1,time=3)
d1 <- max(abs(estimate(gl)$coefmat - estimate(gl2)$coefmat ))
d2 <- max(abs(head(gl$iid)- head(gl2$iid)))
###
d3 <- max(abs(head(gl$iidI)- head(gl2$iidI)))
d4 <- max(abs(gl$coefI-gl2$coefI))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
expect_true( (d3< 0.000000001) & ( d4 < 0.000000001 ))
c(d1,d2,d3,d4)

gl <- binregRatio(Event(time,status)~treatment+z+cluster(nn),data=hf2,cause=1,
		  time=3,type="II")
gl2 <- binregRatio(Event(time,status)~treatment+z,data=hf2,cause=1,time=3,type="II")
estimate(gl)$coefmat 
estimate(gl2)$coefma
d1
d1 <- max(abs(estimate(gl)$coefmat - estimate(gl2)$coefmat ))
d1
expect_true( (d1< 0.000000001) )

## }}}
}
test_binregRatioClust()

test_binregTSRClust <- function() {
## {{{ works binregTSR 

 library(mets)
 ddf <- mets:::gsim(200,covs=1,null=0,cens=1,ce=2)
 datat <- ddf$datat
 datat$id2 <- sample(200,200)[datat$id]
 datat2 <- dsort(datat,~id2+time)
###     
 bb <- binregTSR(Event(entry,time,status)~+1+cluster(id),datat,time=2,cause=c(1),
             cens.code=0,treat.model0=A0.f~+1,treat.model1=A1.f~A0.f,
             augmentR1=~X11+X12+TR,augmentR0=~X01+X02,
             augmentC=~A01+A02+X01+X02+A11t+A12t+X11+X12+TR,
             response.code=2)
###
  bb2 <- binregTSR(Event(entry,time,status)~+1+cluster(id),datat2,time=2,cause=c(1),
             cens.code=0,treat.model0=A0.f~+1,treat.model1=A1.f~A0.f,
             augmentR1=~X11+X12+TR,augmentR0=~X01+X02,
             augmentC=~A01+A02+X01+X02+A11t+A12t+X11+X12+TR,
             response.code=2)
  d1 <- 0
  for (i in 1:4) d1 <- max(d1,max(abs(bb2$riskG[[i]]-bb$riskG[[i]])))
  d2 <- 0
  for (i in 1:4) d2 <- max(d2,max(abs(bb2$riskG.iid[[i]]-bb$riskG.iid[[i]])))
  expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
 print( c(d1,d2))

 bb <- binregTSR(Event(entry,time,status)~+1+cluster(id2),datat,time=2,cause=c(1),
             cens.code=0,treat.model0=A0.f~+1,treat.model1=A1.f~A0.f,
             augmentR1=~X11+X12+TR,augmentR0=~X01+X02,
             augmentC=~A01+A02+X01+X02+A11t+A12t+X11+X12+TR,
             response.code=2)
###
  bb2 <- binregTSR(Event(entry,time,status)~+1+cluster(id2),datat2,time=2,cause=c(1),
             cens.code=0,treat.model0=A0.f~+1,treat.model1=A1.f~A0.f,
             augmentR1=~X11+X12+TR,augmentR0=~X01+X02,
             augmentC=~A01+A02+X01+X02+A11t+A12t+X11+X12+TR,
             response.code=2)
  d1 <- 0
  for (i in 1:4) d1 <- max(d1,max(abs(bb2$riskG[[i]]-bb$riskG[[i]])))
  d2 <- 0
  for (i in 1:4) d2 <- max(d2,max(abs(bb2$riskG.iid[[i]]-bb$riskG.iid[[i]])))
  expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
  print(c(d1,d2))

  ## }}}
}
test_binregTSRClust()

test_phregIPTWClust <- function() {
## {{{ phreg_IPTW

gl <- phreg_IPTW(Event(entry,time,status==1)~treatment+cluster(idn),data=hf, treat.model=treatment~z)
gl2 <-phreg_IPTW(Event(entry,time,status==1)~treatment+cluster(idn),data=hf2, treat.model=treatment~z)
d1 <- max(abs(estimate(gl)$coefmat - estimate(gl2)$coefmat ))
d2 <- max(abs(t(head(iid(gl2))- head(iid(gl)))))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
c(d1,d2)

library(mets)
data <- mets:::simLT(0.7,200,beta=0.3,betac=0,ce=1,betao=0.3)
dfactor(data) <- Z.f~Z
data$id2 <- sample(200,200)
data2 <- data[order(data$id2),]
out <- phreg_IPTW(Surv(time,status)~Z.f+cluster(id2),data=data,treat.model=Z.f~X)
out2 <- phreg_IPTW(Surv(time,status)~Z.f+cluster(id2),data=data2,treat.model=Z.f~X)
d1 <- max(abs(iid(out)-iid(out2)))
d2 <- max(abs(estimate(out)$coefmat-estimate(out2)$coefmat))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))

data(calgb8923)
calgt <- calgb8923
calgt$treatvar <- 1
calgt$id2 <- sample(388,388)[calgt$id]
calgt2 <- dsort(calgt,~id2+time)
###
tm=At.f~factor(Count2)
tm=At.f~factor(Count2)*A0.f
tm=At.f~factor(Count2)+age+sex+wbc
###
gl <- phreg_IPTW(Event(start,time,status==1)~A0+A10+cluster(id2),calgt,treat.model=tm,treat.var="treatvar")
gl2 <- phreg_IPTW(Event(start,time,status==1)~A0+A10+cluster(id2),calgt2,treat.model=tm,treat.var="treatvar")
d1 <- max(abs(estimate(gl2)$coefmat- estimate(gl)$coefmat))
d2 <- max(abs(gl$IID-gl2$IID))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))

## }}}
}
test_phregIPTWClust()

test_phregRCT <- function() {
## {{{ phreg_rct


gl <- phreg_rct(Event(entry,time,status==1)~treatment+cluster(idn),data=hf,augmentR0=~z,augmentC=~z+Count1)
###
gl2 <-phreg_rct(Event(entry,time,status==1)~treatment+cluster(idn),data=hf2,augmentR0=~z,augmentC=~z+Count1)
d1 <- max(abs(gl$coefs-gl2$coefs))
d2 <- max(abs(head(gl$iid[[1]])- head(gl2$iid[[1]])))
d3 <- max(abs(head(gl$iid[[2]])- head(gl2$iid[[2]])))
expect_true( (d1< 0.000000001) )
expect_true( (d2< 0.000000001) & ( d3 < 0.000000001 ))
c(d1,d2,d3)

gl <- phreg_rct(Event(entry,time,status==1)~treatment+cluster(id2),data=hf,augmentR0=~z,augmentC=~z+Count1,return.augmentR0=TRUE)
###
gl2 <-phreg_rct(Event(entry,time,status==1)~treatment+cluster(id2),data=hf2, augmentR0=~z,augmentC=~z+Count1,return.augmentR0=TRUE)
d1 <- max(abs(gl$coefs-gl2$coefs))
###gl$coefs- gl2$coefs
d2 <- max(abs(head(gl$iid[[1]])- head(gl2$iid[[1]])))
d3 <- max(abs(head(gl$iid[[2]])- head(gl2$iid[[2]])))
expect_true( (d1< 0.000000001) )
expect_true( (d2< 0.000000001) & ( d3 < 0.000000001 ))
c(d1,d2,d3)

gl$AugR0-gl2$AugR0
gl$AugCdyn-gl2$AugCdyn

###gl$data.augR0$ea- gl2$data.augR0$ea
###gl$data.augR0$XR- gl2$data.augR0$XR
###gl$data.augR0$XRpi- gl2$data.augR0$XRpi
###gl$data.augR0$p0- gl2$data.augR0$p0
######
###names(gl$data.augR0)

gl <- phreg_rct(Event(entry,time,status==1)~treatment+cluster(id2),data=hf,augmentR0=~z,augmentC=~z+Count1,return.augmentR0=TRUE,estpr=0)
###
gl2 <-phreg_rct(Event(entry,time,status==1)~treatment+cluster(id2),data=hf2, augmentR0=~z,augmentC=~z+Count1,return.augmentR0=TRUE,estpr=0)
d1 <- max(abs(gl$coefs-gl2$coefs))
d2 <- max(abs(head(gl$iid[[1]])- head(gl2$iid[[1]])))
d3 <- max(abs(head(gl$iid[[2]])- head(gl2$iid[[2]])))
expect_true( (d1< 0.000000001) )
expect_true( (d2< 0.000000001) & ( d3 < 0.000000001 ))
c(d1,d2,d3)


data(calgb8923)
calgt <- calgb8923
calgt$treatvar <- 1
set.seed(100)
calgt$id2 <- sample(388,388)[calgt$id]
calgt2 <- dsort(calgt,~id2+time)
###
tm=At.f~factor(Count2)
tm=At.f~factor(Count2)*A0.f
tm=At.f~factor(Count2)+age+sex+wbc
###
gl <- phreg_rct(Event(start,time,status==1)~A0+A10+cluster(id2),calgt,RCT=TRUE,treat.model=tm,treat.var="treatvar",
		augmentR0=~sex+age+wbc,augmentR1=~sex+wbc+TR,augmentC=~sex+age+wbc)
gl2 <- phreg_rct(Event(start,time,status==1)~A0+A10+cluster(id2),calgt2,RCT=TRUE,treat.model=tm,treat.var="treatvar",
		augmentR0=~sex+age+wbc,augmentR1=~sex+wbc+TR,augmentC=~sex+age+wbc)
d1 <- max(abs(gl2$coefs- gl$coefs))
d2 <- 0
for (i in 1:length(gl$iid)) d2 <- max(d2,max(abs(head(gl$iid[[i]])- head(gl2$iid[[i]]))))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
c(d1,d2)
d3 <- max( gl$AugR0-gl2$AugR0, gl$AugR1-gl2$AugR1, gl$AugR01-gl2$AugR01, gl$AugCdyn-gl2$AugCdyn)
expect_true( (d3< 0.000000001) )
###
c(gl$AugR0-gl2$AugR0, gl$AugR1-gl2$AugR1, gl$AugR01-gl2$AugR01, gl$AugCdyn-gl2$AugCdyn)


gl <- phreg_rct(Event(start,time,status==1)~A0+A10+cluster(id2),calgt,RCT=FALSE,treat.model=tm,treat.var="treatvar",
		augmentR0=~sex+age+wbc,augmentR1=~sex+wbc+TR,augmentC=~sex+age+wbc)
gl2 <- phreg_rct(Event(start,time,status==1)~A0+A10+cluster(id2),calgt2,RCT=FALSE,treat.model=tm,treat.var="treatvar",
		augmentR0=~sex+age+wbc,augmentR1=~sex+wbc+TR,augmentC=~sex+age+wbc)
d1 <- max(abs(gl2$coefs- gl$coefs))
d2 <- 0
for (i in 1:length(gl$iid)) d2 <- max(d2,max(abs(head(gl$iid[[i]])- head(gl2$iid[[i]]))))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
c(d1,d2)
d3 <- max( gl$AugR0-gl2$AugR0,gl$AugR1-gl2$AugR1,gl$AugR01-gl2$AugR01,gl$AugCdyn-gl2$AugCdyn)
expect_true( (d3< 0.000000001) )


## }}}
}
test_phregRCT()

test_WAClust <- function() {
## {{{ WA_recurrent

w1 <- WA_recurrent(Event(entry,time,status)~treatment+cluster(idn),data=hf,time=2)
w2 <- WA_recurrent(Event(entry,time,status)~treatment+cluster(idn),data=hf2,time=2)
s1 <- summary(w1)
s2 <- summary(w2)
s1
s2
d1 <- max(abs(s2$rmst$coefmat- s1$rmst$coefmat))
d2 <- max(abs(s2$meanNtD$coefmat- s2$meanNtD$coefmat))
d3 <- max(abs(s2$ratio$coefmat - s1$ratio$coefmat))
d4 <- max(abs(s2$meanpt$coefmat- s1$meanpt$coefmat))
c(d1,d2,d3,d4)
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
expect_true( (d3< 0.000000001) & ( d4 < 0.000000001 ))
c(d1,d2,d3,d4)
###
d1 <- max(abs(w1$ET$riskDR$riskG.iid- w2$ET$riskDR$riskG.iid))
d2 <- max(abs(w1$ET$riskDR$riskDR.iid- w2$ET$riskDR$riskDR.iid))
d3 <- max(abs(w1$ET$riskDR$iid -w2$ET$riskDR$iid))
head(w1$ET$riskDR$iid) 
head(w2$ET$riskDR$iid)
d4 <- max(abs((w1$RAW$iid- w2$RAW$iid)))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
expect_true( (d3< 0.000000001) & ( d4 < 0.000000001 ))
c(d1,d2,d3,d4)


w1 <- WA_recurrent(Event(entry,time,status)~treatment+cluster(idn),data=hf,time=2,augmentC=~z+Count1,augmentR=~z)
w2 <- WA_recurrent(Event(entry,time,status)~treatment+cluster(idn),data=hf2,time=2,augmentC=~z+Count1,augmentR=~z)
###
s1 <- summary(w1)
s2 <- summary(w2)
d1 <- max(abs(s2$rmst$coefmat- s1$rmst$coefmat))
d3 <- max(abs(s2$meanNtD$coefmat- s2$meanNtD$coefmat))
d4 <- max(abs(s2$ratio$coefmat - s1$ratio$coefmat))
d5 <- max(abs(s2$meanpt$coefmat- s1$meanpt$coefmat))
c(d1,d2,d3,d4,d5)
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
expect_true( (d3< 0.000000001) & ( d4 < 0.000000001 ))
expect_true( (d5< 0.000000001) )
c(d1,d2,d3,d4,d5)
###
d4 <- max(abs((w1$RAW$iid- w2$RAW$iid)))
d1 <- max(abs(head(w1$ET$riskDR$riskG.iid)- head(w2$ET$riskDR$riskG.iid)))
d2 <- max(abs(head(w1$ET$riskDR$riskDR.iid)- head(w2$ET$riskDR$riskDR.iid)))
d3 <- max(abs(head(w1$ET$riskDR$iid) - head(w2$ET$riskDR$iid)))
d5 <- max(abs(head(w1$ET$riskDRC$iid) - head(w2$ET$riskDRC$iid)))
expect_true( (d1< 0.000000001) & ( d2 < 0.000000001 ))
expect_true( (d3< 0.000000001) & ( d4 < 0.000000001 ))
expect_true( (d5< 0.000000001) )
c(d1,d2,d3,d4,d5)

## }}}
}
test_WAClust()

