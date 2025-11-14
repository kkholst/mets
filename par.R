library("targeted")
## devtools::load_all("~/Software/targeted")
devtools::load_all("~/Software/mets")

m <- lava::lvm(stop ~ x) |>
  lava::distribution(~ start, value = lava::coxWeibull.lvm(shape=3,scale=5)) |>
  lava::distribution(~ stop, value = lava::coxWeibull.lvm(shape=3,scale=5)) |>  
  transform(~status, value = \(...) TRUE) |>
  distribution(~id, value = lava::Sequence.lvm(integer = TRUE))
d0 <- lava::sim(m, 2e4, p = c("stop~x" = 1))
d <- subset(d0, start < stop)

des <- targeted::design(Event(start, stop, status) ~ x + strata(x>0),
                        data=d, specials=c("strata"))




library(eha)
phreg.par(Surv(stop, status) ~ x, data=d0)
phreg.par(Surv(stop, status) ~ x, data=d)

mets::phreg(Surv(stop, status) ~ x, data=d)
mets::phreg(Surv(start, stop, status) ~ x, data=d)


phreg(Surv(start, stop, status) ~ x, data=d)

eha::phreg(Surv(stop, status)~1,data=d)

weibreg(Surv(stop,status)~x,data=d)

 ## Note in simulation A(t) = lambda*t^scale
 ## but here A(t) (scale*t)^shape, hence
 ## lambda := scale^(1/shape)
 tt <- seq(0,100,length.out=100)
 plot(tt,exp(-(exp(-4)*tt)^exp(.5))
 y <- runif(2e4,0,100)
 (op <- phreg.par(y,TRUE))

 tt <- seq(0,100,length.out=100)
 cc <- survival::coxph(Surv(y,rep(TRUE,length(y)))~1)
 plot(survfit(cc),mark.time=FALSE)
 lines(tt,exp(-(exp(-4)*tt)^exp(.484)),col="red")
(op <- phreg.weibull(d$y,TRUE,cbind(d$x)))
(a <- survival::survreg(Surv(y,status)~1+x,dist="weibull",data=d))
(e <- eha::weibreg(Surv(y,status)~x,data=d))
