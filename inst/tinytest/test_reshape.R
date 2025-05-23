## Reshaping data
library("tinytest")

test_fastreshape1 <- function() {
  m <- lava::lvm()
  lava::regression(m,c(y1,y2,y3)~x) <- c(1,10,100)
  lava::distribution(m,~x) <- f <- function(n,...) rbinom(n,1,0.5)+1
  d <- lava::sim(m,10);
  dd <- fast.reshape(d,var="y")
  d1 <- fast.reshape(dd,id="id")
  expect_true(sum((d[, lava::endogenous(m)]-
                   d1[,lava::endogenous(m)])^2)<1e-20)
  d2 <- fast.reshape(dd,id="id",var="y",num="num")
  expect_true(sum((d-d2[,colnames(d)])^2)<1e-20)
}
test_fastreshape1()

test_fastreshape2 <- function() {
  testdata <- data.frame(hour=c(12,13,14,11,12,14,15,16),id=c(1,1,1,2,2,3,3,3),y=round(rnorm(8),2))
  widetest <- reshape(testdata,v.names="y",idvar="id",direction="wide",timevar="hour")
  wide <- fast.reshape(testdata,varying="y",id="id",num="hour",sep=".")
  expect_equivalent(widetest,wide[,colnames(widetest)])
}
test_fastreshape2()

test_fastrehape_different_data_types <- function() {
  d <- data.frame(time1=c(1:5),
                  time2=c(6.070311,
                          2.026996,
                          7.584480,
                          8.630120,
                          8.193392))
  dd <- mets::fast.reshape(d)
  expect_equivalent(dd[,1],as.vector(t(d)))

  d <- data.frame(time1=c(TRUE,FALSE,TRUE,FALSE,TRUE),
                  time2=c(6.070311,
                          2.026996,
                          7.584480,
                          8.630120,
                          8.193392))
  dd <- mets::fast.reshape(d)
  expect_equivalent(dd[,1],as.vector(t(d)))
}
test_fastrehape_different_data_types()


## fast.reshape(fast.reshape(d,var=c("y","z","w")),id="id",var=c("y","z","w"))
## library(mets)
## x <- matrix(1:10,5,2)
## x[3,2] <- 8
## x[3,2] <- NA
## cluster <- c(1,1,2,2,3)
## x <- cbind(x,cluster)
## x
## ud <- fast.reshape(data.frame(x),"cluster")
## ud
## ###
## out=cluster.index(cluster)
## out
## ###
## ud <- faster.reshape(x,cluster)
## ud
## ud <- faster.reshape(data.frame(x),cluster)
## ud
## ###
## colnames(x) <- c("y1","y2","cluster")
## x
## ud <- fast.reshape(data.frame(x),"cluster")
## ud
## ###
## num <- c(2,1,1,2,2)
## x
## out <- faster.reshape(x,cluster)
## out
## out <- faster.reshape(x,cluster,num=num)
## out

