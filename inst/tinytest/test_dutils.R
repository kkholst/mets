# dutils
library("tinytest")

test_dsort <- function() {
  data("hubble",package="lava")
  expect_equivalent(
              order(dsort(hubble, ~sigma)$sigma),
              seq_len(nrow(hubble)))
}
test_dsort()

test_dsort2 <- function() {
  data("hubble",package="lava")
  h1 <- hubble
  drename(h1,fun=toupper) <- ~.
  expect_equivalent(colnames(h1),
                              toupper(names(hubble)))
}
test_dsort2()

test_daggregate <- function() {
  dd <- data.frame(a = 1:20, b = 20:1,
                   g1 = rep(0:1, 10),
                   g2 = rep(0:1, each = 10),
                   g3 = rbinom(20, 1, 0.5))
  dd$g1[1:2] <- NA
  dd$g2[2:3] <- NA
  dd$a[3:7] <- NA
  ##dcor(dd,.~g1+g2)

  r1 <- dcor(dd, "[acg][12]*$", regex=TRUE, use="pairwise")
  r2 <- cor(dd[,c("a","g1","g2")], use="pairwise")
  expect_equal(norm(r1-r2), 0)
  r1 <- dcor(dd, "[ab]",subset=is.na(g1), regex=TRUE, use="pairwise")
  dd0 <- subset(dd, is.na(g1))
  r2 <- cor(dd[, c("a","b")], use="pairwise")
  expect_equal(norm(r1-r2), 0)
  suppressWarnings(r1 <- dcor(dd, ~.|is.na(g1)|is.na(g2))) # do not warn about zero sd
  dd0 <- subset(dd, is.na(g1) | is.na(g2), select=-c(g1,g2))
  suppressWarnings(r2 <- dcor(dd0)) # do not warn about zero sd
  expect_equal(r1,r2)
  r1 <- dcor(dd, use="pairwise")
  expect_equal(norm(r1-cor(dd, use="pairwise")), 0)
}
test_daggregate()


test_dby <- function() {
  dd <- dby2(iris, . ~ Species, mean, median, REDUCE=T)
  val <- dreshape(dd,varying=list(mean="mean*",median="median*"),dropid=TRUE)
  val$num <- gsub("mean.","",val$num)

  val <- dreshape(dd,varying=c("mean","median"),dropid=TRUE)
  val$num <- gsub("^.","",val$num)
    
  expect_true(ncol(val)==4)
  expect_true(nrow(val)==3*4)
  val0 <- subset(val, Species=="setosa" & num=="Sepal.Width")
  val1 <- subset(iris, Species=="setosa", select="Sepal.Width")[,1]
  expect_true(abs(val0[1,"mean"]-mean(val1))<1e-16)
  expect_true(abs(val0[1,"median"]-median(val1))<1e-16)
}
test_dby()

test_dsample <- function() {
  n <- 10
  d <- data.frame(id=rep(0:1,5),x=seq(0,1,length.out=10),y=seq(1,0,length.out=10),
                  id2=rep(1:5,each=2))
  d1 <- dsample(d, ~id, size=3)
  expect_true(ncol(d1)==5)
  expect_true(nrow(d1)==15)

  d2 <- dsample(d, .~id|x>0.5, size=3)
  expect_true(ncol(d2)==4)
  expect_true(all(d2$x>0))

  d3 <- dsample(d, .~id+id2, size=3)
  expect_true(ncol(d3)==3)

  d4 <- dsample(d, .~id+id2|x>0, size=3)
  expect_true(ncol(d4)==3)
  expect_true(all(d4$x>0))
}
test_dsample()
