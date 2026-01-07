FROM rocker/r-devel-ubsan-clang

RUN RD -e \
	'.libPaths("/usr/lib/R/library"); install.packages(c("devtools", "codetools", "timereg", "prodlim", "cmprsk", "testthat", "icenReg", "optimx", "prodlim", "ucminf", "knitr", "riskRegression", "bookdown", "rmarkdown", "ggplot2", "cowplot", "mvtnorm", "numDeriv", "lava", "Rcpp", "RcppArmadillo"))'

WORKDIR /mets

CMD cd /; RD CMD check --library=/usr/lib/R/library --as-cran mets
