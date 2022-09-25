FROM rocker/r-devel-ubsan-clang

RUN RD -e \
	'.libPaths("/usr/lib/R/library"); install.packages(c("devtools", "timereg", "prodlim", "cmprsk", "testthat", "ucminf", "knitr", "bookdown", "rmarkdown", "ggplot2", "cowplot"))'

RUN RD -e \
	'.libPaths("/usr/lib/R/library"); install.packages(c("mvtnorm", "numDeriv", "lava", "Rcpp", "RcppArmadillo"))'

WORKDIR /mets

CMD cd /; RD CMD check --library=/usr/lib/R/library --as-cran mets
