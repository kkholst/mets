
##' @export
loglikMVN <- function(yl, yu, status, mu, S, thres) {
  .loglikMVN(
          Yl=as.matrix(yl),
          Yu=as.matrix(yu),
          Status=as.integer(status),
          Mu=as.matrix(mu),
          S=as.matrix(S),
          z=NULL,su=NULL,
          Threshold=as.matrix(thres),
          itol=lava::lava.options()$itol
  )
}

##' @export
scoreMVN <- function(y, m, s, dm, ds) {
    .scoreMVN(
          Y=as.matrix(y),
          Mu=as.matrix(m),
          dMu=as.matrix(dm),
          S=as.matrix(s),
          dS=as.matrix(ds),
          itol=lava::lava.options()$itol
      )
}
