##' Scatterplot with contours of the (kernel) estimated density
##'
##' @title Scatter plot function
##' @param data bivariate data to plot (data.frame or matrix with 2 columns)
##' @param density.args argumnets to marginal estimator (`density`
##'  continuous data, `barplot` for categorical )
##' @param kernsmooth.args arguments to 2d-kernel smoother
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param col color of points
##' @param col2 color of contour / density plot
##' @param alpha transparency level of points
##' @param grid should grid be added to the plot
##' @param ... arguments to lower level plot functions
##' @author Klaus Kähler Holst
##' @inheritParams graphics::mosaicplot
##' @export
##' @examples
##' data("twinbmi", package="mets")
##' twinwide <- fast.reshape(twinbmi, id="tvparnr",varying=c("bmi"))
##' datamz <- log(subset(twinwide, zyg=="MZ")[,c("bmi1","bmi2")])
##'
##' # continuous variables
##' plot_twin(datamz)
##'
##' # categorical variables
##' datamz2 <- datamz
##' datamz2[, 1] <- cut(datamz[, 1], 4)
##' datamz2[, 2] <- cut(datamz[, 2], 4)
##' plot_twin(datamz2, color = TRUE)
##'
##' # survival variables
##' cens1 <- rbinom(nrow(datamz), 1, 0.5)
##' cens2 <- rbinom(nrow(datamz), 1, 0.5)
##' datamz2[, 1] <- Event(datamz[, 1], cens1)
##' datamz2[, 2] <- suppressWarnings(Event(datamz[, 2], cens2))
##' plot_twin(datamz2)
##'
##' rm(datamz, datamz2, cens1, cens2)
plot_twin <- function (data,
                       marginal.args = list(),
                       kernsmooth.args = list(),
                       xlab, ylab,
                       col = "black", col2 = "lightblue",
                       alpha = 0.3,
                       grid = TRUE,
                       side.plot = TRUE,
                       ...
                       ) {

  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  if (!ncol(data) > 1) stop("data much contain at least 2 columns")
  data <- na.omit(data[, 1:2])
  if (missing(xlab)) xlab <- colnames(data)[1]
  if (missing(ylab)) ylab <- colnames(data)[2]
  x <- data[, 1, drop = TRUE]
  y <- data[, 2, drop = TRUE]
  if (class(x)[1] != class(y)[1]) {
      stop("Expecting tabular data with 2 columns of the same type")
  }
  if (side.plot) {
      layout(
          matrix(c(2, 4, 1, 3), 2, 2, byrow = TRUE),
          c(3, 1), c(1, 3), TRUE
      )
  }
  par(mar = c(4, 4, 1, 1))
  if (inherits(x, c("character", "factor"))) {
    tt <- prop.table(table(data))
    tt <- tt[, rev(seq_len(ncol(tt)))]
    mosaicplot(tt,
        main = "",
        xlab = xlab,
        ylab = ylab,
        ...
        )
    if (side.plot) {
        do.call(barplot, c(list(table(data[, 1])), marginal.args))
        do.call(barplot, c(list(table(data[, 2]), horiz = TRUE), marginal.args))
    }
    par(def.par)
    return(invisible())
  }
  if (inherits(x, c("Surv", "Event"))) {
    x0 <- x[, 1]
    y0 <- y[, 1]
    rg <- range(x0, y0)
    plot(x0, y0,
        xlim = rg, ylim = rg,
        xlab = "", ylab = "",
        type = "n",
        ...
        )
    if (grid) graphics::grid()
    cens1 <- !x[, 2]
    cens2 <- !y[, 2]
    if (any(cens1)) {
      idx <- which(cens1)
      points(x0[idx], y0[idx],
          pch = 3,
          col = lava::Col(col, alpha)
          )
    }
    if (any(cens2)) {
      idx <- which(cens2)
      points(x0[idx], y0[idx],
          pch = 5,
          col = lava::Col(col, alpha)
          )
    }
    if (any((!cens1) & (!cens2))) {
        idx <- which(!cens1 & !cens2)
        points(x0[idx], y0[idx],
            pch = 16, col = lava::Col(col, alpha)
        )
    }
    if (!side.plot) {
      par(def.par)
      return(invisible())
    }
    p1 <- phreg(t ~ 1, data.frame(t = x))
    p2 <- phreg(t ~ 1, data.frame(t = y))
    chaz1 <- basecumhaz(p1)[[1]]$cumhaz
    chaz2 <- basecumhaz(p2)[[1]]$cumhaz
    chaz1[, 2] <- exp(-chaz1[, 2])
    chaz2[, 2] <- exp(-chaz2[, 2])
    top <- max(chaz2[,2], chaz1[,2])
    ## par(mar = c(3, 3, 1, 1))
    plot(chaz1[, 1], chaz1[, 2],
        axes = FALSE,
        xlim = rg, ylim = c(0, top),
        type = "s",
        xlab = xlab,
        ylab = "Surv. prob."
        )
    if (grid)  grid()
    axis(2)
    axis(1)
    ## par(mar = c(3, 3, 1, 1))
    plot(chaz2[, 2], chaz2[, 1],
        axes = FALSE,
        ylim = rg, , xlim = c(0, top),
        type = "s",
        ylab = ylab,
        xlab = "Surv. prob."
        )
    if (grid)  grid()
    axis(1)
    axis(2)
    plot(1, type = "n", axes = FALSE, xlab="", ylab="")
    legend("left",
           c("no cens", "1. cens", "2. cens"),
           pch = c(16, 3, 5),
           bty = "n"
    )

    par(def.par)
    return(invisible())
  }

  if (is.null(kernsmooth.args$bandwidth)) { # bandwidth
      kernsmooth.args$bandwidth <- as.numeric(apply(data, 2, sd)) *
          nrow(data)^(-1 / 5)
  }
  d1 <- do.call(density, c(list(x), marginal.args))
  d2 <- do.call(density, c(list(y), marginal.args))
  rg <- range(x, y)
  top <- max(c(d1$y, d2$y)) # nolint
  plot(x, y,
       xlim = rg, ylim = rg,
       xlab = "", ylab = "",
       pch = 16, col = lava::Col(col, alpha),
       ...
       )
  graphics::grid()
  est <- do.call(KernSmooth::bkde2D, c(list(data), kernsmooth.args))
  with(est, contour(x1, x2, fhat, add = TRUE, col = col2))
  if (side.plot) {
    with(d1, plot(x, y,
      xlab = "", ylab = "",
      type = "l", lwd = 0.5,
      axes = FALSE, ylim = c(0, top)
    ))
    x <- with(d1, c(x, rev(x)))
    y <- with(d1, c(y, rep(0, length(y))))
    polygon(x, y, lty = 0, col = col2)
    axis(2)
    mtext(xlab, 1, line = 2)
    with(d2, plot(y, x,
       xlab = "", ylab = "",
       type = "l", lwd = 0.5,
       axes = FALSE, xlim = c(0, top)
     ))
    x <- with(d2, c(x, rev(x)))
    y <- with(d2, c(y, rep(0, length(y))))
    polygon(y, x, lty = 0, col = col2)
    axis(1)
    mtext(ylab, side = 2, line = 2, srt = 20)
  }
  par(def.par)
}
