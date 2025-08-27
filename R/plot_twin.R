##' Scatterplot with contours of the (kernel) estimated density
##'
##' @title Scatter plot function
##' @param data bivariate data to plot (data.frame or matrix with 2 columns)
##' @param density.args argumnets to (1d) density estimator
##' @param kernsmooth.args arguments to 2d-kernel smoother
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param col color of points
##' @param col2 color of contour / density plot
##' @param alpha transparency level of points
##' @param grid should grid be added to the plot
##' @author Klaus Kähler Holst
##' @export
plot_twin <- function (data,
                       density.args = list(),
                       kernsmooth.args = list(),
                       xlab, ylab,
                       col = "black", col2 = "lightblue",
                       alpha = 0.3,
                       grid = TRUE
                       ) {

  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  if (!ncol(data) > 1) stop("data much contain at least 2 columns")
  data <- na.omit(data[, 1:2])
  xlab <- colnames(data)[1]
  ylab <- colnames(data)[2]
  if (is.null(kernsmooth.args$bandwidth)) { # bandwidth
    kernsmooth.args$bandwidth <- as.numeric(apply(data, 2, sd)) *
      nrow(data)^(-1 / 5)
  }
  x <- data[, 1, drop = TRUE]
  y <- data[, 2, drop = TRUE]
  d1 <- do.call(density, c(list(x), density.args))
  d2 <- do.call(density, c(list(y), density.args))
  rg <- range(x, y)
  top <- max(c(d1$y, d2$y)) # nolint
  layout(
    matrix(c(2, 0, 1, 3), 2, 2, byrow = TRUE),
    c(3, 1), c(1, 3), TRUE
  )
  par(mar = c(3, 3, 1, 1))
  plot(x, y,
       xlim = rg, ylim = rg,
       xlab = "", ylab = "",
       pch = 16, col = lava::Col(col, alpha)
       )
  graphics::grid()
  est <- do.call(KernSmooth::bkde2D, c(list(data), kernsmooth.args))
  with(est, contour(x1, x2, fhat, add=TRUE, col = col2))
  with(d1, plot(x, y, type = "l", lwd=0.5, axes = FALSE, ylim=c(0, top)))
  x <- with(d1, c(x, rev(x)))
  y <- with(d1, c(y, rep(0, length(y))))
  polygon(x, y, lty = 0, col = col2)
  axis(2)
  mtext(xlab, 1, line = 2)
  with(d2, plot(y, x, type = "l", lwd=0.5, axes = FALSE, xlim=c(0, top)))
  x <- with(d2, c(x, rev(x)))
  y <- with(d2, c(y, rep(0, length(y))))
  polygon(y, x, lty = 0, col = col2)
  axis(1)
  mtext(ylab, side = 2, line = 2, srt = 20)
  par(def.par)
}
