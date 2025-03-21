#' @rdname mex
#' @export
ggplot.mex <- function(data = NULL, mapping,
                       ptcol = "blue",
                       col = "cornflowerblue",
                       fill = "orange",
                       plot. = TRUE,
                       quantiles = seq(0.1, by = 0.2, len = 5),
                       ..., environment) {
  mar <- data[[1]]
  dep <- data[[2]]
  z <- as.matrix(dep$Z[order(mar$transformed[mar$transformed[, dep$which] > dep$dth, dep$which]), ]) # order Z since this is plotted agains p which is an ordered sequence
  colnames(z) <- colnames(dep$Z)
  n <- nrow(z)

  # Determine x-axis values to estimate CE quantiles at along conditioned var
  xmax <- max(mar$data[, dep$which]) # max for conditioned variable
  # marginal parameters
  sig <- coef(mar)[3, dep$which]
  xi <- coef(mar)[4, dep$which]
  marThr <- mar$mth[dep$which] # mth
  marP <- mar$mqu[dep$which] # mqu
  if (xi < 0) upperEnd <- marThr - sig / xi # upper end of domain
  len <- 501
  depThr <- c(quantile(mar$data[, dep$which], dep$dqu)) # dth
  dif <- xmax - depThr
  xlim <- unname(c(depThr - 0.1 * dif, depThr + 1.5 * dif))

  # Determine upper limit of x-axis
  if (xi < 0 && xlim[2] > upperEnd) {
    xlim[2] <- upperEnd
    plim <- 1
  } else {
    plim <- pgpd(xlim[2], sigma = sig, xi = xi, u = marThr)
  }
  # CDF values to plot at
  p <- seq(dep$dqu, 1 - 1 / n, length = n)
  plotp <- seq(dep$dqu, plim, len = len)[-len] # take out largest point to avoid Inf in p2q transform
  # transform to original scale from; these will be x-values in plot
  plotx <- revTransform(
    plotp,
    data = mar$data[, dep$which],
    qu = dep$dqu,
    th = depThr,
    sigma = sig,
    xi = xi
  )
  # convert probs to Laplace scale (these are values to calculate CE line at)
  xq <- dep$margins$p2q(plotp)

  plotfn <- function(i) {
    plots <- vector("list", length = 3)
    dat <- data.frame(p = p, z = z[, i], absz = abs(z[, i] - mean(z[, i])))
    # Plot residuals vs F(p)
    plots[[1]] <- ggplot(dat, aes(p, z)) +
      geom_point(colour = col, alpha = 0.5) +
      labs(
        x = paste("F(", dep$conditioningVariable, ")", sep = ""),
        y = paste("Z   ", colnames(z)[i], " | ", dep$conditioningVariable, sep = "")
      ) +
      geom_smooth(col = ptcol, fill = fill)

    # Plot absolute residuals vs F(p)
    plots[[2]] <- ggplot(dat, aes(p, absz)) +
      geom_point(colour = col, alpha = 0.5) +
      labs(
        x = paste("F(", dep$conditioningVariable, ")", sep = ""),
        y = paste("|Z - mean(Z)|   ", colnames(z)[i], " | ", dep$conditioningVariable, sep = "")
      ) +
      geom_smooth(col = ptcol, fill = fill)

    # conditioned vs conditioning variable to show CE against (backdrop)
    dat <- data.frame(
      x = mar$data[, dep$which],
      y = as.matrix(mar$data[, -dep$which])[, i]
    )

    # pull dependence coefficients and quantiles of residuals
    co <- coef(dep)[, i]
    zq <- quantile(dep$Z[, i], quantiles)
    # calculates regression lines from quantiles of residuals
    yq <- sapply(zq, function(z, co, xq) {
      co["a"] * xq + co["c"] - co["d"] * log(xq) + xq^co["b"] * z
    }, xq, co = co)

    # add to plot
    plots[[3]] <- ggplot(dat, aes(x, y)) +
      geom_point(col = col, alpha = 0.5) +
      labs(x = dep$conditioningVariable, y = colnames(z)[i]) +
      geom_vline(xintercept = depThr) # add dependence threshold

    # add quantile contour lines
    if (is.null(mar$referenceMargin)) {
      trns <- dat$y # conditioning variable on original scale
      # marginal parameters for this variable
      qu <- mar$mqu[-dep$which][i]
      th <- mar$mth[-dep$which][i]
      sigma <- coef(mar)[3, -dep$which][i]
      xi <- coef(mar)[4, -dep$which][i]
    } else {
      IndexDependentVar <- (1:(dim(dep$Z)[2] + 1))[-dep$which][i]
      trns <- mar$referenceMargin$transData[[IndexDependentVar]]
      qu <- mar$referenceMargin$mqu[IndexDependentVar]
      th <- mar$referenceMargin$mth[IndexDependentVar]
      sigma <- exp(mar$referenceMargin$models[[IndexDependentVar]]$coefficients[1])
      xi <- mar$referenceMargin$models[[IndexDependentVar]]$coefficients[2]
    }

    # for regression lines, transform to original scale
    ploty <- apply(dep$margins$q2p(yq), 2, revTransform,
      data = trns,
      qu = qu, th = th, sigma = sigma, xi = xi
    )
    # add each quantile line
    for (j in 1:length(quantiles)) {
      plots[[3]] <- plots[[3]] + geom_line(data = data.frame(x = plotx, y = ploty[, j]), mapping = aes(x, y), linetype = 2, col = ptcol)
    }

    plots
  }

  p <- unlist(lapply(1:ncol(z), plotfn), recursive = FALSE)
  np <- length(p)
  Index <- c(matrix(1:np, ncol = 3, byrow = TRUE)) # to make matrix of plots fill column-wise instead of row-wise
  p <- p[Index]

  # The loess smoother can tend to throw warnings, so suppress
  if (plot.) suppressWarnings(do.call("grid.arrange", c(p, list(nrow = 3))))
  invisible(p)
}
