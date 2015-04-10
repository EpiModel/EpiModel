
#' @title RColorBrewer Color Ramp for EpiModel Plots
#'
#' @description Returns vector of colors consistent with a high-brightness set
#'              of colors from an \code{RColorBrewer} palette.
#'
#' @param plt \code{RColorBrewer} palette from \code{\link{brewer.pal}}
#' @param n Number of colors to return
#' @param delete.lights Delete the lightest colors from the color palette,
#'        helps with plotting in many high-contrast palettes
#'
#' @details
#' \code{RColorBrewer} provides easy access to helpful color palettes, but the
#' built-in palettes are limited to the set of colors in the existing palette.
#' This function expands the palette size to any number of colors by filling
#' in the gaps. Also, colors within the "div" and "seq" set of palettes whose
#' colors are very light (close to white) are deleted by default for better
#' visualization of plots.
#'
#' @return
#' A vector of length equal to \code{n} with a range of color values consistent
#' with an RColorBrewer color palette.
#'
#' @seealso \code{\link{RColorBrewer}}
#' @keywords colorUtils internal
#' @export
#'
#' @examples
#' # Shows a 100-color ramp for 4 RColorBrewer palettes
#' par(mfrow = c(2, 2), mar=c(1, 1, 2, 1))
#' pals <- c("Spectral", "Greys", "Blues", "Set1")
#' for (i in seq_along(pals)) {
#'  plot(1:100, 1:100, type = "n", axes = FALSE, main = pals[i])
#'  abline(v = 1:100, lwd = 6, col = brewer_ramp(100, pals[i]))
#' }
#'
brewer_ramp <- function(n, plt, delete.lights = TRUE){

  if (n < 1) {
    stop("n must be a positive integer", call. = FALSE)
  }

  bpi <- brewer.pal.info
  if (!(plt %in% row.names(bpi))) {
    stop("plt must match an RColorBrewer palette name. See RColorBrewer::brewer.pal.info",
         .call = FALSE)
  }

  pltmax <- bpi[row.names(bpi)==plt, ]$maxcolors
  pltcat <- bpi[row.names(bpi)==plt, ]$category

  if (pltcat == "div") {
    if (delete.lights == TRUE) {
      colors <- brewer.pal(pltmax, plt)[-c(4:7)]
    } else {
      colors <- brewer.pal(pltmax, plt)
    }
  }
  if (pltcat == "qual") {
    colors <- brewer.pal(pltmax, plt)
  }
  if (pltcat == "seq") {
    if (delete.lights == TRUE) {
      colors <- rev(brewer.pal(pltmax, plt)[-c(1:3)])
    } else {
      colors <- rev(brewer.pal(pltmax, plt))
    }
  }
  if (plt == "Set1") {
    colors <- brewer.pal(9, "Set1")[-6]
  }

  pal <- colorRampPalette(colors)

  return(pal(n))
}


#' @title Delete Elements from Attribute List
#'
#' @description Deletes elements from the master attribute list.
#'
#' @param attrList Attribute list.
#' @param ids ID numbers to delete from the list.
#'
#' @export
#' @keywords internal
deleteAttr <- function(attrList, ids) {
  if (length(ids) > 0) {
    attrList <- lapply(attrList, function(x) x[-ids])
  }
  return(attrList)
}


#' @title Obtain Transparent Colors
#'
#' @description Returns an RGB transparent color from any standard R color.
#'
#' @param col Vector consisting of colors, of any of the three kinds of
#'        \code{R} color specifications (named, hexadecimal, or integer; see
#'        \code{\link{col2rgb}}).
#' @param alpha Vector of transparency levels, where 0 is transparent and 1
#'        is opaque.
#' @param invisible Supresses printing of the RGB color.
#'
#' @details
#' The purpose of this function is to facilitate color transparency, which is
#' used widely in \code{EpiModel} plots. This is an internal function that is
#' not ordinarily called by the end-user. This function allows that one of col
#' or alpha may be of length greater than 1.
#'
#' @return
#' A vector of length equal to the input \code{col} vector or the \code{alpha},
#' vector, if one or the other is of length greater than 1, containing the
#' transformed color values in hexidemical format.
#'
#' @seealso \code{\link{rgb}}, \code{\link{col2rgb}}
#'
#' @export
#' @keywords colorUtils internal
#'
#' @examples
#' ## Example 1: Bubble plot with multiple length color vector
#' n <- 25
#' x <- sort(sample(1:200, n))
#' y <- 10 + 2*x + rnorm(n, 0, 10)
#' z <- rpois(n, 10)
#' cols <- transco(c("steelblue", "black"), 0.5)
#' par(mar=c(2, 2, 1, 1))
#' plot(x, y, cex = z/4, pch = 21, col = "black",
#'      bg = cols[1], lwd = 1.2, axes = FALSE,
#'      ylim = c(0, 500), xlim = c(0, 250),
#'      yaxs = "r", xaxs = "r")
#' axis(2, seq(0, 500, 100), col = "white", las = 2,
#'     cex.axis = 0.9, mgp = c(2, 0.5, 0))
#' axis(1, seq(0, 250, 50), cex.axis = 0.9,
#'      mgp = c(2, 0.5, 0))
#' abline(h = seq(100, 500, 100), col = cols[2])
#'
#' ## Example 2: Network plot with multiple length alpha vector
#' net <- network.initialize(500, directed = FALSE)
#' vcol <- transco("firebrick",
#'                 alpha = seq(0, 1, length = network.size(net)))
#' par(mar = c(0, 0, 0, 0))
#' plot(net, vertex.col = vcol, vertex.border = "grey70",
#'      vertex.cex = 1.5, edge.col = "grey50")
#'
transco <- function(col,
                    alpha = 1,
                    invisible = FALSE
                    ) {

  if (length(alpha) > 1 && length(col) > 1) {
    stop("Length of col or length of alpha must be 1")
  }

  if (alpha > 1 || alpha < 0) {
    stop("Specify alpha between 0 and 1")
  }

  newa <- floor(alpha * 255)
  t1 <- col2rgb(col, alpha = FALSE)
  t2 <- rep(NA, length(col))

  if (length(col) > 1) {
    for (i in seq_along(col)) {
      t2[i] <- rgb(t1[1,i], t1[2,i], t1[3,i], newa, maxColorValue = 255)
    }
  }
  if (length(alpha) > 1) {
    for (i in seq_along(alpha)) {
      t2[i] <- rgb(t1[1,1], t1[2,1], t1[3,1], newa[i], maxColorValue = 255)
    }
  }
  if (length(col) == 1 && length(alpha) == 1) {
    t2 <- rgb(t1[1,1], t1[2,1], t1[3,1], newa, maxColorValue = 255)
  }

  if (invisible == TRUE) {
    invisible(t2)
  } else {
    return(t2)
  }
}


#' @title Stable Sampling Function
#'
#' @description Provides a sampling function useful for dynamic simulations, in
#'              which the length of the input vector may be multiple lengths and
#'              the size of the sample may be 0.
#'
#' @param x Either a vector of one or more elements from which to choose, or a
#'        positive integer.
#' @param size Non-negative integer giving the number of items to choose.
#' @param replace Should sampling be with replacement?
#' @param prob Vector of probability weights for obtaining the elements of the
#'        vector being sampled.
#'
#' @export
#' @keywords internal
ssample <- function(x, size, replace = FALSE, prob = NULL) {

  if (length(x) > 1) {
    return(sample(x, size, replace, prob))
  }

  if (length(x) == 1 && size > 0) {
    return(x)
  }

  if (length(x) == 1 && size == 0) {
    return(NULL)
  }

}

