
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
brewer_ramp <- function(n, plt, delete.lights = TRUE) {

  if (n < 1) {
    stop("n must be a positive integer", call. = FALSE)
  }

  bpi <- brewer.pal.info
  if (!(plt %in% row.names(bpi))) {
    stop("plt must match an RColorBrewer palette name. See
         RColorBrewer::brewer.pal.info",
         .call = FALSE)
  }

  pltmax <- bpi[row.names(bpi) == plt, ]$maxcolors
  pltcat <- bpi[row.names(bpi) == plt, ]$category

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

  if (class(attrList) != "list") {
    stop("attrList must be a list", call. = FALSE)
  }

  attr_length <- vapply(attrList, length, numeric(1))
  expected_length <- length(attrList[["active"]])
  wrong_length_attr <- names(attr_length)[attr_length != expected_length]

  if (length(wrong_length_attr > 0)) {
    stop(
      "The following attributes do not have the correct number of elements: \n",
      paste0(wrong_length_attr, collapse = ", "),
      "\n\n", "Check if they are initiated when new nodes are created."
    )
  }

  if (length(ids) > 0) {
    attrList <- lapply(attrList, function(x) x[-ids])
  }
  return(attrList)
}

#' @title Delete Elements from Attribute List
#'
#' @description Deletes elements from the master attribute list.
#'
#' @param dat Master list object containing a sublist of attributes.
#' @param ids ID numbers to delete from the list.
#'
#' @export
#' @keywords internal
delete_attr <- function(dat, ids) {
  attrList <- dat$attr

  if (class(attrList) != "list") {
    stop("dat object does not contain a valid attribute list", call. = FALSE)
  }
  if (length(unique(sapply(attrList, length))) != 1) {
    stop("attribute list must be rectangular (same number of obs per element)")
  }

  if (length(ids) > 0) {
    attrList <- lapply(attrList, function(x) x[-ids])
  }

  dat$attr <- attrList
  return(dat)
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


#' @title Add New Epidemiology Variables
#'
#' @description Inspired by \code{dplyr::mutate}, \code{mutate_epi} adds new
#'              variables to the epidemiological and related variables within
#'              simulated model objects of any class in \code{EpiModel}.
#'
#' @param x An \code{EpiModel} object of class \code{dcm}, \code{icm}, or
#'        \code{netsim}.
#' @param ... Name-value pairs of expressions (see examples below).
#'
#' @export
#'
#' @examples
#' # DCM example
#' param <- param.dcm(inf.prob = 0.2, act.rate = 0.25)
#' init <- init.dcm(s.num = 500, i.num = 1)
#' control <- control.dcm(type = "SI", nsteps = 500)
#' mod1 <- dcm(param, init, control)
#' mod1 <- mutate_epi(mod1, prev = i.num/num)
#' plot(mod1, y = "prev")
#'
#' # Network model example
#' nw <- network_initialize(n = 100)
#' nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#' est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.15)
#' init <- init.net(i.num = 1, i.num.g2 = 0)
#' control <- control.net(type = "SI", nsteps = 10, nsims = 3,
#'                        verbose = FALSE)
#' mod1 <- netsim(est1, param, init, control)
#' mod1
#'
#' # Add the prevalences to the dataset
#' mod1 <- mutate_epi(mod1, i.prev = i.num / num,
#'                          i.prev.g2 = i.num.g2 / num.g2)
#' plot(mod1, y = c("i.prev", "i.prev.g2"), qnts = 0.5, legend = TRUE)
#'
#' # Add incidence rate per 100 person years (assume time step = 1 week)
#' mod1 <- mutate_epi(mod1, ir100 = 5200*(si.flow + si.flow.g2) /
#'                                       (s.num + s.num.g2))
#' as.data.frame(mod1)
#' as.data.frame(mod1, out = "mean")
#'
mutate_epi <- function(x, ...) {

  dt <- lazy_dots(...)
  ndat <- lazy_eval(dt, x$epi)

  not.df <- which(sapply(ndat, class) != "data.frame")
  if (length(not.df) > 0) {
    for (jj in not.df) {
      ndat[jj][[1]] <- data.frame(rep(ndat[jj][[1]],
                                      length.out = x$control$nsteps))
      names(ndat[[jj]]) <- "run1"
    }
  }

  x$epi <- c(x$epi, ndat)
  return(x)

}

#' @title Apportion Least-Remainder Method
#'
#' @description Apportions a vector of values given a specified frequency
#'              distribution of those values such that the length of the output
#'              is robust to rounding and other instabilities.
#'
#' @param vector.length Length for the output vector.
#' @param values Values for the output vector.
#' @param proportions Proportion distribution with one number for each value.
#'        This must sum to 1.
#' @param shuffled If \code{TRUE}, randomly shuffle the order of the vector.
#'
#' @export
#'
apportion_lr <- function(vector.length, values,
                         proportions, shuffled = FALSE) {

  if (vector.length != round(vector.length)) {
    stop("argument vector.length must be a positive integer")
  }
  if (vector.length <= 0) {
    stop("argument vector.length must be a positive integer")
  }
  if (is.vector(values) == FALSE) {
    stop("argument values must be a vector")
  }
  if (!(length(proportions) == length(values) && round(sum(proportions),
                                                       10) == 1) &&
      (!(length(proportions) == length(values) - 1 && round(sum(proportions),
                                                            10) <= 1 &&
         round(sum(proportions), 10) >= 0))) {
    stop("error in proportions length or proportions sum")
  }

  if (length(proportions) == length(values) - 1) {
    proportions <- c(proportions, 1 - round(sum(proportions), 10))
  }
  result <- rep(NA, vector.length)
  exp.nums <- proportions * vector.length
  counts <- floor(exp.nums)
  remainders <- exp.nums - counts
  leftovers <- vector.length - sum(counts)
  if (leftovers > 0) {
    additions <- order(remainders, decreasing = TRUE)[1:leftovers]
    counts[additions]   <- counts[additions] + 1
  }
  result <- rep(values, counts)
  if (shuffled == TRUE) {
    result <- sample(result, length(result))
  }

  return(result)
}


#' @title Message to find in which module a `condition` occured
#'
#' @description Returns vector of colors consistent with a high-brightness set
#'              of colors from an \code{RColorBrewer} palette.
#'
#' @param cond The type of `condition` handled (message, warning, error)
#' @param module The name of the module where the `condition` occured
#' @param at The time step the `condition` occured
#' @param msg The `condition`'s message
#'
#' @return A formatted string describing where and whe the `condition` occured
#'         as well as the `condition`'s message.
#'
#' @keywords internal
netsim_cond_msg <- function(cond, module, at, msg) {
  paste0("\n\tA ", cond, " occured in module '", module, "' at step ", at)
}
