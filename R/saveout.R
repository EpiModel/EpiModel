
#' @title Save dcm Data to Output List Format
#'
#' @description This function transfers the data from the master \code{df}
#'              object to the output \code{out} object at the end of each
#'              simulation in \code{\link{dcm}}.
#'
#' @param all master object in \code{df} simulations.
#' @param s current run number.
#' @param param param list set in \code{\link{param.dcm}}
#' @param control control list set in \code{\link{control.dcm}}.
#' @param out out list passed back in for updating at runs 2+.
#'
#' @keywords internal
#' @export
#'
saveout.dcm <- function(df, s, param, control, out) {

  if (s == 1) {
    out <- list()
    out$param <- param
    out$control <- control
    out$epi <- list()
    for (j in 2:ncol(df)) {
      out$epi[[names(df)[j]]] <- data.frame(df[, j])
    }
  } else {
    for (j in 2:ncol(df)) {
      out$epi[[names(df)[j]]][, s] <- data.frame(df[, j])
    }
  }

  ## Processing for final run
  if (s == control$nruns) {

    # Set column names for varying list elements
    for (s in as.vector(which(lapply(out$epi, class) == "data.frame"))) {
      colnames(out$epi[[s]]) <- paste0("run", 1:control$nruns)
    }

  }

  return(out)
}


#' @title Save icm Data to Output List Format
#'
#' @description This function transfers the data from the master \code{all}
#'              object to the output \code{out} object at the end of each
#'              simulation in \code{\link{icm}}.
#'
#' @param all master object in \code{icm} simulations.
#' @param s current run number.
#' @param out out list passed back in for updating at runs 2+.
#'
#' @keywords internal
#' @export
#'
saveout.icm <- function(all, s, out) {

  if (s == 1) {
    out <- list()
    out$param <- all$param
    out$control <- all$control
    out$epi <- list()
    for (j in 1:length(all$out)) {
      out$epi[[names(all$out)[j]]] <- data.frame(all$out[j])
    }
  } else {
    for (j in 1:length(all$out)) {
      out$epi[[names(all$out)[j]]][, s] <- data.frame(all$out[j])
    }
  }

  ## Processing for final run
  if (s == all$control$nsims) {

    # Remove functions from control list
    ftodel <- grep(".FUN", names(out$control), value = TRUE)
    out$control[ftodel] <- NULL

    # Set column names for varying list elements
    for (i in as.vector(which(lapply(out$epi, class) == "data.frame"))) {
      colnames(out$epi[[i]]) <- paste0("sim", 1:all$control$nsims)
    }

  }

  return(out)
}


#' @title Save netsim Data to Output List Format
#'
#' @description This function transfers the data from the master \code{all}
#'              object to the output \code{out} object at the end of each
#'              simulation in \code{\link{netsim}}.
#'
#' @param all master object in \code{netsim} simulations.
#' @param s current simulation number.
#' @param out out list passed back in for updating at simulations 2+.
#'
#' @keywords internal
#' @export
#'
saveout.net <- function(all, s, out) {

  if (s == 1) {
    out <- list()
    out$param <- all$param
    out$control <- all$control

    out$epi <- list()
    for (j in 1:length(all$out)) {
      out$epi[[names(all$out)[j]]] <- data.frame(all$out[j])
    }

    out$stats <- list()
    if (all$control$save.nwstats == TRUE) {
      out$stats$nwstats <- list(all$stats$nwstats)
    }
    if (all$control$save.transmat == TRUE) {
      if (!is.null(all$stats$transmat)) {
        row.names(all$stats$transmat) <- 1:nrow(all$stats$transmat)
        out$stats$transmat <- list(all$stats$transmat)
      } else {
        out$stats$transmat <- list(data.frame())
      }
    }
    if (all$control$save.network == TRUE) {
      out$network <- list(all$nw)
    }
  }

  if (s > 1) {
    for (j in 1:length(all$out)) {
      out$epi[[names(all$out)[j]]][, s] <- data.frame(all$out[j])
    }

    if (all$control$save.nwstats == TRUE) {
      out$stats$nwstats[[s]] <- all$stats$nwstats
    }
    if (all$control$save.transmat == TRUE) {
      if (!is.null(all$stats$transmat)) {
        row.names(all$stats$transmat) <- 1:nrow(all$stats$transmat)
        out$stats$transmat[[s]] <- all$stats$transmat
      } else {
        out$stats$transmat[[s]] <- data.frame()
      }
    }
    if (all$control$save.network == TRUE) {
      out$network[[s]] <- all$nw
    }
  }

  ## Final processing
  if (s == all$control$nsims) {

    # Set names for out
    simnames <- paste0("sim", 1:all$control$nsims)
    for (i in as.vector(which(lapply(out$epi, class) == "data.frame"))) {
      colnames(out$epi[[i]]) <- simnames
    }
    if (all$control$save.nwstats == TRUE) {
      names(out$stats$nwstats) <- simnames

      # Pull formation terms
      formation <- all$nwparam$formation
      a <- summary(update.formula(all$nw ~ ., formation), at = 1)
      if (ncol(a) == 1) {
        attributes(out$stats$nwstats)$formation.terms <- rownames(a)
      } else {
        attributes(out$stats$nwstats)$formation.terms <- colnames(a)
      }

      # Pull target statistics
      attributes(out$stats$nwstats)$target.stats <- all$nwparam$target.stats

    }
    if (all$control$save.transmat == TRUE) {
      names(out$stats$transmat) <- simnames[1:length(out$stats$transmat)]
    }
    if (all$control$save.network == TRUE) {
      names(out$network) <- simnames
    }

    # Remove functions from control list
    ftodel <- grep(".FUN", names(out$control), value = TRUE)
    out$control[ftodel] <- NULL

    out$temp <- NULL
  }

  return(out)
}


