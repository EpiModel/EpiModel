
#' @title Save dcm Data to Output List Format
#'
#' @description This function transfers the data from the master `df`
#'              object to the output `out` object at the end of each
#'              simulation in `\link{dcm}`.
#'
#' @param df Master object in `\link{dcm}` simulations.
#' @param s Current run number.
#' @param param Param list set in `\link{param.dcm}`
#' @param control Control list set in `\link{control.dcm}`.
#' @param out Out list passed back in for updating at runs 2+.
#'
#' @keywords internal
#' @export
#'
saveout.dcm <- function(df, s, param, control, out = NULL) {

  if (s == 1) {
    out <- list()
    out$param <- param
    out$control <- control
    if (length(control$nsteps) == 1) {
      out$control$timesteps <- seq(1, control$nsteps, control$dt)
    } else {
      out$control$timesteps <- control$nsteps
    }

    out$epi <- list()
    for (j in 2:ncol(df)) {
      out$epi[[names(df)[j]]] <- data.frame(df[, j])
    }
  } else {
    for (j in 2:ncol(df)) {
      out$epi[[names(df)[j]]][, s] <- data.frame(df[, j])
    }
  }

  # Remove NA's from flows by setting last value to penultimate value
  ns <- control$nsteps
  lr.na <- sapply(out$epi, function(x) is.na(x[ns, s]) & !is.na(x[ns - 1, s]))
  wh.lr.na <- as.numeric(which(lr.na == TRUE))
  if (length(wh.lr.na) > 0) {
    for (jj in wh.lr.na) {
      out$epi[[jj]][ns, s] <- out$epi[[jj]][ns - 1, s]
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
#' @description This function transfers the data from the master `all`
#'              object to the output `out` object at the end of each
#'              simulation in `\link{icm}`.
#'
#' @param dat Master object in `icm` simulations.
#' @param s Current run number.
#' @param out Out list passed back in for updating at runs 2+.
#'
#' @keywords internal
#' @export
#'
saveout.icm <- function(dat, s, out = NULL) {

  if (s == 1) {
    out <- list()
    out$param <- dat$param
    out$control <- dat$control
    out$epi <- list()
    for (j in 1:length(dat$epi)) {
      out$epi[[names(dat$epi)[j]]] <- data.frame(dat$epi[j])
    }
  } else {
    for (j in 1:length(dat$epi)) {
      out$epi[[names(dat$epi)[j]]][, s] <- data.frame(dat$epi[j])
    }
  }

  ## Processing for final run
  if (s == dat$control$nsims) {

    # Remove functions from control list
    ftodel <- grep(".FUN", names(out$control), value = TRUE)
    out$control[ftodel] <- NULL

    # Set column names for varying list elements
    for (i in as.vector(which(lapply(out$epi, class) == "data.frame"))) {
      colnames(out$epi[[i]]) <- paste0("sim", 1:dat$control$nsims)
    }

  }

  return(out)
}


#' @title Save netsim Data to Output List Format
#'
#' @description This function transfers the data from the master `dat`
#'              object to the output `out` object at the end of each
#'              simulation in `\link{netsim}`.
#'
#' @param dat Master object in `netsim` simulations.
#' @param s Current simulation number.
#' @param out Out list passed back in for updating at simulations 2+.
#'
#' @keywords internal
#' @export
#'
saveout.net <- function(dat, s, out = NULL) {

  # Counts number of simulated networks
  num.nw <- length(dat$nw)

  if (s == 1) {
    out <- list()
    out$param <- dat$param
    out$control <- dat$control
    out$nwparam <- dat$nwparam
    out$control$num.nw <- num.nw

    out$epi <- list()
    for (j in 1:length(dat$epi)) {
      out$epi[[names(dat$epi)[j]]] <- data.frame(dat$epi[j])
    }

    out$stats <- list()
    if (dat$control$save.nwstats == TRUE) {
      out$stats$nwstats <- list(dat$stats$nwstats)
    }

    if (dat$control$tergmLite == FALSE) {
      if (!is.null(dat$stats$transmat)) {
        row.names(dat$stats$transmat) <- 1:nrow(dat$stats$transmat)
        out$stats$transmat <- list(dat$stats$transmat)
      } else {
        out$stats$transmat <- list(data.frame())
      }
      class(out$stats$transmat) <- c("transmat", class(out$stats$transmat))
      out$network <- list(dat$nw)

    }

    if (!is.null(dat$control$save.other)) {
      for (i in 1:length(dat$control$save.other)) {
        el.name <- dat$control$save.other[i]
        out[[el.name]] <- list(dat[[el.name]])
      }
    }
  }

  if (s > 1) {
    for (j in 1:length(dat$epi)) {
      out$epi[[names(dat$epi)[j]]][, s] <- data.frame(dat$epi[j])
    }

    if (dat$control$save.nwstats == TRUE) {
      out$stats$nwstats[[s]] <- dat$stats$nwstats
    }

    if (dat$control$tergmLite == FALSE) {
      if (!is.null(dat$stats$transmat)) {
        row.names(dat$stats$transmat) <- 1:nrow(dat$stats$transmat)
        out$stats$transmat[[s]] <- dat$stats$transmat
      } else {
        out$stats$transmat[[s]] <- data.frame()
      }
      out$network[[s]] <- dat$nw

    }

    if (!is.null(dat$control$save.other)) {
      for (i in 1:length(dat$control$save.other)) {
        el.name <- dat$control$save.other[i]
        out[[el.name]][[s]] <- dat[[el.name]]
      }
    }
  }

  ## Final processing
  if (s == dat$control$nsims) {

    # Set names for out
    simnames <- paste0("sim", 1:dat$control$nsims)
    for (i in as.vector(which(lapply(out$epi, class) == "data.frame"))) {
      colnames(out$epi[[i]]) <- simnames
    }
    if (dat$control$save.nwstats == TRUE) {
      names(out$stats$nwstats) <- simnames
    }

    if (dat$control$tergmLite == FALSE) {
      names(out$stats$transmat) <- simnames[1:length(out$stats$transmat)]
      names(out$network) <- simnames
    }

    # Remove functions from control list
    ftodel <- grep(".FUN", names(out$control), value = TRUE)
    out$control[ftodel] <- NULL
    out$control$currsim <- NULL
    environment(out$control$nwstats.formula) <- NULL

    if (!("temp" %in% dat$control$save.other)) {
      out$temp <- NULL
    }

  }

  return(out)
}


#' @title Save a list of netsim Data to Output List Format
#'
#' @description This function transfers the data from a list of the master
#'              `dat` objects to the output `out` object at the end of
#'              all simulations in `\link{netsim}`.
#'
#' @param dat_list A list of Master objects in `netsim` simulations.
#'
#' @return
#' A list of class `netsim` with the following elements:
#' \itemize{
#'  \item \strong{param:} the epidemic parameters passed into the model through
#'        `param`, with additional parameters added as necessary.
#'  \item \strong{control:} the control settings passed into the model through
#'        `control`, with additional controls added as necessary.
#'  \item \strong{epi:} a list of data frames, one for each epidemiological
#'        output from the model. Outputs for base models always include the
#'        size of each compartment, as well as flows in, out of, and between
#'        compartments.
#'  \item \strong{stats:} a list containing two sublists, `nwstats` for any
#'        network statistics saved in the simulation, and `transmat` for
#'        the transmission matrix saved in the simulation. See
#'        `\link{control.net}` and the Tutorial for further details.
#'  \item \strong{network:} a list of `networkDynamic` objects,
#'         one for each model simulation.
#' }
#'
#' @keywords internal
#' @export
#'
process_out.net <- function(dat_list) {
  for (s in seq_along(dat_list)) {
    # Set output
    if (s == 1) {
      out <- saveout.net(dat_list[[s]], s)
    } else {
      out <- saveout.net(dat_list[[s]], s, out)
    }
  }
  class(out) <- "netsim"

  return(out)
}
