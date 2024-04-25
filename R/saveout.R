
#' @title Save dcm Data to Output List Format
#'
#' @description This function transfers the data from the main \code{df}
#'              object to the output \code{out} object at the end of each
#'              run in \code{\link{dcm}}.
#'
#' @param df Main object in \code{\link{dcm}} simulations.
#' @param s Current run number.
#' @param param Param list set in \code{\link{param.dcm}}.
#' @param control Control list set in \code{\link{control.dcm}}.
#' @param out Out list passed back in for updating at runs 2+.
#'
#' @return A list with the following elements:
#' \itemize{
#'  \item \strong{param:} the epidemic parameters passed into the model through
#'        \code{\link{param.dcm}}, with additional parameters added as
#'        necessary.
#'  \item \strong{control:} the control settings passed into the model through
#'        \code{\link{control.dcm}}, with additional controls added as
#'        necessary.
#'  \item \strong{epi:} a list of data frames, one for each epidemiological
#'        output from the model.
#' }
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
  ns <- nrow(out$epi[[1]])
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
#' @description This function transfers the data from the main \code{icm_dat}
#'              class data object to the output \code{out} object at the end of
#'              each simulation in \code{\link{icm}}.
#'
#' @inheritParams prevalence.icm
#' @param s Current simulation number.
#' @param out Out list passed back in for updating at simulations 2+.
#'
#' @return
#' A list with the following elements:
#' \itemize{
#'  \item \strong{param:} the epidemic parameters passed into the model through
#'        \code{\link{param.icm}}, with additional parameters added as
#'        necessary.
#'  \item \strong{control:} the control settings passed into the model through
#'        \code{\link{control.icm}}, with additional controls added as
#'        necessary.
#'  \item \strong{epi:} a list of data frames, one for each epidemiological
#'        output from the model.
#' }
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
    for (j in seq_along(dat$epi)) {
      out$epi[[names(dat$epi)[j]]] <- data.frame(dat$epi[j])
    }
  } else {
    for (j in seq_along(dat$epi)) {
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
#' @description This function transfers the data from the main \code{netsim_dat}
#'              object to the output \code{out} object at the end of each
#'              simulation in \code{\link{netsim}}.
#'
#' @inheritParams recovery.net
#' @param s Current simulation number.
#' @param out Out list passed back in for updating at simulations 2+.
#'
#' @return
#' A list with the following elements:
#' \itemize{
#'  \item \strong{param:} the epidemic parameters passed into the model through
#'        \code{\link{param.net}}, with additional parameters added as
#'        necessary.
#'  \item \strong{control:} the control settings passed into the model through
#'        \code{\link{control.net}}, with additional controls added as
#'        necessary.
#'  \item \strong{epi:} a list of data frames, one for each epidemiological
#'        output from the model.
#'  \item \strong{stats:} a list containing two sublists, \code{nwstats} for any
#'        network statistics saved in the simulation, and \code{transmat} for
#'        the transmission matrix saved in the simulation.
#'  \item \strong{network:} a list of \code{networkDynamic} objects,
#'         one for each model simulation.
#' }
#'
#' @keywords internal
#' @export
#'
saveout.net <- function(dat, s, out = NULL) {

  if (s == 1) {
    out <- list()
    out$param <- dat$param
    out$control <- dat$control
    out$nwparam <- dat$nwparam
    out$num.nw <- dat$num.nw

    out$coef.form <- list()
    out$coef.form[[s]] <- lapply(dat$nwparam, `[[`, "coef.form")

    out$epi <- list()
    for (j in seq_along(dat$epi)) {
      out$epi[[names(dat$epi)[j]]] <- data.frame(dat$epi[j])
    }

    out$run <- list()
    out$run[[s]] <- dat$run

    out$attr.history <- list()
    out$attr.history[[s]] <- dat$attr.history

    out$raw.records <- list()
    out$raw.records[[s]] <- dat$raw.records

    out$stats <- list()
    if (dat$control$save.nwstats == TRUE) {
      ## bind rows
      nwstats <- lapply(dat$stats$nwstats, dplyr::bind_rows)
      ## compute ess
      nwstats <- lapply(nwstats, function(y) structure(y, ess = ess(y)))
      ## store as first element in list on output object
      out$stats$nwstats <- list(nwstats)
    }

    if (dat$control$save.transmat == TRUE) {
      if (!is.null(dat$stats$transmat)) {
        transmat <- dplyr::bind_rows(dat$stats$transmat)
        row.names(transmat) <- seq_len(nrow(transmat))
        out$stats$transmat[[s]] <- transmat
      } else {
        out$stats$transmat[[s]] <- dplyr::tibble()
      }
      class(out$stats$transmat) <- c("transmat", class(out$stats$transmat))
    }

    if (dat$control$save.network == TRUE) {
      ## call get_network to use most up-to-date el and attr in tergmLite case
      out$network <- list(lapply(seq_len(dat$num.nw), get_network, x = dat))
    }

    if (!is.null(dat$control$save.other)) {
      for (i in seq_along(dat$control$save.other)) {
        el.name <- dat$control$save.other[i]
        out[[el.name]] <- list(dat[[el.name]])
      }
    }

    if (dat$control$save.diss.stats == TRUE &&
          dat$control$save.network == TRUE &&
          dat$control$tergmLite == FALSE &&
          !is.null(dat$nwparam)) {

      ## for each simulated network, if dissolution model is edges-only, compute diss stats
      out$diss.stats <- list(lapply(seq_len(dat$num.nw), function(network) {
        if (dat$nwparam[[network]]$coef.diss$diss.model.type == "edgesonly") {
          toggles_to_diss_stats(tedgelist_to_toggles(as.data.frame(dat$nw[[network]])),
                                dat$nwparam[[network]]$coef.diss,
                                dat$control$nsteps,
                                dat$nw[[network]])
        } else {
          NULL
        }
      }))
    }
  }

  if (s > 1) {
    if (!is.null(dat$param$random.params.values)) {
      for (nms in names(dat$param$random.params.values)) {
        if (length(dat$param$random.params.values[[nms]]) > 1) {
          if (!is.list(out$param$random.params.values[[nms]])) {
            out$param$random.params.values[[nms]] <- list(
              out$param$random.params.values[[nms]]
            )
          }

          out$param$random.params.values[[nms]] <- c(
            out$param$random.params.values[[nms]],
            list(dat$param$random.params.values[[nms]])
          )

        } else {
          out$param$random.params.values[[nms]] <- c(
            out$param$random.params.values[[nms]],
            dat$param$random.params.values[[nms]]
          )
        }
      }
    }

    out$coef.form[[s]] <- lapply(dat$nwparam, `[[`, "coef.form")

    for (j in seq_along(dat$epi)) {
      out$epi[[names(dat$epi)[j]]][, s] <- data.frame(dat$epi[j])
    }

    out$attr.history[[s]] <- dat$attr.history
    out$raw.records[[s]] <- dat$raw.records

    if (dat$control$save.nwstats == TRUE) {
      ## bind rows
      nwstats <- lapply(dat$stats$nwstats, dplyr::bind_rows)
      ## compute ess
      nwstats <- lapply(nwstats, function(y) structure(y, ess = ess(y)))
      ## store as s'th element in list on output object
      out$stats$nwstats[[s]] <- nwstats
    }

    if (dat$control$save.transmat == TRUE) {
      if (!is.null(dat$stats$transmat)) {
        transmat <- dplyr::bind_rows(dat$stats$transmat)
        row.names(transmat) <- seq_len(nrow(transmat))
        out$stats$transmat[[s]] <- transmat
      } else {
        out$stats$transmat[[s]] <- dplyr::tibble()
      }
    }

    if (dat$control$save.network == TRUE) {
      ## call get_network to use most up-to-date el and attr in tergmLite case
      out$network[[s]] <- lapply(seq_len(dat$num.nw), get_network, x = dat)
    }

    if (!is.null(dat$control$save.other)) {
      for (i in seq_along(dat$control$save.other)) {
        el.name <- dat$control$save.other[i]
        out[[el.name]][[s]] <- dat[[el.name]]
      }
    }

    if (dat$control$save.diss.stats == TRUE &&
          dat$control$save.network == TRUE &&
          dat$control$tergmLite == FALSE &&
          !is.null(dat$nwparam)) {

      ## for each simulated network, if dissolution model is edges-only, compute diss stats
      out$diss.stats[[s]] <- lapply(seq_len(dat$num.nw), function(network) {
        if (dat$nwparam[[network]]$coef.diss$diss.model.type == "edgesonly") {
          toggles_to_diss_stats(tedgelist_to_toggles(as.data.frame(dat$nw[[network]])),
                                dat$nwparam[[network]]$coef.diss,
                                dat$control$nsteps,
                                dat$nw[[network]])
        } else {
          NULL
        }
      })
    }
  }

  ## Final processing
  if (s == dat$control$nsims) {

    # Set names for out
    simnames <- paste0("sim", seq_len(dat$control$nsims))
    for (i in as.vector(which(lapply(out$epi, class) == "data.frame"))) {
      colnames(out$epi[[i]]) <- simnames
    }

    top_lvl_elts <- c("attr.history", ".records")
    for (elt in top_lvl_elts) {
      out[[elt]] <- name_saveout_elts(out[[elt]], elt, simnames)
    }

    out$coef.form <- name_saveout_elts(out$coef.form, "coef.form", simnames)

    if (dat$control$save.nwstats == TRUE) {
      out$stats$nwstats <- name_saveout_elts(out$stats$nwstats,
                                             "stats$nwstats", simnames)
    }

    if (dat$control$save.transmat == TRUE) {
      out$stats$transmat <- name_saveout_elts(out$stats$transmat,
                                              "stats$transmat", simnames)
    }

    if (dat$control$save.network == TRUE) {
      out$network <- name_saveout_elts(out$network, "network", simnames)
    }

    if (dat$control$save.diss.stats == TRUE &&
          dat$control$save.network == TRUE &&
          dat$control$tergmLite == FALSE) {
      out$diss.stats <- name_saveout_elts(out$diss.stats, "diss.stats", simnames)
    }

    for (el.name in dat$control$save.other) {
      out[[el.name]] <- name_saveout_elts(out[[el.name]], el.name, simnames)
    }

    # Remove functions from control list
    ftodel <- grep(".FUN", names(out$control), value = TRUE)
    out$control[ftodel] <- NULL
    out$control$currsim <- NULL
    environment(out$control$nwstats.formula) <- NULL
  }

  return(out)
}


#' @title Save a List of netsim Data to Output List Format
#'
#' @description This function transfers the data from a list of the main
#'              \code{netsim_dat} objects to the output \code{out} object at the
#'              end of all simulations in \code{\link{netsim}}.
#'
#' @param dat_list A list of main \code{netsim_dat} objects in \code{netsim}
#'        simulations.
#'
#' @return
#' A list of class \code{netsim} with the following elements:
#' \itemize{
#'  \item \strong{param:} the epidemic parameters passed into the model through
#'        \code{param}, with additional parameters added as necessary.
#'  \item \strong{control:} the control settings passed into the model through
#'        \code{control}, with additional controls added as necessary.
#'  \item \strong{epi:} a list of data frames, one for each epidemiological
#'        output from the model. Outputs for base models always include the
#'        size of each compartment, as well as flows in, out of, and between
#'        compartments.
#'  \item \strong{stats:} a list containing two sublists, \code{nwstats} for any
#'        network statistics saved in the simulation, and \code{transmat} for
#'        the transmission matrix saved in the simulation. See
#'        \code{\link{control.net}} and the
#'        \href{http://www.epimodel.org/tut.html}{tutorials} for further
#'        details.
#'  \item \strong{network:} a list of \code{networkDynamic} objects,
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

name_saveout_elts <- function(elt, elt_name, simnames) {
  if (length(elt) == length(simnames)) {
    names(elt) <- simnames
  } else if (length(elt) > 0) {
    warning("The number of `", elt_name, "` is not the number of simulations.",
            "Unable to assign each to the correct simulation.")
  }
  return(elt)
}
