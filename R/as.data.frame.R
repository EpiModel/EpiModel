
#' @title Extract Model Data for Deterministic Compartmental Models
#'
#' @description This function extracts a model run from an object of class
#'              \code{dcm} into a data frame using the generic
#'              \code{as.data.frame} function.
#'
#' @param x An \code{EpiModel} object of class \code{\link{dcm}}.
#' @param run Run number for model; used for multiple-run sensitivity models. If
#'        not specified, will output data from all runs in a stacked data frame.
#' @param row.names See \code{\link{as.data.frame.default}}.
#' @param optional See \code{\link{as.data.frame.default}}.
#' @param ...  See \code{\link{as.data.frame.default}}.
#'
#' @details
#' Model output from \code{dcm} simulations are available as a data
#' frame with this helper function. The output data frame will include
#' columns for time, the size of each compartment, the overall population
#' size (the sum of compartment sizes), and the size of each flow.
#'
#' For models with multiple runs (i.e., varying parameters - see example below),
#' the default with the \code{run} parameter not specified will output all runs
#' in a single stacked data frame.
#'
#' @method as.data.frame dcm
#' @keywords extract
#' @export
#'
#' @examples
#' ## Example 1: One-group SIS model with varying act.rate
#' param <- param.dcm(inf.prob = 0.2, act.rate = seq(0.05, 0.5, 0.05),
#'                    rec.rate = 1/50)
#' init <- init.dcm(s.num = 500, i.num = 1)
#' control <- control.dcm(type = "SIS", nsteps = 10)
#' mod1 <- dcm(param, init, control)
#' as.data.frame(mod1)
#' as.data.frame(mod1, run = 1)
#' as.data.frame(mod1, run = 10)
#'
#' ## Example 2: Two-group SIR model with vital dynamics
#' param <- param.dcm(inf.prob = 0.2, inf.prob.g2 = 0.1,
#'                    act.rate = 3, balance = "g1",
#'                    rec.rate = 1/50, rec.rate.g2 = 1/50,
#'                    b.rate = 1/100, b.rate.g2 = NA,
#'                    ds.rate = 1/100, ds.rate.g2 = 1/100,
#'                    di.rate = 1/90, di.rate.g2 = 1/90,
#'                    dr.rate = 1/100, dr.rate.g2 = 1/100)
#' init <- init.dcm(s.num = 500, i.num = 1, r.num = 0,
#'                  s.num.g2 = 500, i.num.g2 = 1, r.num.g2 = 0)
#' control <- control.dcm(type = "SIR", nsteps = 10)
#' mod2 <- dcm(param, init, control)
#' as.data.frame(mod2)
#'
as.data.frame.dcm <- function(x, row.names = NULL, optional = FALSE, run,
                              ...) {

  df <- data.frame(time = x$control$timesteps)
  nruns <- x$control$nruns
  if (missing(run)) {
    run <- 1:nruns
  }

  # Output for models with 1 run
  if (nruns == 1) {
    if (length(run) > 1 || max(run) > 1) {
      stop("Maximum run is 1")
    }
    for (i in seq_along(x$epi)) {
      df[, i + 1] <- x$epi[[i]]
    }
  }

  # Output for models with multiple runs
  if (nruns > 1) {
    if (max(run) > nruns) {
      stop("Maximum run is ", nruns, call. = FALSE)
    }
    for (j in run) {
      if (j == min(run)) {
        for (i in seq_along(x$epi)) {
          df[, i + 1] <- x$epi[[i]][, j]
        }
        df$run <- j
      } else {
        tdf <- data.frame(time = 1:x$control$nsteps)
        for (i in seq_along(x$epi)) {
          tdf[, i + 1] <- x$epi[[i]][, j]
        }
        tdf$run <- j
        df <- rbind(df, tdf)
      }
    }
    df <- df[, c(ncol(df), 1:(ncol(df) - 1))]
  }

  if (nruns > 1) {
    names(df)[3:ncol(df)] <- names(x$epi)
  } else {
    names(df)[2:ncol(df)] <- names(x$epi)
  }

  return(df)
}


#' @title Extract Model Data for Stochastic Models
#'
#' @description This function extracts model simulations for objects of classes
#'              \code{icm} and \code{netsim} into a data frame using
#'              the generic \code{as.data.frame} function.
#'
#' @param x An \code{EpiModel} object of class \code{icm} or \code{netsim}.
#' @param out Data output to data frame: \code{"mean"} for row means across
#'        simulations, \code{"sd"} for row standard deviations across simulations,
#'        \code{"qnt"} for row quantiles at the level specified in \code{qval},
#'        or \code{"vals"} for values from individual simulations.
#' @param sim If \code{out="vals"}, the simulation number to output. If not
#'        specified, then data from all simulations will be output.
#' @param qval Quantile value required when \code{out="qnt"}.
#' @param row.names See \code{\link{as.data.frame.default}}.
#' @param optional See \code{\link{as.data.frame.default}}.
#' @param ...  See \code{\link{as.data.frame.default}}.
#'
#' @details
#' These methods work for both \code{icm} and \code{netsim} class models. The
#' available output includes time-specific means, standard deviations, quantiles,
#' and simulation values (compartment and flow sizes) from these stochastic model
#' classes. Means, standard deviations, and quantiles are calculated by taking the
#' row summary (i.e., each row of data is corresponds to a time step) across all
#' simulations in the model output.
#'
#' @method as.data.frame icm
#' @keywords extract
#' @export
#'
#' @examples
#' ## Stochastic ICM SIS model
#' param <- param.icm(inf.prob = 0.8, act.rate = 2, rec.rate = 0.1)
#' init <- init.icm(s.num = 500, i.num = 1)
#' control <- control.icm(type = "SIS", nsteps = 10,
#'                        nsims = 3, verbose = FALSE)
#' mod <- icm(param, init, control)
#'
#' # Default output all simulation runs, default to all in stacked data.frame
#' as.data.frame(mod)
#' as.data.frame(mod, sim = 2)
#'
#' # Time-specific means across simulations
#' as.data.frame(mod, out = "mean")
#'
#' # Time-specific standard deviations across simulations
#' as.data.frame(mod, out = "sd")
#'
#' # Time-specific quantile values across simulations
#' as.data.frame(mod, out = "qnt", qval = 0.25)
#' as.data.frame(mod, out = "qnt", qval = 0.75)
#'
#' \dontrun{
#' ## Stochastic SI network model
#' nw <- network.initialize(n = 100, directed = FALSE)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#' est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' param <- param.net(inf.prob = 0.5)
#' init <- init.net(i.num = 10)
#' control <- control.net(type = "SI", nsteps = 10, nsims = 3, verbose = FALSE)
#' mod <- netsim(est, param, init, control)
#'
#' # Same data extraction methods as with ICMs
#' as.data.frame(mod)
#' as.data.frame(mod, out = "mean")
#' as.data.frame(mod, out = "sd")
#' as.data.frame(mod, out = "qnt", qval = 0.25)
#' as.data.frame(mod, out = "qnt", qval = 0.75)
#' }
#'
as.data.frame.icm <- function(x, row.names = NULL, optional = FALSE,
                              out = "vals", sim, qval, ...) {

  df <- data.frame(time = 1:x$control$nsteps)
  nsims <- x$control$nsims

  # Will handle both icm and netsim here
  if (class(x) == "icm") {
    groups <- x$param$groups
  }
  if (class(x) == "netsim") {
    groups <- x$param$modes
  }


  if (out == "vals") {

    if (missing(sim)) {
      sim <- "all"
    }

    # Output for models with 1 sim
    if (nsims == 1) {
      for (i in seq_along(x$epi)) {
        df[, i + 1] <- x$epi[[i]]
      }
      df$sim <- 1
      df <- df[, c(ncol(df), 1:(ncol(df) - 1))]
    }

    # Output for models with multiple sims
    if (nsims > 1) {
      if (sim == "all") {
        for (j in 1:nsims) {
          if (j == 1) {
            for (i in seq_along(x$epi)) {
              df[, i + 1] <- x$epi[[i]][, j]
            }
            df$sim <- j
          } else {
            tdf <- data.frame(time = 1:x$control$nsteps)
            for (i in seq_along(x$epi)) {
              tdf[, i + 1] <- x$epi[[i]][, j]
            }
            tdf$sim <- j
            df <- rbind(df, tdf)
          }
        }
        df <- df[, c(ncol(df), 1:(ncol(df) - 1))]
      } else {
        if (sim > nsims) {
          stop(paste("Specify sim between 1 and", nsims))
        }
        for (i in seq_along(x$epi)) {
          df[, i + 1] <- x$epi[[i]][, sim]
        }
        df$sim <- sim
        df <- df[, c(ncol(df), 1:(ncol(df) - 1))]
      }
    }
  }

  ## Output means
  if (out == "mean") {
    if (nsims == 1) {
      for (i in seq_along(x$epi)) {
        df[, i + 1] <- x$epi[[i]]
      }
    }
    if (nsims > 1) {
      for (i in seq_along(x$epi)) {
        df[, i + 1] <- rowMeans(x$epi[[i]], na.rm = TRUE)
      }
    }
  }

  ## Output standard deviations
  if (out == "sd") {
    if (nsims == 1) {
      for (i in seq_along(x$epi)) {
        df[, i + 1] <- 0
      }
    }
    if (nsims > 1) {
      for (i in seq_along(x$epi)) {
        df[, i + 1] <- apply(x$epi[[i]], 1, sd, na.rm = TRUE)
      }
    }
  }

  if (out == "qnt") {
    if (missing(qval) || length(qval) > 1 || (qval > 1 | qval < 0)) {
      stop("Must specify qval as single value between 0 and 1", call. = FALSE)
    }
    for (i in seq_along(x$epi)) {
      df[, i + 1] <- apply(x$epi[[i]], 1, quantile, probs = qval,
                           na.rm = TRUE, names = FALSE)
    }
  }

  if (out == "vals") {
    names(df)[3:ncol(df)] <- names(x$epi)
  } else {
    names(df)[2:ncol(df)] <- names(x$epi)
  }

  return(df)
}


#' @method as.data.frame netsim
#' @export
#' @rdname as.data.frame.icm
as.data.frame.netsim <- function(x, row.names = NULL, optional = FALSE,
                                 out = "vals", sim, ...) {
  if (missing(sim)) {
    sim <- "all"
  }
  df <- as.data.frame.icm(x, row.names = row.names, optional = optional,
                          sim = sim, out = out, ...)
  return(df)

}

#' @title Extract Timed Edgelists netdx Objects
#'
#' @description This function extracts timed edgelists for objects of class
#'              \code{netdx} into a data frame using the generic
#'              \code{as.data.frame} function.
#'
#' @param x An \code{EpiModel} object of class \code{netdx}.
#' @param sim The simulation number to output. If not specified, then data from
#'        all simulations will be output.
#' @param row.names See \code{\link{as.data.frame.default}}.
#' @param optional See \code{\link{as.data.frame.default}}.
#' @param ...  See \code{\link{as.data.frame.default}}.
#'
#' @method as.data.frame netdx
#' @keywords extract
#' @export
#'
#' @examples
#' # Initialize and parameterize the network model
#' nw <- network.initialize(n = 100, directed = FALSE)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#'
#' # Model estimation
#' est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' # Simulate the network with netdx
#' dx <- netdx(est, nsims = 3, nsteps = 10, verbose = FALSE)
#'
#' # Extract data from the first simulation
#' as.data.frame(dx, sim = 1)
#'
#' # Extract data from all simulations
#' as.data.frame(dx)
#'
as.data.frame.netdx <- function(x, row.names = NULL, optional = FALSE,
                                sim, ...) {

  if (missing(sim)) {
    sim <- 1:x$nsims
  }
  if (max(sim) > x$nsims) {
    stop(paste("Maximum sim is", x$nsims), call. = FALSE)
  }

  if (length(sim) == 1) {
    df <- x$edgelist[[sim]]
  } else {
    for (j in sim) {
      if (j == min(sim)) {
        df <- x$edgelist[[j]]
        df$sim <- j
      } else {
        tdf <- x$edgelist[[j]]
        tdf$sim <- j
        df <- rbind(df, tdf)
      }
    }
    df <- df[, c(ncol(df), 1:(ncol(df) - 1))]
  }

  return(df)
}
