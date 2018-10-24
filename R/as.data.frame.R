
#' @title Extract Model Data for Deterministic Compartmental Models
#'
#' @description This function extracts a model run from an object of class
#'              \code{dcm} into a data frame using the generic
#'              \code{as.data.frame} function.
#'
#' @param x An \code{EpiModel} object of class \code{\link{dcm}}.
#' @param run Run number for model; used for multiple-run sensitivity models.
#' @param row.names See \code{\link{as.data.frame.default}}.
#' @param optional See \code{\link{as.data.frame.default}}.
#' @param ...  See \code{\link{as.data.frame.default}}.
#'
#' @details
#' Model output from a \code{dcm} simulation are available as a data
#' frame with this helper function. The output data frame will include
#' columns for time, the size of each compartment, the overall population
#' size (the sum of compartment sizes), and the size of each flow.
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
#' control <- control.dcm(type = "SIS", nsteps = 500)
#' mod1 <- dcm(param, init, control)
#' head(as.data.frame(mod1, run = 1))
#' head(as.data.frame(mod1, run = 10))
#'
#' ## Example 2: Two-group SIR model with vital dynamics
#' param <- param.dcm(inf.prob = 0.2, inf.prob.g2 = 0.1,
#'                    act.rate = 3, balance = "g1",
#'                    rec.rate = 1/50, rec.rate.g2 = 1/50,
#'                    a.rate = 1/100, a.rate.g2 = NA,
#'                    ds.rate = 1/100, ds.rate.g2 = 1/100,
#'                    di.rate = 1/90, di.rate.g2 = 1/90,
#'                    dr.rate = 1/100, dr.rate.g2 = 1/100)
#' init <- init.dcm(s.num = 500, i.num = 1, r.num = 0,
#'                  s.num.g2 = 500, i.num.g2 = 1, r.num.g2 = 0)
#' control <- control.dcm(type = "SIR", nsteps = 500)
#' mod2 <- dcm(param, init, control)
#' head(as.data.frame(mod2))
#' tail(as.data.frame(mod2))
#'
as.data.frame.dcm <- function(x, row.names = NULL, optional = FALSE, run = 1,
                              ...) {

  df <- data.frame(time = x$control$timesteps)
  nruns <- x$control$nruns

  # Output for models with 1 run
  if (nruns == 1) {
    if (run > 1) {
      stop("Specify run = 1")
    }
    for (i in seq_along(x$epi)) {
      df[, i + 1] <- x$epi[[i]]
    }
  }

  # Output for models with multiple runs
  if (nruns > 1) {
    if (run > nruns) {
      stop(paste("Specify run between 1 and", nruns))
    }
    for (i in seq_along(x$epi)) {
      df[, i + 1] <- x$epi[[i]][, run]
    }
  }

  names(df)[2:ncol(df)] <- names(x$epi)

  return(df)
}


#' @title Extract Model Data for Stochastic Models
#'
#' @description This function extracts model simulations for objects of classes
#'              \code{icm} and \code{netsim} into a data frame using
#'              the generic \code{as.data.frame} function.
#'
#' @param x An \code{EpiModel} object of class \code{icm} or \code{netsim}.
#' @param sim If \code{out="vals"}, the simulation number to output, or the default of
#'        \code{out="all"}, which outputs data from all simulations bound together.
#' @param out Data output to data frame: \code{"mean"} for row means across
#'        simulations, \code{"sd"} for row standard deviations across simulations,
#'        \code{"qnt"} for row quantiles at the level specified in \code{qval},
#'        or \code{"vals"} for values from one individuals simulation(s).
#' @param qval Quantile value necessary when \code{out="qnt"}.
#' @param row.names See \code{\link{as.data.frame.default}}.
#' @param optional See \code{\link{as.data.frame.default}}.
#' @param ...  See \code{\link{as.data.frame.default}}.
#'
#' @details
#' These methods work for both \code{icm} and \code{netsim}
#' class models. The available output includes time-specific means,
#' standard deviations, quantiles, and simulation values (compartment and flow
#' sizes from one simulation) from these stochastic model classes. Means and
#' standard deviations are calculated by taking the row summary across all
#' simulations for each time step in the model output.
#'
#' @method as.data.frame icm
#' @keywords extract
#' @export
#'
#' @examples
#' ## Stochastic ICM SIS model with 5 simulations
#' param <- param.icm(inf.prob = 0.8, act.rate = 2, rec.rate = 0.1)
#' init <- init.icm(s.num = 500, i.num = 1)
#' control <- control.icm(type = "SIS", nsteps = 25,
#'                        nsims = 2, verbose = FALSE)
#' mod <- icm(param, init, control)
#'
#' # Default output is mean across simulations
#' as.data.frame(mod)
#'
#' # Standard deviations of simulations
#' as.data.frame(mod, out = "sd")
#'
#' # Quantile values for interquartile interval
#' as.data.frame(mod, out = "qnt", qval = 0.25)
#' as.data.frame(mod, out = "qnt", qval = 0.75)
#'
#' # Individual simulation runs, with default sim="all"
#' as.data.frame(mod, out = "vals")
#' as.data.frame(mod, out = "vals", sim = 2)
#'
#' ## Stochastic SI network model
#' nw <- network.initialize(n = 100, directed = FALSE)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#' est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' param <- param.net(inf.prob = 0.5)
#' init <- init.net(i.num = 10)
#' control <- control.net(type = "SI", nsteps = 10, nsims = 2, verbose = FALSE)
#' mod <- netsim(est, param, init, control)
#'
#' as.data.frame(mod)
#' as.data.frame(mod, out = "vals")
#'
as.data.frame.icm <- function(x, row.names = NULL, optional = FALSE,
                              sim = "all", out = "mean", qval, ...) {

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
                                 sim = "all", out = "mean", ...) {

  df <- as.data.frame.icm(x, row.names = row.names, optional = optional,
                          sim = sim, out = out, ...)
  return(df)

}
