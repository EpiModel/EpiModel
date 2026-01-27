
#' @title Extract Model Data for Deterministic Compartmental Models
#'
#' @description This function extracts a model run from an object of class
#'              `dcm` into a data frame using the generic
#'              `as.data.frame` function.
#'
#' @param x An `EpiModel` object of class [dcm()].
#' @param run Run number for model; used for multiple-run sensitivity models. If
#'        not specified, will output data from all runs in a stacked data frame.
#' @param row.names See [as.data.frame.default()].
#' @param optional See [as.data.frame.default()].
#' @param ...  See [as.data.frame.default()].
#'
#' @details
#' Model output from [dcm()] simulations are available as a data
#' frame with this helper function. The output data frame will include
#' columns for time, the size of each compartment, the overall population
#' size (the sum of compartment sizes), and the size of each flow.
#'
#' For models with multiple runs (i.e., varying parameters - see example below),
#' the default with the `run` parameter not specified will output all runs
#' in a single stacked data frame.
#'
#' @return A data frame containing the data from `x`.
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
#'                    a.rate = 1/100, a.rate.g2 = NA,
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
        tdf <- data.frame(time = seq(1, x$control$nsteps, x$control$dt))
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
#'              `icm` and `netsim` into a data frame using
#'              the generic `as.data.frame` function.
#'
#' @param x An `EpiModel` object of class `icm` or `netsim`.
#' @param out Data output to data frame: `"mean"` for row means across
#'        simulations, `"sd"` for row standard deviations across
#'        simulations, `"qnt"` for row quantiles at the level specified in
#'        `qval`, or `"vals"` for values from individual simulations.
#' @param sim If `out="vals"`, the simulation number to output. If not
#'        specified, then data from all simulations will be output.
#' @param qval Quantile value required when `out="qnt"`.
#' @param row.names See [as.data.frame.default()].
#' @param optional See [as.data.frame.default()].
#' @param repair What to do with epi trackers that are too short. "drop" will
#'        remove theme from the output, "pad" will add `NA` rows to get to the
#'        right size(Default = "drop").
#' @param ...  See [as.data.frame.default()].
#'
#' @details
#' These methods work for both `icm` and `netsim` class models. The
#' available output includes time-specific means, standard deviations,
#' quantiles, and simulation values (compartment and flow sizes) from these
#' stochastic model classes. Means, standard deviations, and quantiles are
#' calculated by taking the row summary (i.e., each row of data is corresponds
#' to a time step) across all simulations in the model output.
#'
#' @return A data frame containing the data from `x`.
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
#' ## Stochastic SI Network Model
#' nw <- network_initialize(n = 100)
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
#' as.data.frame(mod, sim = 2)
#' as.data.frame(mod, out = "mean")
#' as.data.frame(mod, out = "sd")
#' as.data.frame(mod, out = "qnt", qval = 0.25)
#' as.data.frame(mod, out = "qnt", qval = 0.75)
#' }
#'
as.data.frame.icm <- function(x, row.names = NULL, optional = FALSE,
                              out = "vals", sim = NULL, qval = NULL,
                              repair = "drop", ...) {
  nsims <- x$control$nsims
  sim <- if (is.null(sim)) seq_len(nsims) else sim
  if (max(sim) > nsims) {
    stop("Maximum sim is ", nsims, call. = FALSE)
  }

  # Find and repair EPIs with the wrong length
  n_steps <- x$control$nsteps
  epi_lengths <- vapply(x$epi, nrow, 0)
  broken_epi <- names(x$epi)[epi_lengths != n_steps]
  if (length(broken_epi) > 0) {
    if (repair == "drop") {
      warning("Dropping: ", paste0(broken_epi, sep = ", "), " - wrong length")
      x$epi <- x$epi[epi_lengths == n_steps]
    } else if (repair == "pad") {
      warning("Padding with NA: ", paste0(broken_epi, sep = ", "),
              " - wrong length")
      x$epi[broken_epi] <- lapply(x$epi[broken_epi], \(d) {
        n_miss <- n_steps - nrow(d)
        padding_rows <- as.data.frame(matrix(NA, nrow = n_miss, ncol = ncol(d)))
        colnames(padding_rows) <- colnames(d)
        rbind(d, padding_rows)
      })
    } else {
      stop(
        "Some EPIs had the wrong size: ", paste0(broken_epi, sep = ", "),
        "No repair method set, stopping"
      )
    }
  }

  # How to combine the different simulations
  combine_fun <- if (out == "vals") {
    function(epi) Reduce(c, epi[, sim, drop = FALSE])
  } else if (out == "mean") {
    function(epi) rowMeans(epi[, sim, drop = FALSE], na.rm = TRUE)
  } else if (out == "sd") {
    function(epi) apply(epi[, sim, drop = FALSE], 1, sd, na.rm = TRUE)
  } else if (out == "qnt") {
    if (is.null(qval) || length(qval) > 1 || (qval > 1 || qval < 0)) {
      stop("Must specify qval as single value between 0 and 1", call. = FALSE)
    }
    function(epi) {
      apply(epi[, sim, drop = FALSE], 1, quantile,
            probs = qval, na.rm = TRUE, names = FALSE)
    }
  }
  epi_ls <- lapply(x$epi, combine_fun)

  if (out == "vals") { # add columns `time` and `sim`
    common_elts <- list(
      sim = rep(sim, each = n_steps),
      time = rep(seq_len(n_steps), length(sim))
    )
    df <- as.epi.data.frame(as.data.frame(c(common_elts, epi_ls)))
  } else { # add column `time`
    df <- as.data.frame(c(list(time = seq_len(n_steps)), epi_ls))
  }

  return(df)
}


#' @method as.data.frame netsim
#' @export
#' @rdname as.data.frame.icm
as.data.frame.netsim <- function(x, row.names = NULL, optional = FALSE,
                                 out = "vals", sim = NULL, repair = "drop",
                                 ...) {
  as.data.frame.icm(x, row.names = row.names, optional = optional,
                    sim = sim, out = out, repair = repair, ...)
}

#' @title Extract Timed Edgelists for netdx Objects
#'
#' @description This function extracts timed edgelists for objects of class
#'              `netdx` into a data frame using the generic
#'              `as.data.frame` function.
#'
#' @param x An `EpiModel` object of class `netdx`.
#' @param sim The simulation number to output. If not specified, then data from
#'        all simulations will be output.
#' @param row.names See [as.data.frame.default()].
#' @param optional See [as.data.frame.default()].
#' @param ...  See [as.data.frame.default()].
#'
#' @return A data frame containing the data from `x`.
#'
#' @method as.data.frame netdx
#' @keywords extract
#' @export
#'
#' @examples
#' # Initialize and parameterize the network model
#' nw <- network_initialize(n = 100)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#'
#' # Model estimation
#' est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' # Simulate the network with netdx
#' dx <- netdx(est, nsims = 3, nsteps = 10, keep.tedgelist = TRUE,
#'             verbose = FALSE)
#'
#' # Extract data from the first simulation
#' as.data.frame(dx, sim = 1)
#'
#' # Extract data from all simulations
#' as.data.frame(dx)
#'
as.data.frame.netdx <- function(x, row.names = NULL, optional = FALSE,
                                sim, ...) {

  if (is.null(x$tedgelist)) {
    stop("Edgelist not saved in netdx object, check keep.tedgelist parameter",
         call. = FALSE)
  }

  if (missing(sim)) {
    sim <- 1:x$nsims
  }
  if (max(sim) > x$nsims) {
    stop(paste("Maximum sim is", x$nsims), call. = FALSE)
  }

  if (length(sim) == 1) {
    df <- x$tedgelist[[sim]]
  } else {
    for (j in sim) {
      if (j == min(sim)) {
        df <- x$tedgelist[[j]]
        df$sim <- j
      } else {
        tdf <- x$tedgelist[[j]]
        tdf$sim <- j
        df <- rbind(df, tdf)
      }
    }
    df <- df[, c(ncol(df), 1:(ncol(df) - 1))]
  }

  return(df)
}

#' Validate and Convert to epi.data.frame
#'
#' This methods ensures that the `data.frame` is correctly formatted as an
#' `epi.data.frame`
#'
#' @param df A `data.frame` to convert into an `epi.data.frame`
#'
#' @export
as.epi.data.frame <- function(df) {
  if (!inherits(df, "data.frame"))
    stop("Input must be a `data.frame`")
  if (!all(c("time", "sim") %in% names(df)))
    stop("Input must contain a `time` and a `sim` column")
  n_steps <- max(df$time)
  if (!all(by(df$time, df$sim, length) == n_steps))
    stop("Input must have the same number of time step per simulation")
  class(df) <- c("epi.data.frame", class(df))
  return(df)
}