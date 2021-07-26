
#' @title Extract networkDynamic and network Objects from Network Simulations
#'
#' @description Extracts the \code{networkDynamic} object from a either a
#'              network epidemic model object generated with \code{netsim} or a
#'              network diagnostic simulation generated with \code{netdx}, with
#'              the option to collapse the extracted \code{networkDynamic}
#'              object down to a static \code{network} object.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netsim}} or
#'        \code{\link{netdx}}.
#' @param sim Simulation number of extracted network.
#' @param network Network number, for \code{netsim} objects with multiple
#'        overlapping networks (advanced use, and not applicable to \code{netdx}
#'        objects).
#' @param collapse If \code{TRUE}, collapse the \code{networkDynamic} object to
#'        a static \code{network} object at a specified time step.
#' @param at If \code{collapse} is \code{TRUE}, the time step at which the
#'        extracted network should be collapsed.
#' @param ergm.create.nd If \code{TRUE} and \code{x} contains a static ERGM
#'        (i.e., a netest model with duration = 1), then create a
#'        \code{networkDynamic} object from the stored list of static
#'        \code{network} objects; if \code{FALSE}, output the network list
#'        directly.
#'
#' @details
#' This function requires that the \code{networkDynamic} is saved during the
#' network simulation while running either \code{\link{netsim}} or
#' \code{\link{netdx}}. For the former, that is specified with the
#' \code{tergmLite} parameter in \code{\link{control.net}} set to \code{FALSE}.
#' For the latter, that is specified with the \code{keep.tedgelist} parameter
#' directly in \code{\link{netdx}}.
#'
#' @keywords extract
#' @export
#'
#' @examples
#' # Set up network and TERGM formiula
#' nw <- network_initialize(n = 100)
#' nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#'
#' # Estimate the model
#' est <- netest(nw, formation, target.stats, coef.diss)
#'
#' # Run diagnostics, saving the networkDynamic objects
#' dx <- netdx(est, nsteps = 10, nsims = 3, keep.tnetwork = TRUE,
#'      verbose = FALSE)
#'
#' # Extract the network for simulation 2 from dx object
#' get_network(dx, sim = 2)
#'
#' # Extract and collapse the network from simulation 1 at time step 5
#' get_network(dx, collapse = TRUE, at = 5)
#'
#' # Parameterize the epidemic model, and simulate it
#' param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.15)
#' init <- init.net(i.num = 10, i.num.g2 = 10)
#' control <- control.net(type = "SI", nsteps = 10, nsims = 3, verbose = FALSE)
#' mod <- netsim(est, param, init, control)
#'
#' # Extract the network for simulation 2 from mod object
#' get_network(mod, sim = 2)
#'
#' ## Extract and collapse the network from simulation 1 at time step 5
#' get_network(mod, collapse = TRUE, at = 5)
#'
get_network <- function(x, sim = 1, network = 1, collapse = FALSE, at,
                        ergm.create.nd = TRUE) {

  ## Warnings and checks ##
  if (!(class(x) %in% c("netsim", "netdx"))) {
    stop("x must be of class netsim or netdx", call. = FALSE)
  }

  nsims <- ifelse(class(x) == "netsim", x$control$nsims, x$nsims)
  if (length(sim) > 1 || sim > nsims) {
    stop("Specify a single sim between 1 and ", nsims, call. = FALSE)
  }

  if (class(x) == "netsim") {
    if (x$control$tergmLite == TRUE) {
      stop("Network object not saved in netsim object when 'tergmLite = TRUE'.
           Check control.net settings.",
           call. = FALSE)
    }
    if (x$control$tergmLite == TRUE || is.null(x$network)) {
      stop("Network object not saved in netsim object.
            Check control.net settings.", call. = FALSE)
    }
  } else if (class(x) == "netdx") {
    if (is.null(x$network)) {
      stop("Network object not saved in netdx object.
            Check keep.tnetwork parameter", call. = FALSE)
    }
  }

  if (class(x) == "netsim") {
    if (network > x$control$num.nw) {
      stop("Specify network between 1 and ", x$control$num.nw, call. = FALSE)
    }
  }

  nsteps <- ifelse(class(x) == "netsim", x$control$nsteps, x$nsteps)
  if (collapse == TRUE && (missing(at) || at > nsteps)) {
    stop("Specify collapse time step between 1 and ", nsteps,
         call. = FALSE)
  }
  ## Extraction ##
  if (class(x) == "netsim") {
    if (x$control$isTERGM == TRUE) {
      out <- x$network[[sim]][[network]]
    } else if (x$control$isTERGM == FALSE) {
      if (ergm.create.nd == TRUE) {
        out <- suppressMessages(
            networkDynamic(network.list = x$network[[sim]], start = 1,
                           base.net = x$network[[sim]][[1]], create.TEAs = TRUE))
      } else {
        out <- x$network[[sim]]
      }
    }
  } else if (class(x) == "netdx") {
    out <- x$network[[sim]]
  }

  ## Collapsing
  if (collapse == TRUE) {
    out <- network.collapse(out, at = at)
  }

  return(out)
}


#' @title Extract Transmissions Matrix from Network Epidemic Model
#'
#' @description Extracts the matrix of transmission data for each transmission
#'              event that occured within a network epidemic model.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netsim}}.
#' @param sim Simulation number of extracted network.
#'
#' @return
#' A data frame with the following columns
#' \itemize{
#'  \item \strong{at:} the time step at which the transmission occurred.
#'  \item \strong{sus:} the ID number of the susceptible (newly infected) node.
#'  \item \strong{inf:} the ID number of the infecting node.
#'  \item \strong{infDur:} the duration of the infecting node's disease at the
#'        time of the transmission.
#'  \item \strong{transProb:} the probability of transmission per act.
#'  \item \strong{actRate:} the rate of acts per unit time.
#'  \item \strong{finalProb:} the final transmission probability for the
#'        transmission event.
#' }
#'
#' @keywords extract
#' @export
#'
#' @examples
#' ## Simulate SI epidemic on two-group Bernoulli random graph
#' nw <- network_initialize(n = 100)
#' nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#' est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#' param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.15)
#' init <- init.net(i.num = 10, i.num.g2 = 10)
#' control <- control.net(type = "SI", nsteps = 10, nsims = 3, verbose = FALSE)
#' mod <- netsim(est, param, init, control)
#'
#' ## Extract the transmission matrix from simulation 2
#' get_transmat(mod, sim = 2)
#'
get_transmat <- function(x, sim = 1) {

  ## Warnings and checks
  if (class(x) != "netsim") {
    stop("x must be of class netsim", call. = FALSE)
  }

  if (sim > x$control$nsims) {
    stop("Specify sim between 1 and ", x$control$nsims, call. = FALSE)
  }

  if (x$control$save.transmat == FALSE || is.null(x$stats$transmat)) {
    stop("transmat not saved in netsim object, check control.net settings",
         call. = FALSE)
  }

  ## Extraction
  out <- x$stats$transmat[[sim]]
  out <- as.data.frame(out)
  class(out) <- c("transmat", class(out))
  return(out)
}


#' @title Extract Network Statistics from netsim or netdx Object
#'
#' @description Extracts a data frame of network statistics from a network
#'              epidemic model simulated with \code{netsim} or a network
#'              diagnostics object simulated with \code{netdx}.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netsim}} or
#'        \code{\link{netdx}}.
#' @param sim A vector of simulation numbers from the extracted object
#' @param network Network number, for \code{netsim} objects with multiple
#'        overlapping networks (advanced use, and not applicable to \code{netdx}
#'        objects).
#'
#' @keywords extract
#' @export
#'
#' @examples
#' # Two-group Bernoulli random graph TERGM
#' nw <- network_initialize(n = 100)
#' nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#' est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' dx <- netdx(est, nsim = 3, nsteps = 10, verbose = FALSE,
#'             nwstats.formula = ~edges + isolates)
#' get_nwstats(dx)
#' get_nwstats(dx, sim = 1)
#'
#' # SI epidemic model
#' param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.15)
#' init <- init.net(i.num = 10, i.num.g2 = 10)
#' control <- control.net(type = "SI", nsteps = 10, nsims = 3,
#'                        nwstats.formula = ~edges + meandeg + degree(0:5),
#'                        verbose = FALSE)
#' mod <- netsim(est, param, init, control)
#'
#' # Extract the network statistics from all or sets of simulations
#' get_nwstats(mod)
#' get_nwstats(mod, sim = 2)
#' get_nwstats(mod, sim = c(1, 3))
#'
#' # On the fly summary stats
#' summary(get_nwstats(mod))
#' colMeans(get_nwstats(mod))
#'
get_nwstats <- function(x, sim, network = 1) {

  ## Warnings and checks ##
  if (!(class(x) %in% c("netsim", "netdx"))) {
    stop("x must be of class netsim or netdx", call. = FALSE)
  }

  if (class(x) == "netsim") {
    nsims <- x$control$nsims
    nsteps <- x$control$nsteps
  } else {
    if (x$dynamic == TRUE) {
      nsims <- x$nsims
      nsteps <- x$nsteps
    } else {
      nsims <- 1
      nsteps <- x$nsims
    }
  }

  if (missing(sim)) {
    sim <- 1:nsims
  }
  if (max(sim) > nsims) {
    stop("Specify sims less than or equal to ", nsims, call. = FALSE)
  }

  if (class(x) == "netsim") {
    if (x$control$save.nwstats == FALSE || is.null(x$stats$nwstats)) {
      stop("Network statistics not saved in netsim object, check control.net
           settings", call. = FALSE)
    }
    if (network > x$control$num.nw) {
      stop("Specify network between 1 and ", x$control$num.nw, call. = FALSE)
    }
  }

  ## Extraction
  if (class(x) == "netsim") {
    out <- lapply(x$stats$nwstats, function(n) n[[network]])
    out <- out[sim]
  } else if (class(x) == "netdx") {
    out <- x$stats[sim]
  }

  out <- as.data.frame(do.call("rbind", out))
  out$time <- rep(1:min(nsteps, nrow(out)), length(sim))
  out$sim <- rep(sim, each = min(nsteps, nrow(out)))
  row.names(out) <- 1:nrow(out)
  out <- out[, c((ncol(out) - 1):ncol(out), 1:(ncol(out) - 2))]

  if (class(x) == "netdx" && x$dynamic == FALSE) {
    out <- out[, -2]
    names(out)[1] <- "sim"
  }

  return(out)
}


#' @title Extract Network Model Parameters
#'
#' @description Extracts a list of network model parameters saved in the
#'              initialization module.
#'
#' @param x Master data object used in \code{netsim} simulations.
#' @param network Network number, for simulations with multiple networks
#'        representing the population.
#'
#' @keywords extract internal
#' @export
#'
get_nwparam <- function(x, network = 1) {
  x$nwparam[[network]]
}


#' @title Extract Network Simulations
#'
#' @description Subsets the entire \code{netsim} object to a subset of
#'              simulations, essentially functioning like a reverse of
#'              \code{merge}.
#'
#' @param x An object of class \code{netsim}.
#' @param sims A numeric vector of simulation numbers to retain in the output
#'        object, or \code{"mean"} which selects the one simulation with the
#'        value of the variable specified in \code{var} closest to the mean of
#'        \code{var} across all simulations.
#' @param var A character vector of variables to retain from \code{x} if
#'        \code{sims} is a numeric vector, or a single variable name for
#'        selecting the average simulation from the set if \code{sims = "mean"}.
#'
#' @keywords extract
#' @export
#'
#' @examples
#' # Network model estimation
#' nw <- network_initialize(n = 100)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#' est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' # Epidemic model
#' param <- param.net(inf.prob = 0.3)
#' init <- init.net(i.num = 10)
#' control <- control.net(type = "SI", nsteps = 10, nsims = 3, verbose.int = 0)
#' mod1 <- netsim(est1, param, init, control)
#'
#' # Get sim 2
#' s.g2 <- get_sims(mod1, sims = 2)
#'
#' # Get sims 2 and 3 and keep only a subset of variables
#' s.g2.small <- get_sims(mod1, sims = 2:3, var = c("i.num", "si.flow"))
#'
#' # Extract the mean simulation for the variable i.num
#' sim.mean <- get_sims(mod1, sims = "mean", var = "i.num")
#'
get_sims <- function(x, sims, var) {

  if (class(x) != "netsim") {
    stop("x must be of class netsim", call. = FALSE)
  }

  nsims <- x$control$nsims

  if (missing(sims)) {
    stop("Specify sims as a vector of simulations or \"mean\" ", call. = FALSE)
  }
  if (length(sims) == 1 && sims ==
      "mean" && (missing(var) || length(var) > 1)) {
    stop("If sims == 'mean' then var must be a single varible name",
         call. = FALSE)
  }

  if (length(sims) == 1 && sims == "mean") {
    d <- tail(x$epi[[var]], 1)
    md <- mean(as.numeric(d))
    sims <- which.min(abs(d - md))
  }

  if (max(sims) > nsims) {
    stop("Maximum sims value for this object is ", nsims, call. = FALSE)
  }

  delsim <- setdiff(1:nsims, sims)
  out <- x
  if (length(delsim) > 0) {
    for (i in seq_along(out$epi)) {
      out$epi[[i]] <- out$epi[[i]][, -delsim, drop = FALSE]
    }
    if (!is.null(out$network)) {
      out$network[delsim] <- NULL
    }
    if (!is.null(out$stats$nwstats)) {
      out$stats$nwstats[delsim] <- NULL
    }
    if (!is.null(out$stats$transmat)) {
      out$stats$transmat[delsim] <- NULL
    }
    if (!is.null(out$control$save.other)) {
      oname <- out$control$save.other
      for (i in seq_along(oname)) {
        out[[oname[i]]][delsim] <- NULL
      }
    }
  }
  out$control$nsims <- length(sims)

  if (!missing(var)) {
    match.vars <- which(var %in% names(x$epi))
    out$epi <- out$epi[match.vars]
  }

  return(out)
}


#' @title Get Arguments from EpiModel Parameterization Functions
#'
#' @description Returns a list of argument names and values for use for
#'              parameter processing functions.
#'
#' @param formal.args The output of \code{formals(sys.function())}.
#' @param dot.args The output of \code{list(...)}.
#'
#' @export
#'
get_args <- function(formal.args, dot.args) {
  p <- list()
  formal.args[["..."]] <- NULL
  for (arg in names(formal.args)) {
    p[arg] <- list(get(arg, pos = parent.frame()))
  }

  names.dot.args <- names(dot.args)
  if (length(dot.args) > 0) {
    for (i in 1:length(dot.args)) {
      p[[names.dot.args[i]]] <- dot.args[[i]]
    }
  }
  return(p)
}

#' @title Extract the node records from Network Simulations
#'
#' @param sims An \code{EpiModel} object of class \code{netsim}.
#'
#' @return A list of \code{data.frame}s, one for each "measure" recorded in the
#' simulation by the `record_node_value` function.
#'
#' @details
#' See the "Time Varying Attributes" vignette
#'
#' @export
#'
get_node_records <- function(sims) {
  if (!any(class(sims) %in% "netsim")) {
    stop("`sims` must be of class netsim")
  }

  simnames <- names(sims[["node.records"]])

  dfs <- list()

  for (name in simnames) {
    records <- sims[["node.records"]][[name]]
    measures <- lapply(records, function(x) x[["measure"]])
    measure.names <- unique(measures)

    simnum <- as.numeric(sub("[^0-9]*", "", name))

    for (m in measure.names) {
      parts <- Filter(function(x) x[["measure"]] == m, records)
      parts <- lapply(parts, as.data.frame)
      d <- do.call("rbind", parts)
      d[["sim"]] <- simnum
      d <- d[, c("sim", "step", "measure", "uids", "values")]

      if (is.null(dfs[[m]])) {
        dfs[[m]] <- d
      } else {
        dfs[[m]] <- rbind(dfs[[m]], d)
      }
    }
  }

  return(dfs)
}
