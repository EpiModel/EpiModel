
#' @title Extract networkDynamic Object from Network Epidemic Model
#'
#' @description Extracts the networkDynamic object from a network epidemic model
#'              simulated with \code{netsim}, with the option to collapse the
#'              extracted network at a specific time step.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netsim}}.
#' @param sim Simulation number of extracted network.
#' @param network Network number, for simulations with multiple networks
#'        representing the population.
#' @param collapse If \code{TRUE}, collapse the \code{networkDynamic} object to
#'        a static \code{network} object at a specified time step.
#' @param at If \code{collapse} is used, the time step at which the extracted
#'        network should be collapsed.
#'
#' @keywords extract
#' @export
#'
#' @examples
#' \dontrun{
#' ## Simulate SI epidemic on bipartite Bernoulli random graph
#' nw <- network.initialize(n = 100, bipartite = 50, directed = FALSE)
#' formation <- ~ edges
#' target.stats <- 50
#' dissolution <- ~ offset(edges)
#' duration <- 20
#' coef.diss <- dissolution_coefs(dissolution, duration)
#' est <- netest(nw,
#'                formation,
#'                dissolution,
#'                target.stats,
#'                coef.diss,
#'                verbose = FALSE)
#' param <- param.net(inf.prob = 0.3, inf.prob.m2 = 0.15)
#' init <- init.net(i.num = 10, i.num.m2 = 10)
#' control <- control.net(type = "SI", nsteps = 10, nsims = 3,
#'                        verbose = FALSE)
#' mod <- netsim(est, param, init, control)
#'
#' ## Extract the network from simulation 2
#' get_network(mod, sim = 2)
#'
#' ## Extract and collapse the network from simulation 1
#' get_network(mod, collapse = TRUE, at = 5)
#' }
#'
get_network <- function(x, sim = 1, network = 1, collapse = FALSE, at) {

  ## Warnings and checks
  if (class(x) != "netsim") {
    stop("x must be of class netsim", call. = FALSE)
  }

  if (sim > x$control$nsims) {
    stop("Specify sim between 1 and ", x$control$nsims, call. = FALSE)
  }

  if (x$control$save.network == FALSE || is.null(x$network)) {
    stop("Network object not saved in netsim object, check control.net settings", call. = FALSE)
  }

  if (network > x$control$num.nw) {
    stop("Specify network between 1 and ", x$control$num.nw, call. = FALSE)
  }

  if (collapse == TRUE && (missing(at) || at > x$control$nsteps)) {
    stop("Specify collapse time step between 1 and ", x$control$nsteps, call. = FALSE)
  }

  ## Extraction
  if (x$control$num.nw == 1) {
    out <- x$network[[sim]]
  } else {
    out <- x$network[[sim]][[network]]
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
#' A data frame with the following collumns
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
#' \dontrun{
#' ## Simulate SI epidemic on bipartite Bernoulli random graph
#' nw <- network.initialize(n = 100, bipartite = 50, directed = FALSE)
#' formation <- ~ edges
#' target.stats <- 50
#' dissolution <- ~ offset(edges)
#' duration <- 20
#' coef.diss <- dissolution_coefs(dissolution, duration)
#' est <- netest(nw,
#'                formation,
#'                dissolution,
#'                target.stats,
#'                coef.diss,
#'                verbose = FALSE)
#' param <- param.net(inf.prob = 0.3, inf.prob.m2 = 0.15)
#' init <- init.net(i.num = 10, i.num.m2 = 10)
#' control <- control.net(type = "SI", nsteps = 10, nsims = 3,
#'                        verbose = FALSE)
#' mod <- netsim(est, param, init, control)
#'
#' ## Extract the transmission matrix from simulation 2
#' get_transmat(mod, sim = 2)
#' }
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
    stop("transmat not saved in netsim object, check control.net settings", call. = FALSE)
  }


  ## Extraction
  out <- x$stats$transmat[[sim]]
  out <- as.data.frame(out)

  return(out)
}


#' @title Extract Network Statistics from Network Epidemic Model
#'
#' @description Extracts a data frame of network statistics from a network
#'              epidemic model.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netsim}}.
#' @param sim A vector of simulation numbers of extracted network.
#' @param network Network number, for simulations with multiple networks
#'        representing the population.
#'
#' @keywords extract
#' @export
#'
#' @examples
#' \dontrun{
#' ## Simulate SI epidemic on bipartite Bernoulli random graph
#' nw <- network.initialize(n = 100, bipartite = 50, directed = FALSE)
#' formation <- ~ edges
#' target.stats <- 50
#' dissolution <- ~ offset(edges)
#' duration <- 20
#' coef.diss <- dissolution_coefs(dissolution, duration)
#' est <- netest(nw,
#'                formation,
#'                dissolution,
#'                target.stats,
#'                coef.diss,
#'                verbose = FALSE)
#' param <- param.net(inf.prob = 0.3, inf.prob.m2 = 0.15)
#' init <- init.net(i.num = 10, i.num.m2 = 10)
#' control <- control.net(type = "SI", nsteps = 10, nsims = 3,
#'                        nwstats.formula = ~ edges + meandeg + degree(0:5),
#'                        verbose = FALSE)
#' mod <- netsim(est, param, init, control)
#'
#' ## Extract the network statistics from simulation 2
#' get_nwstats(mod)
#' get_nwstats(mod, sim = c(1,3))
#' }
#'
get_nwstats <- function(x, sim, network = 1) {

  ## Warnings and checks
  if (class(x) != "netsim") {
    stop("x must be of class netsim", call. = FALSE)
  }

  if (missing(sim)) {
    sim <- 1:x$control$nsims
  }
  if (max(sim) > x$control$nsims) {
    stop("Specify sims less than or equal to ", x$control$nsims, call. = FALSE)
  }

  if (x$control$save.nwstats == FALSE || is.null(x$stats$nwstats)) {
    stop("Network statistics not saved in netsim object, check control.net settings", call. = FALSE)
  }

  if (network > x$control$num.nw) {
    stop("Specify network between 1 and ", x$control$num.nw, call. = FALSE)
  }

  ## Extraction
  if (x$control$num.nw == 1) {
    out <- x$stats$nwstats[sim]
  } else {
    out <- lapply(x$stats$nwstats, function(n) n[[network]])
    out <- out[sim]
  }

  out <- lapply(out, as.data.frame)

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

  out <- x$nwparam[[network]]
  return(out)
}


#' @title Extract Network Simulations
#'
#' @description Subsets the entire \code{netsim} object to a subset of
#'              simulations, essentially functioning like a reverse of
#'              \code{merge}.
#'
#' @param x An object of class \code{netsim}.
#' @param sims A vector of simulation numbers to retain in the output object,
#'        or \code{"mean"} which automatically selects the one simulation with
#'        the number infected at the final time step closest to the mean across
#'        all simulations.
#'
#' @keywords extract
#' @export
#'
get_sims <- function(x, sims) {

  nsims <- x$control$nsims

  if (missing(sims)) {
    stop("Specify sims as a vector of simulations or \"mean\" ", call. = FALSE)
  }

  if (length(sims) == 1 && sims == "mean") {
    d <- tail(x$epi$i.num, 1)
    md <- mean(as.numeric(d))
    sims <- which.min(abs(d - md))
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
  return(out)
}
