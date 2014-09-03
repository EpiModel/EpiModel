
## Get functions


#' @title Extract networkDynamic Object from Network Epidemic Model
#'
#' @description Extracts the networkDynamic object from a network epidemic model
#'              simulated with \code{netsim}, with the option to collapse the
#'              extracted network at a specific time step.
#'
#' @param x an \code{EpiModel} object of class \code{\link{netsim}}.
#' @param sim simulation number of extracted network.
#' @param collapse if \code{TRUE}, collapse the \code{networkDynamic} object to
#'        a static \code{network} object at a specified time step.
#' @param at if \code{collapse} is used, the time step at which the extracted
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
get_network <- function(x, sim = 1, collapse = FALSE, at) {

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

  if (collapse == TRUE && (missing(at) || at > x$control$nsteps)) {
    stop("Specify collapse time step between 1 and ", x$control$nsteps, call. = FALSE)
  }

  ## Extraction
  out <- x$network[[sim]]


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
#' @param x an \code{EpiModel} object of class \code{\link{netsim}}.
#' @param sim simulation number of extracted network.
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

  return(out)
}


#' @title Extract Network Statistics from Network Epidemic Model
#'
#' @description Extracts a data frame of network statistics from a network
#'              epidemic model.
#'
#' @param x an \code{EpiModel} object of class \code{\link{netsim}}.
#' @param sim simulation number of extracted network.
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
#' get_nwstats(mod, sim = 2)
#' }
#'
get_nwstats <- function(x, sim = 1) {

  ## Warnings and checks
  if (class(x) != "netsim") {
    stop("x must be of class netsim", call. = FALSE)
  }

  if (sim > x$control$nsims) {
    stop("Specify sim between 1 and ", x$control$nsims, call. = FALSE)
  }

  if (x$control$save.nwstats == FALSE || is.null(x$stats$nwstats)) {
    stop("Network statistics not saved in netsim object, check control.net settings", call. = FALSE)
  }


  ## Extraction
  out <- x$stats$nwstats[[sim]]

  return(out)
}
