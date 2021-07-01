
#' @title Initializes EpiModel netsim Object for tergmLite Simulation
#'
#' @param dat A list object containing a \code{networkDynamic} object and other
#'        initialization information passed from \code{netsim}.
#'
#' @details
#' This function is typically used within the initialization modules of
#' \code{EpiModel} to establish the necessary infrastructure needed for 
#' \code{tergmLite} network resimulation.  The example below demonstrates 
#' the specific information returned.
#'
#' @export
#'
#' @return
#' Returns the list object \code{dat} and adds the element \code{el} which is
#' an edgelist representation of the network.  Also converts the \code{nw}
#' element to a \code{networkLite} representation.
#' 
#' @examples
#' \dontrun{
#' library("EpiModel")
#' nw <- network_initialize(100)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#' x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' param <- param.net(inf.prob = 0.3)
#' init <- init.net(i.num = 10)
#' control <- control.net(type = "SI", nsteps = 100, nsims = 5, tergmLite = TRUE)
#'
#' # networkLite representation after initialization
#' dat <- crosscheck.net(x, param, init, control)
#' dat <- initialize.net(x, param, init, control)
#' str(dat, max.level = 1)
#'
#' # Element added is el (edgelist representation of network)...
#' dat$el
#'
#' # ... and nw is now a networkLite
#' dat$nw[[1]]
#' }
#'
init_tergmLite <- function(dat) {

  num_nw <- ifelse(inherits(dat$nw, "network"), 1, length(dat$nw))

  dat$el <- list()
  dat$control$mcmc.control <- list()
  dat$control$nwstats.formulas <- list()
  
  for (i in 1:num_nw) {
    nwp <- dat$nwparam[[i]]
    is_tergm <- all(nwp$coef.diss$duration > 1)

    if(is(dat$nw[[i]], "networkDynamic")) {
      nw <- as.networkLite(network.collapse(dat$nw[[i]], at = 1))
    } else {
      nw <- as.networkLite(dat$nw[[i]])    
    }

    dat$el[[i]] <- as.edgelist(nw)
    attributes(dat$el[[i]])$vnames <- NULL

    nwstats_formula_name <- paste(c("nwstats.formula", if (num_nw > 1) i), collapse = ".")
    nwstats_formula <- NVL(dat$control[[nwstats_formula_name]], trim_env(~.))
    if (identical(nwstats_formula, "formation")) nwstats_formula <- nwp$formation
    dat$control$nwstats.formulas[[i]] <- nwstats_formula

    if (is_tergm) {
      mcmc_control_name <- paste(c("mcmc.control.tergm", if (num_nw > 1) i), collapse = ".")
      dat$control$mcmc.control[[i]] <- check.control.class("simulate.formula.tergm", "init_tergmLite", NVL(dat$control[[mcmc_control_name]], control.simulate.formula.tergm()))

      if (dat$control$tergmLite.track.duration == TRUE) {
        nw %n% "time" <- 1
        nw %n% "lasttoggle" <- cbind(as.edgelist(nw), 1)
      }
    } else {
      mcmc_control_name <- paste(c("mcmc.control.ergm", if (num_nw > 1) i), collapse = ".")
      dat$control$mcmc.control[[i]] <- check.control.class("simulate.formula", "init_tergmLite", NVL(dat$control[[mcmc_control_name]], control.simulate.formula()))
    }

    dat$nw[[i]] <- nw
  }
  
  return(dat)
}
