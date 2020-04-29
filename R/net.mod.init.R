
#' @title Initialization: netsim Module
#'
#' @description This function initializes the master \code{dat} object on which
#'              data are stored, simulates the initial state of the network, and
#'              simulates disease status and other attributes.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netest}}.
#' @param param An \code{EpiModel} object of class \code{\link{param.net}}.
#' @param init An \code{EpiModel} object of class \code{\link{init.net}}.
#' @param control An \code{EpiModel} object of class \code{\link{control.net}}.
#' @param s Simulation number, used for restarting dependent simulations.
#'
#' @export
#' @keywords internal
#'
initialize.net <- function(x, param, init, control, s) {

  if (control$start == 1) {

    # Master Data List --------------------------------------------------------
    dat <- list()
    dat$param <- param
    dat$init <- init
    dat$control <- control

    dat$attr <- list()
    dat$stats <- list()
    dat$temp <- list()


    # Initial Network Simulation ----------------------------------------------
    nw <- simulate(x$fit, basis = x$fit$newnetwork,
                   control = control$set.control.ergm)

    if (control$resimulate.network == TRUE) {
      if (class(x$fit) == "stergm") {
        nw <- network.collapse(nw, at = 1)
      }
      nw <- sim_nets(x, nw, nsteps = 1, control)
    }
    if (control$resimulate.network == FALSE) {
      nw <- sim_nets(x, nw, nsteps = control$nsteps, control)
    }
    nw <- activate.vertices(nw, onset = 1, terminus = Inf)
    dat$nw[[1]] <- nw

    # Network Parameters ------------------------------------------------------
    dat$nwparam <- list(x[-which(names(x) == "fit")])
    groups <- length(unique(get.vertex.attribute(nw, "group")))
    dat <- set_param(dat, "groups", groups)

    # Nodal Attributes --------------------------------------------------------

    # Standard attributes
    num <- network.size(nw)
    dat <- set_attr(dat, "active", rep(1, num), override.length.check = TRUE)
    dat <- set_attr(dat, "entrTime", rep(1, num))
    dat <- set_attr(dat, "exitTime", rep(NA, num))

    ## Pull attr on nw to dat$attr
    dat <- copy_nwattr_to_datattr(dat)

    ## Store current proportions of attr
    nwterms <- get_network_term_attr(nw)
    if (!is.null(nwterms)){
      dat$temp$nwterms <- nwterms
      dat$temp$t1.tab <- get_attr_prop(dat, nwterms)
    }

    # Conversions for tergmLite
    if (control$tergmLite == TRUE) {
      dat <- tergmLite::init_tergmLite(dat)
    }

    ## Infection Status and Time
    dat <- init_status.net(dat)


    # Summary Stats -----------------------------------------------------------
    dat <- do.call(control[["prevalence.FUN"]],list(dat, at = 1))


    # Restart/Reinit Simulations ----------------------------------------------
  } else if (control$start > 1) {
    dat <- list()

    dat$nw <- x$network[[s]]
    dat$param <- x$param
    dat$control <- control
    dat$nwparam <- x$nwparam
    dat$epi <- sapply(x$epi, function(var) var[s])
    names(dat$epi) <- names(x$epi)
    dat$attr <- x$attr[[s]]
    dat$stats <- sapply(x$stats, function(var) var[[s]])
    dat$temp <- list()
  }

  return(dat)
}


#' @title Disease Status Initialization Module for netsim
#'
#' @description This function sets the initial disease status on the
#'              network given the specified initial conditions.
#'
#' @param dat Master list object containing a \code{networkDynamic} object and other
#'        initialization information passed from \code{\link{netsim}}.
#'
#' @details
#' This internal function sets, either randomly or deterministically, the nodes
#' that are infected at the starting time of network simulations, \eqn{t_1}.
#' If the number to be initially infected is passed, this function may set the
#' initial number infected based on the number specified, either as a a set of
#' random draws from a binomial distribution or as the exact number specified. In
#' either case, the specific nodes infected are a random sample from the network.
#' In contrast, a set of specific nodes may be infected by passing the vector to
#' \code{\link{netsim}}.
#'
#' This module sets the time of infection for those nodes set infected
#' at the starting time of network simulations, \eqn{t_1}. For vital
#' dynamics models, the infection time for those nodes is a random draw from an
#' exponential distribution with the rate parameter defined by the \code{di.rate}
#' argument. For models without vital dynamics, the infection time is a random
#' draw from a uniform distribution of integers with a minimum of 1 and a maximum
#' of the number of time steps in the model. In both cases, to set the infection
#' times to be in the past, these times are multiplied by -1, and 2 is added to
#' allow for possible infection times up until step 2, when the disease simulation
#' time loop starts.
#'
#' @seealso This is an initialization module for \code{\link{netsim}}.
#'
#' @export
#' @keywords netMod internal
#'
init_status.net <- function(dat) {

  type <- get_control(dat, "type", override.null.error = TRUE)
  nsteps <- get_control(dat, "nsteps")
  tergmLite <- get_control(dat, "tergmLite")
  vital <- get_param(dat, "vital")
  groups <- get_param(dat, "groups")
  status.vector <- get_init(dat, "status.vector", override.null.error = TRUE)
  if (type %in% c("SIS", "SIR")){
    rec.rate <- get_param(dat, "rec.rate")
  }
  if (vital == TRUE) {
    di.rate <- get_param(dat, "di.rate")
  }

  # Variables ---------------------------------------------------------------
  i.num <- get_init(dat, "i.num", override.null.error = TRUE)
  if (type  == "SIR" & is.null(status.vector)){
    r.num <- get_init(dat, "r.num")
  }

  num <- sum(get_attr(dat, "active") == 1)

  if (groups == 2) {
    group <- get_attr(dat, "group")
    i.num.g2 <- get_init(dat, "i.num.g2")
    if (type  == "SIR" & is.null(status.vector)) {
      r.num.g2 <- get_init(dat, "r.num.g2", override.null.error = TRUE)
    }
  } else {
    group <- rep(1, num)
  }

  statOnNw <- "status" %in% dat$temp$nwterms

  # Status ------------------------------------------------------------------

  ## Status passed on input network
  if (statOnNw == FALSE) {
    if (!is.null(status.vector)) {
      status <- status.vector
    } else {
      status <- rep("s", num)
      status[sample(which(group == 1), size = i.num)] <- "i"
      if (groups == 2) {
        status[sample(which(group == 2), size = i.num.g2)] <- "i"
      }
      if (type == "SIR"  && !is.null(type)) {
        status[sample(which(group == 1 & status == "s"), size = r.num)] <- "r"
        if (groups == 2) {
          status[sample(which(group == 2 & status == "s"), size = r.num.g2)] <- "r"
        }
      }
    }
    dat <- set_attr(dat, "status", status)
  } else {
    status <- get.vertex.attribute(dat$nw[[1]], "status")
    dat <- set_attr(dat, "status", status)
  }


  ## Set up TEA status
  if (tergmLite == FALSE) {
    if (statOnNw == FALSE) {
      nw <- set.vertex.attribute(dat$nw[[1]], "status", status)
      dat <- set_network(dat, nw, 1)
    }
    nw <- activate.vertex.attribute(dat$nw[[1]],
                                    prefix = "testatus",
                                    value = status,
                                    onset = 1,
                                    terminus = Inf)
    dat <- set_network(dat, nw, 1)
  }


  # Infection Time ----------------------------------------------------------
  ## Set up inf.time vector
  if(!is.null(type)){
    idsInf <- which(status == "i")
    infTime <- rep(NA, length(status))
    infTime.vector <- get_init(dat, "infTime.vector", override.null.error = TRUE)

    if (!is.null(infTime.vector)) {
      infTime <- infTime.vector
    } else {
      # If vital dynamics, infTime is a geometric draw over the duration of infection
      if (vital == TRUE && di.rate > 0) {
        if (type == "SI") {
          infTime[idsInf] <- -rgeom(n = length(idsInf), prob = di.rate) + 2
        } else {
          infTime[idsInf] <- -rgeom(n = length(idsInf),
                                    prob = di.rate +
                                      (1 - di.rate)*mean(rec.rate)) + 2
        }
      } else {
        if (type == "SI" || mean(rec.rate) == 0) {
          # if no recovery, infTime a uniform draw over the number of sim time steps
          infTime[idsInf] <- ssample(1:(-nsteps + 2),
                                     length(idsInf), replace = TRUE)
        } else {
          infTime[idsInf] <- -rgeom(n = length(idsInf), prob = mean(rec.rate)) + 2
        }
      }
    }

    dat <- set_attr(dat, "infTime", infTime)
  }

  return(dat)
}
