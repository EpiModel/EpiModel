
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


    # Network Simulation ------------------------------------------------------

    nw <- simulate(x$fit, basis = x$fit$newnetwork,
                   control = control$set.control.ergm)

    if (control$depend == TRUE) {
      if (class(x$fit) == "stergm") {
        nw <- network.collapse(nw, at = 1)
      }
      nw <- sim_nets(x, nw, nsteps = 1, control)
    }
    if (control$depend == FALSE) {
      nw <- sim_nets(x, nw, nsteps = control$nsteps, control)
    }
    nw <- activate.vertices(nw, onset = 1, terminus = Inf)
    dat$nw <- nw


    # Network Parameters ------------------------------------------------------
    dat$nwparam <- list(x[-which(names(x) == "fit")])
    groups <- length(unique(get.vertex.attribute(nw, "group")))
    dat$param$groups <- groups

    # Convert to tergmLite method
    if (control$tgl == TRUE) {
      num <- network.size(nw)
      if (groups == 1){
      dat$attr$group <- rep(1, num)
      } else {
        dat$attr$group <- groupids(dat$nw)
      }
      nw <- tergmLite::init_tergmLite(dat)
      dat$nw <- nw
      dat$param$num <- num
    }


    # Initialization ----------------------------------------------------------

    ## Initialize persistent IDs
    if (control$tgl == FALSE) {
      if (control$use.pids == TRUE) {
        dat$nw <- init_pids(dat$nw, dat$param$groups, dat$control$pid.prefix)
      }
    }

    ## Pull network val to attr
    form <- get_nwparam(dat)$formation
    fterms <- get_formula_term_attr(form, nw)
    dat <- copy_toall_attr(dat, at = 1, fterms)

    ## Store current proportions of attr
    dat$temp$fterms <- fterms
    dat$temp$t1.tab <- get_attr_prop(dat$nw, fterms)


    ## Infection Status and Time Modules
    if (control$tgl == FALSE) {
    dat <- init_status.net(dat)
    } else{
      dat <- init_status.net.tgl(dat)
    }

    ## Get initial prevalence
    dat <- do.call(control[["get_prev.FUN"]],list(dat, at = 1))
  } else {
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

  # Variables ---------------------------------------------------------------
  tea.status <- dat$control$tea.status
  i.num <- dat$init$i.num
  i.num.g2 <- dat$init$i.num.g2
  r.num <- dat$init$r.num
  r.num.g2 <- dat$init$r.num.g2

  status.vector <- dat$init$status.vector
  num <- network.size(dat$nw)
  statOnNw <- "status" %in% dat$temp$fterms

  groups <- dat$param$groups
  if (groups == 1) {
    group <- rep(1, num)
  } else {
    group <- idgroup(dat$nw)
  }

  type <- dat$control$type


  # Status ------------------------------------------------------------------

  ## Status passed on input network
  if (statOnNw == TRUE) {
    status <- get.vertex.attribute(dat$nw, "status")
  } else {
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
  }
  dat$attr$status <- status


  ## Save out other attr
  dat$attr$active <- rep(1, length(status))
  dat$attr$entrTime <- rep(1, length(status))
  dat$attr$exitTime <- rep(NA, length(status))
  if (tea.status == TRUE) {
    dat$nw <- activate.vertex.attribute(dat$nw,
                                        prefix = "testatus",
                                        value = status,
                                        onset = 1,
                                        terminus = Inf)
  }


  # Infection Time ----------------------------------------------------------
  ## Set up inf.time vector
  if(!is.null(type)){
    idsInf <- which(status == "i")
    infTime <- rep(NA, length(status))

    if (!is.null(dat$init$infTime.vector)) {
      infTime <- dat$init$infTime.vector
    } else {
      # If vital dynamics, infTime is a geometric draw over the duration of infection
      if (dat$param$vital == TRUE && dat$param$di.rate > 0) {
        if (dat$control$type == "SI") {
          infTime[idsInf] <- -rgeom(n = length(idsInf), prob = dat$param$di.rate) + 2
        } else {
          infTime[idsInf] <- -rgeom(n = length(idsInf),
                                    prob = dat$param$di.rate +
                                      (1 - dat$param$di.rate)*mean(dat$param$rec.rate)) + 2
        }
      } else {
        if (dat$control$type == "SI" || mean(dat$param$rec.rate) == 0) {
          # if no recovery, infTime a uniform draw over the number of sim time steps
          infTime[idsInf] <- ssample(1:(-dat$control$nsteps + 2),
                                     length(idsInf), replace = TRUE)
        } else {
          infTime[idsInf] <- -rgeom(n = length(idsInf), prob = mean(dat$param$rec.rate)) + 2
        }
      }
    }

    dat$attr$infTime <- infTime
  }

  return(dat)
}


#' @title TergmLite: Disease Status Initialization Module for netsim
#'
#' @description This function sets the initial disease status on the
#'              network given the specified initial conditions.
#'
#' @param dat Master list object containing a \code{networkDynamic} object and other
#'        initialization information passed from \code{\link{netsim}}.
#'
#' @details
#' This internal function sets randomly, the nodes that are infected at the
#' starting time of network simulations, \eqn{t_1}. If the number to be initially
#' infected is passed, this function sets the initial number infected based
#' on the number specified as a a set of random draws from a binomial
#' distribution or as the exact number specified.
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
init_status.net.tgl <- function(dat) {

  # Variables ---------------------------------------------------------------
  i.num <- dat$init$i.num
  i.num.g2 <- dat$init$i.num.g2
  r.num <- dat$init$r.num
  r.num.g2 <- dat$init$r.num.g2

  active <- dat$attr$active
  num <- dat$param$num
  groups <- dat$param$groups
  group <- dat$attr$group

  type <- dat$control$type

  # Status ---------------------------------------------------------------
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

  dat$attr$status <- status

  ## Save out other attr
  dat$attr$active <- rep(1, length(status))
  dat$attr$entrTime <- rep(1, length(status))
  dat$attr$exitTime <- rep(NA, length(status))


  # Infection Time ----------------------------------------------------------
  ## Set up inf.time vector
  if(!is.null(type)){
    idsInf <- which(status == "i")
    infTime <- rep(NA, length(status))

    if (!is.null(dat$init$infTime.vector)) {
      infTime <- dat$init$infTime.vector
    } else {
      # If vital dynamics, infTime is a geometric draw over the duration of infection
      if (dat$param$vital == TRUE && dat$param$di.rate > 0) {
        if (dat$control$type == "SI") {
          infTime[idsInf] <- -rgeom(n = length(idsInf), prob = dat$param$di.rate) + 2
        } else {
          infTime[idsInf] <- -rgeom(n = length(idsInf),
                                    prob = dat$param$di.rate +
                                      (1 - dat$param$di.rate)*mean(dat$param$rec.rate)) + 2
        }
      } else {
        if (dat$control$type == "SI" || mean(dat$param$rec.rate) == 0) {
          # if no recovery, infTime a uniform draw over the number of sim time steps
          infTime[idsInf] <- ssample(1:(-dat$control$nsteps + 2),
                                     length(idsInf), replace = TRUE)
        } else {
          infTime[idsInf] <- -rgeom(n = length(idsInf), prob = mean(dat$param$rec.rate)) + 2
        }
      }
    }

    dat$attr$infTime <- infTime
  }

  return(dat)

}


#' @title Persistent ID Initialization
#'
#' @description This function initializes the persistent IDs for
#'              a \code{networkDynamic} object.
#'
#' @param nw An object of class \code{networkDynamic}.
#' @param prefixes Character string prefix for group-specific ID.
#'
#' @details
#' This function is used for \code{\link{netsim}} simulations
#' over two-group networks for populations with vital dynamics. Persistent IDs are
#' required in this situation because when new nodes are added to the
#' first group in a two-group network, the IDs for the second mode shift
#' upward. Persistent IDs allow for an analysis of disease transmission
#' chains for these simulations. These IDs are also invoked in the
#' \code{\link{arrivals.net}} module when the persistent IDs of incoming nodes
#' must be set.
#'
#' @export
#' @keywords netMod internal
#' @seealso \code{\link{initialize.pids}}
#'
#' @examples
#' # Initialize network with 25 female and 75 male
#' nw <- network.initialize(100)
#' nw <- set.vertex.attribute(nw, "male", rep(c(25, 75), c(25,75)))
#'
#' # Set persistent IDs using the default F/M prefix
#' nw <- init_pids(nw)
#' nw %v% "vertex.names"
#'
#' # Use another prefix combination
#' nw <- init_pids(nw, c("A", "B"))
#' nw %v% "vertex.names"
#'
init_pids <- function(nw, groups, prefixes=c("F", "M")) {

  if (is.null(nw$gal$vertex.pid)) {
    if (groups == 1) {
      nw <- initialize.pids(nw)
    } else {
      t0.pids <- c(paste0(prefixes[1], 1:length(groupids(nw, 1))),
                   paste0(prefixes[2], 1:length(groupids(nw, 2))))

      nw <- set.network.attribute(nw, "vertex.pid", "vertex.names")
      nw <- set.vertex.attribute(nw, "vertex.names", t0.pids)
    }
  }

  return(nw)
}
