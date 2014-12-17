
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
    if (class(x$fit) == "network") {
      nw <- simulate(x$formation,
                     basis = x$fit,
                     coef = x$coef.form.crude,
                     constraints = x$constraints)
    } else {
      nw <- simulate(x$fit)
    }
    modes <- ifelse(nw %n% "bipartite", 2, 1)
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


    # Network Parameters ------------------------------------------------------
    dat$nw <- nw

    dat$nwparam <- list(x[-which(names(x) == "fit")])

    dat$param$modes <- modes


    # Initialization ----------------------------------------------------------

    ## Infection Status and Time Modules
    dat <- init_status.net(dat)


    ## Initialize persistent IDs
    if (modes == 2 & param$vital == TRUE & control$delete.nodes == FALSE) {
      dat$nw <- init_pids(dat$nw, dat$control$pid.prefix)
    }


    ## Pull network val to attr
    form <- get_nwparam(dat)$formation
    fterms <- get_formula_terms(form)
    dat <- copy_toall_attr(dat, at = 1, fterms)


    ## Store current proportions of attr
    dat$temp$t1.tab <- get_attr_prop(dat$nw, fterms)


    ## Get initial prevalence
    dat <- get_prev.net(dat, at = 1)
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
  i.num.m2 <- dat$init$i.num.m2
  r.num <- dat$init$r.num
  r.num.m2 <- dat$init$r.num.m2

  status.vector <- dat$init$status.vector
  status.rand <- dat$init$status.rand
  num <- network.size(dat$nw)
  form <- get_nwparam(dat)$form
  statOnNw <- "status" %in% get_formula_terms(form)

  modes <- dat$param$modes
  if (modes == 1) {
    mode <- rep(1, num)
  } else {
    mode <- idmode(dat$nw)
  }
  nM1 <- sum(mode == 1)
  nM2 <- sum(mode == 2)

  type <- dat$control$type


  # Status ------------------------------------------------------------------

  ## Status passed on input network
  if (statOnNw == TRUE) {
    status <- get.vertex.attribute(dat$nw, "status")
  } else {
    if (!is.null(status.vector)) {
      status <- status.vector
    } else {
      ## Stochastic status
      if (status.rand == TRUE) {
        status <- rep(NA, num)
        if (type == "SIR") {
          status[which(mode == 1)] <- sample(
            x = c("s", "i", "r"),
            size = nM1,
            replace = TRUE,
            prob = c(1-(i.num/nM1)-(r.num/nM1), i.num/nM1, r.num/nM1))
          if (sum(status == "i" & mode == 1) == 0 & i.num > 0) {
            status[sample(which(mode == 1), size = i.num)] <- "i"
          }
          if (sum(status == "r" & mode == 1) == 0 & r.num > 0) {
            status[sample(which(mode == 1), size = r.num)] <- "r"
          }
          if (modes == 2) {
            status[which(mode == 2)] <- sample(
              x = c("s", "i", "r"),
              size = nM2,
              replace = TRUE,
              prob = c(1-(i.num.m2/nM2)-(r.num.m2/nM2), i.num.m2/nM2, r.num.m2/nM2))
            if (sum(status == "i" & mode == 2) == 0 & i.num.m2 > 0) {
              status[sample(which(mode == 2), size = i.num.m2)] <- "i"
            }
            if (sum(status == "r" & mode == 2) == 0 & r.num.m2 > 0) {
              status[sample(which(mode == 2), size = r.num.m2)] <- "r"
            }
          }
        } else {
          status[which(mode == 1)] <- sample(
            x = c("s", "i"),
            size = nM1,
            replace = TRUE,
            prob = c(1-(i.num/nM1), i.num/nM1))
          if (sum(status == "i" & mode == 1) == 0 & i.num > 0) {
            status[sample(which(mode == 1), size = i.num)] <- "i"
          }
          if (modes == 2) {
            status[which(mode == 2)] <- sample(
              x = c("s", "i"),
              size = nM2,
              replace = TRUE,
              prob = c(1-(i.num.m2/nM2), i.num.m2/nM2))
            if (sum(status == "i" & mode == 2) == 0 & i.num.m2 > 0) {
              status[sample(which(mode == 2), size = i.num.m2)] <- "i"
            }
          }
        }
      }

      ## Deterministic status
      if (status.rand == FALSE) {
        status <- rep("s", num)
        status[sample(which(mode == 1), size = i.num)] <- "i"
        if (modes == 2) {
          status[sample(which(mode == 2), size = i.num.m2)] <- "i"
        }
        if (type == "SIR") {
          status[sample(which(mode == 1 & status == "s"), size = r.num)] <- "r"
          if (modes == 2) {
            status[sample(which(mode == 2 & status == "s"), size = r.num.m2)] <- "r"
          }
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
  idsInf <- which(status == "i")
  infTime <- rep(NA, length(status))

  # If vital=TRUE, infTime is a uniform draw over the duration of infection
  if (dat$param$vital == TRUE && dat$param$di.rate > 0) {
    infTime[idsInf] <- -rgeom(n = length(idsInf), prob = dat$param$di.rate)+2
  } else {
    if (dat$control$type == "SI" || mean(dat$param$rec.rate) == 0) {
      # infTime a uniform draw over the number of sim time steps
      infTime[idsInf] <- ssample(1:(-dat$control$nsteps + 2),
                                     length(idsInf), replace = TRUE)
    } else {
      if (modes == 1) {
        infTime[idsInf] <- ssample(1:(-round(1/mean(dat$param$rec.rate)) + 2),
                                   length(idsInf), replace = TRUE)
      }
      if (modes == 2) {
        infM1 <- which(status == "i" & mode == 1)
        infTime[infM1] <- ssample(1:(-round(1/mean(dat$param$rec.rate)) + 2),
                                   length(infM1), replace = TRUE)
        infM2 <- which(status == "i" & mode == 2)
        infTime[infM2] <- ssample(1:(-round(1/mean(dat$param$rec.rate.m2)) + 2),
                                  length(infM2), replace = TRUE)

      }
    }
  }
  dat$attr$infTime <- infTime

  return(dat)
}


#' @title Persistent ID Initialization
#'
#' @description This function initializes the persistent IDs for
#'              a \code{networkDynamic} object.
#'
#' @param nw An object of class \code{networkDynamic}.
#' @param prefixes Character string prefix for mode-specific ID.
#'
#' @details
#' This function is used for \code{\link{netsim}} simulations
#' over bipartite networks for populations with vital dynamics. Persistent IDs are
#' required in this situation because when new nodes are added to the
#' first mode in a bipartite network, the IDs for the second mode shift
#' upward. Persistent IDs allow for an analysis of disease transmission
#' chains for these simulations. These IDs are also invoked in the
#' \code{\link{births.net}} module when the persistent IDs of incoming nodes
#' must be set.
#'
#' @export
#' @keywords netMod internal
#' @seealso \code{\link{initialize.pids}}
#'
#' @examples
#' # Initialize network with 25 female and 75 male
#' nw <- network.initialize(100, bipartite = 25)
#'
#' # Set persistent IDs using the default F/M prefix
#' nw <- init_pids(nw)
#' nw %v% "vertex.names"
#'
#' # Use another prefix combination
#' nw <- init_pids(nw, c("A", "B"))
#' nw %v% "vertex.names"
#'
init_pids <- function(nw, prefixes=c("F", "M")) {

  # Set persistent IDs
  t0.pids <- c(paste0(prefixes[1], 1:length(modeids(nw, 1))),
               paste0(prefixes[2], 1:length(modeids(nw, 2))))

  # Initialize persistent IDs on network
  nw <- set.network.attribute(nw, "vertex.pid", "vertex.names")
  nw <- set.vertex.attribute(nw, "vertex.names", t0.pids)

  return(nw)
}
