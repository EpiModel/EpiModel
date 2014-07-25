
#' @title Initialization: netsim Module
#'
#' @description This function initializes the master \code{all} object on which
#'              data are stored, simulates the initial state of the network, and
#'              simulates disease status and other attributes.
#'
#' @param x an \code{EpiModel} object of class \code{\link{netest}}.
#' @param param an \code{EpiModel} object of class \code{\link{param.net}}.
#' @param init an \code{EpiModel} object of class \code{\link{init.net}}.
#' @param control an \code{EpiModel} object of class \code{\link{control.net}}.
#'
#' @export
#' @keywords internal
#'
initialize.net <- function(x, param, init, control) {

  # Master Data List --------------------------------------------------------
  all <- list()
  all$param <- param
  all$init <- init
  all$control <- control

  all$attr <- list()
  all$stats <- list()
  all$nwparam <- list()
  all$temp <- list()


  # Network Simulation ------------------------------------------------------
  nw <- simulate(x$fit)
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
  all$nw <- nw
  all$nwparam$formation <- x$formation
  all$nwparam$dissolution <- x$dissolution
  all$nwparam$coef.form <- x$coef.form
  all$nwparam$coef.diss <- x$coef.diss
  all$nwparam$constraints <-  x$constraints
  all$nwparam$target.stats <- x$target.stats
  if (is.null(x$constraints)) {
    all$nwparam$constraints <- ~ .
  } else {
    all$nwparam$constraints <- x$constraints
  }
  all$param$modes <- modes


  # Initialization ----------------------------------------------------------

  ## Infection Status and Time Modules
  all <- init_status.net(all)


  ## Initialize persistent IDs
  if (modes == 2 & param$vital == TRUE & control$delete.nodes == FALSE) {
    all$nw <- init_pids(all$nw, all$control$pid.prefix)
  }


  ## Pull network val to attr
  all <- copy_toall_attr(all, at = 1)


  ## Store current proportions of attr
  t <- get_formula_terms(all$nwparam$formation)
  all$temp$t1.tab <- get_attr_prop(all$nw, t)


  ## Get initial prevalence
  all <- get_prev.net(all, at = 1)


  return(all)
}


#' @title Disease Status Initialization Module for netsim
#'
#' @description This function sets the initial disease status on the
#'              network given the specified initial conditions.
#'
#' @param all a list object containing a \code{networkDynamic} object and other
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
init_status.net <- function(all) {

  # Variables ---------------------------------------------------------------
  tea.status <- all$control$tea.status
  i.num <- all$init$i.num
  i.num.m2 <- all$init$i.num.m2
  r.num <- all$init$r.num
  r.num.m2 <- all$init$r.num.m2

  status.vector <- all$init$status.vector
  status.rand <- all$init$status.rand
  num <- network.size(all$nw)
  statOnNw <- "status" %in% get_formula_terms(all$nwparam$formation)

  modes <- all$param$modes
  if (modes == 1) {
    mode <- rep(1, num)
  } else {
    mode <- idmode(all$nw)
  }
  nM1 <- sum(mode == 1)
  nM2 <- sum(mode == 2)

  type <- all$control$type



  # Status ------------------------------------------------------------------

  ## Status passed on input network
  if (statOnNw == TRUE) {
    status <- get.vertex.attribute(all$nw, "status")
  } else {
    if (!is.null(status.vector)) {
      status <- status.vector
    } else {
      ## Stochastic status
      if (status.rand == TRUE) {
        status <- rep(NA, num)
        if (type == "SIR") {
          status[which(mode == 1)] <- sample(
            x = c(0:2),
            size = nM1,
            replace = TRUE,
            prob = c(1-(i.num/nM1)-(r.num/nM1), i.num/nM1, r.num/nM1))
          if (sum(status == 1 & mode == 1) == 0 & i.num > 0) {
            status[sample(which(mode == 1), size = i.num)] <- 1
          }
          if (sum(status == 2 & mode == 1) == 0 & r.num > 0) {
            status[sample(which(mode == 1), size = r.num)] <- 2
          }
          if (modes == 2) {
            status[which(mode == 2)] <- sample(
              x = c(0:2),
              size = nM2,
              replace = TRUE,
              prob = c(1-(i.num.m2/nM2)-(r.num.m2/nM2), i.num.m2/nM2, r.num.m2/nM2))
            if (sum(status == 1 & mode == 2) == 0 & i.num.m2 > 0) {
              status[sample(which(mode == 2), size = i.num.m2)] <- 1
            }
            if (sum(status == 2 & mode == 2) == 0 & r.num.m2 > 0) {
              status[sample(which(mode == 2), size = r.num.m2)] <- 2
            }
          }
        } else {
          status[which(mode == 1)] <- sample(
            x = c(0:1),
            size = nM1,
            replace = TRUE,
            prob = c(1-(i.num/nM1), i.num/nM1))
          if (sum(status == 1 & mode == 1) == 0 & i.num > 0) {
            status[sample(which(mode == 1), size = i.num)] <- 1
          }
          if (modes == 2) {
            status[which(mode == 2)] <- sample(
              x = c(0:1),
              size = nM2,
              replace = TRUE,
              prob = c(1-(i.num.m2/nM2), i.num.m2/nM2))
            if (sum(status == 1 & mode == 2) == 0 & i.num.m2 > 0) {
              status[sample(which(mode == 2), size = i.num.m2)] <- 1
            }
          }
        }
      }

      ## Deterministic status
      if (status.rand == FALSE) {
        status <- rep(0, num)
        status[sample(which(mode == 1), size = i.num)] <- 1
        if (modes == 2) {
          status[sample(which(mode == 2), size = i.num.m2)] <- 1
        }
        if (type == "SIR") {
          status[sample(which(mode == 1 & status == 0), size = r.num)] <- 2
          if (modes == 2) {
            status[sample(which(mode == 2 & status == 0), size = r.num.m2)] <- 2
          }
        }
      }


    }
  }
  all$attr$status <- status

  ## Save out other attr
  all$attr$active <- rep(1, length(status))
  if (tea.status == TRUE) {
    all$nw <- activate.vertex.attribute(all$nw,
                                        prefix = "testatus",
                                        value = status,
                                        onset = 1,
                                        terminus = Inf)
  }


  # Infection Time ----------------------------------------------------------
  ## Set up inf.time vector
  idsInf <- which(status == 1)
  infTime <- rep(NA, length(status))

  # If vital=TRUE, infTime is a uniform draw over the duration of infection
  if (all$param$vital == TRUE && all$param$di.rate > 0) {
    infTime[idsInf] <- -rgeom(n = length(idsInf),
                              prob = all$param$di.rate)+2
  } else {
    if (all$control$type == "SI" || all$param$rec.rate == 0) {
      # infTime a uniform draw over the number of sim time steps
      infTime[idsInf] <- ssample(1:(-all$control$nsteps + 2),
                                     length(idsInf),
                                     replace = TRUE)
    } else {
      infTime[idsInf] <- ssample(1:(-round(1/all$param$rec.rate) + 2),
                                 length(idsInf),
                                 replace = TRUE)
      #TODO: divide this by mode if rec.rate != rec.rate.m2
    }
  }
  all$attr$infTime <- infTime

  return(all)
}



#' @title Persistent ID Initialization
#'
#' @description This function initializes the persistent IDs for
#'              a \code{networkDynamic} object.
#'
#' @param nw an object of class \code{networkDynamic}.
#' @param prefixes character string prefix for mode-specific ID.
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
