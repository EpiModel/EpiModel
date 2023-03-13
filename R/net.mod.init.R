
#' @title Initialization: netsim Module
#'
#' @description This function initializes the main \code{dat} object on which
#'              data are stored, simulates the initial state of the networks,
#'              and simulates disease status and other attributes.
#'
#' @param x If \code{control$start == 1}, either a fitted network model object
#'        of class \code{netest} or a list of such objects. If
#'        \code{control$start > 1}, an object of class \code{netsim}. When
#'        multiple networks are used, the node sets (including network size
#'        and nodal attributes) are assumed to be the same for all networks.
#' @param param An \code{EpiModel} object of class \code{\link{param.net}}.
#' @param init An \code{EpiModel} object of class \code{\link{init.net}}.
#' @param control An \code{EpiModel} object of class \code{\link{control.net}}.
#' @param s Simulation number, used for restarting dependent simulations.
#' @details When re-initializing a simulation, the \code{netsim} object passed
#'          to \code{initialize.net} must contain the elements \code{param},
#'          \code{nwparam}, \code{epi}, \code{attr}, \code{temp},
#'          \code{coef.form}, and \code{num.nw}. If \code{tergmLite == TRUE} it
#'          must also contain the elements \code{el} and \code{net_attr}. If
#'          \code{tergmLite == FALSE} it must also contain the element
#'          \code{network}.
#'
#' @return A \code{dat} main list object.
#'
#' @export
#' @keywords internal
#'
initialize.net <- function(x, param, init, control, s) {

  if (control$start == 1) {

    # Main Data List --------------------------------------------------------
    dat <- create_dat_object(param, init, control)

    # network and stats initialization
    dat <- init_nets(dat, x)

    ## Store current proportions of attr
    if (!is.null(dat$temp$nwterms)) {
      dat$temp$t1.tab <- get_attr_prop(dat, dat$temp$nwterms)
    }

    # simulate first time step
    dat <- sim_nets_t1(dat)
    dat <- summary_nets(dat, at = 1L)

    ## Infection Status and Time
    dat <- init_status.net(dat)

    # Summary Stats -----------------------------------------------------------
    dat <- do.call(control[["prevalence.FUN"]], list(dat, at = 1))

    # Restart/Reinit Simulations ----------------------------------------------
  } else if (control$start > 1) {
    ## check that required names are present
    required_names <- c(
      "param",
      "nwparam",
      "epi",
      "attr",
      "temp",
      "coef.form",
      "num.nw",
      if (control[["tergmLite"]] == TRUE) c("el", "net_attr"),
      if (control[["tergmLite"]] == FALSE) "network"
    )
    missing_names <- setdiff(required_names, names(x))
    if (length(missing_names) > 0) {
      stop("x is missing the following elements required for re-initialization: ",
           paste.and(missing_names), call. = FALSE)
    }

    dat <- create_dat_object(param = x$param, control = control)

    dat$num.nw <- x$num.nw
    if (control[["tergmLite"]] == TRUE) {
      dat$el <- x$el[[s]]
      dat$net_attr <- x$net_attr[[s]]
    }
    if (control[["tergmLite"]] == FALSE) {
      dat$nw <- x$network[[s]]
    }

    # copy if present
    if (length(x[["el.cuml"]]) >= s) {
      dat[["el.cuml"]] <- x[["el.cuml"]][[s]]
    }
    if (length(x[["_last_unique_id"]]) >= s) {
      dat[["_last_unique_id"]] <- x[["_last_unique_id"]][[s]]
    }

    dat$nwparam <- x$nwparam
    for (network in seq_len(dat$num.nw)) {
      dat$nwparam[[network]]$coef.form <- x$coef.form[[s]][[network]]
    }
    dat$epi <- sapply(x$epi, function(var) var[s])
    names(dat$epi) <- names(x$epi)
    dat$attr <- x$attr[[s]]
    dat$temp <- x$temp[[s]]

    dat$stats <- lapply(x$stats, function(var) var[[s]])
    if (get_control(dat, "save.nwstats") == TRUE) {
      dat$stats$nwstats <- lapply(dat$stats$nwstats,
                                  function(oldstats) {
                                    padded_vector(list(oldstats),
                                                  get_control(dat, "nsteps") -
                                                    get_control(dat, "start") + 2L)
                                  })
    }
    if (is.data.frame(dat$stats$transmat)) {
      nsteps <- get_control(dat, "nsteps")
      dat$stats$transmat <- padded_vector(list(dat$stats$transmat), nsteps)
    }
  }

  return(dat)
}


#' @title Disease Status Initialization Module for netsim
#'
#' @description This function sets the initial disease status on the
#'              network given the specified initial conditions.
#'
#' @inheritParams recovery.net
#'
#' @details
#' This internal function sets, either randomly or deterministically, the nodes
#' that are infected at \eqn{t_1}, the starting time of network simulations. If
#' the number to be initially infected is passed, this function sets the initial
#' number infected based on the number specified, either as a set of random
#' draws from a binomial distribution or as the exact number specified. In
#' either case, the specific nodes infected are a random sample from the
#' network. In contrast, a set of specific nodes may be infected by passing a
#' vector containing the status of each node to \code{\link{netsim}}.
#'
#' For the initially infected nodes, this module sets the time of infection as
#' \eqn{t_1}, the starting time of network simulations. For models with vital
#' dynamics, the infection time for those initially infected nodes is a random
#' draw from an exponential distribution with the rate parameter defined by the
#' \code{di.rate} argument. For models without vital dynamics, the infection
#' time is a random draw from a uniform distribution of integers with a minimum
#' of 1 and a maximum of the number of time steps in the model. In both cases,
#' to set the infection times to be in the past, these times are multiplied by
#' -1, and 2 is added to allow for possible infection times up until step 2,
#' when the disease simulation time loop starts.
#'
#' @inherit recovery.net return
#'
#' @seealso This is an initialization module for \code{\link{netsim}}.
#'
#' @export
#' @keywords netMod internal
#'
init_status.net <- function(dat) {

  type <- get_control(dat, "type", override.null.error = TRUE)
  type <- if (is.null(type)) "None" else type

  nsteps <- get_control(dat, "nsteps")
  tergmLite <- get_control(dat, "tergmLite")
  vital <- get_param(dat, "vital")
  groups <- get_param(dat, "groups")
  status.vector <- get_init(dat, "status.vector", override.null.error = TRUE)
  if (type %in% c("SIS", "SIR")) {
    rec.rate <- get_param(dat, "rec.rate")
  }
  if (vital == TRUE) {
    di.rate <- get_param(dat, "di.rate")
  }

  # Variables ---------------------------------------------------------------
  i.num <- get_init(dat, "i.num", override.null.error = TRUE)
  if (type  == "SIR" && is.null(status.vector)) {
    r.num <- get_init(dat, "r.num")
  }

  num <- sum(get_attr(dat, "active") == 1)

  if (groups == 2) {
    group <- get_attr(dat, "group")
    if (!all(group %in% c(1, 2))) {
      stop(
        "When using the `group` attribute, the only authorized values",
        " are 1 and 2.\n",
        "The values found were: ", paste0(unique(group), collapse = ", ")
      )
    }

    i.num.g2 <- get_init(dat, "i.num.g2")
    if (type  == "SIR" && is.null(status.vector)) {
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
      if (type == "SIR") {
        status[sample(which(group == 1 & status == "s"), size = r.num)] <- "r"
        if (groups == 2) {
          status[sample(which(group == 2 & status == "s"),
                        size = r.num.g2)] <- "r"
        }
      }
    }
    dat <- set_attr(dat, "status", status)
  } else {
    status <- get_attr(dat, "status") # already copied in copy_nwattr_to_datattr
  }


  ## Set up TEA status
  if (tergmLite == FALSE) {
    if (statOnNw == FALSE) {
      for (network in seq_len(dat$num.nw)) {
        dat$nw[[network]] <- set_vertex_attribute(dat$nw[[network]],
                                                  "status",
                                                  status)
      }
    }
    for (network in seq_len(dat$num.nw)) {
      dat$nw[[network]] <- activate.vertex.attribute(dat$nw[[network]],
                                                     prefix = "testatus",
                                                     value = status,
                                                     onset = 1,
                                                     terminus = Inf)
    }
  }


  # Infection Time ----------------------------------------------------------
  ## Set up inf.time vector
  if (type == "None") {
    infTime <- rep(NA, num)
    idsInf <- idsInf <- which(status == "i")
    infTime[idsInf] <- 1
    dat <- set_attr(dat, "infTime", infTime)
  } else {
    idsInf <- which(status == "i")
    infTime <- rep(NA, length(status))
    infTime.vector <- get_init(dat, "infTime.vector",
                               override.null.error = TRUE)

    if (!is.null(infTime.vector)) {
      infTime <- infTime.vector
    } else {
      # If vital dynamics, infTime is a geometric draw over the duration of
      # infection
      if (vital == TRUE && di.rate > 0) {
        if (type == "SI") {
          infTime[idsInf] <- -rgeom(n = length(idsInf), prob = di.rate) + 2
        } else {
          infTime[idsInf] <- -rgeom(n = length(idsInf),
                                    prob = di.rate +
                                      (1 - di.rate) * mean(rec.rate)) + 2
        }
      } else {
        if (type == "SI" || mean(rec.rate) == 0) {
          # if no recovery, infTime a uniform draw over the number of sim time
          # steps
          infTime[idsInf] <- ssample(1:(-nsteps + 2),
                                     length(idsInf), replace = TRUE)
        } else {
          infTime[idsInf] <- -rgeom(n = length(idsInf),
                                    prob = mean(rec.rate)) + 2
        }
      }
    }

    dat <- set_attr(dat, "infTime", infTime)
  }

  return(dat)
}

#' @title Network Data and Stats Initialization
#'
#' @description This function initializes the network data and stats on the main
#'              \code{dat} object.
#'
#' @param dat A \code{dat} object obtained from \code{\link{create_dat_object}},
#'        including the \code{control} argument.
#' @param x Either a fitted network model object of class \code{netest}, or a
#'        list of such objects.
#'
#' @return A \code{dat} main list object with network data and stats
#'         initialized.
#'
#' @export
#' @keywords internal
#'
init_nets <- function(dat, x) {
  if (inherits(x, "netest")) {
    x <- list(x)
  }

  ## initialize network data on dat object
  dat$num.nw <- length(x)
  dat$nwparam <- lapply(x, function(y) y[!(names(y) %in% c("fit", "newnetwork"))])
  nws <- lapply(x, `[[`, "newnetwork")
  nw <- nws[[1]]
  if (get_control(dat, "tergmLite") == TRUE) {
    dat$el <- lapply(nws, as.edgelist)
    dat$net_attr <- lapply(nws, get_network_attributes)
  } else {
    dat$nw <- nws
  }

  # Nodal Attributes --------------------------------------------------------

  # Standard attributes
  num <- network.size(nw)
  dat <- append_core_attr(dat, 1, num)

  groups <- length(unique(get_vertex_attribute(nw, "group")))
  dat <- set_param(dat, "groups", groups)

  ## Pull attr on nw to dat$attr
  dat <- copy_nwattr_to_datattr(dat, nw)

  ## record names of relevant vertex attributes
  dat$temp$nwterms <- get_network_term_attr(nw)

  ## initialize stats data structure
  if (get_control(dat, "save.nwstats") == TRUE) {
    if (get_control(dat, "resimulate.network") == TRUE) {
      dat$stats$nwstats <- rep(list(padded_vector(list(), get_control(dat, "nsteps"))),
                               length.out = length(dat$nwparam))
    } else {
      dat$stats$nwstats <- rep(list(list()), length.out = length(dat$nwparam))
    }
  }

  return(dat)
}
