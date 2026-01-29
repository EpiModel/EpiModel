#' @title Initialization: netsim Module
#'
#' @description This function initializes the main `netsim_dat` class data
#'              object on which data are stored, simulates the initial state of
#'              the networks, and simulates disease status and other attributes.
#'
#' @param x If `control$start == 1`, either a fitted network model object
#'        of class `netest` or a list of such objects. If
#'        `control$start > 1`, an object of class `netsim`. When
#'        multiple networks are used, the node sets (including network size
#'        and nodal attributes) are assumed to be the same for all networks.
#' @param param An `EpiModel` object of class [param.net()].
#' @param init An `EpiModel` object of class [init.net()].
#' @param control An `EpiModel` object of class [control.net()].
#' @param s Simulation number, used for restarting dependent simulations.
#' @details When re-initializing a simulation, the `netsim` object passed
#'          to `initialize.net` must contain the elements `param`,
#'          `nwparam`, `epi`, `coef.form`, and `num.nw`.
#'
#' @return A `netsim_dat` class main data object.
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
    if (!is.null(dat$run$nwterms)) {
      dat$run$t1.tab <- get_attr_prop(dat, dat$run$nwterms)
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
    required_names <- c("param", "nwparam", "epi", "run", "coef.form", "num.nw")
    missing_names <- setdiff(required_names, names(x))
    if (length(missing_names) > 0) {
      stop(
        "x is missing the following elements required for re-initialization: ",
        paste.and(missing_names), call. = FALSE
      )
    }

    # recycle sims in the restart object
    # e.g. 5 sim out of a size 3 restart object we will give: 1, 2, 3, 1, 2
    s <- (s - 1) %% length(x$run) + 1

    dat <- create_dat_object(
      param = param,
      control = control,
      run = x$run[[s]]
    )

    missing_params <- setdiff(names(x$param), names(param))
    for (mp in missing_params) {
      dat <- set_param(dat, mp, x$param[[mp]])
    }

    dat$num.nw <- x$num.nw

    dat$nwparam <- x$nwparam
    for (network in seq_len(dat$num.nw)) {
      dat$nwparam[[network]]$coef.form <- x$coef.form[[s]][[network]]
    }
    dat$epi <- sapply(x$epi, function(var) var[s])
    names(dat$epi) <- names(x$epi)

    dat$stats <- lapply(x$stats, function(var) var[[s]])
    if (get_control(dat, "save.nwstats") == TRUE) {
      nsteps <- get_control(dat, "nsteps")
      start <- get_control(dat, "start")
      dat$stats$nwstats <- lapply(dat$stats$nwstats,
        function(oldstats) padded_vector(list(oldstats), nsteps - start + 2L)
      )
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
#' vector containing the status of each node to [netsim()].
#'
#' For the initially infected nodes, this module sets the time of infection as
#' \eqn{t_1}, the starting time of network simulations. For models with vital
#' dynamics, the infection time for those initially infected nodes is a random
#' draw from an exponential distribution with the rate parameter defined by the
#' `di.rate` argument. For models without vital dynamics, the infection
#' time is a random draw from a uniform distribution of integers with a minimum
#' of 1 and a maximum of the number of time steps in the model. In both cases,
#' to set the infection times to be in the past, these times are multiplied by
#' -1, and 2 is added to allow for possible infection times up until step 2,
#' when the disease simulation time loop starts.
#'
#' @inherit recovery.net return
#'
#' @seealso This is an initialization module for [netsim()].
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

  statOnNw <- "status" %in% dat$run$nwterms

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
        dat$run$nw[[network]] <- set_vertex_attribute(dat$run$nw[[network]],
                                                      "status",
                                                      status)
      }
    }
    for (network in seq_len(dat$num.nw)) {
      dat$run$nw[[network]] <- activate.vertex.attribute(dat$run$nw[[network]],
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
#'              `netsim_dat` class data object.
#'
#' @param dat A main data object of class `netsim_dat` obtained from
#'        [create_dat_object()], including the `control`
#'        argument.
#' @param x Either a fitted network model object of class `netest`, or a
#'        list of such objects.
#'
#' @return A `netsim_dat` class main data object with network data and
#'         stats initialized.
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
    dat$run$el <- lapply(nws, as.edgelist)
    dat$run$net_attr <- lapply(nws, get_network_attributes)
  } else {
    dat$run$nw <- nws
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
  dat$run$nwterms <- get_network_term_attr(nw)

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

#' @title Helper to use a `data.frame` to initialize some attributes
#'
#' @description Uses `dat$init$init_attr` to overwrite some attributes of the
#' nodes at initialization
#'
#' @details
#' If an `init_attr` `data.frame` is present in `dat$init`, use it to overwrite
#' the attributes it contains.
#' `init_attr` must have a number of rows equal to the number of nodes in the
#' model as the attributes will be overwritten one to one, ensuring the correct
#' ordering.
#' `init_attr` columns MUST have a corresponding attribute already initialized.
#' See "R/default_attributes.R" for adding new attributes to the model.
#' `init_attr` is removed from `dat$init` at the end of the function to free up
#' its memory.
#'
#' @inheritParams recovery.net
#' @inherit recovery.net return
#'
#' @export
#'
overwrite_attrs <- function(dat) {
  init_attr <- get_init(dat, "init_attr", override.null.error = TRUE)
  if (is.null(init_attr)) {
    return(dat)
  }
  message("init_attr used for initialization of attributes")

  status <- get_attr(dat, "status")
  if (nrow(init_attr) != length(status)) {
    stop("init_attr should contains the same number of nodes as the model")
  }

  new_attrs <- setdiff(names(init_attr), names(dat$attr))
  if (length(new_attrs) > 0) {
    stop(
      "Some attributes in `init_attr` are not present in `dat`: ",
      paste0(new_attrs, collapse = " ,")
    )
  }

  core_attrs <- c("active", "entrTime", "exitTime", "unique_id")
  for (attr_name in setdiff(names(init_attr), core_attrs)) {
    dat <- set_attr(dat, attr_name, init_attr[[attr_name]])
  }

  dat$init$init_attr <- NULL
  dat
}
