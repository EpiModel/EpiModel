
#' @title Initialize Networks Used in netsim
#'
#' @description This function initializes the networks used in
#'              [netsim()]. The initial edge set for a given network
#'              is obtained either from simulating the cross-sectional model
#'              (if `edapprox == TRUE`) or from the `newnetwork`
#'              element of the `netest` object (if
#'              `edapprox == FALSE`). Once the initial edge sets are
#'              determined, the first time step is simulated if
#'              `resimulate.network == TRUE`, and all time steps are
#'              simulated if `resimulate.network == FALSE`. Initializes the
#'              `num(.g2)` epi fields used in
#'              [edges_correct()] for computing edge coefficient
#'              adjustments.
#'
#' @inheritParams recovery.net
#'
#' @inherit recovery.net return
#'
#' @export
#' @keywords netUtils internal
#'
sim_nets_t1 <- function(dat) {
  dat.updates <- NVL(get_control(dat, "dat.updates"), function(dat, ...) dat)

  ## simulate zeroth timestep cross-sectional network
  ## (default to newnetwork, set already in mod.initialize, if full tergm fit)
  dat <- dat.updates(dat = dat, at = 0L, network = 0L)
  for (network in seq_len(dat$num.nw)) {
    nwparam <- get_nwparam(dat, network = network)

    ## (re)construct input network
    nw <- get_network(dat, network)

    ## simulate t0 basis network nw if using edapprox
    if (nwparam$edapprox == TRUE) {
      nw <- simulate(nwparam$formula,
                     coef = nwparam$coef.form.crude,
                     basis = nw,
                     constraints = nwparam$constraints,
                     control = get_network_control(dat, network, "set.control.ergm"),
                     dynamic = FALSE)
    }

    if (get_control(dat, "tergmLite") == TRUE) {
      ## set up time and lasttoggle if tracking duration
      if (get_network_control(dat, network, "tergmLite.track.duration") == TRUE) {
        nw %n% "time" <- 0L
        nw %n% "lasttoggle" <- cbind(as.edgelist(nw), 0L)
      }
    } else {
      ## set up vertex and edge activity in networkDynamic case
      nw <- networkDynamic::as.networkDynamic(nw)
      nw <- networkDynamic::activate.vertices(nw, onset = 0L, terminus = Inf)
      nw <- networkDynamic::activate.edges(nw, onset = 0L, terminus = Inf)
      if (get_control(dat, "resimulate.network") == TRUE) {
        nw %n% "net.obs.period" <- list(observations = list(c(0, 1)),
                                        mode = "discrete",
                                        time.increment = 1,
                                        time.unit = "step")
      }
    }
    ## set network on dat object
    dat <- set_network(dat, network, nw)
    ## update dat object as needed
    dat <- dat.updates(dat = dat, at = 0L, network = network)
  }

  if (get_control(dat, "resimulate.network") == TRUE) {
    nsteps <- 1L
  } else {
    nsteps <- get_control(dat, "nsteps")
  }

  ## initialize num(.g2) run fields
  dat <- update_sim_num(dat)

  ## In tergmLite mode the post-simulate state replaces dat$run$el and the
  ## cross-section edgelist is no longer recoverable. Stash it now so
  ## seed_cumulative_edgelist_t1() can diff t=0 vs t=1 below. See #1016.
  if (get_control(dat, "tergmLite") && get_control(dat, "cumulative.edgelist")) {
    dat$run$el_t0_seed <- lapply(seq_len(dat$num.nw), function(network) {
      dat$run$el[[network]]
    })
  }

  ## simulate first timestep (if resimulate.network == TRUE)
  ## or all timesteps (if resimulate.network == FALSE)
  dat <- dat.updates(dat = dat, at = 1L, network = 0L)
  for (network in seq_len(dat$num.nw)) {
    dat <- simulate_dat(dat, at = 1L, network = network, nsteps = nsteps)
    dat <- dat.updates(dat = dat, at = 1L, network = network)
  }

  dat <- seed_cumulative_edgelist_t1(dat)

  return(dat)
}

#' @title Simulate a Network for a Specified Number of Time Steps
#'
#' @description This function simulates a dynamic network over one or multiple
#'              time steps for TERGMs or one or multiple cross-sectional network
#'              panels for ERGMs, for use in [netsim()] modeling.
#'              Network statistics are also extracted and saved if
#'              `save.nwstats == TRUE` and
#'              `resimulate.network == FALSE`.
#'
#' @inheritParams recovery.net
#' @param network index of the network to simulate
#' @param nsteps number of time steps to simulate
#'
#' @inherit recovery.net return
#'
#' @export
#' @keywords netUtils internal
#'
simulate_dat <- function(dat, at, network = 1L, nsteps = 1L) {
  ## determine formula and coefficients; set discordance_fraction in ergm case
  nwparam <- get_nwparam(dat, network = network)
  simulation_control <- get_network_control(dat, network, "set.control.tergm")
  if (nwparam$coef.diss$duration[1] > 1) {
    formula <- ~Form(nwparam$formation) +
      Persist(nwparam$coef.diss$dissolution)
    coef <- c(nwparam$coef.form, nwparam$coef.diss$coef.adj)
  } else {
    formula <- nwparam$formation
    coef <- nwparam$coef.form
    simulation_control$MCMC.prop.args <- list(discordance_fraction = 0)
  }

  ## determine output type
  if (get_control(dat, "tergmLite") == FALSE) {
    output <- "networkDynamic"
  } else {
    output <- "final"
  }

  ## determine monitor, if needed; note that we only obtain
  ## stats in simulate_dat if resimulate.network == FALSE
  if (get_control(dat, "save.nwstats") == TRUE &&
        get_control(dat, "resimulate.network") == FALSE) {
    monitor <- get_network_control(dat, network, "nwstats.formula")
  } else {
    monitor <- NULL # will be handled by summary_nets, if needed
  }

  if (get_control(dat, "tergmLite") == FALSE &&
        get_control(dat, "resimulate.network") == TRUE) {
    time_offset <- 0L
  } else {
    time_offset <- 1L
  }

  ## always TERGM simulation
  nw <- suppressWarnings(simulate(formula,
                                  coef = coef,
                                  basis = get_network(x = dat, network = network),
                                  constraints = nwparam$constraints,
                                  time.start = at - time_offset,
                                  time.offset = time_offset,
                                  time.slices = nsteps,
                                  output = output,
                                  control = simulation_control,
                                  monitor = monitor,
                                  dynamic = TRUE))

  ## update network (and el, if tergmLite) on the dat object
  dat <- set_network(x = dat, network = network, nw = nw)

  ## if monitor was used, record the results
  if (!is.null(monitor)) {
    new.nwstats <- attributes(nw)$stats
    keep.cols <- which(!duplicated(colnames(new.nwstats)))
    new.nwstats <- new.nwstats[, keep.cols, drop = FALSE]
    dat$stats$nwstats[[network]] <- list(as_tibble(new.nwstats))
  }

  return(dat)
}

#' @title Resimulate Dynamic Network at Time 2+
#'
#' @description This function resimulates the dynamic network in stochastic
#'              network models simulated in [netsim()] with dependence
#'              between the epidemic and demographic processes and the network
#'              structure.
#'
#' @inheritParams recovery.net
#'
#' @inherit recovery.net return
#'
#' @export
#' @keywords netUtils internal
#'
resim_nets <- function(dat, at) {

  # Calculate active attribute
  active <- get_attr(dat, "active")
  idsActive <- which(active == 1)
  anyActive <- ifelse(length(idsActive) > 0, TRUE, FALSE)
  if (dat$param$groups == 2) {
    group <- get_attr(dat, "group")
    groupids.1 <- which(group == 1)
    groupids.2 <- which(group == 2)
    nActiveG1 <- length(intersect(groupids.1, idsActive))
    nActiveG2 <- length(intersect(groupids.2, idsActive))
    anyActive <- ifelse(nActiveG1 > 0 & nActiveG2 > 0, TRUE, FALSE)
  }

  # Network resimulation, with dat.updates interspersed
  if (anyActive == TRUE && get_control(dat, "resimulate.network") == TRUE) {
    ## Edges Correction
    dat <- edges_correct(dat, at)
    dat <- update_sim_num(dat)

    ## network resimulation
    dat.updates <- NVL(get_control(dat, "dat.updates"), function(dat, ...) dat)
    dat <- dat.updates(dat = dat, at = at, network = 0L)
    for (network in seq_len(dat$num.nw)) {
      dat <- simulate_dat(dat = dat, at = at, network = network)
      dat <- dat.updates(dat = dat, at = at, network = network)
    }
  }

  # Cummulative edgelist
  truncate.el.cuml <- get_control(dat, "truncate.el.cuml")
  for (network in seq_len(dat$num.nw)) {
    dat <- update_cumulative_edgelist(dat, network, truncate.el.cuml)
  }

  return(dat)
}


#' @title Nodes Eligible to Form Ties for the Edges Correction
#'
#' @description Returns a logical vector marking the nodes the edges correction
#'              should count, or `NULL` when the `edges.correct.attr` control is
#'              unset and every node counts.
#'
#' @param dat Main `netsim_dat` object containing a `networkDynamic` object and
#'        other initialization information passed from [netsim()].
#' @param active The `active` nodal attribute.
#'
#' @return A logical vector the same length as `active`, or `NULL`.
#'
#' @keywords internal
#'
edges_correct_eligible <- function(dat, active) {
  attr.name <- get_control(dat, "edges.correct.attr",
                           override.null.error = TRUE)
  if (is.null(attr.name)) {
    return(NULL)
  }
  if (length(attr.name) != 1 || !is.character(attr.name)) {
    stop("The `edges.correct.attr` control must be a single attribute name.")
  }
  eligible <- get_attr(dat, attr.name)
  if (length(eligible) != length(active)) {
    stop("The nodal attribute `", attr.name, "` named by the ",
         "`edges.correct.attr` control has length ", length(eligible),
         ", but the `active` attribute has length ", length(active), ".")
  }
  !is.na(eligible) & eligible == 1
}


#' @title Adjustment for the Edges Coefficient with Changing Network Size
#'
#' @description Adjusts the edges coefficient in a dynamic network model
#'              simulated in [netsim()] to preserve the mean
#'              degree of nodes in the network. Requires `at >= 2`.
#'              Maintains the `num(.g2)` epi fields (initialized in
#'              [sim_nets_t1()]) for computing the coefficient
#'              adjustment.
#'
#' @details
#' By default the correction counts every active node. That is correct when
#' every active node can form a tie, which is the usual case.
#'
#' It is not correct when a model carries a subpopulation that stays active but
#' is structurally excluded from the network, for example an age band past a
#' sexual-cessation age whose target statistics are all zero, so that `ergm`
#' pins its terms off and no tie incident to those nodes can form. The
#' population the correction sees then grows while the population that can hold
#' an edge does not, and the whole of the adjustment lands on the nodes that
#' can. The `edges.correct.attr` control names a binary nodal attribute
#' identifying the nodes eligible to form ties, and the correction then counts
#' only those.
#'
#' The size of the error is the size of the excluded share. A model whose
#' excluded band grows to a quarter of the population dilutes mean degree among
#' the remainder by the same quarter, and none of it appears in the network
#' diagnostics, because [netdx()] runs before any node has been excluded.
#'
#' @inheritParams recovery.net
#'
#' @inherit recovery.net return
#'
#' @keywords internal
#' @export
#'
edges_correct <- function(dat, at) {

  resimulate.network <- get_control(dat, "resimulate.network")
  groups <- get_param(dat, "groups")
  active <- get_attr(dat, "active")

  if (resimulate.network == TRUE) {
    eligible <- edges_correct_eligible(dat, active)
    use.eligible <- !is.null(eligible)
    if (!use.eligible) {
      eligible <- TRUE
    }

    if (groups == 1) {
      old.num <- if (use.eligible) dat$run$num.elig else dat$run$num
      new.num <- sum(active == 1 & eligible)
      # A run restarted from output saved before this control existed carries
      # no eligible count. Treating the first step as a no-op is right: there
      # is no previous count to compare against.
      if (is.null(old.num)) {
        old.num <- new.num
      }
      adjustment <- log(old.num) - log(new.num)
    }
    if (groups == 2) {
      group <- get_attr(dat, "group")
      old.num.g1 <- if (use.eligible) dat$run$num.elig else dat$run$num
      old.num.g2 <- if (use.eligible) dat$run$num.elig.g2 else dat$run$num.g2
      new.num.g1 <- sum(active == 1 & group == 1 & eligible)
      new.num.g2 <- sum(active == 1 & group == 2 & eligible)
      if (is.null(old.num.g1)) {
        old.num.g1 <- new.num.g1
      }
      if (is.null(old.num.g2)) {
        old.num.g2 <- new.num.g2
      }
      adjustment <-
        log(2 * old.num.g1 * old.num.g2 / (old.num.g1 + old.num.g2)) -
        log(2 * new.num.g1 * new.num.g2 / (new.num.g1 + new.num.g2))
    }

    # An empty count makes the adjustment infinite, and an infinite adjustment
    # added to the edges coefficient makes every proposed tie be rejected for
    # the rest of the run, which looks like a slow collapse rather than an
    # error. Leaving the coefficients alone and saying so is more useful.
    if (!is.finite(adjustment)) {
      warning("edges_correct() computed a non-finite edges adjustment at ",
              "timestep ", at, " and left the edges coefficients unchanged. ",
              "This means a population count reached zero.")
      return(dat)
    }

    for (network in seq_len(dat$num.nw)) {
      dat$nwparam[[network]]$coef.form[1] <-
        dat$nwparam[[network]]$coef.form[1] + adjustment
    }
  }

  return(dat)
}

# Store the current number of nodes in the network or per network group in
# `dat$run`.
# This is used to adjuste the `edges` coefficient of tergm
update_sim_num <- function(dat) {
  active <- get_attr(dat, "active")
  # NULL unless the `edges.correct.attr` control is set, in which case these
  # are the counts `edges_correct()` reads on the next timestep.
  eligible <- edges_correct_eligible(dat, active)
  if (get_param(dat, "groups") == 1) {
    dat$run$num <- sum(active == 1)
    if (!is.null(eligible)) {
      dat$run$num.elig <- sum(active == 1 & eligible)
    }
  } else {
    group <- get_attr(dat, "group")
    dat$run$num <- sum(active == 1 & group == 1)
    dat$run$num.g2 <- sum(active == 1 & group == 2)
    if (!is.null(eligible)) {
      dat$run$num.elig <- sum(active == 1 & group == 1 & eligible)
      dat$run$num.elig.g2 <- sum(active == 1 & group == 2 & eligible)
    }
  }
  return(dat)
}
