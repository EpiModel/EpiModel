#' @title Merge Data across Stochastic Individual Contact Model Simulations
#'
#' @description Merges epidemiological data from two independent simulations of
#'              stochastic individual contact models from [icm()].
#'
#' @param x An `EpiModel` object of class [icm()].
#' @param y Another `EpiModel` object of class [icm()], with the
#'        identical model parameterization as `x`.
#' @param ...  Additional merge arguments (not used).
#'
#' @details
#' This merge function combines the results of two independent simulations of
#' [icm()] class models, simulated under separate function calls. The
#' model parameterization between the two calls must be exactly the same, except
#' for the number of simulations in each call. This allows for manual
#' parallelization of model simulations.
#'
#' This merge function does not work the same as the default merge, which allows
#' for a combined object where the structure differs between the input elements.
#' Instead, the function checks that objects are identical in model
#' parameterization in every respect (except number of simulations) and binds
#' the results.
#'
#' @return An `EpiModel` object of class [icm()] containing the
#'         data from both `x` and `y`.
#'
#' @method merge icm
#' @keywords extract
#' @export
#'
#' @examples
#' param <- param.icm(inf.prob = 0.2, act.rate = 0.8)
#' init <- init.icm(s.num = 1000, i.num = 100)
#' control <- control.icm(
#'   type = "SI", nsteps = 10,
#'   nsims = 3, verbose = FALSE
#' )
#' x <- icm(param, init, control)
#'
#' control <- control.icm(
#'   type = "SI", nsteps = 10,
#'   nsims = 1, verbose = FALSE
#' )
#' y <- icm(param, init, control)
#'
#' z <- merge(x, y)
#'
#' # Examine separate and merged data
#' as.data.frame(x)
#' as.data.frame(y)
#' as.data.frame(z)
#'
merge.icm <- function(x, y, ...) {
  ## Check structure
  if (length(x) != length(y) || !identical(names(x), names(y))) {
    stop("x and y have different structure")
  }
  if (
    x$control$nsims > 1 &&
      y$control$nsims > 1 &&
      !identical(sapply(x, class), sapply(y, class))
  ) {
    stop("x and y have different structure")
  }

  ## Check params
  if (!identical(x$param, y$param)) {
    stop("x and y have different parameters")
  }

  check_controls <- identical(
    x$control[-which(names(x$control) == "nsims")],
    y$control[-which(names(y$control) == "nsims")]
  )
  if (!check_controls) {
    stop("x and y have different controls")
  }

  z <- x
  new.range <- (x$control$nsims + 1):(x$control$nsims + y$control$nsims)

  # Merge data
  for (i in seq_along(x$epi)) {
    if (x$control$nsims == 1) {
      x$epi[[i]] <- data.frame(x$epi[[i]])
    }
    if (y$control$nsims == 1) {
      y$epi[[i]] <- data.frame(y$epi[[i]])
    }
    z$epi[[i]] <- cbind(x$epi[[i]], y$epi[[i]])
    names(z$epi[[i]])[new.range] <- paste0("sim", new.range)
  }

  z$control$nsims <- max(new.range)

  return(z)
}


#' @title Merge Model Simulations across netsim Objects
#'
#' @description Merges epidemiological data from two independent simulations of
#'              stochastic network models from `netsim`.
#'
#' @param x An `EpiModel` object of class [netsim()].
#' @param y Another `EpiModel` object of class [netsim()],
#'        with the identical model parameterization as `x`.
#' @param keep.transmat If `TRUE`, keep the transmission matrices from the
#'        original `x` and `y` elements. Note: transmission matrices
#'        only saved when (`save.transmat == TRUE`).
#' @param keep.network If `TRUE`, keep the `networkDynamic` objects
#'        from the original `x` and `y` elements. Note: network
#'        only saved when (`tergmLite == FALSE`).
#' @param keep.nwstats If `TRUE`, keep the network statistics (as set by
#'        the `nwstats.formula` parameter in `control.net`) from
#'        the original `x` and `y` elements.
#' @param keep.other If `TRUE`, keep the other simulation elements (as set
#'        by the `save.other` parameter in `control.net`) from the
#'        original `x` and `y` elements. Elements with their own `keep.*`
#'        argument (e.g. `run`) follow that argument instead.
#' @param param.error If `TRUE`, if `x` and `y` have different
#'        params (in [param.net()]) or controls (passed in
#'        [control.net()]) an error will prevent the merge. Use
#'        `FALSE` to override that check.
#' @param keep.diss.stats If `TRUE`, keep `diss.stats` from the
#'        original `x` and `y` objects.
#' @param keep.run If `TRUE`, keep the `run` sublists (as set by the
#'        `save.run` parameter in `control.net`) from the original `x` and
#'        `y` elements. These are required to restart from a simulation of the
#'        merged object in [netsim()], selected with [get_sims()].
#' @param keep.cumulative.edgelist If `TRUE`, keep the cumulative edgelists
#'        (as set by the `save.cumulative.edgelist` parameter in
#'        `control.net`) from the original `x` and `y` elements. `FALSE` by
#'        default, as these grow with the length of the simulation.
#' @param keep.attr.history If `TRUE`, keep the attribute histories (as set by
#'        [record_attr_history()]) from the original `x` and `y` elements.
#' @param keep.raw.records If `TRUE`, keep the raw records (as set by
#'        [record_raw_object()]) from the original `x` and `y` elements.
#' @param ...  Additional merge arguments (not currently used).
#'
#' @details
#' This merge function combines the results of two independent simulations of
#' [netsim()] class models, simulated under separate function calls.
#' The model parameterization between the two calls must be exactly the same,
#' except for the number of simulations in each call. This allows for manual
#' parallelization of model simulations.
#'
#' This merge function does not work the same as the default merge, which allows
#' for a combined object where the structure differs between the input elements.
#' Instead, the function checks that objects are identical in model
#' parameterization in every respect (except number of simulations) and binds
#' the results.
#'
#' After dropping the elements not kept (see the `keep.*` arguments), `x` and
#' `y` must hold the same elements.
#'
#' @return An `EpiModel` object of class [netsim()] containing
#'         the data from both `x` and `y`.
#'
#' @method merge netsim
#' @keywords extract
#' @export
#'
#' @examples
#' \donttest{
#' # Network model
#' nw <- network_initialize(n = 100)
#' coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 10)
#' est <- netest(nw,
#'   formation = ~edges, target.stats = 25,
#'   coef.diss = coef.diss, verbose = FALSE
#' )
#'
#' # Epidemic models
#' param <- param.net(inf.prob = 1)
#' init <- init.net(i.num = 1)
#' control <- control.net(
#'   type = "SI", nsteps = 20, nsims = 2,
#'   save.nwstats = TRUE,
#'   nwstats.formula = ~ edges + degree(0),
#'   verbose = FALSE
#' )
#' x <- netsim(est, param, init, control)
#' y <- netsim(est, param, init, control)
#'
#' # Merging
#' z <- merge(x, y)
#'
#' # Examine separate and merged data
#' as.data.frame(x)
#' as.data.frame(y)
#' as.data.frame(z)
#' }
#'
merge.netsim <- function(
  x,
  y,
  keep.transmat = TRUE,
  keep.network = TRUE,
  keep.nwstats = TRUE,
  keep.other = TRUE,
  param.error = TRUE,
  keep.diss.stats = TRUE,
  keep.run = TRUE,
  keep.cumulative.edgelist = FALSE,
  keep.attr.history = TRUE,
  keep.raw.records = TRUE,
  ...
) {
  x <- trim_netsim(
    x,
    keep.transmat = keep.transmat,
    keep.network = keep.network,
    keep.nwstats = keep.nwstats,
    keep.other = keep.other,
    keep.diss.stats = keep.diss.stats,
    keep.run = keep.run,
    keep.cumulative.edgelist = keep.cumulative.edgelist,
    keep.attr.history = keep.attr.history,
    keep.raw.records = keep.raw.records
  )

  y <- trim_netsim(
    y,
    keep.transmat = keep.transmat,
    keep.network = keep.network,
    keep.nwstats = keep.nwstats,
    keep.other = keep.other,
    keep.diss.stats = keep.diss.stats,
    keep.run = keep.run,
    keep.cumulative.edgelist = keep.cumulative.edgelist,
    keep.attr.history = keep.attr.history,
    keep.raw.records = keep.raw.records
  )

  # TODO: check param / control there ?
  #       return a list of elements (persim, stats, ...)
  similar_components <- check_similar_sim(x, y)

  # Check that `param` and `control` are identical if `param.error == TRUE`
  if (param.error) {
    if (!identical(x$param, y$param)) {
      stop("x and y have different parameters")
    }
    fmla_controls <- c(
      "monitors",
      "nwstats.formula",
      "set.control.tergm",
      "set.control.ergm",
      "dat.updates",
      "future.use.plan"
    )
    check_controls <- identical(
      x$control[setdiff(names(x$control), c("nsims", fmla_controls))],
      y$control[setdiff(names(y$control), c("nsims", fmla_controls))]
    )

    ## handle formulas separately due to environments
    check_controls <- Reduce(
      function(init, elt) {
        init && isTRUE(all.equal(x$control[[elt]], y$control[[elt]]))
      },
      fmla_controls,
      init = check_controls
    )

    if (!check_controls) {
      stop("x and y have different controls")
    }
  }

  # Perform the merging
  out <- x
  out$control$nsims <- sum(similar_components$n_sims)
  newnames <- get_sim_names(out$control$nsims)

  # Merge epi data
  for (i in seq_along(x$epi)) {
    out$epi[[i]] <- cbind(x$epi[[i]], y$epi[[i]])
    names(out$epi[[i]]) <- newnames
  }

  # Merge per simulation elements
  for (elt in similar_components$per_sim) {
    out[[elt]] <- c(x[[elt]], y[[elt]])
    names(out[[elt]]) <- newnames
  }

  for (elt in similar_components$stats) {
    out$stats[[elt]] <- c(x$stats[[elt]], y$stats[[elt]])
    names(out$stats[[elt]]) <- newnames
  }

  return(out)
}

check_similar_sim <- function(sim1, sim2) {
  save.other <- union(sim1$control$save.other, sim2$control$save.other)

  # Same top level names
  if (!setequal(names(sim1), names(sim2))) {
    stop_mismatch("elements", names(sim1), names(sim2), save.other)
  }

  sim1_n_sims <- unique(vapply(sim1$epi, ncol, 1L))
  sim2_n_sims <- unique(vapply(sim2$epi, ncol, 1L))

  # Same per sim elements
  sim1_per_sim <- get_per_sim_element_names(sim1, sim1_n_sims)
  sim2_per_sim <- get_per_sim_element_names(sim2, sim2_n_sims)

  if (!setequal(sim1_per_sim, sim2_per_sim)) {
    stop_mismatch(
      "per-simulation elements", sim1_per_sim, sim2_per_sim, save.other,
      problem = "is not per-simulation in"
    )
  }

  # similar epi
  if (!setequal(names(sim1$epi), names(sim2$epi))) {
    stop("sim1 and sim2 do not save the same epi Trackers")
  }
  sim1_n_steps <- unique(vapply(sim1$epi, nrow, 1L))
  sim2_n_steps <- unique(vapply(sim2$epi, nrow, 1L))
  if (sim1_n_steps != sim2_n_steps) {
    stop("sim1 and sim2 do not have the same number of steps")
  }

  # same `num.nw`
  if (sim1$num.nw != sim2$num.nw) {
    stop("sim1 and sim2 do not have the same number of networks")
  }

  # Do they have `run`
  sim1_has_run <- "run" %in% sim1_per_sim
  sim2_has_run <- "run" %in% sim2_per_sim

  # similar `attrs`
  if (sim1_has_run && sim2_has_run) {
    sim1_attrs_names <- unique(unlist(lapply(sim1$run, \(x) names(x$attr))))
    sim2_attrs_names <- unique(unlist(lapply(sim2$run, \(x) names(x$attr))))
    if (!setequal(sim1_attrs_names, sim2_attrs_names)) {
      stop("sim1 and sim2 do not store the same attributes")
    }
  }

  # nwstats & transmat
  sim1_stats <- get_per_sim_element_names(sim1$stats, sim1_n_sims)
  sim2_stats <- get_per_sim_element_names(sim2$stats, sim2_n_sims)

  if (!setequal(sim1_stats, sim2_stats)) {
    stop_mismatch("stats", sim1_stats, sim2_stats, save.other, prefix = "stats$")
  }

  list(
    n_sims = c(sim1_n_sims, sim2_n_sims),
    per_sim = sim1_per_sim,
    stats = sim1_stats
  )
}

# Name of the `keep.*` argument governing each element, NA if none
keep_arg_name <- function(elts, save.other) {
  own <- c(
    "run", "cumulative.edgelist", "attr.history", "raw.records",
    "network", "diss.stats", "transmat", "nwstats"
  )
  ifelse(
    elts %in% own, paste0("keep.", elts),
    ifelse(elts %in% save.other, "keep.other", NA_character_)
  )
}

# Error listing the elements held by only one of `x` and `y`
stop_mismatch <- function(what, x_elts, y_elts, save.other,
                          problem = "is missing from", prefix = "") {
  elts <- union(setdiff(x_elts, y_elts), setdiff(y_elts, x_elts))
  side <- ifelse(elts %in% x_elts, "y", "x")
  arg <- keep_arg_name(elts, save.other)
  hint <- ifelse(is.na(arg), "", paste0(" (set `", arg, " = FALSE` to drop it)"))
  stop(
    "x and y do not hold the same ", what, ":\n",
    paste0("  - `", prefix, elts, "` ", problem, " ", side, hint, "\n", collapse = ""),
    call. = FALSE
  )
}

trim_netsim <- function(
  sim,
  keep.transmat,
  keep.network,
  keep.nwstats,
  keep.other,
  keep.diss.stats,
  keep.run,
  keep.cumulative.edgelist,
  keep.attr.history,
  keep.raw.records
) {
  top_level <- c(
    "run",
    "cumulative.edgelist",
    "attr.history",
    "raw.records",
    "network",
    "diss.stats"
  )

  # elements with their own `keep.*` argument are not governed by `keep.other`
  other <- setdiff(sim$control$save.other, top_level)
  other_elts <- setNames(rep(keep.other, length(other)), other)

  top_level <- vapply(top_level, \(x) get(paste0("keep.", x)), logical(1))
  top_level <- c(other_elts, top_level)

  for (elt in names(top_level)) {
    if (!top_level[elt]) {
      sim[[elt]] <- NULL
    }
  }

  if (!keep.transmat) {
    sim$stats$transmat <- NULL
  }
  if (!keep.nwstats) {
    sim$stats$nwstats <- NULL
  }
  if (!keep.transmat && !keep.nwstats) {
    sim$stats <- NULL
  }

  sim
}

get_per_sim_element_names <- function(obj, n_sims) {
  sim_names <- get_sim_names(n_sims)
  names(Filter(function(v) setequal(v, sim_names), lapply(obj, names)))
}

get_sim_names <- function(n_sims) {
  paste0("sim", seq_len(n_sims))
}

# TODO: get_sims should not warn on improper `netsim`. If we want this, we
# should have a dedicated function. I think it's a good idea but if we plan to
# update the netsim / dat structures, I would wait until then
# TODO: similarly, caclulating `nsims` with the epi. That breaks if epi's have
# different length. But I would argue it's out of scope and should be fix with a
# general netsim object checker longer term. Current use of `control$nsims` is
# no better as it simply does not check anything
# TODO: explain these differences when push
