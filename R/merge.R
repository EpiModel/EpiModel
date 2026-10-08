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
#'        original `x` and `y` elements.
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
#' @param keep.attr.history If `TRUE`, keep the recorded histories (as set by
#'        [record_attr_history()] and [record_raw_object()]) from the original
#'        `x` and `y` elements. This governs both `attr.history` and
#'        `raw.records`.
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
#' Only the per-simulation elements not kept (see the `keep.*` arguments) may
#' be missing from one of `x` and `y`, e.g. from the result of a previous merge
#' that dropped them. The transmission matrices and network statistics are an
#' exception: they are dropped when missing from one of them.
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
  ...
) {
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
      "dat.updates"
    )
    check_controls  <- identical(
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

  ## Per-simulation elements, and the argument keeping each one. The formation
  ## coefficients are always kept, as they are required along with `run` to
  ## restart from the merged object. `attr.history` and `raw.records` are the
  ## two halves of the same recording facility and are saved together by
  ## `saveout.net`. `keep.other` has the final say on the elements requested
  ## through `save.other`, including the ones above (e.g. `run`).
  keep_args <- c(
    network = "keep.network",
    diss.stats = "keep.diss.stats",
    run = "keep.run",
    cumulative.edgelist = "keep.cumulative.edgelist",
    attr.history = "keep.attr.history",
    raw.records = "keep.attr.history"
  )
  keep_args[union(x$control$save.other, y$control$save.other)] <- "keep.other"
  keep_elts <- c(
    coef.form = TRUE,
    vapply(keep_args, get, logical(1), envir = environment())
  )
  keep_stats_elts <- c(nwstats = keep.nwstats, transmat = keep.transmat)

  ## Check structure. Only the per-simulation elements not kept may differ,
  ## e.g. when `x` is the result of a previous merge that dropped them.
  kept_elts <- setdiff(union(names(x), names(y)), names(keep_elts)[!keep_elts])
  in_x <- kept_elts %in% names(x)
  in_y <- kept_elts %in% names(y)
  if (!all(in_x & in_y)) {
    missing_elts <- kept_elts[!(in_x & in_y)]
    drop_hints <- ifelse(
      missing_elts %in% names(keep_args),
      paste0(" (set `", keep_args[missing_elts], " = FALSE` to drop it)"),
      ""
    )
    stop(
      "x and y have different structure, the elements kept in the merge ",
      "must be in both:\n",
      paste0(
        "  - `", missing_elts, "` is missing from ",
        ifelse(in_x[!(in_x & in_y)], "y", "x"), drop_hints, "\n",
        collapse = ""
      )
    )
  }
  if (x$control$nsims > 1 && y$control$nsims > 1) {
    x_classes <- vapply(x[kept_elts], function(i) class(i)[1], character(1))
    y_classes <- vapply(y[kept_elts], function(i) class(i)[1], character(1))
    diff_classes <- x_classes != y_classes
    if (any(diff_classes)) {
      stop(
        "x and y have different structure, the elements kept in the merge ",
        "must be of the same class in both:\n",
        paste0(
          "  - `", kept_elts[diff_classes], "`: `", x_classes[diff_classes],
          "` in x, `", y_classes[diff_classes], "` in y\n",
          collapse = ""
        )
      )
    }
  }

  # Perform the merging
  out <- x
  out$control$nsims <- as.integer(x$control$nsims + y$control$nsims)
  newnames <- paste0("sim", seq_len(out$control$nsims))

  # Merge epi data
  for (i in seq_along(x$epi)) {
    out$epi[[i]] <- cbind(x$epi[[i]], y$epi[[i]])
    names(out$epi[[i]]) <- newnames
  }

  ## Per-simulation elements: the kept ones are bound, the others dropped. A
  ## kept element held by neither side, or empty on both (e.g. the
  ## `attr.history` of restart points built by `make_restart_point()`), is
  ## left as is.
  for (elt in names(keep_elts)) {
    if (!keep_elts[[elt]]) {
      out[[elt]] <- NULL
    } else if (length(x[[elt]]) > 0 || length(y[[elt]]) > 0) {
      out[[elt]] <- c(x[[elt]], y[[elt]])
      names(out[[elt]]) <- newnames
    }
  }

  for (elt in names(keep_stats_elts)) {
    if (
      keep_stats_elts[[elt]] &&
        length(x$stats[[elt]]) > 0 &&
        length(y$stats[[elt]]) > 0
    ) {
      out$stats[[elt]] <- c(x$stats[[elt]], y$stats[[elt]])
      names(out$stats[[elt]]) <- newnames
    } else {
      out$stats[[elt]] <- NULL
    }
  }

  return(out)
}
