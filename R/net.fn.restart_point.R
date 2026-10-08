#' Make a Lightweight Restart Point From a `netsim` Object with tergmLite
#'
#' Extract the elements required for re-initializing a `netsim` simulation from
#' a completed simulation. This function also resets the Unique IDs and Time
#' values to reduce the size of the simulation. This function only works for
#' simulations where `control$tergmLite = TRUE`
#'
#' @param sim_obj a `netsim` object from an ended `netsim` call.
#' @param sim_num the number of the simulation to extract from the `netsim`
#'        object (default = 1).
#' @param keep_steps The number of simulation steps to keep from the previous
#'        run. By default only keep one but more is possible if some
#'        back-history is wanted.
#' @param time_attrs a `character` vector containing the names of the attributes
#'        that are expressed in time-steps. These will be offset so the last
#'        step in the original simulation become the step 1 (default) in the new
#'        ones. If no such attributes exist, pass `c()`.
#'
#' @details
#' The restart point created always contains a single simulation and drops the
#' `attr.history` and the `raw.records` from the initial simulation.
#'
#' `sim_obj` does not need to hold the parameters: when restarting, [netsim()]
#' takes them from its `param` argument alone, which must hold all the
#' parameters of the model.
#'
#' The epi trackers, cumulative edgelists, transmission matrix and `nwstats` are
#' truncated to only contain the last `keep_steps` entries. The unique IDs in
#' the transmission matrix are re-based like the ones of the attributes.
#'
#' Warning: the `time_attrs` argument is mandatory. Almost all simulations worth
#' restarting have such attributes (e.g. time.of.hiv.infection). If no such
#' argument exists, passing `c()` will allow the function to run while ensuring
#' that this was done on purpose.
#'
#' When restarting from the output of this function, it is suggested to express
#' the time steps in a relative manner in `control.net`:
#' ```
#' control.net(
#'   start = restart_point$control$nsteps + 1,
#'   nsteps = restart_point$control$nsteps + 1 + 104
#' )
#' ```
#'
#' @examples
#' \dontrun{
#' # With  pre-existing `sim`, `param` and `init` object (see `netsim`)
#'
#' # List all attributes that store a time step
#' time_attrs <- c(
#'   "inf.time",
#'   "stage.time",
#'   "aids.time",
#'   "prep.start.last"
#' )
#' # Make a restart point from simulation 1, re-run for 10 more timesteps
#' x <- make_restart_point(sim, time_attrs, sim_num = 1, keep_steps = 1)
#' control <- control_msm(
#'   start = x$control$nsteps + 1,
#'   nsteps = x$control$nsteps + 1 + 10
#' )
#' sim <- netsim(x, param, init, control)
#' }
#'
#' @return a trimmed `netsim` object with only one simulation that is ready to
#'         be used as a restart point.
#'
#' @export
make_restart_point <- function(
  sim_obj,
  time_attrs,
  sim_num = 1,
  keep_steps = 1
) {
  if (!inherits(sim_obj, c("netsim"))) {
    stop("`sim_obj` must be an object of class `netsim`")
  }
  required_names <- c(
    "control",
    "nwparam",
    "epi",
    "run",
    "coef.form",
    "num.nw"
  )
  missing_names <- setdiff(required_names, names(sim_obj))
  if (length(missing_names) > 0) {
    stop(
      "`sim_obj` is missing the following elements required for",
      " re-initialization: ",
      paste.and(missing_names)
    )
  }

  nsims <- sim_obj$control$nsims
  if (length(sim_num) != 1 || !sim_num %in% seq_len(nsims)) {
    stop(
      "`sim_num` must be a single number >= 1 and <= `control$nsims` (",
      nsims, ")"
    )
  }

  if (!sim_obj$control$tergmLite) {
    stop("Only `netsim` object with `tergmLite == TRUE` are supported")
  }

  # Select the simulation of interest, that renames the selected sim: `sim1`
  x <- get_sims(sim_obj, sims = sim_num)
  n_steps <- x$control$nsteps
  run_ls <- x$run[[1]]

  # Keep only the last `keep_steps` rows of each epi
  if (keep_steps < 1 || keep_steps > n_steps) {
    stop("`keep_steps` must be >= 1 and <= `sim_obj$control$nsteps`")
  }
  keep_rows <- (n_steps - keep_steps + 1):n_steps
  x$epi <- lapply(x$epi, function(r) r[keep_rows, , drop = FALSE])

  # Time correction
  time_offset <- n_steps - keep_steps
  x$control$start <- 1
  x$control$nsteps <- keep_steps

  # Time attributes - offset so last step is now `keep_steps`
  time_attrs <- union(c("entrTime", "exitTime"), time_attrs)
  missing_attrs <- setdiff(time_attrs, names(run_ls$attr))
  if (length(missing_attrs) > 0) {
    stop(
      "Some time attributes are not present in the attributes list: ",
      paste.and(missing_attrs)
    )
  }
  run_ls$attr[time_attrs] <- lapply(
    run_ls$attr[time_attrs],
    function(v) v - time_offset
  )

  # Fix UIDs
  uid_offset <- min(run_ls$attr$unique_id) - 1
  run_ls$attr$unique_id <- run_ls$attr$unique_id - uid_offset
  run_ls$last_unique_id <- run_ls$last_unique_id - uid_offset

  # Cumulative Edgelist - fix time and UIDs
  run_ls$el_cuml_cur <- lapply(
    run_ls$el_cuml_cur,
    function(el) {
      el$head <- el$head - uid_offset
      el$tail <- el$tail - uid_offset
      el$start <- el$start - time_offset
      el
    }
  )
  # For Historical one - truncate to 1 (only edges in the kept history)
  run_ls$el_cuml_hist <- lapply(
    run_ls$el_cuml_hist,
    function(el) {
      el$head <- el$head - uid_offset
      el$tail <- el$tail - uid_offset
      el$start <- el$start - time_offset
      el$stop <- el$stop - time_offset
      el[el$stop >= 1, , drop = FALSE]
    }
  )

  # the edgelist stores the name of the vertices. We don't use it with
  # `tergmLite` and it takes a lot of space
  run_ls$el <- lapply(run_ls$el, function(x) {
    attr(x, "vnames") <- NULL
    x
  })

  x$run[[1]] <- run_ls

  # If transmat was saved, trim it, offset the `at` column and the UIDs (a
  # simulation without transmissions holds an empty `tibble`, left as is)
  if (!is.null(x$stats$transmat) && nrow(x$stats$transmat[[1]]) > 0) {
    tsmt <- x$stats$transmat[[1]]
    tsmt$at <- tsmt$at - time_offset
    for (uid_col in intersect(c("sus", "inf"), names(tsmt))) {
      tsmt[[uid_col]] <- tsmt[[uid_col]] - uid_offset
    }
    x$stats$transmat[[1]] <- tsmt[tsmt$at > 0, , drop = FALSE]
  }

  # If `nwstats` are saved, keep only the last rows
  if (!is.null(x$stats$nwstats)) {
    x$stats$nwstats[[1]] <- lapply(
      x$stats$nwstats[[1]],
      function(d) d[keep_rows, , drop = FALSE]
    )
  }

  # Output ---------------------------------------------------------------------
  x$attr.history <- list()
  x$raw.records <- list()
  return(x)
}
