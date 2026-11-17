# File names of the `<n>.rds` elements of a restart pool directory, in order
#' @noRd
get_restart_pool_files <- function(pool_dir) {
  pool <- list.files(pool_dir, pattern = "^[0-9]+\\.rds$")
  pool <- pool[!dir.exists(file.path(pool_dir, pool))]
  pool[order(as.integer(sub("\\.rds$", "", pool)))]
}

# The `pool_index.csv` of the restart pool at `pool_dir`. A directory is a
# restart pool if it holds one, with a `pool_num` column numbering the
# elements `1` to `N` and an `epimodel_version` column
#' @noRd
read_restart_pool_index <- function(pool_dir) {
  index_path <- file.path(pool_dir, "pool_index.csv")
  if (!file.exists(index_path)) {
    stop(
      "\"", pool_dir, "\" is not a restart pool: it has no `pool_index.csv`. ",
      "Make restart pools with `make_restart_pool`"
    )
  }
  pool_index <- utils::read.csv(index_path)
  missing_cols <- setdiff(c("pool_num", "epimodel_version"), names(pool_index))
  if (length(missing_cols) > 0) {
    stop(
      "\"", index_path, "\" has no ",
      paste0("`", missing_cols, "`", collapse = " or "),
      " column. Remake the pool with `make_restart_pool`"
    )
  }
  if (nrow(pool_index) == 0 ||
        !identical(pool_index$pool_num, seq_len(nrow(pool_index)))) {
    stop(
      "The `pool_num` column of \"", index_path, "\" must number the pool ",
      "elements `1` to `N`"
    )
  }
  pool_index$epimodel_version <- as.character(pool_index$epimodel_version)
  pool_index
}

# Stop unless `version`, the EpiModel version a restart point was made with,
# has the major.minor version of the installed EpiModel: `netsim` does not
# restart from other ones. `what` names the restart point in the error.
# Shared by the restart pool functions and by `initialize.net`.
#' @noRd
check_restart_point_version <- function(version, what) {
  installed <- utils::packageVersion("EpiModel")
  if (is.null(version)) version <- NA_character_
  version <- numeric_version(as.character(version), strict = FALSE)
  if (is.na(version) || format(version[, 1:2]) != format(installed[, 1:2])) {
    stop(
      what, " was made with EpiModel ",
      if (is.na(version)) "< 2.7 (no version recorded)" else paste0("v", version),
      " but EpiModel v", installed, " is installed: `netsim` cannot restart ",
      "from it. Please rebuild the restart object with the current version"
    )
  }
  invisible(TRUE)
}

# The hidden sibling directory `make_restart_pool` writes the pool at
# `pool_dir` in, before renaming it to `pool_dir`
#' @noRd
get_partial_pool_dir <- function(pool_dir) {
  file.path(dirname(pool_dir), paste0(".", basename(pool_dir), ".partial"))
}

#' Make a Restart Pool From Simulations Run With Scenarios
#'
#' Creates a directory of single simulation restart points (`1.rds`, `2.rds`,
#' ...) from the batch files written by [netsim_scenarios()]. Such a directory
#' can be passed as `path_to_x` to the scenario runners (see
#' [netsim_path_wrapper()]).
#'
#' Pool element `i` is made from row `i` of `keep_sims`, with
#' [make_restart_point()] (keeping only the last time step). Each batch file is
#' read once, and each element is written to the pool as soon as it is made, so
#' only one batch file's elements are in memory at a time.
#' The pool directory also receives a `pool_index.csv` file
#' with one row per pool element: `pool_num`, `scenario_name`, `batch_number`,
#' `sim_number`, `file_path` (the batch file it was made from) and
#' `epimodel_version` (the EpiModel version recorded in the element, which
#' `netsim` requires to have the installed major.minor version to restart from
#' it). A directory is used as a restart pool only if it holds this file.
#'
#' The pool is written in a hidden `.<pool name>.partial` directory next to
#' `pool_dir`, renamed to `pool_dir` once complete: `pool_dir` never holds a
#' partial pool, even if R is killed. The partial directory is removed if an
#' error occurs, or at the next call if R was killed.
#'
#' The EpiModelHPC package provides `slurmworkflow` step templates that run
#' [netsim_scenarios()] and this function on a SLURM cluster.
#'
#' @param pool_dir Path to the restart pool directory to create. It must not
#'   exist.
#' @param scenarios_dir The directory containing the
#'   `sim__<scenario>__<batch>.rds` files.
#' @param scenario_name The name of the scenario to make the pool from.
#' @param keep_sims A `data.frame` with a `batch_number` and a `sim_number`
#'   column, one row per simulation to put in the pool. `sim_number` is the
#'   position of the simulation within its batch file (the `sim_number`
#'   column of [merge_netsim_scenarios_tibble()], not the `sim` one).
#' @param time_attrs A character vector of the attributes holding time steps,
#'   passed to [make_restart_point()].
#'
#' @return The pool index `data.frame` (invisibly).
#'
#' @seealso [validate_restart_pool()], [netsim_path_wrapper()]
#'
#' @examples
#' \dontrun{
#' d_sim <- readRDS("data/run/calib/merged_tibbles/df__default.rds")
#' keep_sims <- dplyr::distinct(d_sim, batch_number, sim_number)[1:32, ]
#' make_restart_pool(
#'   pool_dir = "data/run/restart_pool",
#'   scenarios_dir = "data/run/calib",
#'   scenario_name = "default",
#'   keep_sims = keep_sims,
#'   time_attrs = c("infTime")
#' )
#' validate_restart_pool("data/run/restart_pool")
#' }
#'
#' @export
make_restart_pool <- function(
  pool_dir,
  scenarios_dir,
  scenario_name,
  keep_sims,
  time_attrs
) {
  if (dir.exists(pool_dir)) {
    stop("A restart pool already exists at: \"", pool_dir, "\"")
  }
  if (!all(c("batch_number", "sim_number") %in% names(keep_sims))) {
    stop("`keep_sims` must have a `batch_number` and a `sim_number` column")
  }
  if (nrow(keep_sims) == 0) {
    stop("`keep_sims` has no rows")
  }
  if (anyDuplicated(keep_sims[c("batch_number", "sim_number")]) > 0) {
    stop("`keep_sims` contains the same simulation more than once")
  }

  b_infos <- get_scenarios_batches_infos(scenarios_dir)
  b_infos <- b_infos[b_infos$scenario_name == scenario_name, , drop = FALSE]
  missing_batches <- setdiff(keep_sims$batch_number, b_infos$batch_number)
  if (length(missing_batches) > 0) {
    stop(
      "No simulation file for scenario \"", scenario_name, "\" and batch(es) ",
      paste(missing_batches, collapse = ", "), " in \"", scenarios_dir, "\""
    )
  }

  partial_dir <- get_partial_pool_dir(pool_dir)
  if (dir.exists(partial_dir)) {
    message(
      "Removing the partial restart pool left by an interrupted run at: \n\"",
      partial_dir, "\""
    )
    unlink(partial_dir, recursive = TRUE)
  }
  on.exit(
    if (dir.exists(partial_dir)) unlink(partial_dir, recursive = TRUE),
    add = TRUE
  )
  dir.create(partial_dir, recursive = TRUE)

  file_paths <- character(nrow(keep_sims))
  epimodel_versions <- character(nrow(keep_sims))
  for (batch in unique(keep_sims$batch_number)) {
    rows <- which(keep_sims$batch_number == batch)
    batch_path <- b_infos$file_path[b_infos$batch_number == batch]
    sim_obj <- readRDS(batch_path)
    out_of_range <- setdiff(
      keep_sims$sim_number[rows],
      seq_len(sim_obj$control$nsims)
    )
    if (length(out_of_range) > 0) {
      stop(
        "`sim_number` ", paste(out_of_range, collapse = ", "),
        " out of range for batch ", batch, " (", sim_obj$control$nsims,
        " simulations)"
      )
    }
    for (i in rows) {
      restart_point <- make_restart_point(
        sim_obj,
        time_attrs = time_attrs,
        sim_num = keep_sims$sim_number[i]
      )
      x_ver <- attr(restart_point, "version")
      epimodel_versions[i] <-
        if (is.null(x_ver)) NA_character_ else as.character(x_ver)
      saveRDS(restart_point, file.path(partial_dir, paste0(i, ".rds")))
    }
    file_paths[rows] <- batch_path
    rm(sim_obj)
  }

  pool_index <- data.frame(
    pool_num = seq_len(nrow(keep_sims)),
    scenario_name = scenario_name,
    batch_number = keep_sims$batch_number,
    sim_number = keep_sims$sim_number,
    file_path = as.character(file_paths),
    epimodel_version = epimodel_versions
  )
  utils::write.csv(
    pool_index,
    file.path(partial_dir, "pool_index.csv"),
    row.names = FALSE
  )
  # a single `rename(2)`: `partial_dir` is on the same filesystem
  if (!file.rename(partial_dir, pool_dir)) {
    stop(
      "Could not rename the partial restart pool at \"", partial_dir,
      "\" to \"", pool_dir, "\""
    )
  }

  invisible(pool_index)
}

#' Validate a Restart Pool
#'
#' Checks that a restart pool (a directory made by [make_restart_pool()],
#' holding a `pool_index.csv`) can be used by [netsim_path_wrapper()]. This
#' function is not called by the simulation functions: call it after making a
#' pool, or before starting a long run from it.
#'
#' Errors are raised when:
#' - `path` is a directory without a `pool_index.csv`, or its `pool_index.csv`
#'   has no `pool_num` or `epimodel_version` column, or does not number the
#'   elements `1` to `N`;
#' - an element listed in `pool_index.csv` has no `<pool_num>.rds` file;
#' - an element is not a single simulation `netsim` restart point with the
#'   elements needed to restart (`param`, `control`, `nwparam`, `epi`, `run`,
#'   `coef.form` and `num.nw`);
#' - an element was made with another major.minor version of EpiModel than
#'   the one installed, according to `pool_index.csv` (checked by
#'   [netsim_path_wrapper()] before running) or to the element itself (checked
#'   by `netsim`). Run this function with the EpiModel version the simulations
#'   will use (e.g. in the HPC `renv`);
#' - the elements do not all hold the same number of time steps: `netsim`
#'   requires `control$start` to be that number + 1 for each of them.
#'
#' Warnings are raised when:
#' - the pool directory holds `<n>.rds` files not listed in `pool_index.csv`,
#'   which are never used;
#' - `path` is a single file holding several simulations:
#'   [netsim_path_wrapper()] restarts every replicate from its first
#'   simulation.
#'
#' @param path Path to a restart pool directory, or to a single restart point
#'   file.
#'
#' @return A `data.frame` with one row per pool element (invisibly):
#'   `pool_num`, `file_path`, `n_runs`, `nsteps` (the number of time steps
#'   kept in the element, `1` for a pool made by [make_restart_pool()]) and
#'   `epimodel_version`.
#'
#' @seealso [make_restart_pool()], [netsim_path_wrapper()]
#'
#' @export
validate_restart_pool <- function(path) {
  is_pool_dir <- dir.exists(path)
  if (is_pool_dir) {
    pool_index <- read_restart_pool_index(path)
    # The check `netsim_path_wrapper` makes before running
    for (n in pool_index$pool_num) {
      check_restart_point_version(
        pool_index$epimodel_version[n],
        paste0("Element ", n, " of the restart pool at \"", path, "\"")
      )
    }
    pool <- paste0(pool_index$pool_num, ".rds")
    file_paths <- file.path(path, pool)
    missing_paths <- file_paths[!file.exists(file_paths)]
    if (length(missing_paths) > 0) {
      stop(
        "Missing restart pool element(s): ",
        paste0("\"", missing_paths, "\"", collapse = ", ")
      )
    }
    unlisted <- setdiff(get_restart_pool_files(path), pool)
    if (length(unlisted) > 0) {
      warning(
        "\"", path, "\" holds element files not listed in its ",
        "`pool_index.csv`, which are never used: ",
        paste(unlisted, collapse = ", ")
      )
    }
  } else if (file.exists(path)) {
    file_paths <- path
  } else {
    stop("`path` points to neither a file nor a directory: \"", path, "\"")
  }

  required_elts <- c(
    "param", "control", "nwparam", "epi", "run", "coef.form", "num.nw"
  )
  infos <- vector("list", length(file_paths))
  for (i in seq_along(file_paths)) {
    x <- readRDS(file_paths[i])
    if (!inherits(x, "netsim")) {
      stop("\"", file_paths[i], "\" is not a `netsim` object")
    }
    missing_elts <- setdiff(required_elts, names(x))
    if (length(missing_elts) > 0) {
      stop(
        "\"", file_paths[i], "\" is missing the restart element(s): ",
        paste(missing_elts, collapse = ", ")
      )
    }
    # The check `netsim` makes when it starts
    check_restart_point_version(
      attr(x, "version"),
      paste0("\"", file_paths[i], "\"")
    )
    n_runs <- length(x$run)
    if (is_pool_dir && n_runs != 1) {
      stop(
        "\"", file_paths[i], "\" holds ", n_runs, " simulations. ",
        "Restart pool elements must hold a single one"
      )
    }

    infos[[i]] <- data.frame(
      pool_num = i,
      file_path = as.character(file_paths[i]),
      n_runs = n_runs,
      nsteps = x$control$nsteps,
      epimodel_version = as.character(attr(x, "version"))
    )
  }
  infos <- do.call(rbind, infos)

  if (length(unique(infos$nsteps)) > 1) {
    elts_by_steps <- split(infos$pool_num, infos$nsteps)
    stop(
      "The pool elements must all hold the same number of time steps, as ",
      "`netsim` requires `control$start` to be that number + 1. Found: ",
      paste0(
        names(elts_by_steps), " step(s) for element(s) ",
        vapply(elts_by_steps, paste, character(1), collapse = ", "),
        collapse = "; "
      )
    )
  }

  if (!is_pool_dir && infos$n_runs > 1) {
    warning(
      "\"", path, "\" holds ", infos$n_runs, " simulations. ",
      "`netsim_path_wrapper` restarts every replicate from the first one. ",
      "Use `make_restart_pool` to make a pool of these simulations."
    )
  }

  invisible(infos)
}

#' Run Single Simulations From a Restart Pool or a Single Object
#'
#' Runs one [netsim()] call with a single simulation for each element of
#' `sim_nums`, in parallel with `future.apply::future_lapply` (using the
#' current `future` plan).
#'
#' If `path_to_x` is a restart pool directory (one holding a `pool_index.csv`,
#' see [make_restart_pool()]) of `N` elements, simulation `k` starts from pool
#' element `(k - 1) %% N + 1`: the pool is reused in order when there are more
#' simulations than pool elements; the path of the restart point used is saved
#' in `$restart_point`. Before any simulation starts, the elements used are
#' checked to have been made with the installed major.minor version of
#' EpiModel, according to `pool_index.csv`. `netsim` checks the version of
#' each restart point again when it starts.
#' If `path_to_x` is a file, every simulation
#' starts from it (from its first simulation if it holds several).
#'
#' `control$nsims`, `control$ncores` and `control$future.use.plan` are
#' overridden. As `netsim` does for several simulations, `control$verbose` is
#' set to `FALSE` when `sim_nums` has more than one element.
#' If `control$.checkpoint.dir` is set, each simulation
#' checkpoints in its own `<.checkpoint.dir>/sim_<k>` sub directory.
#'
#' Use [validate_restart_pool()] to check a pool before using it.
#'
#' @param path_to_x Path to a restart pool directory (holding a
#'   `pool_index.csv`), or to a single file saved with `saveRDS` (a fitted
#'   network model or a restart point, see the `x` argument to [netsim()]).
#' @param sim_nums The simulation numbers to run, used to pick the pool
#'   elements.
#' @inheritParams netsim
#'
#' @return A list of `netsim` objects with one simulation each, in the order
#'   of `sim_nums`.
#'
#' @seealso [make_restart_pool()], [validate_restart_pool()]
#'
#' @export
netsim_path_wrapper <- function(path_to_x, param, init, control, sim_nums) {
  control$nsims <- 1
  control$ncores <- 1
  control$future.use.plan <- FALSE
  # `netsim` turns `verbose` off for several simulations. Here each one is its
  # own `netsim` call, so do it for the whole set
  if (length(sim_nums) > 1) control$verbose <- FALSE

  if (dir.exists(path_to_x)) {
    # case: a restart pool, its size is given by its `pool_index.csv`.
    # `validate_restart_pool` makes the full checks
    pool_index <- read_restart_pool_index(path_to_x)
    start_nums <- (sim_nums - 1) %% nrow(pool_index) + 1
    # First check of the EpiModel version of the elements used, before
    # starting any simulation. `netsim` checks it again when it starts
    for (n in unique(start_nums)) {
      check_restart_point_version(
        pool_index$epimodel_version[match(n, pool_index$pool_num)],
        paste0("Element ", n, " of the restart pool at \"", path_to_x, "\"")
      )
    }
    orig_paths <- as.character(file.path(path_to_x, paste0(start_nums, ".rds")))
  } else if (file.exists(path_to_x)) {
    # case: it's a file
    orig_paths <- rep(path_to_x, length(sim_nums))
  } else {
    stop("`path_to_x` points to neither of file or directory.")
  }

  checkpoint_dir <- control[[".checkpoint.dir"]]
  # Each `netsim` call is sequential (`ncores = 1`), so it never makes the
  # #1104 check, but the simulations run in the current plan's workers
  warn_worker_globals("netsim_path_wrapper", param = param, init = init, control = control)
  sim_list <- future.apply::future_lapply(
    seq_along(orig_paths),
    function(i) {
      if (!is.null(checkpoint_dir)) {
        control[[".checkpoint.dir"]] <-
          paste0(checkpoint_dir, "/sim_", sim_nums[i])
      }
      sim <- netsim(readRDS(orig_paths[i]), param, init, control)

      # If `sim` is a `netsim` - raw.output = FALSE
      if (inherits(sim, "netsim")) {
        # `merge.netsim` requires identical controls
        sim$control[[".checkpoint.dir"]] <- checkpoint_dir
        sim$restart_point <- list(sim1 = orig_paths[i])
      }
      sim
    },
    future.seed = TRUE
  )
  sim_list
}
