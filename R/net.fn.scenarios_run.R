#' Run EpiModel Network Simulations With Scenarios
#'
#' This function will run `n_rep` replications of each scenarios in
#' the `scenarios_list`. It runs them as multiple batches of up to
#' `n_cores` simulations at a time. The simfiles are then stored in the
#' `output_dir` folder and are named using the following pattern:
#' "sim__name_of_scenario__2.rds". Where the last number is the batch number
#' for this particular scenario. Each scenario is therefore run over
#' `ceiling(n_rep / n_cores)` batches.
#'
#' The batches run one after the other in the current R session, while the
#' simulations of a batch run in parallel, one `netsim` call per simulation
#' (see [netsim_path_wrapper()]), by default with
#' `future::plan("multisession", workers = n_cores)`. Each `netsim` call still
#' runs in a fresh `future` process, but the packages listed in `libraries` are
#' attached to the current session, not to the workers.
#' Using `control$future.use.plan <- future::tweak(<your plan>)` will bypass
#' the `multisession` setting.
#'
#' The EpiModelHPC package provides `slurmworkflow` step templates that run the
#' same code on a SLURM cluster.
#'
#' @param scenarios_list A list of scenarios to be run. Produced by the
#'   [create_scenario_list()] function
#' @param ... for compatibility reasons
#'
#' @inheritParams .run_batch
#' @inheritSection .run_batch Checkpointing
#'
#' @return `NULL`, invisibly. The simulations are saved in `output_dir`.
#'
#' @export
netsim_scenarios <- function(path_to_x, param, init, control,
                             scenarios_list, n_rep, n_cores,
                             output_dir, libraries = NULL, ...) {
  p_list <- .scenarios_setup(
    path_to_x = path_to_x,
    param = param,
    init = init,
    control = control,
    scenarios_list = scenarios_list,
    n_rep = n_rep,
    n_cores = n_cores,
    output_dir = output_dir,
    libraries = libraries
  )

  for (i in seq_along(p_list$scenarios_list)) {
    args <- c(
      list(
        scenario = p_list$scenarios_list[[i]],
        batch_num = p_list$batchs_list[[i]]
      ),
      p_list$MoreArgs
    )
    do.call(.run_batch, args)
  }

  invisible(NULL)
}

#' Set Up the Arguments to Run a Batch of Scenario Simulations
#'
#' Exported for the EpiModelHPC `slurmworkflow` step templates, which run the
#' same code on a SLURM cluster; it is not meant to be called directly.
#'
#' Expands the scenarios over their batches and builds the `MoreArgs` list
#' passed to [.run_batch()] by [netsim_scenarios()], checking the
#' `control` settings the scenario runners support.
#'
#' @inheritParams netsim_scenarios
#'
#' @return a list of arguments for `.run_batch`
#'
#' @keywords internal
#' @export
.scenarios_setup <- function(path_to_x, param, init, control,
                             scenarios_list, n_rep, n_cores,
                             output_dir, libraries) {
  if (is.null(scenarios_list)) {
    scenarios_list <- data.frame(.at = 0, .scenario.id = "empty_scenario")
    scenarios_list <- create_scenario_list(scenarios_list)
  }

  n_batch <- ceiling(n_rep / n_cores)
  batchs_list <- rep(seq_len(n_batch), length(scenarios_list))
  scenarios_list <- rep(scenarios_list, each = n_batch)

  if (isTRUE(control$future.use.plan)) {
    stop(
      "The scenario runners do not accept `control$future.use.plan = TRUE`. ",
      "Use `control$future.use.plan <- future::tweak(<your plan>)` to set up ",
      "a custom plan."
    )
  }

  if (isTRUE(control$raw.output)) {
    stop(
      "`raw.output = TRUE` is not supported by the scenario runners: each ",
      "simulation of a batch is run as its own `netsim()` call and the batch ",
      "is rebuilt with `merge.netsim()`, which needs `netsim` objects. ",
      "Set `control$raw.output <- FALSE`."
    )
  }

  list(
    scenarios_list = scenarios_list,
    batchs_list = batchs_list,
    MoreArgs = list(
      path_to_x = path_to_x,
      param = param,
      init = init,
      control = control,
      libraries = libraries,
      output_dir = output_dir,
      n_batch = n_batch,
      n_rep = n_rep,
      n_cores = n_cores
    )
  )
}

#' Build the `future` Plan to Run the Simulations of a Batch
#'
#' Exported for the EpiModelHPC `slurmworkflow` step templates, which run the
#' same code on a SLURM cluster; it is not meant to be called directly.
#'
#' Only builds the plan: `with(future::plan(<plan>), local = TRUE)` must be
#' called by the function running the simulations, as it resets the plan when
#' its caller exits.
#'
#' `strategy` is a function rather than a name, as `future::tweak()` looks a
#' name up from the calling frame (this namespace, then the search path) before
#' the `future` namespace.
#'
#' @param control A `control.net` list. Its `future.use.plan` element is used
#'   as the plan when it holds one set with `future::tweak()`.
#' @param strategy A `future` backend function, e.g. `future::multicore`.
#' @param n_cores The number of workers to use when `control$future.use.plan`
#'   holds no plan.
#'
#' @return A `future` plan: `control$future.use.plan` when it is a
#'   `future::tweak()` plan, else `future::tweak(strategy, workers = n_cores)`.
#'
#' @keywords internal
#' @export
.batch_plan <- function(control, strategy, n_cores) {
  if (inherits(control$future.use.plan, c("tweaked", "future"))) {
    control$future.use.plan
  } else {
    future::tweak(strategy, workers = n_cores)
  }
}

#' Run and Save One Batch of `netsim` Simulations
#'
#' Exported for the EpiModelHPC `slurmworkflow` step templates, which run the
#' same code on a SLURM cluster; it is not meant to be called directly.
#'
#' This inner function is called by [netsim_scenarios()]. It runs the
#' simulations of batch `batch_num` with [netsim_path_wrapper()], one `netsim`
#' call per simulation, merges them and saves them as
#' "sim__<scenario id>__<batch_num>.rds" in `output_dir`.
#' The simulations of batch `b` are numbered `(b - 1) * n_cores + 1` to
#' `b * n_cores` over all the batches of a scenario. With a restart pool of
#' size `N`, simulation `k` starts from pool element `(k - 1) %% N + 1`, the
#' same for every scenario.
#' The simulations are run with the `future::tweak()` plan in
#' `control$future.use.plan`, or with
#' `future::plan("multisession", workers = n_cores)` if there is none. The
#' `slurmworkflow` step templates set it to `multicore` before calling this
#' function.
#'
#' @param path_to_x Path to a fitted network model or a restart point saved
#'   with `saveRDS` (See the `x` argument to the [netsim()] function), or to a
#'   restart pool directory (see [make_restart_pool()]).
#' @param scenario A single "`EpiModel` scenario" to be used in the simulation
#' @param batch_num The batch number, calculated from the number of replications
#'   and CPUs required.
#' @param n_batch The number of batches to be run `ceiling(n_rep / n_cores)`.
#' @param n_rep The number of replication to be run for each scenario.
#' @param n_cores The number of CPUs on which the simulations will be run.
#' @param output_dir The folder where the simulation files are to be stored.
#' @param libraries A character vector containing the name of the libraries
#'   required for the model to run. (e.g. EpiModelHIV or EpiModelCOVID)
#' @inheritParams netsim
#'
#' @section Checkpointing:
#' This function takes care of editing `.checkpoint.dir` to create unique sub
#' directories for each scenario, batch and simulation:
#' "<.checkpoint.dir>/sim__<scenario id>__<batch_num>/sim_<k>". The
#' `control.net` way of setting up checkpoints can be used
#' transparently. Unless `.checkpoint.keep` is `TRUE`, the batch directory is
#' removed once the batch file is saved.
#'
#' @keywords internal
#' @export
.run_batch <- function(scenario, batch_num,
                       path_to_x, param, init, control,
                       libraries, output_dir,
                       n_batch, n_rep, n_cores) {
  start_time <- Sys.time()
  lapply(libraries, function(l) library(l, character.only = TRUE))

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  sim_nums_offset <- (batch_num - 1) * n_cores
  # On last batch, adjust the number of simulation to be run
  if (batch_num == n_batch) n_cores <- n_rep - n_cores * (n_batch - 1)
  sim_nums <- sim_nums_offset + seq_len(n_cores)

  param_sc <- use_scenario(param, scenario)

  if (!is.null(control[[".checkpoint.dir"]])) {
    control[[".checkpoint.dir"]] <- paste0(
      control[[".checkpoint.dir"]], "/sim__", scenario[["id"]], "__", batch_num
    )
  }

  # The step templates set a `multicore` plan before reaching here. After the
  # last batch adjustment of `n_cores`, to not start idle workers
  sim_plan <- .batch_plan(control, future::multisession, n_cores)

  message("Starting simulation for scenario: ", scenario[["id"]])
  message("Batch number: ", batch_num, " / ", n_batch)
  with(future::plan(sim_plan), local = TRUE)
  # `keep.cumulative.edgelist` is the only `merge.netsim` default that drops an
  # element a single `netsim` call with several simulations keeps
  sim <- netsim_path_wrapper(
    path_to_x = path_to_x,
    param = param_sc,
    init = init,
    control = control,
    sim_nums = sim_nums
  ) |>
    Reduce(f = \(x, y) merge(x, y, keep.cumulative.edgelist = TRUE))

  file_name <- paste0("sim__", scenario[["id"]], "__", batch_num, ".rds")
  message("Saving simulation in file: ", file_name)
  saveRDS(sim, file.path(output_dir, file_name))

  # Each `netsim` call only clears its own `sim_<k>` checkpoint directory
  batch_ckpt_dir <- control[[".checkpoint.dir"]]
  if (!is.null(batch_ckpt_dir) && !isTRUE(control[[".checkpoint.keep"]])) {
    unlink(batch_ckpt_dir, recursive = TRUE)
  }

  message("Done in: ", format(Sys.time() - start_time))
  invisible(NULL)
}

#' Helper function to access the file name elements of scenarios
#'
#' This function returns the list of simulation files and the corresponding
#' scenario name and batch number present in a given directory. It is meant to
#' be used after [netsim_scenarios()]. Only files named
#' `sim__<scenario>__<batch>.rds` are listed: other files, and directories, are
#' ignored.
#'
#' @param scenario_dir the directory where `netsim_scenarios` saved it's
#'   simulations.
#'
#' @return a `tibble` with three columns: `file_path` - the full paths of
#' the simulation file, `scenario_name` the associated scenario name,
#' `batch_number` the associated batch number. Sorted by scenario name, then
#' batch number (`2` before `10`).
#'
#' @export
get_scenarios_batches_infos <- function(scenario_dir) {
  if (!dir.exists(scenario_dir)) {
    stop(
      "\"", scenario_dir, "\" is not an existing directory: check the path ",
      "where `netsim_scenarios` saved the simulations"
    )
  }
  file_names <- list.files(scenario_dir, pattern = "^sim__.+__[0-9]+\\.rds$")
  file_names <- file_names[!dir.exists(file.path(scenario_dir, file_names))]

  parts <- regmatches(
    file_names,
    regexec("^sim__(.+)__([0-9]+)\\.rds$", file_names)
  )

  dplyr::tibble(
    file_path = file.path(scenario_dir, file_names),
    scenario_name = vapply(parts, `[`, character(1), 2),
    batch_number = as.integer(vapply(parts, `[`, character(1), 3))
  ) |>
    dplyr::arrange(.data$scenario_name, .data$batch_number)
}

#' Create a Single Sim File per Scenarios Using the Files From
#' `netsim_scenarios`
#'
#' Merges the batch files of each scenario into a single `netsim` object,
#' saved as `merged__<scenario>.rds` in `output_dir`.
#'
#' @param sim_dir The folder where the simulation files are to be stored.
#' @param output_dir The folder where the merged files will be stored.
#' @param truncate.at Time step at which to left-truncate the time series.
#'
#' @inheritParams merge.netsim
#'
#' @return `NULL`, invisibly. The merged files are saved in `output_dir`.
#'
#' @export
merge_netsim_scenarios <- function(sim_dir, output_dir,
                                   keep.transmat = TRUE, keep.network = TRUE,
                                   keep.nwstats = TRUE, keep.other = TRUE,
                                   param.error = FALSE, keep.diss.stats = TRUE,
                                   truncate.at = NULL, keep.run = TRUE,
                                   keep.cumulative.edgelist = FALSE,
                                   keep.attr.history = TRUE) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  batches_infos <- get_scenarios_batches_infos(sim_dir)

  future.apply::future_lapply(
    unique(batches_infos$scenario_name),
    function(scenario) {
      scenario_infos <- dplyr::filter(
        batches_infos,
        .data$scenario_name == scenario
      )
      file_paths <- scenario_infos$file_path
      merged <- NULL
      for (j in seq_along(file_paths)) {
        current <- readRDS(file_paths[j])
        if (!is.null(truncate.at)) {
          current <- truncate_sim(current, truncate.at)
        }

        if (j == 1) {
          merged <- current
        } else {
          merged <- merge(
            merged,
            current,
            keep.transmat = keep.transmat,
            keep.network = keep.network,
            keep.nwstats = keep.nwstats,
            keep.other = keep.other,
            param.error = param.error,
            keep.diss.stats = keep.diss.stats,
            keep.run = keep.run,
            keep.cumulative.edgelist = keep.cumulative.edgelist,
            keep.attr.history = keep.attr.history
          )
        }
      }

      # saved once per scenario, after its last batch is merged
      saveRDS(
        merged,
        file.path(output_dir, paste0("merged__", scenario, ".rds"))
      )
      invisible(NULL)
    }
  )

  invisible(NULL)
}

#' Create a Single Sim File per Scenarios Using the Files From
#' `netsim_scenarios`
#'
#' Converts the batch files of each scenario into a single `data.frame` of
#' simulation output, saved as `df__<scenario>.rds` in `output_dir`. The
#' `batch_number`, `sim` and `time` columns always come first.
#'
#' @param steps_to_keep Numbers of time steps add the end of the simulation to
#'   keep in the `data.frame`.
#' @param cols <tidy-select> columns to keep in the `data.frame`. By default all
#'    columns are kept. And in any case, the `batch_number`, `sim` and `time`
#'    are always kept.
#'
#' @inheritParams merge_netsim_scenarios
#'
#' @return `NULL`, invisibly. The merged files are saved in `output_dir`.
#'
#' @export
merge_netsim_scenarios_tibble <- function(sim_dir, output_dir, steps_to_keep,
                                          cols = dplyr::everything()) {
  expr <- rlang::enquo(cols)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  batches_infos <- get_scenarios_batches_infos(sim_dir)

  for (scenario in unique(batches_infos$scenario_name)) {
    scenario_infos <- dplyr::filter(
      batches_infos,
      .data$scenario_name == scenario
    )

    df_list <- future.apply::future_lapply(
      seq_len(nrow(scenario_infos)),
      function(i) {
        sc_inf <- scenario_infos[i, ]
        d <- readRDS(sc_inf$file_path) |>
          dplyr::as_tibble() |>
          dplyr::filter(.data$time >= max(.data$time) - steps_to_keep)

        d_fix <- dplyr::select(d, "sim", "time")
        d_var <- dplyr::select(d, -c("sim", "time"))

        d_var <- dplyr::select(d_var, !!expr)

        dplyr::bind_cols(d_fix, d_var) |>
          dplyr::mutate(batch_number = as.integer(sc_inf$batch_number)) |>
          dplyr::select("batch_number", "sim", "time", dplyr::everything())
      }
    )
    df_sc <- dplyr::bind_rows(df_list) |>
      dplyr::mutate(
        tmp_max_sim = max(.data$sim),
        sim_number = .data$sim,
        sim = .data$tmp_max_sim * (.data$batch_number - 1) + .data$sim
      ) |>
      dplyr::select(-c("tmp_max_sim"))
    saveRDS(df_sc, file.path(output_dir, paste0("df__", scenario, ".rds")))
  }

  invisible(NULL)
}

#' Helper function to access the infos on merged scenarios `data.frame`
#'
#' This function returns the list of scenario tibble files and the corresponding
#' scenario name present in a given directory. It is meant to
#' be used after [merge_netsim_scenarios_tibble()]. Only files named
#' `df__<scenario>.rds` are listed: other files, and directories, are ignored.
#'
#' @param scenario_dir the directory where `merge_netsim_scenarios_tibble` saved
#' the merged tibbles.
#'
#' @return a `tibble` with two columns: `file_path` - the full path of
#' the scenario tibble file and `scenario_name` the associated scenario name.
#'
#' @export
get_scenarios_tibble_infos <- function(scenario_dir) {
  if (!dir.exists(scenario_dir)) {
    stop(
      "\"", scenario_dir, "\" is not an existing directory: check the path ",
      "where `merge_netsim_scenarios_tibble` saved the tibbles"
    )
  }
  file_names <- list.files(scenario_dir, pattern = "^df__.+\\.rds$")
  file_names <- file_names[!dir.exists(file.path(scenario_dir, file_names))]

  parts <- regmatches(
    file_names,
    regexec("^df__(.+)\\.rds$", file_names)
  )

  dplyr::tibble(
    file_path = file.path(scenario_dir, file_names),
    scenario_name = vapply(parts, `[`, character(1), 2)
  )
}
