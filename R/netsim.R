
#' @title Stochastic Network Models
#'
#' @description Simulates stochastic network epidemic models for infectious
#'              disease. This is typically the final step in the network
#'              modeling pipeline, after network estimation with [`netest`]
#'              and model diagnostics with [`netdx`].
#'
#' @param x If `control$start == 1`, either a fitted network model object
#'        of class `netest` or a list of such objects. If
#'        `control$start > 1`, an object of class `netsim`. When
#'        multiple networks are used (multi-layer models), pass a list of
#'        `netest` objects, one per network layer; the node sets (including
#'        network size and nodal attributes) are assumed to be the same for
#'        all networks.
#' @param param Model parameters, as an object of class [`param.net`].
#'        Includes transmission probability (`inf.prob`), act rate
#'        (`act.rate`), recovery rate (`rec.rate`), and demographic rates
#'        for models with vital dynamics. Custom parameters for extended
#'        models may also be passed through `param.net`.
#' @param init Initial conditions, as an object of class [`init.net`].
#'        Specifies the initial number of infected nodes (`i.num`) and,
#'        for SIR models, recovered nodes (`r.num`). For two-group models,
#'        the corresponding `.g2` parameters are also required.
#' @param control Control settings, as an object of class [`control.net`].
#'        Key settings include `type` (disease model: `"SI"`, `"SIR"`, or
#'        `"SIS"`), `nsteps` (number of time steps), `nsims` (number of
#'        simulations), `tergmLite` (lightweight mode for performance), and
#'        `resimulate.network` (required for models with vital dynamics).
#'        For extended models, custom module functions are also passed here.
#'
#' @details
#' Stochastic network models explicitly represent phenomena within and across
#' edges (pairs of nodes that remain connected) over time. This enables edges to
#' have duration, allowing for repeated transmission-related acts within the
#' same dyad, specification of edge formation and dissolution rates, control
#' over the temporal sequencing of multiple edges, and specification of
#' network-level features. A detailed description of these models, along with
#' examples, is found in the [Network Modeling for Epidemics](https://epimodel.github.io/sismid/)
#' course materials.
#'
#' The `netsim` function performs modeling of both the base model types
#' and original models. Base model types include one-group and two-group models
#' with disease types for Susceptible-Infected (SI), Susceptible-Infected-Recovered (SIR),
#' and Susceptible-Infected-Susceptible (SIS).
#'
#' Original models may be parameterized by writing new process modules that
#' either take the place of existing modules (for example, disease recovery), or
#' supplement the set of existing processes with a new one contained in a new
#' module. This functionality is documented in the
#' [Extending EpiModel](https://epimodel.github.io/sismid/9_extending/mod9-Intro.html)
#' section of the [Network Modeling for Epidemics](https://epimodel.github.io/sismid/)
#' course materials. The list of modules within `netsim` available for
#' modification is listed in [modules.net()].
#'
#' @section Module Pipeline:
#' At each time step, `netsim` executes a series of modules in sequence.
#' For base models, the default pipeline is:
#'
#'  1. **Network resimulation** ([`resim_nets`]): updates the network
#'     structure (if `resimulate.network = TRUE`).
#'  2. **Infection** ([`infection.net`]): simulates disease transmission
#'     across discordant edges (where one partner is susceptible and the
#'     other is infected).
#'  3. **Recovery** ([`recovery.net`]): simulates recovery from infection
#'     (SIR and SIS models only).
#'  4. **Departures** ([`departures.net`]): simulates node exits from the
#'     network (if vital dynamics are enabled).
#'  5. **Arrivals** ([`arrivals.net`]): simulates new node entries into the
#'     network (if vital dynamics are enabled).
#'  6. **Prevalence** ([`prevalence.net`]): calculates summary statistics.
#'
#' See [modules.net()] for full details on each module.
#'
#' @section Performance and tergmLite:
#' Setting `tergmLite = TRUE` in [`control.net`] uses a lightweight network
#' representation (`networkLite`) that is substantially faster than the full
#' `networkDynamic` representation, often by a factor of 20--50x. This is
#' recommended for large networks or when running many simulations. The
#' trade-off is that full `networkDynamic` objects are not available for
#' post-hoc analysis; use the cumulative edgelist
#' (`cumulative.edgelist = TRUE`) instead to track partnership histories.
#'
#' The `resimulate.network` control must be set to `TRUE` when demographic
#' processes (arrivals and departures) change the network composition over
#' time. Without it, the network structure evolves independently of the
#' epidemic and demographic dynamics.
#'
#' @section Multi-Network Models:
#' For models with multiple overlapping network layers (e.g., sexual and
#' needle-sharing networks), pass a list of [`netest`] objects to the `x`
#' argument, one per network layer. Each layer has its own
#' formation/dissolution dynamics but shares the same node set. See the
#' `multilayer` documentation and `test-multinets.R` for examples.
#'
#' @section Restarting and Checkpointing:
#' Simulations can be checkpointed and restarted if interrupted. Set
#' `.checkpoint.steps` and `.checkpoint.dir` in [`control.net`] to enable
#' automatic checkpointing. To restart a simulation from a prior `netsim`
#' output, pass the `netsim` object as `x` and set `control$start` to one
#' greater than the final time step of the prior simulation. See the
#' Checkpointing Simulations section of [`control.net`] for full details.
#'
#' @return
#' A list of class `netsim` with the following elements:
#'
#'  * **param:** the epidemic parameters passed into the model through
#'        `param`, with additional parameters added as necessary.
#'  * **control:** the control settings passed into the model through
#'        `control`, with additional controls added as necessary.
#'  * **epi:** a list of data frames, one for each epidemiological
#'        output from the model. Outputs for base models always include the
#'        size of each compartment, as well as flows in, out of, and between
#'        compartments.
#'  * **stats:** a list containing two sublists, `nwstats` for any
#'        network statistics saved in the simulation, and `transmat` for
#'        the transmission matrix saved in the simulation. See
#'        [control.net()] for further
#'        details.
#'  * **network:** a list of lists of `networkDynamic` or
#'        `networkLite` objects, with one list of objects for each model
#'        simulation.
#'
#' If `control$raw.output == TRUE`: A list of the raw (pre-processed)
#' `netsim_dat` objects, for use in simulation continuation.
#'
#' The `epi` data can be extracted as a data frame with
#' [as.data.frame.netsim()], with options for per-simulation values
#' (`out = "vals"`), means (`out = "mean"`), standard deviations
#' (`out = "sd"`), or quantiles (`out = "qnt"`).
#'
#' @references
#' Jenness SM, Goodreau SM and Morris M. EpiModel: An R Package for Mathematical
#' Modeling of Infectious Disease over Networks. Journal of Statistical
#' Software. 2018; 84(8): 1-47.
#'
#' @seealso Estimate the network model with [`netest`] and diagnose model
#'          fit with [`netdx`] before running simulations. Extract model
#'          results with [as.data.frame.netsim()]. Summarize the
#'          time-specific model results with [summary.netsim()]. Plot the
#'          model results with [plot.netsim()]. Extract the network with
#'          [`get_network`], the transmission matrix with [`get_transmat`],
#'          and derive new epi statistics with [`mutate_epi`].
#'
#' @keywords model
#' @export
#'
#' @examples
#' \dontrun{
#' ## Example 1: SI Model without Network Feedback
#' # Network model estimation
#' nw <- network_initialize(n = 100)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#' est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' # Epidemic model
#' param <- param.net(inf.prob = 0.3)
#' init <- init.net(i.num = 10)
#' control <- control.net(type = "SI", nsteps = 100, nsims = 5, verbose.int = 0)
#' mod1 <- netsim(est1, param, init, control)
#'
#' # Print, plot, and summarize the results
#' mod1
#' plot(mod1)
#' summary(mod1, at = 50)
#'
#' ## Example 2: SIR Model with Network Feedback (Two-Group)
#' nw <- network_initialize(n = 100)
#' nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
#' formation <- ~edges
#' target.stats <- 50
#'
#' # Recalculate dissolution coefficient with departure rate
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20,
#'                                d.rate = 0.0021)
#'
#' # Reestimate the model with new coefficient
#' est2 <- netest(nw, formation, target.stats, coef.diss)
#'
#' # Set parameters to include demographic rates
#' param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.15,
#'                    rec.rate = 0.02, rec.rate.g2 = 0.02,
#'                    a.rate = 0.002, a.rate.g2 = NA,
#'                    ds.rate = 0.001, ds.rate.g2 = 0.001,
#'                    di.rate = 0.001, di.rate.g2 = 0.001,
#'                    dr.rate = 0.001, dr.rate.g2 = 0.001)
#' init <- init.net(i.num = 10, i.num.g2 = 10,
#'                  r.num = 0, r.num.g2 = 0)
#' control <- control.net(type = "SIR", nsteps = 100, nsims = 5,
#'                        resimulate.network = TRUE, tergmLite = TRUE)
#'
#' # Simulate the model with new network fit
#' mod2 <- netsim(est2, param, init, control)
#'
#' # Print, plot, and summarize the results
#' mod2
#' plot(mod2)
#' summary(mod2, at = 40)
#'
#' ## Example 3: Post-Simulation Analysis
#' # Extract epi data as a data frame
#' as.data.frame(mod1)
#' as.data.frame(mod1, out = "mean")
#'
#' # Extract the transmission matrix from simulation 1
#' get_transmat(mod1, sim = 1)
#'
#' # Derive new epi statistics and plot them
#' mod1 <- mutate_epi(mod1, prev = i.num / num)
#' plot(mod1, y = "prev")
#' }
#'
netsim <- function(x, param, init, control) {
  crosscheck.net(x, param, init, control)
  control <- netsim_validate_control(control)
  if (!is.null(control[["verbose.FUN"]]))
    control[["verbose.FUN"]](control, type = "startup")

  if (control$ncores == 1 && isFALSE(control$future.use.plan)) {
    dat_list <- lapply(
      seq_len(control$nsims),
      function(s) {
        netsim_initialize(x, param, init, control, s)
      }
    )
    dat_list <- Map(netsim_run, dat = dat_list, s = seq_along(dat_list))
  } else {
    if (inherits(control$future.use.plan, c("tweaked", "future"))) {
      with(future::plan(control$future.use.plan), local = TRUE)
    } else if (control$ncores > 1 && isFALSE(control$future.use.plan)) {
      with(future::plan("multisession", workers = control$ncores), local = TRUE)
    } # else if (isTRUE(control$future.use.plan)) - use plan defined by user
    dat_list <- future.apply::future_lapply(
      seq_len(control$nsims),
      function(s) {
        dat <- netsim_initialize(x, param, init, control, s)
        netsim_run(dat, s)
      },
      future.seed = TRUE
    )
  }

  out <- if (control$raw.output) dat_list else process_out.net(dat_list)

  if (!control[[".checkpoint.keep"]])
    netsim_clear_checkpoint(control)

  return(out)
}

# Ensure default values for control.net
#
netsim_validate_control <- function(control) {
  # Controls to be set to TRUE or FALSE if missing
  control_default_bool <- list(
    "TRUE" = c(
      ".checkpoint.compress"
    ),
    "FALSE" = c(
      "resimulate.network",
      "raw.output",
      "verbose",
      ".checkpoint.keep",
      ".traceback.on.error",
      ".dump.frame.on.error",
      "cumulative.edgelist",
      "future.use.plan",
      "tergmLite",
      "save.network",
      "save.diss.stats",
      "save.nwstats",
      "save.transmat",
      "save.network",
      "save.cumulative.edgelist",
      "save.run"
    )
  )

  if (is.null(control$save.other))
    control$save.other <- character(0)

  for (val in names(control_default_bool)) {
    for (flag in control_default_bool[[val]])
      if (is.null(control[[flag]])) control[[flag]] <- as.logical(val)
  }

  # truncate the cumulative edgelists to keep only active partnerships
  if (is.null(control$truncate.el.cuml))
    control$truncate.el.cuml <- 0

  if (is.null(control$start))
    control$start <- 1

  if (control$nsims == 1) {
    control$ncores <- 1
  } else {
    control$ncores <- min(parallel::detectCores(), control$ncores)
    control$verbose <- FALSE
  }

  control[[".checkpointed"]] <- netsim_is_checkpointed(control)
  if (!control[[".checkpointed"]]) {
    control[[".checkpoint.steps"]] <- Inf
  } else {
    control[[".checkpoint.steps"]] <- as.integer(control[[".checkpoint.steps"]])
  }

  if (control$nsteps < 1 || control$nsteps == Inf) {
    stop("`control$nsteps` must be positive and not infinite.")
  }

  return(control)
}

# Create the `dat` object or load it if checkpointed
#
netsim_initialize <- function(x, param, init, control, s = 1) {
  if (netsim_is_resume_checkpoint(control, s)) {
    dat <- netsim_load_checkpoint(control, s)
  } else {
    param <- generate_random_params(param, verbose = FALSE)
    dat <- control[["initialize.FUN"]](x, param, init, control, s)
    dat <- make_module_list(dat)
    if (get_control(dat, "start") != 1) {
      dat <- set_current_timestep(dat, get_control(dat, "start") - 1L)
    }
    if (get_control(dat, ".checkpointed"))
      netsim_save_checkpoint(dat, s)

    check_end_horizon_control(dat)
  }

  return(dat)
}

# Run a simulation from a initialized dat object
#
netsim_run <- function(dat, s = 1) {
  while (get_current_timestep(dat) < get_control(dat, "nsteps")) {
    steps_to_run <- get_control(dat, "nsteps") - get_current_timestep(dat)
    steps_to_run <- min(steps_to_run, get_control(dat, ".checkpoint.steps"))
    dat <- netsim_run_nsteps(dat, steps_to_run, s)

    if (get_control(dat, ".checkpointed"))
      netsim_save_checkpoint(dat, s)
  }

  return(dat)
}

# Run N simulation steps on a `dat` object
#
netsim_run_nsteps <- function(dat, nsteps, s) {
  for (n in seq_len(nsteps)) {
    dat <- increment_timestep(dat)
    dat <- netsim_run_modules(dat, s)
  }
  return(dat)
}

# Run all the modules of a `dat` object once
#
netsim_run_modules <- function(dat, s) {
  at <- get_current_timestep(dat)

  dat <- withCallingHandlers(
    expr = {
      current_mod <- "epimodel.internal"
      # Applies updaters, if any
      dat <- input_updater(dat)
      dat <- trigger_end_horizon(dat)

      modules <- get_modules(dat)

      for (i in seq_along(modules)) {
        current_mod <- names(modules)[i]
        dat <- modules[[i]](dat, at)
      }

      current_mod <- "epimodel.internal"
      # Run the user-provided trackers, if any
      dat <- epi_trackers(dat)

      ## Verbose module
      if (!is.null(get_control(dat, "verbose.FUN"))) {
        current_mod <- "verbose.FUN"
        do.call(get_control(dat, "verbose.FUN"),
                list(dat, type = "progress", s, at))
      }

      dat
    },
    message = function(e) message(netsim_cond_msg("MESSAGE", current_mod, at)),
    warning = function(e) message(netsim_cond_msg("WARNING", current_mod, at)),
    error = function(e) {
      message(netsim_cond_msg("ERROR", current_mod, at))
      netsim_error_logger(dat, s)
    }
  )

  return(dat)
}
