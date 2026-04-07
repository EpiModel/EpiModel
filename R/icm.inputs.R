
#' @title Epidemic Parameters for Stochastic Individual Contact Models
#'
#' @description Sets the epidemic parameters for stochastic individual contact
#'              models simulated with `icm`.
#'
#' @inheritParams param.dcm
#'
#' @details
#' `param.icm` sets the epidemic parameters for the stochastic individual
#' contact models simulated with the [icm()] function. Models
#' may use the base types, for which these parameters are used, or new process
#' modules which may use these parameters (but not necessarily).
#'
#' For base models, the model specification will be chosen as a result of
#' the model parameters entered here and the control settings in
#' [control.icm()]. One-group and two-group models are available,
#' where the former assumes a homogeneous mixing in the population and the
#' latter assumes some form of heterogeneous mixing between two distinct
#' partitions in the population (e.g., men and women). Specifying any group two
#' parameters (those with a `.g2`) implies the simulation of a two-group
#' model. All the parameters for a desired model type must be specified, even if
#' they are zero.
#'
#' @section Act Balancing:
#' In two-group models, a balance between the number of acts for group 1 members
#' and those for group 2 members must be maintained. With purely heterogeneous
#' mixing, the product of one group size and act rate must equal the product of
#' the other group size and act rate: \eqn{N_1 \alpha_1 = N_2 \alpha_2}, where
#' \eqn{N_i} is the group size and \eqn{\alpha_i} the group-specific act rate
#' at time \eqn{t}. The `balance` parameter here specifies which group's
#' act rate should control the others with respect to balancing.
#'
#' @return An `EpiModel` object of class `param.icm`.
#'
#' @seealso Use [init.icm()] to specify the initial conditions and
#'          [control.icm()] to specify the control settings. Run the
#'          parameterized model with [icm()].
#'
#' @keywords parameterization
#'
#' @export
#'
param.icm <- function(inf.prob, inter.eff, inter.start, act.rate, rec.rate,
                      a.rate, ds.rate, di.rate, dr.rate, inf.prob.g2,
                      act.rate.g2, rec.rate.g2, a.rate.g2, ds.rate.g2,
                      di.rate.g2, dr.rate.g2, balance, ...) {

  # Get arguments
  p <- list()
  formal.args <- formals(sys.function())
  formal.args[["..."]] <- NULL
  for (arg in names(formal.args)) {
    if (as.logical(mget(arg) != "")) {
      p[arg] <- list(get(arg))
    }
  }
  dot.args <- list(...)
  names.dot.args <- names(dot.args)
  if (length(dot.args) > 0) {
    for (i in seq_along(dot.args)) {
      p[[names.dot.args[i]]] <- dot.args[[i]]
    }
  }

  if ("b.rate" %in% names.dot.args) {
    stop("The b.rate parameter has been removed. Use a.rate instead.",
         call. = FALSE)
  }
  if ("b.rate.g2" %in% names.dot.args) {
    stop("The b.rate.g2 parameter has been removed. Use a.rate.g2 instead.",
         call. = FALSE)
  }

  ## Defaults and checks
  if (is.null(p$act.rate)) {
    p$act.rate <- 1
  }
  p$vital <- ifelse(!is.null(p$a.rate) | !is.null(p$ds.rate) |
                      !is.null(p$di.rate) | !is.null(p$dr.rate), TRUE, FALSE)

  p$groups <- ifelse(any(grepl(".g2", names(p))) == TRUE, 2, 1)

  if (p$groups == 2 && (is.null(p$balance) ||
                          !(p$balance %in% c("g1", "g2")))) {
    stop("Specify balance=\"g1\" or balance=\"g2\" with 2-group models")
  }

  if (!is.null(p$inter.eff) && is.null(p$inter.start)) {
    p$inter.start <- 1
  }

  ## Output
  class(p) <- c("param.icm", "list")
  return(p)
}


#' @title Initial Conditions for Stochastic Individual Contact Models
#'
#' @description Sets the initial conditions for stochastic individual contact
#'              models simulated with `icm`.
#'
#' @param s.num Number of initial susceptible persons. For two-group models,
#'        this is the number of initial group 1 susceptible persons.
#' @param i.num Number of initial infected persons. For two-group models, this
#'        is the number of initial group 1 infected persons.
#' @param r.num Number of initial recovered persons. For two-group models, this
#'        is the number of initial group 1 recovered persons. This parameter is
#'        only used for the `SIR` model type.
#' @param s.num.g2 Number of initial susceptible persons in group 2. This
#'        parameter is only used for two-group models.
#' @param i.num.g2 Number of initial infected persons in group 2. This parameter
#'        is only used for two-group models.
#' @param r.num.g2 Number of initial recovered persons in group 2. This
#'        parameter is only used for two-group `SIR` models.
#' @param ... Additional initial conditions passed to model.
#'
#' @details
#' The initial conditions for a model solved with [icm()] should be
#' input into the `init.icm` function. This function handles initial
#' conditions for both base models and original models using new modules.
#'
#' @return An `EpiModel` object of class `init.icm`.
#'
#' @seealso Use [param.icm()] to specify model parameters and
#'          [control.icm()] to specify the control settings. Run the
#'          parameterized model with [icm()].
#'
#' @keywords parameterization
#'
#' @export
#'
init.icm <- function(s.num, i.num, r.num,
                     s.num.g2, i.num.g2, r.num.g2, ...) {

  # Get arguments
  p <- list()
  formal.args <- formals(sys.function())
  formal.args[["..."]] <- NULL
  for (arg in names(formal.args)) {
    if (as.logical(mget(arg) != "")) {
      p[arg] <- list(get(arg))
    }
  }
  dot.args <- list(...)
  names.dot.args <- names(dot.args)
  if (length(dot.args) > 0) {
    for (i in seq_along(dot.args)) {
      p[[names.dot.args[i]]] <- dot.args[[i]]
    }
  }


  ## Output
  class(p) <- c("init.icm", "list")
  return(p)
}


#' @title Control Settings for Stochastic Individual Contact Models
#'
#' @description Sets the controls for stochastic individual contact models
#'              simulated with [icm()].
#'
#' @param type Disease type to be modeled, with the choice of `"SI"` for
#'        Susceptible-Infected diseases, `"SIR"` for
#'        Susceptible-Infected-Recovered diseases, and `"SIS"` for
#'        Susceptible-Infected-Susceptible diseases.
#' @param nsteps Number of time steps to solve the model over. This must be a
#'        positive integer.
#' @param nsims Number of simulations to run.
#' @param verbose If `TRUE`, print model progress to the console.
#' @param verbose.int Time step interval for printing progress to console, where
#'        0 (the default) prints completion status of entire simulation and
#'        positive integer `x` prints progress after every `x` time
#'        steps.
#'
#' @details
#' `control.icm` sets the required control settings for any stochastic
#' individual contact model solved with the [icm()] function. ICM simulations
#' use the built-in SI, SIR, and SIS disease types only. The `type` argument
#' is required and has no default. For custom or extension epidemic models, use
#' the network model class via [control.net()] instead.
#'
#' @return An `EpiModel` object of class `control.icm`.
#'
#' @seealso Use [param.icm()] to specify model parameters and
#'          [init.icm()] to specify the initial conditions. Run the
#'          parameterized model with [icm()].
#'
#' @keywords parameterization
#'
#' @export
#'
control.icm <- function(type, nsteps, nsims = 1,
                        verbose = FALSE, verbose.int = 0) {

  # Get arguments
  p <- list()
  formal.args <- formals(sys.function())
  for (arg in names(formal.args)) {
    if (as.logical(mget(arg) != "")) {
      p[arg] <- list(get(arg))
    }
  }

  ## Defaults and checks
  if (is.null(p$type) || !(p$type %in% c("SI", "SIS", "SIR"))) {
    stop("Specify type as \"SI\", \"SIS\", or \"SIR\" ", call. = FALSE)
  }
  if (is.null(p$nsteps)) {
    stop("Specify nsteps", call. = FALSE)
  }


  ## Output
  p <- set.control.class("control.icm", p)
  return(p)
}


#' @title Cross Checking of Inputs for Stochastic Individual Contact Models
#'
#' @description This function checks that the three parameter lists from
#'              [param.icm()], [init.icm()], and
#'              [control.icm()] are consistent.
#'
#' @param param An `EpiModel` object of class [param.icm()].
#' @param init An `EpiModel` object of class [init.icm()].
#' @param control An `EpiModel` object of class [control.icm()].
#'
#' @return
#' This function returns no objects.
#'
#' @export
#' @keywords internal
#'
crosscheck.icm <- function(param, init, control) {
  check.control.class("icm", "EpiModel crosscheck.icm")

  ## Main class check
  if (!inherits(param, "param.icm")) {
    stop("param must be an object of class param.icm", call. = FALSE)
  }
  if (!inherits(init, "init.icm")) {
    stop("init must be an object of class init.icm", call. = FALSE)
  }
  if (!inherits(control, "control.icm")) {
    stop("control must be an object of class control.icm", call. = FALSE)
  }

  ## Check that rec.rate is supplied for SIR models
  if (control$type %in% c("SIR", "SIS")) {
    if (is.null(param$rec.rate)) {
      stop("Specify rec.rate in param.icm", call. = FALSE)
    }
    if (param$groups == 2 && is.null(param$rec.rate.g2)) {
      stop("Specify rec.rate.g2 in param.icm", call. = FALSE)
    }
  }

  ## Check that parameters and init are supplied for SIR models
  if (control$type == "SIR") {
    if (is.null(init$r.num)) {
      stop("Specify r.num in init.icm", call. = FALSE)
    }
    if (param$groups == 2 && is.null(init$r.num.g2)) {
      stop("Specify r.num.g2 in init.icm", call. = FALSE)
    }
  }

  ## Check that groups implied by init and params are consistent
  if (any(grepl(".g2", names(init))) == TRUE) {
    init.groups <- 2
  } else {
    init.groups <- 1
  }
  if (param$groups == 2 && init.groups == 1) {
    stop("Group 2 parameters specified in param.icm, but missing group 2, ",
         "initial states in init.icm", call. = FALSE)
  }
  if (param$groups == 1 && init.groups == 2) {
    stop("Group 2 initial stats specified in init.icm, but missing group 2 ",
         "parameters in param.icm", call. = FALSE)
  }

  ## Deprecated parameters
  if (!is.null(param$trans.rate)) {
    stop("The trans.rate parameter is deprecated. Use the inf.prob ",
         "parameter instead.", call. = FALSE)
  }
  if (!is.null(param$trans.rate.g2)) {
    stop("The trans.rate.g2 parameter is deprecated. Use the inf.prob.g2 ",
         "parameter instead.", call. = FALSE)
  }

  ## Assign built-in modules based on group parameter
  ## initialize.icm handles both 1-group and 2-group (no .bip variant)
  control[["initialize.FUN"]] <- initialize.icm
  bi.mods <- c("infection.FUN", "recovery.FUN",
               "departures.FUN", "arrivals.FUN", "prevalence.FUN")
  if (param$groups == 1) {
    for (mod in bi.mods) {
      control[[mod]] <- get(gsub(".FUN", ".icm", mod))
    }
  } else {
    for (mod in bi.mods) {
      control[[mod]] <- get(gsub(".FUN", ".icm.bip", mod))
    }
  }

  ## In-place assignment to update param and control
  on.exit(assign("param", param, pos = parent.frame()))
  on.exit(assign("control", control, pos = parent.frame()), add = TRUE)
}
