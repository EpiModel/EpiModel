
#' @title Epidemic Parameters for Stochastic Individual Contact Models
#'
#' @description Sets the epidemic parameters for stochastic individual contact
#'              models simulated with `icm`.
#'
#' @inheritParams param.dcm
#'
#' @details
#' `param.icm` sets the epidemic parameters for the stochastic individual
#' contact models simulated with the `\link{icm}` function. Models
#' may use the base types, for which these parameters are used, or new process
#' modules which may use these parameters (but not necessarily). A detailed
#' description of ICM parameterization for base models is found in the
#' \href{http://statnet.github.io/tut/BasicICMs.html}{Basic ICMs} tutorial.
#'
#' For base models, the model specification will be chosen as a result of
#' the model parameters entered here and the control settings in
#' `\link{control.icm}`. One-group and two-group models are available, where
#' the former assumes a homogenous mixing in the population and the latter
#' assumes a purely heterogenous mixing between two distinct partitions in the
#' population (e.g., men and women). Specifying any group two parameters (those
#' with a `.g2`) implies the simulation of a two-group model. All the
#' parameters for a desired model type must be specified, even if they are zero.
#'
#' @section Act Balancing:
#' In two-group models, a balance between the number of acts for group 1 members
#' and those for group 2 members must be maintained. With purely heterogenous
#' mixing, the product of one group size and act rate must equal the product of
#' the other group size and act rate: \eqn{N_1 \alpha_1 = N_2 \alpha_2}, where
#' \eqn{N_i} is the group size and \eqn{\alpha_i} the group-specific act rates
#' at time \eqn{t}. The `balance` parameter here specifies which group's
#' act rate should control the others with respect to balancing. See the
#' \href{http://statnet.github.io/tut/BasicDCMs.html}{Basic DCMs} tutorial.
#'
#' @section New Modules:
#' To build original models outside of the base models, new process modules
#' may be constructed to replace the existing modules or to supplement the existing
#' set. These are passed into the control settings in `\link{control.icm}`.
#' New modules may use either the existing model parameters named here, an
#' original set of parameters, or a combination of both. The `...` allows
#' the user to pass an arbitrary set of original model parameters into
#' `param.icm`. Whereas there are strict checks with default modules for
#' parameter validity, these checks are the user's responsibility with new modules.
#'
#' @seealso Use `\link{init.icm}` to specify the initial conditions and
#'          `\link{control.icm}` to specify the control settings. Run the
#'          parameterized model with `\link{icm}`.
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
    for (i in 1:length(dot.args)) {
      p[[names.dot.args[i]]] <- dot.args[[i]]
    }
  }

  if ("b.rate" %in% names.dot.args) {
    p$a.rate <- dot.args$b.rate
    message("EpiModel 1.7.0 onward renamed the birth rate parameter b.rate to a.rate. ",
            "See documentation for details.")
  }
  if ("b.rate.g2" %in% names.dot.args) {
    p$a.rate.g2 <- dot.args$b.rate.g2
    message("EpiModel 1.7.0 onward renamed the birth rate parameter b.rate.g2 to a.rate.g2. ",
            "See documentation for details.")
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
#' @param s.num Number of initial susceptible. For two-group models, this is
#'        the number of initial group 1 susceptible.
#' @param i.num Number of initial infected. For two-group models, this is the
#'        number of initial group 1 infected.
#' @param r.num Number of initial recovered. For two-group models, this is the
#'        number of initial group 1 recovered. This parameter is only used for
#'        the `SIR` model type.
#' @param s.num.g2 Number of initial susceptible in group 2. This parameter is
#'        only used for two-group models.
#' @param i.num.g2 Number of initial infected in group 2. This parameter is only
#'        used for two-group models.
#' @param r.num.g2 Number of initial recovered in group 2. This parameter is
#'        only used for two-group `SIR` models.
#' @param ... Additional initial conditions passed to model.
#'
#' @details
#' The initial conditions for a model solved with `\link{icm}` should be
#' input into the `init.icm` function. This function handles initial
#' conditions for both base models and original models using new modules. For
#' an overview of initial conditions for base ICM class models, consult the
#' \href{http://statnet.github.io/tut/BasicICMs.html}{Basic ICMs} tutorial.
#'
#' @seealso Use `\link{param.icm}` to specify model parameters and
#'          `\link{control.icm}` to specify the control settings. Run the
#'          parameterized model with `\link{icm}`.
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
    for (i in 1:length(dot.args)) {
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
#'              simulated with `\link{icm}`.
#'
#' @param type Disease type to be modeled, with the choice of `"SI"` for
#'        Susceptible-Infected diseases, `"SIR"` for
#'        Susceptible-Infected-Recovered diseases, and `"SIS"` for
#'        Susceptible-Infected-Susceptible diseases.
#' @param nsteps Number of time steps to solve the model over. This must be a
#'        positive integer.
#' @param nsims Number of simulations to run.
#' @param initialize.FUN Module to initialize the model at the outset, with the
#'        default function of `\link{initialize.icm}`.
#' @param infection.FUN Module to simulate disease infection, with the default
#'        function of `\link{infection.icm}`.
#' @param recovery.FUN Module to simulate disease recovery, with the default
#'        function of `\link{recovery.icm}`.
#' @param departures.FUN Module to simulate departures or exits, with the default
#'        function of `\link{departures.icm}`.
#' @param arrivals.FUN Module to simulate arrivals or entries, with the default
#'        function of `\link{arrivals.icm}`.
#' @param prevalence.FUN Module to calculate disease prevalence at each time step,
#'        with the default function of `\link{prevalence.icm}`.
#' @param verbose If `TRUE`, print model progress to the console.
#' @param verbose.int Time step interval for printing progress to console, where
#'        0 (the default) prints completion status of entire simulation and
#'        positive integer `x` prints progress after each `x` time
#'        steps.
#' @param skip.check If `TRUE`, skips the default error checking for the
#'        structure and consistency of the parameter values, initial conditions,
#'        and control settings before running base epidemic models. Setting
#'        this to `FALSE` is recommended when running models with new modules
#'        specified.
#' @param ... Additional control settings passed to model.
#'
#' @details
#' `control.icm` sets the required control settings for any stochastic
#' individual contact model solved with the `\link{icm}` function. Controls
#' are required for both base model types and when passing original process
#' modules. For an overview of control settings for base ICM class models,
#' consult the \href{http://statnet.github.io/tut/BasicICMs.html}{Basic ICMs}
#' tutorial. For all base models, the `type` argument is a necessary
#' parameter and it has no default.
#'
#' @section New Modules:
#' Base ICM models use a set of module functions that specify
#' how the individual agents in the population are subjected to infection, recovery,
#' demographics, and other processes. Core modules are those listed in the
#' `.FUN` arguments. For each module, there is a default function used in
#' the simulation. The default infection module, for example, is contained in
#' the `\link{infection.icm}` function.
#'
#' For original models, one may substitute replacement module functions for any of
#' the default functions. New modules may be added to the workflow at each time
#' step by passing a module function via the `...` argument.
#'
#' @seealso Use `\link{param.icm}` to specify model parameters and
#'          `\link{init.icm}` to specify the initial conditions. Run the
#'          parameterized model with `\link{icm}`.
#'
#' @keywords parameterization
#'
#' @export
#'
control.icm <- function(type, nsteps, nsims = 1,
                        initialize.FUN = initialize.icm,
                        infection.FUN = NULL, recovery.FUN = NULL,
                        departures.FUN = NULL, arrivals.FUN = NULL,
                        prevalence.FUN = NULL, verbose = FALSE,
                        verbose.int = 0, skip.check = FALSE, ...) {

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
    for (i in 1:length(dot.args)) {
      p[[names.dot.args[i]]] <- dot.args[[i]]
    }
  }

  if ("births.FUN" %in% names(dot.args)) {
    p$arrivals.FUN <- dot.args$births.FUN
    p$births.FUN <- dot.args$births.FUN <- NULL
    message("EpiModel 1.7.0 onward renamed the birth function births.FUN to arrivals.FUN. See documentation for details.")
  }
  if ("deaths.FUN" %in% names(dot.args)) {
    p$departures.FUN <- dot.args$deaths.FUN
    p$deaths.FUN <- dot.args$deaths.FUN <- NULL
    message("EpiModel 1.7.0 onward renamed the death function deaths.FUN to departures.FUN. See documentation for details.")
  }


  ## Module classification
  p$bi.mods <- grep(".FUN", names(formal.args), value = TRUE)
  p$user.mods <- grep(".FUN", names(dot.args), value = TRUE)


  ## Defaults and checks
  if (is.null(p$type) | !(p$type %in% c("SI", "SIS", "SIR"))) {
    stop("Specify type as \"SI\", \"SIS\", or \"SIR\" ", call. = FALSE)
  }
  if (is.null(p$nsteps)) {
    stop("Specify nsteps", call. = FALSE)
  }


  ## Output
  class(p) <- c("control.icm", "list")
  return(p)
}


#' @title Cross Checking of Inputs for Stochastic Individual Contact Models
#'
#' @description This function checks that the three parameter lists from
#'              `\link{param.icm}`, `\link{init.icm}`, and
#'              `\link{control.icm}` are consistent.
#'
#' @param param An `EpiModel` object of class `\link{param.icm}`.
#' @param init An `EpiModel` object of class `\link{init.icm}`.
#' @param control An `EpiModel` object of class `\link{control.icm}`.
#'
#' @return
#' This function returns no objects.
#'
#' @export
#' @keywords internal
#'
crosscheck.icm <- function(param, init, control) {

  ## Main class check
  if (!inherits(param, "param.icm")) {
    stop("param must an object of class param.icm", call. = FALSE)
  }
  if (!inherits(init, "init.icm")) {
    stop("init must an object of class init.icm", call. = FALSE)
  }
  if (!inherits(control, "control.icm")) {
    stop("control must an object of class control.icm", call. = FALSE)
  }

  if (control$skip.check == FALSE) {

    ## Check that rec.rate is supplied for SIR models
    if (control$type %in% c("SIR", "SIS")) {
      if (is.null(param$rec.rate)) {
        stop("Specify rec.rate in param.icm", call. = FALSE)
      }
      if (param$groups == 2 & is.null(param$rec.rate.g2)) {
        stop("Specify rec.rate.g2 in param.icm", call. = FALSE)
      }
    }


    ## Check that paramets and init are supplied for SIR models
    if (control$type == "SIR") {
      if (is.null(init$r.num)) {
        stop("Specify r.num in init.icm", call. = FALSE)
      }
      if (param$groups == 2 & is.null(init$r.num.g2)) {
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
      stop("Group 2 parameters specified in param.dcm, but missing group 2, ",
           "initial states in init.icm", call. = FALSE)
    }
    if (param$groups == 1 && init.groups == 2) {
      stop("Group 2 initial stats specified in init.dcm, but missing group 2 ",
           "parameters in param.icm", call. = FALSE)
    }

    ## Deprecated parameters
    bim <- grep(".FUN", names(formals(control.icm)), value = TRUE)
    um <- which(grepl(".FUN", names(control)) & !(names(control) %in% bim))
    if (length(um) == 0 && !is.null(control$type)) {
      if (!is.null(param$trans.rate)) {
        stop("The trans.rate parameter is deprecated. Use the inf.prob ",
             "parameter instead.", call. = FALSE)
      }
      if (!is.null(param$trans.rate.g2)) {
        stop("The trans.rate.g2 parameter is deprecated. Use the inf.prob.g2 ",
             "parameter instead.", call. = FALSE)
      }
    }

  }


  ## Assign modules based on group parameter
  if (!is.null(control$type)) {
    def <- grep(".FUN",names(control))
    args <- names(control)[def]

    if (param$groups == 1) {
      for (i in 1:length(args)) {
        if (is.null(control[[args[i]]])) {
          temp <- get(gsub(".FUN",".icm",args[i]))
          control[[args[i]]] <- temp
        }
      }
      message("Default modules set to appropriate two-group functions. See documentation ",
              "for details.")
    }
    else {
      for (i in 1:length(args)) {
        if (is.null(control[[args[i]]])) {
          temp <- get(gsub(".FUN",".icm.bip",args[i]))
          control[[args[i]]] <- temp
        }
      }
      message("Default modules set to appropriate two-group functions. See documentation ",
              "for details.")
    }
  }

  ## In-place assignment to update param and control
  on.exit(assign("param", param, pos = parent.frame()))
  on.exit(assign("control", control, pos = parent.frame()), add = TRUE)
}
