
#' @title Epidemic Parameters for Deterministic Compartmental Models
#'
#' @description Sets the epidemic parameters for deterministic compartmental
#'              models simulated with \code{dcm}.
#'
#' @param inf.prob Probability of infection per transmissible act between
#'        a susceptible and an infected person. In two-group models, this is the
#'        probability of infection for the group 1 members.
#' @param inter.eff Efficacy of an intervention which affects the per-act
#'        probability of infection. Efficacy is defined as 1 - the relative
#'        hazard of infection given exposure to the intervention, compared to no
#'        exposure.
#' @param inter.start Time step at which the intervention starts, between 1 and
#'        the number of time steps specified in the model. This will default to
#'        1 if the \code{inter.eff} is defined but this parameter is not.
#' @param act.rate Average number of transmissible acts per person per unit time.
#'        For two-group models, this is the number of acts per group 1 persons
#'        per unit time; a balance between the acts in groups 1 and 2 is necessary,
#'        and set using the \code{balance} parameter (see details).
#' @param rec.rate Average rate of recovery with immunity (in \code{SIR} models)
#'        or re-susceptibility (in \code{SIS} models). The recovery rate is the
#'        reciprocal of the disease duration. For two-group models, this is the
#'        recovery rate for group 1 persons only. This parameter is only used for
#'        \code{SIR} and \code{SIS} models.
#' @param a.rate Arrival or entry rate. For one-group models, the arrival rate is the
#'        rate of new arrivals per person per unit time. For two-group models, the
#'        arrival rate may be parameterized as a rate per group 1 person time (with
#'        group 1 persons representing females), and with the \code{a.rate.g2}
#'        rate set as described below.
#' @param ds.rate Departure or exit rate for susceptible. For two-group models, it
#'        is the rate for the group 1 susceptible only.
#' @param di.rate Departure or exit rate for infected. For two-group models, it is
#'        the rate for the group 1 infected only.
#' @param dr.rate Departure or exit rate for recovered. For two-group models, it is
#'        the rate for the group 1 recovered only. This parameter is only used for
#'        \code{SIR} models.
#' @param inf.prob.g2 Probability of infection per transmissible act
#'        between a susceptible group 2 person and an infected group 1 person.
#'        It is the probability of infection to group 2 members.
#' @param act.rate.g2 Average number of transmissible acts per group 2 person per
#'        unit time; a balance between the acts in groups 1 and 2 is necessary,
#'        and set using the \code{balance} parameter (see details).
#' @param rec.rate.g2 Average rate of recovery with immunity (in \code{SIR} models)
#'        or re-susceptibility (in \code{SIS} models) for group 2 persons. This
#'        parameter is only used for two-group \code{SIR} and \code{SIS} models.
#' @param a.rate.g2 Arrival or entry rate for group 2. This may either be specified
#'        numerically as the rate of new arrivals per group 2 persons per unit time,
#'        or as \code{NA} in which case the group 1 rate, \code{a.rate}, governs
#'        the group 2 rate. The latter is used when, for example, the first group
#'        is conceptualized as female, and the female population size determines
#'        the arrival rate. Such arrivals are evenly allocated between the two groups.
#' @param ds.rate.g2 Departure or exit rate for group 2 susceptible.
#' @param di.rate.g2 Departure or exit rate for group 2 infected.
#' @param dr.rate.g2 Departure or exit rate for group 2 recovered. This parameter is
#'        only used for \code{SIR} model types.
#' @param balance For two-group models, balance the \code{act.rate} to the rate
#'        set for group 1 (with \code{balance="g1"}) or group 2 (with
#'        \code{balance="g2"}). See details.
#' @param ... Additional arguments passed to model.
#'
#' @details
#' \code{param.dcm} sets the epidemic parameters for deterministic compartmental
#' models solved with the \code{\link{dcm}} function. The models may use the
#' base types, for which these parameters are used, or original model
#' specifications for which these parameters may be used (but not necessarily).
#' A detailed description of DCM parameterization for base models is found
#' in the \href{http://epimodel.org/tut.html}{Basic DCMs} tutorial.
#'
#' For base models, the model specification will be selected as a function
#' of the model parameters entered here and the control settings in
#' \code{\link{control.dcm}}. One-group and two-group models are available, where
#' the former assumes a homogeneous mixing in the population and the latter
#' assumes a purely heterogeneous mixing between two distinct partitions in the
#' population (e.g., men and women). Specifying any group two parameters (those
#' with a \code{.g2}) implies the simulation of a two-group model. All the
#' parameters for a desired model type must be specified, even if they are zero.
#'
#' @section Act Balancing:
#' In two-group models, a balance between the number of acts for group 1 members
#' and those for group 2 members must be maintained. With purely heterogeneous
#' mixing, the product of one group size and act rate must equal the product of
#' the other group size and act rate: \eqn{N_1 \alpha_1 = N_2 \alpha_2}, where
#' \eqn{N_i} is the group size and \eqn{\alpha_i} the group-specific act rates
#' at time \eqn{t}. The \code{balance} parameter here specifies which group's
#' act rate should control the others with respect to balancing. See the
#' \href{http://epimodel.org/tut.html}{Basic DCMs} tutorial for further details.
#'
#' @section Sensitivity Analyses:
#' \code{dcm} has been designed to easily run DCM sensitivity analyses, where a
#' series of models varying one or more of the model parameters is run. This is
#' possible by setting any parameter as a vector of length greater than one. See
#' both the example below and the \href{http://epimodel.org/tut.html}{Basic DCMs}
#' tutorial.
#'
#' @section New Model Types:
#' To build original model specifications outside of the base models, start
#' by consulting the \href{http://epimodel.org/tut.html}{New DCMs with EpiModel}
#' tutorial. Briefly, an original model may use either the existing model
#' parameters named here, an original set of parameters, or a combination of
#' both. The \code{...} argument allows the user to pass an arbitrary set of new
#' model parameters into \code{param.dcm}. Whereas there are strict checks for
#' base models that the model parameters are valid, parameter validity is the
#' user responsibility with these original models.
#'
#' @seealso Use \code{\link{init.dcm}} to specify the initial conditions and
#'          \code{\link{control.dcm}} to specify the control settings. Run the
#'          parameterized model with \code{\link{dcm}}.
#'
#' @keywords parameterization
#'
#' @export
#'
param.dcm <- function(inf.prob, inter.eff, inter.start, act.rate, rec.rate,
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
    message("EpiModel 1.7.0 onward renamed the birth rate parameter b.rate to a.rate. ",
            "See documentation for details.")
  }

  if (!is.null(p$inter.eff) && is.null(p$inter.start)) {
    p$inter.start <- 1
  }

  class(p) <- c("param.dcm", "list")
  return(p)
}


#' @title Initial Conditions for Deterministic Compartmental Models
#'
#' @description Sets the initial conditions for deterministic compartmental
#'              models simulated with \code{dcm}.
#'
#' @param s.num Number of initial susceptible. For two-group models, this is
#'        the number of initial group 1 susceptible.
#' @param i.num Number of initial infected. For two-group models, this is the
#'        number of initial group 1 infected.
#' @param r.num Number of initial recovered. For two-group models, this is the
#'        number of initial group 1 recovered. This parameter is only used for
#'        the \code{SIR} model type.
#' @param s.num.g2 Number of initial susceptible in group 2. This parameter is
#'        only used for two-group models.
#' @param i.num.g2 Number of initial infected in group 2. This parameter is only
#'        used for two-group models.
#' @param r.num.g2 Number of initial recovered in group 2. This parameter is
#'        only used for two-group \code{SIR} models.
#' @param ... Additional initial conditions passed to model.
#'
#' @details
#' The initial conditions for a model solved with \code{\link{dcm}} should be
#' input into the \code{init.dcm} function. This function handles initial
#' conditions for both base model types and original models. For an overview
#' of initial conditions for base DCM class models, consult the
#' \href{http://epimodel.org/tut.html}{Basic DCMs} tutorial.
#'
#' Original models may use the parameter names listed as arguments here, a new
#' set of names, or a combination of both. With new models, initial conditions
#' must be input in the same order that the solved derivatives from the model
#' are output. More details on this requirement are outlined in the
#' \href{http://epimodel.org/tut.html}{Solving New DCMs} tutorial.
#'
#' @seealso Use \code{\link{param.dcm}} to specify model parameters and
#'          \code{\link{control.dcm}} to specify the control settings. Run the
#'          parameterized model with \code{\link{dcm}}.
#'
#' @keywords parameterization
#'
#' @export
#'
init.dcm <- function(s.num, i.num, r.num, s.num.g2, i.num.g2, r.num.g2,
                     ...) {

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

  # Reorder arguments
  mc <- names(as.list(sys.call()[-1]))
  out.p <- list()
  for (i in seq_along(mc)) {
    out.p[[mc[i]]] <- p[[mc[i]]]
  }

  ## Output
  class(out.p) <- c("init.dcm", "list")
  return(out.p)
}


#' @title Control Settings for Deterministic Compartmental Models
#'
#' @description Sets the controls for deterministic compartmental models
#'              simulated with \code{\link{dcm}}.
#'
#' @param type Disease type to be modeled, with the choice of \code{"SI"} for
#'        Susceptible-Infected diseases, \code{"SIR"} for
#'        Susceptible-Infected-Recovered diseases, and \code{"SIS"} for
#'        Susceptible-Infected-Susceptible diseases.
#' @param nsteps Number of time steps to solve the model over or vector of times
#'        to solve the model over. If the number of time steps, then this must be
#'        a positive integer of length 1.
#' @param dt Time unit for model solutions, with the default of 1. Model
#'        solutions for fractional time steps may be obtained by setting this to a
#'        number between 0 and 1.
#' @param odemethod Ordinary differential equation (ODE) integration method, with
#'        the default of the "Runge-Kutta 4" method (see \code{\link{ode}} for
#'        other options).
#' @param dede If \code{TRUE}, use the delayed differential equation solver,
#'        which allows for time-lagged variables.
#' @param new.mod If not running an base model type, a function with a new
#'        model to be simulated (see details).
#' @param sens.param If \code{TRUE}, evaluate arguments in parameters with length
#'        greater than 1 as sensitivity analyses, with one model run per value of
#'        the parameter. If \code{FALSE}, one model will be run with parameters
#'        of arbitrary length.
#' @param print.mod If \code{TRUE}, print the model form to the console.
#' @param verbose If \code{TRUE}, print model progress to the console.
#' @param ... additional control settings passed to model.
#'
#' @details
#' \code{control.dcm} sets the required control settings for any deterministic
#' compartmental models solved with the \code{\link{dcm}} function. Controls are
#' required for both base model types and original models. For an overview of
#' control settings for base DCM class models, consult the
#' \href{http://epimodel.org/tut.html}{Basic DCMs} tutorial.
#' For all base models, the \code{type} argument is a necessary parameter
#' and it has no default.
#'
#' @section New Model Functions:
#' The form of the model function for base models may be displayed with the
#' \code{print.mod} argument set to \code{TRUE}. In this case, the model will not
#' be run. These model forms may be used as templates to write original model
#' functions.
#'
#' These new models may be input and solved with \code{\link{dcm}} using the
#' \code{new.mod} argument, which requires as input a model function. Details and
#' examples are found in the \href{http://epimodel.org/tut.html}{New DCMs}
#' tutorial.
#'
#' @seealso Use \code{\link{param.dcm}} to specify model parameters and
#'          \code{\link{init.dcm}} to specify the initial conditions. Run the
#'          parameterized model with \code{\link{dcm}}.
#'
#' @keywords parameterization
#'
#' @export
#'
control.dcm <- function(type, nsteps, dt = 1, odemethod = "rk4",
                        dede = FALSE, new.mod = NULL, sens.param = TRUE,
                        print.mod = FALSE, verbose = FALSE, ...) {



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
  if (!is.null(p$new.mod)) {
    p$new.mod.name <- as.list(match.call())$new.mod
  }

  ## Defaults and checks
  if (is.null(p$nsteps)) {
    stop("Specify nsteps")
  }


  # Check type for base models
  if (is.null(p$new.mod)) {
    if (is.null(p$type) || !(p$type %in% c("SI", "SIS", "SIR"))) {
      stop("Specify type as \"SI\", \"SIS\", or \"SIR\" ", call. = FALSE)
    }
  }

  ## Output
  class(p) <- c("control.dcm", "list")
  return(p)
}


#' @title Cross Checking of Inputs for Deterministic Compartmental Models
#'
#' @description This function checks that the three parameter lists from
#'              `\link{param.dcm}`, `\link{init.dcm}`, and
#'              `\link{control.dcm}` are consistent.
#'
#' @param param An `EpiModel` object of class `\link{param.dcm}`.
#' @param init An `EpiModel` object of class `\link{init.dcm}`.
#' @param control An `EpiModel` object of class `\link{control.dcm}`.
#'
#' @return
#' This function returns no objects.
#'
#' @export
#' @keywords internal
#'
crosscheck.dcm <- function(param, init, control) {

  # Main class check --------------------------------------------------------
  if (!inherits(param, "param.dcm")) {
    stop("param must an object of class param.dcm", call. = FALSE)
  }
  if (!inherits(init, "init.dcm")) {
    stop("init must an object of class init.dcm", call. = FALSE)
  }
  if (!inherits(control, "control.dcm")) {
    stop("control must an object of class control.dcm", call. = FALSE)
  }

  # Parameter checks for base models ----------------------------------
  if (is.null(control$new.mod)) {

    ## Defaults
    if (is.null(param$act.rate)) {
      param$act.rate <- 1
    }
    if (is.null(param$vital)) {
      if (!is.null(param$a.rate) |
          !is.null(param$ds.rate) |
          !is.null(param$di.rate) |
          !is.null(param$dr.rate)) {
        param$vital <- TRUE
      } else {
        param$vital <- FALSE
      }
    }

    if (any(grepl(".g2", names(param))) == TRUE) {
      param$groups <- 2
    } else {
      param$groups <- 1
    }

    if (param$groups == 2 && (is.null(param$balance) ||
                              !(param$balance %in% c("g1", "g2")))) {
      stop("Specify balance=\"g1\" or balance=\"g2\" with 2-group models",
           call. = FALSE)
    }

    ## Error checks
    # Specify inf.prob
    if (is.null(param$inf.prob)) {
      stop("Specify inf.prob in param.dcm", call. = FALSE)
    }

    # Check that rec.rate is supplied for SIR models
    if (control$type %in% c("SIR", "SIS") & is.null(param$rec.rate)) {
      stop("Specify rec.rate in param.dcm", call. = FALSE)
      if (param$groups == 2 & is.null(param$rec.rate.g2)) {
        stop("Specify rec.rate.g2 in param.dcm", call. = FALSE)
      }
    }

    # Check that r.num is supplied for SIR models
    if (control$type == "SIR" & is.null(init$r.num)) {
      stop("Specify r.num in init.dcm", call. = FALSE)
      if (param$groups == 2 & is.null(init$r.num.g2)) {
        stop("Specify r.num.g2 in init.dcm", call. = FALSE)
      }
    }

    # Check that groups implied by init and params are consistent
    if (any(grepl(".g2", names(init))) == TRUE) {
      init.groups <- 2
    } else {
      init.groups <- 1
    }
    if (param$groups == 2 && init.groups == 1) {
      stop("Group 2 parameters specified in param.dcm,
           \rbut missing group 2 initial states in init.dcm",
           call. = FALSE)
    }
    if (param$groups == 1 && init.groups == 2) {
      stop("Group 2 initial stats specified in init.dcm,
           but missing group 2 parameters in param.dcm",
           call. = FALSE)
    }

    # Over-specified initial conditions
    if (control$type != "SIR" & any(c("r.num", "r.num.g2") %in% names(init))) {
      stop("Specified initial number recovered for non-SIR model",
           call. = FALSE)
    }

    # Deprecated parameters
    if (!is.null(param$trans.rate)) {
      stop("The trans.rate parameter is deprecated. Use the inf.prob parameter instead.",
           call. = FALSE)
    }
    if (!is.null(param$trans.rate.g2)) {
      stop("The trans.rate.g2 parameter is deprecated. Use the inf.prob.g2 parameter instead.",
           call. = FALSE)
    }
  }

  on.exit(assign("param", param, pos = parent.frame()))
}
