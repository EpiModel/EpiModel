
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
#' in the
#' \href{http://statnet.github.io/tut/BasicDCMs.html}{Basic DCMs} tutorial.
#'
#' For base models, the model specification will be selected as a function
#' of the model parameters entered here and the control settings in
#' \code{\link{control.dcm}}. One-group and two-group models are available, where
#' the former assumes a homogenous mixing in the population and the latter
#' assumes a purely heterogenous mixing between two distinct partitions in the
#' population (e.g., men and women). Specifying any group two parameters (those
#' with a \code{.g2}) implies the simulation of a two-group model. All the
#' parameters for a desired model type must be specified, even if they are zero.
#'
#' @section Act Balancing:
#' In two-group models, a balance between the number of acts for group 1 members
#' and those for group 2 members must be maintained. With purely heterogenous
#' mixing, the product of one group size and act rate must equal the product of
#' the other group size and act rate: \eqn{N_1 \alpha_1 = N_2 \alpha_2}, where
#' \eqn{N_i} is the group size and \eqn{\alpha_i} the group-specific act rates
#' at time \eqn{t}. The \code{balance} parameter here specifies which group's
#' act rate should control the others with respect to balancing. See the
#' \href{http://statnet.github.io/tut/BasicDCMs.html}{Basic DCMs} tutorial
#' for further details.
#'
#' @section Sensitivity Analyses:
#' \code{dcm} has been designed to easily run DCM sensitivity analyses, where a
#' series of models varying one or more of the model parameters is run. This is
#' possible by setting any parameter as a vector of length greater than one. See
#' both the example below and the
#' \href{http://statnet.github.io/tut/BasicDCMs.html}{Basic DCMs} tutorial.
#'
#' @section New Model Types:
#' To build original model specifications outside of the base models, start
#' by consulting the \href{http://statnet.github.io/tut/NewDCMs.html}{Solving
#' New DCMs with EpiModel} tutorial. Briefly, an original model may use either
#' the existing model parameters named here, an original set of parameters, or
#' a combination of both. The \code{...} argument allows the user to pass an
#' arbitrary set of new model parameters into \code{param.dcm}. Whereas there are
#' strict checks for base models that the model parameters are valid,
#' parameter validity is the user's responsibility with these original models.
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
    a.rate <- dot.args$b.rate
    message("EpiModel 1.7.0 onward renamed the birth rate parameter b.rate to a.rate. See documentation for details.")
  }

  if ("b.flow" %in% names.dot.args) {
    a.flow <- dot.args$b.flow
    message("EpiModel 1.7.0 onward renamed the birth summary statistic b.flow to a.flow. See documentation for details.")
  }

  if (!is.null(p$inter.eff) && is.null(p$inter.start)) {
    p$inter.start <- 1
  }

  class(p) <- c("param.dcm", "list")
  return(p)
}


#' @title Epidemic Parameters for Stochastic Individual Contact Models
#'
#' @description Sets the epidemic parameters for stochastic individual contact
#'              models simulated with \code{icm}.
#'
#' @inheritParams param.dcm
#'
#' @details
#' \code{param.icm} sets the epidemic parameters for the stochastic individual
#' contact models simulated with the \code{\link{icm}} function. Models
#' may use the base types, for which these parameters are used, or new process
#' modules which may use these parameters (but not necessarily). A detailed
#' description of ICM parameterization for base models is found in the
#' \href{http://statnet.github.io/tut/BasicICMs.html}{Basic ICMs} tutorial.
#'
#' For base models, the model specification will be chosen as a result of
#' the model parameters entered here and the control settings in
#' \code{\link{control.icm}}. One-group and two-group models are available, where
#' the former assumes a homogenous mixing in the population and the latter
#' assumes a purely heterogenous mixing between two distinct partitions in the
#' population (e.g., men and women). Specifying any group two parameters (those
#' with a \code{.g2}) implies the simulation of a two-group model. All the
#' parameters for a desired model type must be specified, even if they are zero.
#'
#' @section Act Balancing:
#' In two-group models, a balance between the number of acts for group 1 members
#' and those for group 2 members must be maintained. With purely heterogenous
#' mixing, the product of one group size and act rate must equal the product of
#' the other group size and act rate: \eqn{N_1 \alpha_1 = N_2 \alpha_2}, where
#' \eqn{N_i} is the group size and \eqn{\alpha_i} the group-specific act rates
#' at time \eqn{t}. The \code{balance} parameter here specifies which group's
#' act rate should control the others with respect to balancing. See the
#' \href{http://statnet.github.io/tut/BasicDCMs.html}{Basic DCMs} tutorial.
#'
#' @section New Modules:
#' To build original models outside of the base models, new process modules
#' may be constructed to replace the existing modules or to supplement the existing
#' set. These are passed into the control settings in \code{\link{control.icm}}.
#' New modules may use either the existing model parameters named here, an
#' original set of parameters, or a combination of both. The \code{...} allows
#' the user to pass an arbitrary set of original model parameters into
#' \code{param.icm}. Whereas there are strict checks with default modules for
#' parameter validity, these checks are the user's responsibility with new modules.
#'
#' @seealso Use \code{\link{init.icm}} to specify the initial conditions and
#'          \code{\link{control.icm}} to specify the control settings. Run the
#'          parameterized model with \code{\link{icm}}.
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
    a.rate <- dot.args$b.rate
    message("EpiModel 1.7.0 onward renamed the birth rate parameter b.rate to a.rate. See documentation for details.")
  }

  if ("b.flow" %in% names.dot.args) {
    a.flow <- dot.args$b.flow
    message("EpiModel 1.7.0 onward renamed the birth summary statistic b.flow to a.flow. See documentation for details.")
  }

  if ("b.rand" %in% names.dot.args) {
    a.rand <- dot.args$b.rand
    message("EpiModel 1.7.0 onward renamed the stochastic birth model flag b.rand to a.rand. See documentation for details.")
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


#' @title Epidemic Parameters for Stochastic Network Models
#'
#' @description Sets the epidemic parameters for stochastic network models
#'              simulated with \code{\link{netsim}}.
#'
#' @param inf.prob Probability of infection per transmissible act between
#'        a susceptible and an infected person. In bipartite models, this is the
#'        probability of infection to the mode 1 nodes. This may also be a vector
#'        of probabilities, with each element corresponding to the probability
#'        in that time step of infection (see Time-Varying Parameters below).
#' @param inter.eff Efficacy of an intervention which affects the per-act
#'        probability of infection. Efficacy is defined as 1 - the relative
#'        hazard of infection given exposure to the intervention, compared to no
#'        exposure.
#' @param inter.start Time step at which the intervention starts, between 1 and
#'        the number of time steps specified in the model. This will default to
#'        1 if the \code{inter.eff} is defined but this parameter is not.
#' @param act.rate Average number of transmissible acts \emph{per partnership}
#'        per unit time (see \code{act.rate} Parameter below). This may also be a vector
#'        of rates, with each element corresponding to the rate in in that time
#'        step of infection (see Time-Varying Parameters below).
#' @param rec.rate Average rate of recovery with immunity (in \code{SIR} models)
#'        or re-susceptibility (in \code{SIS} models). The recovery rate is the
#'        reciprocal of the disease duration. For bipartite models, this is the
#'        recovery rate for mode 1 persons only. This parameter is only used for
#'        \code{SIR} and \code{SIS} models. This may also be a vector
#'        of rates, with each element corresponding to the rate in that time step
#'        of infection (see Time-Varying Parameters below).
#' @param a.rate Arrival or entry rate. For one-mode models, the arrival rate is the
#'        rate of new arrivals per person per unit time. For bipartite models, the
#'        arrival rate may be parameterized as a rate per mode 1 person time (with
#'        mode 1 persons representing females), and with the \code{a.rate.g2}
#'        rate set as described below.
#' @param ds.rate Departure or exit rate for susceptible. For bipartite models, it
#'        is the rate for the mode 1 susceptible only.
#' @param di.rate Departure or exit rate for infected. For bipartite models, it is
#'        the rate for the mode 1 infected only.
#' @param dr.rate Departure or exit rate for recovered. For bipartite models, it is
#'        the rate for the mode 1 recovered only. This parameter is only used for
#'        \code{SIR} models.
#' @param inf.prob.m2 Probability of transmission given a transmissible act
#'        between a susceptible mode 2 person and an infected mode 1 person.
#'        It is the probability of transmission to mode 2 members.
#' @param rec.rate.m2 Average rate of recovery with immunity (in \code{SIR} models)
#'        or re-susceptibility (in \code{SIS} models) for mode 2 persons. This
#'        parameter is only used for bipartite \code{SIR} and \code{SIS} models.
#' @param a.rate.m2 Arrival or entry rate for mode 2. This may either be specified
#'        numerically as the rate of new arrivals per mode 2 persons per unit time,
#'        or as \code{NA} in which case the mode 1 rate, \code{a.rate}, governs
#'        the mode 2 rate. The latter is used when, for example, the first mode
#'        is conceptualized as female, and the female population size determines
#'        the arrival rate. Such arrivalss are evenly allocated between the two modes.
#' @param ds.rate.m2 Departure or exit rate for mode 2 susceptible.
#' @param di.rate.m2 Departure or exit rate for mode 2 infected.
#' @param dr.rate.m2 Departure or exit rate for mode 2 recovered. This parameter is
#'        only used for \code{SIR} model types.
#' @param ... Additional arguments passed to model.
#'
#' @details
#' \code{param.net} sets the epidemic parameters for the stochastic network
#' models simulated with the \code{\link{netsim}} function. Models
#' may use the base types, for which these parameters are used, or new process
#' modules which may use these parameters (but not necessarily). A detailed
#' description of network model parameterization for base models is found in
#' the \href{http://statnet.github.io/tut/BasicNet.html}{Basic Network Models}
#' tutorial.
#'
#' For base models, the model specification will be chosen as a result of
#' the model parameters entered here and the control settings in
#' \code{\link{control.net}}. One-mode and two-mode models are available, where
#' the the latter assumes a heterogenous mixing between two distinct partitions
#' in the population (e.g., men and women). Specifying any bipartite parameters
#' (those with a \code{.m2}) implies the simulation of a bipartite model. All the
#' parameters for a desired model type must be specified, even if they are zero.
#'
#' @section The \code{act.rate} Parameter:
#' A key difference between these network models and DCM/ICM classes is the
#' treatment of transmission events. With DCM and ICM, contacts or partnerships are
#' mathematically instantaneous events: they have no duration in time, and thus
#' no changes may occur within them over time. In contrast, network models allow
#' for partnership durations defined by the dynamic network model, summarized in
#' the model dissolution coefficients calculated in \code{\link{dissolution_coefs}}.
#' Therefore, the \code{act.rate} parameter has a different interpretation here,
#' where it is the number of transmissible acts \emph{per partnership} per unit
#' time.
#'
#' @section Time-Varying Parameters:
#' The \code{inf.prob}, \code{act.rate}, \code{rec.rate} arguments (and their
#' \code{.m2} companions) may be specified as time-varying parameters by passing
#' in a vector of probabilities or rates, respectively. The value in each
#' position on the vector then corresponds to the probability or rate at that
#' discrete time step for the infected partner. For example, an \code{inf.prob}
#' of \code{c(0.5, 0.5, 0.1)} would simulate a 0.5 transmission probability for
#' the first two time steps of a person's infection, followed by a 0.1 for the
#' third time step. If the infected person has not recovered or exited the
#' population by the fourth time step, the third element in the vector will carry
#' forward until one of those events occurs or the simulation ends. For further examples,
#' see the NME tutorial,
#' \href{https://statnet.github.io/nme/d3-s6.html}{Time-Varying Biology & Behavior}.
#'
#' @section New Modules:
#' To build original models outside of the base models, new process modules
#' may be constructed to replace the existing modules or to supplement the existing
#' set. These are passed into the control settings in \code{\link{control.net}}.
#' New modules may use either the existing model parameters named here, an
#' original set of parameters, or a combination of both. The \code{...} allows
#' the user to pass an arbitrary set of original model parameters into
#' \code{param.net}. Whereas there are strict checks with default modules for
#' parameter validity, these checks are the user's responsibility with new modules.
#'
#' @seealso Use \code{\link{init.net}} to specify the initial conditions and
#'          \code{\link{control.net}} to specify the control settings. Run the
#'          parameterized model with \code{\link{netsim}}.
#'
#' @keywords parameterization
#'
#' @export
#'
param.net <- function(inf.prob, inter.eff, inter.start, act.rate, rec.rate,
                      a.rate, ds.rate, di.rate, dr.rate, inf.prob.m2,
                      rec.rate.m2, a.rate.m2, ds.rate.m2, di.rate.m2,
                      dr.rate.m2, ...) {

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

  ##Updated parameter names
  if ("b.rate" %in% names.dot.args) {
    a.rate <- dot.args$b.rate
    message("EpiModel 1.7.0 onward renamed the birth rate parameter b.rate to a.rate. See documentation for details.")
  }

  ## Defaults and checks
  if (missing(act.rate)) {
    p$act.rate <- 1
  }
  p$vital <- ifelse(!missing(a.rate) | !missing(ds.rate) |
                    !missing(di.rate) | !missing(dr.rate), TRUE, FALSE)
  if ("act.rate.m2" %in% names.dot.args) {
    warning("act.rate.m2 parameter was entered. If using built-in models, only act.rate parameter will apply.",
            call. = FALSE)
  }


  if (!is.null(p$inter.eff) && is.null(p$inter.start)) {
    p$inter.start <- 1
  }

  ## Output
  class(p) <- c("param.net", "list")
  return(p)
}

