
#' @title Epidemic Parameters for Stochastic Network Models
#'
#' @description Sets the epidemic parameters for stochastic network models
#'              simulated with \code{\link{netsim}}.
#'
#' @param inf.prob Probability of infection per transmissible act between
#'        a susceptible and an infected person. In two-group models, this is the
#'        probability of infection to the group 1 nodes. This may also be a vector
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
#'        reciprocal of the disease duration. For two-group models, this is the
#'        recovery rate for mode 1 persons only. This parameter is only used for
#'        \code{SIR} and \code{SIS} models. This may also be a vector
#'        of rates, with each element corresponding to the rate in that time step
#'        of infection (see Time-Varying Parameters below).
#' @param a.rate Arrival or entry rate. For one-mode models, the arrival rate is the
#'        rate of new arrivals per person per unit time. For two-group models, the
#'        arrival rate may be parameterized as a rate per mode 1 person time (with
#'        group 1 persons representing females), and with the \code{a.rate.g2}
#'        rate set as described below.
#' @param ds.rate Departure or exit rate for susceptible. For two-group models, it
#'        is the rate for the group 1 susceptible only.
#' @param di.rate Departure or exit rate for infected. For two-group models, it is
#'        the rate for the group 1 infected only.
#' @param dr.rate Departure or exit rate for recovered. For two-group models, it is
#'        the rate for the group 1 recovered only. This parameter is only used for
#'        \code{SIR} models.
#' @param inf.prob.g2 Probability of transmission given a transmissible act
#'        between a susceptible group 2 person and an infected group 1 person.
#'        It is the probability of transmission to group 2 members.
#' @param rec.rate.g2 Average rate of recovery with immunity (in \code{SIR} models)
#'        or re-susceptibility (in \code{SIS} models) for group 2 persons. This
#'        parameter is only used for two-group \code{SIR} and \code{SIS} models.
#' @param a.rate.g2 Arrival or entry rate for group 2. This may either be specified
#'        numerically as the rate of new arrivals per group 2 persons per unit time,
#'        or as \code{NA} in which case the mode 1 rate, \code{a.rate}, governs
#'        the group 2 rate. The latter is used when, for example, the first group
#'        is conceptualized as female, and the female population size determines
#'        the arrival rate. Such arrivalss are evenly allocated between the two groups.
#' @param ds.rate.g2 Departure or exit rate for group 2 susceptible.
#' @param di.rate.g2 Departure or exit rate for group 2 infected.
#' @param dr.rate.g2 Departure or exit rate for group 2 recovered. This parameter is
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
#' \code{\link{control.net}}. One-group and two-group models are available, where
#' the latter assumes a heterogenous mixing between two distinct partitions
#' in the population (e.g., men and women). Specifying any two-group parameters
#' (those with a \code{.g2}) implies the simulation of a two-group model. All the
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
#' \code{.g2} companions) may be specified as time-varying parameters by passing
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
                      a.rate, ds.rate, di.rate, dr.rate, inf.prob.g2,
                      rec.rate.g2, a.rate.g2, ds.rate.g2, di.rate.g2,
                      dr.rate.g2, ...) {

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
  if (missing(act.rate)) {
    p$act.rate <- 1
  }
  p$vital <- ifelse(!missing(a.rate) | !missing(ds.rate) |
                      !missing(di.rate) | !missing(dr.rate), TRUE, FALSE)
  if ("act.rate.g2" %in% names.dot.args) {
    warning("act.rate.g2 parameter was entered. ",
            "If using built-in models, only act.rate parameter will apply.",
            call. = FALSE)
  }


  if (!is.null(p$inter.eff) && is.null(p$inter.start)) {
    p$inter.start <- 1
  }

  ## Output
  class(p) <- c("param.net", "list")
  return(p)
}


#' @title Initial Conditions for Stochastic Network Models
#'
#' @description Sets the initial conditions for stochastic network models
#'              simulated with \code{netsim}.
#'
#' @param i.num Number of initial infected. For two-group models, this is the
#'        number of initial group 1 infected.
#' @param r.num Number of initial recovered. For two-group models, this is the
#'        number of initial group 1 recovered. This parameter is only used for
#'        the \code{SIR} model type.
#' @param i.num.g2 Number of initial infected in group 2. This parameter is only
#'        used for two-group models.
#' @param r.num.g2 Number of initial recovered in group 2. This parameter is
#'        only used for two-group \code{SIR} models.
#' @param status.vector A vector of length equal to the size of the input network,
#'        containing the status of each node. Setting status here overrides any
#'        inputs passed in the \code{.num} arguments.
#' @param infTime.vector A vector of length equal to the size of the input network,
#'        containing the (historical) time of infection for each of those nodes
#'        with a current status of \code{"i"}. Can only be used if \code{status.vector}
#'        is used, and must contain \code{NA} values for any nodes whose status
#'        is not \code{"i"}.
#' @param ... Additional initial conditions passed to model.
#'
#' @details
#' The initial conditions for a model solved with \code{\link{netsim}} should be
#' input into the \code{init.net} function. This function handles initial
#' conditions for both base models and new modules. For an overview of
#' specifying initial conditions across a variety of base network models,
#' consult the \href{http://statnet.github.io/tut/BasicNet.html}{Basic Network
#' Models} tutorial.
#'
#' @seealso Use \code{\link{param.net}} to specify model parameters and
#'          \code{\link{control.net}} to specify the control settings. Run the
#'          parameterized model with \code{\link{netsim}}.
#'
#' @keywords parameterization
#'
#' @export
#'
#' @examples
#' # Example of using status.vector and infTime.vector together
#' n <- 100
#' status <- sample(c("s", "i"), size = n, replace = TRUE, prob = c(0.8, 0.2))
#' infTime <- rep(NA, n)
#' infTime[which(status == "i")] <- -rgeom(sum(status == "i"), prob = 0.01) + 2
#'
#' init.net(status.vector = status, infTime.vector = infTime)
#'
init.net <- function(i.num, r.num, i.num.g2, r.num.g2,
                     status.vector, infTime.vector, ...) {

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

  ## Defaults and checks
  if (!is.null(p$i.num) & !is.null(p$status.vector)) {
    stop("Use i.num OR status.vector to set initial infected")
  }
  if (!is.null(p$infTime.vector) & is.null(p$status.vector)) {
    stop("infTime.vector may only be used if status.vector is used")
  }
  if (!is.null(p$infTime.vector) & length(p$infTime.vector) != length(p$status.vector)) {
    stop("Length of infTime.vector must match length of status.vector")
  }


  ## Output
  class(p) <- c("init.net", "list")
  return(p)
}


#' @title Control Settings for Stochastic Network Models
#'
#' @description Sets the controls for stochastic network models simulated with
#'              \code{\link{netsim}}.
#'
#' @param type Disease type to be modeled, with the choice of \code{"SI"} for
#'        Susceptible-Infected diseases, \code{"SIR"} for
#'        Susceptible-Infected-Recovered diseases, and \code{"SIS"} for
#'        Susceptible-Infected-Susceptible diseases.
#' @param nsteps Number of time steps to simulate the model over. This must be a
#'        positive integer that is equal to the final step of a simulation. If
#'        simulation is restarted with \code{start} argument, this number must
#'        be at least one greater than that argument's value.
#' @param nsims The total number of disease simulations.
#' @param ncores Number of processor cores to run multiple simulations
#'        on, using the \code{foreach} and \code{doParallel} implementations.
#' @param start For dependent simulations, time point to start up simulation.
#'        For restarted simulations, this must be one greater than the final time
#'        step in the prior simulation and must be less than the value in
#'        \code{nsteps}.
#' @param resimulate.network If \code{TRUE}, resimulate the network at each time step. This
#'        occurs by default with two varieties of dependent models: if there are
#'        any vital dynamic parameters in the model (or if non-standard arrival or
#'        departures modules are passed into \code{control.net}), or if the network model
#'        formation formula includes the "status" attribute.
#' @param rec.rand If \code{TRUE}, use a stochastic recovery model, with the
#'        number of recovered at each time step a function of random draws from
#'        a binomial distribution with the probability equal to \code{rec.rate}.
#'        If \code{FALSE}, then a deterministic rounded count of the expectation
#'        implied by that rate.
#' @param a.rand If \code{TRUE}, use a stochastic arrival model, with the
#'        number of arrivals at each time step a function of random draws from a
#'        binomial distribution with the probability equal to the governing arrival
#'        rates. If \code{FALSE}, then a deterministic rounded count of the
#'        expectation implied by those rates.
#' @param d.rand If \code{TRUE}, use a stochastic departure model, with the number of
#'        departures at each time step a function of random draws from a binomial
#'        distribution with the probability equal to the governing departure rates.
#'        If \code{FALSE}, then a deterministic rounded count of the expectation
#'        implied by those rates.
#' @param tergmLite Logical indicating usage of either \code{tergm} (\code{tergmLite = FALSE}),
#'        or \code{tergmLite} (\code{tergmLite = TRUE}). Default of \code{FALSE}.
#' @param attr.rules A list containing the  rules for setting the attributes of
#'        incoming nodes, with one list element per attribute to be set (see
#'        details below).
#' @param epi.by A character vector of length 1 containing a nodal attribute for
#'        which subgroup epidemic prevalences should be calculated. This nodal
#'        attribute must be contained in the network model formation formula,
#'        otherwise it is ignored.
#' @param initialize.FUN Module to initialize the model at time 1, with the
#'        default function of \code{\link{initialize.net}}.
#' @param departures.FUN Module to simulate departure or exit, with the default function
#'        of \code{\link{departures.net}}.
#' @param arrivals.FUN Module to simulate arrivals or entries, with the default
#'        function of \code{\link{arrivals.net}}.
#' @param recovery.FUN Module to simulate disease recovery, with the default
#'        function of \code{\link{recovery.net}}.
#' @param resim_nets.FUN Module to resimulate the network at each time step,
#'        with the default function of \code{\link{resim_nets}}.
#' @param infection.FUN Module to simulate disease infection, with the default
#'        function of \code{\link{infection.net}}.
#' @param nwupdate.FUN Module to handle updating of network structure and nodal
#'        attributes due to exogenous epidemic model processes, with the default
#'        function of \code{\link{nwupdate.net}}.
#' @param prevalence.FUN Module to calculate disease prevalence at each time step,
#'        with the default function of \code{\link{prevalence.net}}.
#' @param verbose.FUN Module to print simulation progress to screen, with the
#'        default function of \code{\link{verbose.net}}.
#' @param module.order A character vector of module names that lists modules the
#'        order in which they should be evaluated within each time step. If
#'        \code{NULL}, the modules will be evaluated as follows: first any
#'        new modules supplied through \code{...} in the order in which they are
#'        listed, then the built-in modules in their order of the function listing.
#'        The \code{initialize.FUN} will always be run first and the
#'        \code{verbose.FUN} always last.
#' @param set.control.ergm Control arguments passed to simulate.ergm. See the
#'        help file for \code{\link{netdx}} for details and examples on specifying
#'        this parameter.
#' @param set.control.stergm Control arguments passed to simulate.stergm. See the
#'        help file for \code{\link{netdx}} for details and examples on specifying
#'        this parameter.
#' @param save.nwstats If \code{TRUE}, save network statistics in a data frame.
#'        The statistics to be saved are specified in the \code{nwstats.formula}
#'        argument.
#' @param nwstats.formula A right-hand sided ERGM formula that includes network
#'        statistics of interest, with the default to the formation formula terms.
#' @param save.transmat If \code{TRUE}, save a transmission matrix for each
#'        simulation. This object contains one row for each transmission event
#'        (see \code{\link{discord_edgelist}}).
#' @param save.network If \code{TRUE}, save a \code{networkDynamic} object
#'        containing full edge history for each simulation.
#' @param save.other A vector of elements on the \code{dat} master data list
#'        to save out after each simulation. One example for base models is
#'        the attribute list, "attr", at the final time step.
#' @param verbose If \code{TRUE}, print model progress to the console.
#' @param verbose.int Time step interval for printing progress to console, where
#'        0 prints completion status of entire simulation and positive integer
#'        \code{x} prints progress after each \code{x} time steps. The default
#'        is to print progress after each time step.
#' @param skip.check If \code{TRUE}, skips the default error checking for the
#'        structure and consistency of the parameter values, initial conditions,
#'        and control settings before running base epidemic models. Setting
#'        this to \code{FALSE} is recommended when running models with new modules
#'        specified.
#' @param raw_output If \code{TRUE}, \code{netsim} will output a list of nestsim
#'        Data (one per simulation) instead of a formatted \code{netsim} object.
#' @param ... Additional control settings passed to model.
#'
#' @details
#' \code{control.net} sets the required control settings for any network model
#' solved with the \code{\link{netsim}} function. Controls are required for both
#' base model types and when passing original process modules. For an overview
#' of control settings for base models, consult the
#' \href{http://statnet.github.io/tut/BasicNet.html}{Basic Network Models} tutorial.
#' For all base models, the \code{type} argument is a necessary parameter
#' and it has no default.
#'
#' @section The attr.rules Argument:
#' The \code{attr.rules} parameter is used to specify the rules for how nodal
#' attribute values for incoming nodes should be set. These rules are only
#' necessary for models in which there are incoming nodes (i.e., arrivals) and
#' there is also a nodal attribute in the network model formation formula set in
#' \code{\link{netest}}. There are three rules available for each attribute
#' value:
#' \itemize{
#'  \item \strong{"current":} new nodes will be assigned this attribute in
#'        proportion to the distribution of that attribute among existing nodes
#'        at that current time step.
#'  \item \strong{"t1":} new nodes will be assigned this attribute in proportion
#'        to the distribution of that attribute among nodes at time 1 (that is,
#'        the proportions set in the original network for \code{\link{netest}}).
#'  \item \strong{<Value>:} all new nodes will be assigned this specific value,
#'        with no variation.
#' }
#' For example, the rules list
#' \code{attr.rules = list(race = "t1", sex = "current", status = "s")}
#' specifies how the race, sex, and status attributes should be set for incoming
#' nodes. By default, the rule is "current" for all attributes except status,
#' in which case it is "s" (that is, all incoming nodes are susceptible).
#'
#' @section New Modules:
#' Base network models use a set of module functions that specify how the
#' individual nodes in the network are subjected to infection, recovery,
#' demographics, and other processes. Core modules are those listed in the
#' \code{.FUN} arguments. For each module, there is a default function used in
#' the simulation. The default infection module, for example, is contained in
#' the \code{\link{infection.net}} function.
#'
#' For original models, one may substitute replacement module functions for any of
#' the default functions. New modules may be added to the workflow at each time
#' step by passing a module function via the \code{...} argument. Consult the
#' \href{http://statnet.github.io/tut/NewNet.html}{New Network Models} tutorial.
#' One may remove existing modules, such as \code{arrivals.FUN}, from the workflow
#' by setting the parameter value for that argument to \code{NULL}.
#'
#' @seealso Use \code{\link{param.net}} to specify model parameters and
#'          \code{\link{init.net}} to specify the initial conditions. Run the
#'          parameterized model with \code{\link{netsim}}.
#'
#' @keywords parameterization
#'
#' @export
#'
control.net <- function(type,
                        nsteps,
                        start = 1,
                        nsims = 1,
                        ncores = 1,
                        resimulate.network,
                        rec.rand = TRUE,
                        a.rand = TRUE,
                        d.rand = TRUE,
                        tergmLite = FALSE,
                        attr.rules,
                        epi.by,
                        initialize.FUN = initialize.net,
                        resim_nets.FUN = resim_nets,
                        infection.FUN = NULL,
                        recovery.FUN = NULL,
                        departures.FUN = NULL,
                        arrivals.FUN = NULL,
                        nwupdate.FUN = nwupdate.net,
                        prevalence.FUN = NULL,
                        verbose.FUN = verbose.net,
                        module.order = NULL,
                        set.control.ergm,
                        set.control.stergm,
                        save.nwstats = TRUE,
                        nwstats.formula = "formation",
                        save.transmat = TRUE,
                        save.network = TRUE,
                        save.other,
                        verbose = TRUE,
                        verbose.int = 1,
                        skip.check = FALSE,
                        raw_output = FALSE,
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
  bi.mods <- grep(".FUN", names(formal.args), value = TRUE)
  bi.nms <- bi.mods
  p$bi.mods <- bi.mods
  p$user.mods <- grep(".FUN", names(dot.args), value = TRUE)


  if (missing(resimulate.network)) {
    arg.list <- as.list(match.call())
    if ((!is.null(arg.list$departures.FUN) && arg.list$departures.FUN != "departures.net") |
        (!is.null(arg.list$arrivals.FUN) && arg.list$arrivals.FUN != "arrivals.net")) {
      p$resimulate.network <- TRUE
    }
  }

  # Using tergmLite --> resimulate.network = TRUE
  if (tergmLite == TRUE) {
    p$resimulate.network <- TRUE
  }

  # Temporary until we develop a nwstats fix for tergmLite
  if (tergmLite == TRUE) {
    p$save.nwstats <- FALSE
  }

  ## Defaults and checks

  #Check whether any base modules have been redefined by user (note: must come after above)
  bi.nms <- bi.nms[-which(bi.nms %in% c("initialize.FUN", "resim_nets.FUN", "verbose.FUN", "nwupdate.FUN"))]
  if (length(bi.nms) > 0){
    flag1 <- logical()
    for (args in 1:length(bi.nms)) {
      if (!(is.null(p[[bi.nms[args]]])) ) {
        temp1 <- get(gsub(".FUN",".net",bi.nms[args]))
        temp2 <- p[[bi.nms[args]]]
        flag1[args] <- identical(temp1,temp2)
      }
    }

    if (!is.null(p$type) && sum(flag1, na.rm = TRUE) != length(flag1)) {
      stop("Control parameter 'type' must be null if any user defined base modules are present")
    }
  }

  if (!is.null(p$type) && length(p$user.mods) > 0) {
    stop("Control parameter 'type' must be null if any user specified modules are present")
  }

  if (is.null(p$nsteps)) {
    stop("Specify nsteps")
  }

  if (missing(attr.rules)) {
    p$attr.rules <- list()
  }

  if (!is.null(p$epi.by)) {
    if (length(p$epi.by) > 1) {
      stop("Length of epi.by currently limited to 1")
    } else {
      p$epi.by <- epi.by
    }
  }

  if (is.null(p$set.control.stergm)) {
    p$set.control.stergm <- control.simulate.network(MCMC.burnin.min = 1000)
  }
  if (is.null(p$set.control.ergm)) {
    p$set.control.ergm <- control.simulate.ergm(MCMC.burnin = 2e5)
  }

  if (is.null(p$type)) {
    names <- unlist(lapply(sys.call()[-1], as.character))
    pos <- which(names(names) %in% grep(".FUN", names(names), value = TRUE))
    p$f.names <- as.vector(names)[pos]
    p$f.args  <- grep(".FUN", names(names), value = TRUE)
  }

  if (p$type != "SIR" && !is.null(p$type)) {
    p$f.names <- c("arrivals.FUN", "departures.FUN", "infection.FUN", "prevalence.FUN")
    p$f.args  <- c("arrivals.net", "departures.net", "infection.net", "prevalence.net")
  }

  if (p$type %in% c("SIR", "SIS") && !is.null(p$type)) {
    p$f.names <- c("arrivals.FUN", "departures.FUN", "infection.FUN",
                   "recovery.FUN", "prevalence.FUN")
    p$f.args  <- c("arrivals.net", "departures.net", "infection.net",
                   "recovery.net", "prevalence.net")
  }

  ## Output
  class(p) <- c("control.net", "list")
  return(p)
}


#' @title Cross Checking of Inputs for Stochastic Network Models
#'
#' @description This function checks that the estimation object from
#'              \code{\link{netest}} and the three parameter lists from
#'              \code{\link{param.net}}, \code{\link{init.net}}, and
#'              \code{\link{control.net}} are consistent.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netest}}.
#' @param param An \code{EpiModel} object of class \code{\link{param.net}}.
#' @param init An \code{EpiModel} object of class \code{\link{init.net}}.
#' @param control An \code{EpiModel} object of class \code{\link{control.net}}.
#'
#' @return
#' This function returns no objects.
#'
#' @export
#' @keywords internal
#'
crosscheck.net <- function(x, param, init, control) {

  if (!is.null(control$type) && length(control$user.mods) == 0) {

    if (control$start == 1 && control$skip.check == FALSE) {

      # Main class check --------------------------------------------------------
      if (class(x) != "netest" && class(x) != "netsim") {
        stop("x must be either an object of class netest or class netsim",
             call. = FALSE)
      }
      if (!inherits(param, "param.net")) {
        stop("param must an object of class param.net", call. = FALSE)
      }
      if (!inherits(init, "init.net")) {
        stop("init must an object of class init.net", call. = FALSE)
      }
      if (!inherits(control, "control.net")) {
        stop("control must an object of class control.net", call. = FALSE)
      }

      nw <- x$fit$newnetwork

      # Defaults ----------------------------------------------------------------

      # Is status in network formation formula?
      statOnNw <- ("status" %in% get_formula_term_attr(x$formation, nw))

      # Set dependent modeling defaults if vital or status on nw
      if (is.null(control$resimulate.network)) {
        if (param$vital == TRUE | statOnNw == TRUE) {
          control$resimulate.network <- TRUE
        } else {
          control$resimulate.network <- FALSE
        }
      }

      nGroups <- length(unique(get.vertex.attribute(nw, "group")))
      nGroups <- ifelse(nGroups == 2, 2, 1)

      if (nGroups == 2 & is.null(control$pid.prefix)) {
        control$pid.prefix <- c("g1.", "g2.")
      }

      if (statOnNw == TRUE && is.null(control$attr.rules$status)) {
        control$attr.rules$status <- "s"
      }



      # Checks ------------------------------------------------------------------

      # Check that prevalence in NW attr status and initial conditions match
      if (statOnNw == TRUE) {
        nw1 <- sum(get.vertex.attribute(nw, "status") == 1)
        init1 <- sum(unlist(init[grep("i.num", names(init), value = TRUE)]))
        if ("i.num" %in% names(init) && nw1 != init1) {
          warning("Overriding init infected settings with network status attribute",
                  call. = FALSE, immediate. = TRUE)
          if (interactive()) Sys.sleep(4)
        }
      }

      # If status not in formation formula but set on original network, state that it
      #   will be ignored
      if (statOnNw == FALSE & "status" %in% names(nw$val[[1]])) {
        warning("Overriding status vertex attribute on network with init.net conditions",
                call. = FALSE, immediate. = TRUE)
        if (interactive()) Sys.sleep(4)
      }


      # Check consistency of status vector to network structure
      if (!is.null(init$status.vector)) {
        if (length(init$status.vector) != network.size(nw)) {
          stop("Length of status.vector is unequal to size of initial network")
        }
        svals <- sort(unique(init$status.vector))
        if (control$type == "SIR") {
          if (any(svals %in% c("s", "i", "r") == FALSE)) {
            stop("status.vector contains values other than \"s\", \"i\", and \"r\" ",
                 call. = FALSE)
          }
        } else {
          if (any(svals %in% c("s", "i") == FALSE)) {
            stop("status.vector contains values other than \"s\" and \"i\" ",
                 call. = FALSE)
          }
        }
      }

      # Two-group model checks for inital conditions
      if (nGroups == 2 & is.null(init$i.num.g2) &
          is.null(init$status.vector) & statOnNw == FALSE) {
        stop("Specify i.num.g2 for two-group model simulations", call. = FALSE)
      }

      # Recovery rate and initial recovered checks
      if (control$type %in% c("SIR", "SIS")) {
        if (is.null(param$rec.rate)) {
          stop("Specify rec.rate in param.net", call. = FALSE)
        }
        if (nGroups == 2 & is.null(param$rec.rate.g2)) {
          stop("Specify rec.rate.g2 in param.net", call. = FALSE)
        }
      }
      if (control$type == "SIR") {
        if (is.null(init$r.num) & is.null(init$status.vector) & statOnNw == FALSE) {
          stop("Specify r.num in init.net", call. = FALSE)
        }
        if (nGroups == 2 & is.null(init$r.num.g2) & is.null(init$status.vector) &
            statOnNw == FALSE) {
          stop("Specify r.num.g2 in init.net", call. = FALSE)
        }
      }

      # Check demographic parameters for two-group models
      if (nGroups == 2 & param$vital == TRUE) {
        if (is.null(param$a.rate.g2)) {
          stop("Specify a.rate.g2 in param.net", call. = FALSE)
        }
        if (is.null(param$ds.rate.g2)) {
          stop("Specify ds.rate.g2 in param.net", call. = FALSE)
        }
        if (is.null(param$di.rate.g2)) {
          stop("Specify di.rate.g2 in param.net", call. = FALSE)
        }
        if (control$type == "SIR") {
          if (is.null(param$dr.rate.g2)) {
            stop("Specify dr.rate.g2 in param.net", call. = FALSE)
          }
        }
      }


      ## Deprecated parameters
      bim <- grep(".FUN", names(formals(control.net)), value = TRUE)
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

    if (control$start > 1) {

      control$resimulate.network <- TRUE

      if (control$skip.check == FALSE) {
        if (class(x) != "netsim") {
          stop("x must be a netsim object if control setting start > 1",
               call. = FALSE)
        }
        if (is.null(x$attr)) {
          stop("x must contain attr to restart simulation, see save.other ",
               "control setting", call. = FALSE)
        }
        if (is.null(x$network)) {
          stop("x must contain network object to restart simulation, ",
               "see save.network control setting", call. = FALSE)
        }
        if (control$nsteps < control$start) {
          stop("control setting nsteps must be > control setting start in ",
               "restarted simulations", call. = FALSE)
        }
        if (control$start > x$control$nsteps + 1) {
          stop("control setting start must be 1 greater than the nsteps in the ",
               "prior simulation", call. = FALSE)
        }


      }

    }


    ## Assign modules based on group parameter
    if (!is.null(control$type)) {
      def <- grep(".FUN",names(control))
      args <- names(control)[def]
      flag <- length(grep(".g2",names(param)))

      if (flag == 0) {
        for (i in 1:length(args)) {
          if (is.null(control[[args[i]]])) {
            temp <- get(gsub(".FUN",".net",args[i]))
            control[[args[i]]] <- temp
          }
        }
      }
      else {
        for (i in 1:length(args)) {
          if (is.null(control[[args[i]]])) {
            temp <- get(gsub(".FUN",".2g.net",args[i]))
            control[[args[i]]] <- temp
          }
        }
      }
    }

    ## In-place assignment to update param and control
    assign("param", param, pos = parent.frame())
    assign("control", control, pos = parent.frame())
  }


  if (!is.null(control$type) && length(control$user.mods) > 0) {
    stop("Control setting 'type' must be NULL if any user-specified modules specified.",
         call. = FALSE)
  }

  if (is.null(control$type)) {
    control$type <- "SI"
  }

  if (is.null(control$type) && length(grep("rec", names(param))) != 0){
    control$type <- "SIR"
  }

  ## In-place assignment to update param and control
  assign("param", param, pos = parent.frame())
  assign("control", control, pos = parent.frame())
}
