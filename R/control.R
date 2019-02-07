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
#' \href{http://statnet.github.io/tut/BasicDCMs.html}{Basic DCMs} tutorial.
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
#' examples are found in the \href{http://statnet.github.io/tut/NewDCMs.html}{Solving
#' New DCMs} tutorial.
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


#' @title Control Settings for Stochastic Individual Contact Models
#'
#' @description Sets the controls for stochastic individual contact models
#'              simulated with \code{\link{icm}}.
#'
#' @param type Disease type to be modeled, with the choice of \code{"SI"} for
#'        Susceptible-Infected diseases, \code{"SIR"} for
#'        Susceptible-Infected-Recovered diseases, and \code{"SIS"} for
#'        Susceptible-Infected-Susceptible diseases.
#' @param nsteps Number of time steps to solve the model over. This must be a
#'        positive integer.
#' @param nsims Number of simulations to run.
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
#' @param initialize.FUN Module to initialize the model at the outset, with the
#'        default function of \code{\link{initialize.icm}}.
#' @param infection.FUN Module to simulate disease infection, with the default
#'        function of \code{\link{infection.icm}}.
#' @param recovery.FUN Module to simulate disease recovery, with the default
#'        function of \code{\link{recovery.icm}}.
#' @param departures.FUN Module to simulate departures or exits, with the default
#'        function of \code{\link{departures.icm}}.
#' @param arrivals.FUN Module to simulate arrivals or entries, with the default
#'        function of \code{\link{arrivals.icm}}.
#' @param get_prev.FUN Module to calculate disease prevalence at each time step,
#'        with the default function of \code{\link{get_prev.icm}}.
#' @param verbose If \code{TRUE}, print model progress to the console.
#' @param verbose.int Time step interval for printing progress to console, where
#'        0 (the default) prints completion status of entire simulation and
#'        positive integer \code{x} prints progress after each \code{x} time
#'        steps.
#' @param skip.check If \code{TRUE}, skips the default error checking for the
#'        structure and consistency of the parameter values, initial conditions,
#'        and control settings before running base epidemic models. Setting
#'        this to \code{FALSE} is recommended when running models with new modules
#'        specified.
#' @param ... Additional control settings passed to model.
#'
#' @details
#' \code{control.icm} sets the required control settings for any stochastic
#' individual contact model solved with the \code{\link{icm}} function. Controls
#' are required for both base model types and when passing original process
#' modules. For an overview of control settings for base ICM class models,
#' consult the \href{http://statnet.github.io/tut/BasicICMs.html}{Basic ICMs}
#' tutorial. For all base models, the \code{type} argument is a necessary
#' parameter and it has no default.
#'
#' @section New Modules:
#' Base ICM models use a set of module functions that specify
#' how the individual agents in the population are subjected to infection, recovery,
#' demographics, and other processes. Core modules are those listed in the
#' \code{.FUN} arguments. For each module, there is a default function used in
#' the simulation. The default infection module, for example, is contained in
#' the \code{\link{infection.icm}} function.
#'
#' For original models, one may substitute replacement module functions for any of
#' the default functions. New modules may be added to the workflow at each time
#' step by passing a module function via the \code{...} argument.
#'
#' @seealso Use \code{\link{param.icm}} to specify model parameters and
#'          \code{\link{init.icm}} to specify the initial conditions. Run the
#'          parameterized model with \code{\link{icm}}.
#'
#' @keywords parameterization
#'
#' @export
#'
control.icm <- function(type, nsteps, nsims = 1, rec.rand = TRUE, a.rand = TRUE,
                        d.rand = TRUE, initialize.FUN = initialize.icm,
                        infection.FUN = NULL, recovery.FUN = NULL,
                        departures.FUN = NULL, arrivals.FUN = NULL,
                        get_prev.FUN = NULL, verbose = FALSE,
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
  
  
  ## Module classification
  p$bi.mods <- grep(".FUN", names(formal.args), value = TRUE)
  p$user.mods <- grep(".FUN", names.dot.args, value = TRUE)
  
  
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
#' @param depend If \code{TRUE}, resimulate the network at each time step. This
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
#' @param tea.status If \code{TRUE}, use a temporally extended attribute (TEA)
#'        to store disease status. A TEA is needed for plotting static networks
#'        at different time steps and for animating dynamic networks with evolving
#'        status. TEAs are computationally inefficient for large simulations and
#'        should be toggled off in those cases. This argument automatically set
#'        to \code{FALSE} if \code{delete.nodes=TRUE}.
#' @param attr.rules A list containing the  rules for setting the attributes of
#'        incoming nodes, with one list element per attribute to be set (see
#'        details below).
#' @param epi.by A character vector of length 1 containing a nodal attribute for
#'        which subgroup epidemic prevalences should be calculated. This nodal
#'        attribute must be contained in the network model formation formula,
#'        otherwise it is ignored.
#' @param use.pids If \code{TRUE}, use persistent ids for vertices; otherwise,
#'        numeric ids will be recycled in models with vital dynamics. For one-mode
#'        simulations, this will be a random hexidecimal value; for bipartite
#'        simulations, it will be based on \code{pid.prefix}.
#' @param pid.prefix For bipartite network simulations with vital dynamics,
#'        a character vector of length 2 containing the prefixes, with the
#'        default of \code{c("F", "M")}.
#' @param initialize.FUN Module to initialize the model at time 1, with the
#'        default function of \code{\link{initialize.net}}.
#' @param departures.FUN Module to simulate departure or exit, with the default function
#'        of \code{\link{departures.net}}.
#' @param arrivals.FUN Module to simulate arrivals or entries, with the default
#'        function of \code{\link{arrivals.net}}.
#' @param recovery.FUN Module to simulate disease recovery, with the default
#'        function of \code{\link{recovery.net}}.
#' @param edges_correct.FUN Module to adjust the edges coefficient in response
#'        to changes to the population size, with the default function of
#'        \code{\link{edges_correct}} that preserves mean degree.
#' @param resim_nets.FUN Module to resimulate the network at each time step,
#'        with the default function of \code{\link{resim_nets}}.
#' @param infection.FUN Module to simulate disease infection, with the default
#'        function of \code{\link{infection.net}}.
#' @param get_prev.FUN Module to calculate disease prevalence at each time step,
#'        with the default function of \code{\link{get_prev.net}}.
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
#' @param delete.nodes If \code{TRUE}, delete inactive nodes from the network
#'        after each time step, otherwise deactivate them but keep them in the
#'        network object. Deleting nodes increases computational efficiency in
#'        large network simulations.
#' @param save.transmat If \code{TRUE}, save a transmission matrix for each
#'        simulation. This object contains one row for each transmission event
#'        (see \code{\link{discord_edgelist}}).
#' @param save.network If \code{TRUE}, save a \code{networkDynamic} object
#'        containing full edge history for each simulation. If \code{delete.nodes}
#'        is set to \code{TRUE}, this will only contain a static network with the
#'        edge configuration at the final time step of each simulation.
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
control.net <- function(type, nsteps, start = 1, nsims = 1, ncores = 1,
                        depend, rec.rand = TRUE,
                        a.rand = TRUE, d.rand = TRUE, tea.status = TRUE,
                        attr.rules, epi.by, use.pids = TRUE, pid.prefix,
                        initialize.FUN = initialize.net, departures.FUN = NULL,
                        arrivals.FUN = NULL, recovery.FUN = NULL,
                        edges_correct.FUN = edges_correct,
                        resim_nets.FUN = resim_nets, infection.FUN = NULL,
                        get_prev.FUN = NULL, verbose.FUN = verbose.net,
                        module.order = NULL, set.control.ergm, set.control.stergm,
                        save.nwstats = TRUE, nwstats.formula = "formation",
                        delete.nodes = FALSE, save.transmat = TRUE,
                        save.network = TRUE, save.other, verbose = TRUE,
                        verbose.int = 1, skip.check = FALSE, ...) {
  
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
  
  ## Module classification
  bi.mods <- grep(".FUN", names(formal.args), value = TRUE)
  p$bi.mods <- bi.mods
  p$user.mods <- grep(".FUN", names.dot.args, value = TRUE)
  
  
  if (missing(depend)) {
    arg.list <- as.list(match.call())
    if ((!is.null(arg.list$departures.FUN) && arg.list$departures.FUN != "departures.net") |
        (!is.null(arg.list$arrivals.FUN) && arg.list$arrivals.FUN != "arrivals.net")) {
      p$depend <- TRUE
    }
  }
  
  ## Defaults and checks

  #Check whether any base modules have been redefined by user (note: must come after above)
  bi.nms <- p$bi.mods
  bi.nms <- bi.nms[-which(bi.nms %in% c("initialize.FUN", "edges_correct.FUN", 
                                        "resim_nets.FUN", "verbose.FUN"))]
  flag1 <- logical()
  for (args in 1:length(bi.nms)) {
    if (!(is.null(p[[bi.nms[args]]])) ) {
      temp1 <- get(gsub(".FUN",".net",bi.nms[args]))
      temp2 <- p[[bi.nms[args]]]
      flag1[args] <- identical(temp1,temp2)
    }
  }
  
  if (!is.null(p$type) && sum(flag1) != length(flag1)) {
    stop("Control parameter 'type' must be null if any user defined base modules are present")
  }
  
  
  ##Base modes: adjusting for recovery module
  #Note: likely not the best way to do this as references obj. out of scope
  
  flag.rec <- sum(grep("rec", names(param)) > 0)
  
  if ( flag.rec == TRUE) {
    flag2 <- logical()
    for (args in 1:length(bi.nms)) {
      flag2[args] <- ifelse(is.null(p[[bi.nms[args]]]), TRUE, FALSE)
    }
  
    if (is.null(p$type) && sum(flag2) != 0) {
      stop(paste("If control parameter 'type' is not specified, user must specify all base modules. 
           Missing: ",paste(bi.nms[flag2==TRUE], collapse=' '),". See modules.net for details."),
           call. = FALSE)
    }
  }
  else {
    bi.nms <- bi.nms[-which(bi.nms %in% c("recovery.FUN"))]
    flag2 <- logical()
    for (args in 1:length(bi.nms)) {
      flag2[args] <- ifelse(is.null(p[[bi.nms[args]]]), TRUE, FALSE)
    }
    
    if (is.null(p$type) && sum(flag2) != 0) {
      stop(paste("If control parameter 'type' is not specified, user must specify all base modules. 
           Missing: ",paste(bi.nms[flag2==TRUE], collapse=' '),". See modules.net for details."),
           call. = FALSE)
    }
  }

  
  if (!is.null(p$type) && length(p$user.mods) > 0) {
    stop("Control setting 'type' must be NULL if any user-specified modules specified.",
         calll. = FALSE)
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
  
  if (p$delete.nodes == TRUE) {
    p$tea.status <- FALSE
  }
  
  ## Output
  class(p) <- c("control.net", "list")
  return(p)
  }
