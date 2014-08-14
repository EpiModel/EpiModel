
#' @title Control Settings for Deterministic Compartmental Models
#'
#' @description Sets the controls for deterministic compartmental models
#'              simulated with \code{\link{dcm}}.
#'
#' @param type disease type to be modeled, with the choice of \code{"SI"} for
#'        Susceptible-Infected diseases, \code{"SIR"} for
#'        Susceptible-Infected-Recovered diseases, and \code{"SIS"} for
#'        Susceptible-Infected-Susceptible diseases.
#' @param nsteps number of time steps to solve the model over. This must be a
#'        positive integer.
#' @param dt time unit for model solutions, with the default of 1. Model
#'        solutions for fractional time steps may be obtained by setting this to a
#'        number between 0 and 1.
#' @param odemethod ordinary differential equation (ODE) integration method, with
#'        the default of the "Runge-Kutta 4" method (see \code{\link{ode}} for
#'        other options).
#' @param new.mod if not running a built-in model type, a function with a new
#'        model to be simulated (see details).
#' @param print.mod if \code{TRUE}, print the model form to the console.
#' @param verbose if \code{TRUE}, print model progress to the console.
#' @param ... additional control settings passed to model.
#'
#' @details
#' \code{control.dcm} sets the required control settings for any deterministic
#' compartmental models solved with the \code{\link{dcm}} function. Controls are
#' required for both built-in model types and original models. For an overview of
#' control settings for built-in DCM class models, consult the
#' \href{http://statnet.org/EpiModel/vignette/Tutorial.pdf}{EpiModel Tutorial}.
#' For all built-in models, the \code{type} argument is a necessary parameter
#' and it has no default.
#'
#' @section New Model Functions:
#' The form of the model function for built-in models may be displayed with the
#' \code{print.mod} argument set to \code{TRUE}. In this case, the model will not
#' be run. These model forms may be used as templates to write original model
#' functions.
#'
#' These new models may be input and solved with \code{\link{dcm}} using the
#' \code{new.mod} argument, which requires as input a model function. Details and
#' examples are found in the \href{http://statnet.org/EpiModel/vignette/NewDCMs.html}{Solving
#' New DCMs with EpiModel} tutorial.
#'
#' @seealso Use \code{\link{param.dcm}} to specify model parameters and
#'          \code{\link{init.dcm}} to specify the initial conditions. Run the
#'          parameterized model with \code{\link{dcm}}.
#'
#' @keywords parameterization
#'
#' @export
#'
control.dcm <- function(type,
                        nsteps,
                        dt,
                        odemethod,
                        new.mod,
                        print.mod,
                        verbose,
                        ...) {

  ## Pull parameters
  out <- as.list(match.call(expand.dots = TRUE)[-1])


  ## Split lists
  out <- split_list(out)


  ## Eval args
  out <- eval_list(out)


  ## Defaults and checks
  if (is.null(out$nsteps)) {
    stop("Specify nsteps")
  }

  # Respecify dt as sequence vector
  if (is.null(out$dt)) {
    dt <- 1
  }
  out$dt <- seq(1, nsteps, dt)

  if (is.null(out$odemethod)) {
    out$odemethod <- "rk4"
  }
  if (is.null(out$new.mod)) {
    out$new.mod <- NULL
  }
  if (is.null(out$print.mod)) {
    out$print.mod <- FALSE
  }
  if (is.null(out$verbose)) {
    out$verbose <- TRUE
  }

  # Check type for integrated models
  if (is.null(out$new.mod)) {
    if (is.null(out$type) || !(out$type %in% c("SI", "SIS", "SIR"))) {
      stop("Specify type as \"SI\", \"SIS\", or \"SIR\" ")
    }
  }

  ## Output
  class(out) <- "control.dcm"
  return(out)
}


#' @title Control Settings for Stochastic Individual Contact Models
#'
#' @description Sets the controls for stochastic individual contact models
#'              simulated with \code{\link{icm}}.
#'
#' @param type disease type to be modeled, with the choice of \code{"SI"} for
#'        Susceptible-Infected diseases, \code{"SIR"} for
#'        Susceptible-Infected-Recovered diseases, and \code{"SIS"} for
#'        Susceptible-Infected-Susceptible diseases.
#' @param nsteps number of time steps to solve the model over. This must be a
#'        positive integer.
#' @param nsims number of simulations to run.
#' @param rec.rand if \code{TRUE}, use a stochastic recovery model, with the
#'        number of recovered at each time step a function of random draws from
#'        a binomial distribution with the probability equal to \code{rec.rate}.
#'        If \code{FALSE}, then a deterministic rounded count of the expectation
#'        implied by that rate.
#' @param b.rand if \code{TRUE}, use a stochastic birth model, with the
#'        number of births at each time step a function of random draws from a
#'        binomial distribution with the probability equal to the governing birth
#'        rates. If \code{FALSE}, then a deterministic rounded count of the
#'        expectation implied by those rates.
#' @param d.rand if \code{TRUE}, use a stochastic death model, with the number of
#'        deaths at each time step a function of random draws from a binomial
#'        distribution with the probability equal to the governing death rates.
#'        If \code{FALSE}, then a deterministic rounded count of the expectation
#'        implied by those rates.
#' @param initialize.FUN module to initialize the model at the outset, with the
#'        default function of \code{\link{initialize.icm}}.
#' @param infection.FUN module to simulate disease infection, with the default
#'        function of \code{\link{infection.icm}}.
#' @param recovery.FUN module to simulate disease recovery, with the default
#'        function of \code{\link{recovery.icm}}.
#' @param deaths.FUN module to simulate deaths or exits, with the default
#'        function of \code{\link{deaths.icm}}.
#' @param births.FUN module to simulate births or entries, with the default
#'        function of \code{\link{births.icm}}.
#' @param get_prev.FUN module to calculate disease prevalence at each time step,
#'        with the default function of \code{\link{get_prev.icm}}.
#' @param verbose if \code{TRUE}, print model progress to the console.
#' @param verbose.int time step interval for printing progress to console, where
#'        0 (the default) prints completion status of entire simulation and
#'        positive integer \code{x} prints progress after each \code{x} time
#'        steps.
#' @param ... additional control settings passed to model.
#'
#' @details
#' \code{control.icm} sets the required control settings for any stochastic
#' individual contact model solved with the \code{\link{icm}} function. Controls
#' are required for both built-in model types and when passing original process
#' modules. For an overview of control settings for built-in ICM class models,
#' consult the \href{http://statnet.org/EpiModel/vignette/Tutorial.pdf}{EpiModel
#' Tutorial}. For all built-in models, the \code{type} argument is a necessary
#' parameter and it has no default.
#'
#' @section New Modules:
#' Built-in ICM models use a set of module functions that specify
#' how the individual agents in the population are subjected to infection, recovery,
#' demographics, and other processes. Core modules are those listed in the
#' \code{.FUN} arguments. For each module, there is a default function used in
#' the simulation. The default infection module, for example, is contained in
#' the \code{\link{infection.icm}} function.
#'
#' For original models, one may substitute replacement module functions for any
#' the default functions. New modules may be added to the workflow at each time
#' step by passing a module function via the \code{...} argument. Further details
#' and examples of passing new modules to \code{icm} are found in the
#' \href{http://statnet.org/EpiModel/vignette/NewICMs.html}{Solving New ICMs with
#' EpiModel} tutorial.
#'
#' @seealso Use \code{\link{param.icm}} to specify model parameters and
#'          \code{\link{init.icm}} to specify the initial conditions. Run the
#'          parameterized model with \code{\link{icm}}.
#'
#' @keywords parameterization
#'
#' @export
#'
control.icm <- function(type,
                        nsteps,
                        nsims,
                        rec.rand,
                        b.rand,
                        d.rand,
                        initialize.FUN,
                        infection.FUN,
                        recovery.FUN,
                        deaths.FUN,
                        births.FUN,
                        get_prev.FUN,
                        verbose,
                        verbose.int,
                        ...) {

  ## Pull parameters
  out <- as.list(match.call(expand.dots = TRUE)[-1])


  ## Split lists
  out <- split_list(out)


  ## Eval args
  out <- eval_list(out)


  ## Defaults and checks
  if (is.null(out$type) | !(out$type %in% c("SI", "SIS", "SIR"))) {
    stop("Specify type as \"SI\", \"SIS\", or \"SIR\" ", call. = FALSE)
  }
  if (is.null(out$nsteps)) {
    stop("Specify nsteps", call. = FALSE)
  }
  if (is.null(out$nsims)) {
    out$nsims <- 1
  }
  if (is.null(out$rec.rand)) {
    out$rec.rand <- TRUE
  }
  if (is.null(out$b.rand)) {
    out$b.rand <- TRUE
  }
  if (is.null(out$d.rand)) {
    out$d.rand <- TRUE
  }
  if (is.null(out$initialize.FUN)) {
    out$initialize.FUN <- initialize.icm
  }
  if (is.null(out$infection.FUN)) {
    out$infection.FUN <- infection.icm
  }
  if (is.null(out$recovery.FUN)) {
    out$recovery.FUN <- recovery.icm
  }
  if (is.null(out$deaths.FUN)) {
    out$deaths.FUN <- deaths.icm
  }
  if (is.null(out$births.FUN)) {
    out$births.FUN <- births.icm
  }
  if (is.null(out$get_prev.FUN)) {
    out$get_prev.FUN <- get_prev.icm
  }
  if (missing(verbose)) {
    out$verbose <- TRUE
  }
  if (is.null(out$verbose.int)) {
    out$verbose.int <- 0
  }


  ## Output
  class(out) <- "control.icm"
  return(out)
}


#' @title Control Settings for Stochastic Network Models
#'
#' @description Sets the controls for stochastic network models simulated with
#'              \code{\link{netsim}}.
#'
#' @param type disease type to be modeled, with the choice of \code{"SI"} for
#'        Susceptible-Infected diseases, \code{"SIR"} for
#'        Susceptible-Infected-Recovered diseases, and \code{"SIS"} for
#'        Susceptible-Infected-Susceptible diseases.
#' @param nsteps number of time steps to simulate the model over. This must be a
#'        positive integer.
#' @param nsims the total number of disease simulations.
#' @param depend if \code{TRUE}, resimulate the network at each time step. This
#'        occurs by default with two varieties of dependent models: if there are
#'        any vital dynamic parameters in the model, or if the network model
#'        formation formula includes the "status" attribute.
#' @param rec.rand if \code{TRUE}, use a stochastic recovery model, with the
#'        number of recovered at each time step a function of random draws from
#'        a binomial distribution with the probability equal to \code{rec.rate}.
#'        If \code{FALSE}, then a deterministic rounded count of the expectation
#'        implied by that rate.
#' @param b.rand if \code{TRUE}, use a stochastic birth model, with the
#'        number of births at each time step a function of random draws from a
#'        binomial distribution with the probability equal to the governing birth
#'        rates. If \code{FALSE}, then a deterministic rounded count of the
#'        expectation implied by those rates.
#' @param d.rand if \code{TRUE}, use a stochastic death model, with the number of
#'        deaths at each time step a function of random draws from a binomial
#'        distribution with the probability equal to the governing death rates.
#'        If \code{FALSE}, then a deterministic rounded count of the expectation
#'        implied by those rates.
#' @param tea.status if \code{TRUE}, use a temporally extended attribute (TEA)
#'        to store disease status. A TEA is needed for plotting static networks
#'        at different time steps and for animating dynamic networks with evolving
#'        status. TEAs are computationally inefficient for large simulations and
#'        should be toggled off in those cases.
#' @param attr.rules a list containing the  rules for setting the attributes of
#'        incoming nodes, with one list element per attribute to be set (see
#'        details below).
#' @param epi.by a character vector of length 1 containing a nodal attribute for
#'        which subgroup epidemic prevalences should be calculated. This nodal
#'        attribute must be contained in the network model formation formula,
#'        otherwise it is ignored.
#' @param pid.prefix a character vector of length 2 containing the prefixes for
#'        persistent ids initialized in bipartite networks with vital dynamics,
#'        with the default of \code{c("F", "M")}.
#' @param initialize.FUN module to initialize the model at time 1, with the
#'        default function of \code{\link{initialize.net}}.
#' @param infection.FUN module to simulate disease infection, with the default
#'        function of \code{\link{infection.net}}.
#' @param recovery.FUN module to simulate disease recovery, with the default
#'        function of \code{\link{recovery.net}}.
#' @param deaths.FUN module to simulate death or exit, with the default function
#'        of \code{\link{deaths.net}}.
#' @param births.FUN module to simulate births or entries, with the default
#'        function of \code{\link{births.net}}.
#' @param resim_nets.FUN module to resimulate the network at each time step,
#'        with the default function of \code{\link{resim_nets}}
#' @param get_prev.FUN module to calculate disease prevalence at each time step,
#'        with the default function of \code{\link{get_prev.net}}.
#' @param set.control.stergm control arguments passed to simulate.stergm. See the
#'        help file for \code{\link{netdx}} for details and examples on specifying
#'        this parameter.
#' @param save.nwstats if \code{TRUE}, save network statistics in a data frame.
#'        The statistics to be saved are specified in the \code{nwstats.formula}
#'        argument.
#' @param nwstats.formula a right-hand sided ERGM formula that includes network
#'        statistics of interest, with the default to the formation formula terms.
#' @param delete.nodes if \code{TRUE}, delete inactive nodes from the network
#'        after each time step, otherwise deactivate them but keep them in the
#'        network object. Deleting nodes increases computational efficiency in
#'        large network simulations.
#' @param save.transmat if \code{TRUE}, save a transmission matrix for each
#'        simulation. This object contains one row for each transmission event
#'        (see \code{\link{discord_edgelist}}).
#' @param save.network if \code{TRUE}, save a \code{networkDynamic} object
#'        containing full edge history for each simulation. If \code{delete.nodes}
#'        is set to \code{TRUE}, this will only contain a static network with the
#'        edge configuration at the final time step of each simulation.
#' @param verbose if \code{TRUE}, print model progress to the console.
#' @param verbose.int time step interval for printing progress to console, where
#'        0 prints completion status of entire simulation and positive integer
#'        \code{x} prints progress after each \code{x} time steps. The default
#'        is to print progress after each time step.
#' @param ... additional control settings passed to model.
#'
#' @details
#' \code{control.net} sets the required control settings for any network model
#' solved with the \code{\link{netsim}} function. Controls are required for both
#' built-in model types and when passing original process modules. For an overview
#' of control settings for built-in network models, consult the
#' \href{http://statnet.org/EpiModel/vignette/Tutorial.pdf}{EpiModel Tutorial}.
#' For all built-in models, the \code{type} argument is a necessary parameter
#' and it has no default.
#'
#' @section The attr.rules Argument:
#' The \code{attr.rules} parameter is used to specify the rules for how nodal
#' attribute values for incoming nodes should be set. These rules are only
#' necessary for models in which there are incoming nodes (i.e., births) and also
#' there is a nodal attribute in the network model formation formula set in
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
#' \code{attr.rules = list(race = "t1", sex = "current", status = 0)}
#' specifies how the race, sex, and status attributes should be set for incoming
#' nodes. By default, the rule is "current" for all attributes except status,
#' in which case it is 0 (that is, all incoming nodes are susceptible).
#'
#' @section New Modules:
#' Built-in network models use a set of module functions that specify how the
#' individual nodes in the network are subjected to infection, recovery,
#' demographics, and other processes. Core modules are those listed in the
#' \code{.FUN} arguments. For each module, there is a default function used in
#' the simulation. The default infection module, for example, is contained in
#' the \code{\link{infection.net}} function.
#'
#' For original models, one may substitute replacement module functions for any
#' the default functions. New modules may be added to the workflow at each time
#' step by passing a module function via the \code{...} argument. This functionality
#' is new to \code{EpiModel} and further documentation will be posted in tutorial
#' vignettes at the \href{http://statnet.org/trac/wiki/EpiModel}{EpiModel website}.
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
                        nsims,
                        depend,
                        rec.rand,
                        b.rand,
                        d.rand,
                        tea.status,
                        attr.rules,
                        epi.by,
                        pid.prefix,
                        initialize.FUN,
                        infection.FUN,
                        recovery.FUN,
                        deaths.FUN,
                        births.FUN,
                        resim_nets.FUN,
                        get_prev.FUN,
                        set.control.stergm,
                        save.nwstats,
                        nwstats.formula,
                        delete.nodes,
                        save.transmat,
                        save.network,
                        verbose,
                        verbose.int,
                        ...) {

  ## Pull parameters
  out <- as.list(match.call(expand.dots = TRUE)[-1])


  ## Split lists
  out <- split_list(out, exclude = "attr.rules")


  ## Eval args
  out <- eval_list(out)


  ## Defaults and checks
  if (is.null(out$type) | !(out$type %in% c("SI", "SIS", "SIR"))) {
    stop("Specify type as \"SI\", \"SIS\", or \"SIR\" ")
  }
  if (is.null(out$nsteps)) {
    stop("Specify nsteps")
  }
  if (is.null(out$nsims)) {
    out$nsims <- 1
  }
  if (is.null(out$rec.rand)) {
    out$rec.rand <- TRUE
  }
  if (is.null(out$b.rand)) {
    out$b.rand <- TRUE
  }
  if (is.null(out$d.rand)) {
    out$d.rand <- TRUE
  }
  if (is.null(out$tea.status)) {
    out$tea.status <- TRUE
  }
  if (missing(attr.rules)) {
    out$attr.rules <- list()
  }
  if (!is.null(out$epi.by)) {
    if (length(out$epi.by) > 1) {
      stop("Length of epi.by currently limited to 1")
    } else {
      out$epi.by <- epi.by
    }
  }
  if (is.null(out$initialize.FUN)) {
    out$initialize.FUN <- initialize.net
  }
  if (is.null(out$infection.FUN)) {
    out$infection.FUN <- infection.net
  }
  if (is.null(out$recovery.FUN)) {
    out$recovery.FUN <- recovery.net
  }
  if (is.null(out$deaths.FUN)) {
    out$deaths.FUN <- deaths.net
  }
  if (is.null(out$births.FUN)) {
    out$births.FUN <- births.net
  }
  if (is.null(out$resim_nets.FUN)) {
    out$resim_nets.FUN <- resim_nets
  }
  if (is.null(out$get_prev.FUN)) {
    out$get_prev.FUN <- get_prev.net
  }
  if (is.null(out$set.control.stergm)) {
    out$set.control.stergm <- control.simulate.network(MCMC.burnin.min = 1000)
  }
  if (is.null(out$save.nwstats)) {
    out$save.nwstats <- TRUE
  }
  if (is.null(out$delete.nodes)) {
    out$delete.nodes <- FALSE
  }
  if (out$delete.nodes == TRUE) {
    out$tea.status <- FALSE
  }
  if (is.null(out$save.transmat)) {
    out$save.transmat <- TRUE
  }
  if (is.null(out$save.network)) {
    out$save.network <- TRUE
  }
  if (missing(verbose)) {
    out$verbose <- TRUE
  }
  if (is.null(out$verbose.int)) {
    out$verbose.int <- 1
  }


  ## Output
  class(out) <- "control.net"
  return(out)
}
