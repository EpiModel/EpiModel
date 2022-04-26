#' @title Epidemic Parameters for Stochastic Network Models
#'
#' @description Sets the epidemic parameters for stochastic network models
#'              simulated with \code{\link{netsim}}.
#'
#' @param inf.prob Probability of infection per transmissible act between
#'        a susceptible and an infected person. In two-group models, this is the
#'        probability of infection to the group 1 nodes. This may also be a
#'        vector of probabilities, with each element corresponding to the
#'        probability in that time step of infection (see Time-Varying
#'        Parameters below).
#' @param inter.eff Efficacy of an intervention which affects the per-act
#'        probability of infection. Efficacy is defined as 1 - the relative
#'        hazard of infection given exposure to the intervention, compared to no
#'        exposure.
#' @param inter.start Time step at which the intervention starts, between 1 and
#'        the number of time steps specified in the model. This will default to
#'        1 if \code{inter.eff} is defined but this parameter is not.
#' @param act.rate Average number of transmissible acts \emph{per partnership}
#'        per unit time (see \code{act.rate} Parameter below). This may also be
#'        a vector of rates, with each element corresponding to the rate in
#'        that time step of infection (see Time-Varying Parameters below).
#' @param rec.rate Average rate of recovery with immunity (in \code{SIR} models)
#'        or re-susceptibility (in \code{SIS} models). The recovery rate is the
#'        reciprocal of the disease duration. For two-group models, this is the
#'        recovery rate for group 1 persons only. This parameter is only used
#'        for \code{SIR} and \code{SIS} models. This may also be a vector
#'        of rates, with each element corresponding to the rate in that time
#'        step of infection (see Time-Varying Parameters below).
#' @param a.rate Arrival or entry rate. For one-group models, the arrival rate
#'        is the rate of new arrivals per person per unit time. For two-group
#'        models, the arrival rate is parameterized as a rate per group 1
#'        person per unit time, with the \code{a.rate.g2} rate set as described
#'        below.
#' @param ds.rate Departure or exit rate for susceptible persons. For two-group
#'        models, it is the rate for group 1 susceptible persons only.
#' @param di.rate Departure or exit rate for infected persons. For two-group
#'        models, it is the rate for group 1 infected persons only.
#' @param dr.rate Departure or exit rate for recovered persons. For two-group
#'        models, it is the rate for group 1 recovered persons only. This
#'        parameter is only used for \code{SIR} models.
#' @param inf.prob.g2 Probability of transmission given a transmissible act
#'        between a susceptible group 2 person and an infected group 1 person.
#'        It is the probability of transmission to group 2 members.
#' @param rec.rate.g2 Average rate of recovery with immunity (in \code{SIR}
#'        models) or re-susceptibility (in \code{SIS} models) for group 2
#'        persons. This parameter is only used for two-group \code{SIR} and
#'        \code{SIS} models.
#' @param a.rate.g2 Arrival or entry rate for group 2. This may either be
#'        specified numerically as the rate of new arrivals per group 2 person
#'        per unit time, or as \code{NA}, in which case the group 1 rate,
#'        \code{a.rate}, governs the group 2 rate. The latter is used when, for
#'        example, the first group is conceptualized as female, and the female
#'        population size determines the arrival rate. Such arrivals are evenly
#'        allocated between the two groups.
#' @param ds.rate.g2 Departure or exit rate for group 2 susceptible persons.
#' @param di.rate.g2 Departure or exit rate for group 2 infected persons.
#' @param dr.rate.g2 Departure or exit rate for group 2 recovered persons. This
#'        parameter is only used for \code{SIR} model types.
#'
#' @param ... Additional arguments passed to model.
#'
#' @details
#' \code{param.net} sets the epidemic parameters for the stochastic network
#' models simulated with the \code{\link{netsim}} function. Models
#' may use the base types, for which these parameters are used, or new process
#' modules which may use these parameters (but not necessarily). A detailed
#' description of network model parameterization for base models is found in
#' the \href{http://www.epimodel.org/tut.html}{Basic Network Models} tutorial.
#'
#' For base models, the model specification will be chosen as a result of
#' the model parameters entered here and the control settings in
#' \code{\link{control.net}}. One-group and two-group models are available,
#' where the latter assumes a heterogeneous mixing between two distinct
#' partitions in the population (e.g., men and women). Specifying any two-group
#' parameters (those with a \code{.g2}) implies the simulation of a two-group
#' model. All the parameters for a desired model type must be specified, even if
#' they are zero.
#'
#' @section The \code{act.rate} Parameter:
#' A key difference between these network models and DCM/ICM classes is the
#' treatment of transmission events. With DCM and ICM, contacts or partnerships
#' are mathematically instantaneous events: they have no duration in time, and
#' thus no changes may occur within them over time. In contrast, network models
#' allow for partnership durations defined by the dynamic network model,
#' summarized in the model dissolution coefficients calculated in
#' \code{\link{dissolution_coefs}}. Therefore, the \code{act.rate} parameter has
#' a different interpretation here, where it is the number of transmissible acts
#' \emph{per partnership} per unit time.
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
#' population by the fourth time step, the third element in the vector will
#' carry forward until one of those events occurs or the simulation ends. For
#' further examples, see the \href{https://statnet.org/nme/}{NME Course
#' Tutorials}.
#'
#' @section Random Parameters:
#' In addition to deterministic parameters in either fixed or time-varying
#' varieties above, one may also include a generator for random parameters.
#' These might include a vector of potential parameter values or a statistical
#' distribution definition; in either case, one draw from the generator would
#' be completed per individual simulation. This is possible by passing a list
#' named \code{random.params} into \code{param.net}, with each element of
#' \code{random.params} a named generator function. See the help page and
#' examples in \code{\link{generate_random_params}}. A simple factory function
#' for sampling is provided with \code{\link{param_random}} but any function
#' will do.
#'
#' @section New Modules:
#' To build original models outside of the base models, new process modules
#' may be constructed to replace the existing modules or to supplement the
#' existing set. These are passed into the control settings in
#' \code{\link{control.net}}. New modules may use either the existing model
#' parameters named here, an original set of parameters, or a combination of
#' both. The \code{...} allows the user to pass an arbitrary set of original
#' model parameters into \code{param.net}. Whereas there are strict checks with
#' default modules for parameter validity, these checks are the user's
#' responsibility with new modules.
#'
#' @return An \code{EpiModel} object of class \code{param.net}.
#'
#' @seealso Use \code{\link{init.net}} to specify the initial conditions and
#'          \code{\link{control.net}} to specify the control settings. Run the
#'          parameterized model with \code{\link{netsim}}.
#'
#' @keywords parameterization
#'
#' @export
#'
#' @examples
#' ## Example SIR model parameterization with fixed and random parameters
#' # Network model estimation
#' nw <- network_initialize(n = 100)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#' est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' # Random epidemic parameter list (here act.rate values are sampled uniformly
#' # with helper function param_random, and inf.prob follows a general Beta
#' # distribution with the parameters shown below)
#' my_randoms <- list(
#'   act.rate = param_random(1:3),
#'   inf.prob = function() rbeta(1, 1, 2)
#' )
#'
#' # Parameters, initial conditions, and control settings
#' param <- param.net(rec.rate = 0.02, random.params = my_randoms)
#'
#' # Printing parameters shows both fixed and and random parameter functions
#' param
#'
#' # Set initial conditions and controls
#' init <- init.net(i.num = 10, r.num = 0)
#' control <- control.net(type = "SIR", nsteps = 10, nsims = 3, verbose = FALSE)
#'
#' # Simulate the model
#' sim <- netsim(est, param, init, control)
#'
#' # Printing the sim object shows the randomly drawn values for each simulation
#' sim
#'
#' # Parameter sets can be extracted with:
#' get_param_set(sim)
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
    for (i in seq_along(dot.args)) {
      p[[names.dot.args[i]]] <- dot.args[[i]]
    }
  }

  ## random.params checks
   if ("random.params" %in% names.dot.args) {
     for (nm in names(p[["random.params"]])) {
       if (nm %in% names(p)) {
        warning(
          "The parameter `", nm, "` is defined twice, once as fixed",
          " and once as a random parameter.\n Only the random parameter",
          " definition will be used."
        )
       }
     }
   }

  ## Defaults and Checks
  if ("b.rate" %in% names.dot.args) {
    p[["a.rate"]] <- dot.args[["b.rate"]]
    stop("EpiModel 1.7.0 onward renamed the birth rate parameter b.rate
         to a.rate. ", "See documentation for details.",
         call. = FALSE)
  }
  if ("b.rate.g2" %in% names.dot.args) {
    p[["a.rate.g2"]] <- dot.args[["b.rate.g2"]]
    stop("EpiModel 1.7.0 onward renamed the birth rate parameter b.rate.g2 to
         a.rate.g2. ", "See documentation for details.",
         call. = FALSE)
  }
  # Check for mode to group suffix change
  m2.flag <- grep(".m2", names(p))
  if (length(m2.flag) > 0) {
    names(p) <- gsub(".m2", ".g2", names(p))
    warning("EpiModel 2.0+ has updated parameter suffixes. ",
            "All .m2 parameters changed to .g2. See documentation.")

  }
  if (missing(act.rate)) {
    p[["act.rate"]] <- 1
  }
  p[["vital"]] <- ifelse(!missing(a.rate) | !missing(ds.rate) |
                      !missing(di.rate) | !missing(dr.rate), TRUE, FALSE)
  if ("act.rate.g2" %in% names.dot.args) {
    warning("act.rate.g2 parameter was entered. ",
            "If using built-in models, only act.rate parameter will apply.",
            call. = FALSE)
  }

  if (!is.null(p[["inter.eff"]]) && is.null(p[["inter.start"]])) {
    p[["inter.start"]] <- 1
  }

  ## Output
  class(p) <- c("param.net", "list")
  return(p)
}

#' @title Update Model Parameters for Stochastic Network Models
#'
#' @description Updates epidemic model parameters originally set with
#'              \code{\link{param.net}} and adds new parameters.
#'
#' @param param Object of class \code{param.net}, output from function of same
#'              name.
#' @param new.param.list Named list of new parameters to add to original
#'        parameters.
#'
#' @details
#' This function can update any original parameters specified with
#' \code{\link{param.net}} and add new parameters. This function would be used
#' if the inputs to \code{\link{param.net}} were a long list of fixed model
#' parameters that needed supplemental replacements or additions for particular
#' model runs (e.g., changing an intervention efficacy parameter but leaving all
#' other parameters fixed).
#'
#' The \code{new.param.list} object should be a named list object containing
#' named parameters matching those already in \code{x} (in which case those
#' original parameter values will be replaced) or not matching (in which case
#' new parameters will be added to \code{param}).
#'
#' @return
#' An updated list object of class \code{param.net}, which can be passed to the
#' EpiModel function \code{\link{netsim}}.
#'
#' @examples
#' x <- param.net(inf.prob = 0.5, act.rate = 2)
#' y <- list(inf.prob = 0.75, dx.rate = 0.2)
#' z <- update_params(x, y)
#' print(z)
#'
#' @export
#'
update_params <- function(param, new.param.list) {

  if (!inherits(param, "param.net")) {
    stop("x should be object of class param.net")
  }
  if (class(new.param.list) != "list") {
    stop("new.param.list should be object of class list")
  }

  for (ii in seq_along(new.param.list)) {
    param[[names(new.param.list)[ii]]] <- new.param.list[[ii]]
  }

  return(param)
}


#' @title Create a Value Sampler for Random Parameters
#'
#' @description This function returns a 0 argument function that can be used as
#'   a generator function in the \code{random.params} argument of the
#'   \code{\link{param.net}} function.
#'
#' @param values A vector of values to sample from.
#' @param prob A vector of weights to use during sampling. If \code{NULL},
#'        all values have the same probability of being picked
#'        (default = \code{NULL}).
#'
#' @return A 0 argument generator function to sample one of the values from the
#' \code{values} vector.
#'
#' @seealso \code{\link{param.net}} and \code{\link{generate_random_params}}
#' @export
#'
#' @examples
#' # Define function with equal sampling probability
#' a <- param_random(1:5)
#' a()
#'
#' # Define function with unequal sampling probability
#' b <- param_random(1:5, prob = c(0.1, 0.1, 0.1, 0.1, 0.6))
#' b()
#'
param_random <- function(values, prob = NULL) {
  if (!is.null(prob) && length(prob) != length(values)) {
    stop("incorrect number of probabilites")
  }

  f <- function() {
    return(sample(x = values, size = 1, prob = prob, replace = TRUE))
  }

  return(f)
}


#' @title Generate Values for Random Parameters
#'
#' @description This function uses the generative functions in the
#'              \code{random.params} list to create values for the parameters.
#'
#' @param param The \code{param} argument received by the \code{netsim}
#'              functions.
#' @param verbose Should the function output the generated values
#'                (default = FALSE)?
#'
#' @return A fully instantiated \code{param} list.
#'

#' @section \code{random.params}:
#' The \code{random.params} argument to the \code{\link{param.net}} function
#' must be a named list of functions that each return a value that can be used
#' as the argument with the same name. In the example below, \code{param_random}
#' is a function factory provided by EpiModel for \code{act.rate} and
#' for \code{tx.halt.part.prob} we provide bespoke functions. A function factory
#' is a function that returns a new function
#' (see https://adv-r.hadley.nz/function-factories.html).
#'
#' @section Generator Functions:
#' The functions used inside \code{random_params} must be 0 argument functions
#' returning a valid value for the parameter with the same name.
#'
#' @section \code{param_random_set}:
#' The \code{random_params} list can optionally contain a
#' \code{param_random_set} element. It must be a \code{data.frame} of possible
#' values to be used as parameters.
#'
#' The column names must correspond either to:
#' the name of one parameter, if this parameter is of size 1; or the name of one
#' parameter with "_1", "_2", etc. appended, with the number representing the
#' position of the value, if this parameter is of size > 1. This means that the
#' parameter names cannot contain any underscores "_" if you intend to use
#' \code{param_random_set}.
#'
#' The point of the \code{param.random.set} \code{data.frame} is to allow the
#' random parameters to be correlated. To achieve this, a whole row of the
#' \code{data.frame} is selected for each simulation.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' ## Example with only the generator function
#'
#' # Define random parameter list
#' my_randoms <- list(
#'   act.rate = param_random(c(0.25, 0.5, 0.75)),
#'   tx.prob = function() rbeta(1, 1, 2),
#'   stratified.test.rate = function() c(
#'     rnorm(1, 0.05, 0.01),
#'     rnorm(1, 0.15, 0.03),
#'     rnorm(1, 0.25, 0.05)
#'   )
#' )
#'
#' # Parameter model with fixed and random parameters
#' param <- param.net(inf.prob = 0.3, random.params = my_randoms)
#'
#' # Below, `tx.prob` is set first to 0.3 then assigned a random value using
#' # the function from `my_randoms`. A warning notifying of this overwrite is
#' # therefore produced.
#' param <- param.net(tx.prob = 0.3, random.params = my_randoms)
#'
#'
#' # Parameters are drawn automatically in netsim by calling the function
#' # within netsim_loop. Demonstrating draws here but this is not used by
#' # end user.
#' paramDraw <- generate_random_params(param, verbose = TRUE)
#' paramDraw
#'
#'
#' ## Addition of the `param.random.set` `data.frame`
#'
#' # This function will generate sets of correlated parameters
#'  generate_correlated_params <- function() {
#'    param.unique <- runif(1)
#'    param.set.1 <- param.unique + runif(2)
#'    param.set.2 <- param.unique * rnorm(3)
#'
#'    return(list(param.unique, param.set.1, param.set.2))
#'  }
#'
#'  # Data.frame set of random parameters :
#'  correlated_params <- t(replicate(10, unlist(generate_correlated_params())))
#'  correlated_params <- as.data.frame(correlated_params)
#'  colnames(correlated_params) <- c(
#'    "param.unique",
#'    "param.set.1_1", "param.set.1_2",
#'    "param.set.2_1", "param.set.2_2", "param.set.2_3"
#'  )
#'
#' # Define random parameter list with the `param.random.set` element
#' my_randoms <- list(
#'   act.rate = param_random(c(0.25, 0.5, 0.75)),
#'   param.random.set = correlated_params
#' )
#'
#' # Parameter model with fixed and random parameters
#' param <- param.net(inf.prob = 0.3, random.params = my_randoms)
#'
#' # Parameters are drawn automatically in netsim by calling the function
#' # within netsim_loop. Demonstrating draws here but this is not used by
#' # end user.
#' paramDraw <- generate_random_params(param, verbose = TRUE)
#' paramDraw
#'
#' }
generate_random_params <- function(param, verbose = FALSE) {
  if (is.null(param[["random.params"]]) ||
      length(param[["random.params"]]) == 0) {
    return(param)
  } else {
    random.params <- param[["random.params"]]
  }

  if (!is.list(random.params)) {
    stop("`random.params` must be named list of functions")
  }

  rng_names <- names(random.params)
  if (any(rng_names == "")) {
    stop("all elements of `random.params` must be named")
  }

  rng_values <- list()

  if ("param.random.set" %in% rng_names) {
    # Take `param.random.set` out of the `random.params` list
    param.random.set <- random.params[["param.random.set"]]
    random.params[["param.random.set"]] <- NULL
    rng_names <- names(random.params)

    if (!is.data.frame(param.random.set)) {
      stop("`param.random.set` must be a data.frame")
    }

    # Check the format of the names
    set.elements <- names(param.random.set)
    correct_format <- grepl("^[a-zA-Z0-9.]*(_[0-9]+)?$", set.elements)
    if (!all(correct_format)) {
      stop("The following column names in `param.random.set` are malformed: \n",
        paste0(set.elements[!correct_format], collapse = ", "), "\n\n",
        "you can check the names with ",
        '`grepl("^[a-zA-Z0-9.]*(_[0-9]+)?$", your.names)` \n',
        "Example: 'unique.param', 'param.set_1', 'param.set_2'"
      )
    }

    # Construct a `data.frame` matching the names with the parameters
    set.elements <- names(param.random.set)
    set.elements.split <- data.frame(do.call(
      rbind,
      strsplit(set.elements, "_")
    ))

    nums <- grepl("^[0-9]+$", set.elements.split[, 2])
    set.elements.split[, 2] <- as.numeric(
      ifelse(nums, set.elements.split[, 2], "1")
    )

    set.elements <- cbind(set.elements, set.elements.split)
    colnames(set.elements) <- c("name", "param", "position")

    # Pick one row of the `data.frame`
    sampled.row <- sample.int(nrow(param.random.set), 1)
    param.random.set <- param.random.set[sampled.row, ]

    # Set the new values in `rng_values`
    for (i in seq_len(nrow(set.elements))) {
      value <- param.random.set[1, set.elements[i, "name"][[1]]]
      parameter <- set.elements[i, "param"][[1]]
      position <- set.elements[i, "position"][[1]]

      rng_values[[parameter]][position] <- value
    }
  }

  if (!all(vapply(random.params, is.function, TRUE))) {
    stop("all elements of `random.params` must be functions \n",
         "(Except 'param.random.set')")
  }

  duplicated.rng <- names(rng_values) %in% rng_names
  if (any(duplicated.rng)) {
    warning("Some parameters are set to be randomly assigned twice: \n",
            paste0(rng_names[duplicated.rng], collapse = ", "), "\n\n",
            "The version from a generator function will be used")
  }

  rng_values[rng_names] <- lapply(random.params, do.call, args = list())
  for (nm in rng_names) {
    param[nm] <- rng_values[nm]
  }

  param[["random.params.values"]] <- rng_values

  if (verbose == TRUE) {
    msg <-
     "The following values were randomly generated for the given parameters: \n"
    msg <- c(msg, paste0("`", names(rng_values), "`: ", rng_values, "\n"))
    message(msg)
  }

  return(param)
}


#' @title Initial Conditions for Stochastic Network Models
#'
#' @description Sets the initial conditions for stochastic network models
#'              simulated with \code{netsim}.
#'
#' @param i.num Number of initial infected persons. For two-group models, this
#'        is the number of initial group 1 infected persons.
#' @param r.num Number of initial recovered persons. For two-group models, this
#'        is the number of initial group 1 recovered persons. This parameter is
#'        only used for the \code{SIR} model type.
#' @param i.num.g2 Number of initial infected persons in group 2. This parameter
#'        is only used for two-group models.
#' @param r.num.g2 Number of initial recovered persons in group 2. This
#'        parameter is only used for two-group \code{SIR} models.
#' @param status.vector A vector of length equal to the size of the input
#'        network, containing the status of each node. Setting status here
#'        overrides any inputs passed in the \code{.num} arguments.
#' @param infTime.vector A vector of length equal to the size of the input
#'        network, containing the (historical) time of infection for each of
#'        those nodes with a current status of \code{"i"}. Can only be used if
#'        \code{status.vector} is used, and must contain \code{NA} values for
#'        any nodes whose status is not \code{"i"}.
#' @param ... Additional initial conditions passed to model.
#'
#' @details
#' The initial conditions for a model solved with \code{\link{netsim}} should be
#' input into the \code{init.net} function. This function handles initial
#' conditions for both base models and new modules. For an overview of
#' specifying initial conditions across a variety of base network models,
#' consult the \href{http://www.epimodel.org/tut.html}{Basic Network Models}
#' tutorials.
#'
#' @return An \code{EpiModel} object of class \code{init.net}.
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
    for (i in seq_along(dot.args)) {
      p[[names.dot.args[i]]] <- dot.args[[i]]
    }
  }

  ## Defaults and checks

  # Checks that .m2 syntax is change to .g2
  m2.init <- grep("m2", names(p), value = TRUE)
  if (length(m2.init) > 0) {
    names(p) <- gsub(".m2", ".g2", names(p))
    warning("EpiModel 2.0+ updated initial condition suffixes. ",
            "All .m2 initial conditions changed to .g2. See documentation.")
  }
  if (!is.null(p[["i.num"]]) & !is.null(p[["status.vector"]])) {
    stop("Use i.num OR status.vector to set initial infected")
  }
  if (!is.null(p[["infTime.vector"]]) & is.null(p[["status.vector"]])) {
    stop("infTime.vector may only be used if status.vector is used")
  }
  if (!is.null(p[["infTime.vector"]]) &
      length(p[["infTime.vector"]]) != length(p[["status.vector"]])) {
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
#' @param start For models with network resimulation, time point to start up
#'        simulation. For restarted simulations, this must be one greater than
#'        the final time step in the prior simulation and must be less than the
#'        value in \code{nsteps}.
#' @param resimulate.network If \code{TRUE}, resimulate the network at each time
#'        step. This is required when the epidemic or demographic processes
#'        impact the network structure (e.g., vital dynamics).
#' @param tergmLite Logical indicating usage of either \code{tergm}
#'        (\code{tergmLite = FALSE}), or \code{tergmLite}
#'        (\code{tergmLite = TRUE}). Default of \code{FALSE}. (See the
#'        \href{https://statnet.org/tut/EpiModel2.html#tergmLite}{EpiModel 2.0
#'        migration document} for details on \code{tergmLite}.)
#' @param cumulative.edgelist If \code{TRUE}, calculates a cumulative edgelist
#'        within the network simulation module. This is used when tergmLite is
#'        used and the entire networkDynamic object is not used.
#' @param truncate.el.cuml Number of time steps of the cumulative edgelist to
#'        retain. See help file for \code{\link{update_cumulative_edgelist}} for
#'        options.
#' @param attr.rules A list containing the  rules for setting the attributes of
#'        incoming nodes, with one list element per attribute to be set (see
#'        details below).
#' @param epi.by A character vector of length 1 containing a nodal attribute for
#'        which subgroup epidemic prevalences should be calculated. This nodal
#'        attribute must be contained in the network model formation formula,
#'        otherwise it is ignored.
#' @param initialize.FUN Module to initialize the model at time 1, with the
#'        default function of \code{\link{initialize.net}}.
#' @param departures.FUN Module to simulate departure or exit, with the default
#'        function of \code{\link{departures.net}}.
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
#' @param prevalence.FUN Module to calculate disease prevalence at each time
#'        step, with the default function of \code{\link{prevalence.net}}.
#' @param verbose.FUN Module to print simulation progress to screen, with the
#'        default function of \code{\link{verbose.net}}.
#' @param module.order A character vector of module names that lists modules in
#'        the order in which they should be evaluated within each time step. If
#'        \code{NULL}, the modules will be evaluated as follows: first any
#'        new modules supplied through \code{...} in the order in which they are
#'        listed, then the built-in modules in the order in which they are
#'        listed as arguments above. \code{initialize.FUN} will
#'        always be run first and \code{verbose.FUN} will always be run last.
#' @param save.nwstats If \code{TRUE}, save network statistics in a data frame.
#'        The statistics to be saved are specified in the \code{nwstats.formula}
#'        argument.
#' @param nwstats.formula A right-hand sided ERGM formula that includes network
#'        statistics of interest, with the default to the formation formula
#'        terms.
#' @param save.network If \code{TRUE}, networkDynamic/network object is
#'        saved at simulation end. Implicit reset to \code{FALSE} if
#'        \code{tergmLite = TRUE} (because network history is not saved with
#'        tergmLite).
#' @param save.transmat If \code{TRUE}, complete transmission matrix is saved at
#'        simulation end. Default of \code{TRUE}.
#' @param save.other A character vector of elements on the \code{dat} main data
#'        list to save out after each simulation. One example for base models is
#'        the attribute list, \code{"attr"}, at the final time step.
#' @param verbose If \code{TRUE}, print model progress to the console.
#' @param verbose.int Time step interval for printing progress to console, where
#'        0 prints completion status of entire simulation and positive integer
#'        \code{x} prints progress after every \code{x} time steps. The default
#'        is to print progress after each time step.
#' @param skip.check If \code{TRUE}, skips the default error checking for the
#'        structure and consistency of the parameter values, initial conditions,
#'        and control settings before running base epidemic models. Setting
#'        this to \code{FALSE} is recommended when running models with new
#'        modules specified.
#' @param raw.output If \code{TRUE}, \code{netsim} will output a list of netsim
#'        data (one per simulation) instead of a formatted \code{netsim} object.
#' @param tergmLite.track.duration If \code{TRUE}, track duration information
#'        for models in \code{tergmLite} simulations.
#' @param set.control.ergm Control arguments passed to \code{ergm}'s
#'        \code{simulate_formula.network}.  In \code{netsim}, this is only used
#'        when initializing the network with \code{edapprox = TRUE}; all other
#'        simulations in \code{netsim} use \code{tergm}.
#' @param set.control.tergm Control arguments passed to \code{tergm}'s
#'        \code{simulate_formula.network}. See the help file for 
#'        \code{\link{netdx}} for details and examples on specifying this 
#'        parameter.
#' @param set.control.stergm Deprecated; use \code{set.control.tergm} instead.
#' @param ... Additional control settings passed to model.
#'
#' @details
#' \code{control.net} sets the required control settings for any network model
#' solved with the \code{\link{netsim}} function. Controls are required for both
#' base model types and when passing original process modules. For an overview
#' of control settings for base models, consult the
#' \href{http://www.epimodel.org/tut.html}{Basic Network Models} tutorials.
#' For all base models, the \code{type} argument is a necessary parameter
#' and it has no default.
#'
#' @section The \code{attr.rules} Argument:
#' The \code{attr.rules} parameter is used to specify the rules for how nodal
#' attribute values for incoming nodes should be set. These rules are only
#' necessary for models in which there are incoming nodes (i.e., arrivals).
#' There are three rules available for each attribute value:
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
#' For original models, one may substitute replacement module functions for any
#' of the default functions. New modules may be added to the workflow at each
#' time step by passing a module function via the \code{...} argument. Consult
#' the \href{http://www.epimodel.org/tut.html}{New Network Models} tutorials.
#' One may remove existing modules, such as \code{arrivals.FUN}, from the
#' workflow by setting the parameter value for that argument to \code{NULL}.
#'
#' @return An \code{EpiModel} object of class \code{control.net}.
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
                        resimulate.network = FALSE,
                        tergmLite = FALSE,
                        cumulative.edgelist = FALSE,
                        truncate.el.cuml = 0,
                        attr.rules,
                        epi.by,
                        initialize.FUN = initialize.net,
                        resim_nets.FUN = resim_nets,
                        infection.FUN = NULL,
                        recovery.FUN = NULL,
                        departures.FUN = NULL,
                        arrivals.FUN = NULL,
                        nwupdate.FUN = nwupdate.net,
                        prevalence.FUN = prevalence.net,
                        verbose.FUN = verbose.net,
                        module.order = NULL,
                        save.nwstats = TRUE,
                        nwstats.formula = "formation",
                        save.transmat = TRUE,
                        save.network,
                        save.other,
                        verbose = TRUE,
                        verbose.int = 1,
                        skip.check = FALSE,
                        raw.output = FALSE,
                        tergmLite.track.duration = FALSE,
                        set.control.ergm = control.simulate.formula(
                          MCMC.burnin = 2e5),
                        set.control.stergm = NULL,
                        set.control.tergm = control.simulate.formula.tergm(),
                        ...) {
  if (!missing(set.control.stergm)) {
    warning("set.control.stergm is deprecated and will be removed in a future 
             version; use set.control.tergm instead.")
  }

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

  if ("births.FUN" %in% names(dot.args)) {
    p[["arrivals.FUN"]] <- dot.args[["births.FUN"]]
    p[["births.FUN"]] <- dot.args[["births.FUN"]] <- NULL
    stop("EpiModel 1.7.0+ renamed the birth function births.FUN to
         arrivals.FUN. ", "See documentation for details.", call. = FALSE)
  }
  if ("deaths.FUN" %in% names(dot.args)) {
    p[["departures.FUN"]] <- dot.args[["deaths.FUN"]]
    p[["deaths.FUN"]] <- dot.args[["deaths.FUN"]] <- NULL
    stop("EpiModel 1.7.0+ renamed the death function deaths.FUN to
         departures.FUN. ", "See documentation for details.", call. = FALSE)
  }

  if ("depend" %in% names(dot.args)) {
    p[["resimulate.network"]] <- dot.args[["depend"]]
    stop("EpiModel >= 2.0 has renamed the control.net setting `depend` to ",
         "`resimulate.network`. Update your code accordingly.",
         call. = FALSE)
  }

  ## Module classification
  bi.mods <- grep(".FUN", names(formal.args), value = TRUE)
  p[["bi.mods"]] <- character()
  bi.nms <- bi.mods
  index <- 1
  if (is.null(p[["type"]])) {
    for (i in seq_along(bi.mods)) {
      if (!is.null(p[[bi.mods[i]]])) {
        p[["bi.mods"]][index] <- bi.mods[i]
        index <- index + 1
      }
    }
  } else{
    p[["bi.mods"]] <- bi.mods
  }
  p[["user.mods"]] <- grep(".FUN", names(dot.args), value = TRUE)
  p[["f.names"]] <- c(p[["bi.mods"]], p[["user.mods"]])

  ## Defaults and checks

  #Check whether any base modules have been redefined by user (note: must come
  # after above)
  bi.nms <- bi.nms[-which(bi.nms %in% c("initialize.FUN", "resim_nets.FUN",
                                        "verbose.FUN", "nwupdate.FUN",
                                        "prevalence.FUN"))]
  if (length(bi.nms) > 0) {
    flag1 <- logical()
    for (args in seq_along(bi.nms)) {
      if (!(is.null(p[[bi.nms[args]]]))) {
        temp1 <- get(gsub(".FUN", ".net", bi.nms[args]))
        temp2 <- p[[bi.nms[args]]]
        flag1[args] <- identical(temp1, temp2)
      }
    }

    if (!is.null(p[["type"]]) && sum(flag1, na.rm = TRUE) != length(flag1)) {
      stop("Control parameter 'type' must be null if any user defined base
           modules are present", call. = FALSE)
    }
  }

  if (!is.null(p[["type"]]) && length(p[["user.mods"]]) > 0) {
    stop("Control parameter 'type' must be null if any user specified modules
         are present", call. = FALSE)
  }

  if (is.null(p[["nsteps"]])) {
    stop("Specify nsteps", call. = FALSE)
  }

  if (missing(attr.rules)) {
    p[["attr.rules"]] <- list()
  }

  if (!is.null(p[["epi.by"]])) {
    if (length(p[["epi.by"]]) > 1) {
      stop("Length of epi.by currently limited to 1")
    } else {
      p[["epi.by"]] <- epi.by
    }
  }

  if (is.null(p[["save.network"]])) {
    p[["save.network"]] <- TRUE
  }
  if (p[["tergmLite"]] == TRUE) {
    p[["save.network"]] <- FALSE
  }
  if (p[["tergmLite"]] == TRUE & p[["resimulate.network"]] == FALSE) {
    message("Because tergmLite = TRUE, resetting resimulate.network = TRUE",
            call. = FALSE)
    p[["resimulate.network"]] <- TRUE
  }

  ## Output
  p <- set.control.class("control.net", p)
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
  check.control.class("net", "EpiModel crosscheck.net")

  if (!is.null(control[["type"]]) && length(control[["user.mods"]]) == 0) {

    if (control[["start"]] == 1 && control[["skip.check"]] == FALSE) {

      # Main class check ----------------------------------------------------
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

      # Pull network object from netest object
      if (x[["edapprox"]] == TRUE) {
        nw <- x[["fit"]][["newnetwork"]]
      } else {
        nw <- x[["fit"]][["network"]]      
      }
      
      # Defaults ------------------------------------------------------------

      # Is status in network formation formula?
      statOnNw <- get_vertex_attribute(nw, "status") %in% c("s", "i", "r")
      statOnNw <- ifelse(sum(statOnNw) > 0, TRUE, FALSE)

      nGroups <- length(unique(get_vertex_attribute(nw, "group")))
      nGroups <- ifelse(nGroups == 2, 2, 1)

      if (nGroups == 2 & is.null(control[["pid.prefix"]])) {
        control[["pid.prefix"]] <- c("g1.", "g2.")
      }

      if (statOnNw == TRUE && is.null(control[["attr.rules"]][["status"]])) {
        control[["attr.rules"]][["status"]] <- "s"
      }


      # Checks -------------------------------------------------------------

      # Check that prevalence in NW attr status and initial conditions match
      if (statOnNw == TRUE) {
        nw1 <- sum(get_vertex_attribute(nw, "status") == 1)
        init1 <- sum(unlist(init[grep("i.num", names(init), value = TRUE)]))
        if ("i.num" %in% names(init) && nw1 != init1) {
          warning("Overriding init infected settings with network
                  status attribute", call. = FALSE, immediate. = TRUE)
          if (interactive()) Sys.sleep(4)
        }
      }

      # Check consistency of status vector to network structure
      if (!is.null(init[["status.vector"]])) {
        if (length(init[["status.vector"]]) != network.size(nw)) {
          stop("Length of status.vector is unequal to size of initial network")
        }
        svals <- sort(unique(init[["status.vector"]]))
        if (control[["type"]] == "SIR") {
          if (any(svals %in% c("s", "i", "r") == FALSE)) {
            stop("status.vector contains values other than \"s\", \"i\",
                 and \"r\" ",
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
      if (nGroups == 2 & is.null(init[["i.num.g2"]]) &
          is.null(init[["status.vector"]]) & statOnNw == FALSE) {
        stop("Specify i.num.g2 for two-group model simulations", call. = FALSE)
      }

      # Recovery rate and initial recovered checks
      if (control[["type"]] %in% c("SIR", "SIS")) {
        if (is.null(param[["rec.rate"]])) {
          stop("Specify rec.rate in param.net", call. = FALSE)
        }
        if (nGroups == 2 & is.null(param[["rec.rate.g2"]])) {
          stop("Specify rec.rate.g2 in param.net", call. = FALSE)
        }
      }
      if (control[["type"]] == "SIR") {
        if (is.null(init[["r.num"]]) & is.null(init[["status.vector"]]) &
            statOnNw == FALSE) {
          stop("Specify r.num in init.net", call. = FALSE)
        }
        if (nGroups == 2 & is.null(init[["r.num.g2"]]) &
            is.null(init[["status.vector"]]) & statOnNw == FALSE) {
          stop("Specify r.num.g2 in init.net", call. = FALSE)
        }
      }

      # Check demographic parameters for two-group models
      if (nGroups == 2 & param[["vital"]] == TRUE) {
        if (is.null(param[["a.rate.g2"]])) {
          stop("Specify a.rate.g2 in param.net", call. = FALSE)
        }
        if (is.null(param[["ds.rate.g2"]])) {
          stop("Specify ds.rate.g2 in param.net", call. = FALSE)
        }
        if (is.null(param[["di.rate.g2"]])) {
          stop("Specify di.rate.g2 in param.net", call. = FALSE)
        }
        if (control[["type"]] == "SIR") {
          if (is.null(param[["dr.rate.g2"]])) {
            stop("Specify dr.rate.g2 in param.net", call. = FALSE)
          }
        }
      }

    }

    if (control[["start"]] > 1) {

      control[["resimulate.network"]] <- TRUE

      if (control[["skip.check"]] == FALSE) {
        if (class(x) != "netsim") {
          stop("x must be a netsim object if control setting start > 1",
               call. = FALSE)
        }
        if (is.null(x[["attr"]])) {
          stop("x must contain attr to restart simulation, see save.other ",
               "control setting", call. = FALSE)
        }
        if (is.null(x[["network"]])) {
          stop("x must contain network object to restart simulation, ",
               call. = FALSE)
        }
        if (control[["nsteps"]] < control[["start"]]) {
          stop("control setting nsteps must be > control setting start in ",
               "restarted simulations", call. = FALSE)
        }
        if (control[["start"]] > x[["control"]][["nsteps"]] + 1) {
          stop("control setting start must be 1 greater than the nsteps in
               the ", "prior simulation", call. = FALSE)
        }


      }

    }


    ## Assign modules based on group parameter
    if (!is.null(control[["type"]])) {
      def <- grep(".FUN", names(control))
      args <- names(control)[def]
      flag <- length(grep(".g2", names(param)))

      if (flag == 0) {
        for (i in seq_along(args)) {
          if (is.null(control[[args[i]]])) {
            temp <- get(gsub(".FUN", ".net", args[i]))
            control[[args[i]]] <- temp
            control[["f.args"]][i] <- gsub(".FUN", ".net", args[i])
          }
        }
      }
      else {
        for (i in seq_along(args)) {
          if (is.null(control[[args[i]]])) {
            temp <- get(gsub(".FUN", ".2g.net", args[i]))
            control[[args[i]]] <- temp
            control[["f.args"]][i] <- gsub(".FUN", ".2g.net", args[i])
          }
        }
      }
    }
  }

  if (!is.null(control[["type"]]) && length(control[["user.mods"]]) > 0) {
    stop("Control setting 'type' must be NULL if any user-specified modules
         specified.", call. = FALSE)
  }

  ## In-place assignment to update param and control
  assign("param", param, pos = parent.frame())
  assign("control", control, pos = parent.frame())
}

#' Parameters List for Stochastic Network Models from a Formatted Data Frame
#'
#' Sets the epidemic parameters for stochastic network models with
#' \code{\link{netsim}} using a specially formatted data frame of parameters.
#'
#' @param long.param.df A \code{data.frame} of parameters. See details for the
#'   expected format.
#'
#' @return A list object of class \code{param.net}, which can be passed to
#'   EpiModel function \code{\link{netsim}}.
#'
#' @details
#' The \code{long.param.df} first 3 columns must be:
#' 1. 'param': the name of the parameter. If it is a multi dimension parameter,
#'    prepend the name, with "_1", "_2", etc for the position in the vector.
#' 2. 'value': the value for the parameter (or the Nth position if
#' multidimensional)
#' 3. 'type': a character string containing either "numeric", "logical" or
#'    "character". The type the value will be cast to.
#'
#' Appart from these 3 columns, the \code{data.frame} can contain any number
#' of other columns. A typical use case would be to have a "details" and
#' "source" columns to document where these parameters come from.
#'
#' @export
param.net_from_table <- function(long.param.df) {
  # Checks
  if (!all(c("param", "value", "type") %in% names(long.param.df))) {
    stop(
      "The `data.frame` must contain the following 3 columns:\n",
      "'param', 'value'", " and 'type"
    )
  }
  if (!all(long.param.df[["type"]] %in% c("numeric", "logical", "character"))) {
    stop("The `type` column must contain only 'numeric', 'logical' or",
         " 'character'")
  }
  check_params_names(long.param.df[["param"]])

  # To flat params
  flat.params <- Map(
    function(g, x) g(x),
    g = lapply(paste0("as.", long.param.df[["type"]]), get),
    x = long.param.df[["value"]]
  )
  names(flat.params) <- long.param.df[["param"]]

  # To param.list
  param <- unflatten_params(flat.params)
  class(param) <- c("param.net", "list")

  return(param)
}
