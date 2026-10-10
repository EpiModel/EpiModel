#' @title Parameters That Vary With the Duration of Infection
#'
#' @description Marks transmission probabilities, act rates, or recovery rates
#'              as varying with the time since infection, for the `inf.prob`,
#'              `inf.prob.g2`, `act.rate`, `rec.rate`, and `rec.rate.g2`
#'              arguments of [param.net()]. The values are given either by
#'              stage, with the length of each stage in time steps, or one per
#'              time step since infection. The built-in [netsim()] modules
#'              apply the value at each infected node's current duration of
#'              infection.
#'
#' @param x A numeric vector of values. With `durations`, one value per stage
#'        of infection, optionally named for the stages (such as
#'        `c(acute = 0.25, chronic = 0.02)`). Without `durations`, one value per
#'        time step since infection: the first element applies in the first
#'        time step of an infection, the second in the second, and so on, and
#'        the last element carries forward for the rest of the infection.
#' @param durations Optional numeric vector of the same length as `x`, giving
#'        the length of each stage in time steps. Every stage but the last must
#'        last a whole number of time steps, one or more. The last stage lasts
#'        until the infection ends, so its duration must be `Inf`; for a value
#'        that changes after a finite last stage, add a stage with the value
#'        that follows.
#'
#' @details
#' A duration of infection is the number of time steps since the infected
#' node was infected, counted from 1: a node infected at time step 10 is in
#' its first step of infection at time steps 10 and 11, its second at time
#' step 12, and so on. For transmission, the duration is that of the infected
#' partner; for recovery, that of the recovering node.
#'
#' The two forms describe the same thing.
#' `by_infection_duration(c(acute = 0.25, chronic = 0.02), durations = c(10, Inf))`
#' gives 0.25 in the first ten time steps of an infection and 0.02 from the
#' eleventh on, which can represent an acute stage of high infectiousness;
#' `by_infection_duration(c(rep(0.25, 10), 0.02))` is the same parameter given
#' per time step. The stage form suits parameters reported by stage of
#' infection; the per-step form suits a profile computed elsewhere, such as an
#' infectiousness curve.
#'
#' This is variation over the course of each infection, not over calendar
#' time: to change a parameter at a given time step of the simulation, use
#' scenarios or parameter updaters (see `vignette("model-parameters", package
#' = "EpiModel")`).
#'
#' In multi-layer models, an entry of a [multilayer()] parameter may itself be
#' a `by_infection_duration` object:
#' `inf.prob = multilayer(by_infection_duration(c(0.5, 0.1)), 0.05)`.
#' Arithmetic on the object, such as `x * 2`, scales the values and keeps the
#' stages.
#'
#' With the built-in modules, a plain vector passed to these arguments stops
#' [netsim()] with an error, since a plain vector carries no record of what its
#' positions mean; custom modules keep their own reading of plain vectors, such
#' as values by group or by layer. Parameter tables cannot represent these
#' objects yet, so [param.net_to_table()] stops on them.
#'
#' @return An object of class `by_infection_duration`: the numeric vector of
#'         values `x`, with the stage durations in its `durations` attribute.
#'         Without `durations`, every stage lasts one time step but the last.
#'
#' @seealso [param.net()], [multilayer()].
#'
#' @export
#'
#' @examples
#' # An acute stage of ten time steps with a higher transmission probability
#' inf.prob <- by_infection_duration(c(acute = 0.25, chronic = 0.02),
#'                                   durations = c(10, Inf))
#' inf.prob
#'
#' # Recovery impossible in the first 20 time steps of an infection, then
#' # certain
#' by_infection_duration(c(0, 1), durations = c(20, Inf))
#'
#' # One value per time step, here an infectiousness profile
#' by_infection_duration(round(dgamma(1:15, shape = 3, rate = 0.5), 3))
#'
#' param <- param.net(inf.prob = inf.prob, act.rate = 1)
#' param
#'
#' \dontrun{
#' nw <- network_initialize(n = 100)
#' est <- netest(nw, formation = ~edges, target.stats = 50,
#'               coef.diss = dissolution_coefs(~offset(edges), 10),
#'               verbose = FALSE)
#' sim <- netsim(est, param, init.net(i.num = 10),
#'               control.net(type = "SI", nsteps = 25, nsims = 1,
#'                           verbose = FALSE))
#' tm <- get_transmat(sim)
#' table(tm$infDur, tm$transProb)
#' }
#'
by_infection_duration <- function(x, durations = NULL) {
  if (!is.numeric(x) || length(x) == 0 || anyNA(x)) {
    stop("`x` must be a numeric vector of length one or more, without ",
         "missing values.", call. = FALSE)
  }
  n <- length(x)
  if (is.null(durations)) {
    durations <- c(rep(1, n - 1), Inf)
  } else {
    if (!is.numeric(durations) || length(durations) != n || anyNA(durations)) {
      stop("`durations` must be a numeric vector with one duration per ",
           "value of `x`.", call. = FALSE)
    }
    if (is.finite(durations[n])) {
      stop("The last stage lasts until the infection ends, so its duration ",
           "must be Inf. For a value that changes after a finite stage, add ",
           "a stage with the value that follows, as in durations = c(10, Inf).",
           call. = FALSE)
    }
    finite <- durations[-n]
    if (any(!is.finite(finite) | finite < 1 | finite != round(finite))) {
      stop("Every stage but the last must last a whole number of time steps, ",
           "one or more.", call. = FALSE)
    }
  }
  structure(as.numeric(unclass(x)), names = names(x),
            durations = as.numeric(durations), class = "by_infection_duration")
}

# The stage durations of a by_infection_duration() object, defaulting to one
# time step per stage for an object without the attribute.
stage_durations <- function(x) {
  d <- attr(x, "durations")
  if (is.null(d)) {
    d <- c(rep(1, length(x) - 1), Inf)
  }
  d
}

# Whether a by_infection_duration() object is in the per-step form: every
# stage but the last lasts one time step.
is_per_step <- function(x) {
  d <- stage_durations(x)
  all(d[-length(d)] == 1)
}

#' @export
print.by_infection_duration <- function(x, ...) {
  values <- as.numeric(unclass(x))
  if (is_per_step(x)) {
    cat("Values by duration of infection (time steps since infection),",
        "the last carried forward:\n")
    print(values, ...)
    return(invisible(x))
  }
  d <- stage_durations(x)
  first <- cumsum(c(1, d[-length(d)]))
  last <- first + d - 1
  steps <- ifelse(is.finite(last), paste0(first, "-", last), paste0(first, "+"))
  steps[d == 1] <- as.character(first[d == 1])
  stage <- if (is.null(names(x))) as.character(seq_along(values)) else names(x)
  cat("Values by stage of infection (time steps since infection):\n")
  print(data.frame(Stage = stage, Steps = steps, Value = values,
                   stringsAsFactors = FALSE),
        row.names = FALSE, ...)
  invisible(x)
}

#' @export
format.by_infection_duration <- function(x, ...) {
  values <- structure(as.numeric(unclass(x)), names = names(x))
  out <- paste0("by_infection_duration(",
                paste(deparse(values), collapse = ""))
  if (!is_per_step(x)) {
    out <- paste0(out, ", durations = ",
                  paste(deparse(stage_durations(x)), collapse = ""))
  }
  paste0(out, ")")
}

# The values of a by_infection_duration() object as one value per time step
# since infection, through the end of its last finite stage plus one step for
# the last stage. Initialization uses its mean to backdate infection times.
as_duration_vector <- function(x) {
  if (!inherits(x, "by_infection_duration")) {
    return(x)
  }
  d <- stage_durations(x)
  d[length(d)] <- 1
  rep(as.numeric(unclass(x)), d)
}

# The value of a parameter for each infected node given its duration of
# infection: a scalar applies at every duration, and a by_infection_duration()
# object gives the value of the stage the duration falls in, the last stage
# lasting until the infection ends. `name` and `layer` label the error for a
# plain vector.
infection_duration_value <- function(x, infDur, name, layer = NULL) {
  if (inherits(x, "by_infection_duration")) {
    values <- as.numeric(unclass(x))
    ends <- cumsum(stage_durations(x))
    ends <- ends[-length(ends)]
    return(values[1 + findInterval(infDur - 1, ends)])
  }
  if (length(x) != 1) {
    stop_plain_duration_vector(name, x, layer)
  }
  rep(x, length(infDur))
}

# As infection_duration_value(), for a parameter that may also vary by layer
# through multilayer(): each discordant edge gets the entry for its layer.
layer_duration_value <- function(x, network, infDur, name) {
  if (!inherits(x, "multilayer")) {
    return(infection_duration_value(x, infDur, name))
  }
  out <- numeric(length(network))
  for (k in unique(network)) {
    idx <- which(network == k)
    out[idx] <- infection_duration_value(x[[k]], infDur[idx], name, layer = k)
  }
  out
}

stop_plain_duration_vector <- function(name, x, layer = NULL) {
  stop("`", name, "`", if (!is.null(layer)) paste0(" (layer ", layer, ")"),
       " is a plain vector of length ", length(x), ". The built-in modules ",
       "no longer read a plain vector as varying with the duration of ",
       "infection; wrap it in by_infection_duration() to give one value per ",
       "time step since infection.", call. = FALSE)
}

# Check the parameters that the built-in netsim modules read by duration of
# infection: each must be a single value or a by_infection_duration() object,
# and inf.prob, inf.prob.g2, and act.rate may also be multilayer() objects of
# those. Called by crosscheck.net() only when a built-in model type is used,
# since custom modules give plain vectors their own meaning.
check_duration_params <- function(param) {
  for (name in c("inf.prob", "inf.prob.g2", "act.rate", "rec.rate",
                 "rec.rate.g2")) {
    x <- param[[name]]
    if (is.null(x)) {
      next
    }
    if (inherits(x, "multilayer")) {
      if (name %in% c("rec.rate", "rec.rate.g2")) {
        stop("`", name, "` cannot be a multilayer() object: recovery does ",
             "not depend on the network layer.", call. = FALSE)
      }
      for (k in seq_along(x)) {
        check_duration_value(x[[k]], name, layer = k)
      }
    } else {
      check_duration_value(x, name)
    }
  }
  invisible(TRUE)
}

check_duration_value <- function(x, name, layer = NULL) {
  if (inherits(x, "by_infection_duration")) {
    return(invisible(TRUE))
  }
  if (!is.numeric(x)) {
    stop("`", name, "`", if (!is.null(layer)) paste0(" (layer ", layer, ")"),
         " must be numeric.", call. = FALSE)
  }
  if (length(x) != 1) {
    stop_plain_duration_vector(name, x, layer)
  }
  invisible(TRUE)
}

# by_infection_duration() objects apply to netsim only; DCM and ICM read
# vector parameters differently (sensitivity runs for DCMs).
stop_if_by_infection_duration <- function(param, model) {
  bid <- names(param)[vapply(param, inherits, logical(1),
                             "by_infection_duration")]
  if (length(bid) > 0) {
    stop("by_infection_duration() parameters apply to network models ",
         "simulated with netsim() only, not to ", model, " models: `",
         paste(bid, collapse = "`, `"), "`.", call. = FALSE)
  }
  invisible(TRUE)
}
