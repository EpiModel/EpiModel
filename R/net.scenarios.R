#' Make a list of EpiModel scenarios from a data.frame of scenarios
#'
#' An EpiModel scenario allows one or multiple set of parameters to be applied
#' to a model a predefined timesteps. They are usually used by a researcher who
#' wants to model counterfactuals using a pre calibrated model.
#'
#' @param scenarios.df a \code{data.frame}
#'
#' @return a list of EpiModel scenarios
#'
#' @section scenarios.df:
#' The \code{scenarios.df} is a \code{data.frame} of values to be used as
#' parameters.
#'
#' It must contain a ".at" column, specifying when the changes should occur.
#' It requires the "updater" module of EpiModel. *See, vignette*. If the ".at"
#' value of a row is less than two, the changes will be applied to the
#' parameter list iteself. The second mandatory column is ".scenario.id". It
#' is used to distinguish the different scenarios. If multiple rows share the
#' same ".scenario.id", the resulting scenario will contain one updater per row.
#' This permits modifying parameters at multiple points in time. (e.g. an
#' intervention limited in time).
#'
#' The other column names must correspond either to:
#' the name of one parameter if this parameter is of size 1 or the name of the
#' parameter with "_1", "_2", "_N" with the second part being the position of
#' the value for a parameter of size > 1. This means that the parameter names
#' cannot contain any underscore "_". (e.g "a.rate", "d.rate_1", "d.rate_2")
#'
#' @export
create_scenario_list <- function(scenarios.df) {
  check_scenarios_df(scenarios.df)
  scenarios.names <- unique(scenarios.df[[".scenario.id"]])

  scenarios <- lapply(
    scenarios.names,
    function(id) make_scenario(
      scenarios.df[scenarios.df[[".scenario.id"]] == id, ])
  )

  names(scenarios) <- scenarios.names

  return(scenarios)
}

#' a scenario is a list with the following elements:
#'   - `id`: the identifier of the scenario
#'   - `.param.updater.list`: as in the updater module
#'
#' when running "use_scenario", if any updater has an `at` value < 2, it is
#' applied immediatly (so before netsim)
#'
#' @noRd
make_scenario <- function(scenario.rows) {
  scenario <- list(
    id = scenario.rows[[".scenario.id"]][1],
    .param.updater.list = vector(mode = "list", length = nrow(scenario.rows))
  )

  elements.at <- scenario.rows[[".at"]]
  # ensures that `scenario.rows` is a `tibble`
  scenario.rows <- dplyr::select(scenario.rows, - c(".at", ".scenario.id"))

  for (i in seq_along(elements.at)) {
    scenario[[".param.updater.list"]][[i]] <- list(
      at = elements.at[[i]],
      param = unflatten_params(scenario.rows[i, ])
    )
  }

  return(scenario)
}

#' Apply a scenario object to a param.net object
#'
#' @param scenario a scenario object usually created from a \code{data.frame} of
#' scenarios using the \code{create_scenario_list} function. See the vignette
#' "network-model-scenarios".
#'
#' @section scenario:
#' A scenario is a list containing an "id" field, the name of the scenario and
#' a ".param.updater.list" containing a list of updaters that modifies the
#' parameters of the model at given time steps. If a scenario contains a
#' parameter not defined in the \code{param} object, an error will be produced.
#' See the vignette "model-parameters" for the technical detail of their
#' implementation.
#'
#' @inheritParams update_params
#' @inherit update_params return
#'
#' @export
use_scenario <- function(param, scenario) {
  scenario.params <- unique(unlist(lapply(
    scenario[[".param.updater.list"]],
    function(element) names(element[["param"]]
  ))))

  undef.params <- setdiff(scenario.params, names(param))
  if (length(undef.params) > 0) {
    stop("Some of the scenario parameters are not defined in `param`: \n'",
         paste0(undef.params, collapse = "', '"), "'")
  }

  elements.at <- vapply(
    scenario[[".param.updater.list"]],
    function(element) element[["at"]],
    numeric(1)
  )

  for (i in which(elements.at < 2)) {
    param <- update_params(
      param,
      scenario[[".param.updater.list"]][[i]][["param"]])
  }

  param[[".param.updater.list"]] <- c(
    param[[".param.updater.list"]],
    scenario[[".param.updater.list"]][elements.at >= 2]
  )

  param[[".scenario.id"]] <- scenario[["id"]]

  return(param)
}

#' Helper function validating the format of a `scenarios.df`
#' @noRd
check_scenarios_df <- function(scenarios.df) {
  checks <- c(
    all(c(".scenario.id", ".at") %in% names(scenarios.df)),
    all(as.integer(scenarios.df[[".at"]]) == scenarios.df[[".at"]])
  )

  if (! all(checks)) {
    stop(
      "A `data.frame` of scenarios must have a '.scenario.id' column \n",
      "and a '.at' column containing integers."
    )
  }
}

#' helper function to make a ragged param list into a flat one
#' @noRd
flatten_params <- function(params) {
  params <- remove_special_params(params)
  params.length <- vapply(params, length, 0)
  params.n <- sum(params.length)
  params.flat <- vector(mode = "list", length = params.n)

  i <- 1
  n <- 1
  while (n <= params.n) {
    l <- params.length[i]
    cur.param <- as.list(params[[i]])
    cur.name <- names(params)[i]
    if (l > 1)
      cur.name <- paste0(cur.name, "_", seq_len(l))
    params.flat[n:(n + l - 1)] <- cur.param
    names(params.flat)[n:(n + l - 1)] <- cur.name

    n <- n + l
    i <- i + 1
  }
  return(params.flat)
}

#' list the "special parameters" from a param list. They include some EpiModel
#' internals as well as all parameters starting with "."
#' @noRd
list_special_params <- function(params) { builtin.special.params <- c(
    "random.params",
    "random.params.values"
  )

  builtin.special.params <- intersect(builtin.special.params, names(params))
  dot.special.params <- names(params)[grep("^\\.", names(params))]

  return(unique(c(builtin.special.params, dot.special.params)))
}

#' helper function to remove the "special parameters" from a param list.
#' @noRd
remove_special_params <- function(params) {
  special.params.names <- list_special_params(params)
  params[!names(params) %in% special.params.names]
}

#' helper function to make a flat param list into a ragged one
#' @noRd
unflatten_params <- function(params.flat) {
  check_params_names(names(params.flat))
  check_params_flat(params.flat)

  set.elements.names <- vapply(
    names(params.flat),
    function(x) sub("_.*$", "", x),
    ""
  )

  params.names <- unique(set.elements.names)
  params <- vector(mode = "list", length = length(params.names))
  names(params) <- params.names

  for (i in seq_along(set.elements.names)) {
    nme <- set.elements.names[[i]]
    params[[nme]] <- c(params[[nme]], params.flat[[i]])
  }

  return(params)
}

#' helper function to check if a list of flat parameters is actually flat
#' @noRd
check_params_flat <- function(params.flat) {
  params.length <- vapply(params.flat, length, 0)
  params.list <- vapply(params.flat, is.list, TRUE)
  if (any(params.length != 1) || any(params.list)) {
    stop("A flat parameter list should contain only length 1 non list elements")
  }
  invisible(TRUE)
}

#' helper function to check the correctness of the flat parameters names
#' @noRd
check_params_names <- function(params.names) {
  params.pattern <- "^[[:alpha:]][[:alnum:].]*(_[1-9]+)?$"
  correct.format <- grepl(params.pattern, params.names)

  if (!all(correct.format)) {
    stop("The following flat parameter names are malformed: \n`",
      paste0(params.names[!correct.format], collapse = "`, `"), "`\n\n",
      "you can check the names with ",
      '`grepl("', params.pattern, '", your.names)` \n',
      "Example: 'unique.param', 'param.set_1', 'param.set_2'"
    )
  }

  special.params.names <- list_special_params(params.names)
  if (length(special.params.names) != 0) {
    stop("The following special parameter names are not allowed: \n`",
         paste0(special.params.names, collapse = "`, `"), "`\n\n"
    )
  }

  invisible(return(TRUE))
}
