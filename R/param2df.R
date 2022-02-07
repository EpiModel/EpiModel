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
make_scenarios_list <- function(scenarios.df) {
  check_scenarios_df(scenarios.df)

  scenarios <- lapply(
    unique(scenario.df[[".scenario.id"]]),
    function(id) make_scenario(dplyr::filter(scenario.df, .scenario.id == id)),
  )

  return(scenarios)
}

#' a scenario is a list with the following elements:
#'   - `id`: the identifier of the scenario
#'   - `param.updater.list`: as in the updater module
#'
#' when running "use_scenario", if any updater has an `at` value < 2, it is
#' applied immediatly (so before netsim)
#'
#' @noRd
make_scenario <- function(scenario.rows) {
  scenario <- list(
    id = scenario.rows[1, ".scenario.id"],
    param.updater.list = vector(mode = "list", length = nrow(scenario.rows))
  )

  elements.at <- scenario.rows[[".at"]]
  # ensures that `scenario.rows` is a `tibble`
  scenario.rows <- dplyr::select(scenario.rows, - c(".at", ".scenario.id"))

  for (i in seq_along(elements.at)) {
    scenario[["param.updater.list"]][[i]] <- list(
      at = elements.at[[i]],
      params = unflatten_params(scenario.rows[i, ])
    )
  }

  return(scenario)
}

#' Apply a scenario object to a param.net object
#'
#' @section scenario:
#' Can be made by make_scenarios_list from a scenarios.df
use_scenario <- function(param, scenario) {
  elements.at <- vapply(
    scenario[["param.updater.list"]],
    function(element) element[["at"]],
    numeric(1)
  )

  for (i in which(elements.at < 2)) {
    param <- update_params(
      param,
      scenario[["param.updater.list"]][[i]][["params"]])
  }

  param[["param.updater.list"]] <- c(
    param[["param.updater.list"]],
    scenario[["param.updater.list"]][elements.at >= 2]
  )

  param[[".scenario.id"]] <- scenario[["id"]]

  # print a message describing the scenario

  return(param)
}

#' Helper function validating the format of a `scenarios.df`
#' @noRd
check_scenarios_df <- function(scenarios.df) {
  checks <- c(
    all(c(".scenario.name", ".at") %in% names(scenarios.df)),
    all(as.integer(scenarios.df[[".at"]]) == scenarios.df[[".at"]])
  )

  if (! all(checks)) {
    stop(
      "A `data.frame` of scenarios must have a '.scenario.id' column",
      "and an '.at' column containing integers."
    )
  }
}

#' Make a list of scenarios from a data.frame of scenarios
#'
#' An EpiModel scenario is a list of parameters to be used by a model. They are
#' usually used by a researcher who wants to model counterfactuals using a pre
#' calibrated model.
#'
#' @param scenarios.df a \code{data.frame}
#'
#' @return a list of parameter list
#'
#' @section scenarios.df:
#' The \code{scenarios.df} is a \code{data.frame} of values to be used as
#' parameters.
#'
#' The column names must correspond either to:
#' the name of one parameter if this parameter is of size 1 or the name of the
#' parameter with "_1", "_2", "_N" with the second part being the position of
#' the value for a parameter of size > 1. This means that the parameter names
#' cannot contain any underscore "_". (e.g "a.rate", "d.rate_1", "d.rate_2")
#'
#' @section Scenario Names:
#' If \code{scenarios.df} contains a column ".scenario.id", it will be used to
#' set the names of the scenarios. Otherwise, they will be named sequentially.
#'
#' @section When Should a Scenario Apply:
#' If \code{scenarios.df} contains a column ".at", it will be used to define
#' when the parameters defined in the scenario should be used. If the column is
#' not present or if ".at" is less than 2, the parameters will be changed before
#' the initialization of the model. Otherwise they will occurs at the beginning
#' of the specified timestep.
#'
#' @export
make_scenarios_list <- function(scenarios.df) {
  if (!is.null(scenarios.df[[".scenario.id"]])) {
    scenarios.names <- scenarios.df[[".scenario.id"]]
    scenarios.df[[".scenario.id"]] <- NULL
  } else {
    scenarios.names <- seq_len(nrow(scenarios.df))
  }

  scenarios.list <- lapply(
    seq_len(nrow(scenarios.df)),
    function(row) unflatten_params(scenarios.df[row, ])
  )
  names(scenarios.list) <- scenarios.names

  return(scenarios.list)
}

#' Make a list of scenarios updaters from a data.frame of scenarios
#'
#' An EpiModel scenario updater is a list of parameters to be used by a model.
#' They are usually used by a researcher who wants to model counterfactuals
#' using a pre calibrated model.
#'
#' @param scenarios.df a \code{data.frame}
#'
#' @return a list of parameter list
#'
#' @section scenarios.df:
#' The \code{scenarios.df} is a \code{data.frame} of values to be used as
#' parameters.
#'
#' The column names must correspond either to:
#' the name of one parameter if this parameter is of size 1 or the name of the
#' parameter with "_1", "_2", "_N" with the second part being the position of
#' the value for a parameter of size > 1. This means that the parameter names
#' cannot contain any underscore "_". (e.g "a.rate", "d.rate_1", "d.rate_2")
#'
#' @section :
#' If \code{scenarios.df} contains a column ".scenario.id", it will be used to
#' set the names of the scenarios. Otherwise, they will be named sequentially.
#'
#' @export
make_scenarios_updaters <- function(scenarios.df, at) {
  scenarios.list <- make_scenarios_list(scenarios.df)
  updaters.list <- lapply(
    scenarios.list,
    function(scenario) list(list(at = at, param = scenario))
  )
  names(updaters.list) <- names(scenarios.list)
  return(updaters.list)
}

make_scenarios_df <- function(scenarios.list) {
  scenarios.rows <- lapply(scenarios.list, flatten_params)
  scenarios.df <- dplyr::bind_rows(scenarios.rows)
  scenarios.df <- dplyr::mutate(
    scenarios.df,
    .scenario.id = names(scenarios.list)
  )
  dplyr::select(scenarios.df, .scenario.id, dplyr::everything())
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
list_special_params <- function(params) {
  builtin.special.params <- c(
    "param.updater.list",
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

# update params to create scenarios
#
# params_orig <- function() {
#   list(
#     aasd.ad = sample(seq_len(9), 1),
#     d.adk = sample(seq_len(9), 3, replace = TRUE),
#     eas.d = sample(letters, 1, replace = TRUE),
#     fd.dd = sample(LETTERS, 4, replace = TRUE)
#   )
# }
#
# scenarios.df <- dplyr::bind_rows(
#   lapply(1:5, function(x) tibble::as_tibble(flatten_params(params_orig())))
# )
# scenarios.df[[".scenario.id"]] <- replicate(
#   5,
#   paste0(sample(LETTERS, 10, replace = TRUE), collapse = "")
# )
#
# scenarios.list <- make_scenarios_list(scenarios.df)
# scenarios.updaters <- make_scenarios_updaters(scenarios.df, 100)
#
# param_df <- read.csv("/home/adrien/Documents/Projects/BigNets/data/input/calibration.csv")
# make_scenarios_list(param_df)
