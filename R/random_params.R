#' Create a value sampler for random parameters definition
#'
#' This function return a 0 argument function that can be used as a generator
#' function in the `random_params` argument of the `param.net` function.
#'
#' @param values a vector of values to sample from
#' @param prob a vector of weights to use during sampling, if `NULL`,
#' all values have the same probability of being picked (default = `NULL`)
#' @return one of the values from the `values` vector
#'
#' @export
param_random <- function(values, prob = NULL) {
  if (!is.null(prob) && length(prob) != length(values)) {
    stop("incorrect number of probabilites")
  }

  f <- function() {
    return(sample(x = values, size = 1, prob = weights, replace = TRUE))
  }

  return(f)
}

#' Generate values for the random parameters
#'
#' This function uses the generative function in `random.params` to create
#' values for the parameters
#'
#' @param param the `param` argument received by the `netsim` functions
#' @return a fully instanciated `param` list
#'
#' @section `random_params`:
#' The `random_params` argument to the `param.net` function must be a named list
#' of functions that return a values that can be used as the argument with the
#' same name. In the example below, `param_random` is a function factory
#' provided by EpiModel for `acte.rate` and `tx.halt.part.prob` we provide
#' bespoke functions
#'
#' @section Generator Functions:
#' The function used inside `random_parameters` must be 0 argument functions
#' returning a valid value for the parameter with the same name.
#'
#' @examples
#' \dontrun{
#' param <- param.net(
#'   inf.prob = 0.3,
#'   random_params = list(
#'     act.rate = param_random(c(0.25, 0.5, 0.75), prob = c(0.1, 0.2, 0.7)),
#'     tx.halt.part.prob = function() rbeta(1, 1, 2)
#'   )
#' )
#'
#' param <- generate_random_params(param)
#' }
#'
generate_random_params <- function(param) {
  if (is.null(param$random.params) || length(param$random.params) == 0) {
    return(param)
  }

  if (!is.list(param$random.params)) {
    stop("`random.params` must be named list of functions")
  }

  rng_names <- names(param$random.params)
  if (any(rng_names == "")) {
    stop("all elements of `random.params` must be named")
  }

  if (! all(vapply(param$random.params, is.function, TRUE))) {
    stop("all elements of `random.params` must be functions")
  }

  param[rng_names] <- lapply(param$random.params, do.call, args = list())

  return(param)
}
