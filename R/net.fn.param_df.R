#' Get the `param.df` formated type of a value
#' @noRd
param.df_type <- function(x) {
  if (is.logical(x)) {
    "logical"
  } else if (is.numeric(x)) {
    "numeric"
  } else {
    "character"
  }
}

#' Coerce a value to a `param.df` type
#' @noRd
as_param.df_type <- function(x, type) {
  if (type == "logical") {
    as.logical(x)
  } else if (type == "numeric") {
    as.numeric(x)
  } else {
    as.character(x)
  }
}

#' List the recognized `param.df` types
#' @noRd
list_param.df_types <- function() {
  c("numeric", "logical", "character")
}

#' @title Parameters List for Stochastic Network Models from a Formatted
#'        Data Frame
#'
#' @description Sets the epidemic parameters for stochastic network models with
#'              [netsim()] using a specially formatted data frame of
#'              parameters.
#'
#' @param long.param.df A `data.frame` of parameters. See details for the
#'                      expected format.
#'
#' @return A list object of class `param.net`, which can be passed to
#'         [netsim()].
#'
#' @section long.param.df:
#' It is possible to set input parameters using a specifically formatted
#' `data.frame` object. The first 3 columns of this `data.frame` must
#' be:
#'
#'  * `param`: The name of the parameter. If this is a non-scalar
#'    parameter (a vector of length > 1), end the parameter name with the
#'    position on the vector (e.g., `"p_1"`, `"p_2"`, ...).
#'  * `value`: the value for the parameter (or the value of the
#'    parameter in the Nth position if non-scalar).
#'  * `type`: a character string containing either `"numeric"`,
#'    `"logical"`, or `"character"` to define the parameter object
#'    class.
#'
#'
#' In addition to these 3 columns, the `data.frame` can contain any number
#' of other columns, such as `details` or `source` columns to document
#' parameter meta-data. However, these extra columns will not be used by
#' EpiModel.
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
  if (!all(long.param.df$type %in% list_param.df_types())) {
    stop("The `type` column must contain only: ",
         paste0(list_param.df_types(), collapse = ", "))
  }
  check_params_names(long.param.df$param)

  duplicated_params <- duplicated(long.param.df$param)
  duplicated_params <- long.param.df$param[duplicated_params]
  if (length(duplicated_params) > 0) {
    stop("The following parameters are duplicated: `",
         paste0(duplicated_params, collapse = "`, `"), "`")
  }

  # To flat params
  flat.params <- Map(as_param.df_type, long.param.df$value, long.param.df$type)
  names(flat.params) <- long.param.df$param

  # To param.list
  param <- unflatten_params(flat.params)
  class(param) <- c("param.net", "list")

  return(param)
}


#' Coerce a list of parameters to a `long.param.df`
#'
#' @param params A list of parameters to be formatted into a `long.param.df`
#'
#' @inheritSection param.net_from_table long.param.df
#'
#' @details
#' Parameters given per network layer with [multilayer()] or by duration of
#' infection with [by_infection_duration()] are not supported and stop with
#' an error. A parameter table records each element of a vector as its own row
#' but not what the vector means, so such a value would read back as a plain
#' vector, which the built-in modules do not accept. Remove these parameters
#' before converting and set them again after reading the table back.
#'
#' @return  A `data.frame` of parameters.
#'
#' @export
param.net_to_table <- function(params) {
  # a multilayer or by_infection_duration value would flatten into a plain
  # vector, which loses what its positions mean
  structured <- vapply(params, inherits, logical(1),
                       c("multilayer", "by_infection_duration"))
  if (any(structured)) {
    stop("param.net_to_table() does not support multilayer() or ",
         "by_infection_duration() parameters: `",
         paste(names(params)[structured], collapse = "`, `"), "`. A ",
         "parameter table cannot record what the positions of a vector ",
         "mean, so these would read back as plain vectors. Remove them before ",
         "converting and set them again after reading the table back.",
         call. = FALSE)
  }
  flat_params <- flatten_params(params)
  check_params_names(names(flat_params))
  dplyr::tibble(
    param = names(flat_params),
    value = vapply(flat_params, as.character, ""),
    type = vapply(flat_params, param.df_type, "")
  )
}