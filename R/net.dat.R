#' @title Get the nodale attributes of a Master list object in network models
#'
#' @description Helper function to access the nodale attributes list of the
#'              Master list object passed by modules during the course of a
#'              network simulation
#'
#' @param attr_names the name of the attribute(s) to get (default = NULL)
#' @param list_out Should the output be a list containing the element(s)?
#'                 See **Value**. (default = FALSE)
#' @inheritParams nwupdate.net
#'
#' @return If `attr_names` is null, the full nodal attribute list. If `attr_names`
#'         is of length one and `list_out == FALSE` (default), the `attr_names`
#'         element of the nodal attribute list. If `list_out == TRUE`, a list
#'         containing the `attr_names` element(s)  of the nodale attribute list.
#'
#' @examples
#'
#' ```
#' get_attr(dat)
#' get_attr(dat, "age")
#' get_attr(dat, "age", list_out = TRUE)
#' get_attr(dat, c("age", "status"), list_out = TRUE)
#' ```
#'
#' @export
get_attr <- function(dat, attr_names = NULL, list_out = FALSE) {
  return(get_dat_elt(dat, "attr", attr_names, list_out))
}

#' @title Get the parameters of a Master list object in network models
#'
#' @description Helper function to access the parameters list of the
#'              Master list object passed by modules during the course of a
#'              network simulation
#'
#' @param param_names the name(s) of the parameter(s) to get (default = NULL)
#' @param list_out Should the output be a list containing the element(s)?
#'                 See **Value**. (default = FALSE)
#' @inheritParams nwupdate.net
#'
#' @return
#'   If `param_names` is null, the full parameters list.
#'   If `param_names` is of length one and `list_out == FALSE` (default),
#'   the `param_names` element of the nodal param_ibute list.
#'   If `list_out == TRUE`, a list containing the `param_name` element(s)
#'   of the paramets list.
#'
#' @examples
#'
#' ```
#' get_param(dat)
#' get_param(dat, "inf.prob")
#' get_param(dat, "inf.prob", list_out = TRUE)
#' get_param(dat, c("inf.prob", "act.rate"), list_out = TRUE)
#' ```
#'
#' @export
get_param <- function(dat, param_names = NULL, list_out = FALSE) {
  return(get_dat_elt(dat, "param", param_names, list_out))
}

#' @title Get the Epidemic Outputs of a Master list object in network models
#'
#' @description Helper function to access the Epidemic Outputs list of the
#'              Master list object passed by modules during the course of a
#'              network simulation
#'
#' @param epi_names the name(s) of the output(s) to get (default = NULL)
#' @param list_out Should the output be a list containing the element(s)?
#'                 See **Value**. (default = FALSE)
#' @inheritParams nwupdate.net
#'
#' @return
#'   If `epi_names` is null, the full Epidemic Outputs list.
#'   If `epi_names` is of length one and `list_out == FALSE` (default),
#'   the `epi_names` elements of the Output list.
#'   If `list_out == TRUE`, a list containing the `out_name` element(s)
#'   of the parameters list.
#'
#' @examples
#'
#' ```
#' get_epi(dat)
#' get_epi(dat, "i")
#' get_epi(dat, "i", list_out = TRUE)
#' get_epi(dat, c("i", "s"), list_out = TRUE)
#' ```
#'
#' @export
get_epi <- function(dat, epi_names = NULL, list_out = FALSE) {
  return(get_dat_elt(dat, "epi", epi_names, list_out))
}

#' @title internal helper function to get elements of a Master list object
#'
#' @param elt the name of the sublist to access
#' @param attrnames the name(s) of the `elt` list elements to get
#' @param list_out Should the output be a list containing the element(s)?
#'                 See **Value**.
#' @inheritParams nwupdate.net
#'
#' @return
#'   If `attrname` is null, the full `elt` list.
#'   If `attrname` is of length one and `list_out == FALSE` (default),
#'   the `attrname` elements of the `elt` list.
#'   If `list_out == TRUE`, a list containing the `attrname` element(s)
#'   of the `elt` list.
#'
#' @examples
#'
#' ```
#' get_dat_elt(dat, "attr")
#' get_dat_elt(dat, "attr", "age")
#' get_dat_elt(dat, "attr", "age", list_out = TRUE)
#' get_dat_elt(dat, "attr", c("age", "status"), list_out = TRUE)
#' ```
#'
#' @keywords Internal
get_dat_elt <- function(dat, elt, attrnames, list_out) {
 if (is.null(attrnames)) {
    return(dat[[elt]])
  } else if (list_out) {
    return(dat[[elt]][attrnames])
  } else {
    return(dat[[elt]][[attrnames]])
  }
}

#' @title Edit an element of the nodale attributes of a Master list object
#'
#' @description Helper function to modify the nodale attribute list of the
#'              Master list object passed by modules during the course of a
#'              network simulation
#'
#' @param attr_name the name of the attribute to modify
#' @param value the new value for `attr_name`
#' @inheritParams nwupdate.net
#'
#' @return The modified Master list object
#'
#' @examples
#'
#' ```
#' dat <- set_attr(dat, "age", new_ages)
#' ```
#'
#' @keywords Internal
set_attr <- function(dat, attr_name, value) {
  return(set_dat_elt(dat, "attr", attr_name, value))
}

#' @title Edit an element of the parameters of a Master list object
#'
#' @description Helper function to modify the parameters list of the
#'              Master list object passed by modules during the course of a
#'              network simulation
#'
#' @param param_name the name of the parameter to modify
#' @param value the new value for `param_name`
#' @inheritParams nwupdate.net
#'
#' @return The modified Master list object
#'
#' @examples
#'
#' ```
#' dat <- set_param(dat, "act.rate", new_rate)
#' ```
#'
#' @keywords Internal
set_param <- function(dat, param_name, value) {
  return(set_dat_elt(dat, "param", param_name, value))
}

#' @title Edit an Epidemic Outcome of the parameters of a Master list object
#'
#' @description Helper function to modify the Epidemic Outcome list of the
#'              Master list object passed by modules during the course of a
#'              network simulation
#'
#' @param epi_name the name of the outcome to modify
#' @param value the new value for `epi_name`
#' @inheritParams nwupdate.net
#'
#' @return The modified Master list object
#'
#' @examples
#'
#' ```
#' dat <- set_epi(dat, "i", new_rate)
#' ```
#'
#' @keywords Internal
set_epi <- function(dat, epi_name, value) {
  return(set_dat_elt(dat, "epi", epi_name, value))
}

#' @title internal helper function to edit elements of a Master list object
#'
#' @param elt the name of the sublist to access
#' @param attrname the name of the output(s) to get (default = NULL)
#' @inheritParams nwupdate.net
#'
#' @return The modified Master list object
#'
#' @examples
#'
#' ```
#' dat <- set_dat_elt(dat, "attr", "age", new_ages)
#' ```
#'
#' @keywords Internal
set_dat_elt <- function(dat, elt, attrname, value) {
  dat[[elt]][[attrname]] <- value

  return(dat)
}
