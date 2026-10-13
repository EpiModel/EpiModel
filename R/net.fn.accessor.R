#' @title Functions to Access and Edit the Main netsim_dat Object in Network Models
#'
#' @description These `get_`, `set_`, `append_`, and `add_`
#'              functions allow a safe and efficient way to retrieve and mutate
#'              the main `netsim_dat` class object of network models
#'              (typical variable name `dat`). They are intended for use
#'              *inside module functions* that run during a [netsim()]
#'              simulation, not for editing the user-facing `param.net`,
#'              `init.net`, or `control.net` inputs before a run or the
#'              output object returned by [netsim()]. See the **Intended
#'              Usage** section below for details and alternatives.
#'
#' The accessors are documented on three pages, by the part of `dat` they act
#' on:
#'
#'  * **Nodal attributes**, with one value per node ([get_attr()]):
#'    `get_attr()`, `get_attr_list()`, `set_attr()`, `add_attr()`,
#'    `append_attr()`, `append_core_attr()`, `remove_node_attr()`, and
#'    `get_core_attributes()`.
#'  * **Epidemic outputs**, with one value per time step ([get_epi()]):
#'    `get_epi()`, `get_epi_list()`, `set_epi()`, and `add_epi()`.
#'  * **Parameters, control settings, and initial conditions**
#'    ([get_param()]): `get_param()`, `get_param_list()`, `set_param()`,
#'    `add_param()`, the same four functions for `control` and `init`, and
#'    `get_network_control()`.
#'
#' Nodes are referred to by their positional identifiers; see
#' [get_unique_ids()] and [get_posit_ids()] to convert them to and from
#' their unique identifiers.
#'
#' @section Intended Usage:
#' These accessors operate on the live `netsim_dat` object that EpiModel
#' constructs internally when a simulation starts and passes through each
#' module at every time step. The expected calling context is **inside a
#' user-supplied module function** (e.g., a custom `infection.FUN`,
#' `recovery.FUN`, or arrivals/departures module) registered with
#' [control.net()].
#'
#' They are *not* a general-purpose tool for editing the user-facing input
#' or output objects of [netsim()]:
#'
#'  * Calling them on a [param.net()], [init.net()], or [control.net()]
#'    object before a run will appear to succeed (those objects are also
#'    list-like) but produces an object that is not a valid simulation
#'    input.
#'  * Calling them on the object returned by [netsim()] modifies only the
#'    saved input slot, not anything that would change the result of a
#'    re-run. Note also that [netsim()] strips the `*.FUN` entries from
#'    its returned `control` slot, so feeding that slot back into a new
#'    simulation will fail.
#'
#' For modifications outside a running simulation, use the following
#' instead:
#'
#'  * **Editing parameters before a run:** [update_params()] for a
#'    `param.net` object, or direct list assignment
#'    (e.g., `p$inf.prob <- 0.5`).
#'  * **Editing init or control before a run:** direct list assignment
#'    (e.g., `ctrl$nsteps <- 1000`), or rebuild with a fresh call to the
#'    [init.net()] / [control.net()] constructor.
#'  * **Scheduled mid-simulation changes:** the `.param.updater.list`
#'    argument to [param.net()] and the `.control.updater.list` argument
#'    to [control.net()]. See the "model-parameters" vignette for the
#'    underlying updater module and the scenario API built on top of it.
#'
#' @section Mutability:
#' The `set_`, `append_`, and `add_` functions DO NOT modify the
#' `netsim_dat` object in place. The result must be assigned back to
#' `dat` in order to be registered: `dat <- set_*(dat, item, value)`.
#'
#' @section `set_` and `append_` vs `add_`:
#' The `set_` and `append_` functions edit a pre-existing element or
#' create a new one if it does not exist already by calling the `add_`
#' functions internally.
#'
#' @seealso [update_params()] for editing a `param.net` object outside a
#'   simulation; [param.net()], [init.net()], [control.net()] for input
#'   constructors; [netsim()] for running a simulation.
#'
#' @name net-accessor
NULL

raw_get_attr_list <- function(dat) return(dat$run$attr)
raw_set_attr_list <- function(dat, attr_list) {
  dat$run$attr <- attr_list
  return(dat)
}
raw_set_attr <- function(dat, item, val) {
  dat$run$attr[[item]] <- val
  return(dat)
}

# Nodal attributes -------------------------------------------------------------

#' @rdname get_attr
#' @order 2
#' @export
get_attr_list <- function(dat, item = NULL) {
  attr_list <- raw_get_attr_list(dat)
  if (is.null(item)) {
    return(attr_list)
  }

  missing_item <- setdiff(item, names(attr_list))
  if (length(missing_item) > 0L) {
    stop("There are no attributes called `",
         paste(missing_item, collapse = ", "),
         "` in the attributes list of the main list object (dat)")
  }

  out <- attr_list[item]

  return(out)
}

#' @title Get and Set Nodal Attributes in Network Models
#'
#' @description Retrieve, create, and edit the nodal attributes of the main
#'              `netsim_dat` object of network models, from inside module
#'              functions. Each nodal attribute is a vector with one value per
#'              node, indexed by the positional identifiers of the nodes
#'              (`posit_ids`).
#'
#'  * `get_attr()` returns one attribute, for all nodes or for the nodes in
#'    `posit_ids`. `get_attr_list()` returns several attributes as a named
#'    list.
#'  * `set_attr()` replaces the values of an attribute, for all nodes or for
#'    the nodes in `posit_ids`, and creates the attribute if it does not
#'    exist. `add_attr()` creates an attribute set to `NA` for all nodes.
#'  * `append_core_attr()` adds `n.new` nodes and sets their core
#'    attributes. `append_attr()` then adds the values of each other
#'    attribute for the same nodes.
#'  * `remove_node_attr()` removes the nodes in `posit_ids` from every
#'    attribute.
#'  * `get_core_attributes()` lists the core attributes with their types.
#'
#' See [net-accessor] for an overview of the accessors.
#'
#' @inheritParams recovery.net
#' @param item Name of the nodal attribute. For `get_attr_list()`, a
#'        character vector of names, or `NULL` (default) for all attributes.
#' @param posit_ids Positional identifiers of the nodes. For `get_attr()` and
#'        `set_attr()`, the nodes to get or edit, with `NULL` (default) for all
#'        nodes; for `remove_node_attr()`, the nodes to remove.
#' @param value New values of the attribute. For `set_attr()`, one value per
#'        node, or with `posit_ids`, one value per selected node or a single
#'        value for all of them. For `append_attr()`, one value per new node
#'        or a single value for all of them.
#' @param override.null.error If `TRUE`, return `NULL` instead of an error
#'        when `item` does not exist (default = `FALSE`).
#' @param override.length.check If `TRUE`, `set_attr()` accepts a `value`
#'        whose length differs from the number of nodes (default = `FALSE`).
#'        Every attribute must again hold one value per node by the end of
#'        the time step.
#' @param n.new Number of new nodes: for `append_core_attr()`, the nodes to
#'        add; for `append_attr()`, the number of values to append at the end
#'        of `item`.
#' @param at Current time step, recorded as the entry time (`entrTime`) of
#'        the new nodes.
#'
#' @return For `get_attr()`, the attribute vector, or its values for the
#'         nodes in `posit_ids`; for `get_attr_list()`, a named list of
#'         attribute vectors; for `get_core_attributes()`, a named list with a
#'         value of the correct type for each core attribute. The other
#'         functions return the updated `netsim_dat` object.
#'
#' @section Core Attributes:
#' The `append_core_attr` function initializes the attributes necessary for
#' EpiModel to work (the four core attributes are: "active", "unique_id",
#' "entrTime", and "exitTime"). These attributes are used in the initialization
#' phase of the simulation, to create the nodes (see
#' [initialize.net()]); and also used when adding nodes during the
#' simulation (see [arrivals.net()]).
#'
#' @section One Value Per Node:
#' Every nodal attribute holds exactly one value per node, and [netsim()]
#' checks this at the end of each time step. When `set_attr` and `append_attr`
#' create an attribute, the nodes they do not write to are filled with `NA` to
#' preserve that: `set_attr` with `posit_ids` fills the unselected nodes, and
#' `append_attr` fills the nodes that predate the ones being appended. This
#' means `append_core_attr` must be called before `append_attr` when adding new
#' nodes, as it is what makes those nodes known to the rest of the object.
#'
#' @inheritSection net-accessor Intended Usage
#' @inheritSection net-accessor Mutability
#'
#' @seealso [net-accessor] for an overview of the accessors, and [get_epi()]
#'   and [get_param()] for the other groups; [arrive_nodes()] and
#'   [depart_nodes()] for adding and removing nodes from the networks;
#'   [get_unique_ids()] and [get_posit_ids()] for converting node identifiers.
#'
#' @examples
#' dat <- create_dat_object(control = list(nsteps = 150))
#' dat <- append_core_attr(dat, at = 1, n.new = 10)
#'
#' dat <- add_attr(dat, "age")
#' dat <- set_attr(dat, "age", runif(10, 15, 65))
#' dat <- set_attr(dat, "status", rep(c("s", "i"), c(8, 2)))
#' get_attr_list(dat)
#' get_attr_list(dat, c("age", "active"))
#' get_attr(dat, "status")
#' get_attr(dat, "status", posit_ids = c(1, 9))
#'
#' # Edit selected nodes
#' dat <- set_attr(dat, "status", "i", posit_ids = 3)
#' get_attr(dat, "status")
#'
#' # Add 5 nodes: core attributes first, then the others
#' dat <- append_core_attr(dat, at = 2, n.new = 5)
#' dat <- append_attr(dat, "status", "s", 5)
#' dat <- append_attr(dat, "age", 15, 5)
#' get_attr(dat, "age")
#' get_attr(dat, "entrTime")
#'
#' @order 1
#' @export
get_attr <- function(dat, item, posit_ids = NULL, override.null.error = FALSE) {
  attr_list <- raw_get_attr_list(dat)

  if (!item %in% names(attr_list)) {
    if (!override.null.error) {
      stop("There is no attribute called `", item,
           "` in the attributes list of the main list object (dat)")
    }
    return(NULL)
  }

  if (is.null(posit_ids)) {
    return(attr_list[[item]])
  }

  assert_valid_posit_ids(dat, posit_ids)
  attr <- attr_list[[item]]
  return(attr[posit_ids])
}

#' @rdname get_attr
#' @order 4
#' @export
add_attr <- function(dat, item) {
  attr_list <- get_attr_list(dat)
  if (item %in% names(attr_list)) {
    stop("Cannot create the attribute '", item, "': exists already")
  }
  dat <- raw_set_attr(dat, item, rep(NA, length(attr_list$active)))
  return(dat)
}

#' @rdname get_attr
#' @order 3
#' @export
set_attr <- function(dat, item, value, posit_ids = NULL,
                     override.length.check = FALSE) {
  attr_list <- raw_get_attr_list(dat)
  if (!item %in% names(attr_list)) {
    dat <- add_attr(dat, item)
    # re-read, so that the `posit_ids` branch below writes into the full-length
    # NA vector just created rather than into the pre-creation snapshot
    attr_list <- raw_get_attr_list(dat)
  }

  if (is.null(posit_ids)) {
    if (!override.length.check &&
          length(value) != length(attr_list$active)) {
      stop(
        "When trying to edit the ", `item`, " nodal attribute: ",
        "The size of the `value` vector is not equal to the number of nodes in",
        " the network. \n",
        "Expected: ", length(attr_list$active), "\n",
        "Given: ", length(value)
      )
    }
    dat <- raw_set_attr(dat, item, value)
  } else {
    assert_valid_posit_ids(dat, posit_ids)
    if (length(value) != 1L && length(value) != length(posit_ids)) {
      stop(
        "When trying to edit the `", item, "` nodal attribute: ",
        "The size of the `value` vector is not equal to the number of nodes ",
        "selected by the `posit_ids` vector nor of length 1. \n",
        "Expected: ", length(attr_list$active[posit_ids]), " or 1 \n",
        "Given: ", length(value)
      )
    }
    attr <- attr_list[[item]]
    attr[posit_ids] <- value
    dat <- raw_set_attr(dat, item, attr)
  }

  return(dat)
}

#' @rdname get_attr
#' @order 5
#' @export
append_attr <- function(dat, item, value, n.new) {
  if (!is.numeric(n.new) || n.new < 0L) {
    stop("`n.new` must be numeric and greater than or equal to zero.")
  }

  if (length(value) == 1L) {
    new_vals <- rep(value, n.new)
  } else if (length(value) == n.new) {
    new_vals <- value
  } else {
    stop("`value` must be of length one or `n.new`.")
  }

  old_vals <- get_attr(dat, item, override.null.error = TRUE)
  if (is.null(old_vals) && item != "active") {
    # the attribute does not exist yet: the nodes that predate the ones being
    # appended get NA, so that the result is still one value per node.
    # `append_core_attr` appends to `active` first, so `active` already covers
    # the new nodes here. `active` itself is exempt, as it is what defines the
    # number of nodes and creating it is how the first nodes come to exist
    n_active <- length(get_posit_ids(dat))
    old_vals <- rep(NA, max(n_active - n.new, 0L))
  }
  dat <- set_attr(dat, item, c(old_vals, new_vals),
                  override.length.check = TRUE)

  return(dat)
}

#' @rdname get_attr
#' @order 7
#' @export
remove_node_attr <- function(dat, posit_ids) {
  if (is.null(posit_ids)) {
    return(dat)
  }
  assert_valid_posit_ids(dat, posit_ids)
  attr_list <- raw_get_attr_list(dat)
  attr_list <- lapply(attr_list, function(x) x[-posit_ids])
  dat <- raw_set_attr_list(dat, attr_list)
  return(dat)
}

# Epidemic outputs -------------------------------------------------------------

#' @rdname get_epi
#' @order 2
#' @export
get_epi_list <- function(dat, item = NULL) {
  if (is.null(item)) {
    out <- dat$epi
  } else {
    missing_item <- setdiff(item, names(dat$epi))
    if (length(missing_item) > 0L) {
      stop("There is no epi output called `",
           paste(missing_item, collapse = ", "),
           "` in the epi output list of the main list object (dat)")
    }
    out <- dat$epi[item]
  }
  return(out)
}

#' @title Get and Set Epidemic Outputs in Network Models
#'
#' @description Retrieve, create, and edit the epidemic outputs of the main
#'              `netsim_dat` object of network models, from inside module
#'              functions. Each epidemic output is a vector with one value per
#'              time step, which becomes a column of the data frame returned by
#'              [as.data.frame.netsim()] after the simulation.
#'
#'  * `get_epi()` returns one output, for all time steps or for the time
#'    steps in `at`. `get_epi_list()` returns several outputs as a named
#'    list.
#'  * `set_epi()` records `value` at one time step, the current one by
#'    default, and creates the output if it does not exist. `add_epi()`
#'    creates an output set to `NA` at every time step.
#'
#' Unlike the other `set_` functions, `set_epi()` takes a time step `at`,
#' which comes last and must be passed by name.
#'
#' See [net-accessor] for an overview of the accessors.
#'
#' @inheritParams get_attr
#' @param item Name of the epidemic output. For `get_epi_list()`, a
#'        character vector of names, or `NULL` (default) for all outputs.
#' @param value Value of the output at time step `at`.
#' @param at For `get_epi()`, the time steps at which to return the output,
#'        with `NULL` (default) for all time steps. For `set_epi()`, the time
#'        step at which to record `value`, which defaults to the current time
#'        step ([get_current_timestep()]) when `NULL`.
#' @param ... Must be empty. It makes `at` a named-only argument of
#'        `set_epi()`, so that calls written for the `set_epi(dat, item, at,
#'        value)` signature of EpiModel 2.6 and earlier error instead of
#'        silently swapping `at` and `value`.
#'
#' @return For `get_epi()`, the output vector, or its values at the time
#'         steps in `at`; for `get_epi_list()`, a named list of output
#'         vectors. The other functions return the updated `netsim_dat`
#'         object.
#'
#' @inheritSection net-accessor Intended Usage
#' @inheritSection net-accessor Mutability
#'
#' @seealso [net-accessor] for an overview of the accessors, and [get_attr()]
#'   and [get_param()] for the other groups; [epi_trackers()] for outputs
#'   computed by tracker functions at each time step; [mutate_epi()] for
#'   outputs derived after the simulation.
#'
#' @examples
#' dat <- create_dat_object(control = list(nsteps = 5))
#' dat <- set_current_timestep(dat, 1)
#'
#' dat <- add_epi(dat, "i.num")
#' dat <- set_epi(dat, "i.num", 10)
#' dat <- set_epi(dat, "s.num", 90)
#' get_epi_list(dat)
#'
#' # Record at a time step other than the current one
#' dat <- set_epi(dat, "i.num", 12, at = 2)
#' get_epi_list(dat, c("i.num", "s.num"))
#' get_epi(dat, "i.num")
#' get_epi(dat, "i.num", at = c(1, 2))
#'
#' @order 1
#' @export
get_epi <- function(dat, item, at = NULL, override.null.error = FALSE) {
  if (!item %in% names(dat$epi)) {
    if (!override.null.error) {
      stop("There is no epi output called `", item,
           "` in the epi output list of the main list object (dat)")
    }
    return(NULL)
  }

  if (is.null(at)) {
    return(dat$epi[[item]])
  }

  assert_valid_time_steps(dat, at)
  return(dat$epi[[item]][at])
}

#' @rdname get_epi
#' @order 4
#' @export
add_epi <- function(dat, item) {
  if (item %in% names(dat$epi)) {
    stop("Cannot create the epi output, ", item, ": exists already")
  }
  dat$epi[[item]] <- padded_vector(NA_integer_, get_control(dat, "nsteps"))
  return(dat)
}

#' @rdname get_epi
#' @order 3
#' @export
set_epi <- function(dat, item, value, ..., at = NULL) {
  if (...length() > 0L) {
    stop_old_signature(
      "set_epi",
      new_usage = "set_epi(dat, item, value, at = NULL)",
      old_usage = "set_epi(dat, item, at, value)",
      fix = "set_epi(dat, item, value)",
      dots_names = ...names()
    )
  }
  check_record_name(item, "item")
  at <- resolve_record_at(dat, at)

  if (!item %in% names(dat$epi)) {
    dat <- add_epi(dat, item)
  }

  # ensure that the vector is of size `nsteps`, right padded with NA
  if (at > length(dat$epi[[item]])) {
    nsteps <- get_control(dat, "nsteps")
    dat$epi[[item]] <- padded_vector(dat$epi[[item]], nsteps)
  }

  dat$epi[[item]][at] <- value
  return(dat)
}

# Parameters, controls, and initial conditions ---------------------------------

#' @rdname get_param
#' @order 2
#' @export
get_param_list <- function(dat, item = NULL) {
  if (is.null(item)) {
    return(dat$param)
  }

  missing_item <- setdiff(item, names(dat$param))
  if (length(missing_item) > 0L) {
    stop("There are no parameters called `",
         paste(missing_item, collapse = ", "),
         "` in the parameter list of the main list object (dat)")
  }
  return(dat$param[item])
}

#' @title Get and Set Parameters, Controls, and Initial Conditions in Network
#'        Models
#'
#' @description Retrieve, create, and edit the parameters, control settings,
#'              and initial conditions of the main `netsim_dat` object of
#'              network models, from inside module functions. They start as
#'              the [param.net()], [control.net()], and [init.net()] inputs of
#'              [netsim()], and the three share one set of accessors:
#'
#'  * `get_param()`, `get_control()`, and `get_init()` return one element.
#'    `get_param_list()`, `get_control_list()`, and `get_init_list()` return
#'    several elements as a named list.
#'  * `set_param()`, `set_control()`, and `set_init()` replace the value of an
#'    element, and create the element if it does not exist. `add_param()`,
#'    `add_control()`, and `add_init()` create an element set to `NA`.
#'  * `get_network_control()` returns the entry for one network layer of a
#'    control setting given as a [multilayer()] object.
#'
#' Unlike [set_attr()] and [set_epi()], these `set_` functions take neither
#' nodes nor a time step: the new value is used from then on, until it is set
#' again.
#'
#' See [net-accessor] for an overview of the accessors.
#'
#' @inheritParams get_attr
#' @param item Name of the parameter, control setting, or initial condition.
#'        For the `get_*_list()` functions, a character vector of names, or
#'        `NULL` (default) for all elements.
#' @param value New value of the element.
#' @param network Index of the network layer.
#'
#' @return For `get_param()`, `get_control()`, and `get_init()`, the value of
#'         the element; for `get_network_control()`, its entry for layer
#'         `network`; for the `get_*_list()` functions, a named list of
#'         elements. The other functions return the updated `netsim_dat`
#'         object.
#'
#' @inheritSection net-accessor Intended Usage
#' @inheritSection net-accessor Mutability
#'
#' @seealso [net-accessor] for an overview of the accessors, and [get_attr()]
#'   and [get_epi()] for the other groups; [update_params()] for editing a
#'   `param.net` object outside a simulation.
#'
#' @examples
#' dat <- create_dat_object(control = list(nsteps = 150))
#'
#' dat <- add_param(dat, "x")
#' get_param(dat, "x")
#' dat <- set_param(dat, "x", 0.4)
#' dat <- set_param(dat, "y", 0.8)
#' get_param_list(dat)
#' get_param(dat, "x")
#'
#' # The control and init accessors work the same way
#' dat <- add_control(dat, "x")
#' dat <- set_control(dat, "x", 0.4)
#' get_control_list(dat, c("nsteps", "x"))
#' get_control(dat, "nsteps")
#'
#' dat <- add_init(dat, "x")
#' dat <- set_init(dat, "x", 0.4)
#' get_init_list(dat)
#' get_init(dat, "x")
#'
#' # A control setting with one entry per network layer
#' dat <- set_control(dat, "z", multilayer(0.1, 0.2))
#' get_network_control(dat, network = 2, "z")
#'
#' @order 1
#' @export
get_param <- function(dat, item, override.null.error = FALSE) {
  if (!item %in% names(dat$param)) {
    if (!override.null.error) {
      stop("There is no parameter called `", item,
           "` in the parameter list of the main list object (dat)")
    }
    return(NULL)
  }

  return(dat$param[[item]])
}

#' @rdname get_param
#' @order 4
#' @export
add_param <- function(dat, item) {
  if (item %in% names(dat$param)) {
    stop("Cannot create the parameter, ", item, ": exists already")
  }
  dat$param[[item]] <- NA
  return(dat)
}

#' @rdname get_param
#' @order 3
#' @export
set_param <- function(dat, item, value) {
  if (!item %in% names(dat$param)) {
    dat <- add_param(dat, item)
  }
  dat$param[[item]] <- value
  return(dat)
}

#' @rdname get_param
#' @order 6
#' @export
get_control_list <- function(dat, item = NULL) {
  if (is.null(item)) {
    return(dat$control)
  }

  missing_item <- setdiff(item, names(dat$control))
  if (length(missing_item) > 0L) {
    stop("There is no control value called `",
         paste(missing_item, collapse = ", "),
         "` in the control list of the main list object (dat)")
  }

  return(dat$control[item])
}

#' @rdname get_param
#' @order 5
#' @export
get_control <- function(dat, item, override.null.error = FALSE) {
  if (!item %in% names(dat$control)) {
    if (!override.null.error) {
      stop("There is no control value called `", item,
           "` in the control list of the main list object (dat)")
    }
    return(NULL)
  }
  return(dat$control[[item]])
}

#' @rdname get_param
#' @order 7
#' @export
get_network_control <- function(dat, network, item, override.null.error = FALSE) {
  if (!item %in% names(dat$control)) {
    if (!override.null.error) {
      stop("There is no control value called `", item,
           "` in the control list of the main list object (dat)")
    }
    return(NULL)
  }
  out <- dat$control[[item]]

  if (!is(out, "multilayer")) {
    stop("Control value `", item, "` accessed through `get_network_control` ",
         "is not of class `multilayer`.")
  }
  return(out[[network]])
}

#' @rdname get_param
#' @order 9
#' @export
add_control <- function(dat, item) {
  if (item %in% names(dat$control)) {
    stop("Cannot create the control value, ", item,
         ": exists already")
  }
  dat$control[[item]] <- NA
  return(dat)
}

#' @rdname get_param
#' @order 8
#' @export
set_control <- function(dat, item, value) {
  if (!item %in% names(dat$control)) {
    dat <- add_control(dat, item)
  }
  dat$control[[item]] <- value
  return(dat)
}

#' @rdname get_param
#' @order 11
#' @export
get_init_list <- function(dat, item = NULL) {
  if (is.null(item)) {
    return(dat$init)
  }

  missing_item <- setdiff(item, names(dat$init))
  if (length(missing_item) > 0L) {
    stop("There is no init value called `",
         paste(missing_item, collapse = ", "),
         "` in the init list of the main list object (dat)")
  }
  return(dat$init[item])
}

#' @rdname get_param
#' @order 10
#' @export
get_init <- function(dat, item, override.null.error = FALSE) {
  if (!item %in% names(dat$init)) {
    if (!override.null.error) {
      stop("There is no init value called `", item,
           "` in the init list of the main list object (dat)")
    }
    return(NULL)
  }
  return(dat$init[[item]])
}

#' @rdname get_param
#' @order 13
#' @export
add_init <- function(dat, item) {
  if (item %in% names(dat$init)) {
    stop("Cannot create the init value, ", item,
         ": exists already")
  }
  dat$init[[item]] <- NA
  return(dat)
}

#' @rdname get_param
#' @order 12
#' @export
set_init <- function(dat, item, value) {
  if (!item %in% names(dat$init)) {
    dat <- add_init(dat, item)
  }
  dat$init[[item]] <- value
  return(dat)
}

# Core Attributes --------------------------------------------------------------

#' @rdname get_attr
#' @order 8
#' @export
get_core_attributes <- function() {
  list(
    active = 1L,
    entrTime = 1L,
    exitTime =  NA_integer_,
    unique_id = 1L
  )
}

#' @rdname get_attr
#' @order 6
#' @export
append_core_attr <- function(dat, at, n.new) {
  dat <- append_attr(dat, "active", 1L, n.new)
  dat <- append_attr(dat, "entrTime", as.integer(at), n.new)
  dat <- append_attr(dat, "exitTime", NA_integer_, n.new)
  dat <- update_unique_ids(dat, n.new)
  return(dat)
}

#' @title Create the Unique Identifiers for New Nodes
#'
#' @description This function is called by [append_core_attr()] and
#' appends new `unique_ids` to the created nodes. It also keeps track of
#' the already used `unique_ids` with the `dat$run$last_unique_id`
#' variable.
#'
#' @inheritParams recovery.net
#' @param n.new The number of new nodes to give `unique_ids` to.
#'
#' @inherit recovery.net return
#'
#' @keywords internal
update_unique_ids <- function(dat, n.new) {
  last_unique_id <- dat$run$last_unique_id
  next_unique_ids <- seq_len(n.new) + last_unique_id
  dat$run$last_unique_id <- last_unique_id + as.integer(n.new)
  dat <- append_attr(dat, "unique_id", next_unique_ids, n.new)
  return(dat)
}

#' @title Check that All Attributes in the Main `netsim_dat` Object are of Equal
#'        Length
#'
#' @description Every nodal attribute must hold exactly one value per node, the
#'              number of nodes being the length of the `active` attribute.
#'              This is checked at the end of each time step of [netsim()].
#'              An attribute of the wrong length is usually created by an
#'              [append_attr()] call that does not cover every new node, or by
#'              direct assignment into `dat$run$attr` in place of the accessors.
#'
#' @inheritParams recovery.net
#'
#' @return invisible(TRUE) if everything is correct; an error if not.
#'
#' @keywords internal
check_attr_lengths <- function(dat) {
  attr_list <- raw_get_attr_list(dat)
  attr_lengths <- vapply(attr_list, length, integer(1))

  if (!"active" %in% names(attr_lengths)) {
    stop("The `active` nodal attribute is missing from the main list object ",
         "(dat), so the number of nodes is undefined")
  }

  expected_length <- attr_lengths[["active"]]
  wrong_lengths <- attr_lengths[attr_lengths != expected_length]

  if (length(wrong_lengths) > 0L) {
    stop(
      "Some nodal attributes are not of the correct length \n",
      "Expected length (the length of `active`): ", expected_length, "\n",
      "Wrong length attributes: \n",
      paste0("`", names(wrong_lengths), "`: ", wrong_lengths, "\n",
             collapse = "")
    )
  }

  return(invisible(TRUE))
}

# Unique / Positional identifier converters ------------------------------------

#' @title Convert Unique Identifiers to/from Positional Identifiers
#'
#' @description EpiModel refers to its nodes either by positional identifiers
#'              (`posit_ids`), which describe the position of a node in the
#'              `attr` vector, or by unique identifiers
#'              (`unique_ids`), which allow references to nodes even after
#'              they are deactivated.
#'
#' @section All elements:
#'   When `unique_ids` or `posit_ids` is NULL (default)
#'   the full list of positional IDs or unique IDs is returned.
#'
#' @section Deactivated nodes:
#'   When providing `unique_ids` of deactivated nodes to
#'   `get_posit_ids`, `NA`s are returned instead and a warning is
#'   produced.
#'
#' @inheritParams recovery.net
#' @param unique_ids A vector of node unique identifiers (default = NULL).
#' @param posit_ids A vector of node positional identifiers (default = NULL).
#' @return A vector of unique or positional identifiers.
#'
#' @name unique_id-tools
NULL

#' @rdname unique_id-tools
#' @export
get_unique_ids <- function(dat, posit_ids = NULL) {
  if (is.null(posit_ids)) {
    return(get_attr(dat, "unique_id"))
  }
  unique_ids <- get_attr(dat, "unique_id", posit_ids = posit_ids)
  return(unique_ids)
}

#' @rdname unique_id-tools
#' @export
get_posit_ids <- function(dat, unique_ids = NULL) {
  if (is.null(unique_ids)) {
    return(seq_along(get_attr(dat, "active")))
  }
  posit_ids <- base::match(unique_ids, get_attr(dat, "unique_id"))
  if (any(is.na(posit_ids))) {
    warning(
      "While converting `unique_ids` to `posit_ids`, some `unique_ids`",
      " correspond to deactivated nodes and NAs were produced"
    )
  }
  return(posit_ids)
}

#' @title Are These Nodes Active (Positional IDs)
#'
#' @inheritParams recovery.net
#' @param posit_ids A vector of node positional identifiers.
#'
#' @return A logical vector with TRUE if the node is still active and FALSE
#' otherwise.
#'
#' @export
is_active_posit_ids <- function(dat, posit_ids) {
  active <- get_attr(dat, "active")
  return(active[posit_ids] %in% 1L)
}

#' @title Are These Nodes Active (Unique IDs)
#'
#' @inheritParams recovery.net
#' @param unique_ids A vector of node unique identifiers.
#'
#' @return A logical vector with TRUE if the node is still active and FALSE
#' otherwise.
#'
#' @export
is_active_unique_ids <- function(dat, unique_ids) {
  suppressWarnings({
    posit_ids <- get_posit_ids(dat, unique_ids)
  })
  return(is_active_posit_ids(dat, posit_ids))
}

#' @title Grow a Vector to a Given Size, Padding it With Empty Elements
#'
#' @description Grow a vector to a given size, padding it with `NULL` if `orig` is a `list`
#' and with `NA` otherwise
#'
#' @param orig A vector to grow.
#' @param size The final size of the vector.
#'
#' @return
#' A vector of size `size` padded with `NULL`s or `NA`s at the end.
#'
#' @export
padded_vector <- function(orig, size) {
  if (is.list(orig)) {
    out <- c(orig, vector(mode = "list", size - length(orig)))
  } else {
    out <- c(orig, rep(NA, size - length(orig)))
  }
  return(out)
}

#' Check if `posit_ids` are corrects (numeric and <= NNodes)
#' @noRd
assert_valid_posit_ids <- function(dat, posit_ids) {
  if (!is.numeric(posit_ids))
    stop("`posit_ids` must be `numeric`")
  if (any(posit_ids > length(get_attr(dat, "active"))))
    stop("Some `posit_ids` are larger than the number of nodes ",
         " in the network")
  invisible(posit_ids)
}

#' Check if `time_steps` are corrects (numeric and <= nsteps)
#' @noRd
assert_valid_time_steps <- function(dat, time_steps) {
  if (!is.numeric(time_steps))
    stop("time steps must be `numeric`")
  if ((any(time_steps > get_control(dat, "nsteps"))))
    stop("Some time steps are larger than the number of
          steps planned for for the simulation (control$nsteps)")
  invisible(time_steps)
}
