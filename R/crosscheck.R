
#' @title Cross Checking of Inputs for Deterministic Compartmental Models
#'
#' @description This function checks that the three parameter lists from
#'              \code{\link{param.dcm}}, \code{\link{init.dcm}}, and
#'              \code{\link{control.dcm}} are consistent.
#'
#' @param param An \code{EpiModel} object of class \code{\link{param.dcm}}.
#' @param init An \code{EpiModel} object of class \code{\link{init.dcm}}.
#' @param control An \code{EpiModel} object of class \code{\link{control.dcm}}.
#'
#' @return
#' This function returns no objects.
#'
#' @export
#' @keywords internal
#'
crosscheck.dcm <- function(param, init, control) {

  # Main class check --------------------------------------------------------
  if (!inherits(param, "param.dcm")) {
    stop("param must an object of class param.dcm", call. = FALSE)
  }
  if (!inherits(init, "init.dcm")) {
    stop("init must an object of class init.dcm", call. = FALSE)
  }
  if (!inherits(control, "control.dcm")) {
    stop("control must an object of class control.dcm", call. = FALSE)
  }

  # Parameter checks for base models ----------------------------------
  if (is.null(control$new.mod)) {

    ## Defaults
    if (is.null(param$act.rate)) {
      param$act.rate <- 1
    }
    if (is.null(param$vital)) {
      if (!is.null(param$a.rate) |
          !is.null(param$ds.rate) |
          !is.null(param$di.rate) |
          !is.null(param$dr.rate)) {
        param$vital <- TRUE
      } else {
        param$vital <- FALSE
      }
    }

    if (any(grepl(".g2", names(param))) == TRUE) {
      param$groups <- 2
    } else {
      param$groups <- 1
    }

    if (param$groups == 2 && (is.null(param$balance) ||
                                !(param$balance %in% c("g1", "g2")))) {
      stop("Specify balance=\"g1\" or balance=\"g2\" with 2-group models",
           call. = FALSE)
    }

    ## Error checks
    # Specify inf.prob
    if (is.null(param$inf.prob)) {
      stop("Specify inf.prob in param.dcm", call. = FALSE)
    }

    # Check that rec.rate is supplied for SIR models
    if (control$type %in% c("SIR", "SIS") & is.null(param$rec.rate)) {
      stop("Specify rec.rate in param.dcm", call. = FALSE)
      if (param$groups == 2 & is.null(param$rec.rate.g2)) {
        stop("Specify rec.rate.g2 in param.dcm", call. = FALSE)
      }
    }

    # Check that r.num is supplied for SIR models
    if (control$type == "SIR" & is.null(init$r.num)) {
      stop("Specify r.num in init.dcm", call. = FALSE)
      if (param$groups == 2 & is.null(init$r.num.g2)) {
        stop("Specify r.num.g2 in init.dcm", call. = FALSE)
      }
    }

    # Check that groups implied by init and params are consistent
    if (any(grepl(".g2", names(init))) == TRUE) {
      init.groups <- 2
    } else {
      init.groups <- 1
    }
    if (param$groups == 2 && init.groups == 1) {
      stop("Group 2 parameters specified in param.dcm,
           \rbut missing group 2 initial states in init.dcm",
           call. = FALSE)
    }
    if (param$groups == 1 && init.groups == 2) {
      stop("Group 2 initial stats specified in init.dcm,
           but missing group 2 parameters in param.dcm",
           call. = FALSE)
    }

    # Over-specified initial conditions
    if (control$type != "SIR" & any(c("r.num", "r.num.g2") %in% names(init))) {
      stop("Specified initial number recovered for non-SIR model",
           call. = FALSE)
    }

    # Deprecated parameters
    if (!is.null(param$trans.rate)) {
      stop("The trans.rate parameter is deprecated. Use the inf.prob parameter instead.",
           call. = FALSE)
    }
    if (!is.null(param$trans.rate.g2)) {
      stop("The trans.rate.g2 parameter is deprecated. Use the inf.prob.g2 parameter instead.",
           call. = FALSE)
    }
  }

  on.exit(assign("param", param, pos = parent.frame()))
}


#' @title Cross Checking of Inputs for Stochastic Individual Contact Models
#'
#' @description This function checks that the three parameter lists from
#'              \code{\link{param.icm}}, \code{\link{init.icm}}, and
#'              \code{\link{control.icm}} are consistent.
#'
#' @param param An \code{EpiModel} object of class \code{\link{param.icm}}.
#' @param init An \code{EpiModel} object of class \code{\link{init.icm}}.
#' @param control An \code{EpiModel} object of class \code{\link{control.icm}}.
#'
#' @return
#' This function returns no objects.
#'
#' @export
#' @keywords internal
#'
crosscheck.icm <- function(param, init, control) {

  ## Main class check
  if (!inherits(param, "param.icm")) {
    stop("param must an object of class param.icm", call. = FALSE)
  }
  if (!inherits(init, "init.icm")) {
    stop("init must an object of class init.icm", call. = FALSE)
  }
  if (!inherits(control, "control.icm")) {
    stop("control must an object of class control.icm", call. = FALSE)
  }

  if (control$skip.check == FALSE) {

    ## Check that rec.rate is supplied for SIR models
    if (control$type %in% c("SIR", "SIS")) {
      if (is.null(param$rec.rate)) {
        stop("Specify rec.rate in param.icm", call. = FALSE)
      }
      if (param$groups == 2 & is.null(param$rec.rate.g2)) {
        stop("Specify rec.rate.g2 in param.icm", call. = FALSE)
      }
    }


    ## Check that paramets and init are supplied for SIR models
    if (control$type == "SIR") {
      if (is.null(init$r.num)) {
        stop("Specify r.num in init.icm", call. = FALSE)
      }
      if (param$groups == 2 & is.null(init$r.num.g2)) {
        stop("Specify r.num.g2 in init.icm", call. = FALSE)
      }
    }

    ## Check that groups implied by init and params are consistent
    if (any(grepl(".g2", names(init))) == TRUE) {
      init.groups <- 2
    } else {
      init.groups <- 1
    }
    if (param$groups == 2 && init.groups == 1) {
      stop("Group 2 parameters specified in param.dcm, but missing group 2, ",
           "initial states in init.icm", call. = FALSE)
    }
    if (param$groups == 1 && init.groups == 2) {
      stop("Group 2 initial stats specified in init.dcm, but missing group 2 ",
           "parameters in param.icm", call. = FALSE)
    }

    ## Deprecated parameters
    bim <- grep(".FUN", names(formals(control.icm)), value = TRUE)
    um <- which(grepl(".FUN", names(control)) & !(names(control) %in% bim))
    if (length(um) == 0 && !is.null(control$type)) {
      if (!is.null(param$trans.rate)) {
        stop("The trans.rate parameter is deprecated. Use the inf.prob ",
             "parameter instead.", call. = FALSE)
      }
      if (!is.null(param$trans.rate.g2)) {
        stop("The trans.rate.g2 parameter is deprecated. Use the inf.prob.g2 ",
             "parameter instead.", call. = FALSE)
      }
    }

  }


  ## Assign modules based on group parameter
  if (!is.null(control$type)) {
    def <- grep(".FUN",names(control))
    args <- names(control)[def]

    if (param$groups == 1) {
      for (i in 1:length(args)) {
        if (is.null(control[[args[i]]])) {
          temp <- get(gsub(".FUN",".icm",args[i]))
          control[[args[i]]] <- temp
        }
      }
      message("Default modules set to appropriate two-group functions. See documentation ",
              "for details.")
    }
    else {
      for (i in 1:length(args)) {
        if (is.null(control[[args[i]]])) {
          temp <- get(gsub(".FUN",".icm.bip",args[i]))
          control[[args[i]]] <- temp
        }
      }
     message("Default modules set to appropriate two-group functions. See documentation ",
             "for details.")
    }
  }

  ## In-place assignment to update param and control
  on.exit(assign("param", param, pos = parent.frame()))
  on.exit(assign("control", control, pos = parent.frame()), add = TRUE)
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

  if (!is.null(control$type) && length(control$user.mods) == 0) {

  if (control$start == 1 && control$skip.check == FALSE) {

    # Main class check --------------------------------------------------------
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

    nw <- x$fit$newnetwork

    # Defaults ----------------------------------------------------------------

    # Is status in network formation formula?
    statOnNw <- ("status" %in% get_formula_term_attr(x$formation, nw))

    # Set dependent modeling defaults if vital or status on nw
    if (is.null(control$depend)) {
      if (param$vital == TRUE | statOnNw == TRUE) {
        control$depend <- TRUE
      } else {
        control$depend <- FALSE
      }
    }

    bip <- length(unique(get.vertex.attribute(nw, "group")))
    bip <- ifelse(bip == 2, TRUE, FALSE)

    if (bip == TRUE & is.null(control$pid.prefix)) {
      control$pid.prefix <- c("F", "M")
    }

    if (statOnNw == TRUE && is.null(control$attr.rules$status)) {
      control$attr.rules$status <- "s"
    }

    # Using tergmLite; then 'depend = TRUE':

    if (control$tgl == TRUE) {
      control$depend <- TRUE
    }


    # Checks ------------------------------------------------------------------

    # Check that prevalence in NW attr status and initial conditions match
    if (statOnNw == TRUE) {
      nw1 <- sum(get.vertex.attribute(nw, "status") == 1)
      init1 <- sum(unlist(init[grep("i.num", names(init), value = TRUE)]))
      if ("i.num" %in% names(init) && nw1 != init1) {
        warning("Overriding init infected settings with network status attribute",
                call. = FALSE, immediate. = TRUE)
        if (interactive()) Sys.sleep(4)
      }
    }

    # If status not in formation formula but set on original network, state that it
    #   will be ignored
    if (statOnNw == FALSE & "status" %in% names(nw$val[[1]])) {
      warning("Overriding status vertex attribute on network with init.net conditions",
              call. = FALSE, immediate. = TRUE)
      if (interactive()) Sys.sleep(4)
    }


    # Check consistency of status vector to network structure
    if (!is.null(init$status.vector)) {
      if (length(init$status.vector) != network.size(nw)) {
        stop("Length of status.vector is unequal to size of initial network")
      }
      svals <- sort(unique(init$status.vector))
      if (control$type == "SIR") {
        if (any(svals %in% c("s", "i", "r") == FALSE)) {
          stop("status.vector contains values other than \"s\", \"i\", and \"r\" ",
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
    if (bip == TRUE & is.null(init$i.num.g2) &
          is.null(init$status.vector) & statOnNw == FALSE) {
      stop("Specify i.num.g2 for two-group model simulations", call. = FALSE)
    }

    # Recovery rate and initial recovered checks
    if (control$type %in% c("SIR", "SIS")) {
      if (is.null(param$rec.rate)) {
        stop("Specify rec.rate in param.net", call. = FALSE)
      }
      if (bip == TRUE & is.null(param$rec.rate.g2)) {
        stop("Specify rec.rate.g2 in param.net", call. = FALSE)
      }
    }
    if (control$type == "SIR") {
      if (is.null(init$r.num) & is.null(init$status.vector) & statOnNw == FALSE) {
        stop("Specify r.num in init.net", call. = FALSE)
      }
      if (bip == TRUE & is.null(init$r.num.g2) & is.null(init$status.vector) &
          statOnNw == FALSE) {
        stop("Specify r.num.g2 in init.net", call. = FALSE)
      }
    }

    # Check demographic parameters for two-group models
    if (bip == TRUE & param$vital == TRUE) {
      if (is.null(param$a.rate.g2)) {
        stop("Specify a.rate.g2 in param.net", call. = FALSE)
      }
      if (is.null(param$ds.rate.g2)) {
        stop("Specify ds.rate.g2 in param.net", call. = FALSE)
      }
      if (is.null(param$di.rate.g2)) {
        stop("Specify di.rate.g2 in param.net", call. = FALSE)
      }
      if (control$type == "SIR") {
        if (is.null(param$dr.rate.g2)) {
          stop("Specify dr.rate.g2 in param.net", call. = FALSE)
        }
      }
    }


    ## Deprecated parameters
    bim <- grep(".FUN", names(formals(control.net)), value = TRUE)
    um <- which(grepl(".FUN", names(control)) & !(names(control) %in% bim))
    if (length(um) == 0 && !is.null(control$type)) {
      if (!is.null(param$trans.rate)) {
        stop("The trans.rate parameter is deprecated. Use the inf.prob ",
             "parameter instead.", call. = FALSE)
      }
      if (!is.null(param$trans.rate.g2)) {
        stop("The trans.rate.g2 parameter is deprecated. Use the inf.prob.g2 ",
             "parameter instead.", call. = FALSE)
      }
    }

  }

  if (control$start > 1) {

    control$depend <- TRUE

    if (control$skip.check == FALSE) {
      if (class(x) != "netsim") {
        stop("x must be a netsim object if control setting start > 1",
             call. = FALSE)
      }
      if (is.null(x$attr)) {
        stop("x must contain attr to restart simulation, see save.other ",
             "control setting", call. = FALSE)
      }
      if (is.null(x$network)) {
        stop("x must contain network object to restart simulation, ",
             "see save.network control setting", call. = FALSE)
      }
      if (control$nsteps < control$start) {
        stop("control setting nsteps must be > control setting start in ",
             "restarted simulations", call. = FALSE)
      }
      if (control$start > x$control$nsteps + 1) {
        stop("control setting start must be 1 greater than the nsteps in the ",
             "prior simulation", call. = FALSE)
      }


    }

  }


  ## Assign modules based on group parameter
  if (!is.null(control$type)) {
    def <- grep(".FUN",names(control))
    args <- names(control)[def]
    flag <- length(grep(".g2",names(param)))

    if (flag == 0) {
      for (i in 1:length(args)) {
        if (is.null(control[[args[i]]])) {
          temp <- get(gsub(".FUN",".net",args[i]))
          control[[args[i]]] <- temp
        }
      }
    }
    else {
      for (i in 1:length(args)) {
        if (is.null(control[[args[i]]])) {
          temp <- get(gsub(".FUN",".net.grp",args[i]))
          control[[args[i]]] <- temp
        }
      }
    }
  }

  ## In-place assignment to update param and control
  assign("param", param, pos = parent.frame())
  assign("control", control, pos = parent.frame())
  }


  if (!is.null(control$type) && length(control$user.mods) > 0) {
    stop("Control setting 'type' must be NULL if any user-specified modules specified.",
         call. = FALSE)
  }

  if (is.null(control$type)) {
    control$type <- "SI"
  }

  if (is.null(control$type) && length(grep("rec", names(param))) != 0){
    control$type <- "SIR"
  }

  ## In-place assignment to update param and control
  assign("param", param, pos = parent.frame())
  assign("control", control, pos = parent.frame())
}
