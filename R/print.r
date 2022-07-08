
#' @export
print.dcm <- function(x, ...) {

  # New model
  new.mod <- ifelse(!is.null(x$control$new.mod), TRUE, FALSE)

  cat("EpiModel Simulation")
  cat("\n=======================")
  cat("\nModel class:", class(x))

  cat("\n\nSimulation Summary")
  cat("\n-----------------------")
  if (new.mod == FALSE) {
    cat("\nModel type:", x$control$type)
  }
  cat("\nNo. runs:", x$control$nruns)
  cat("\nNo. time steps:", x$control$nsteps)
  if (new.mod == FALSE) {
    cat("\nNo. groups:", x$param$groups)
  }

  cat("\n\nModel Parameters")
  cat("\n-----------------------\n")
  pToPrint <- which(!(names(x$param) %in% c("groups", "vital")))
  for (i in pToPrint) {
    cat(names(x$param)[i], "=", x$param[[i]], fill = 60)
  }

  cat("\nModel Output")
  cat("\n-----------------------")
  cat("\nVariables:", names(x$epi), fill = 60)

  invisible()
}


#' @export
print.icm <- function(x, ...) {

  cat("EpiModel Simulation")
  cat("\n=======================")
  cat("\nModel class:", class(x))

  cat("\n\nSimulation Summary")
  cat("\n-----------------------")
  cat("\nModel type:", x$control$type)
  cat("\nNo. simulations:", x$control$nsims)
  cat("\nNo. time steps:", x$control$nsteps)
  cat("\nNo. groups:", x$param$groups)

  cat("\n\nModel Parameters")
  cat("\n-----------------------\n")
  pToPrint <- which(!(names(x$param) %in% c("groups", "vital")))
  for (i in pToPrint) {
    cat(names(x$param)[i], "=", x$param[[i]], fill = 60)
  }

  cat("\nModel Output")
  cat("\n-----------------------")
  cat("\nVariables:", names(x$epi), fill = 60)

  invisible()
}


#' @export
print.netest <- function(x, digits = 3, ...) {

  cat("EpiModel Network Estimation")
  cat("\n=======================")
  cat("\nModel class:", class(x))
  estmeth <- ifelse(x$edapprox == TRUE, "ERGM with Edges Approximation",
                    "Full STERGM Fit")
  cat(paste("\nEstimation Method:", estmeth))

  cat("\n\nModel Form")
  cat("\n-----------------------")
  cat("\nFormation: "); print(x$formation)
  cat("Target Statistics:", x$target.stats)
  cat("\nConstraints: "); cat(paste0(as.character(x$constraints)[1],
                                     as.character(x$constraints)[2]))
  cat("\n\nDissolution: "); cat(as.character(x$coef.diss$dissolution), sep = "")
  cat("\nTarget Statistics:", x$coef.diss$duration)

  invisible()
}

#' @rdname print.netdx
#' @title Utility Function for Printing netdx Object
#' @description Prints basic information and statistics from a \code{netdx}
#'              object.
#' @param x an object of class \code{netdx}
#' @param digits number of digits to print in statistics tables
#' @param ... additional arguments (currently ignored)
#' @details
#' Given a \code{netdx} object, \code{print.netdx} prints the diagnostic method
#' (static/dynamic), number of simulations, and (if dynamic) the number of time
#' steps per simulation used in generating the \code{netdx} object, as well as
#' printing the formation statistics table and (if present) the duration and
#' dissolution statistics tables.  The statistics tables are interpreted as
#' follows.
#'
#' Each row has the name of a particular network statistic.  In the formation
#' table, these correspond to actual network statistics in the obvious way.
#' In the duration and dissolution tables, these correspond to dissolution
#' model dyad types: in a homogeneous dissolution model, all dyads are of the
#' \code{edges} type; in a heterogeneous dissolution model, a dyad with a
#' nonzero \code{nodematch} or \code{nodemix} change statistic in the
#' dissolution model has type equal to that statistic, and has type equal to
#' \code{edges} otherwise.  The statistics of interest for the duration and
#' dissolution tables are, respectively, the mean age of extant edges and the
#' edge dissolution rate, broken down by dissolution model dyad type.  (The
#' current convention is to treat the mean age and dissolution rate for a
#' particular dissolution dyad type as 0 on time steps with no edges of that
#' type; this behavior may be changed in the future.)
#' 
#' The columns are named \code{Target}, \code{Sim Mean}, \code{Pct Diff}, 
#' \code{Sim SE}, \code{Z Score}, \code{SD(1-Sim Mean)}, and 
#' \code{SD(Statistic)}.  The \code{Sim Mean} column refers to the mean 
#' statistic value, across all time steps in all simulations in the dynamic
#' case, and across all sampled networks in all simulations in the static case.
#' The \code{Sim SE} column refers to the standard error in the mean, estimated
#' using \code{\link[=effectiveSize]{coda::effectiveSize}}.  The \code{Target} 
#' column indicates the target value (if present) for the statistic, and the 
#' \code{Pct Diff} column gives \code{(Sim Mean - Target)/Target} when 
#' \code{Target} is present.  The \code{Z Score} column gives 
#' \code{(Sim Mean - Target)/(Sim SE)}.  The \code{SD(1-Sim Mean)} column gives
#' the empirical standard deviation across simulations of the mean statistic
#' value within simulation, and \code{SD(Statistic)} gives the empirical 
#' standard deviation of the statistic value across all the simulated data.
#' @export
print.netdx <- function(x, digits = 3, ...) {

  cat("EpiModel Network Diagnostics")
  cat("\n=======================")
  dxmethod <- ifelse(x$dynamic == TRUE, "Dynamic", "Static")
  cat("\nDiagnostic Method:", dxmethod)
  cat("\nSimulations:", x$nsims)
  if (x$dynamic == TRUE) {
    cat("\nTime Steps per Sim:", x$nsteps)
  }

  cat("\n\nFormation Diagnostics")
  cat("\n----------------------- \n")
  print_nwstats_table(x$stats.table.formation, digits)

  if (x$dynamic == TRUE & !is.null(x$stats.table.duration)) {
    cat("\nDuration Diagnostics")
    cat("\n----------------------- \n")
    print_nwstats_table(x$stats.table.duration, digits)
  }
  if (x$dynamic == TRUE & !is.null(x$stats.table.dissolution)) {
    cat("\nDissolution Diagnostics")
    cat("\n----------------------- \n")
    print_nwstats_table(x$stats.table.dissolution, digits)
  }
  # TODO Remove nodefactor in future release.
  if (x$coef.diss$diss.model.type == "nodefactor") {
    cat("----------------------- \n")
    cat("* Duration and dissolution results are averaged over for dissolution
        models containing a nodefactor term.")
  }
  invisible()
}


#' @export
print.netsim <- function(x, nwstats = TRUE, digits = 3, network = 1, ...) {

  nsims <- x$control$nsims
  if (nsims == 1) {
    simnames <- "sim1"
  }
  if (nsims == 2) {
    simnames <- "sim1 sim2"
  }
  if (nsims > 2) {
    simnames <- paste0("sim1 ... sim", nsims)
  }

  cat("EpiModel Simulation")
  cat("\n=======================")
  cat("\nModel class:", class(x))

  cat("\n\nSimulation Summary")
  cat("\n-----------------------")
  cat("\nModel type:", x$control$type)
  cat("\nNo. simulations:", nsims)
  cat("\nNo. time steps:", x$control$nsteps)
  cat("\nNo. NW groups:", x$param$groups)

  # Parameters
  cat("\n\n")
  print(x$param)

  if (is.null(x$control$type)) {
    cat("\nModel Functions")
    cat("\n-----------------------\n")
    for (i in seq_along(x$control$f.names)) {
      cat(x$control$f.names[i], "\n")
    }
  }

  cat("\nModel Output")
  cat("\n-----------------------")
  cat("\nVariables:", names(x$epi), fill = 60)
  if (!(is.null(x$network))) {
    cat("Networks:", simnames)
  }
  if (!(is.null(x$stats$transmat))) {
    if (!is.null(x$network)) {
      cat("\nTransmissions:", simnames)
    } else {
      cat("Transmissions:", simnames)
    }
  }
  if (!is.null(x$control$save.other)) {
    names_present <- intersect(x$control$save.other, names(x))
    if (length(names_present) > 0) {
      cat("\nOther Elements:", names_present)
    }
  }
  cat("")

  if (nwstats && !is.null(x$stats$nwstats)) {
    stats <- lapply(x$stats$nwstats, 
                    function(stats_list) stats_list[[network]])
    nsims <- x$control$nsims

    target.stats <- x$nwparam[[network]]$target.stats
    ts.attr.names <- x$nwparam[[network]]$target.stats.names
    if (length(ts.attr.names) != length(target.stats)) {
      target.stats <- target.stats[which(target.stats > 0)]
    }
    ts.out <- data.frame(names = ts.attr.names,
      targets = target.stats)

    stats.table.formation <- make_formation_table(stats, ts.out)

    cat("\n\nFormation Diagnostics")
    cat("\n----------------------- \n")
    print_nwstats_table(stats.table.formation, digits = digits)
    cat("\n")

    cat("\nDissolution Diagnostics")
    cat("\n----------------------- \n")

    if (x$control$save.diss.stats &&
        x$control$save.network &&
        ! x$control$tergmLite &&
        ! is.null(x$diss.stats) &&
        x$nwparam[[network]]$coef.diss$dissolution == ~ offset(edges)) {

      dissolution.stats <- make_dissolution_stats(
        lapply(seq_len(x$control$nsims), function(sim) x$diss.stats[[sim]][[network]]),
        x$nwparam[[network]]$coef.diss,
        x$control$nsteps,
        verbose = FALSE
      )
      
      print_nwstats_table(dissolution.stats$stats.table.dissolution, digits)

    } else {
      cat("Not available when:")
      cat("\n- `control$tergmLite == TRUE`")
      cat("\n- `control$save.network == FALSE`")
      cat("\n- `control$save.diss.stats == FALSE`")
      cat("\n- dissolution formula is not `~ offset(edges)`")
      cat("\n- `keep.diss.stats == FALSE` (if merging)")
      cat("\n")
    }
  }

  cat("\n")
  invisible()
}


#' @export
print.disscoef <- function(x, ...) {

  cat("Dissolution Coefficients")
  cat("\n=======================")
  cat("\nDissolution Model: "); cat(as.character(x$dissolution), sep = "")
  cat("\nTarget Statistics:", x$duration)
  cat("\nCrude Coefficient:", x$coef.crude)
  cat("\nMortality/Exit Rate:", x$d.rate)
  cat("\nAdjusted Coefficient:", x$coef.adj)
  cat("\n")

  invisible()
}

#' @title Format One Parameter for Printing with the \code{print.param.xxx}
#'        Functions
#'
#' @param param_name The name of the parameter to print.
#' @param param_value The value of the parameter to print.
#'
#' @keywords internal
format_param <- function(param_name, param_value) {
  if (is.numeric(param_value) && length(param_value) > 10) {
    cat(param_name, "=", param_value[1:10], "...", fill = 80)
  } else if (is.data.frame(param_value)) {
    cat(param_name, "= <data.frame>\n")
  } else if (is.list(param_value)) {
    cat(param_name, "= <list>\n")
  } else if (inherits(param_value, "lm")) {
    cat(param_name, "= <lm/glm>\n")
  } else {
    cat(param_name, "=", param_value, fill = 80)
  }
}

#' @export
print.param.dcm <- function(x, ...) {

  pToPrint <- seq_along(x)

  cat("DCM Parameters")
  cat("\n===========================\n")
  for (i in pToPrint) {
    format_param(names(x)[i], x[[i]])
  }

  invisible()
}

#' @export
print.param.icm <- function(x, ...) {

  pToPrint <- which(!(names(x) %in% c("vital")))

  cat("ICM Parameters")
  cat("\n===========================\n")
  for (i in pToPrint) {
    format_param(names(x)[i], x[[i]])
  }

  invisible()
}

#' @export
print.param.net <- function(x, ...) {

  randoms <- c("random.params", "random.params.values")
  pToPrint <- which(!(names(x) %in% c("vital", randoms)))
  rng_values <- list()
  rng_defs <- NULL

  if (all(randoms %in% names(x))) {
    rng_values <- x$random.params.values
    pToPrint <- pToPrint[! names(x)[pToPrint] %in% names(rng_values)]
  } else if (randoms[1] %in% names(x)) {
    rng_defs <- names(x[[randoms[1]]])
    pToPrint <- pToPrint[! names(x)[pToPrint] %in% rng_defs]
  }

  cat("Fixed Parameters")
    cat("\n---------------------------\n")
  for (i in pToPrint) {
    format_param(names(x)[i], x[[i]])
  }

  if (!is.null(rng_defs)) {
    cat("\nRandom Parameters")
    cat("\n(Not drawn yet)")
    cat("\n---------------------------\n")
    for (prm in rng_defs) {
      if (prm == "param.random.set") {
        cat(prm, "= <data.frame> ( dimensions:",
            dim(x$random.param$param.random.set), ")\n")
      } else {
        cat(prm, "= <function>\n")
      }
    }
  }

  if (length(rng_values) > 0) {
    cat("\nRandom Parameters")
    cat("\n---------------------------\n")
    for (i in seq_along(rng_values)) {
      format_param(names(rng_values)[i], rng_values[[i]])
    }
  }

  invisible()
}

#' @export
print.init.dcm <- function(x, ...) {

  pToPrint <- seq_along(x)

  cat("DCM Initial Conditions")
  cat("\n===========================\n")
  for (i in pToPrint) {
    if (inherits(x[[i]], c("integer", "numeric")) && length(x[[i]]) > 10) {
      cat(names(x)[i], "=", x[[i]][1:5], "...", fill = 80)
    } else if (inherits(x[[i]], "data.frame")) {
      cat(names(x)[i], "= <data.frame>\n")
    } else if (inherits(x[[i]], "list")) {
      cat(names(x)[i], "= <list>\n")
    } else {
      cat(names(x)[i], "=", x[[i]], fill = 80)
    }
  }

  invisible()
}

#' @export
print.init.icm <- function(x, ...) {

  pToPrint <- seq_along(x)

  cat("ICM Initial Conditions")
  cat("\n===========================\n")
  for (i in pToPrint) {
    if (inherits(x[[i]], c("integer", "numeric")) && length(x[[i]]) > 10) {
      cat(names(x)[i], "=", x[[i]][1:5], "...", fill = 80)
    } else if (inherits(x[[i]], "data.frame")) {
      cat(names(x)[i], "= <data.frame>\n")
    } else if (inherits(x[[i]], "list")) {
      cat(names(x)[i], "= <list>\n")
    } else {
      cat(names(x)[i], "=", x[[i]], fill = 80)
    }
  }

  invisible()
}

#' @export
print.init.net <- function(x, ...) {

  pToPrint <- seq_along(x)

  cat("Network Model Initial Conditions")
  cat("\n=================================\n")
  for (i in pToPrint) {
    if (inherits(x[[i]], c("integer", "numeric")) && length(x[[i]]) > 10) {
      cat(names(x)[i], "=", x[[i]][1:5], "...", fill = 80)
    } else if (inherits(x[[i]], "data.frame")) {
      cat(names(x)[i], "= <data.frame>\n")
    } else if (inherits(x[[i]], "list")) {
      cat(names(x)[i], "= <list>\n")
    } else {
      cat(names(x)[i], "=", x[[i]], fill = 80)
    }
  }

  invisible()
}

#' @export
print.control.dcm <- function(x, ...) {

  pToPrint <- seq_along(names(x))
  pToPrint <- pToPrint[-which(names(x) == "new.mod")]
  if (!is.null(x$new.mod)) {
    names(x)[which(names(x) == "new.mod.name")] <- "new.mod"
  }

  cat("DCM Control Settings")
  cat("\n===========================\n")
  for (i in pToPrint) {
    cat(names(x)[i], "=", x[[i]], fill = 80)
  }

  invisible()
}

#' @export
print.control.icm <- function(x, ...) {

  pToPrint <- which(!grepl(".FUN", names(x)) &
                      !(names(x) %in% c("bi.mods", "user.mods")))

  cat("ICM Control Settings")
  cat("\n===========================\n")
  for (i in pToPrint) {
    cat(names(x)[i], "=", x[[i]], fill = 80)
  }
  cat("Base Modules:", x$bi.mods, fill = 80)
  if (length(x$user.mods) > 0) {
    cat("Extension Modules:", x$user.mods, fill = 80)
  }

  invisible()
}

#' @export
print.control.net <- function(x, ...) {

  pToPrint <- which(
    !grepl(".FUN", names(x)) &
    names(x) != "f.args" &
    names(x) != "f.names" &
    names(x) != "set.control.stergm" &
    names(x) != "set.control.tergm" &
    names(x) != "set.control.ergm" &
    !grepl("^mcmc\\.control", names(x)) &
    !(names(x) %in% c("bi.mods", "user.mods"))
  )



  cat("Network Model Control Settings")
  cat("\n===============================\n")
  for (i in pToPrint) {
    if (inherits(x[[i]], "formula")) {
      cat(names(x)[i], "= "); cat(paste0(as.character(x[[i]])[1],
                                         as.character(x[[i]])[2]), "\n")
    } else if (inherits(x[[i]], "data.frame")) {
      cat(names(x)[i], "= <data.frame>\n")
    } else if (inherits(x[[i]], "list")) {
      cat(names(x)[i], "= <list>\n")
    } else {
      cat(names(x)[i], "=", x[[i]], fill = 80)
    }
  }

  if (!is.null(x$module.order)) {
    funToPrint <- x$module.order
  } else {
    funToPrint <- names(x)[grep(".FUN", names(x))]
    funToPrint <- funToPrint[!funToPrint %in% c("initialize.FUN",
                                                "verbose.FUN")]
  }

  cat("Dynamic Modules:", funToPrint)
  cat("\n")

  invisible()
}

#' @title Print Helper For Network Stats Tables
#'
#' @param nwtable A formation or dissolution statistics \code{data.frame}.
#' @param digits Argument to be passed to \code{round}.
#'
#' @keywords internal
print_nwstats_table <- function(nwtable, digits) {
  print(as.data.frame(round(as.matrix(nwtable), digits = digits)))
}
