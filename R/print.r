
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

  if (x$dynamic == TRUE & !is.null(x$stats.table.dissolution)) {
    cat("\nDissolution Diagnostics")
    cat("\n----------------------- \n")
    if (x$coef.diss$dissolution == ~ offset(edges)) {
      print_nwstats_table(x$stats.table.dissolution, digits)
      if (x$coef.diss$model.type == "hetero") {
        cat("----------------------- \n")
        cat("* Heterogeneous dissolution model results averaged over")
      }
    } else {
      cat("Not available when:")
      cat("\n- dissolution formula is not `~ offset(edges)`")
      cat("\n")
    }
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
      cat("Transsmissions:", simnames)
    }
  }
  if (!is.null(x$control$save.other)) {
    cat("\nOther Elements:", x$control$save.other)
  }
  cat("")

  if (nwstats) {
    stats <- x$stats$nwstats
    nsims <- x$control$nsims

    target.stats <- x$nwparam[[network]]$target.stats
    ts.attr.names <- x$nwparam[[network]]$target.stats.names
    if (length(ts.attr.names) != length(target.stats)) {
      target.stats <- target.stats[which(target.stats > 0)]
    }
    ts.out <- data.frame(names = ts.attr.names,
      targets = target.stats)

    merged.stats <- Reduce(
      function(a, x) rbind(a, x[[network]]),
      stats, init = c()
    )

    stats.table.formation <- make_formation_table(merged.stats, ts.out)

    cat("\n\nFormation Diagnostics")
    cat("\n----------------------- \n")
    print_nwstats_table(stats.table.formation, digits = digits)
    cat("\n")

    cat("\nDissolution Diagnostics")
    cat("\n----------------------- \n")

    if (x$control$save.network &&
        ! x$control$tergmLite &&
        x$nwparam[[network]]$coef.diss$dissolution == ~ offset(edges)) {
      diag.sim <- lapply(
        seq_len(x$control$nsims),
        get_network, network = network, x = x
      )
      sim.df <- lapply(diag.sim, as.data.frame)

      dissolution.stats <- make_dissolution_stats(
        sim.df,
        x$nwparam[[network]]$coef.diss,
        x$control$nsteps,
        verbose = FALSE
      )

      print_nwstats_table(dissolution.stats$stats.table.dissolution, digits)

      if (x$nwparam[[network]]$coef.diss$model.type == "hetero") {
        cat("----------------------- \n")
        cat("* Heterogeneous dissolution model results averaged over")
      }
    } else {
      cat("Not available when:")
      cat("\n- `control$tergmLite == TRUE`")
      cat("\n- `control$save.network == FALSE`")
      cat("\n- dissolution formula is not `~ offset(edges)`")
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

#' Format one parameter for printing with the `print.param.xxx` functions
#'
#' @param param_name The name of the parameter to print
#' @param param_value The value of the parameter to print
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
      if (prm == "param_random_set") {
        cat(prm, "= <data.frame> ( dimensions:", dim(x$random.param$param_random_set), ")\n")
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
    if (class(x[[i]]) %in% c("integer", "numeric") && length(x[[i]]) > 10) {
      cat(names(x)[i], "=", x[[i]][1:5], "...", fill = 80)
    } else if (class(x[[i]]) == "data.frame") {
      cat(names(x)[i], "= <data.frame>\n")
    } else if (class(x[[i]]) == "list") {
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
    if (class(x[[i]]) %in% c("integer", "numeric") && length(x[[i]]) > 10) {
      cat(names(x)[i], "=", x[[i]][1:5], "...", fill = 80)
    } else if (class(x[[i]]) == "data.frame") {
      cat(names(x)[i], "= <data.frame>\n")
    } else if (class(x[[i]]) == "list") {
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
    if (class(x[[i]]) %in% c("integer", "numeric") && length(x[[i]]) > 10) {
      cat(names(x)[i], "=", x[[i]][1:5], "...", fill = 80)
    } else if (class(x[[i]]) == "data.frame") {
      cat(names(x)[i], "= <data.frame>\n")
    } else if (class(x[[i]]) == "list") {
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
    names(x) != "set.control.ergm" &
    !grepl("^mcmc\\.control", names(x)) &
    !(names(x) %in% c("bi.mods", "user.mods"))
  )



  cat("Network Model Control Settings")
  cat("\n===============================\n")
  for (i in pToPrint) {
    if (class(x[[i]]) == "formula") {
      cat(names(x)[i], "= "); cat(paste0(as.character(x[[i]])[1],
                                         as.character(x[[i]])[2]), "\n")
    } else if (class(x[[i]]) == "data.frame") {
      cat(names(x)[i], "= <data.frame>\n")
    } else if (class(x[[i]]) == "list") {
      cat(names(x)[i], "= <list>\n")
    } else {
      cat(names(x)[i], "=", x[[i]], fill = 80)
    }
  }

  funToPrint <- names(x)[grep(".FUN", names(x))]
  funToPrint <- funToPrint[-which(funToPrint %in% c("initialize.FUN",
                                                    "verbose.FUN"))]
  if (is.null(x$module.order)) {
    cat("Dynamic Modules:", funToPrint)
  } else {
    order <- unlist(lapply(funToPrint, function(y) which(y == x$module.order)))
    funToPrint.mo <- funToPrint[order]
    funtoPrint.nmo <- funToPrint[-which(funToPrint %in% x$module.order)]
    cat("Dynamic Modules:", funToPrint)
  }

  invisible()
}

#' @title Print Helper For Network Stats Tables
#'
#' @param nwtable a formation or dissolution statistics \code{data.frame}
#' @param digits argument to be passed to \code{round}
#'
#' @keywords internal
print_nwstats_table <- function(nwtable, digits) {
  print(as.data.frame(round(as.matrix(nwtable), digits = digits)))
}
