
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
  if (new.mod == FALSE) {
    cat("\nCompartments:", names(x$epi)[grep("num", names(x$epi))], fill = 60)
    cat("Flows:", names(x$epi)[grep("flow", names(x$epi))], fill = 60)
  } else {
    cat("\nAll Output:", names(x$epi), fill = 60)
  }

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
  cat("\nCompartments:", names(x$epi)[grep("num", names(x$epi))], fill = 60)
  cat("Flows:", names(x$epi)[grep("flow", names(x$epi))], fill = 60)

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
  print(round(x$stats.table.formation, digits = digits))

  if (x$dynamic == TRUE) {
    cat("\nDissolution Diagnostics")
    cat("\n----------------------- \n")
    print(round(x$stats.table.dissolution, digits = digits))
    if (x$coef.diss$model.type == "hetero") {
      cat("----------------------- \n")
      cat("* Heterogeneous dissolution model results averaged over")
    }
  }

  invisible()
}


#' @export
print.netsim <- function(x, ...) {

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
  cat("\nNo. NW modes:", x$param$modes)

  cat("\n\nModel Parameters")
  cat("\n-----------------------\n")
  pToPrint <- which(!(names(x$param) %in% c("modes", "vital")))
  for (i in pToPrint) {
    if (class(x$param[[i]]) == "numeric" && length(x$param[[i]]) > 5) {
      cat(names(x$param)[i], "=", x$param[[i]][1:3], "...", fill = 80)
    } else if (class(x$param[[i]]) == "data.frame") {
      cat(names(x$param)[i], "= <data.frame>\n")
    } else if (class(x$param[[i]]) == "list") {
      cat(names(x$param)[i], "= <list>\n")
    } else {
      cat(names(x$param)[i], "=", x$param[[i]], fill = 80)
    }
  }

  cat("\nModel Output")
  cat("\n-----------------------")
  cat("\nCompartments:", names(x$epi)[grep("num", names(x$epi))], fill = 60)
  cat("Flows:", names(x$epi)[grep("flow", names(x$epi))], fill = 60)
  othOut <- names(x$epi)[-c(grep("num", names(x$epi)), grep("flow", names(x$epi)))]
  if (length(othOut) > 0) {
    cat("Other Output:", othOut, fill = 60)
  }
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


#' @export
print.param.dcm <- function(x, ...) {

  pToPrint <- seq_along(x)

  cat("DCM Parameters")
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
print.param.icm <- function(x, ...) {

  pToPrint <- which(!(names(x) %in% c("vital")))

  cat("ICM Parameters")
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
print.param.net <- function(x, ...) {

  pToPrint <- which(!(names(x) %in% c("vital")))

  cat("Network Model Parameters")
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
  if (is.null(x$new.mod)) {
    pToPrint <- pToPrint[-which(names(x) == "new.mod")]
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
  cat("Built-In Modules:", x$bi.mods, fill = 80)
  if (length(x$user.mods) > 0) {
    cat("User Modules:", x$user.mods, fill = 80)
  }

  invisible()
}

#' @export
print.control.net <- function(x, ...) {

  pToPrint <- which(!grepl(".FUN", names(x)) &
                    names(x) != "set.control.stergm" &
                    !(names(x) %in% c("bi.mods", "user.mods")))

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
  cat("Built-In Modules:", x$bi.mods, fill = 80)
  if (length(x$user.mods) > 0) {
    cat("User Modules:", x$user.mods, fill = 80)
  }

  invisible()
}
