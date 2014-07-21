
#' @export
print.dcm <- function(x, ...) {

  # New model
  new.mod <- ifelse(!is.null(x$control$new.mod), TRUE, FALSE)

  cat("EpiModel Object")
  cat("\n=======================")
  cat("\nModel class:", class(x))

  cat("\n\nSimulation Summary")
  cat("\n-----------------------")
  if (new.mod == FALSE) {
    cat("\nModel type:", x$control$type)
  }
  cat("\nNo. runs:", x$control$nruns)
  cat("\nNo. time steps:", max(x$control$dt))
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

  cat("EpiModel Object")
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
print.netest <- function(x, digits=3, ...) {

  cat("EpiModel Object")
  cat("\n=======================")
  cat("\nModel class:", class(x))
  estmeth <- ifelse(x$edapprox == TRUE, "ERGM with Edges Approximation",
                                        "Full STERGM Fit")
  cat(paste("\nEstimation Method:", estmeth))

  cat("\n\nModel Form")
  cat("\n-----------------------")
  cat("\nFormation: "); print(x$formation)
  cat("Formation Targets:", x$target.stats)
  cat("\nDissolution: "); print(x$dissolution)
  cat("Edge Duration Target:", x$coef.diss$duration)
  cat("\nConstraints: "); cat(paste0(as.character(x$constraints)[1],
                                   as.character(x$constraints)[2]))

  invisible()
}


#' @export
print.netdx <- function(x, digits = 3, ...) {

  cat("EpiModel Network Diagnostics")
  cat("\n=======================")
  cat("\nNo. Simulations:", x$nsims)
  cat("\nNo. Time Steps:", x$nsteps)

  cat("\n\nFormation Diagnostics")
  cat("\n----------------------- \n")
  print(round(x$stats.table.formation, digits = digits))

  cat("\nDuration Diagnostics")
  cat("\n----------------------- \n")
  print(round(x$stats.table.duration, digits = digits))

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

  cat("EpiModel Object")
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
    cat(names(x$param)[i], "=", x$param[[i]], fill = 60)
  }

  cat("\nModel Output")
  cat("\n-----------------------")
  cat("\nCompartments:", names(x$epi)[grep("num", names(x$epi))])
  cat("\nFlows:", names(x$epi)[grep("flow", names(x$epi))])
  if (!(is.null(x$network))) {
    cat("\nNetworks:", simnames)
  }
  if (!(is.null(x$stats$transmat))) {
    cat("\nTransmissions:", simnames)
  }
  cat("")

  invisible()
}


#' @export
print.disscoef <- function(x, ...) {

  cat("Dissolution Coefficients")
  cat("\n=======================")
  cat("\nDissolution Model: "); ; print(x$dissolution)
  cat("Edge Duration:", x$duration)
  cat("\nAdjusted Coefficient:", x$coef.adj)
  cat("\nCrude Coefficient:", x$coef.crude)

  invisible()
}


#' @export
print.param.dcm <- function(x, ...) {

  pToPrint <- seq_along(x)

  cat("DCM Parameters")
  cat("\n===========================\n")
  for (i in pToPrint) {
    if (class(x[[i]]) == "numeric" && length(x[[i]]) > 10) {
      cat(names(x)[i], "=", x[[i]][1:3], "...", fill = 80)
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
    if (class(x[[i]]) == "numeric" && length(x[[i]]) > 5) {
      cat(names(x)[i], "=", x[[i]][1:3], "...", fill = 80)
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
    if (class(x[[i]]) == "numeric" && length(x[[i]]) > 5) {
      cat(names(x)[i], "=", x[[i]][1:3], "...", fill = 80)
    } else {
      cat(names(x)[i], "=", x[[i]], fill = 80)
    }
  }

  invisible()
}

#' @export
print.init.dcm <- function(x, ...) {

  cat("DCM Initial Conditions")
  cat("\n===========================\n")
  for (i in seq_along(x)) {
    cat(names(x)[i], "=", x[[i]], fill = 80)
  }

  invisible()
}

#' @export
print.init.icm <- function(x, ...) {

  cat("ICM Initial Conditions")
  cat("\n===========================\n")
  for (i in seq_along(x)) {
    cat(names(x)[i], "=", x[[i]], fill = 80)
  }

  invisible()
}

#' @export
print.init.net <- function(x, ...) {

  cat("Network Model Initial Conditions")
  cat("\n=================================\n")
  for (i in seq_along(x)) {
    cat(names(x)[i], "=", x[[i]], fill = 80)
  }

  invisible()
}

#' @export
print.control.dcm <- function(x, ...) {

  pToPrint <- which(!(names(x) %in% c("dt")))

  cat("DCM Control Settings")
  cat("\n===========================\n")
  for (i in pToPrint) {
    cat(names(x)[i], "=", x[[i]], fill = 80)
  }

  invisible()
}

#' @export
print.control.icm <- function(x, ...) {

  pToPrint <- which(!grepl(".FUN", names(x)))

  cat("ICM Control Settings")
  cat("\n===========================\n")
  for (i in pToPrint) {
    cat(names(x)[i], "=", x[[i]], fill = 80)
  }

  invisible()
}

#' @export
print.control.net <- function(x, ...) {

  pToPrint <- which(!grepl(".FUN", names(x)) &
                    names(x) != "set.control.stergm")

  cat("Network Model Control Settings")
  cat("\n===============================\n")
  for (i in pToPrint) {
    if (class(x[[i]]) == "formula") {
      cat(names(x)[i], "= "); cat(paste0(as.character(x[[i]])[1],
                                         as.character(x[[i]])[2]), "\n")
    } else if (class(x[[i]]) == "data.frame") {
      cat(names(x)[i], "= <data.frame>")
    } else {
      cat(names(x)[i], "=", x[[i]], fill = 80)
    }
  }

  invisible()
}
