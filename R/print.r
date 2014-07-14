
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
  cat(paste("\nEsimation Method:", estmeth))

  cat("\n\nERGM Model Form")
  cat("\n-----------------------")
  cat("\nFormation: "); print(x$formation)
  cat("Dissolution: "); print(x$dissolution)
  cat("Constraints: "); cat(paste0(as.character(x$constraints)[1],
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
  cat("\nCompartments:", names(x$epi)[grep("num", names(x$epi))], fill = 60)
  cat("Flows:", names(x$epi)[grep("flow", names(x$epi))], fill = 60)
  if (!(is.null(x$network))) {
    cat("Networks:", simnames)
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

  cat("Model Parameters for dcm")
  cat("\n===========================\n")
  for (i in seq_along(x)) {
    cat(names(x)[i], "=", x[[i]], fill = 80)
  }

  invisible()
}

#' @export
print.param.icm <- function(x, ...) {

  pToPrint <- which(!(names(x) %in% c("vital")))

  cat("Model Parameters for icm")
  cat("\n===========================\n")
  for (i in pToPrint) {
    cat(names(x)[i], "=", x[[i]], fill = 80)
  }

  invisible()
}

#' @export
print.param.net <- function(x, ...) {

  pToPrint <- which(!(names(x) %in% c("vital")))

  cat("Model Parameters for net")
  cat("\n===========================\n")
  for (i in pToPrint) {
    cat(names(x)[i], "=", x[[i]], fill = 80)
  }

  invisible()
}


#' @export
print.init.dcm <- function(x, ...) {

  cat("Initial Conditions for dcm")
  cat("\n===========================\n")
  for (i in seq_along(x)) {
    cat(names(x)[i], "=", x[[i]], fill = 80)
  }

  invisible()
}

#' @export
print.init.icm <- function(x, ...) {

  cat("Initial Conditions for icm")
  cat("\n===========================\n")
  for (i in seq_along(x)) {
    cat(names(x)[i], "=", x[[i]], fill = 80)
  }

  invisible()
}

#' @export
print.init.net <- function(x, ...) {

  cat("Initial Conditions for net")
  cat("\n===========================\n")
  for (i in seq_along(x)) {
    cat(names(x)[i], "=", x[[i]], fill = 80)
  }

  invisible()
}

#' @export
print.control.dcm <- function(x, ...) {

  pToPrint <- which(!(names(x) %in% c("dt")))

  cat("Control Settings for dcm")
  cat("\n===========================\n")
  for (i in pToPrint) {
    cat(names(x)[i], "=", x[[i]], fill = 80)
  }

  invisible()
}

#' @export
print.control.icm <- function(x, ...) {

  pToPrint <- which(!grepl(".FUN", names(x)))

  cat("Control Settings for dcm")
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

  cat("Control Settings for net")
  cat("\n===========================\n")
  for (i in pToPrint) {
    cat(names(x)[i], "=", x[[i]], fill = 80)
  }

  invisible()
}
