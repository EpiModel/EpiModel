
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
  print(round(x$stats.table.formation, digits = digits))

  if (x$dynamic == TRUE & !is.null(x$stats.table.dissolution)) {
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
print.netsim <- function(x, formation.stats = FALSE, ...) {

  nsims <- x$control$nsims
  if (nsims == 1) {
    simnames <- "sim1"
  }
  if (nsims == 2) {
    simnames <- "sim1 s.g2"
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

  cat("\n\nModel Parameters")
  cat("\n-----------------------\n")
  pToPrint <- which(!(names(x$param) %in% c("groups", "vital")))
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

  cat("\nModel Functions")
  cat("\n-----------------------\n")
  for (i in 1:length(x$control$f.args)) {
    cat(x$control$f.names[i])
  }
  #cat("\n")

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

  if (formation.stats == TRUE) {

    stats <- x$stats$nwstats
    nsims <- x$control$nsims

    ## Merged stats across all simulations
    if (nsims > 1) {
      merged.stats <- matrix(NA, nrow = nrow(stats[[1]]) * nsims,
                             ncol = ncol(stats[[1]]))
      for (i in 1:ncol(stats[[1]])) {
        merged.stats[,i] <- as.numeric(sapply(stats, function(x) c(x[,i])))
      }
      colnames(merged.stats) <- colnames(stats[[1]])
    } else {
      merged.stats <- stats[[1]]
    }

    ## Calculate mean/sd from merged stats
    stats.means <- colMeans(merged.stats)
    stats.sd <- apply(merged.stats, 2, sd)
    stats.table <- data.frame(sorder = 1:length(names(stats.means)),
                              names = names(stats.means),
                              stats.means, stats.sd)

    ## Get stats from for target statistics
    ts.attr.names <- x$nwparam[[1]]$target.stats.names
    target.stats <- x$nwparam[[1]]$target.stats
    if (length(ts.attr.names) != length(target.stats)) {
      target.stats <- target.stats[which(target.stats > 0)]
    }
    ts.out <- data.frame(names = ts.attr.names,
                         targets = target.stats)

    ## Create stats.formation table for output
    stats.table <- merge(ts.out, stats.table, all = TRUE)
    stats.table <- stats.table[do.call("order",
                                       stats.table[, "sorder", drop = FALSE]), , drop = FALSE]
    rownames(stats.table) <- stats.table$names

    stats.table$reldiff <- 100 * (stats.table$stats.means - stats.table$targets) /
      stats.table$targets
    stats.table.formation <- stats.table[, c(2, 4, 6, 5)]
    colnames(stats.table.formation) <- c("Target", "Sim Mean", "Pct Diff", "Sim SD")

    cat("\n\nFormation Diagnostics")
    cat("\n----------------------- \n")
    print(round(stats.table.formation, digits = 2))
    cat("")
  }

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
    } else if (inherits(x[[i]], "data.frame")) {
      cat(names(x)[i], "= <data.frame>\n")
    } else if (inherits(x[[i]], "list")) {
      cat(names(x)[i], "= <list>\n")
    } else if (inherits(x[[i]], "lm")) {
      cat(names(x)[i], "= <lm/glm>\n")
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

  pToPrint <- which(!grepl(".FUN", names(x)) &
                      names(x) != "f.args" &
                      names(x) != "f.names" &
                      names(x) != "set.control.stergm" &
                      names(x) != "set.control.ergm" &
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

  funToPrint <- names(x)[grep(".FUN", names(x))]
  funToPrint <- funToPrint[-which(funToPrint %in% c("initialize.FUN",
                                                    "verbose.FUN"))]
  if(is.null(x$module.order)) {
    cat("Dynamic Modules:", funToPrint)
  } else {
    order <- unlist(lapply(funToPrint, function(y) which(y == x$module.order)))
    funToPrint.mo <- funToPrint[order]
    funtoPrint.nmo <- funToPrint[-which(funToPrint %in% x$module.order)]
    cat("Dynamic Modules:", funToPrint)
  }

  invisible()
}
