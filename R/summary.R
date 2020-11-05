
#' @title Summary Model Statistics
#'
#' @description Extracts and prints model statistics solved with \code{dcm}.
#'
#' @param object An \code{EpiModel} object of class \code{dcm}.
#' @param at Time step for model statistics.
#' @param run Model run number, for \code{dcm} class models with multiple runs
#'        (sensitivity analyses).
#' @param digits Number of significant digits to print.
#' @param ... Additional summary function arguments (not used).
#'
#' @details
#' Summary statistics for the main epidemiological outcomes (state and
#' transition size and prevalence) from an \code{dcm} model. Time-specific
#' summary measures are provided, so it is necessary to input a time of
#' interest. Formultiple-run models (sensitivity analyses), input a model run
#' number. See examples below.
#'
#' @seealso \code{\link{dcm}}
#'
#' @method summary dcm
#' @keywords extract
#' @export
#'
#' @examples
#' ## Deterministic SIR model with varying act.rate
#' param <- param.dcm(inf.prob = 0.2, act.rate = 2:4, rec.rate = 1/3,
#'                    a.rate = 0.011, ds.rate = 0.01,
#'                    di.rate = 0.03, dr.rate = 0.01)
#' init <- init.dcm(s.num = 1000, i.num = 1, r.num = 0)
#' control <- control.dcm(type = "SIR", nsteps = 50)
#' mod <- dcm(param, init, control)
#' summary(mod, at = 25, run = 1)
#' summary(mod, at = 25, run = 3)
#' summary(mod, at = 26, run = 3)
#'
summary.dcm <- function(object, at, run = 1, digits = 3, ...) {

  nruns <- object$control$nruns
  type <- object$control$type
  groups <- object$param$groups
  vital <- object$param$vital
  nsteps <- object$control$nsteps

  if (!is.null(object$control$new.mod)) {
    stop("summary method not available for new model types in dcm")
  }

  df <- as.data.frame(object, run = run)

  if (missing(at) || (at > nsteps | at < 1)) {
    stop("Specify at between 1 and ", nsteps)
  }
  df <- df[df$time == at, ]

  ## Prevalence calculations
  df$s.prev <- df$s.num / df$num
  df$i.prev <- df$i.num / df$num
  if (type == "SIR") {
    df$r.prev <- df$r.num / df$num
  }
  if (groups == 2) {
    df$s.prev.g2 <- df$s.num.g2 / df$num.g2
    df$i.prev.g2 <- df$i.num.g2 / df$num.g2
    if (type == "SIR") {
      df$r.prev.g2 <- df$r.num.g2 / df$num.g2
    }
  }

  if (type == "SI") {
    stats <- with(df, c(s.num, s.prev,
                        i.num, i.prev,
                        num, 1,
                        si.flow, NA))
    mat <- matrix(stats, byrow = TRUE, nrow = length(stats) / 2)
    rownames(mat) <- c("Suscept.", "Infect.", "Total", "S -> I")
    if (vital == TRUE) {
      stats <- with(df, c(a.flow, NA,
                          ds.flow, NA,
                          di.flow, NA))
      mat <- rbind(mat, matrix(stats, byrow = TRUE, nrow = length(stats) / 2))
      rownames(mat)[rownames(mat) == ""] <- c("Arrival ->", "S Departure ->",
                                              "I Departure ->")
    }
    if (groups == 2) {
      stats <- with(df, c(s.num.g2, s.prev.g2,
                          i.num.g2, i.prev.g2,
                          num.g2, 1,
                          si.flow.g2, NA))
      mat.g2 <- matrix(stats, byrow = TRUE, nrow = length(stats) / 2)
      if (vital == TRUE) {
        stats <- with(df, c(a.flow.g2, NA,
                            ds.flow.g2, NA,
                            di.flow.g2, NA))
        mat.g2 <- rbind(mat.g2, matrix(stats, byrow = TRUE,
                                       nrow = length(stats) / 2))
      }
      mat <- cbind(mat, mat.g2)
    }
  }
  if (type == "SIR") {
    stats <- with(df, c(s.num, s.prev,
                        i.num, i.prev,
                        r.num, r.prev,
                        num, 1,
                        si.flow, NA,
                        ir.flow, NA))
    mat <- matrix(stats, byrow = TRUE, nrow = length(stats) / 2)
    rownames(mat) <- c("Suscept.", "Infect.", "Recov.", "Total",
                       "S -> I", "I -> R")
    if (vital == TRUE) {
      stats <- c(df$a.flow, NA,
                 df$ds.flow, NA,
                 df$di.flow, NA,
                 df$dr.flow, NA)
      mat <- rbind(mat, matrix(stats, byrow = TRUE, nrow = length(stats) / 2))
      rownames(mat)[rownames(mat) == ""] <- c("Arrival ->",
                                              "S Departure ->",
                                              "I Departure ->",
                                              "R Departure ->")
    }
    if (groups == 2) {
      stats <- with(df, c(s.num.g2, s.prev.g2,
                          i.num.g2, i.prev.g2,
                          r.num.g2, r.prev.g2,
                          num.g2, 1,
                          si.flow.g2, NA,
                          ir.flow.g2, NA))
      mat.g2 <- matrix(stats, byrow = TRUE, nrow = length(stats) / 2)
      if (vital == TRUE) {
        stats <- with(df, c(a.flow.g2, NA,
                            ds.flow.g2, NA,
                            di.flow.g2, NA,
                            dr.flow.g2, NA))
        mat.g2 <- rbind(mat.g2, matrix(stats, byrow = TRUE,
                                       nrow = length(stats) / 2))
      }
      mat <- cbind(mat, mat.g2)
    }
  }
  if (type == "SIS") {
    stats <- with(df, c(s.num, s.prev,
                        i.num, i.prev,
                        num, 1,
                        si.flow, NA,
                        is.flow, NA))
    mat <- matrix(stats, byrow = TRUE, nrow = length(stats) / 2)
    rownames(mat) <- c("Suscept.", "Infect.", "Total", "S -> I", "I -> S")
    if (vital == TRUE) {
      stats <- c(df$a.flow, NA,
                 df$ds.flow, NA,
                 df$di.flow, NA)
      mat <- rbind(mat, matrix(stats, byrow = TRUE, nrow = length(stats) / 2))
      rownames(mat)[rownames(mat) == ""] <- c("Arrival ->", "S Departure ->",
                                              "I Departure ->")
    }
    if (groups == 2) {
      stats <- with(df, c(s.num.g2, s.prev.g2,
                          i.num.g2, i.prev.g2,
                          num.g2, 1,
                          si.flow.g2, NA,
                          is.flow.g2, NA))
      mat.g2 <- matrix(stats, byrow = TRUE, nrow = length(stats) / 2)
      if (vital == TRUE) {
        stats <- with(df, c(a.flow.g2, NA,
                            ds.flow.g2, NA,
                            di.flow.g2, NA))
        mat.g2 <- rbind(mat.g2, matrix(stats, byrow = TRUE,
                                       nrow = length(stats) / 2))
      }
      mat <- cbind(mat, mat.g2)
    }
  }

  if (groups == 1) {
    colnames(mat) <- c("n", "pct")
  }
  if (groups == 2) {
    colnames(mat) <- c("n:g1", "pct:g1", "n:g2", "pct:g2")
  }

  mat <- round(mat, digits)

  # print it
  cat("EpiModel Summary")
  cat("\n=======================")
  cat("\nModel class:", class(object))

  cat("\n\nSimulation Summary")
  cat("\n-----------------------")
  cat("\nModel type:", type)
  cat("\nNo. runs:", nruns)
  cat("\nNo. time steps:", object$nsteps)
  cat("\nNo. groups:", groups)

  if (groups == 1) {
    statsep <- paste(rep("-", 30), collapse = "")
  }
  if (groups == 2) {
    statsep <- paste(rep("-", 60), collapse = "")
  }

  cat("\n\nModel Statistics\n")
  cat(statsep)
  cat("\nTime:", at)
  cat("\t Run:", run, "\n")
  cat(statsep, "\n")
  print(mat, print.gap = 2)
  cat(statsep, "\n")

}


#' @title Summary Model Statistics
#'
#' @description Extracts and prints model statistics simulated with \code{icm}.
#'
#' @param object An \code{EpiModel} object of class \code{icm}.
#' @param at Time step for model statistics.
#' @param digits Number of significant digits to print.
#' @param ... Additional summary function arguments.
#'
#' @details
#' Summary statistics for the main epidemiological outcomes (state and
#' transition size and prevalence) from an \code{icm} model. Time-specific
#' summary measures are provided, so it is necessary to input a time of
#' interest.
#'
#' @seealso \code{\link{icm}}
#'
#' @method summary icm
#' @keywords extract
#' @export
#'
#' @examples
#' ## Stochastic ICM SI model with 3 simulations
#' param <- param.icm(inf.prob = 0.2, act.rate = 1)
#' init <- init.icm(s.num = 500, i.num = 1)
#' control <- control.icm(type = "SI", nsteps = 50,
#'                        nsims = 5, verbose = FALSE)
#' mod <- icm(param, init, control)
#' summary(mod, at = 25)
#' summary(mod, at = 50)
#'
summary.icm <- function(object, at, digits = 3, ...) {

  nsims <- object$control$nsims
  type <- object$control$type
  groups <- object$param$groups
  vital <- object$param$vital
  nsteps <- object$control$nsteps

  if (missing(at) || (at > nsteps | at < 1)) {
    stop("Specify a time step between 1 and ", nsteps)
  }

  df.mn <- as.data.frame(object, out = "mean")
  df.mn <- df.mn[df.mn$time == at, ]

  df.sd <- as.data.frame(object, out = "sd")
  df.sd <- df.sd[df.sd$time == at, ]

  if (type == "SI") {

    ## Prevalence calcs
    s.prev <- df.mn$s.num / df.mn$num
    i.prev <- df.mn$i.num / df.mn$num
    if (groups == 2) {
      s.prev.g2 <- df.mn$s.num.g2 / df.mn$num.g2
      i.prev.g2 <- df.mn$i.num.g2 / df.mn$num.g2
    }

    ## Group 1 stats
    stats <- c(df.mn$s.num, df.sd$s.num, s.prev,
               df.mn$i.num, df.sd$i.num, i.prev,
               df.mn$num, df.sd$num, 1,
               df.mn$si.flow, df.sd$si.flow, NA)
    mat <- matrix(stats, byrow = TRUE, nrow = length(stats) / 3)
    rownames(mat) <- c("Suscept.", "Infect.", "Total", "S -> I")
    if (vital == TRUE) {
      stats <- c(df.mn$a.flow, df.sd$a.flow, NA,
                 df.mn$ds.flow, df.sd$ds.flow, NA,
                 df.mn$di.flow, df.sd$di.flow, NA)
      mat <- rbind(mat, matrix(stats, byrow = TRUE, nrow = length(stats) / 3))
      rownames(mat)[rownames(mat) == ""] <- c("Arrival ->", "S Departure ->",
                                              "I Departure ->")
    }

    ## Group 2 stats
    if (groups == 2) {
      stats <- c(df.mn$s.num.g2, df.sd$s.num.g2, s.prev.g2,
                 df.mn$i.num.g2, df.sd$i.num.g2, i.prev.g2,
                 df.mn$num.g2, df.sd$num.g2, 1,
                 df.mn$si.flow.g2, df.sd$si.flow.g2, NA)
      mat.g2 <- matrix(stats, byrow = TRUE, nrow = length(stats) / 3)
      if (vital == TRUE) {
        stats <- c(df.mn$a.flow.g2, df.sd$a.flow.g2, NA,
                   df.mn$ds.flow.g2, df.sd$ds.flow.g2, NA,
                   df.mn$di.flow.g2, df.sd$di.flow.g2, NA)
        mat.g2 <- rbind(mat.g2, matrix(stats, byrow = TRUE,
                                       nrow = length(stats) / 3))
      }
      mat <- cbind(mat, mat.g2)
    }

  } # end SI summary


  if (type == "SIR") {

    ## Prevalence calcs
    s.prev <- df.mn$s.num / df.mn$num
    i.prev <- df.mn$i.num / df.mn$num
    r.prev <- df.mn$r.num / df.mn$num
    if (groups == 2) {
      s.prev.g2 <- df.mn$s.num.g2 / df.mn$num.g2
      i.prev.g2 <- df.mn$i.num.g2 / df.mn$num.g2
      r.prev.g2 <- df.mn$r.num.g2 / df.mn$num.g2
    }

    ## Group 1 stats
    stats <- c(df.mn$s.num, df.sd$s.num, s.prev,
               df.mn$i.num, df.sd$i.num, i.prev,
               df.mn$r.num, df.sd$r.num, r.prev,
               df.mn$num, df.sd$num, 1,
               df.mn$si.flow, df.sd$si.flow, NA,
               df.mn$ir.flow, df.sd$ir.flow, NA)
    mat <- matrix(stats, byrow = TRUE, nrow = length(stats) / 3)
    rownames(mat) <- c("Suscept.", "Infect.", "Recov.", "Total",
                       "S -> I", "I -> R")
    if (vital == TRUE) {
      stats <- c(df.mn$a.flow, df.sd$a.flow, NA,
                 df.mn$ds.flow, df.sd$ds.flow, NA,
                 df.mn$di.flow, df.sd$di.flow, NA,
                 df.mn$dr.flow, df.sd$dr.flow, NA)
      mat <- rbind(mat, matrix(stats, byrow = TRUE, nrow = length(stats) / 3))
      rownames(mat)[rownames(mat) == ""] <- c("Arrival ->",
                                              "S Departure ->",
                                              "I Departure ->",
                                              "R Departure ->")
    }

    ## Group 2 stats
    if (groups == 2) {
      stats <- c(df.mn$s.num.g2, df.sd$s.num.g2, s.prev.g2,
                 df.mn$i.num.g2, df.sd$i.num.g2, i.prev.g2,
                 df.mn$r.num.g2, df.sd$r.num.g2, r.prev.g2,
                 df.mn$num.g2, df.sd$num.g2, 1,
                 df.mn$si.flow.g2, df.sd$si.flow.g2, NA,
                 df.mn$ir.flow.g2, df.sd$ir.flow.g2, NA)
      mat.g2 <- matrix(stats, byrow = TRUE, nrow = length(stats) / 3)
      if (vital == TRUE) {
        stats <- c(df.mn$a.flow.g2, df.sd$a.flow.g2, NA,
                   df.mn$ds.flow.g2, df.sd$ds.flow.g2, NA,
                   df.mn$di.flow.g2, df.sd$di.flow.g2, NA,
                   df.mn$dr.flow.g2, df.sd$dr.flow.g2, NA)
        mat.g2 <- rbind(mat.g2, matrix(stats, byrow = TRUE,
                                       nrow = length(stats) / 3))
      }
      mat <- cbind(mat, mat.g2)
    }
  } # end SIR summary


  if (type == "SIS") {

    ## Prevalence calcs
    s.prev <- df.mn$s.num / df.mn$num
    i.prev <- df.mn$i.num / df.mn$num
    if (groups == 2) {
      s.prev.g2 <- df.mn$s.num.g2 / df.mn$num.g2
      i.prev.g2 <- df.mn$i.num.g2 / df.mn$num.g2
    }

    ## Group 1 stats
    stats <- c(df.mn$s.num, df.sd$s.num, s.prev,
               df.mn$i.num, df.sd$i.num, i.prev,
               df.mn$num, df.sd$num, 1,
               df.mn$si.flow, df.sd$si.flow, NA,
               df.mn$is.flow, df.sd$is.flow, NA)
    mat <- matrix(stats, byrow = TRUE, nrow = length(stats) / 3)
    rownames(mat) <- c("Suscept.", "Infect.", "Total", "S -> I", "I -> S")
    if (vital == TRUE) {
      stats <- c(df.mn$a.flow, df.sd$a.flow, NA,
                 df.mn$ds.flow, df.sd$ds.flow, NA,
                 df.mn$di.flow, df.sd$di.flow, NA)
      mat <- rbind(mat, matrix(stats, byrow = TRUE, nrow = length(stats) / 3))
      rownames(mat)[rownames(mat) == ""] <- c("Arrival ->", "S Departure ->",
                                            "I Departure ->")
    }

    ## Group 2 stats
    if (groups == 2) {
      stats <- c(df.mn$s.num.g2, df.sd$s.num.g2, s.prev.g2,
                 df.mn$i.num.g2, df.sd$i.num.g2, i.prev.g2,
                 df.mn$num.g2, df.sd$num.g2, 1,
                 df.mn$si.flow.g2, df.sd$si.flow.g2, NA,
                 df.mn$is.flow.g2, df.sd$is.flow.g2, NA)
      mat.g2 <- matrix(stats, byrow = TRUE, nrow = length(stats) / 3)
      if (vital == TRUE) {
        stats <- c(df.mn$a.flow.g2, df.sd$a.flow.g2, NA,
                   df.mn$ds.flow.g2, df.sd$ds.flow.g2, NA,
                   df.mn$di.flow.g2, df.sd$di.flow.g2, NA)
        mat.g2 <- rbind(mat.g2, matrix(stats, byrow = TRUE,
                                       nrow = length(stats) / 3))
      }
      mat <- cbind(mat, mat.g2)
    }
  } # end SIS summary

  if (groups == 1) colnames(mat) <- c("mean", "sd", "pct")
  if (groups == 2) colnames(mat) <- c("mean:g1", "sd:g1", "pct:g1",
                                      "mean:g2", "sd:g2", "pct:g2")
  mat <- round(mat, digits)

  ## Print it
  cat("EpiModel Summary")
  cat("\n=======================")
  cat("\nModel class:", class(object))

  cat("\n\nSimulation Details")
  cat("\n-----------------------")
  cat("\nModel type:", type)
  cat("\nNo. simulations:", nsims)
  cat("\nNo. time steps:", nsteps)
  cat("\nNo. groups:", groups)

  if (groups == 1) {
    statsep <- paste(rep("-", 30), collapse = "")
  }
  if (groups == 2) {
    statsep <- paste(rep("-", 60), collapse = "")
  }
  cat("\n\nModel Statistics\n")
  cat(statsep)
  cat("\nTime:", at, "\n")
  cat(statsep, "\n")
  print(mat, print.gap = 2)
  cat(statsep, "\n")

}


#' @title Summary Model Statistics
#'
#' @description Extracts and prints model statistics simulated with
#' \code{netsim}.
#'
#' @param object An \code{EpiModel} object of class \code{netsim}.
#' @param at Time step for model statistics.
#' @param digits Number of significant digits to print.
#' @param ... Additional summary function arguments.
#'
#' @details
#' Summary statistics for the main epidemiological outcomes (state and
#' transition size and prevalence) from an \code{netsim} model. Time-specific
#' summary measures are provided, so it is necessary to input a time of
#' interest.
#'
#' @seealso \code{\link{netsim}}
#'
#' @method summary netsim
#' @keywords extract
#' @export
#'
#' @examples
#' \dontrun{
#' ## SI Model without Network Feedback
#' # Initialize network and set network model parameters
#' nw <- network_initialize(n = 100)
#' nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#'
#' # Estimate the ERGM models (see help for netest)
#' est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' # Parameters, initial conditions, and controls for model
#' param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.15)
#' init <- init.net(i.num = 10, i.num.g2 = 10)
#' control <- control.net(type = "SI", nsteps = 100, nsims = 5, verbose.int = 0)
#'
#' # Run the model simulation
#' mod <- netsim(est1, param, init, control)
#'
#' summary(mod, at = 1)
#' summary(mod, at = 50)
#' summary(mod, at = 100)
#' }
#'
summary.netsim <- function(object, at, digits = 3, ...) {

  nsims <- object$control$nsims
  type <- object$control$type
  groups <- object$param$groups
  vital <- object$param$vital
  nsteps <- object$control$nsteps

  if (missing(at) || (at > nsteps | at < 1)) {
    stop("Specify at between 1 and ", nsteps)
  }

  df.mn <- as.data.frame(object, out = "mean")
  df.mn <- df.mn[df.mn$time == at, ]

  df.sd <- as.data.frame(object, out = "sd")
  df.sd <- df.sd[df.sd$time == at, ]

  if (type == "SI") {

    ## Prevalence calcs
    s.prev <- df.mn$s.num / df.mn$num
    i.prev <- df.mn$i.num / df.mn$num
    if (groups == 2) {
      s.prev.g2 <- df.mn$s.num.g2 / df.mn$num.g2
      i.prev.g2 <- df.mn$i.num.g2 / df.mn$num.g2
    }

    ## Group 1 stats
    stats <- c(df.mn$s.num, df.sd$s.num, s.prev,
               df.mn$i.num, df.sd$i.num, i.prev,
               df.mn$num, df.sd$num, 1,
               df.mn$si.flow, df.sd$si.flow, NA)
    mat <- matrix(stats, byrow = TRUE, nrow = length(stats) / 3)
    rownames(mat) <- c("Suscept.", "Infect.", "Total", "S -> I")
    if (vital == TRUE) {
      stats <- c(df.mn$a.flow, df.sd$a.flow, NA,
                 df.mn$ds.flow, df.sd$ds.flow, NA,
                 df.mn$di.flow, df.sd$di.flow, NA)
      mat <- rbind(mat, matrix(stats, byrow = TRUE, nrow = length(stats) / 3))
      rownames(mat)[rownames(mat) == ""] <- c("Arrival ->", "S Departure ->",
                                              "I Departure ->")
    }

    ## Group 2 stats
    if (groups == 2) {
      stats <- c(df.mn$s.num.g2, df.sd$s.num.g2, s.prev.g2,
                 df.mn$i.num.g2, df.sd$i.num.g2, i.prev.g2,
                 df.mn$num.g2, df.sd$num.g2, 1,
                 df.mn$si.flow.g2, df.sd$si.flow.g2, NA)
      mat.g2 <- matrix(stats, byrow = TRUE, nrow = 4)
      if (vital == TRUE) {
        stats <- c(df.mn$a.flow.g2, df.sd$a.flow.g2, NA,
                   df.mn$ds.flow.g2, df.sd$ds.flow.g2, NA,
                   df.mn$di.flow.g2, df.sd$di.flow.g2, NA)
        mat.g2 <- rbind(mat.g2, matrix(stats, byrow = TRUE,
                                       nrow = length(stats) / 3))
      }
      mat <- cbind(mat, mat.g2)
    }

  } ## end SI summary


  if (type == "SIR") {

    ## Prevalence calcs
    s.prev <- df.mn$s.num / df.mn$num
    i.prev <- df.mn$i.num / df.mn$num
    r.prev <- df.mn$r.num / df.mn$num
    if (groups == 2) {
      s.prev.g2 <- df.mn$s.num.g2 / df.mn$num.g2
      i.prev.g2 <- df.mn$i.num.g2 / df.mn$num.g2
      r.prev.g2 <- df.mn$r.num.g2 / df.mn$num.g2
    }

    ## Group 1 stats
    stats <- c(df.mn$s.num, df.sd$s.num, s.prev,
               df.mn$i.num, df.sd$i.num, i.prev,
               df.mn$r.num, df.sd$r.num, r.prev,
               df.mn$num, df.sd$num, 1,
               df.mn$si.flow, df.sd$si.flow, NA,
               df.mn$ir.flow, df.sd$ir.flow, NA)
    mat <- matrix(stats, byrow = TRUE, nrow = length(stats) / 3)
    rownames(mat) <- c("Suscept.", "Infect.", "Recov.", "Total",
                       "S -> I", "I -> R")
    if (vital == TRUE) {
      stats <- c(df.mn$a.flow, df.sd$a.flow, NA,
                 df.mn$ds.flow, df.sd$ds.flow, NA,
                 df.mn$di.flow, df.sd$di.flow, NA,
                 df.mn$dr.flow, df.sd$dr.flow, NA)
      mat <- rbind(mat, matrix(stats, byrow = TRUE, nrow = length(stats) / 3))
      rownames(mat)[rownames(mat) == ""] <- c("Arrival ->",
                                              "S Departure ->",
                                              "I Departure ->",
                                              "R Departure ->")
    }

    ## Group 2 stats
    if (groups == 2) {
      stats <- c(df.mn$s.num.g2, df.sd$s.num.g2, s.prev.g2,
                 df.mn$i.num.g2, df.sd$i.num.g2, i.prev.g2,
                 df.mn$r.num.g2, df.sd$r.num.g2, r.prev.g2,
                 df.mn$num.g2, df.sd$num.g2, 1,
                 df.mn$si.flow.g2, df.sd$si.flow.g2, NA,
                 df.mn$ir.flow.g2, df.sd$ir.flow.g2, NA)
      mat.g2 <- matrix(stats, byrow = TRUE, nrow = length(stats) / 3)
      if (vital == TRUE) {
        stats <- c(df.mn$a.flow.g2, df.sd$a.flow.g2, NA,
                   df.mn$ds.flow.g2, df.sd$ds.flow.g2, NA,
                   df.mn$di.flow.g2, df.sd$di.flow.g2, NA,
                   df.mn$dr.flow.g2, df.sd$dr.flow.g2, NA)
        mat.g2 <- rbind(mat.g2, matrix(stats, byrow = TRUE,
                                       nrow = length(stats) / 3))
      }
      mat <- cbind(mat, mat.g2)
    }
  }


  if (type == "SIS") {

    ## Prevalence calcs
    s.prev <- df.mn$s.num / df.mn$num
    i.prev <- df.mn$i.num / df.mn$num
    if (groups == 2) {
      s.prev.g2 <- df.mn$s.num.g2 / df.mn$num.g2
      i.prev.g2 <- df.mn$i.num.g2 / df.mn$num.g2
    }

    ## Group 1 stats
    stats <- c(df.mn$s.num, df.sd$s.num, s.prev,
               df.mn$i.num, df.sd$i.num, i.prev,
               df.mn$num, df.sd$num, 1,
               df.mn$si.flow, df.sd$si.flow, NA,
               df.mn$is.flow, df.sd$is.flow, NA)
    mat <- matrix(stats, byrow = TRUE, nrow = length(stats) / 3)
    rownames(mat) <- c("Suscept.", "Infect.", "Total",
                       "S -> I", "I -> S")
    if (vital == TRUE) {
      stats <- c(df.mn$a.flow, df.sd$a.flow, NA,
                 df.mn$ds.flow, df.sd$ds.flow, NA,
                 df.mn$di.flow, df.sd$di.flow, NA)
      mat <- rbind(mat, matrix(stats, byrow = TRUE, nrow = length(stats) / 3))
      rownames(mat)[rownames(mat) == ""] <- c("Arrival ->", "S Departure ->",
                                              "I Departure ->")
    }

    ## Group 2 stats
    if (groups == 2) {
      stats <- c(df.mn$s.num.g2, df.sd$s.num.g2, s.prev.g2,
                 df.mn$i.num.g2, df.sd$i.num.g2, i.prev.g2,
                 df.mn$num.g2, df.sd$num.g2, 1,
                 df.mn$si.flow.g2, df.sd$si.flow.g2, NA,
                 df.mn$is.flow.g2, df.sd$is.flow.g2, NA)
      mat.g2 <- matrix(stats, byrow = TRUE, nrow = length(stats) / 3)
      if (vital == TRUE) {
        stats <- c(df.mn$a.flow.g2, df.sd$a.flow.g2, NA,
                   df.mn$ds.flow.g2, df.sd$ds.flow.g2, NA,
                   df.mn$di.flow.g2, df.sd$di.flow.g2, NA)
        mat.g2 <- rbind(mat.g2, matrix(stats, byrow = TRUE,
                                       nrow = length(stats) / 3))
      }
      mat <- cbind(mat, mat.g2)
    }
  }

  if (groups == 1) {
    colnames(mat) <- c("mean", "sd", "pct")
  }
  if (groups == 2) {
    colnames(mat) <- c("mean:g1", "sd:g1", "pct:g1",
                       "mean:g2", "sd:g2", "pct:g2")
  }
  mat <- round(mat, digits)

  ## Print it
  cat("\nEpiModel Summary")
  cat("\n=======================")
  cat("\nModel class:", class(object))

  cat("\n\nSimulation Details")
  cat("\n-----------------------")
  cat("\nModel type:", type)
  cat("\nNo. simulations:", nsims)
  cat("\nNo. time steps:", nsteps)
  cat("\nNo. NW groups:", groups)

  if (groups == 1) {
    statsep <- paste(rep("-", 30), collapse = "")
  }
  if (groups == 2) {
    statsep <- paste(rep("-", 60), collapse = "")
  }

  cat("\n\nModel Statistics\n")
  cat(statsep)
  cat("\nTime:", at, "\n")
  cat(statsep, "\n")
  print(mat, print.gap = 2)
  cat(statsep, "\n")

}


#' @title Summary for Network Model Fit
#'
#' @description Prints the summary model fit statistics for an ERGM or STERGM
#'fit.
#'
#' @param object An \code{EpiModel} object of class \code{netest}.
#' @param ... Additional summary function arguments.
#'
#' @method summary netest
#' @keywords internal
#' @export
#'
#' @details
#' This function is simply a wrapper function for \code{summary.ergm} and
#' \code{summary.stergm}. Additionally, if the edges dissolution approximation
#' was used to fit the temporal ERGM, then the dissolution coefficient
#' information will be printed.
#'
summary.netest <- function(object, ...) {

  print(summary(object$fit, ...))

  if (object$edapprox == TRUE) {
    cat("\n")
    print(object$coef.diss)
  }

}
