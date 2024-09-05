#' @title Plot Compartment Diagram for Epidemic Models
#'
#' @description Plots a compartment flow diagram for deterministic compartmental
#'              models, stochastic individual contact models, and stochastic
#'              network models.
#'
#' @param x An \code{EpiModel} object of class \code{dcm}, \code{icm}, or
#'        \code{netsim}.
#' @param at Time step for model statistics.
#' @param digits Number of significant digits to print.
#' @param ... Additional arguments passed to plot (not currently used).
#'
#' @details
#' The \code{comp_plot} function provides a visual summary of an epidemic model
#' at a specific time step. The information contained in \code{comp_plot} is the
#' same as in the \code{summary} functions for a model, but presented
#' graphically as a compartment flow diagram.
#'
#' For \code{dcm} class plots, specify the model run number if the model
#' contains multiple runs, as in a sensitivity analysis. For \code{icm} and
#' \code{netsim} class plots, the \code{run} argument is not used; the plots
#' show the means and standard deviations across simulations at the specified
#' time step.
#'
#' These plots are currently limited to one-group models for each of the three
#' model classes. That functionality may be expanded in future software
#' releases.
#'
#' @export
#' @keywords plot
#'
#' @examples
#' ## Example 1: DCM SIR model with varying act.rate
#' param <- param.dcm(inf.prob = 0.2, act.rate = 5:7,
#'                    rec.rate = 1/3, a.rate = 1/90, ds.rate = 1/100,
#'                    di.rate = 1/35, dr.rate = 1/100)
#' init <- init.dcm(s.num = 1000, i.num = 1, r.num = 0)
#' control <- control.dcm(type = "SIR", nsteps = 25, verbose = FALSE)
#' mod1 <- dcm(param, init, control)
#' comp_plot(mod1, at = 25, run = 3)
#'
#' ## Example 2: ICM SIR model with 3 simulations
#' param <- param.icm(inf.prob = 0.2, act.rate = 3, rec.rate = 1/50,
#'                    a.rate = 1/100, ds.rate = 1/100,
#'                    di.rate = 1/90, dr.rate = 1/100)
#' init <- init.icm(s.num = 500, i.num = 1, r.num = 0)
#' control <- control.icm(type = "SIR", nsteps = 25,
#'                        nsims = 3, verbose = FALSE)
#' mod2 <- icm(param, init, control)
#' comp_plot(mod2, at = 25, digits = 1)
#'
comp_plot <- function(x, at, digits, ...) {
  UseMethod("comp_plot")
}

#' @method comp_plot netsim
#' @rdname comp_plot
#' @export
comp_plot.netsim <- function(x, at = 1, digits = 3, ...) {

  comp_plot.icm(x = x, at = at, digits = digits, ...)

}

#' @method comp_plot icm
#' @rdname comp_plot
#' @export
comp_plot.icm <- function(x, at = 1, digits = 3, ...) {

  # Variables
  nsteps <- x$control$nsteps
  dis.type <- x$control$type
  vital <- x$param$vital

  # Standardize groups
  if (inherits(x, "icm")) {
    groups <- x$param$groups
  }
  if (inherits(x, "netsim")) {
    groups <- x$param$groups
  }
  if (groups != 1) {
    stop("Only 1-group models currently supported",
         call. = FALSE)
  }

  # Time
  if (at > nsteps || at < 1) {
    stop("Specify a timestep between 1 and ", nsteps,
         call. = FALSE)
  }

  ## Dataframe subsets for plots
  df.mn <- as.data.frame(x, out = "mean")
  df.mn <- round(df.mn[at == df.mn$time, ], digits)
  df.sd <- as.data.frame(x, out = "sd")
  df.sd <- round(df.sd[at == df.sd$time, ], digits)

  ## Change graphical parameters
  ops <- list(mar = par()$mar, mfrow = par()$mfrow, mgp = par()$mgp)
  par(mar = c(0, 0, 2, 0))
  options(scipen = 10)

  ## Main Plot
  plot(0:100, 0:100, type = "n", axes = FALSE)
  title(main = paste(dis.type, "Model Diagram"))
  mtext(paste0("Simulation means(sd) | time=", at),
        side = 3, cex = 0.8, line = -1)

  ## 1. SI Model
  if (dis.type == "SI" && groups == 1) {
    mbox(22, 40, "Susceptible", paste0(df.mn$s.num, "(", df.sd$s.num, ")"))
    mbox(57, 40, "Infected", paste0(df.mn$i.num, "(", df.sd$i.num, ")"))
    harrow(22, 40, "si.flow", df.mn$si.flow, dir = "right")
    if (vital == TRUE) {
      varrow(22, 40, "ds.flow", df.mn$ds.flow, dir = "out")
      varrow(57, 40, "di.flow", df.mn$di.flow, dir = "out")
      varrow(22, 40, "a.flow", df.mn$a.flow, dir = "in")
    }
  }

  ## 2. SIR Model
  if (dis.type == "SIR" && groups == 1) {
    mbox(5, 40, "Susceptible", paste0(df.mn$s.num, "(", df.sd$s.num, ")"))
    mbox(40, 40, "Infected", paste0(df.mn$i.num, "(", df.sd$i.num, ")"))
    mbox(75, 40, "Recovered", paste0(df.mn$r.num, "(", df.sd$r.num, ")"))
    harrow(5, 40, "si.flow", df.mn$si.flow, dir = "right")
    harrow(40, 40, "ir.flow", df.mn$ir.flow, dir = "right")
    if (vital == TRUE) {
      varrow(5, 40, "ds.flow", df.mn$ds.flow, dir = "out")
      varrow(40, 40, "di.flow", df.mn$di.flow, dir = "out")
      varrow(75, 40, "dr.flow", df.mn$dr.flow, dir = "out")
      varrow(5, 40, "a.flow", df.mn$a.flow, dir = "in")
    }
  }

  ## 3. SIS Model
  if (dis.type == "SIS" && groups == 1) {
    mbox(22, 40, "Susceptible", paste0(df.mn$s.num, "(", df.sd$s.num, ")"))
    mbox(57, 40, "Infected", paste0(df.mn$i.num, "(", df.sd$i.num, ")"))
    harrow(22, 40, "si.flow", df.mn$si.flow, dir = "right")
    harrow(22, 40, "is.flow", df.mn$is.flow, dir = "left")
    if (vital == TRUE) {
      varrow(22, 40, "ds.flow", df.mn$ds.flow, dir = "out")
      varrow(57, 40, "di.flow", df.mn$di.flow, dir = "out")
      varrow(22, 40, "a.flow", df.mn$a.flow, dir = "in")
    }
  }

  # Reset graphical parameters
  on.exit(par(ops))
}

#' @param run Model run number, for \code{dcm} class models with multiple runs
#'        (sensitivity analyses).
#' @method comp_plot dcm
#' @rdname comp_plot
#' @export
comp_plot.dcm <- function(x, at = 1, digits = 3, run = 1, ...) {


  ## Variables
  nsteps <- x$control$nsteps
  dis.type <- x$control$type
  groups <- x$param$groups
  vital <- x$param$vital

  ## Errors
  if (groups != 1) {
    stop("Only 1-group dcm models currently supported",
         call. = FALSE)
  }

  ## Time
  if (at > nsteps || at < 1) {
    stop("Specify a time step between 1 and ", nsteps)
  }
  intime <- at
  at <- which(x$control$timesteps == intime)

  ## Dataframe subsets
  df <- as.data.frame(x, run = run)
  df <- round(df[at, ], digits)

  ## Change graphical parameters
  ops <- list(mar = par()$mar, mfrow = par()$mfrow, mgp = par()$mgp)
  par(mar = c(0, 0, 2, 0))
  options(scipen = 10)

  ## Main Plot
  plot(0:100, 0:100, type = "n", axes = FALSE)
  title(main = paste(dis.type, "Model Diagram"))
  mtext(paste0("time=", intime, "  |  run=", run),
        side = 3, cex = 0.8, line = -1)

  ## 1. SI Model
  if (dis.type == "SI") {
    mbox(22, 40, "Susceptible", df$s.num)
    mbox(57, 40, "Infected", df$i.num)
    harrow(22, 40, "si.flow", df$si.flow, dir = "right")
    if (vital == TRUE) {
      varrow(22, 40, "ds.flow", df$ds.flow, dir = "out")
      varrow(57, 40, "di.flow", df$di.flow, dir = "out")
      varrow(22, 40, "a.flow", df$a.flow, dir = "in")
    }
  }

  ## 2. SIR Model
  if (dis.type == "SIR") {
    mbox(5, 40, "Susceptible", df$s.num)
    mbox(40, 40, "Infected", df$i.num)
    mbox(75, 40, "Recovered", df$r.num)
    harrow(5, 40, "si.flow", df$si.flow, dir = "right")
    harrow(40, 40, "ir.flow", df$ir.flow, dir = "right")
    if (vital == TRUE) {
      varrow(5, 40, "ds.flow", df$ds.flow, dir = "out")
      varrow(40, 40, "di.flow", df$di.flow, dir = "out")
      varrow(75, 40, "dr.flow", df$dr.flow, dir = "out")
      varrow(5, 40, "a.flow", df$a.flow, dir = "in")
    }
  }

  ## 3. SIS Model
  if (dis.type == "SIS") {
    mbox(22, 40, "Susceptible", df$s.num)
    mbox(57, 40, "Infected", df$i.num)
    harrow(22, 40, "si.flow", df$si.flow, dir = "right")
    harrow(22, 40, "is.flow", df$is.flow, dir = "left")
    if (vital == TRUE) {
      varrow(22, 40, "ds.flow", df$ds.flow, dir = "out")
      varrow(57, 40, "di.flow", df$di.flow, dir = "out")
      varrow(22, 40, "a.flow", df$a.flow, dir = "in")
    }
  }

  # Reset graphical parameters
  on.exit(par(ops))
}

## comp_plot helper utilities ##
#  Text box
mbox <- function(x, y, title, val) {
  polygon(c(x, x + 20, x + 20, x), c(y, y, y + 20, y + 20))
  text(x + 10, y + 10, paste(title, "\n n=", val, sep = ""), cex = 0.9)
}
#  Horizontal arrow
harrow <- function(xbox, ybox, title, val, dir) {
  if (dir == "right") {
    arrows(xbox + 20, ybox + 12, xbox + 35, lwd = 2, length = 0.15)
    text(xbox + 27.5, ybox + 17, paste(title, val, sep = "="), cex = 0.8)
  }
  if (dir == "left") {
    arrows(xbox + 20 + 15, ybox + 5, xbox + 20, lwd = 2, length = 0.15)
    text(xbox + 27.5, ybox + 2, paste(title, val, sep = "="), cex = 0.8)
  }
}
#  Vertical arrow
varrow <- function(xbox, ybox, title, val, dir) {
  if (dir == "out") {
    arrows(xbox + 10, ybox, xbox + 10, ybox - 25, lwd = 2, length = 0.15)
    text(xbox + 10, ybox - 12.5, paste(title, val, sep = "="),
         cex = 0.8, pos = 4)
  }
  if (dir == "in") {
    arrows(xbox + 10, ybox + 45, xbox + 10, ybox + 20, lwd = 2, length = 0.15)
    text(xbox + 10, ybox + 32.5, paste(title, val, sep = "="),
         cex = 0.8, pos = 4)
  }
}
