
verbose.dcm <- function(x, type, s) {

  if (type == "startup") {
    if (x$verbose == TRUE & x$nruns > 1) {
      cat("===============================")
      cat("\nStarting", x$nruns, "dcm Simulations")
      cat("\n===============================\n")
    }
  }

  if (type == "progress") {
    if (x$verbose == TRUE & x$nruns > 1) {
      cat("Run = ", s, "/", x$nruns, "\n", sep="")
    }
  }

}

#' @title Progress Print Module for Stochastic Individual Contact Models
#'
#' @description This function prints progress from stochastic individual contact
#'              models simulated with \code{icm} to the console.
#'
#' @param x if the \code{type} is "startup", then an object of class
#'        \code{control.icm}, otherwise the all master data object in \code{icm}
#'        simulations.
#' @param type progress type, either of "startup" for starting messages before
#'        all simulations, or "progress" for time step specific messages.
#' @param s current simulation number, if type is "progress".
#' @param at current time step, if type is "progress".
#'
#' @export
#' @keywords internal
#'
verbose.icm <- function(x, type, s, at) {

  if (type == "startup") {
    if (x$verbose == TRUE) {
      if (x$nsims == 1) {
        cat("===============================")
        cat("\nStarting 1 icm Simulation")
        cat("\n===============================\n")
      } else {
        cat("===============================")
        cat("\nStarting", x$nsims, "icm Simulations")
        cat("\n===============================")
      }
    }
  }

  if (type == "progress") {
    if (x$control$verbose == TRUE) {
      if (x$control$verbose.int == 0 && at == x$control$nsteps) {
        cat("\nSim = ", s, "/", x$control$nsims, sep="")
      }
      if (x$control$verbose.int > 0 && (at %% x$control$verbose.int == 0)) {
        cat("\014")
        cat("\nNetwork Disease Simulation")
        cat("\n----------------------------")
        cat("\nSimulation: ", s, "/", x$control$nsims, sep="")
        cat("\nTimestep: ", at, "/", x$control$nsteps, sep="")
        if (x$param$groups == 1) {
          cat("\nIncidence:", x$epi$si.flow[at])
        }
        if (x$param$groups == 2) {
          cat("\nIncidence:", x$epi$si.flow[at] + x$epi$si.flow.g2[at])
        }
        if (x$control$type == "SIR") {
          if (x$param$groups == 1) {
            cat("\nRecoveries:", x$epi$ir.flow[at])
          }
          if (x$param$groups == 2) {
            cat("\nRecoveries:", x$epi$ir.flow[at] +
                                 x$epi$ir.flow.g2[at])
          }
        }
        if (x$control$type == "SIS") {
          if (x$param$groups == 1) {
            cat("\nRecoveries:", x$epi$is.flow[at])
          }
          if (x$param$groups == 2) {
            cat("\nRecoveries:", x$epi$is.flow[at] +
                                 x$epi$is.flow.g2[at])
          }
        }
        if (x$param$groups == 1) {
          cat("\nPrevalence:", x$epi$i.num[at])
        }
        if (x$param$groups == 2) {
          cat("\nPrevalence:", x$epi$i.num[at] + x$epi$i.num.g2[at])
        }
        if (x$control$type %in% c("SI", "SIS")) {
          if (x$param$groups == 1) {
            cat("\nPopulation:", x$epi$s.num[at] + x$epi$i.num[at])
          }
          if (x$param$groups == 2) {
            cat("\nPopulation:", x$epi$s.num[at] + x$epi$s.num.g2[at]+
                                 x$epi$i.num[at] + x$epi$i.num.g2[at])
          }
        }
        if (x$control$type == "SIR") {
          if (x$param$groups == 1) {
            cat("\nPopulation:", x$epi$s.num[at] +
                                 x$epi$i.num[at] +
                                 x$epi$r.num[at])
          }
          if (x$param$groups == 2) {
            cat("\nPopulation:", x$epi$s.num[at] +
                                 x$epi$i.num[at] +
                                 x$epi$r.num[at] +
                                 x$epi$s.num.g2[at] +
                                 x$epi$i.num.g2[at] +
                                 x$epi$r.num.g2[at])
          }
        }
        if (x$param$vital == TRUE) {
          if (x$param$groups == 1) {
            cat("\nBirths:", x$epi$b.flow[at])
            cat("\nDeaths, susceptibles:", x$epi$ds.flow[at])
            cat("\nDeaths, infecteds:", x$epi$di.flow[at])
            if (x$control$type == "SIR") {
              cat("\nDeaths, recovered:", x$epi$dr.flow[at])
            }
          }
          if (x$param$groups == 2) {
            cat("\nBirths:", x$epi$b.flow[at] + x$epi$b.flow.g2[at])
            cat("\nDeaths, susceptible:", x$epi$ds.flow[at] +
                                          x$epi$ds.flow.g2[at])
            cat("\nDeaths, infected:", x$epi$di.flow[at] +
                                       x$epi$di.flow.g2[at])
            if (x$control$type == "SIR") {
              cat("\nDeaths, recovered:", x$epi$dr.flow[at] +
                                          x$epi$dr.flow.g2[at])
            }
          }
        }
        cat("\n----------------------------")
      }
    }
  }

}


#' @title Progress Print Module for Stochastic Network Models
#'
#' @description This function prints progress from stochastic network models
#'              simulated with \code{netsim} to the console.
#'
#' @param x if the \code{type} is "startup", then an object of class
#'        \code{control.net}, otherwise the all master data object in \code{netsim}
#'        simulations.
#' @param type progress type, either of "startup" for starting messages before
#'        all simulations, or "progress" for time step specific messages.
#' @param s current simulation number, if type is "progress"
#' @param at current time step, if type is "progress"
#'
#' @export
#' @keywords internal
#'
verbose.net <- function(x, type, s, at) {

  if (type == "startup") {
    if (x$verbose == TRUE) {
      if (x$nsims == 1) {
        cat("===================================")
        cat("\nStarting 1 Epidemic Simulation")
        cat("\n===================================")
      } else {
        cat("===================================")
        cat("\nStarting", x$nsims, "Epidemic Simulations")
        cat("\n===================================")
      }
    }
  }

  if (type == "progress") {
    if (x$control$verbose == TRUE) {
      if (x$control$verbose.int == 0 && at == x$control$nsteps) {
        cat("\nSim = ", s, "/", x$control$nsims, sep="")
      }
      if (x$control$verbose.int > 0 && (at %% x$control$verbose.int == 0)) {
        cat("\014")
        cat("\nNetwork Disease Simulation")
        cat("\n----------------------------")
        cat("\nSimulation: ", s, "/", x$control$nsims, sep="")
        cat("\nTimestep: ", at, "/", x$control$nsteps, sep="")
        if (x$param$modes == 1) {
          cat("\nIncidence:", x$epi$si.flow[at])
        }
        if (x$param$modes == 2) {
          cat("\nIncidence:", x$epi$si.flow[at] + x$epi$si.flow.m2[at])
        }
        if (x$control$type == "SIR") {
          if (x$param$modes == 1) {
            cat("\nRecoveries:", x$epi$ir.flow[at])
          }
          if (x$param$modes == 2) {
            cat("\nRecoveries:", x$epi$ir.flow[at] +
                                 x$epi$ir.flow.m2[at])
          }
        }
        if (x$control$type == "SIS") {
          if (x$param$modes == 1) {
            cat("\nRecoveries:", x$epi$is.flow[at])
          }
          if (x$param$modes == 2) {
            cat("\nRecoveries:", x$epi$is.flow[at] +
                                 x$epi$is.flow.m2[at])
          }
        }
        if (x$param$modes == 1) {
          cat("\nPrevalence:", x$epi$i.num[at])
        }
        if (x$param$modes == 2) {
          cat("\nPrevalence:", x$epi$i.num[at] + x$epi$i.num.m2[at])
        }
        if (x$control$type %in% c("SI", "SIS")) {
          if (x$param$modes == 1) {
            cat("\nPopulation:", x$epi$s.num[at] + x$epi$i.num[at])
          }
          if (x$param$modes == 2) {
            cat("\nPopulation:", x$epi$s.num[at] + x$epi$s.num.m2[at] +
                                 x$epi$i.num[at] + x$epi$i.num.m2[at])
          }
        }
        if (x$control$type == "SIR") {
          if (x$param$modes == 1) {
            cat("\nPopulation:", x$epi$s.num[at] +
                                 x$epi$i.num[at] +
                                 x$epi$r.num[at])
          }
          if (x$param$modes == 2) {
            cat("\nPopulation:", x$epi$s.num[at] +
                                 x$epi$i.num[at] +
                                 x$epi$r.num[at] +
                                 x$epi$s.num.m2[at] +
                                 x$epi$i.num.m2[at] +
                                 x$epi$r.num.m2[at])
          }
        }
        if (x$param$vital == TRUE) {
          if (x$param$modes == 1) {
            cat("\nBirths:", x$epi$b.flow[at])
            cat("\nDeaths, susceptibles:", x$epi$ds.flow[at])
            cat("\nDeaths, infecteds:", x$epi$di.flow[at])
            if (x$control$type == "SIR") {
              cat("\nDeaths, recovered:", x$epi$dr.flow[at])
            }
          }
          if (x$param$modes == 2) {
            cat("\nBirths:", x$epi$b.flow[at] + x$epi$b.flow.m2[at])
            cat("\nDeaths, susceptible:", x$epi$ds.flow[at] +
                                          x$epi$ds.flow.m2[at])
            cat("\nDeaths, infected:", x$epi$di.flow[at] +
                                       x$epi$di.flow.m2[at])
            if (x$control$type == "SIR") {
              cat("\nDeaths, recovered:", x$epi$dr.flow[at] +
                                          x$epi$dr.flow.m2[at])
            }
          }
        }
        cat("\n----------------------------")
      }
    }
  }

}

