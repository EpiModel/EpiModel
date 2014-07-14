
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
          cat("\nIncidence:", x$out$si.flow[at])
        }
        if (x$param$groups == 2) {
          cat("\nIncidence:", x$out$si.flow[at] + x$out$si.flow.g2[at])
        }
        if (x$control$type == "SIR") {
          if (x$param$groups == 1) {
            cat("\nRecoveries:", x$out$ir.flow[at])
          }
          if (x$param$groups == 2) {
            cat("\nRecoveries:", x$out$ir.flow[at] +
                                 x$out$ir.flow.g2[at])
          }
        }
        if (x$control$type == "SIS") {
          if (x$param$groups == 1) {
            cat("\nRecoveries:", x$out$is.flow[at])
          }
          if (x$param$groups == 2) {
            cat("\nRecoveries:", x$out$is.flow[at] +
                                 x$out$is.flow.g2[at])
          }
        }
        if (x$param$groups == 1) {
          cat("\nPrevalence:", x$out$i.num[at])
        }
        if (x$param$groups == 2) {
          cat("\nPrevalence:", x$out$i.num[at] + x$out$i.num.g2[at])
        }
        if (x$control$type %in% c("SI", "SIS")) {
          if (x$param$groups == 1) {
            cat("\nPopulation:", x$out$s.num[at] + x$out$i.num[at])
          }
          if (x$param$groups == 2) {
            cat("\nPopulation:", x$out$s.num[at] + x$out$s.num.g2[at]+
                                 x$out$i.num[at] + x$out$i.num.g2[at])
          }
        }
        if (x$control$type == "SIR") {
          if (x$param$groups == 1) {
            cat("\nPopulation:", x$out$s.num[at] +
                                 x$out$i.num[at] +
                                 x$out$r.num[at])
          }
          if (x$param$groups == 2) {
            cat("\nPopulation:", x$out$s.num[at] +
                                 x$out$i.num[at] +
                                 x$out$r.num[at] +
                                 x$out$s.num.g2[at] +
                                 x$out$i.num.g2[at] +
                                 x$out$r.num.g2[at])
          }
        }
        if (x$param$vital == TRUE) {
          if (x$param$groups == 1) {
            cat("\nBirths:", x$out$b.flow[at])
            cat("\nDeaths, susceptibles:", x$out$ds.flow[at])
            cat("\nDeaths, infecteds:", x$out$di.flow[at])
            if (x$control$type == "SIR") {
              cat("\nDeaths, recovered:", x$out$dr.flow[at])
            }
          }
          if (x$param$groups == 2) {
            cat("\nBirths:", x$out$b.flow[at] + x$out$b.flow.g2[at])
            cat("\nDeaths, susceptible:", x$out$ds.flow[at] +
                                          x$out$ds.flow.g2[at])
            cat("\nDeaths, infected:", x$out$di.flow[at] +
                                       x$out$di.flow.g2[at])
            if (x$control$type == "SIR") {
              cat("\nDeaths, recovered:", x$out$dr.flow[at] +
                                          x$out$dr.flow.g2[at])
            }
          }
        }
        cat("\n----------------------------")
      }
    }
  }

}


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
          cat("\nIncidence:", x$out$si.flow[at])
        }
        if (x$param$modes == 2) {
          cat("\nIncidence:", x$out$si.flow[at] + x$out$si.flow.m2[at])
        }
        if (x$control$type == "SIR") {
          if (x$param$modes == 1) {
            cat("\nRecoveries:", x$out$ir.flow[at])
          }
          if (x$param$modes == 2) {
            cat("\nRecoveries:", x$out$ir.flow[at] +
                                 x$out$ir.flow.m2[at])
          }
        }
        if (x$control$type == "SIS") {
          if (x$param$modes == 1) {
            cat("\nRecoveries:", x$out$is.flow[at])
          }
          if (x$param$modes == 2) {
            cat("\nRecoveries:", x$out$is.flow[at] +
                                 x$out$is.flow.m2[at])
          }
        }
        if (x$param$modes == 1) {
          cat("\nPrevalence:", x$out$i.num[at])
        }
        if (x$param$modes == 2) {
          cat("\nPrevalence:", x$out$i.num[at] + x$out$i.num.m2[at])
        }
        if (x$control$type %in% c("SI", "SIS")) {
          if (x$param$modes == 1) {
            cat("\nPopulation:", x$out$s.num[at] + x$out$i.num[at])
          }
          if (x$param$modes == 2) {
            cat("\nPopulation:", x$out$s.num[at] + x$out$s.num.m2[at] +
                                 x$out$i.num[at] + x$out$i.num.m2[at])
          }
        }
        if (x$control$type == "SIR") {
          if (x$param$modes == 1) {
            cat("\nPopulation:", x$out$s.num[at] +
                                 x$out$i.num[at] +
                                 x$out$r.num[at])
          }
          if (x$param$modes == 2) {
            cat("\nPopulation:", x$out$s.num[at] +
                                 x$out$i.num[at] +
                                 x$out$r.num[at] +
                                 x$out$s.num.m2[at] +
                                 x$out$i.num.m2[at] +
                                 x$out$r.num.m2[at])
          }
        }
        if (x$param$vital == TRUE) {
          if (x$param$modes == 1) {
            cat("\nBirths:", x$out$b.flow[at])
            cat("\nDeaths, susceptibles:", x$out$ds.flow[at])
            cat("\nDeaths, infecteds:", x$out$di.flow[at])
            if (x$control$type == "SIR") {
              cat("\nDeaths, recovered:", x$out$dr.flow[at])
            }
          }
          if (x$param$modes == 2) {
            cat("\nBirths:", x$out$b.flow[at] + x$out$b.flow.m2[at])
            cat("\nDeaths, susceptible:", x$out$ds.flow[at] +
                                          x$out$ds.flow.m2[at])
            cat("\nDeaths, infected:", x$out$di.flow[at] +
                                       x$out$di.flow.m2[at])
            if (x$control$type == "SIR") {
              cat("\nDeaths, recovered:", x$out$dr.flow[at] +
                                          x$out$dr.flow.m2[at])
            }
          }
        }
        cat("\n----------------------------")
      }
    }
  }

}

