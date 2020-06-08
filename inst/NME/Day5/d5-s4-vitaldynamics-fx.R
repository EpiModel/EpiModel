
##
## Simple SI Model with Vital Dynamics
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: August 2018
##


# New Aging Module --------------------------------------------------------

aging <- function(dat, at) {

  ## Initialization of Age
  if (at == 2) {
    # Pull age from the fitted network model
    dat$attr$age <- get.vertex.attribute(dat$nw, "age")
  }

  # Update age on attr and also the network
  dat$attr$age <- dat$attr$age + 1/52
  dat$nw <- set.vertex.attribute(dat$nw, "age", dat$attr$age)

  ## Summary statistics ##
  dat$epi$meanAge[at] <- mean(dat$attr$age, na.rm = TRUE)

  return(dat)
}


# Update Death Module -----------------------------------------------------

dfunc <- function(dat, at) {

  ## Attributes ##
  active <- dat$attr$active
  age <- dat$attr$age
  status <- dat$attr$status

  ## Parameters ##
  mort.rates <- dat$param$mortality.rates
  mort.dis.mult <- dat$param$mortality.disease.mult

  ## Query alive ##
  idsElig <- which(active == 1)
  nElig <- length(idsElig)
  nDeaths <- 0

  if (nElig > 0) {

    ## Calculate age-specific death rates for each eligible node ##
    ## Everyone older than 85 gets the final mortality
    whole_ages_of_elig <- pmin(ceiling(age[idsElig]), 86)
    death_rates_of_elig <- mort.rates[whole_ages_of_elig]

    ## Multiply death rates for diseased persons
    idsElig.inf <- which(status[idsElig] == "i")
    death_rates_of_elig[idsElig.inf] <- death_rates_of_elig[idsElig.inf] * mort.dis.mult

    ## Simulate mortality process
    vecDeaths <- which(rbinom(nElig, 1, death_rates_of_elig) == 1)
    idsDeaths <- idsElig[vecDeaths]
    nDeaths <- length(idsDeaths)

    ## Update nodal attributes on attr and networkDynamic object ##
    if (nDeaths > 0) {
      dat$attr$active[idsDeaths] <- 0
      dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                    v = idsDeaths, deactivate.edges = TRUE)
    }
  }

  ## Summary statistics ##
  dat$epi$d.flow[at] <- nDeaths

  return(dat)
}


# Updated Birth Module ----------------------------------------------------

bfunc <- function(dat, at) {

  ## Parameters ##
  n <- network.size(dat$nw)
  b.rate <- dat$param$birth.rate

  ## Process ##
  nBirthsExp <- n * b.rate
  nBirths <- rpois(1, nBirthsExp)

  if (nBirths > 0) {
    dat$nw <- add.vertices(dat$nw, nv = nBirths)
    newNodes <- (n + 1):(n + nBirths)
    dat$nw <- activate.vertices(dat$nw, onset = at, terminus = Inf, v = newNodes)
  }

  # Update attributes
  if (nBirths > 0) {
    dat$attr$active <- c(dat$attr$active, rep(1, nBirths))
    dat$attr$status <- c(dat$attr$status, rep("s", nBirths))
    dat$attr$infTime <- c(dat$attr$infTime, rep(NA, nBirths))
    dat$attr$entrTime <- c(dat$attr$entrTime, rep(at, nBirths))
    dat$attr$exitTime <- c(dat$attr$exitTime, rep(NA, nBirths))

    # Updated age must go on both attr list and network b/c it's in the ERGM
    dat$attr$age <- c(dat$attr$age, rep(0, nBirths))
    dat$nw <- set.vertex.attribute(dat$nw, "age", 0, newNodes)
  }

  ## Summary statistics ##
  dat$epi$b.flow[at] <- nBirths

  return(dat)
}
