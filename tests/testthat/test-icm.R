context("ICM standard models")

# SI Models ---------------------------------------------------------------

test_that("SI, 1G, CL: 1 sim", {
  param <- param.icm(inf.prob = 0.2, act.rate = 0.25)
  init <- init.icm(s.num = 500, i.num = 1)
  control <- control.icm(type = "SI", nsteps = 5, nsims = 1, verbose = FALSE)
  x <- icm(param, init, control)
  expect_is(x, "icm")
})

test_that("SI, 1G, CL: 2 sims", {
  param <- param.icm(inf.prob = 0.2, act.rate = 0.25)
  init <- init.icm(s.num = 500, i.num = 1)
  control <- control.icm(type = "SI", nsteps = 5, nsims = 2, verbose = FALSE)
  x <- icm(param, init, control)
  expect_is(x, "icm")
})

test_that("SI, 1G, OP: 1 sim", {
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25, a.rate = 1 / 100,
                     ds.rate = 1 / 100, di.rate = 1 / 90)
  init <- init.icm(s.num = 500, i.num = 1)
  control <- control.icm(type = "SI", nsteps = 5, nsims = 1, verbose = FALSE)
  x <- icm(param, init, control)
  expect_is(x, "icm")
})

test_that("SI, 2G, OP: 2 sims", {
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25, inf.prob.g2 = 0.1, balance = "g1",
                     a.rate = 1 / 100, a.rate.g2 = NA, ds.rate = 1 / 100,
                     ds.rate.g2 = 1 / 100, di.rate = 1 / 90,
                     di.rate.g2 = 1 / 90)
  init <- init.icm(s.num = 500, i.num = 1, s.num.g2 = 500, i.num.g2 = 0)
  control <- control.icm(type = "SI", nsteps = 5, nsims = 2, verbose = FALSE)
  x <- icm(param, init, control)
  expect_is(x, "icm")
})


# SIR Models --------------------------------------------------------------

test_that("SIR, 1G, CL: 1 sim", {
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25, rec.rate = 1 / 50)
  init <- init.icm(s.num = 500, i.num = 1, r.num = 0)
  control <- control.icm(type = "SIR", nsteps = 5, nsims = 1, verbose = FALSE)
  x <- icm(param, init, control)
  expect_is(x, "icm")
})

test_that("SIR, 1G, OP: 1 sim", {
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 2, rec.rate = 1 / 50, a.rate = 1 / 100,
                     ds.rate = 1 / 100, di.rate = 1 / 90, dr.rate = 1 / 100)
  init <- init.icm(s.num = 500, i.num = 1, r.num = 0)
  control <- control.icm(type = "SIR", nsteps = 5, nsims = 1, verbose = FALSE)
  x <- icm(param, init, control)
  expect_is(x, "icm")
})

test_that("SIR, 2G, OP: 2 sims", {
  param <- param.icm(inf.prob = 0.2, inf.prob.g2 = 0.1,
                     act.rate = 1, balance = "g1",
                     rec.rate = 1 / 50, rec.rate.g2 = 1 / 50,
                     a.rate = 1 / 100, a.rate.g2 = NA,
                     ds.rate = 1 / 100, ds.rate.g2 = 1 / 100,
                     di.rate = 1 / 90, di.rate.g2 = 1 / 90,
                     dr.rate = 1 / 100, dr.rate.g2 = 1 / 100)
  init <- init.icm(s.num = 500, i.num = 1, r.num = 0,
                   s.num.g2 = 500, i.num.g2 = 1, r.num.g2 = 0)
  control <- control.icm(type = "SIR", nsteps = 5, nsims = 2, verbose = FALSE)
  x <- icm(param, init, control)
  expect_is(x, "icm")
})



# SIS Models --------------------------------------------------------------

test_that("SIS, 1G, CL: 1 sim", {
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25, rec.rate = 1 / 50)
  init <- init.icm(s.num = 500, i.num = 1)
  control <- control.icm(type = "SIS", nsteps = 5, nsims = 1, verbose = FALSE)
  x <- icm(param, init, control)
  expect_is(x, "icm")
})

test_that("SIS, 2G, CL: 2 sims", {
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25, inf.prob.g2 = 0.1, balance = "g1",
                     rec.rate = 1 / 100, rec.rate.g2 = 1 / 100)
  init <- init.icm(s.num = 500, i.num = 1,
                   s.num.g2 = 500, i.num.g2 = 0)
  control <- control.icm(type = "SIS", nsteps = 5, nsims = 2, verbose = FALSE)
  x <- icm(param, init, control)
  expect_is(x, "icm")
})

test_that("SIS, 1G, OP: 1 sim", {
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.5, rec.rate = 1 / 50,
                     a.rate = 1 / 100, ds.rate = 1 / 100, di.rate = 1 / 90)
  init <- init.icm(s.num = 500, i.num = 1)
  control <- control.icm(type = "SIS", nsteps = 10, nsims = 1, verbose = FALSE)
  x <- icm(param, init, control)
  expect_is(x, "icm")
})



# Other -------------------------------------------------------------------

test_that("Act rate balance works when specified to g2", {
  param <- param.icm(inf.prob = 0.2,
                     inf.prob.g2 = 0.1,
                     act.rate.g2 = 0.5,
                     balance = "g2",
                     rec.rate = 1 / 50,
                     rec.rate.g2 = 1 / 50,
                     a.rate = 1 / 100,
                     a.rate.g2 = NA,
                     ds.rate = 1 / 100,
                     ds.rate.g2 = 1 / 100,
                     di.rate = 1 / 90,
                     di.rate.g2 = 1 / 90)
  init <- init.icm(s.num = 50,
                   i.num = 10,
                   s.num.g2 = 100,
                   i.num.g2 = 1)
  control <- control.icm(type = "SIS",
                         nsteps = 10,
                         nsims = 1,
                         verbose = FALSE)
  x <- icm(param, init, control)
  expect_is(x, "icm")
})

context("New ICM models")

test_that("ICM new modules", {

  ## New Aging Module
  aging <- function(dat, at) {

    if (at == 2) {
      n <- sum(dat$attr$active == 1)
      dat$attr$age <- sample(18:49, n, replace = TRUE)
    } else {
      dat$attr$age <- dat$attr$age + 1 / 12
    }

    return(dat)
  }

  ## Replacement Departure Module
  dfunc <- function(dat, at) {

    idsElig <- which(dat$attr$active == 1)
    nElig <- length(idsElig)
    if (nElig > 0) {
      ages <- dat$attr$age[idsElig]
      life.expt <- dat$param$life.expt
      departure.rates <- pmin(1, 1 / (life.expt * 12 - ages * 12))
      vecDepartures <- which(rbinom(nElig, 1, departure.rates) == 1)
      idsDepartures <- idsElig[vecDepartures]
      nDepartures <- length(idsDepartures)
      if (nDepartures > 0) {
        dat$attr$active[idsDepartures] <- 0
      }
    } else {
      nDepartures <- 0
    }

    if (at == 2) {
      dat$epi$d.flow <- c(0, nDepartures)
    } else {
      dat$epi$d.flow[at] <- nDepartures
    }

    return(dat)
  }


  ## Replacement Arrival Module
  afunc <- function(dat, at) {

    growth.rate <- dat$param$growth.rate
    exptPopSize <- dat$epi$num[1] * (1 + growth.rate * at)

    numNeeded <- exptPopSize - sum(dat$attr$active == 1)
    if (numNeeded > 0) {
      nArrivals <- rpois(1, numNeeded)
    } else {
      nArrivals <- 0
    }

    dat$attr$active <- c(dat$attr$active, rep(1, nArrivals))
    dat$attr$group <- c(dat$attr$group, rep(1, nArrivals))
    dat$attr$status <- c(dat$attr$status, rep("s", nArrivals))
    dat$attr$infTime <- c(dat$attr$infTime, rep(NA, nArrivals))
    dat$attr$age <- c(dat$attr$age, rep(18, nArrivals))

    if (at == 2) {
      dat$epi$a.flow <- c(0, nArrivals)
    } else {
      dat$epi$a.flow[at] <- nArrivals
    }

    return(dat)
  }


  param <- param.icm(inf.prob = 0.15,
                     act.rate = 0.15,
                     growth.rate = 0.01 / 12,
                     life.expt = 70)
  init <- init.icm(s.num = 1000, i.num = 100)
  control <- control.icm(type = "SI", nsims = 2, nsteps = 250,
                         departures.FUN = dfunc, arrivals.FUN = afunc,
                         aging.FUN = aging,
                         verbose = FALSE)



  mod <- icm(param, init, control)
  expect_is(mod, "icm")

})
