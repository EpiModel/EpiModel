context("DCM delayed differential equation models")

test_that("Delayed differntial equation models function", {

  lagsi <- function(t, t0, parms) {
    with(as.list(c(t0, parms)), {

      num <- s.num + i.num

      if (t < 2) {
        i.num.foi <- 0
      } else {
        i.num.foi <- lagvalue(t - 1, 2)
      }

      lambda <- inf.prob * act.rate * i.num.foi / num

      dS <- -lambda * s.num
      dI <- lambda * s.num

      list(c(dS, dI),
           num = num,
           si.flow = lambda * s.num,
           i.num.foi = i.num.foi)
    })
  }

  param <- param.dcm(inf.prob = 0.5, act.rate = 0.25)
  init <- init.dcm(s.num = 100, i.num = 10)
  control <- control.dcm(nsteps = 25, new.mod = lagsi, dede = TRUE)
  mod <- dcm(param, init, control)

  expect_is(mod, "dcm")
  expect_true(mod$control$dede)

})
