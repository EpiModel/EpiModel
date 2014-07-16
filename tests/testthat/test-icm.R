context("icm")

test_that("Act rate balance works when specified to g2", {
  param <- param.icm(trans.rate = 0.2,
                     trans.rate.g2 = 0.1,
                     act.rate.g2 = 0.5,
                     balance = "g2",
                     rec.rate = 1/50,
                     rec.rate.g2 = 1/50,
                     b.rate = 1/100,
                     b.rate.g2 = NA,
                     ds.rate = 1/100,
                     ds.rate.g2 = 1/100,
                     di.rate = 1/90,
                     di.rate.g2 = 1/90)
  init <- init.icm(s.num = 50,
                   i.num = 10,
                   s.num.g2 = 100,
                   i.num.g2 = 1)
  control <- control.icm(type = "SIS",
                         nsteps = 25,
                         nsims = 1,
                         verbose = FALSE)
  x <- icm(param, init, control)
  expect_is(x, "icm")
})
