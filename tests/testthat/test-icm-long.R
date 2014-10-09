context("ICM extended models")

# SI Models ---------------------------------------------------------------

test_that("SI, 1G, CL: 1 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25)
  init <- init.icm(s.num = 500,
                   i.num = 1)
  control <- control.icm(type = "SI",
                         nsteps = 500,
                         nsims = 1,
                         verbose = FALSE)
  x <- icm(param, init, control)
  expect_equal(max(x$epi$i.num), 501)
  expect_equal(max(x$epi$s.num+x$epi$i.num), 501)
  expect_is(as.data.frame(x), "data.frame")
  test_icm(x)
})

test_that("SI, 1G, CL: 5 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25)
  init <- init.icm(s.num = 500,
                   i.num = 1)
  control <- control.icm(type = "SI",
                         nsteps = 500,
                         nsims = 5,
                         verbose = FALSE)
  x <- icm(param, init, control)
  expect_equal(max(x$epi$i.num), 501)
  expect_equal(max(x$epi$s.num+x$epi$i.num), 501)
  test_icm(x)
})

test_that("SI, 1G, CL: inf.prob = 0, status.rand = FALSE", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0,
                     act.rate = 0.25)
  init <- init.icm(s.num = 500,
                   i.num = 1,
                   status.rand = FALSE)
  control <- control.icm(type = "SI",
                         nsteps = 25,
                         nsims = 25,
                         verbose = FALSE)
  x <- icm(param, init, control)
  expect_equal(max(x$epi$i.num), 1)
  expect_equal(min(x$epi$i.num), 1)
  test_icm(x)
})

test_that("SI, 2G, CL: 1 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.5,
                     act.rate = 0.25,
                     inf.prob.g2 = 0.1,
                     balance = "g1")
  init <- init.icm(s.num = 500,
                   i.num = 1,
                   s.num.g2 = 500,
                   i.num.g2 = 0)
  control <- control.icm(type = "SI",
                         nsteps = 500,
                         nsims = 1,
                         verbose = FALSE)
  x <- icm(param, init, control)
  expect_equal(min(x$epi$i.num), 1)
  expect_equal(max(x$epi$num), 501)
  test_icm(x)
})

test_that("SI, 2G, CL: 5 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25,
                     inf.prob.g2 = 0.1,
                     balance = "g1")
  init <- init.icm(s.num = 500,
                   i.num = 1,
                   s.num.g2 = 500,
                   i.num.g2 = 0)
  control <- control.icm(type = "SI",
                         nsteps = 500,
                         nsims = 5,
                         verbose = FALSE)
  x <- icm(param, init, control)
  expect_equal(max(x$epi$i.num), 501)
  expect_equal(max(x$epi$s.num+x$epi$i.num), 501)
  test_icm(x)
})

test_that("SI, 1G, OP: 1 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25,
                     b.rate = 1/100,
                     ds.rate = 1/100,
                     di.rate = 1/90)
  init <- init.icm(s.num = 500,
                   i.num = 1)
  control <- control.icm(type = "SI",
                         nsteps = 500,
                         nsims = 1,
                         verbose = FALSE)
  x <- icm(param, init, control)
  test_icm(x)
})

test_that("SI, 1G, OP: 5 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25,
                     b.rate = 1/100,
                     ds.rate = 1/100,
                     di.rate = 1/90)
  init <- init.icm(s.num = 500,
                   i.num = 1)
  control <- control.icm(type = "SI",
                         nsteps = 500,
                         nsims = 5,
                         verbose = FALSE)
  x <- icm(param, init, control)
  test_icm(x)
})

test_that("SI, 2G, OP: 1 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25,
                     inf.prob.g2 = 0.1,
                     balance = "g1",
                     b.rate = 1/100,
                     b.rate.g2 = NA,
                     ds.rate = 1/100,
                     ds.rate.g2 = 1/100,
                     di.rate = 1/90,
                     di.rate.g2 = 1/90)
  init <- init.icm(s.num = 500,
                   i.num = 1,
                   s.num.g2 = 500,
                   i.num.g2 = 0)
  control <- control.icm(type = "SI",
                         nsteps = 500,
                         nsims = 1,
                         verbose = FALSE)
  x <- icm(param, init, control)
  test_icm(x)
})

test_that("SI, 2G, OP: 5 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25,
                     inf.prob.g2 = 0.1,
                     balance = "g1",
                     b.rate = 1/100,
                     b.rate.g2 = NA,
                     ds.rate = 1/100,
                     ds.rate.g2 = 1/100,
                     di.rate = 1/90,
                     di.rate.g2 = 1/90)
  init <- init.icm(s.num = 500,
                   i.num = 1,
                   s.num.g2 = 500,
                   i.num.g2 = 0)
  control <- control.icm(type = "SI",
                         nsteps = 500,
                         nsims = 5,
                         verbose = FALSE)
  x <- icm(param, init, control)
  test_icm(x)
})


# SIR Models --------------------------------------------------------------

test_that("SIR, 1G, CL: 1 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25,
                     rec.rate = 1/50)
  init <- init.icm(s.num = 500,
                   i.num = 1,
                   r.num = 0)
  control <- control.icm(type = "SIR",
                         nsteps = 500,
                         nsims = 1,
                         verbose = FALSE)
  x <- icm(param, init, control)
  expect_equal(max(x$epi$s.num+x$epi$i.num+x$epi$r.num), 501)
  test_icm(x)
})

test_that("SIR, 1G, CL: 5 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25,
                     rec.rate = 1/50)
  init <- init.icm(s.num = 500,
                   i.num = 1,
                   r.num = 0)
  control <- control.icm(type = "SIR",
                         nsteps = 500,
                         nsims = 5,
                         verbose = FALSE)
  x <- icm(param, init, control)
  expect_equal(max(x$epi$s.num+x$epi$i.num+x$epi$r.num), 501)
  test_icm(x)
})

test_that("SIR, 1G, CL: inf.prob = 0, status.rand = FALSE", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0,
                     act.rate = 0.25,
                     rec.rate = 0)
  init <- init.icm(s.num = 500,
                   i.num = 1,
                   r.num = 0,
                   status.rand = FALSE)
  control <- control.icm(type = "SIR",
                         nsteps = 25,
                         nsims = 25,
                         verbose = FALSE)
  x <- icm(param, init, control)
  expect_equal(max(x$epi$i.num), 1)
  expect_equal(min(x$epi$i.num), 1)
  test_icm(x)
})

test_that("SIR, 2G, CL: 1 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25,
                     inf.prob.g2 = 0.1,
                     balance = "g1",
                     rec.rate = 1/100,
                     rec.rate.g2 = 1/100)
  init <- init.icm(s.num = 500, i.num = 1, r.num = 0,
                   s.num.g2 = 500, i.num.g2 = 0, r.num.g2 = 0)
  control <- control.icm(type = "SIR",
                         nsims = 1,
                         nsteps = 1000,
                         verbose = FALSE)
  x <- icm(param, init, control)
  expect_equal(max(x$epi$s.num+x$epi$i.num+x$epi$r.num+
                     x$epi$s.num.g2+x$epi$i.num.g2+x$epi$r.num.g2), 1001)
  test_icm(x)
})

test_that("SIR, 2G, CL: 5 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25,
                     inf.prob.g2 = 0.1,
                     balance = "g1",
                     rec.rate = 1/100,
                     rec.rate.g2 = 1/100)
  init <- init.icm(s.num = 500, i.num = 1, r.num = 0,
                   s.num.g2 = 500, i.num.g2 = 0, r.num.g2 = 0)
  control <- control.icm(type = "SIR",
                         nsims = 5,
                         nsteps = 1000,
                         verbose = FALSE)
  x <- icm(param, init, control)
  expect_equal(max(x$epi$s.num+x$epi$i.num+x$epi$r.num+
                     x$epi$s.num.g2+x$epi$i.num.g2+x$epi$r.num.g2), 1001)
  test_icm(x)
})

test_that("SIR, 1G, OP: 1 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 3,
                     rec.rate = 1/50,
                     b.rate = 1/100,
                     ds.rate = 1/100,
                     di.rate = 1/90,
                     dr.rate = 1/100)
  init <- init.icm(s.num = 500,
                   i.num = 1,
                   r.num = 0)
  control <- control.icm(type = "SIR",
                         nsteps = 500,
                         nsims = 1,
                         verbose = FALSE)
  x <- icm(param, init, control)
  test_icm(x)
})

test_that("SIR, 1G, OP: 5 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 3,
                     rec.rate = 1/50,
                     b.rate = 1/100,
                     ds.rate = 1/100,
                     di.rate = 1/90,
                     dr.rate = 1/100)
  init <- init.icm(s.num = 500,
                   i.num = 1,
                   r.num = 0)
  control <- control.icm(type = "SIR",
                         nsteps = 500,
                         nsims = 5,
                         verbose = FALSE)
  x <- icm(param, init, control)
  test_icm(x)
})

test_that("SIR, 2G, OP: 1 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.2,
                     inf.prob.g2 = 0.1,
                     act.rate = 3,
                     balance = "g1",
                     rec.rate = 1/50,
                     rec.rate.g2 = 1/50,
                     b.rate = 1/100,
                     b.rate.g2 = NA,
                     ds.rate = 1/100,
                     ds.rate.g2 = 1/100,
                     di.rate = 1/90,
                     di.rate.g2 = 1/90,
                     dr.rate = 1/100,
                     dr.rate.g2 = 1/100)
  init <- init.icm(s.num = 500, i.num = 1, r.num = 0,
                   s.num.g2 = 500, i.num.g2 = 0, r.num.g2 = 0)
  control <- control.icm(type = "SIR",
                         nsteps = 500,
                         nsims = 1,
                         verbose = FALSE)
  x <- icm(param, init, control)
  test_icm(x)
})

test_that("SIR, 2G, OP: 5 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.2,
                     inf.prob.g2 = 0.1,
                     act.rate = 3,
                     balance = "g1",
                     rec.rate = 1/50,
                     rec.rate.g2 = 1/50,
                     b.rate = 1/100,
                     b.rate.g2 = NA,
                     ds.rate = 1/100,
                     ds.rate.g2 = 1/100,
                     di.rate = 1/90,
                     di.rate.g2 = 1/90,
                     dr.rate = 1/100,
                     dr.rate.g2 = 1/100)
  init <- init.icm(s.num = 500, i.num = 1, r.num = 0,
                   s.num.g2 = 500, i.num.g2 = 0, r.num.g2 = 0)
  control <- control.icm(type = "SIR",
                         nsteps = 500,
                         nsims = 5,
                         verbose = FALSE)
  x <- icm(param, init, control)
  test_icm(x)
})


# SIS Models --------------------------------------------------------------

test_that("SIS, 1G, CL: 1 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25,
                     rec.rate = 1/50)
  init <- init.icm(s.num = 500,
                   i.num = 1)
  control <- control.icm(type = "SIS",
                         nsteps = 500,
                         nsims = 1,
                         verbose = FALSE)
  x <- icm(param, init, control)
  expect_equal(max(x$epi$s.num+x$epi$i.num), 501)
  test_icm(x)
})

test_that("SIS, 1G, CL: 5 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25,
                     rec.rate = 1/50)
  init <- init.icm(s.num = 500,
                   i.num = 1)
  control <- control.icm(type = "SIS",
                         nsteps = 500,
                         nsims = 5,
                         verbose = FALSE)
  x <- icm(param, init, control)
  expect_equal(max(x$epi$s.num+x$epi$i.num), 501)
  test_icm(x)
})

test_that("SIS, 1G, CL: inf.prob = 0, status.rand = FALSE", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0,
                     act.rate = 0.25,
                     rec.rate = 0)
  init <- init.icm(s.num = 500,
                   i.num = 1,
                   status.rand = FALSE)
  control <- control.icm(type = "SIS",
                         nsteps = 25,
                         nsims = 25,
                         verbose = FALSE)
  x <- icm(param, init, control)
  test_icm(x)
})

test_that("SIS, 2G, CL: 1 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25,
                     inf.prob.g2 = 0.1,
                     balance = "g1",
                     rec.rate = 1/100,
                     rec.rate.g2 = 1/100)
  init <- init.icm(s.num = 500,
                   i.num = 1,
                   s.num.g2 = 500,
                   i.num.g2 = 0)
  control <- control.icm(type = "SIS",
                         nsteps = 1000,
                         nsims = 1,
                         verbose = FALSE)
  x <- icm(param, init, control)
  expect_equal(max(x$epi$s.num+x$epi$i.num), 501)
  expect_equal(max(x$epi$s.num+x$epi$i.num+x$epi$s.num.g2+x$epi$i.num.g2), 1001)
  test_icm(x)
})

test_that("SIS, 2G, CL: 5 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25,
                     inf.prob.g2 = 0.1,
                     balance = "g1",
                     rec.rate = 1/100,
                     rec.rate.g2 = 1/100)
  init <- init.icm(s.num = 500,
                   i.num = 1,
                   s.num.g2 = 500,
                   i.num.g2 = 0)
  control <- control.icm(type = "SIS",
                         nsteps = 1000,
                         nsims = 5,
                         verbose = FALSE)
  x <- icm(param, init, control)
  expect_equal(max(x$epi$s.num+x$epi$i.num), 501)
  expect_equal(max(x$epi$s.num+x$epi$i.num+x$epi$s.num.g2+x$epi$i.num.g2), 1001)
  test_icm(x)
})

test_that("SIS, 1G, OP: 1 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.5,
                     rec.rate = 1/50,
                     b.rate = 1/100,
                     ds.rate = 1/100,
                     di.rate = 1/90)
  init <- init.icm(s.num = 500,
                   i.num = 1)
  control <- control.icm(type = "SIS",
                         nsteps = 500,
                         nsims = 1,
                         verbose = FALSE)
  x <- icm(param, init, control)
  test_icm(x)
})

test_that("SIS, 1G, OP: 5 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.5,
                     rec.rate = 1/50,
                     b.rate = 1/100,
                     ds.rate = 1/100,
                     di.rate = 1/90)
  init <- init.icm(s.num = 500,
                   i.num = 1)
  control <- control.icm(type = "SIS",
                         nsteps = 500,
                         nsims = 5,
                         verbose = FALSE)
  x <- icm(param, init, control)
  test_icm(x)
})

test_that("SIS, 2G, OP: 1 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.2,
                     inf.prob.g2 = 0.1,
                     act.rate = 0.5,
                     balance = "g1",
                     rec.rate = 1/50,
                     rec.rate.g2 = 1/50,
                     b.rate = 1/100,
                     b.rate.g2 = NA,
                     ds.rate = 1/100,
                     ds.rate.g2 = 1/100,
                     di.rate = 1/90,
                     di.rate.g2 = 1/90)
  init <- init.icm(s.num = 500,
                   i.num = 1,
                   s.num.g2 = 500,
                   i.num.g2 = 1)
  control <- control.icm(type = "SIS",
                         nsteps = 500,
                         nsims = 1,
                         verbose = FALSE)
  x <- icm(param, init, control)
  test_icm(x)
})

test_that("SIS, 2G, OP: 5 sim", {
  skip_on_cran()
  param <- param.icm(inf.prob = 0.2,
                     inf.prob.g2 = 0.1,
                     act.rate = 0.5,
                     balance = "g1",
                     rec.rate = 1/50,
                     rec.rate.g2 = 1/50,
                     b.rate = 1/100,
                     b.rate.g2 = NA,
                     ds.rate = 1/100,
                     ds.rate.g2 = 1/100,
                     di.rate = 1/90,
                     di.rate.g2 = 1/90)
  init <- init.icm(s.num = 500,
                   i.num = 1,
                   s.num.g2 = 500,
                   i.num.g2 = 1)
  control <- control.icm(type = "SIS",
                         nsteps = 500,
                         nsims = 5,
                         verbose = FALSE)
  x <- icm(param, init, control)
  test_icm(x)
})
