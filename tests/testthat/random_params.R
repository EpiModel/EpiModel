context("Random parameters generators")

test_that("Random parameters generators", {

  my_randoms <- list(
    act.rate = param_random(c(0.25, 0.5, 0.75)),
    tx.halt.part.prob = function() rbeta(1, 1, 2),
    hiv.test.rate = function() c(
      rnorm(1, 0.015, 0.01),
      rnorm(1, 0.010, 0.01),
      rnorm(1, 0.020, 0.01)
    )
  )

  param <- param.net(inf.prob = 0.3, random.params = my_randoms)
  expect_message(generate_random_params(param, verbose = TRUE))
  expect_silent(generate_random_params(param, verbose = FALSE))

  param <- param.net(inf.prob = 0.3, acte.rate = 0.1)
  expect_equal(generate_random_params(param), param)

  param <- param.net(inf.prob = 0.3, random.params = list())
  expect_equal(generate_random_params(param), param)

  param <- param.net(inf.prob = 0.3, random.params = 4)
  expect_error(generate_random_params(param))

  param <- param.net(inf.prob = 0.3, random.params = list(1))
  expect_error(generate_random_params(param))
})
