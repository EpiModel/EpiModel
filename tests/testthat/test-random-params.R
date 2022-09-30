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

  expect_warning(param <- param.net(
      inf.prob = 0.3,
      act.rate = 0.3,
      random.params = my_randoms)
  )
  expect_message(generate_random_params(param, verbose = TRUE))
  expect_silent(generate_random_params(param, verbose = FALSE))

  param <- param.net(inf.prob = 0.3, act.rate = 0.1)
  expect_equal(generate_random_params(param), param)

  param <- param.net(inf.prob = 0.3, random.params = list())
  expect_equal(generate_random_params(param), param)

  param <- param.net(inf.prob = 0.3, random.params = 4)
  expect_error(generate_random_params(param))

  param <- param.net(inf.prob = 0.3, random.params = list(1))
  expect_error(generate_random_params(param))


  generate_correlated_params <- function() {
    param.unique <- runif(1)
    param.set.1 <- param.unique + runif(2)
    param.set.2 <- param.unique * rnorm(3)

    return(list(param.unique, param.set.1, param.set.2))
  }

  # Data.frame set of random parameters :
  correlated_params <- t(replicate(10, unlist(generate_correlated_params())))
  correlated_params <- as.data.frame(correlated_params)
  colnames(correlated_params) <- c(
    "param.unique",
    "param.set.1_1", "param.set.1_2",
    "param.set.2_1", "param.set.2_2", "param.set.2_3"
  )

  randoms <- c(my_randoms, list(param.random.set = correlated_params))
  param <- param.net(inf.prob = 0.3, random.params = randoms)
  expect_silent(generate_random_params(param))

  # duplicated `act.rate` random definition
  colnames(correlated_params) <- c(
    "act.rate",
    "param.set.1_1", "param.set.1_2",
    "param.set.2_1", "param.set.2_2", "param.set.2_3"
  )
  randoms <- c(my_randoms, list(param.random.set = correlated_params))
  expect_warning(
    param <- param.net(inf.prob = 0.3, act.rate = 0.1, random.params = randoms)
  )
  expect_warning(generate_random_params(param))

  # malformed name "param_set.1_1"
  colnames(correlated_params) <- c(
    "act.rate",
    "param_set.1_1", "param.set.1_2",
    "param.set.2_1", "param.set.2_2", "param.set.2_3"
  )
  randoms <- c(my_randoms, list(param.random.set = correlated_params))
  param <- param.net(inf.prob = 0.3, random.params = randoms)
  expect_error(generate_random_params(param))

  # param.random.set not a data.frame
  randoms <- c(my_randoms, list(param.random.set = list()))
  param <- param.net(inf.prob = 0.3, random.params = randoms)
  expect_error(generate_random_params(param))
})

