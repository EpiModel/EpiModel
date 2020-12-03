# devtools::load_all()

#  my_randoms <- list(
#    act.rate = param_random(c(0.25, 0.5, 0.75)),
#    tx.halt.part.prob = function() rbeta(1, 1, 2),
#    hiv.test.rate = function() c(
#      rnorm(1, 0.015, 0.01),
#      rnorm(1, 0.010, 0.01),
#      rnorm(1, 0.020, 0.01)
#    )
#  )

# param <- param.net( inf.prob = 0.3, random.params = my_randoms)
# param <- generate_random_params(param, TRUE)
# param

# param <- param.net(inf.prob = 0.3, acte.rate = 0.1)
# param <- generate_random_params(param)
# param

# param <- param.net(inf.prob = 0.3, random.params = list())
# param <- generate_random_params(param)
# param

# param <- param.net(inf.prob = 0.3, random.params = 4)
# param <- generate_random_params(param)
# param

# param <- param.net(inf.prob = 0.3, random.params = list(1))
# param <- generate_random_params(param)
# param


