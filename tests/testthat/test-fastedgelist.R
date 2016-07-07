

# tests for evaluating the fast.edgelist mode and ensuring that its output matches
context('fast.edgelist mode')

# make sure the mode switch works
test_that('mode switch works',{
  nw <- network.initialize(n = 100, directed = FALSE)
  formation <- ~edges
  target.stats <- 50
  
  # Example 2: Dependent SIR Model
  # Recalculate dissolution coefficient with death rate
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20,
                                 d.rate = 0.0021)
  
  # Reestimate the model with new coefficient
  set.seed(1)
  est2 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
  
  # Reset parameters to include demographic rates
  param <- param.net(inf.prob = 0.3, 
                     rec.rate = 0.02, 
                     b.rate = 0.002, 
                     ds.rate = 0.001, 
                     di.rate = 0.001, 
                     dr.rate = 0.001)
  init<- init.net(i.num = 10,
                   r.num = 0)
  
  control_old <- control.net(type = "SIR", nsteps = 5, nsims = 1,
                             tea.status = FALSE,
                             save.network=TRUE,
                             use.pids = FALSE,
                             fast.edgelist=FALSE,
                             verbose=FALSE)
  set.seed(1)
  simold <- netsim(est2, param, init, control_old)
  
  control_new <- control.net(type = "SIR", nsteps = 5, nsims = 1,
                             tea.status = FALSE,
                             save.network=TRUE,
                             use.pids = FALSE,
                             fast.edgelist=TRUE,
                             verbose=FALSE)
  set.seed(1)
  simnew <- netsim(est2, param, init, control_new)
  

  
  # since they were both run with the same seed, expect the transmats to be identical
  expect_equal(simold$stats$transmat,simnew$stats$transmat)
})

# TODO: add tests that non-supported models flagged correctly


# ideally, if we set the seed the same, we should get the same set of toggles via either simulate 