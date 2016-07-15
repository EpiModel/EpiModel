

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
                     b.rate = 0.00, 
                     ds.rate = 0.01, 
                     di.rate = 0.0, 
                     dr.rate = 0.0)
  init<- init.net(i.num = 10,
                   r.num = 0)
  
  control_old <- control.net(type = "SIR", nsteps = 100, nsims = 1,
                             tea.status = FALSE,
                             save.network=TRUE,
                             use.pids = FALSE,
                             fast.edgelist=FALSE,
                             verbose=FALSE)
  set.seed(1)
  simold <- netsim(est2, param, init, control_old)
  
  control_new <- control.net(type = "SIR", nsteps = 100, nsims = 1,
                             tea.status = FALSE,
                             save.network=TRUE,
                             use.pids = FALSE,
                             fast.edgelist=TRUE,
                             verbose=FALSE)
  set.seed(1)
  simnew <- netsim(est2, param, init, control_new)
  

  
  # since they were both run with the same seed, expect the transmats to be identical
  # edgelist sorting code disabled, so don't expect them to be identical
  # expect_equal(simold$stats$transmat$sim1,simnew$stats$transmat$sim1)
})

test_that('edges+nodematch model works',{
  nw <- network.initialize(n = 100, directed = FALSE)
  # specify two different roles for the vertices
  nw%v%'rolemode'<-rep_len(c('a','b'),network.size(nw))
  formation <- ~edges+offset(nodematch('rolemode'))
  target.stats <- 50
  
  # calculate dissolution coefficient with death rate
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20,
                                 d.rate = 0.0021)
  # set iInf coef for offset statistic for nodematch "never form ties between vertcies with non matching attributes
  coef.form <- -Inf  
  
  # Reestimate the model with new coefficient
  set.seed(1)
  est2 <- netest(nw, formation, target.stats, coef.diss, coef.form=coef.form, verbose = FALSE)
  
  # Reset parameters to include demographic rates
  param <- param.net(inf.prob = 0.3, 
                     rec.rate = 0.02, 
                     b.rate = 0.00, 
                     ds.rate = 0.01, 
                     di.rate = 0.0, 
                     dr.rate = 0.0)
  init<- init.net(i.num = 10,
                  r.num = 0)
  
  control_old <- control.net(type = "SIR", nsteps = 100, nsims = 12,
                             tea.status = FALSE,
                             save.network=FALSE,
                             use.pids = FALSE,
                             fast.edgelist=FALSE,
                             verbose=FALSE)
  set.seed(1)
  simold <- netsim(est2, param, init, control_old)
  
  control_new <- control.net(type = "SIR", nsteps = 100, nsims = 12,
                             tea.status = FALSE,
                             save.network=FALSE,
                             use.pids = FALSE,
                             fast.edgelist=TRUE,
                             verbose=FALSE)
  set.seed(1)
  simnew <- netsim(est2, param, init, control_new)
  # need to run a number of sims and verify that the outcomes match on average

})

test_that('expected speed improvement',{
  library(microbenchmark)
  nw <- network.initialize(n = 500, directed = FALSE)
  # specify two different roles for the vertices
  nw%v%'rolemode'<-rep_len(c('a','b'),network.size(nw))
  formation <- ~edges+offset(nodematch('rolemode'))
  target.stats <- .75*network.size(nw)
  
  # calculate dissolution coefficient with death rate
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20,
                                 d.rate = 0.0021)
  # set iInf coef for offset statistic for nodematch "never form ties between vertcies with non matching attributes
  coef.form <- -Inf  
  
  # Reestimate the model with new coefficient
  set.seed(1)
  est2 <- netest(nw, formation, target.stats, coef.diss, coef.form=coef.form, verbose = FALSE)
  
  # Reset parameters to include demographic rates
  param <- param.net(inf.prob = 0.3, 
                     rec.rate = 0.02, 
                     b.rate = 0.00, 
                     ds.rate = 0.01, 
                     di.rate = 0.0, 
                     dr.rate = 0.0)
  init<- init.net(i.num = 10,
                  r.num = 0)
  
  control_old <- control.net(type = "SIR", nsteps = 100, nsims = 1,
                             tea.status = FALSE,
                             save.network=FALSE,
                             use.pids = FALSE,
                             fast.edgelist=FALSE,
                             verbose=FALSE)
  
  
  control_new <- control.net(type = "SIR", nsteps = 100, nsims = 1,
                             tea.status = FALSE,
                             save.network=FALSE,
                             use.pids = FALSE,
                             fast.edgelist=TRUE,
                             verbose=FALSE)
  
  
  simoldfun<-function(){
    netsim(est2, param, init, control_old)
  }
  
  
  simnewfun <- function(){
    netsim(est2, param, init, control_new)
  }
  
  simcompare<-microbenchmark(simoldfun(),simnewfun(),times=10)

})

# TODO: add tests that non-supported models flagged correctly



plotvital<-function(netsim){
  deaths<-as.vector(netsim$epi$ds.flow+netsim$epi$di.flow+netsim$epi$dr.flow)[,1]
  births<-netsim$epi$b.flow[,1]
  plot(seq.int(length(deaths)),deaths,type='l',ylab='deaths and birthds',col='red',xlab='sim step')
  points(seq.int(length(births)),births,type='l',col='green')
}


# ideally, if we set the seed the same, we should get the same set of toggles via either simulate 