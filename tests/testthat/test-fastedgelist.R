

# tests for evaluating the fast.edgelist mode and ensuring that its output matches
context('fast.edgelist mode')

# number of sims to run for time comparisons set this to a high value (200)
# if actually comparying output, low value if just checking that the terms run
NUMSIMS <-1  

# number of timesteps to run each model.  100 for realistic output, 5 to just make sure they run
NUMSTEPS <-5

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
  
  control_old <- control.net(type = "SIR", nsteps = NUMSTEPS, nsims = NUMSIMS,
                             tea.status = FALSE,
                             save.network=TRUE,
                             use.pids = FALSE,
                             fast.edgelist=FALSE,
                             verbose=FALSE)
  set.seed(1)
  simold <- netsim(est2, param, init, control_old)
  
  control_new <- control.net(type = "SIR", nsteps = NUMSTEPS, nsims = NUMSIMS,
                             tea.status = FALSE,
                             save.network=TRUE,
                             use.pids = FALSE,
                             fast.edgelist=TRUE,
                             verbose=FALSE)
  set.seed(1)
  simnew <- netsim(est2, param, init, control_new)
  

  
  # edgelist sorting code disabled, so don't expect transmats to be identical
  # expect_equal(simold$stats$transmat$sim1,simnew$stats$transmat$sim1)
})

test_that('edges+nodematch model works',{
  sims = NUMSIMS
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
  
  control_old <- control.net(type = "SIR", nsteps = NUMSTEPS, nsims = sims,
                             tea.status = FALSE,
                             save.network=FALSE,
                             save.transmat = FALSE,
                             use.pids = FALSE,
                             fast.edgelist=FALSE,
                             verbose=TRUE)
  set.seed(1)
  simold <- netsim(est2, param, init, control_old)
  
  control_new <- control.net(type = "SIR", nsteps = NUMSTEPS, nsims = sims,
                             tea.status = FALSE,
                             save.network=FALSE,
                             save.transmat = FALSE,
                             use.pids = FALSE,
                             fast.edgelist=TRUE,
                             verbose=TRUE)
  set.seed(1)
  simnew <- netsim(est2, param, init, control_new)
  # need to run a number of sims and verify that the outcomes match on average
  
  pdf('edgelistPrevTimeseriesComparison.pdf')
  plot(simold, y = 'i.num', qnts = 0.5,main=paste('i.num comparison for',sims, 'net and edgelist models'))
  plot(simnew, y = 'i.num', qnts = 0.5,add=TRUE)
  dev.off()
  
  pdf('edgelistFinalPrevBoxplot.pdf')
  boxplot(x=list(net=as.numeric(simold$epi$i.num[NUMSTEPS,]),
                 edgelist=as.numeric(simnew$epi$i.num[NUMSTEPS,])),
          main='comparison final i.num, 500 sims, 100 steps')
  dev.off()
  IQR(as.numeric(simold$epi$i.num[NUMSTEPS,]))
  IQR(as.numeric(simnew$epi$i.num[NUMSTEPS,]))

})

test_that('edges+nodematch(diff=TRUE) model works',{
  sims = NUMSIMS
  nw <- network.initialize(n = 100, directed = FALSE)
  # specify two different roles for the vertices
  nw%v%'rolemode'<-rep_len(c('a','b','c'),network.size(nw))
  foo <- TRUE  # this is to test that evaluating formula args in calling environment works
               
  valRange <- 2:3
  # formation <- ~edges+offset(nodematch('rolemode',diff=foo,keep=valRange))  # this seems to run fine in manual testing, but not inside testthat, different environment evaluation?
  formation <- ~edges+offset(nodematch('rolemode',diff=TRUE,keep=valRange))
  target.stats <- 50
  
  # calculate dissolution coefficient with death rate
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20,
                                 d.rate = 0.0021)
  # set iInf coef for offset statistic for nodematch "never form ties between vertcies with non matching attributes
  coef.form <- c(-Inf,-Inf)  
  
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
  
  control_old <- control.net(type = "SIR", nsteps = NUMSTEPS, nsims = sims,
                             tea.status = FALSE,
                             save.network=FALSE,
                             save.transmat = FALSE,
                             use.pids = FALSE,
                             fast.edgelist=FALSE,
                             verbose=TRUE)
  set.seed(1)
  simold <- netsim(est2, param, init, control_old)
  
  control_new <- control.net(type = "SIR", nsteps = NUMSTEPS, nsims = sims,
                             tea.status = FALSE,
                             save.network=FALSE,
                             save.transmat = FALSE,
                             use.pids = FALSE,
                             fast.edgelist=TRUE,
                             verbose=TRUE)
  set.seed(1)
  simnew <- netsim(est2, param, init, control_new)
  # need to run a number of sims and verify that the outcomes match on average
  
  pdf('edgelistPrevNodematchDiffTimeseriesComparison.pdf')
  plot(simold, y = 'i.num', qnts = 0.5,main=paste('i.num comparison for',sims, 'net and edgelist models'))
  plot(simnew, y = 'i.num', qnts = 0.5,add=TRUE)
  dev.off()
  
  pdf('edgelistFinalPrevNodematchDiffBoxplot.pdf')
  boxplot(x=list(net=as.numeric(simold$epi$i.num[NUMSTEPS,]),
                 edgelist=as.numeric(simnew$epi$i.num[NUMSTEPS,])),
          main=paste('comparison final i.num,',sims,'sims, ',NUMSTEPS,' steps'))
  dev.off()
  IQR(as.numeric(simold$epi$i.num[NUMSTEPS,]))
  IQR(as.numeric(simnew$epi$i.num[NUMSTEPS,]))
  
})
library(microbenchmark)

# test_that('expected speed improvement',{
#   nw <- network.initialize(n = 500, directed = FALSE)
#   # specify two different roles for the vertices
#   nw%v%'rolemode'<-rep_len(c('a','b'),network.size(nw))
#   formation <- ~edges+offset(nodematch('rolemode'))
#   target.stats <- .75*network.size(nw)
#   
#   # calculate dissolution coefficient with death rate
#   coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20,
#                                  d.rate = 0.0021)
#   # set iInf coef for offset statistic for nodematch "never form ties between vertcies with non matching attributes
#   coef.form <- -Inf  
#   
#   # Reestimate the model with new coefficient
#   set.seed(1)
#   est2 <- netest(nw, formation, target.stats, coef.diss, coef.form=coef.form, verbose = FALSE)
#   
#   # Reset parameters to include demographic rates
#   param <- param.net(inf.prob = 0.3, 
#                      rec.rate = 0.02, 
#                      b.rate = 0.00, 
#                      ds.rate = 0.01, 
#                      di.rate = 0.0, 
#                      dr.rate = 0.0)
#   init<- init.net(i.num = 10,
#                   r.num = 0)
#   
#   control_old <- control.net(type = "SIR", nsteps = 100, nsims = 1,
#                              tea.status = FALSE,
#                              save.network=FALSE,
#                              save.transmat = FALSE,
#                              use.pids = FALSE,
#                              fast.edgelist=FALSE,
#                              verbose=FALSE)
#   
#   
#   control_new <- control.net(type = "SIR", nsteps = 100, nsims = 1,
#                              tea.status = FALSE,
#                              save.network=FALSE,
#                              save.transmat = FALSE,
#                              use.pids = FALSE,
#                              fast.edgelist=TRUE,
#                              verbose=FALSE)
#   
#   
#   simoldfun<-function(){
#     netsim(est2, param, init, control_old)
#   }
#   
#   
#   simnewfun <- function(){
#     netsim(est2, param, init, control_new)
#   }
#   
#   simcompare<-microbenchmark(simoldfun(),simnewfun(),times=10)
# 
# })

# TODO: add tests that non-supported models flagged correctly
test_that('non-supported models give error',{
  sims = NUMSIMS
  nw <- network.initialize(n = 100, directed = FALSE)
  formation <- ~edges+offset(density)
  target.stats <- c(30)
  
  # calculate dissolution coefficient with death rate
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20,
                                 d.rate = 0.0021)
  coef.form <- c(30)  
  # Reestimate the model with new coefficient
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
  
  control_old <- control.net(type = "SIR", nsteps = 2, nsims = sims,
                             tea.status = FALSE,
                             save.network=FALSE,
                             save.transmat = FALSE,
                             use.pids = FALSE,
                             fast.edgelist=TRUE,
                             verbose=TRUE)
  expect_error(netsim(est2, param, init, control_old),regexp = 'fast_edgelist mode does not know how to update the term')
})

compareNetVsElmodel <- function(est,modelName,param,init,numSteps=100,numSims=1){
  
  control_noEl <- control.net(type = "SIR", nsteps = numSteps, nsims = numSims,
                              tea.status = FALSE,
                              save.network=FALSE,
                              save.transmat = FALSE,
                              use.pids = FALSE,
                              fast.edgelist=FALSE,
                              verbose=TRUE)
  message('starting network-based sims')
  simNoEl<- netsim(est, param, init, control_noEl)
  control_el <- control.net(type = "SIR", nsteps = numSteps, nsims = numSims,
                            tea.status = FALSE,
                            save.network=FALSE,
                            save.transmat = FALSE,
                            use.pids = FALSE,
                            fast.edgelist=TRUE,
                            verbose=TRUE)
  message('starting edgelist-based sims')
  simEl <- netsim(est, param, init, control_el)
  message('comparing results')
  # generate plots in external file
  
  pdf(file = paste('comparisonsFor',modelName,'.pdf',sep=''))
  
  # model summary info
  plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "",main='model params:')
  text(3,9,labels= est$formation)
  text(3,8,labels= paste(est$target.stats))
  text(rep(8,9),y = 9:1,labels = paste(names(simEl$param),simEl$param,sep='='))
  
  plot(simNoEl, y = 'i.num', qnts = 0.5,
       mean.col='blue',
       main=paste(modelName,'i.num for',numSims, 'net and el models'))
  plot(simEl, y = 'i.num', qnts = 0.5,add=TRUE,
       mean.col='green' )
  legend(0,1,c('net','el'),fill = c('blue','green'))
  
  boxplot(x=list(net=as.numeric(simNoEl$epi$i.num[numSteps,]),
                 edgelist=as.numeric(simEl$epi$i.num[numSteps,])),
          main=paste(modelName,'final i.num distribution,',numSims,'sims, ',numSteps,' steps'))
  
  dev.off()
  
  invisible(list(simNoEl,simEl))
}

test_that('edges+concurrent model works',{
  sims = NUMSIMS
  nw <- network.initialize(n = 100, directed = FALSE)
  formation <- ~edges+concurrent
  target.stats <- c(50,35)
  
  # calculate dissolution coefficient with death rate
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20,
                                 d.rate = 0.0021)
  
  # Reestimate the model with new coefficient
  set.seed(1)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
  
  models <- compareNetVsElmodel(est, 'Concurrent',param,init,numSteps = NUMSTEPS,numSims = NUMSIMS)
  
})

test_that('edges+nodefactor model works',{
  sims = NUMSIMS
  nw <- network.initialize(n = 100, directed = FALSE)
  # specify two different roles for the vertices
  nw%v%'rolemode'<-rep_len(c('a','b','c'),network.size(nw))
  formation <- ~edges+nodefactor('rolemode')
  target.stats <- c(50,30,30)
  
  # calculate dissolution coefficient with death rate
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20,
                                 d.rate = 0.0021)
  
  # Reestimate the model with new coefficient
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
  
  models <- compareNetVsElmodel(est, 'Nodefactor',param,init,numSteps = NUMSTEPS,numSims = NUMSIMS)
  
})





test_that('edges+absdiff model works',{
  nw <- network.initialize(n = 100, directed = FALSE)
  # specify random numeric value for vertices
  nw%v%'numParam'<-runif(network.size(nw))
  formation <- ~edges+absdiff('numParam')
  target.stats <- c(50,10)
  param <- param.net(inf.prob = 0.3, 
                     rec.rate = 0.02, 
                     b.rate = 0.00, # <- need to leave births at zero because have not defined module to give new values
                     ds.rate = 0.01, 
                     di.rate = 0.01, 
                     dr.rate = 0.01)
  init<- init.net(i.num = 10,
                  r.num = 0)
  
  # calculate dissolution coefficient with death rate
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 2,
                                 d.rate = 0.0021)
  
  # estimate the model with new coefficient
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
  
  models <- compareNetVsElmodel(est, 'Absdiff',param,init,numSteps = NUMSTEPS,numSims = NUMSIMS)
})

# test_that('edges+nodecov model works',{
#   nw <- network.initialize(n = 100, directed = FALSE)
#   # specify random numeric value for vertices
#   nw%v%'numParam'<-runif(network.size(nw))
#   formation <- ~edges+nodecov('numParam',transform=sum,transformname='Sum')
#   target.stats <- c(50,10)
#   param <- param.net(inf.prob = 0.3, 
#                      rec.rate = 0.02, 
#                      b.rate = 0.00, # <- need to leave births at zero because have not defined module to give new values
#                      ds.rate = 0.01, 
#                      di.rate = 0.01, 
#                      dr.rate = 0.01)
#   init<- init.net(i.num = 10,
#                   r.num = 0)
#   
#   # calculate dissolution coefficient with death rate
#   coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 2,
#                                  d.rate = 0.0021)
#   
#   # estimate the model with new coefficient
#   est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#   
#   models <- compareNetVsElmodel(est, 'Nodecov',param,init,numSteps = NUMSTEPS,numSims = NUMSIMS)
# })

test_that('edges+nodecov model works',{
  data(florentine)

  formation <- ~edges+nodecov('wealth')
  target.stats <- c(20,2168)
  param <- param.net(inf.prob = 0.3, 
                     rec.rate = 0.02, 
                     b.rate = 0.00, # <- need to leave births at zero because have not defined module to give new values
                     ds.rate = 0.01, 
                     di.rate = 0.01, 
                     dr.rate = 0.01)
  init<- init.net(i.num = 10,
                  r.num = 0)
  
  # calculate dissolution coefficient with death rate
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 2,
                                 d.rate = 0.0021)
  
  # estimate the model with new coefficient
  est <- netest(flomarriage, formation, target.stats, coef.diss, verbose = FALSE)
  
  models <- compareNetVsElmodel(est, 'Nodecov',param,init,numSteps = NUMSTEPS,numSims = NUMSIMS)
})

test_that('edges+nodemix model works',{
  nw <- network.initialize(n = 100, directed = FALSE)
  # specify random categorical value for vertices
  nw%v%'rolemode'<-rep_len(c('a','b','c'),network.size(nw))
  formation <- ~edges+nodemix('rolemode',base=1)
  target.stats <- c(50,12,8,8,9,8)
  param <- param.net(inf.prob = 0.3, 
                     rec.rate = 0.02, 
                     b.rate = 0.00, # <- need to leave births at zero because have not defined module to give new values
                     ds.rate = 0.01, 
                     di.rate = 0.01, 
                     dr.rate = 0.01)
  init<- init.net(i.num = 10,
                  r.num = 0)
  
  # calculate dissolution coefficient with death rate
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 2,
                                 d.rate = 0.0021)
  
  # estimate the model with new coefficient
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
  
  models <- compareNetVsElmodel(est, 'Nodmix',param,init,numSteps = NUMSTEPS,numSims = NUMSIMS)
})

test_that('edges+degree model works',{
  nw <- network.initialize(n = 100, directed = FALSE)
  formation <- ~edges+degree(1:3)
  target.stats <- c(50,44,14,8)
  param <- param.net(inf.prob = 0.3, 
                     rec.rate = 0.02, 
                     b.rate = 0.00, # <- need to leave births at zero because have not defined module to give new values
                     ds.rate = 0.01, 
                     di.rate = 0.01, 
                     dr.rate = 0.01)
  init<- init.net(i.num = 10,
                  r.num = 0)
  
  # calculate dissolution coefficient with death rate
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 2,
                                 d.rate = 0.0021)
  
  # estimate the model with new coefficient
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
  
  models <- compareNetVsElmodel(est, 'Degree',param,init,numSteps = NUMSTEPS,numSims = NUMSIMS)
})


plotvital<-function(netsim){
  deaths<-as.vector(netsim$epi$ds.flow+netsim$epi$di.flow+netsim$epi$dr.flow)[,1]
  births<-netsim$epi$b.flow[,1]
  plot(seq.int(length(deaths)),deaths,type='l',ylab='deaths and birthds',col='red',xlab='sim step')
  points(seq.int(length(births)),births,type='l',col='green')
}

