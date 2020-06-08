
library(EpiModel)

nsteps <- 100
nsims = 1e3

my.control <- control.icm(type='SIR', nsteps=nsteps, nsims=nsims, verbose=TRUE)
my.param <- param.icm(rec.rate=0.1, act.rate = 0.2, inf.prob = 1)
my.init <- init.icm(s.num = 9, i.num = 1, r.num=0)
aaa <- icm(param = my.param, control = my.control, init = my.init)
aaa <- as.data.frame(aaa)
aaa <- aaa[aaa$time==100,]
hist(aaa$s.num)


