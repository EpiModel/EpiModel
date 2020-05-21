## ----setup, include = FALSE----------------------------------------------
require(knitr)
require(EpiModel)
opts_chunk$set(comment = NA, message = FALSE, tidy = FALSE)

## ----dcmSi1--------------------------------------------------------------
param <- param.dcm(inf.prob = 0.2, act.rate = 0.25)
init <- init.dcm(s.num = 500, i.num = 1)
control <- control.dcm(type = "SI", nsteps = 500)

## ----dcmSi2--------------------------------------------------------------
mod <- dcm(param, init, control)

## ----dcmSiPrint----------------------------------------------------------
mod

## ----dcmSiPlot-----------------------------------------------------------
plot(mod)

## ----dcmSiSumm-----------------------------------------------------------
summary(mod, at = 150)

## ----dcmSir--------------------------------------------------------------
param <- param.dcm(inf.prob = 0.2, act.rate = 1, rec.rate = 1/20,
                   a.rate = 1/95, ds.rate = 1/100, di.rate = 1/80, dr.rate = 1/100)
init <- init.dcm(s.num = 1000, i.num = 1, r.num = 0)
control <- control.dcm(type = "SIR", nsteps = 500, dt = 0.5)
mod <- dcm(param, init, control)

## ----dcmSirPlot, fig.height=4--------------------------------------------
par(mar = c(3.2, 3, 2, 1), mgp = c(2, 1, 0), mfrow = c(1, 2))
plot(mod, popfrac = FALSE, alpha = 0.5,
     lwd = 4, main = "Compartment Sizes")
plot(mod, y = "si.flow", lwd = 4, col = "firebrick",
     main = "Disease Incidence", legend = "n")

## ----dcmSirCPlot---------------------------------------------------------
par(mfrow = c(1, 1))
comp_plot(mod, at = 50, digits = 1)

## ----dcmSis, results='hide'----------------------------------------------
param <- param.dcm(inf.prob = 0.2, act.rate = seq(0.25, 0.5, 0.05), rec.rate = 0.02)
init <- init.dcm(s.num = 500, i.num = 1)
control <- control.dcm(type = "SIS", nsteps = 350)
mod <- dcm(param, init, control)

## ----dcmSisPrint---------------------------------------------------------
mod

## ----dcmSisHead----------------------------------------------------------
head(as.data.frame(mod, run = 5))

## ----dcmSisPlot, fig.height=4--------------------------------------------
par(mfrow = c(1,2), mar = c(3.2,3,2.5,1))
plot(mod, alpha = 1, main = "Disease Prevalence")
plot(mod, y = "si.flow", col = "Greens", alpha = 0.8, main = "Disease Incidence")

## ----dcmSisPlotOpts, eval=FALSE------------------------------------------
plot(mod, col = "black")
plot(mod, col = 1:6)
plot(mod, col = c("black", "red", "blue", "green", "purple", "pink"))
plot(mod, col = rainbow(3), lty = rep(1:2, each = 3), legend = "full")

## ----dcmSisSensOpts------------------------------------------------------
act.rates <- c(0.2, 0.2, 0.4, 0.4, 0.6, 0.6)
inf.probs <- c(0.1, 0.2, 0.1, 0.2, 0.1, 0.2)
param <- param.dcm(inf.prob = inf.probs, act.rate = act.rates,
                   rec.rate = 0.02)
mod <- dcm(param, init, control)
plot(mod)

## ----dcmSi2g, results="hide"---------------------------------------------
param <- param.dcm(inf.prob = 0.4,  inf.prob.g2 = 0.1, act.rate = 0.25, balance = "g1",
                   a.rate = 1/100, a.rate.g2 = NA, ds.rate = 1/100, ds.rate.g2 = 1/100,
                   di.rate = 1/50, di.rate.g2 = 1/50)
init <- init.dcm(s.num = 500, i.num = 1, s.num.g2 = 500, i.num.g2 = 0)
control <- control.dcm(type = "SI", nsteps = 500)
mod <- dcm(param, init, control)

## ----dcmSi2gPrint--------------------------------------------------------
mod

## ----dcmSi2gPlot---------------------------------------------------------
plot(mod)

