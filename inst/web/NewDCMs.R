## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(comment = NA)

## ----load, message = FALSE-----------------------------------------------
library(EpiModel)

## ----introParams---------------------------------------------------------
param <- param.dcm(inf.prob = 0.5, act.rate = 0.25)
init <- init.dcm(s.num = 500, i.num = 1)
control <- control.dcm(type = "SI", nsteps = 100)
mod <- dcm(param, init, control)
plot(mod)

## ----modFunc1------------------------------------------------------------
control <- control.dcm(type = "SI", nsteps = 100, print.mod = TRUE)
mod <- dcm(param, init, control)

## ----dcmModsHelp, eval = FALSE-------------------------------------------
?dcm.mods

## ----printSImodEx--------------------------------------------------------
print(mod_SI_1g_cl)

## ------------------------------------------------------------------------
SEIR <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {

    # Population size
    num <- s.num + e.num + i.num + r.num

    # Effective contact rate and FOI from a rearrangement of Beta * c * D
    ce <- R0 / i.dur
    lambda <- ce * i.num / num

    dS <- -lambda * s.num
    dE <- lambda * s.num - (1 / e.dur) * e.num
    dI <- (1 / e.dur) * e.num - (1 - cfr) * (1 / i.dur) * i.num - cfr *
      (1 / i.dur) * i.num
    dR <- (1 - cfr) * (1 / i.dur) * i.num

    # Compartments and flows are part of the derivative vector
    # Other calculations to be output are outside the vector, but within the
    # containing list
    list(c(dS, dE, dI, dR,
           se.flow = lambda * s.num,
           ei.flow = (1 / e.dur) * e.num,
           ir.flow = (1 - cfr) * (1 / i.dur) * i.num,
           d.flow = cfr * (1 / i.dur) * i.num),
         num = num,
         i.prev = i.num / num,
         ei.prev = (e.num + i.num) / num)
  })
}

## ---- results = "hide"---------------------------------------------------
param <- param.dcm(R0 = 1.9, e.dur = 10, i.dur = 14, cfr = c(0.5, 0.7, 0.9))
init <- init.dcm(s.num = 1e6, e.num = 10, i.num = 0, r.num = 0,
                 se.flow = 0, ei.flow = 0, ir.flow = 0, d.flow = 0)
control <- control.dcm(nsteps = 500, dt = 1, new.mod = SEIR)
mod <- dcm(param, init, control)

## ------------------------------------------------------------------------
mod

## ------------------------------------------------------------------------
par(mfrow = c(1, 2))
plot(mod, y = "i.num", run = 2, main = "Prevalence")
plot(mod, y = "se.flow", run = 2, main = "Incidence")

## ------------------------------------------------------------------------
par(mfrow = c(1, 2))
plot(mod, y = "i.num", main = "Number Infected")
plot(mod, y = "i.prev", main = "Percent Infected", ylim = c(0, 0.5),
     legend = "full")

## ----qmodFunc------------------------------------------------------------
Qmod <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {

    ## Dynamic Calculations ##

    # Population size and prevalence
    h.num <- sh.num + ih.num
    l.num <- sl.num + il.num
    num <- h.num + l.num
    prev <- (ih.num + il.num) / num

    # Contact rates for high specified as a function of
    #   mean and low rates
    c.high <- (c.mean * num - c.low * l.num) / h.num

    # Mixing matrix calculations based on variable Q statistic
    g.hh <- ((c.high * h.num) + (Q * c.low * l.num)) /
            ((c.high * h.num) + (c.low * l.num))
    g.lh <- 1 - g.hh
    g.hl <- (1 - g.hh) * ((c.high * h.num) / (c.low * l.num))
    g.ll <- 1 - g.hl

    # Probability that contact is infected based on mixing probabilities
    p.high <- (g.hh * ih.num / h.num) + (g.lh * il.num / l.num)
    p.low <- (g.ll * il.num / l.num) + (g.hl * ih.num / h.num)

    # Force of infection for high and low groups
    lambda.high <- rho * c.high * p.high
    lambda.low <- rho * c.low * p.low

    ## Derivatives ##
    dS.high <- -lambda.high * sh.num + nu * ih.num
    dI.high <- lambda.high * sh.num - nu * ih.num

    dS.low <- -lambda.low * sl.num + nu * il.num
    dI.low <- lambda.low * sl.num - nu * nil.num

    ## Output ##
    list(c(dS.high, dI.high, dS.low, dI.low),
         num = num, prev = prev)
  })
}

## ----qmodParams----------------------------------------------------------
param <- param.dcm(c.mean = 2, c.low = 1.4, rho = 0.75, nu = 6,
                   Q = c(-0.45, -0.33, 0, 0.5, 1))

## ----qmodInits-----------------------------------------------------------
init <- init.dcm(sh.num = 2e7 * 0.02, ih.num = 1,
                 sl.num = 2e7 * 0.98, il.num = 1)

## ----qmodControls--------------------------------------------------------
control <- control.dcm(nsteps = 25, dt = 0.02, new.mod = Qmod)

## ----qmodRunMod, results = "hide"----------------------------------------
mod <- dcm(param, init, control)

## ----qmodPrint-----------------------------------------------------------
mod

## ----qmodADF-------------------------------------------------------------
head(as.data.frame(mod))
head(as.data.frame(mod, run = 5))

## ----qmodPlot, fig.width = 9.5-------------------------------------------
par(mfrow = c(1, 2))
plot(mod, y = "ih.num", legend = "full", main = "Infected High")
plot(mod, y = "il.num", legend = "full", main = "Infected Low")

## ----qmodPlot2-----------------------------------------------------------
par(mfrow = c(1, 1))
plot(mod, y = "prev", ylim = c(0, 0.05), legend = "full",
     main = "Overall Prevalence")
