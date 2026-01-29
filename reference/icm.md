# Stochastic Individual Contact Models

Simulates stochastic individual contact epidemic models for infectious
disease.

## Usage

``` r
icm(param, init, control)
```

## Arguments

- param:

  Model parameters, as an object of class
  [`param.icm()`](http://epimodel.github.io/EpiModel/reference/param.icm.md).

- init:

  Initial conditions, as an object of class
  [`init.icm()`](http://epimodel.github.io/EpiModel/reference/init.icm.md).

- control:

  Control settings, as an object of class
  [`control.icm()`](http://epimodel.github.io/EpiModel/reference/control.icm.md).

## Value

A list of class `icm` with the following elements:

- **param:** the epidemic parameters passed into the model through
  `param`, with additional parameters added as necessary.

- **control:** the control settings passed into the model through
  `control`, with additional controls added as necessary.

- **epi:** a list of data frames, one for each epidemiological output
  from the model. Outputs for base models always include the size of
  each compartment, as well as flows in, out of, and between
  compartments.

## Details

Individual contact models are intended to be the stochastic
microsimulation analogs to deterministic compartmental models. ICMs
simulate disease spread on individual agents in discrete time as a
function of processes with stochastic variation. The stochasticity is
inherent in all transition processes: infection, recovery, and
demographics.

The `icm` function performs modeling of both the base model types and
original models. Base model types include one-group and two-group models
with disease types for Susceptible-Infected (SI),
Susceptible-Infected-Recovered (SIR), and
Susceptible-Infected-Susceptible (SIS). Original models may be built by
writing new process modules that either take the place of existing
modules (for example, disease recovery), or supplement the set of
existing processes with a new one contained in an original module.

## See also

Extract the model results with
[`as.data.frame.icm()`](http://epimodel.github.io/EpiModel/reference/as.data.frame.icm.md).
Summarize the time-specific model results with
[`summary.icm()`](http://epimodel.github.io/EpiModel/reference/summary.icm.md).
Plot the model results with
[`plot.icm()`](http://epimodel.github.io/EpiModel/reference/plot.icm.md).
Plot a compartment flow diagram with
[`comp_plot()`](http://epimodel.github.io/EpiModel/reference/comp_plot.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## Example 1: SI Model
param <- param.icm(inf.prob = 0.2, act.rate = 0.25)
init <- init.icm(s.num = 500, i.num = 1)
control <- control.icm(type = "SI", nsteps = 500, nsims = 10)
mod1 <- icm(param, init, control)
mod1
plot(mod1)

## Example 2: SIR Model
param <- param.icm(inf.prob = 0.2, act.rate = 0.25, rec.rate = 1/50)
init <- init.icm(s.num = 500, i.num = 1, r.num = 0)
control <- control.icm(type = "SIR", nsteps = 500, nsims = 10)
mod2 <- icm(param, init, control)
mod2
plot(mod2)

## Example 3: SIS Model
param <- param.icm(inf.prob = 0.2, act.rate = 0.25, rec.rate = 1/50)
init <- init.icm(s.num = 500, i.num = 1)
control <- control.icm(type = "SIS", nsteps = 500, nsims = 10)
mod3 <- icm(param, init, control)
mod3
plot(mod3)

## Example 4: SI Model with Vital Dynamics (Two-Group)
param <- param.icm(inf.prob = 0.4,  inf.prob.g2 = 0.1,
                   act.rate = 0.25, balance = "g1",
                   a.rate = 1/100, a.rate.g2 = NA,
                   ds.rate = 1/100, ds.rate.g2 = 1/100,
                   di.rate = 1/50, di.rate.g2 = 1/50)
init <- init.icm(s.num = 500, i.num = 1,
                 s.num.g2 = 500, i.num.g2 = 0)
control <- control.icm(type = "SI", nsteps = 500, nsims = 10)
mod4 <- icm(param, init, control)
mod4
plot(mod4)
} # }
```
