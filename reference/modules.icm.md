# Modules for Stochastic Individual Contact Models

Stochastic individual contact models of infectious disease simulate
epidemics in which contacts between individuals are instantaneous events
in discrete time. They are intended to be the stochastic microsimulation
analogs to deterministic compartmental models.

The [`icm()`](https://epimodel.github.io/EpiModel/reference/icm.md)
function handles both the simulation tasks. Within this function are a
series of modules that initialize the simulation and then simulate new
infections, recoveries, and vital dynamics at each time step. A module
also handles the basic bookkeeping calculations for disease prevalence.

This help page presents a brief overview of the built-in module
functions in the order in which they are used within
[`icm()`](https://epimodel.github.io/EpiModel/reference/icm.md). ICMs
support only the built-in SI, SIR, and SIS disease types. For custom or
extension epidemic models, use the network model class via
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
and
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
instead.

## Initialization Module

This function sets up agent attributes, like disease status, on the
network at the starting time step of disease simulation, \\t_1\\. For
multiple-simulation function calls, these are reset at the beginning of
each simulation.

- [`initialize.icm()`](https://epimodel.github.io/EpiModel/reference/initialize.icm.md):
  sets which agents are initially infected, through the initial
  conditions passed in
  [`init.icm()`](https://epimodel.github.io/EpiModel/reference/init.icm.md).

## Disease Status Modification Modules

The main disease simulation occurs at each time step given the current
state of the population at that step. Infection of agents is simulated
as a function of disease parameters and population composition. Recovery
of agents is likewise simulated with respect to infected nodes. These
functions also analyze the flows for summary measures such as disease
incidence.

- [`infection.icm()`](https://epimodel.github.io/EpiModel/reference/infection.icm.md):
  randomly draws an edgelist given the parameters, subsets the list for
  discordant pairs, and simulates transmission on those discordant pairs
  through a series of draws from a binomial distribution.

- [`recovery.icm()`](https://epimodel.github.io/EpiModel/reference/recovery.icm.md):
  simulates recovery from infection either to a lifelong immune state
  (for SIR models) or back to the susceptible state (for SIS models), as
  a function of the recovery rate specified in the `rec.rate` parameter.
  The recovery rate may vary for two-group models.

## Demographic Modules

Vital dynamics such as arrival and departure processes are simulated at
each time step to update entries into and exits from the population.
These are used in open-population ICMs.

- [`departures.icm()`](https://epimodel.github.io/EpiModel/reference/departures.icm.md):
  randomly simulates departures or exits for agents given the departure
  rate specified in the disease-state and group-specific departure
  parameters in
  [`param.icm()`](https://epimodel.github.io/EpiModel/reference/param.icm.md).
  This involves deactivating agents from the population, but their
  historical data is preserved in the simulation.

- [`arrivals.icm()`](https://epimodel.github.io/EpiModel/reference/arrivals.icm.md):
  randomly simulates new arrivals into the population given the current
  population size and the arrival rate parameters. This involves adding
  new agents into the population.

## Bookkeeping Module

Simulations require bookkeeping at each time step to calculate the
summary epidemiological statistics used in the model output analysis.

- [`prevalence.icm()`](https://epimodel.github.io/EpiModel/reference/prevalence.icm.md):
  calculates the number in each disease state (susceptible, infected,
  recovered) at each time step for those active agents in the
  population.
