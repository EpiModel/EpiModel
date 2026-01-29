# Modules for Stochastic Individual Contact Models

Stochastic individual contact models of infectious disease simulate
epidemics in which contacts between individuals are instantaneous events
in discrete time. They are intended to be the stochastic microsimulation
analogs to deterministic compartmental models.

The [`icm()`](http://epimodel.github.io/EpiModel/reference/icm.md)
function handles both the simulation tasks. Within this function are a
series of modules that initialize the simulation and then simulate new
infections, recoveries, and vital dynamics at each time step. A module
also handles the basic bookkeeping calculations for disease prevalence.

Writing original ICMs will require modifying the existing modules or
adding new modules to the workflow in
[`icm()`](http://epimodel.github.io/EpiModel/reference/icm.md). The
existing modules may be used as a template for replacement or new
modules.

This help page presents a brief overview of the module functions in the
order in which they are used within
[`icm()`](http://epimodel.github.io/EpiModel/reference/icm.md), in order
to help guide users in writing their own module functions. These module
functions are not shown on the help index since they are not called
directly by the end-user. To understand these functions in more detail,
review the separate help pages listed below.

## Initialization Module

This function sets up agent attributes, like disease status, on the
network at the starting time step of disease simulation, \\t_1\\. For
multiple-simulation function calls, these are reset at the beginning of
each simulation.

- [`initialize.icm()`](http://epimodel.github.io/EpiModel/reference/initialize.icm.md):
  sets which agents are initially infected, through the initial
  conditions passed in
  [`init.icm()`](http://epimodel.github.io/EpiModel/reference/init.icm.md).

## Disease Status Modification Modules

The main disease simulation occurs at each time step given the current
state of the population at that step. Infection of agents is simulated
as a function of disease parameters and population composition. Recovery
of agents is likewise simulated with respect to infected nodes. These
functions also analyze the flows for summary measures such as disease
incidence.

- [`infection.icm()`](http://epimodel.github.io/EpiModel/reference/infection.icm.md):
  randomly draws an edgelist given the parameters, subsets the list for
  discordant pairs, and simulates transmission on those discordant pairs
  through a series of draws from a binomial distribution.

- [`recovery.icm()`](http://epimodel.github.io/EpiModel/reference/recovery.icm.md):
  simulates recovery from infection either to a lifelong immune state
  (for SIR models) or back to the susceptible state (for SIS models), as
  a function of the recovery rate specified in the `rec.rate` parameter.
  The recovery rate may vary for two-group models.

## Demographic Modules

Vital dynamics such as arrival and departure processes are simulated at
each time step to update entries into and exits from the population.
These are used in open-population ICMs.

- [`departures.icm()`](http://epimodel.github.io/EpiModel/reference/departures.icm.md):
  randomly simulates departures or exits for agents given the departure
  rate specified in the disease-state and group-specific departure
  parameters in
  [`param.icm()`](http://epimodel.github.io/EpiModel/reference/param.icm.md).
  This involves deactivating agents from the population, but their
  historical data is preserved in the simulation.

- [`arrivals.icm()`](http://epimodel.github.io/EpiModel/reference/arrivals.icm.md):
  randomly simulates new arrivals into the population given the current
  population size and the arrival rate parameters. This involves adding
  new agents into the population.

## Bookkeeping Module

Simulations require bookkeeping at each time step to calculate the
summary epidemiological statistics used in the model output analysis.

- [`prevalence.icm()`](http://epimodel.github.io/EpiModel/reference/prevalence.icm.md):
  calculates the number in each disease state (susceptible, infected,
  recovered) at each time step for those active agents in the
  population.
