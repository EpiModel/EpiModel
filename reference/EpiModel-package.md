# Mathematical Modeling of Infectious Disease Dynamics

|           |            |
|-----------|------------|
| Package:  | EpiModel   |
| Type:     | Package    |
| Version:  | 2.5.5      |
| Date:     | 2026-02-11 |
| License:  | GPL-3      |
| LazyLoad: | yes        |

## Details

The EpiModel software package provides tools for building, solving, and
visualizing mathematical models of infectious disease dynamics. These
tools allow users to simulate epidemic models in multiple frameworks for
both pedagogical purposes ("base models") and novel research purposes
("extension models").

## Model Classes and Infectious Disease Types

EpiModel provides functionality for three classes of epidemic models:

- **Deterministic Compartmental Models:** these continuous-time models
  are solved using ordinary differential equations. EpiModel allows for
  easy specification of sensitivity analyses to compare multiple
  scenarios of the same model across different parameter values.

- **Stochastic Individual Contact Models:** a novel class of
  individual-based, microsimulation models that were developed to add
  random variation in all components of the transmission system, from
  infection to recovery to vital dynamics (arrivals and departures).

- **Stochastic Network Models:** with the underlying statistical
  framework of temporal exponential random graph models (ERGMs) recently
  developed in the **Statnet** suite of software in R, network models
  over epidemics simulate edge (e.g., partnership) formation and
  dissolution stochastically according to a specified statistical model,
  with disease spread across that network.

EpiModel supports three infectious disease types to be run across all of
the three classes.

- **Susceptible-Infectious (SI):** a two-state disease in which there is
  life-long infection without recovery. HIV/AIDS is one example,
  although for this case it is common to model infection stages as
  separate compartments.

- **Susceptible-Infectious-Recovered (SIR):** a three-stage disease in
  which one has life-long recovery with immunity after infection.
  Measles is one example, but modern models for the disease also require
  consideration of vaccination patterns in the population.

- **Susceptible-Infectious-Susceptible (SIS):** a two-stage disease in
  which one may transition back and forth from the susceptible to
  infected states throughout life. Examples include bacterial sexually
  transmitted diseases like gonorrhea.

These basic disease types may be extended in any arbitrarily complex way
to simulate specific diseases for research questions.

## Model Parameterization and Simulation

EpiModel uses three model setup functions for each model class to input
the necessary parameters, initial conditions, and control settings:

- [`param.dcm()`](http://epimodel.github.io/EpiModel/reference/param.dcm.md),
  [`param.icm()`](http://epimodel.github.io/EpiModel/reference/param.icm.md),
  and
  [`param.net()`](http://epimodel.github.io/EpiModel/reference/param.net.md)
  are used to input epidemic parameters for each of the three model
  classes. Parameters include the rate of contacts or acts between
  actors, the probability of transmission per contact, and recovery and
  demographic rates for models that include those transitions.

- [`init.dcm()`](http://epimodel.github.io/EpiModel/reference/init.dcm.md),
  [`init.icm()`](http://epimodel.github.io/EpiModel/reference/init.icm.md),
  and
  [`init.net()`](http://epimodel.github.io/EpiModel/reference/init.net.md)
  are used to input the initial conditions for each class. The main
  conditions are limited to the numbers or, if applicable, the specific
  agents in the population who are infected or recovered at the
  simulation outset.

- [`control.dcm()`](http://epimodel.github.io/EpiModel/reference/control.dcm.md),
  [`control.icm()`](http://epimodel.github.io/EpiModel/reference/control.icm.md),
  and
  [`control.net()`](http://epimodel.github.io/EpiModel/reference/control.net.md)
  are used to specify the remaining control settings for each
  simulation. The core controls for base model types include the disease
  type, number of time steps, and number of simulations. Controls are
  also used to input new model functions (for DCMs) and new model
  modules (for ICMs and network models) to allow the user to simulate
  fully original epidemic models in EpiModel. See the documentation for
  the specific control functions help pages.

With the models parameterized, the functions for simulating epidemic
models are:

- [`dcm()`](http://epimodel.github.io/EpiModel/reference/dcm.md) for
  deterministic compartmental models.

- [`icm()`](http://epimodel.github.io/EpiModel/reference/icm.md) for
  individual contact models.

- Network models are simulated in a three-step process:

  1.  [`netest()`](http://epimodel.github.io/EpiModel/reference/netest.md)
      estimates the statistical model for the network structure itself
      (i.e., how partnerships form and dissolve over time given the
      parameterization of those processes). This function is a wrapper
      around the `ergm` and `tergm` functions in the `ergm` and `tergm`
      packages. The current statistical framework for model simulation
      is called "egocentric inference": target statistics summarizing
      these formation and dissolution processes collected from an
      egocentric sample of the population.

  2.  [`netdx()`](http://epimodel.github.io/EpiModel/reference/netdx.md)
      runs diagnostics on the dynamic model fit by simulating the base
      network over time to ensure the model fits the targets for
      formation and dissolution.

  3.  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md)
      simulates the stochastic network epidemic models, with a given
      network model fit in
      [`netest()`](http://epimodel.github.io/EpiModel/reference/netest.md).
      Here the function requires this model fit object along with the
      parameters, initial conditions, and control settings as defined
      above.

## References

The EpiModel website is at <https://www.epimodel.org/>, and the source
code is at <https://github.com/EpiModel/EpiModel>. Bug reports and
feature requests are welcome.

Our primary methods paper on EpiModel is published in the **Journal of
Statistical Software**. If you use EpiModel for any research or teaching
purposes, please cite this reference:

Jenness SM, Goodreau SM, and Morris M. EpiModel: An R Package for
Mathematical Modeling of Infectious Disease over Networks. Journal of
Statistical Software. 2018; 84(8): 1-47.
[doi:10.18637/jss.v084.i08](https://doi.org/10.18637/jss.v084.i08) .

We have also developed two extension packages for modeling specific
disease dynamics. For HIV and bacterial sexually transmitted infections,
we have developed `EpiModelHIV`, which is available on Github at
<https://github.com/EpiModel/EpiModelHIV>. For COVID-19, we have
developed `EpiModelCOVID`, which is available at
<https://github.com/EpiModel/EpiModelCOVID>.

## See also

Useful links:

- <https://www.epimodel.org/>

- <https://epimodel.github.io/EpiModel/>

- Report bugs at <https://github.com/EpiModel/EpiModel/issues/>

## Author

**Maintainer**: Samuel Jenness <samuel.m.jenness@emory.edu>

Authors:

- Steven M. Goodreau <goodreau@uw.edu>

- Martina Morris <morrism@uw.edu>

- Adrien Le Guillou <contact@aleguillou.org>

- Chad Klumb <cklumb@gmail.com>

Other contributors:

- Skye Bender-deMoll <skyebend@uw.edu> \[contributor\]
