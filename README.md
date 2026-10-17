# EpiModel

<!-- badges: start -->
[![CRAN Version](https://www.r-pkg.org/badges/version/EpiModel)](https://CRAN.R-project.org/package=EpiModel) [![Downloads](https://cranlogs.r-pkg.org/badges/EpiModel?color=blue)](https://CRAN.R-project.org/package=EpiModel) [![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/EpiModel?color=blue)](https://CRAN.R-project.org/package=EpiModel) [![R-CMD-check](https://github.com/EpiModel/EpiModel/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/EpiModel/EpiModel/actions/workflows/check-standard.yaml)
<!-- badges: end -->

<picture>
  <source srcset="https://www.epimodel.org/img/network-logo-static.svg" media="(prefers-reduced-motion: reduce)">
  <img src="https://www.epimodel.org/img/network-logo.svg" align="right" width="230" alt="Animated simulation of an infection spreading across a dynamic contact network, simulated with EpiModel">
</picture>

**Mathematical modeling of infectious disease dynamics in R.**

EpiModel simulates epidemics with three classes of models, from deterministic compartmental models to stochastic network models in which partnerships form, persist, and dissolve according to rules estimated with temporal exponential-family random graph models (TERGMs) from the [Statnet](https://statnet.org) suite of R packages. Built-in SI, SIR, and SIS models run out of the box, and an extension API supports custom disease states, demography, and interventions for new research questions.

| Model class | Simulate with | Contact structure |
|---|---|---|
| Deterministic compartmental models | `dcm()` | Population-level ODEs with random mixing |
| Stochastic individual contact models | `icm()` | Individuals with random mixing |
| Stochastic network models | `netest()`, `netdx()`, `netsim()` | Dynamic partnerships estimated from network data |

All three share one interface: set inputs with `param.*()`, `init.*()`, and `control.*()`, then simulate.

EpiModel is the core package of a broader platform for infectious disease modeling. Visit **[epimodel.org](https://www.epimodel.org/)** for the extension packages, the learning pathway, and the team.

## Installation

The current release version can be found on [CRAN](https://CRAN.R-project.org/package=EpiModel) and installed with:

``` r
install.packages("EpiModel", dependencies = TRUE)
```

To install the development version, use the [remotes](https://github.com/r-lib/remotes) package:

``` r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("EpiModel/EpiModel")
```

## Example

A network model in three steps: estimate a dynamic network from target statistics, check that simulations from it reproduce those targets, then simulate an epidemic over it. This runs in under a minute on a laptop.

``` r
library(EpiModel)

# 1. Estimate: 500 people, 175 partnerships, 110 with 2+ partners, none with 4+
nw <- network_initialize(n = 500)
est <- netest(nw,
  formation = ~edges + concurrent + degrange(from = 4),
  target.stats = c(175, 110, 0),
  coef.diss = dissolution_coefs(~offset(edges), duration = 50)  # mean duration
)

# 2. Diagnose: simulated networks should match the targets
dx <- netdx(est, nsims = 5, nsteps = 500)
plot(dx)

# 3. Simulate an SIS epidemic over the dynamic network
param <- param.net(inf.prob = 0.4, act.rate = 2, rec.rate = 0.05)
init <- init.net(i.num = 10)
control <- control.net(type = "SIS", nsims = 5, nsteps = 500)
sim <- netsim(est, param, init, control)
plot(sim)
```

The [Intro to EpiModel](https://epimodel.github.io/EpiModel/articles/Intro.html) vignette walks through each step and the analysis of the output.

## Learning EpiModel

<table class="readme-resources">
<tr>
<td align="center" valign="top" width="33%">
<a href="https://epimodel.github.io/EpiModel/"><img src="https://www.epimodel.org/img/thumb-package.png" height="110" alt="Stylized EpiModel R code on a dark background"></a><br>
<b><a href="https://epimodel.github.io/EpiModel/">Package Documentation</a></b><br>
Function reference, vignettes, and getting started guides
</td>
<td align="center" valign="top" width="33%">
<a href="https://epimodel.github.io/sismid/"><img src="https://www.epimodel.org/img/thumb-training.png" height="110" alt="Prevalence over time in a stochastic network SIS model"></a><br>
<b><a href="https://epimodel.github.io/sismid/">Training Materials</a></b><br>
Network Modeling for Epidemics (NME), our open course at SISMID
</td>
<td align="center" valign="top" width="33%">
<a href="https://epimodel.github.io/EpiModel-Gallery/"><img src="https://www.epimodel.org/img/thumb-gallery.png" height="110" alt="SEIR and SIS compartmental model flow diagrams"></a><br>
<b><a href="https://epimodel.github.io/EpiModel-Gallery/">EpiModel Gallery</a></b><br>
Worked examples and templates for extending network models
</td>
</tr>
</table>

The main methods paper in the *Journal of Statistical Software* describes the modeling framework and software design: [doi:10.18637/jss.v084.i08](https://doi.org/10.18637/jss.v084.i08). The [learning pathway](https://www.epimodel.org/) on epimodel.org orders these resources from beginning to advanced.

## Getting Help

Users are encouraged to use [GitHub Issues](https://github.com/EpiModel/EpiModel/issues) to ask questions (both technical coding questions and conceptual modeling questions), report bugs, and request new features.

## Citation

If using EpiModel for teaching or research, please include a citation to our main methods paper:

> Jenness SM, Goodreau SM and Morris M. EpiModel: An R Package for Mathematical Modeling of Infectious Disease over Networks. *Journal of Statistical Software.* 2018; 84(8): 1-47. doi: 10.18637/jss.v084.i08

Please also [send us an email](mailto:samuel.m.jenness@emory.edu?Subject=We%20Used%20EpiModel%20in%20Our%20Study!) if you have used EpiModel in your work so we can add it to our records. EpiModel and its extension packages have been used in 130+ published studies; see the [full list](https://github.com/EpiModel/EpiModel/wiki/EpiModel-in-the-Scientific-Literature).

## Authors and Funding

EpiModel is developed by [Samuel Jenness](https://www.epimodel.org/team.html), [Steven Goodreau](https://www.epimodel.org/team.html), [Martina Morris](https://www.epimodel.org/team.html), and [Adrien Le Guillou](https://www.epimodel.org/team.html). Additional contributors are listed on the [contributors](https://github.com/EpiModel/EpiModel/graphs/contributors) page.

Development has been supported by grants from the National Institutes of Health (NIH) and the Centers for Disease Control and Prevention (CDC). Details are on the [project page](https://www.epimodel.org/about.html).

## License

These materials are distributed under the GPL-3 license, with the copyright and attribution requirements listed in the [LICENSE](https://github.com/EpiModel/EpiModel/blob/main/LICENSE.md).
