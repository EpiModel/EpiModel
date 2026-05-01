# EpiModel

[![CRAN
Version](https://www.r-pkg.org/badges/version/EpiModel)](https://CRAN.R-project.org/package=EpiModel)
[![Downloads](http://cranlogs.r-pkg.org/badges/EpiModel?color=blue)](https://CRAN.R-project.org/package=EpiModel)
[![Total
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/EpiModel?color=blue)](https://CRAN.R-project.org/package=EpiModel)
[![Build
Status](https://github.com/EpiModel/EpiModel/workflows/R-CMD-check/badge.svg)](https://github.com/EpiModel/EpiModel/actions)

  
![](https://epimodel.github.io/sismid/images/movie.gif)

EpiModel is an R package for simulating mathematical models of
infectious disease dynamics. Epidemic model classes include
deterministic compartmental models, stochastic individual-contact
models, and stochastic network models. Network models use the robust
statistical methods of exponential-family random graph models (ERGMs)
from the [Statnet](https://statnet.org) suite of R packages. Standard
templates for epidemic modeling include SI, SIR, and SIS disease types.
EpiModel features an API for extending these templates to address novel
research aims.

EpiModel is part of a broader platform for infectious disease modeling.
Visit **[epimodel.org](https://www.epimodel.org/)** for the full
ecosystem, learning pathway, team, and project information.

### Installation

The current release version can be found on
[CRAN](https://CRAN.R-project.org/package=EpiModel) and installed with:

``` r

install.packages("EpiModel", dependencies = TRUE)
```

To install the development version, use the
[remotes](https://github.com/r-lib/remotes) package:

``` r

if (!require("remotes")) install.packages("remotes")
remotes::install_github("EpiModel/EpiModel")
```

### Documentation

**Package Documentation.** Full function reference, vignettes, and
getting started guides are on the [pkgdown
site](https://epimodel.github.io/EpiModel/).

**Methods Paper.** The main methods paper in the *Journal of Statistical
Software* describes the modeling framework and software design:
[doi:10.18637/jss.v084.i08](https://doi.org/10.18637/jss.v084.i08).

**Training Course.** [Network Modeling for Epidemics
(NME)](https://epimodel.github.io/sismid/) is our course at the Summer
Institute in Statistics and Modeling in Infectious Diseases (SISMID).
Materials are open-source and updated annually.

**EpiModel Gallery.** The [EpiModel
Gallery](https://epimodel.github.io/EpiModel-Gallery/) provides worked
examples and templates for extending EpiModel’s stochastic network
models with custom epidemic modules.

### Getting Help

Users are encouraged to use [GitHub
Issues](https://github.com/EpiModel/EpiModel/issues) to ask questions
(both technical coding questions and conceptual modeling questions),
report bugs, and request new features.

### Citation

If using EpiModel for teaching or research, please include a citation to
our main methods paper:

> Jenness SM, Goodreau SM and Morris M. EpiModel: An R Package for
> Mathematical Modeling of Infectious Disease over Networks. *Journal of
> Statistical Software.* 2018; 84(8): 1-47. doi: 10.18637/jss.v084.i08

Please also [send us an
email](mailto:samuel.m.jenness@emory.edu?Subject=We%20Used%20EpiModel%20in%20Our%20Study!)
if you have used EpiModel in your work so we can add it to our records.

### Authors

EpiModel is developed by [Samuel
Jenness](https://www.epimodel.org/team.html), [Steven
Goodreau](https://www.epimodel.org/team.html), [Martina
Morris](https://www.epimodel.org/team.html), and [Adrien Le
Guillou](https://www.epimodel.org/team.html). Additional contributors
are listed on the
[contributors](https://github.com/EpiModel/EpiModel/graphs/contributors)
page. See the full team at
[epimodel.org/team](https://www.epimodel.org/team.html).

### Funding and Literature

The development of EpiModel has been supported by grants from the
National Institutes of Health (NIH) and the Centers for Disease Control
and Prevention (CDC). Details are available on the [project
page](https://www.epimodel.org/about.html).

EpiModel and its extension packages have been used in 125+ published
studies. See the [full
list](https://github.com/EpiModel/EpiModel/wiki/EpiModel-in-the-Scientific-Literature).

### Copyright

These materials are distributed under the GPL-3 license, with the
copyright and attribution requirements listed in the
[LICENSE](https://github.com/EpiModel/EpiModel/blob/main/LICENSE.md).
