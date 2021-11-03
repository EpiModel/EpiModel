# EpiModel

[![Version](http://img.shields.io/badge/Version-2.1.0-orange.svg?style=flat)](https://github.com/EpiModel/EpiModel/releases/tag/v2.1.0) [![](http://cranlogs.r-pkg.org/badges/EpiModel?color=blue)](http://cran.rstudio.com/web/packages/EpiModel/index.html) [![](http://cranlogs.r-pkg.org/badges/grand-total/EpiModel?color=blue)](http://cran.rstudio.com/web/packages/EpiModel/index.html) [![Build Status](https://github.com/EpiModel/EpiModel/workflows/R-CMD-check/badge.svg)](https://github.com/EpiModel/EpiModel/actions) [![Methods](https://img.shields.io/badge/docs-Methods-943ad8.svg)](http://doi.org/10.18637/jss.v084.i08)

<br> <img src="http://www.epimodel.org/movie.gif" align="right"/>

Tools for simulating mathematical models of infectious disease dynamics. Epidemic model classes include deterministic compartmental models, stochastic individual-contact models, and stochastic network models. Network models use the robust statistical methods of exponential-family random graph models (ERGMs) from the Statnet suite of software packages in R. Standard templates for epidemic modeling include SI, SIR, and SIS disease types. EpiModel features an easy API for extending these templates to address novel scientific research aims.

### Lead Authors

|                                                               |                                         |                          |
|---------------------------------------------------------------|-----------------------------------------|--------------------------|
| [Samuel M. Jenness](http://samueljenness.org/)                | Department of Epidemiology              | Emory University         |
| [Steven M. Goodreau](http://faculty.washington.edu/goodreau/) | Department of Anthropology              | University of Washington |
| [Martina Morris](http://faculty.washington.edu/morrism/)      | Departments of Statistics and Sociology | University of Washington |
| [Adrien Le Guillou](http://samueljenness.org/team.html)       | Department of Epidemiology              | Emory University         |

Additional contributors to EpiModel are listed on the [contributors](https://github.com/EpiModel/EpiModel/graphs/contributors) page.

### Installation

The current release version can be found on <a href="http://cran.r-project.org/web/packages/EpiModel/index.html" target="_blank">CRAN</a> and installed with:

``` r
install.packages("EpiModel", dependencies = TRUE)
```

To install this development version, use the <a href="https://github.com/r-lib/remotes" target="_blank">remotes package</a>:

``` r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("EpiModel/EpiModel")
```

### EpiModel 2.0

In July 2020, we released EpiModel version 2.0. This major software package update incorporates a substantial redesign of many elements of the EpiModel infrastructure and application programming interface (API). We anticipate that there will be some minor backwards incompatibilities with any EpiModel code developed with versions 1.x. There is a EpiModel 2.0 migration document available in the <a href="http://www.epimodel.org/tut.html" target="_blank"> Tutorials</a> page.

### Documentation and Support

**Website.** The main website for EpiModel, with tutorials and other supporting files is <a href="http://epimodel.org/" target="_blank">http://epimodel.org/</a>.

**Methods Paper.** A good place to start learning about EpiModel is the main methods paper published in the *Journal of Statistical Software.* It is available at <a href="http://doi.org/10.18637/jss.v084.i08" target="_blank">http://doi.org/10.18637/jss.v084.i08</a>.

**Summer Course.** Network Modeling for Epidemics is our annual 5-day course at the University of Washington where we teach the statistical theory, software tools, and applied modeling methods using EpiModel. <a href="http://statnet.github.io/nme/" target="_blank">Our course materials</a> are open-source and updated annually around the time of the course.

**Getting Help.** Users are encouraged to use Github issues on this repository as a place to ask questions (both technical coding quesitons and conceptual modeling questions), report bugs, and request new features & functionality.

### The EpiModel Gallery

We recently started a new <a href="https://github.com/EpiModel/EpiModel-Gallery" target="_blank">EpiModel Gallery</a> that contains templates of extensions to EpiModel, for now focused on network-based mathematical models. We will be continuing to add new examples the gallery, and encourage users to either file requests for new examples or contribute them following our guidelines.

### Citation

If using EpiModel for teaching or research, please include a citation our main methods paper:

> Jenness SM, Goodreau SM and Morris M. EpiModel: An R Package for Mathematical Modeling of Infectious Disease over Networks. *Journal of Statistical Software.* 2018; 84(8): 1-47. doi: 10.18637/jss.v084.i08

Please also <a href="mailto:samuel.m.jenness@emory.edu?Subject=We Used EpiModel in Our Study!" target="_top">send us an email </a> if you have used EpiModel in your work so we can add the citation below.

### Funding

The primary support for the development of these software tools and statistical methods has been by two National Institutes of Health (NIH) grants. Our applied research projects using EpiModel have received funding from the NIH and Centers for Disease Control and Prevention (CDC). Our team also receives institutional support through center-level NIH grants. A full list of our funding support can be found [here](https://github.com/EpiModel/EpiModel/wiki/EpiModel-Funding).

### EpiModel in the Scientific Literature

EpiModel and its [extension packages](https://github.com/EpiModel/EpiModelHIV) have been used in the following scientific journal articles. A list of these articles can be accessed in a [wiki page](https://github.com/EpiModel/EpiModel/wiki/EpiModel-in-the-Scientific-Literature) or on [Zotero](https://www.zotero.org/groups/2486200/epimodel_literature/library). (If you are aware of others, send us an email at [samuel.m.jenness\@emory.edu](mailto:samuel.m.jenness@emory.edu) to be included in this list.)

### Copyright

These materials are distributed under the GPL-3 license, with the following copyright and attribution requirements listed in the [LICENSE](https://github.com/EpiModel/EpiModel/blob/main/LICENSE.md) document above.
