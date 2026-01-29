# EpiModel

[![CRAN
Version](http://img.shields.io/badge/Version-2.5.0-orange.svg?style=flat)](https://github.com/EpiModel/EpiModel/releases/tag/v2.5.0)
[![](http://cranlogs.r-pkg.org/badges/EpiModel?color=blue)](https://CRAN.R-project.org/package=EpiModel)
[![](http://cranlogs.r-pkg.org/badges/grand-total/EpiModel?color=blue)](https://CRAN.R-project.org/package=EpiModel)
[![Build
Status](https://github.com/EpiModel/EpiModel/workflows/R-CMD-check/badge.svg)](https://github.com/EpiModel/EpiModel/actions)
[![Methods](https://img.shields.io/badge/docs-Methods-943ad8.svg)](https://doi.org/10.18637/jss.v084.i08)

  
![](https://www.epimodel.org/movie.gif)

Tools for simulating mathematical models of infectious disease dynamics.
Epidemic model classes include deterministic compartmental models,
stochastic individual-contact models, and stochastic network models.
Network models use the robust statistical methods of exponential-family
random graph models (ERGMs) from the Statnet suite of software packages
in R. Standard templates for epidemic modeling include SI, SIR, and SIS
disease types. EpiModel features an easy application programming
interface (API) for extending these templates to address novel
scientific research aims.

### Lead Authors

|                                                               |                                              |                          |
|---------------------------------------------------------------|----------------------------------------------|--------------------------|
| [Samuel M. Jenness](http://samueljenness.org/)                | Department of Epidemiology                   | Emory University         |
| [Steven M. Goodreau](http://faculty.washington.edu/goodreau/) | Department of Anthropology                   | University of Washington |
| [Martina Morris](http://faculty.washington.edu/morrism/)      | Departments of Statistics and Sociology      | University of Washington |
| [Adrien Le Guillou](http://samueljenness.org/team.md)         | Department of Epidemiology                   | Emory University         |
| Chad Klumb                                                    | Center for Studies in Demography and Ecology | University of Washington |

Additional contributors to EpiModel are listed on the
[contributors](https://github.com/EpiModel/EpiModel/graphs/contributors)
page.

### Installation

The current release version can be found on
[CRAN](https://CRAN.R-project.org/package=EpiModel) and installed with:

``` r
install.packages("EpiModel", dependencies = TRUE)
```

To install this development version, use the [remotes
package](https://github.com/r-lib/remotes):

``` r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("EpiModel/EpiModel")
```

### Documentation and Support

**Website.** The main website for EpiModel, with tutorials and other
supporting files, is here:
[http://www.epimodel.org/](https://www.epimodel.org/).

**Methods Paper.** A good place to start learning about EpiModel is the
main methods paper published in the *Journal of Statistical Software.*
It is available at:
[https://doi.org/10.18637/jss.v084.i08](https://doi.org/10.18637/jss.v084.i08).

**Summer Course.** Network Modeling for Epidemics is our annual summer
workshop/course where we teach the statistical theory, software tools,
and applied modeling methods using EpiModel. It is now part of the
Summer Institute in Statistics and Modeling in Infectious Diseases
(SISMID) at Emory University. [Our course
materials](https://epimodel.github.io/sismid/) are open-source and
updated annually around the time of the course.

**Getting Help.** Users are encouraged to use Github issues on this
repository as a place to ask questions (both technical coding questions
and conceptual modeling questions), report bugs, and request new
features & functionality. Broader modeling questions can be posted on
the Discussions board here.

### The EpiModel Gallery

The [EpiModel Gallery](https://github.com/EpiModel/EpiModel-Gallery)
contains templates of extensions to EpiModel, for now focused on
network-based mathematical modeling class. We will be continuing to add
new examples the gallery, and encourage users to either file requests
for new examples or else to contribute them directly.

### Citation

If using EpiModel for teaching or research, please include a citation
our main methods paper:

> Jenness SM, Goodreau SM and Morris M. EpiModel: An R Package for
> Mathematical Modeling of Infectious Disease over Networks. *Journal of
> Statistical Software.* 2018; 84(8): 1-47. doi: 10.18637/jss.v084.i08

Please also [send us an
email](mailto:samuel.m.jenness@emory.edu?Subject=We%20Used%20EpiModel%20in%20Our%20Study!)
if you have used EpiModel in your work so we can add the citation below.

### Funding

The primary support for the development of these software tools and
statistical methods has been by two National Institutes of Health (NIH)
grants. Our applied research projects using EpiModel have received
funding from the NIH and Centers for Disease Control and Prevention
(CDC). Our team also receives institutional support through center-level
NIH grants. A full list of our funding support can be found
[here](https://github.com/EpiModel/EpiModel/wiki/EpiModel-Funding).

### EpiModel in the Scientific Literature

EpiModel and its [extension
packages](https://github.com/EpiModel/EpiModelHIV) have been used in the
following scientific journal articles. A list of these articles can be
accessed in a [wiki
page](https://github.com/EpiModel/EpiModel/wiki/EpiModel-in-the-Scientific-Literature)
or on
[Zotero](https://www.zotero.org/groups/2486200/epimodel_literature/library).
(If you are aware of others, send us an email at
<samuel.m.jenness@emory.edu> to be included in this list.)

### Copyright

These materials are distributed under the GPL-3 license, with the
following copyright and attribution requirements listed in the
[LICENSE](https://github.com/EpiModel/EpiModel/blob/main/LICENSE.md)
document above.
