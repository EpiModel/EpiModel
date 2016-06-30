EpiModel
===============

[![Version](http://img.shields.io/badge/Version-1.2.6-orange.svg?style=flat)](https://github.com/statnet/EpiModel/releases/tag/v1.2.6)
[![](http://cranlogs.r-pkg.org/badges/EpiModel?color=yellow)](http://cran.rstudio.com/web/packages/EpiModel/index.html)
<a href='https://travis-ci.org/statnet/EpiModel' target="_blank"><img src='http://img.shields.io/travis/statnet/EpiModel/master.svg?style=flat' alt='Build Status' /></a>
<a href='https://coveralls.io/r/statnet/EpiModel?branch=master' target="_blank"><img src='https://coveralls.io/repos/statnet/EpiModel/badge.svg?branch=master' alt='Coverage Status' /></a>
<a href='http://dx.doi.org/10.5281/zenodo.16767' target="_blank"><img src='http://img.shields.io/badge/DOI-10.5281%2Fzenodo.16767-blue.svg?style=flat' alt='DOI' /></a>

<br>
<img align="right" src="http://www.epimodel.org/movie.gif">

EpiModel: tools for simulating mathematical models of infectious disease. Epidemic model classes include deterministic compartmental models, stochastic individual contact models, and stochastic network models. Disease types include SI, SIR, and SIS epidemics with and without demography, with tools available for expansion to model complex epidemic processes.


#### Installation
The current release version can be found on <a href="http://cran.r-project.org/web/packages/EpiModel/index.html" target="_blank">CRAN</a> and installed with:
```r
install.packages("EpiModel")
```

To install this development version, use the <a href="https://github.com/hadley/devtools" target="_blank">devtools package</a>:
```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("statnet/EpiModel")
```

#### Lead Authors
<table>
  <tr>
    <td><a href="http://samueljenness.org/" target="_blank">Samuel M. Jenness</a></th>
    <td>Department of Epidemiology</th>
    <td>Emory University</th>
  </tr>
  <tr>
    <td><a href="http://faculty.washington.edu/goodreau/" target="_blank">Steven M. Goodreau</a></td>
    <td>Department of Anthropology</td>
    <td>University of Washington</td>
  </tr>
  <tr>
    <td><a href="http://faculty.washington.edu/morrism/" target="_blank">Martina Morris</a></td>
    <td>Departments of Statistics and Sociology</td>
    <td>University of Washington</td>
  </tr>
</table>


#### Documentation
The main website for EpiModel, with tutorials and other supporting files is <a href="http://epimodel.org/" target="_blank">http://epimodel.org/</a>. Users are encouraged to join the <a href="http://mailman11.u.washington.edu/mailman/listinfo/epimodel" target="_blank">email list for EpiModel</a> as a place to ask questions, report bugs, and tell us about your research using these tools.

A good place to start learning about EpiModel is the main vignette, currently under review, but available in pre-press form <a href="http://statnet.github.io/tut/EpiModelVignette.pdf" target="_blank">here!</a>

#### Citation
If using EpiModel for teaching or research, please include a citation:
> Jenness SM, Goodreau SM, Morris M (2016). *EpiModel: Mathematical Modeling of Infectious Disease.* R Package Version 1.2.6. URL: http://epimodel.org/. DOI: 10.5281/zenodo.16767.

#### Funding
Development of this software is supported by the following grants from the National Institutes of Health: R01HD68395 (NICHD), T32HD007543 (NICHD), and R24HD042828 (NICHD).

#### Copyright
These materials are distributed under the GPL-3 license, with the following copyright and attribution requirements listed <a href="http://statnet.csde.washington.edu/attribution.shtml" target="_blank">here</a>.
