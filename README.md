EpiModel
===============

[![Version](http://img.shields.io/badge/Version-1.1.1-orange.svg?style=flat)](https://github.com/statnet/EpiModel/releases/tag/v1.1.1)
[![Build Status](http://img.shields.io/travis/statnet/EpiModel/master.svg?style=flat)](https://travis-ci.org/statnet/EpiModel)
[![DOI](http://img.shields.io/badge/DOI-10.5281%2Fzenodo.12524-blue.svg?style=flat)](http://dx.doi.org/10.5281/zenodo.12524)

EpiModel: tools for simulating mathematical models of infectious disease. Epidemic model classes include deterministic compartmental models, stochastic individual contact models, and stochastic network models. Disease types include SI, SIR, and SIS epidemics with and without demography, with tools available for expansion to model complex epidemic processes.


### Installation
The current release version can be found on <a href="http://cran.r-project.org/web/packages/EpiModel/index.html" target="_blank">CRAN</a> and installed with:
```r
install.packages("EpiModel")
```

To install this development version, use the <a href="https://github.com/hadley/devtools" target="_blank">devtools package</a>:
```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("statnet/EpiModel")
```

### Authors
<table>
  <tr>
    <td><a href="http://samueljenness.org/" target="_blank">Samuel M. Jenness</a></th>
    <td>Department of Epidemiology</th>
    <td>University of Washington</th>
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


### Documentation
The main website for EpiModel, with tutorials and other supporting files, is <a href="http://statnet.github.io/EpiModel/" target="_blank">here</a>. Users are encouraged to join the <a href="http://mailman11.u.washington.edu/mailman/listinfo/epimodel" target="_blank">email list for EpiModel</a> as a place to ask questions, report bugs, and tell us about your research using these tools.

### Other 

#### Citation
If using EpiModel for teaching or research, please include a citation:
> Jenness SM, Goodreau SM, Morris M (2014). *EpiModel: Mathematical Modeling of Infectious Disease.* R Package Version 1.1.1. URL: http://epimodel.org/. DOI: 10.5281/zenodo.12524.

#### Funding
This software is supported by grant number R01HD68395 from the National Institute of Child Health and Human Development.

#### Copyright
These materials are distributed under the GPL-3 license, with the following copyright and attribution requirements listed <a href="http://statnet.csde.washington.edu/attribution.shtml" target="_blank">here</a>.
