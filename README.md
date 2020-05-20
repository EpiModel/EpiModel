EpiModel
===============

[![Version](http://img.shields.io/badge/Version-1.8.0-orange.svg?style=flat)](https://github.com/statnet/EpiModel/releases/tag/v1.8.0)
[![](http://cranlogs.r-pkg.org/badges/EpiModel?color=yellow)](http://cran.rstudio.com/web/packages/EpiModel/index.html)
[![](http://cranlogs.r-pkg.org/badges/grand-total/EpiModel?color=blue)](http://cran.rstudio.com/web/packages/EpiModel/index.html)
[![Build Status](https://travis-ci.org/statnet/EpiModel.svg?branch=master)](https://travis-ci.org/statnet/EpiModel)
[![Methods](https://img.shields.io/badge/docs-Methods-943ad8.svg)](http://doi.org/10.18637/jss.v084.i08)

<br>
<img align="right" src="http://www.epimodel.org/movie.gif">

Tools for simulating mathematical models of infectious disease dynamics. Epidemic model classes include deterministic compartmental models, stochastic individual-contact models, and stochastic network models. Network models use the robust statistical methods of exponential-family random graph models (ERGMs) from the Statnet suite of software packages in R. Standard templates for epidemic modeling include SI, SIR, and SIS disease types. EpiModel features an easy API for extending these templates to address novel scientific research aims.

### Lead Authors
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

Additional contributors to EpiModel are listed on the [contributors](https://github.com/statnet/EpiModel/graphs/contributors) page.


### Installation
The current release version can be found on <a href="http://cran.r-project.org/web/packages/EpiModel/index.html" target="_blank">CRAN</a> and installed with:
```r
install.packages("EpiModel", dependencies = TRUE)
```

To install this development version, use the <a href="https://github.com/r-lib/remotes" target="_blank">remotes package</a>:
```r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("statnet/EpiModel")
```

### Documentation and Support

**Website.** The main website for EpiModel, with tutorials and other supporting files is <a href="http://epimodel.org/" target="_blank">http://epimodel.org/</a>.

**Methods Paper.** A good place to start learning about EpiModel is the main methods paper published in the *Journal of Statistical Software.* It is available at <a href="http://doi.org/10.18637/jss.v084.i08" target="_blank">http://doi.org/10.18637/jss.v084.i08</a>.

**Summer Course.** Network Modeling for Epidemics is our annual 5-day course at the University of Washington where we teach the statistical theory, software tools, and applied modeling methods using EpiModel. <a href="http://statnet.github.io/nme/" target="_blank">Our course materials</a> are open-source and updated annually around the time of the course.

**Email listserv.** Users are encouraged to join the <a href="http://mailman11.u.washington.edu/mailman/listinfo/epimodel" target="_blank">email list for EpiModel</a> as a place to ask questions, report bugs, and tell us about your research using these tools.


### The EpiModel Gallery
We recently started a new <a href="https://github.com/statnet/EpiModel-Gallery" target="_blank">EpiModel Gallery</a> that contains templates of extensions to EpiModel, for now focused on network-based mathematical models. We will be continuing to add new examples the gallery, and encourage users to either file requests for new examples or contribute them following our guidelines.

### Citation
If using EpiModel for teaching or research, please include a citation our main methods paper:

> Jenness SM, Goodreau SM and Morris M. EpiModel: An R Package for Mathematical Modeling of Infectious Disease over Networks. *Journal of Statistical Software.* 2018; 84(8): 1-47. doi: 10.18637/jss.v084.i08

Please also <a href="mailto:samuel.m.jenness@emory.edu?Subject=We Used EpiModel in Our Study!" target="_top">send us an email </a> if you have used EpiModel in your work so we can add the citation below.

### Funding
The primary support for the development of these software tools and statistical methods has been by  two National Institutes of Health (NIH) grants:

* [NIH R01 AI138783](https://projectreporter.nih.gov/project_info_description.cfm?aid=9623724): EpiModel 2.0: Integrated Network Models for HIV/STI Prevention Science (PI: Samuel Jenness)
* [NIH R01 HD68395](https://projectreporter.nih.gov/project_info_description.cfm?aid=8841605): Statistical Methods for Network Epidemiology (PI: Martina Morris)

Our applied research projects using EpiModel have received funding from the NIH and Centers for Disease Control and Prevention (CDC):

* [NIH R21 MH112449](https://projectreporter.nih.gov/project_info_description.cfm?aid=9271672): Modeling Antiretroviral-Based Prevention among MSM in the US (PI: Samuel Jenness)
* [NIH R21 HD075662](https://projectreporter.nih.gov/project_info_description.cfm?aid=8601779): Using Sexual Network Transmission Models to Explain HIV Disparities Between Black and White MSM (PI: Steven Goodreau)
* [NIH R01 AI108490](https://projectreporter.nih.gov/project_info_description.cfm?aid=9024415): Integrated Bio-Social Models for HIV Epidemiology (MPIs: Steven Goodreau, Joshua Herbeck, and John Mittler)
* [CDC U38 PS004646](https://projectreporter.nih.gov/project_info_details.cfm?aid=8926715): Enhancing Models of HIV, Viral Hepatitis, STIs, and Tuberculosis to Inform and Improve Public Health Impact (PI: Patrick Sullivan)

Our team also receives institutional support through the following center-level NIH grants:

* [NIH P30 AI050409](https://projectreporter.nih.gov/project_info_description.cfm?aid=9120767): Center for AIDS Research at Emory University (MPIs: Carlos del Rio and James Curran)
* [NIH P30 AI027757](https://projectreporter.nih.gov/project_info_description.cfm?aid=9069392): Center for AIDS Research at the University of Washington (PI: King Holmes)


### EpiModel in the Scientific Literature

EpiModel and its [extension packages](https://github.com/statnet/EpiModelHIV) have been used in the following scientific journal articles. A list of these articles can be accessed in a [wiki page](https://github.com/statnet/EpiModel/wiki/EpiModel-in-the-Scientific-Literature) or on [Zotero](https://www.zotero.org/groups/2486200/epimodel_literature/library). (If you are aware of others, send us an email at samuel.m.jenness@emory.edu to be included in this list.)


### Copyright
These materials are distributed under the GPL-3 license, with the following copyright and attribution requirements listed in the [LICENSE](https://github.com/statnet/EpiModel/blob/master/LICENSE.md) document above.
