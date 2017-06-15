EpiModel
===============

[![Version](http://img.shields.io/badge/Version-1.5.0-orange.svg?style=flat)](https://github.com/statnet/EpiModel/releases/tag/v1.5.0)
[![](http://cranlogs.r-pkg.org/badges/EpiModel?color=yellow)](http://cran.rstudio.com/web/packages/EpiModel/index.html)
[![](http://cranlogs.r-pkg.org/badges/grand-total/EpiModel?color=yellow)](http://cran.rstudio.com/web/packages/EpiModel/index.html)
[![Build Status](https://travis-ci.org/statnet/EpiModel.svg?branch=master)](https://travis-ci.org/statnet/EpiModel)
[![Vignette](https://img.shields.io/badge/docs-Vignette-943ad8.svg)](http://statnet.github.io/tut/EpiModelVignette.pdf)
<a href='http://dx.doi.org/10.5281/zenodo.16767' target="_blank"><img src='http://img.shields.io/badge/DOI-10.5281%2Fzenodo.16767-blue.svg?style=flat' alt='DOI' /></a>

<br>
<img align="right" src="http://www.epimodel.org/movie.gif">

EpiModel: tools for simulating mathematical models of infectious disease. Epidemic model classes include deterministic compartmental models, stochastic individual contact models, and stochastic network models. Disease types include SI, SIR, and SIS epidemics with and without demography, with tools available for expansion to model complex epidemic processes.


#### Installation
The current release version can be found on <a href="http://cran.r-project.org/web/packages/EpiModel/index.html" target="_blank">CRAN</a> and installed with:
```r
install.packages("EpiModel", dependencies = TRUE)
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

A good place to start learning about EpiModel is the main vignette, currently under review, but <a href="http://statnet.github.io/tut/EpiModelVignette.pdf" target="_blank">currently available in pre-press form here!</a>

#### Citation
If using EpiModel for teaching or research, please include a citation:
> Jenness SM, Goodreau SM, Morris M (2017). *EpiModel: Mathematical Modeling of Infectious Disease.* R Package Version 1.5.0. URL: http://epimodel.org/. DOI: 10.5281/zenodo.16767.

Please also send us an email if you have used EpiModel in your work.

#### Citing Articles

EpiModel has been used in the following scientific articles:

1. Jenness SM, Goodreau SM, Morris M, Cassels S. Effectiveness of Combination Packages for HIV-1 Prevention in Sub-Saharan Africa Depends on Partnership Network Structure. _Sexually Transmitted Infections._ 2016; 92(8): 619-624. [LINK](http://sti.bmj.com/content/early/2016/06/09/sextrans-2015-052476.abstract)

2. Jenness SM, Goodreau SM, Rosenberg E, Beylerian EN, Hoover KW, Smith DK, Sullivan P. Impact of CDC’s HIV Preexposure Prophylaxis Guidelines among MSM in the United States. _Journal of Infectious Diseases._ 2016; 214(12): 1800-1807. [LINK](http://jid.oxfordjournals.org/content/early/2016/07/12/infdis.jiw223.full)

3. Jenness SM, Sharma A, Goodreau SM, Rosenberg ES, Weiss KM, Hoover KW, Smith DK, Sullivan P. Individual HIV Risk versus Population Impact of Risk Compensation after HIV Preexposure Prophylaxis Initiation among Men Who Have Sex with Men. _PLoS One._ 2017; 12(1): e0169484. [LINK](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0169484)

4. Goodreau SM, Rosenberg ES, Jenness SM, Luisi N, Stansfield SE, Millett G, Sullivan P. Sources of Racial Disparities in HIV Prevalence among Men Who Have Sex with Men in Atlanta, GA: A Modeling Study. _Lancet HIV._ Epub ahead of print. DOI: 10.1016/ S2352-3018(17)30067. [LINK](https://www.ncbi.nlm.nih.gov/pubmed/28431923)

5. Jenness SM, Weiss KM, Goodreau SM, Rosenberg E, Gift T, Chesson H, Hoover KW, Smith DK, Liu AY, Sullivan P. Incidence of Gonorrhea and Chlamydia Following HIV Preexposure Prophylaxis among Men Who Have Sex with Men: A Modeling Study. _Clinical Infectious Diseases._ Epub ahead of print. DOI: 10.1093/cid/cix439. [LINK](https://academic.oup.com/cid/article-lookup/doi/10.1093/cid/cix439)

#### Funding
Development of this software has been supported by the following grants from the National Institutes of Health (NIH) and Centers for Disease Control and Prevention (CDC): 

* [NIH R01HD68395](https://projectreporter.nih.gov/project_info_description.cfm?aid=8841605): Statistical Methods for Network Epidemiology (PI: Martina Morris)
* [NIH R21HD075662](https://projectreporter.nih.gov/project_info_description.cfm?aid=8601779): Using Sexual Network Transmission Models to Explain HIV Disparities Between Black and White MSM (PI: Steven Goodreau)
* [NIH R01AI108490](https://projectreporter.nih.gov/project_info_description.cfm?aid=9024415): Integrated Bio-Social Models for HIV Epidemiology (MPIs: Steven Goodreau, Joshua Herbeck, and John Mittler)
* [NIH R21MH112449](https://projectreporter.nih.gov/project_info_description.cfm?aid=9271672): Modeling Antiretroviral-Based Prevention among MSM in the US (PI: Samuel Jenness)
* [NIH P30AI050409](https://projectreporter.nih.gov/project_info_description.cfm?aid=9120767): Center for AIDS Research at Emory University (MPIs: Carlos del Rio and James Curran)
* [NIH P30AI027757](https://projectreporter.nih.gov/project_info_description.cfm?aid=9069392): Center for AIDS Research at the University of Washington (PI: King Holmes)
* [CDC U38PS004646](https://projectreporter.nih.gov/project_info_details.cfm?aid=8926715): Enhancing Models of HIV, Viral Hepatitis, STIs, and Tuberculosis to Inform and Improve Public Health Impact (PI: Patrick Sullivan)
* [NIH P30DA027828](https://projectreporter.nih.gov/project_info_description.cfm?aid=9204281): Center for Prevention Implementation Methodology for Drug Abuse and HIV (Ce-PIM) (PI: Henricks Brown and Brian Mustanski)

#### Copyright
These materials are distributed under the GPL-3 license, with the following copyright and attribution requirements listed <a href="http://statnet.csde.washington.edu/attribution.shtml" target="_blank">here</a>.
