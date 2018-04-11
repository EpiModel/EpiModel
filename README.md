EpiModel
===============

[![Version](http://img.shields.io/badge/Version-1.6.5-orange.svg?style=flat)](https://github.com/statnet/EpiModel/releases/tag/v1.6.5)
[![](http://cranlogs.r-pkg.org/badges/EpiModel?color=yellow)](http://cran.rstudio.com/web/packages/EpiModel/index.html)
[![](http://cranlogs.r-pkg.org/badges/grand-total/EpiModel?color=yellow)](http://cran.rstudio.com/web/packages/EpiModel/index.html)
[![Build Status](https://travis-ci.org/statnet/EpiModel.svg?branch=master)](https://travis-ci.org/statnet/EpiModel)
[![Methods](https://img.shields.io/badge/docs-Methods-943ad8.svg)](https://www.biorxiv.org/content/early/2017/11/03/213009)

<br>
<img align="right" src="http://www.epimodel.org/movie.gif">

Tools for simulating mathematical models of infectious disease dynamics. Epidemic model classes include deterministic compartmental models, stochastic agent-based models, and stochastic network models. Network models use the robust statistical methods of exponential-family random graph models (ERGMs) from the Statnet suite of software packages in R. Standard templates for epidemic modeling include SI, SIR, and SIS disease types. EpiModel features an easy API for extending these templates to address novel scientific research aims.


### Installation
The current release version can be found on <a href="http://cran.r-project.org/web/packages/EpiModel/index.html" target="_blank">CRAN</a> and installed with:
```r
install.packages("EpiModel", dependencies = TRUE)
```

To install this development version, use the <a href="https://github.com/hadley/devtools" target="_blank">devtools package</a>:
```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("statnet/EpiModel")
```

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


### Documentation
The main website for EpiModel, with tutorials and other supporting files is <a href="http://epimodel.org/" target="_blank">http://epimodel.org/</a>. Users are encouraged to join the <a href="http://mailman11.u.washington.edu/mailman/listinfo/epimodel" target="_blank">email list for EpiModel</a> as a place to ask questions, report bugs, and tell us about your research using these tools.

A good place to start learning about EpiModel is the main methods paper. It is currently in-press at the *Journal of Statistical Software* but available as <a href="https://www.biorxiv.org/content/early/2017/11/03/213009" target="_blank">pre-press at bioRxiv here</a>.

### Citation
If using EpiModel for teaching or research, please include a citation:
> Jenness SM, Goodreau SM, Morris M. EpiModel: An R Package for Mathematical Modeling of Infectious Disease over Networks. *bioRxiv* 2017;  213009. DOI: https://doi.org/10.1101/213009.

Please also send us an email if you have used EpiModel in your work so we can add the citation below.

### Applied Uses of EpiModel in the Scientific Literature

EpiModel has been used in the following scientific journal articles.

#### HIV and Other Sexually Transmitted Infections

1. Delaney KP, Rosenberg ES, Kramer MR, Waller LA, Sullivan PS. Optimizing Human Immunodeficiency Virus Testing Interventions for Men Who Have Sex With Men in the United States: A Modeling Study. _Open Forum Infect Dis._ 2015;2(4): ofv153. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/26613096)

2. Jenness SM, Goodreau SM, Morris M, Cassels S. Effectiveness of Combination Packages for HIV-1 Prevention in Sub-Saharan Africa Depends on Partnership Network Structure. _Sexually Transmitted Infections._ 2016; 92(8): 619-624. [[LINK]](http://sti.bmj.com/content/early/2016/06/09/sextrans-2015-052476.abstract)

3. Jenness SM, Goodreau SM, Rosenberg E, Beylerian EN, Hoover KW, Smith DK, Sullivan P. Impact of CDC’s HIV Preexposure Prophylaxis Guidelines among MSM in the United States. _Journal of Infectious Diseases._ 2016; 214(12): 1800-1807. [[LINK]](http://jid.oxfordjournals.org/content/early/2016/07/12/infdis.jiw223.full)

4. Jenness SM, Sharma A, Goodreau SM, Rosenberg ES, Weiss KM, Hoover KW, Smith DK, Sullivan P. Individual HIV Risk versus Population Impact of Risk Compensation after HIV Preexposure Prophylaxis Initiation among Men Who Have Sex with Men. _PLoS One._ 2017; 12(1): e0169484. [[LINK]](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0169484)

5. Goodreau SM, Rosenberg ES, Jenness SM, Luisi N, Stansfield SE, Millett G, Sullivan P. Sources of Racial Disparities in HIV Prevalence among Men Who Have Sex with Men in Atlanta, GA: A Modeling Study. _Lancet HIV._ 2017; 4(7):e311-e320. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/28431923)

6. Jenness SM, Weiss KM, Goodreau SM, Rosenberg E, Gift T, Chesson H, Hoover KW, Smith DK, Liu AY, Sullivan P. Incidence of Gonorrhea and Chlamydia Following HIV Preexposure Prophylaxis among Men Who Have Sex with Men: A Modeling Study. _Clinical Infectious Diseases._ 2017; 65(5): 712-718. [[LINK]](https://academic.oup.com/cid/article-lookup/doi/10.1093/cid/cix439)

7. Vandormael A, Dobra A, Bärnighausen T, de Oliveira T, Tanser F. Incidence rate estimation, periodic testing and the limitations of the mid-point imputation approach. _International Journal of Epidemiology._ Epub 2017 Aug 9. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/29024978)

8. Goodreau SM, Hamilton DT, Jenness SM, Sullivan PS, Valencia RK, Wang LY, Dunville RL, Barrios LC, Rosenberg ES. Targeting Human Immunodeficiency Virus Pre-Exposure Prophylaxis to Adolescent Sexual Minority Males in Higher Prevalence Areas of the United States: A Modeling Study. _J Adolesc Health._ Epub 2017. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/29248392)

9. Herbeck JT, Peebles K, Edlefsen PT, Rolland M, Murphy JT, Gottlieb GS, Abernethy N, Mullins JI, Mittler JE, Goodreau SM. HIV population-level adaptation can rapidly diminish the impact of a partially effective vaccine. _Vaccine._ 2018;36(4): 514-520. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/29241646)


#### Other

1. Ezenwa VO, Archie EA, Craft ME, Hawley DM, Martin LB, Moore J, White L. Host behaviour-parasite feedback: an essential link between animal behaviour and disease ecology. Proc Biol Sci. 2016; 283(1828). [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/27053751)

2. Goldstein ND, Eppes SC, Mackley A, Tuttle D, Paul DA. A Network Model of Hand Hygiene: How Good Is Good Enough to Stop the Spread of MRSA? _Infect Control Hosp Epidemiol._ 2017:1-8. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/28656884)

3. White LA, Forester JD, Craft ME. Covariation between the physiological and behavioral components of pathogen transmission: Host heterogeneity determines epidemic outcomes. _Oikos._ Epub 2017 November. [[LINK]](http://onlinelibrary.wiley.com/doi/10.1111/oik.04527/full).

4. Haeussler K, Hout AV, Baio G. A dynamic Bayesian Markov model for health economic evaluations of interventions against infectious diseases. _arXiv._ arXiv:1512.06881. [[LINK]](https://arxiv.org/abs/1512.06881). 

5. Webber QM, Brigham RM, Park AD, Gillam EH, O’Shea TJ, Willis CK. Social network characteristics and predicted pathogen transmission in summer colonies of female big brown bats (Eptesicus fuscus). _Behavioral Ecology and Sociobiology._ 2016;70(5): 701-12. [[LINK]](https://link.springer.com/article/10.1007/s00265-016-2093-3). 

6. Robinson SJ, Barbieri MM, Murphy S, Baker JD, Harting AL, Craft ME, Littnan CL. Model recommendations meet management reality: implementation and evaluation of a network-informed vaccination effort for endangered Hawaiian monk seals. _Proceeding of the Royal Society B._ 2018; DOI: 10.1098/rspb.2017.1899. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/29321294).

### Funding
Development of this software has been supported by the following grants from the National Institutes of Health (NIH) and Centers for Disease Control and Prevention (CDC): 

* [NIH R01HD68395](https://projectreporter.nih.gov/project_info_description.cfm?aid=8841605): Statistical Methods for Network Epidemiology (PI: Martina Morris)
* [NIH R21HD075662](https://projectreporter.nih.gov/project_info_description.cfm?aid=8601779): Using Sexual Network Transmission Models to Explain HIV Disparities Between Black and White MSM (PI: Steven Goodreau)
* [NIH R01AI108490](https://projectreporter.nih.gov/project_info_description.cfm?aid=9024415): Integrated Bio-Social Models for HIV Epidemiology (MPIs: Steven Goodreau, Joshua Herbeck, and John Mittler)
* [NIH R21MH112449](https://projectreporter.nih.gov/project_info_description.cfm?aid=9271672): Modeling Antiretroviral-Based Prevention among MSM in the US (PI: Samuel Jenness)
* [NIH P30AI050409](https://projectreporter.nih.gov/project_info_description.cfm?aid=9120767): Center for AIDS Research at Emory University (MPIs: Carlos del Rio and James Curran)
* [NIH P30AI027757](https://projectreporter.nih.gov/project_info_description.cfm?aid=9069392): Center for AIDS Research at the University of Washington (PI: King Holmes)
* [CDC U38PS004646](https://projectreporter.nih.gov/project_info_details.cfm?aid=8926715): Enhancing Models of HIV, Viral Hepatitis, STIs, and Tuberculosis to Inform and Improve Public Health Impact (PI: Patrick Sullivan)
* [NIH P30DA027828](https://projectreporter.nih.gov/project_info_description.cfm?aid=9204281): Center for Prevention Implementation Methodology for Drug Abuse and HIV (Ce-PIM) (PI: Henricks Brown and Brian Mustanski)

### Copyright
These materials are distributed under the GPL-3 license, with the following copyright and attribution requirements listed <a href="http://statnet.csde.washington.edu/attribution.shtml" target="_blank">here</a>.
