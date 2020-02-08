EpiModel
===============

[![Version](http://img.shields.io/badge/Version-1.7.5-orange.svg?style=flat)](https://github.com/statnet/EpiModel/releases/tag/v1.7.5)
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
* [NIH P30 DA027828](https://projectreporter.nih.gov/project_info_description.cfm?aid=9204281): Center for Prevention Implementation Methodology for Drug Abuse and HIV (Ce-PIM) (PI: Henricks Brown and Brian Mustanski)



### EpiModel in the Scientific Literature

EpiModel and its [extension packages](https://github.com/statnet/EpiModelHIV) have been used in the following scientific journal articles. (If you are aware of others, send us an email at samuel.m.jenness@emory.edu to be included in this list.)

#### HIV and Sexually Transmitted Infections

1. Delaney KP, Rosenberg ES, Kramer MR, Waller LA, Sullivan PS. Optimizing Human Immunodeficiency Virus Testing Interventions for Men Who Have Sex With Men in the United States: A Modeling Study. _Open Forum Infect Dis._ 2015;2(4): ofv153. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/26613096)

2. Jenness SM, Goodreau SM, Morris M, Cassels S. Effectiveness of Combination Packages for HIV-1 Prevention in Sub-Saharan Africa Depends on Partnership Network Structure. _Sexually Transmitted Infections._ 2016; 92(8): 619-624. [[LINK]](http://sti.bmj.com/content/early/2016/06/09/sextrans-2015-052476.abstract)

3. Jenness SM, Goodreau SM, Rosenberg E, Beylerian EN, Hoover KW, Smith DK, Sullivan P. Impact of CDC’s HIV Preexposure Prophylaxis Guidelines among MSM in the United States. _Journal of Infectious Diseases._ 2016; 214(12): 1800-1807. [[LINK]](http://jid.oxfordjournals.org/content/early/2016/07/12/infdis.jiw223.full)

4. Jenness SM, Sharma A, Goodreau SM, Rosenberg ES, Weiss KM, Hoover KW, Smith DK, Sullivan P. Individual HIV Risk versus Population Impact of Risk Compensation after HIV Preexposure Prophylaxis Initiation among Men Who Have Sex with Men. _PLoS One._ 2017; 12(1): e0169484. [[LINK]](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0169484)

5. Goodreau SM, Rosenberg ES, Jenness SM, Luisi N, Stansfield SE, Millett G, Sullivan P. Sources of Racial Disparities in HIV Prevalence among Men Who Have Sex with Men in Atlanta, GA: A Modeling Study. _Lancet HIV._ 2017; 4(7):e311-e320. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/28431923)

6. Jenness SM, Weiss KM, Goodreau SM, Rosenberg E, Gift T, Chesson H, Hoover KW, Smith DK, Liu AY, Sullivan P. Incidence of Gonorrhea and Chlamydia Following HIV Preexposure Prophylaxis among Men Who Have Sex with Men: A Modeling Study. _Clinical Infectious Diseases._ 2017; 65(5): 712-718. [[LINK]](https://academic.oup.com/cid/article-lookup/doi/10.1093/cid/cix439)

7. Vandormael A, Dobra A, Bärnighausen T, de Oliveira T, Tanser F. Incidence rate estimation, periodic testing and the limitations of the mid-point imputation approach. _International Journal of Epidemiology._ 2018; 47(1): 236-245. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/29024978)

8. Goodreau SM, Hamilton DT, Jenness SM, Sullivan PS, Valencia RK, Wang LY, Dunville RL, Barrios LC, Rosenberg ES. Targeting Human Immunodeficiency Virus Pre-Exposure Prophylaxis to Adolescent Sexual Minority Males in Higher Prevalence Areas of the United States: A Modeling Study. _J Adolesc Health._ 2018; 62(3): 311-319. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/29248392)

9. Herbeck JT, Peebles K, Edlefsen PT, Rolland M, Murphy JT, Gottlieb GS, Abernethy N, Mullins JI, Mittler JE, Goodreau SM. HIV population-level adaptation can rapidly diminish the impact of a partially effective vaccine. _Vaccine._ 2018;36(4): 514-520. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/29241646)

10. Luo W, Katz DA, Hamilton DT, McKenney J, Jenness SM, Goodreau SM, Stekler JD, Rosenberg ES, Sullivan P, Cassels S. Development of an Agent-Based Model to Investigate the Impact of HIV Self-Testing Programs for Men Who Have Sex with Men in Atlanta and Seattle. _Journal of Medical Internet Research Public Health Surveillance._ 2018; 4(2): e58. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/29959112)

11. Jenness SM, Maloney K, Smith SK, Hoover KW, Rosenberg ES, Goodreau SM, Weiss KM, Liu AY, Rao D, Sullivan PS. Addressing Gaps in HIV Preexposure Prophylaxis Care to Reduce Racial Disparities in HIV Incidence in the United States. _American Journal of Epidemiology._ 2019; 188(4): 743-752. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/30312365)

12. Stansfield SE, Mittler JE, Gottlieb GS, Murphy JT, Hamilton DT, Detels R, Wolinsky SM, Jacobson LP, Margolick JB, Rinaldo CR, Herbeck JT, Goodreau SM. Sexual Role and HIV-1 Set Point Viral Load among Men who Have Sex with Men. _Epidemics._ 2019; 26: 68-76. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/30193771)

13. Hamilton DT, Goodreau SM, Jenness SM, Sullivan PS, Wang LY, Dunville RL, Barrios LC, Rosenberg ES. Potential Impact of HIV Preexposure Prophylaxis Among Black and White Adolescent Sexual Minority Males: A Modeling Study. _American Journal of Public Health._ 2018; 108(S4): S284–S291. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/30383415)

14. Goodreau SM, Stansfield SE, Murphy JT, Peebles KC, Gottlieb GS, Abernethy NF, Herbeck JT, Mittler JE. Relational concurrency, stages of infection, and the evolution of HIV set point viral load. _Virus Evolution._ 2018; 4(2): vey032. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/30483403)

15. Goldstein ND, LeVasseur MT, Tran NK, Purtle J, Welles SL, Eppes SC. Modeling HPV vaccination scale-up among urban young men who have sex with men in the context of HIV. _Vaccine._ 2019; 37(29): 3883-3891. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/31155416)

16. Jones J, Weiss K, Mermin J, Dietz P, Rosenberg ES, Gift TL, Chesson H, Sullivan PS, Lyles C, Bernstein KT, Jenness SM. Proportion of Incident Human Immunodeficiency Virus Cases Among Men Who Have Sex With Men Attributable to Gonorrhea and Chlamydia: A Modeling Analysis. _Sexually Transmitted Diseases._ 2019; 46(6): 357-63. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/31095100)

17. Hamilton DT, Rosenberg ES, Jenness SM, Sullivan PS, Wang LY, Dunville RL, Barrios LC, Aslam M, Goodreau SM. Modeling the Joint Effects of Adolescent and Adult PrEP for Sexual Minority Males in the United States. _PloS One._ 2019; 14(5): e0217315. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/31116802)

18. Weiss KM, Jones J, Katz DA, Gift TL, Bernstein K, Workowski K, Rosenberg E, Jenness SM. Epidemiological Impact of Expedited Partner Therapy for Men Who Have Sex with Men: A Modeling Study. _Sexually Transmitted Diseases._ 2019; 46(11): 697–705. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/31644497)

19. Goldstein ND, LeVasseur MT, Tran NK, Purtle J, Welles SL, Eppes SC. Modeling HPV vaccination scale-up among urban young men who have sex with men in the context of HIV. _Vaccine._ 2019; 37(29): 3883-3891. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/31155416)

20. Weiss KM, Jones JS, Anderson EJ, Gift T, Chesson H, Bernstein K, Workowski K, Tuite A, Rosenberg ES, Sullivan PS, Jenness SM. Optimizing Coverage versus Frequency for Sexually Transmitted Infection Screening of Men Who Have Sex with Men. _Open Forum Infectious Diseases._ 2019; 6(10): ofz405. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/31667198)

21. Maloney KM, Driggers RA, Sarkar S, Anderson E, Malik A, Jenness SM. Expanded Choices with Greater HIV Prevention Benefits: A Mathematical Model of Long-Acting Injectable and Daily-Oral Pre-Exposure Prophylaxis. _medRxiv._ 2019; DOI: 10.1101/19012443. [[LINK]](https://doi.org/10.1101/19012443)

22. Wang LY, Hamilton DT, Rosenberg ES, Aslam MV, Sullivan PS, Katz DA, Dunville RL, Barrios LC, Goodreau SM. Cost-Effectiveness of Pre-Exposure Prophylaxis Among Adolescent Sexual Minority Males. _J Adolesc Health._ 2019; pii: S1054-139X(19)30415-X. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/31757626)

23. Mittler JE, Murphy JT, Stansfield SE, Peebles KC, Gottlieb GS, Abernethy NF, Reid MC, Goodreau SM, Herbeck JT. Large benefits to youth-focused HIV treatment-as-prevention efforts in generalized heterosexual populations: An agent-based simulation model. PLoS Comput Biol. 2019 Dec 17;15(12):e1007561. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/31846456)

#### Other Infectious Diseases and Contagious Processes

1. Ezenwa VO, Archie EA, Craft ME, Hawley DM, Martin LB, Moore J, White L. Host behaviour-parasite feedback: an essential link between animal behaviour and disease ecology. _Proc Biol Sci._ 2016; 283(1828). [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/27053751)

2. Webber QM, Brigham RM, Park AD, Gillam EH, O’Shea TJ, Willis CK. Social network characteristics and predicted pathogen transmission in summer colonies of female big brown bats (Eptesicus fuscus). _Behavioral Ecology and Sociobiology._ 2016;70(5): 701-12. [[LINK]](https://link.springer.com/article/10.1007/s00265-016-2093-3). 

3. Goldstein ND, Eppes SC, Mackley A, Tuttle D, Paul DA. A Network Model of Hand Hygiene: How Good Is Good Enough to Stop the Spread of MRSA? _Infect Control Hosp Epidemiol._ 2017; 38(8): 945-52. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/28656884)

4. White LA, Forester JD, Craft ME. Covariation between the physiological and behavioral components of pathogen transmission: Host heterogeneity determines epidemic outcomes. _Oikos._ 2018; 127(4): 538-52. [[LINK]](http://onlinelibrary.wiley.com/doi/10.1111/oik.04527/full).

5. Robinson SJ, Barbieri MM, Murphy S, Baker JD, Harting AL, Craft ME, Littnan CL. Model recommendations meet management reality: implementation and evaluation of a network-informed vaccination effort for endangered Hawaiian monk seals. _Proceeding of the Royal Society B._ 2018; 285(1870): 20171899. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/29321294).

6. Goldstein ND, Jenness SM, Tuttle D, Power M, Paul DA, Eppes SC. Evaluating a neonatal intensive care unit HRSA surveillance programme using agent-based network modeling. _Journal of Hospital Infection._ 2018; 100(3): 337-43.  [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/29751022)

7. Haeussler K, Hout AV, Baio G. A dynamic Bayesian Markov model for health economic evaluations of interventions against infectious diseases. _arXiv._ arXiv:1512.06881. [[LINK]](https://arxiv.org/abs/1512.06881). 

8. Amirpour Haredasht S, Tavornpanich S, Jansen MD, Lyngstad TM, Yatabe T, Brun E, Martínez-López B. A stochastic network-based model to simulate the spread of pancreas disease (PD) in the Norwegian salmon industry based on the observed vessel movements and seaway distance between marine farms. _Prev Vet Med._ 2019; 167: 174-181. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/30055856)

9. Wilson-Aggarwal JK, Ozella L, Tizzoni M, Cattuto C, Swan GJ, Moundai T, Silk MJ, Zingeser JA, McDonald RA. High-resolution contact networks of free-ranging domestic dogs Canis familiaris and implications for transmission of infection. _PLoS Neglected Tropical Diseases._ 2019; 13(7): e0007565. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/31306425)

10. Baker E, Challenor P, Eames M. Predicting the Output From a Stochastic Computer Model When a Deterministic Approximation is Available. _arXiv._ 2019; 1902.01290. [[LINK]](https://arxiv.org/abs/1902.01290)

11. Milwid RM, O'Sullivan TL, Poljak Z, Laskowski M, Greer AL. Comparing the effects of non-homogenous mixing patterns on epidemiological outcomes in equine populations: A mathematical modelling study. _Sci Rep._ 2019; 9(1): 3227. [[LINK]](https://www.ncbi.nlm.nih.gov/pubmed/30824806)

12. Everton SF, Schroeder R. Plagues, Pagans, and Christians: Differential Survival, Social Networks, and the Rise of Christianity. _Journal for the Scientific Study of Religion._ 2019; DOI: 10.1111/jssr.12631. [[LINK]](https://doi.org/10.1111/jssr.12631)

13. Amusan O, Thompson AF, Aderinola TB, Alese BK. Modelling Malicious Attack in Social Networks. _Network and Communication Technologies._ 2020; 5(1): 37-43. [[LINK]](http://www.ccsenet.org/journal/index.php/nct/article/view/0/41983)


### Copyright
These materials are distributed under the GPL-3 license, with the following copyright and attribution requirements listed in the [LICENSE](https://github.com/statnet/EpiModel/blob/master/LICENSE.md) document above.
