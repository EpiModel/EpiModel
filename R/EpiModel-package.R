
#' Mathematical Modeling of Infectious Disease Dynamics
#'
#' @details
#' EpiModel provides tools for building, simulating, and analyzing mathematical
#' models of infectious disease dynamics. Models may be built for pedagogical
#' purposes using the built-in disease types, or for research purposes using the
#' extension API to define arbitrary disease states and processes.
#'
#' The distinguishing capability of EpiModel is stochastic network modeling,
#' in which partnerships form, persist, and dissolve according to a statistical
#' model estimated from network data. Simpler compartmental and individual-based
#' classes are also provided, both for teaching and for comparison against the
#' network results.
#'
#' @section Model Classes and Infectious Disease Types:
#' EpiModel provides functionality for three classes of epidemic models:
#'
#'  * **Deterministic Compartmental Models (DCMs):** continuous-time models
#'        solved with ordinary differential equations, under an assumption of
#'        random mixing within the population. Sensitivity analyses across
#'        parameter values or initial conditions are specified by passing
#'        vectors to [param.dcm()] or [init.dcm()].
#'  * **Stochastic Individual Contact Models (ICMs):** discrete-time,
#'        individual-based analogs of the DCM class, which add random variation
#'        to each component of the transmission system, from infection to
#'        recovery to vital dynamics (arrivals and departures). Contacts are
#'        randomly mixed, with no persistent partnership structure.
#'  * **Stochastic Network Models:** discrete-time, individual-based models in
#'        which the contact structure is represented explicitly as a dynamic
#'        network. Partnership formation and dissolution are simulated from a
#'        temporal exponential-family random graph model (TERGM), using the
#'        statistical framework implemented in the **Statnet** suite of R
#'        packages. Disease then spreads across the edges of that evolving
#'        network. This class represents features that random mixing cannot,
#'        including degree distributions, partnership duration, concurrency,
#'        assortative mixing, and feedback between network structure and
#'        epidemic dynamics.
#'
#'
#' EpiModel supports three built-in infectious disease types, available across
#' all three model classes, in one-group and two-group configurations and with
#' or without demography:
#'
#'   * **Susceptible-Infectious (SI):** a two-state disease with life-long
#'         infection and no recovery. HIV is one example, although for that case
#'         it is common to model infection stages as separate compartments.
#'   * **Susceptible-Infectious-Recovered (SIR):** a three-state disease with
#'         life-long immunity following infection. Measles is one example,
#'         although modern models for measles also require consideration of
#'         vaccination patterns in the population.
#'   * **Susceptible-Infectious-Susceptible (SIS):** a two-state disease in
#'         which one may transition back and forth between the susceptible and
#'         infected states throughout life. Examples include bacterial sexually
#'         transmitted infections such as gonorrhea.
#'
#' These built-in types are starting points. The extension API described below
#' allows them to be replaced with arbitrarily complex disease systems.
#'
#' @section Model Parameterization and Simulation:
#' Each model class uses three setup functions to collect the epidemic
#' parameters, initial conditions, and control settings:
#'
#'  * [param.dcm()], [param.icm()], and [param.net()] set the epidemic
#'        parameters, such as the rate of contacts or acts between agents, the
#'        probability of transmission per act, and the recovery and demographic
#'        rates for models that include those transitions.
#'  * [init.dcm()], [init.icm()], and [init.net()] set the initial
#'        conditions, which for built-in types are the numbers or, where
#'        applicable, the specific agents infected or recovered at the outset.
#'  * [control.dcm()], [control.icm()], and [control.net()] set the remaining
#'        simulation controls, including the disease type, the number of time
#'        steps, and the number of simulations or runs.
#'
#'
#' With the model parameterized, the simulation functions are:
#'
#' * [dcm()] for deterministic compartmental models.
#' * [icm()] for stochastic individual contact models.
#' * Network models are simulated in a three-step process:
#'    1. [netest()] estimates the statistical model for the network structure
#'       itself, that is, how partnerships form and dissolve over time. This
#'       function wraps the estimation routines in the `ergm` and `tergm`
#'       packages. The typical input is a set of target statistics summarizing
#'       the formation process together with a mean partnership duration,
#'       an approach known as egocentric inference because those targets may be
#'       calculated from an egocentric sample of the population. Dissolution
#'       coefficients are calculated with [dissolution_coefs()].
#'    2. [netdx()] runs diagnostics on the fitted model by simulating the
#'       network forward in time without disease, to check that the simulated
#'       network recovers the formation and dissolution targets.
#'    3. [netsim()] simulates the stochastic network epidemic model, given the
#'       fitted network model from [netest()] along with the parameters,
#'       initial conditions, and controls defined above.
#'
#'
#' Network models may also use multiple network layers over a shared node set,
#' for example separate sexual and needle-sharing contact networks. Each layer
#' is estimated separately and the fits are passed together to [netsim()], with
#' [multilayer()] specifying any controls that vary by layer.
#'
#' @section Extending the Built-In Models:
#' Research applications typically require disease states or processes beyond
#' the built-in types. The extension API supports this for the network class:
#' setting `type = NULL` in [control.net()] and supplying module functions
#' (for example, `infection.FUN = my_infection_module`) replaces or supplements
#' the built-in modules. A module is a function with the signature
#' `function(dat, at)` that reads from and returns the `netsim_dat` object
#' holding all simulation state.
#'
#' Modules should use the accessor functions rather than manipulating that
#' object directly. These include [get_attr()] and [set_attr()] for nodal
#' attributes, [get_param()] for parameters, [get_epi()] and [set_epi()] for
#' epidemic summary statistics, [append_core_attr()] and [append_attr()] for
#' new arrivals, [discord_edgelist()] for serodiscordant pairs, and
#' [set_transmat()] for recording transmission events.
#'
#' Extension support is specific to the network class. ICMs support the
#' built-in SI, SIR, and SIS types only; users needing custom individual-based
#' dynamics should use [netsim()]. DCMs accept an original system of
#' differential equations through the `new.mod` argument to [control.dcm()].
#'
#' @section Working with Simulation Output:
#' All three model classes return objects with `print`, `summary`, `plot`, and
#' `as.data.frame` methods, so the same analysis workflow applies throughout.
#' Simulation output may be converted to a long data frame with
#' [as.data.frame.dcm()], [as.data.frame.icm()], or [as.data.frame.netsim()],
#' extended with derived statistics using [mutate_epi()], subset with
#' [get_sims()], and combined across runs with [merge.icm()] or
#' [merge.netsim()].
#'
#' Network models additionally record the simulated networks, retrieved with
#' [get_network()]; the network statistics, retrieved with [get_nwstats()];
#' and the transmission tree, retrieved with [get_transmat()] and convertible
#' to a phylogeny with [as.phylo.transmat()]. Partnership histories may be
#' tracked over time with the cumulative edgelist and queried with
#' [get_partners()], [get_forward_reachable()], and [get_backward_reachable()].
#'
#' @section Learning Resources:
#' The package vignettes cover getting started (`Intro`) and then working with
#' model parameters, network objects, and nodal attributes and summary
#' statistics in extension models.
#'
#' Beyond the package documentation:
#'
#'  * The **Network Modeling for Epidemics** course materials provide a full
#'        treatment of network model theory, estimation, diagnostics, and
#'        epidemic simulation, followed by the extension API and advanced
#'        topics such as multi-layer networks and epidemics over fully observed
#'        networks, at <https://epimodel.github.io/sismid/>.
#'  * The **EpiModel Gallery** provides worked extension models, from adding a
#'        single compartment through multi-stage disease models with
#'        interventions and cost-effectiveness analysis, at
#'        <https://epimodel.github.io/EpiModel-Gallery/>.
#'  * The **EpiModel website** collects these resources and a suggested
#'        learning pathway, at <https://www.epimodel.org/>.
#'
#' @references
#' The EpiModel website is at <https://www.epimodel.org/>, and the source code
#' is at <https://github.com/EpiModel/EpiModel>. Bug reports and feature
#' requests are welcome.
#'
#' Our primary methods paper on EpiModel is published in the **Journal of
#' Statistical Software**. If you use EpiModel for any research or teaching
#' purposes, please cite this reference:
#'
#' Jenness SM, Goodreau SM, and Morris M. EpiModel: An R Package for
#' Mathematical Modeling of Infectious Disease over Networks. Journal of
#' Statistical Software. 2018; 84(8): 1-47. \doi{10.18637/jss.v084.i08}.
#'
#' EpiModel is the foundation for a set of extension packages maintained by our
#' group. `EpiModelHIV` (<https://github.com/EpiModel/EpiModelHIV>) models HIV
#' and bacterial sexually transmitted infections, and `EpiModelCOVID`
#' (<https://github.com/EpiModel/EpiModelCOVID>) models SARS-CoV-2. Simulation
#' on high-performance computing clusters is supported by `EpiModelHPC`
#' (<https://github.com/EpiModel/EpiModelHPC>).
#'
#' @name EpiModel-package
#' @aliases EpiModel
#' @import ergm network networkLite networkDynamic tergm ggplot2
#' @importFrom deSolve dede ode
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom graphics abline arrows boxplot legend lines mtext par plot points polygon text title
#' @importFrom grDevices col2rgb colorRampPalette rgb adjustcolor
#' @importFrom stats complete.cases quantile rbinom rgeom sd setNames simulate
#' @importFrom stats supsmu terms.formula update dgeom coef na.omit
#' @importFrom utils head tail packageVersion
#' @importFrom ape as.phylo collapse.singles
#' @importFrom lazyeval lazy_dots lazy_eval
#' @importFrom statnet.common trim_env set.control.class check.control.class NVL paste.and .Deprecate_once
#' @importFrom methods is
#' @importFrom tibble as_tibble tibble
#' @importFrom coda effectiveSize
#' @importFrom rlang .data
#'
#' @useDynLib EpiModel, .registration = TRUE
#'
#' @keywords package
#'
"_PACKAGE"
