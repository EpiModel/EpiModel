
#' @title Modules for Stochastic Network Models
#'
#' @description
#' Stochastic network models of infectious disease in EpiModel require
#' statistical modeling of networks, simulation of those networks forward
#' through time, and simulation of epidemic dynamics on top of those evolving
#' networks. The [netsim()] function handles both the network and
#' epidemic simulation tasks. Within this function are a series of modules that
#' initialize the simulation and then simulate new infections, recoveries, and
#' demographics on the network. Modules also handle the resimulation of the
#' network and some bookkeeping calculations for disease prevalence.
#'
#' Writing original network models that expand upon our "base" model set will
#' require modifying the existing modules or adding new modules to the workflow
#' in [netsim()]. The existing modules may be used as a template for
#' replacement or new modules.
#'
#' This help page provides an orientation to these module functions, in the
#' order in which they are used within [netsim()], to help guide users
#' in writing their own functions. These module functions are not shown
#' on the help index since they are not called directly by the end-user. To
#' understand these functions in more detail, review the separate help pages
#' listed below.
#'
#' @section Initialization Module:
#' This function sets up nodal attributes, like disease status, on the network
#' at the starting time step of disease simulation, \eqn{t_1}. For
#' multiple-simulation function calls, these are reset at the beginning of each
#' individual simulation.
#'
#'  * [initialize.net()]: sets up the main `netsim_dat` data
#'        structure used in the simulation, initializes which nodes are infected
#'        (via the initial conditions passed in [init.net()]), and
#'        simulates a first time step of the networks given the network model
#'        fit from [netest()].
#'
#'
#' @section Disease Status Modification Modules:
#' The main disease simulation occurs at each time step given the current state
#' of the network at that step. Infection of nodes is simulated as a function of
#' attributes of the nodes and the edges. Recovery of nodes is likewise
#' simulated as a function of nodal attributes of those infected nodes. These
#' functions also calculate summary flow measures such as disease incidence.
#'
#'  * [infection.net()]: simulates disease transmission given an
#'        edgelist of discordant partnerships by calculating the relevant
#'        transmission and act rates for each edge, and then updating the nodal
#'        attributes and summary statistics.
#'  * [recovery.net()]: simulates recovery from infection either
#'        to a lifelong immune state (for SIR models) or back to the susceptible
#'        state (for SIS models), as a function of the recovery rate parameters
#'        specified in [param.net()].
#'
#'
#' @section Demographic Modules:
#' Demographics such as arrival and departure processes are simulated at each
#' time step to update entries into and exits from the network. These are used
#' in epidemic models with network feedback, in which the network is resimulated
#' at each time step to account for the nodal changes affecting the edges.
#'
#'  * [departures.net()]: randomly simulates departure for nodes
#'        given their disease status (susceptible, infected, recovered), and
#'        their group-specific departure rates specified in
#'        [param.net()]. Departures involve deactivating nodes.
#'  * [arrivals.net()]: randomly simulates new arrivals into the
#'        network given the current population size and the arrival rate
#'        specified in the `a.rate` parameters. This involves adding new
#'        nodes into the network.
#'
#'
#' @section Network Resimulation Module:
#' In dependent network models, the network object is resimulated at each time
#' step to account for changes in the size of the network (changed through
#' entries and exits), and the disease status of the nodes.
#'
#'  * [resim_nets()]: resimulates the network object one time step
#'        forward given the set of formation and dissolution coefficients
#'        estimated in [netest()].
#'
#'
#' @section Bookkeeping Module:
#' Network simulations require bookkeeping at each time step to calculate the
#' summary epidemiological statistics used in the model output analysis.
#'
#'  * [prevalence.net()]: calculates the number in each disease
#'        state (susceptible, infected, recovered) at each time step for those
#'        active nodes in the network. If the `epi.by` control is used, it
#'        calculates these statistics by a set of specified nodal attributes.
#'  * [verbose.net()]: summarizes the current state of the
#'        simulation and prints this to the console.
#'
#'
#' @section One- & Two-Group Modules:
#' If epidemic `type` is supplied within [control.net()],
#' EpiModel defaults each of the base epidemic and demographic modules described
#' above (arrivals.FUN, departures.FUN, infection.FUN, recovery.FUN) to the
#' correct .net function based on variables passed to [param.net()]
#' (e.g. num.g2, denoting population size of group two, would select the
#' two-group variants of the aforementioned modules). Two-group modules are
#' denoted by a .2g affix (e.g., recovery.2g.net)
#'
#'
#' @name modules.net
#' @aliases modules.net
#'
NULL