# Disease Status Initialization Module for netsim

This function sets the initial disease status on the network given the
specified initial conditions.

## Usage

``` r
init_status.net(dat)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

## Value

The updated `netsim_dat` main list object.

## Details

This internal function sets, either randomly or deterministically, the
nodes that are infected at \\t_1\\, the starting time of network
simulations. If the number to be initially infected is passed, this
function sets the initial number infected based on the number specified,
either as a set of random draws from a binomial distribution or as the
exact number specified. In either case, the specific nodes infected are
a random sample from the network. In contrast, a set of specific nodes
may be infected by passing a vector containing the status of each node
to [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

For the initially infected nodes, this module sets the time of infection
as \\t_1\\, the starting time of network simulations. For models with
vital dynamics, the infection time for those initially infected nodes is
a random draw from an exponential distribution with the rate parameter
defined by the `di.rate` argument. For models without vital dynamics,
the infection time is a random draw from a uniform distribution of
integers with a minimum of 1 and a maximum of the number of time steps
in the model. In both cases, to set the infection times to be in the
past, these times are multiplied by -1, and 2 is added to allow for
possible infection times up until step 2, when the disease simulation
time loop starts.

## See also

This is an initialization module for
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).
